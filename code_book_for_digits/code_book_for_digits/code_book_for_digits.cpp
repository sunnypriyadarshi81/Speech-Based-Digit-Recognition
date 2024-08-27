// code_book_for_digits.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "math.h"
#include "string.h"
#include "float.h"
#include "stdlib.h"

#define p 12									//NO of cepstrum coefficient 
#define DELTA 0.0001							//used in K-means algo for (current_distortion - prev_distortion)>DELTA
#define EPS 0.03								//used for spliting of codebook vectors (1+EPS) and (1-EPS) 	
#define M 32									//codebook size  (MxP) i.e 32x12
#define N 5										//number of states
#define T 160									//size of observation sequence


int flag=0,change=0,type=1;						//for selecting option 
int u_ind_row=0;								//used in durbins for row index of universe_matrix
int cnt_universe=0;								//size of universe
char fileName[100],universe[100],fileNo[100];

//durbins variables
long double alpha_durbins[p+1][p+1],E[p+1],k[p+1],C[p+1];
long double sine_raised_window[12];				//Liftering window used to  control the sensitivity and it amplifies the middle cepstrum values

							/* cnt_universe used in durbins were as size_of_uni used in read_universe() function */
//LBG_algo
int index_array[M]={0},minimum_dist=0,size_of_uni=0;
long double universe_matrix[30000][12], codebook[M+1][12],copyCodeBook[M+1][12],mappingbook[M+1][12];
long double tokhuraWt[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
int CBsize = 0;


//HMM algo
long double A[N+1][N+1],B[N+1][M+1],pi[N+1] = {0, 1.0, 0, 0, 0, 0},pstar=0,prev_pstar=0;
long double alpha[T+1][N+1], beta[T+1][N+1],gamma[T+1][N+1],Delta[T+1][N+1];
int ob[T+1],Psi[T+1][N+1],Q[T+1];				// psi[] tracke the max probability state sequence and from this matrix Q[] stores the Best sequence path 
long double Abar[N+1][N+1] = {0};
long double Bbar[N+1][M+1] = {0};
long double Pibar[N+1] = {0, 1.0, 0, 0, 0, 0};
long double xi[T+1][N+1][M+1],prob=0;
long double const threshold = 1e-30;		   //Min threshold to be assigned to zero values in matrix B

//for finding observation sequence
long double single_utter_cis[T+1][p+1]={0};
int index_utter=0;

//training functions
void training_model();
void read_all_hmm_files();
void create_observation_sequence(FILE *ptr,int digitNo, int utterNo);
void store_model(FILE *ptr_model_A,FILE *ptr_model_B);
long double avg_A[N+1][N+1],avg_B[N+1][M+1];

// Testing function
void testing_model();
void load_A_and_B(int k);
//   HMM Functions
void print_all_values();
void calc_xi();
void cal_gamma();
void cal_alpha();
void cal_beta();
void viterbi();
void reevaluate_model_parameters();
void load_calculated_model();

void durbins(long double R[],FILE *pt)
{
	
	E[0] = R[0];
	for(int i=1;i<=p;i++)
    {
        k[i] = R[i];
        for(int j=1;j<i;j++)
        {
            k[i]=k[i]-alpha_durbins[j][i-1]*R[i-j];
        }	
        
        k[i] = k[i]/E[i-1];
        alpha_durbins[i][i] = k[i];
        
        for(int j=1;j<i;j++)
        {
            alpha_durbins[j][i] = alpha_durbins[j][i-1] - alpha_durbins[i-j][i-1]*k[i];
        }
        E[i] = (1 - k[i]*k[i])*E[i-1];
    }

	
	for(int i=1;i<=12;i++)
	{
		sine_raised_window[i-1] = (1+ 6*sin((3.141*(long double)i)/12));						//raised_sine_window
	}

	for(int m=1;m<=p;m++)
	{
			C[m] = 0;
			for(int k=1;k<m;k++)
			{
				C[m]+= (long double)k/(long double)m*C[k]*alpha_durbins[m-k][p];
			}
		C[m]+=alpha_durbins[m][p];
		fprintf(pt,"%lf ",C[m]*sine_raised_window[m-1]);

		universe_matrix[u_ind_row][m-1] 			= C[m]*sine_raised_window[m-1];

		if(index_utter<160)
		{
			single_utter_cis[index_utter][m-1]   = C[m]*sine_raised_window[m-1];
		}
		//if(flag1<=5)printf("%lf  ",C[m]*sine_raised_window[m-1]);
	}
	u_ind_row++;

	if(index_utter<160)index_utter++;
	fprintf(pt,"\n");

	cnt_universe++;

	return;
}

//read universe from the TEXT file
void read_universe()
{
	FILE* fp;
	if(change==0 )fp = fopen("universe.txt","r");
	else fp = fopen("new_universe.txt","r");
	while(!feof(fp))
	{
		size_of_uni++;
		for(int j=0;j<12;j++)
		{
			long double x;
			fscanf(fp,"%lf",&x);
			universe_matrix[size_of_uni][j]=x;
		}
	}
}

long double tokhuraDistance(long double curr_region[],long double curr_universe[])
{
	long double ans =0;
	for(int i=0;i<12;i++)
	{
		long double d = curr_region[i] - curr_universe[i];
		ans += (tokhuraWt[i]*d*d);
	}
	return ans;
}

// function which updates the clusters centroid
void updateCodeBook_Kmeans()										
{
	for(int i=0;i<CBsize;i++)
	{
		for(int j=0;j<12;j++)
		{
				codebook[i][j] = mappingbook[i][j]/((long double)index_array[i]);
		}
	}
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<p;j++)mappingbook[i][j]=0;
	}
	for(int k=0;k<M;k++)index_array[k]=0;
	
	return;
}

// function to assign clusters to vectors of universe, it uses tokhura distance
void assign_cluster()
{
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<p;j++)mappingbook[i][j]=0;
	}
	for(int k=0;k<M;k++)index_array[k]=0;
	minimum_dist=0;
	int index=0;
	
	for(int i=0;i<size_of_uni;i++)													             //looping over all vectors of univers 
	{
		long double MinDis = LDBL_MAX,distance;
		
		long double current_universe[12];
		long double dis_array[M];
		for(int a=0;a<M;a++)dis_array[a]=0;

		for(int k=0;k<12;k++)current_universe[k] = universe_matrix[i][k];
		for(int j=0;j<CBsize;j++)
		{
			long double current_region[12];

			for(int y=0;y<12;y++)current_region[y] = codebook[j][y];

			dis_array[j] = tokhuraDistance(current_region,current_universe);
		}
		for(int k=0;k<CBsize;k++)
		{
			if(MinDis>dis_array[k])
			{
				MinDis = dis_array[k];
				index = k;
			}
		}
		//printf(" min distance taken = %lf\n",MinDis);
		minimum_dist=minimum_dist + MinDis;														//for calculating distortion
		for(int w=0;w<12;w++)
		{
			mappingbook[index][w]+=universe_matrix[i][w];
			
		}
		index_array[index]++;
	}

	return;
}


void k_means_algo()
{
	 
	assign_cluster();

	long double curr_dist = minimum_dist/((long double)size_of_uni);
	int round =1;
	long double old_dis =0.0;

	printf("Cycle %d has total distortion: %.8lf\n", round, curr_dist);
	while(abs(curr_dist- old_dis)> DELTA)
	{
		updateCodeBook_Kmeans();
		
		assign_cluster();
		old_dis = curr_dist;

		curr_dist = minimum_dist/((long double)size_of_uni);
		round++;

		printf("Cycle %d has total distortion: %.8lf\n", round, curr_dist);

	}
	return;
}


void updateCodebook_LBG()																	// function which updates the clusters centroid
{
	int k=0;
	for(int i=0; i<M; i++){
		
		for(int j=0;j<12;j++)
		{
			codebook[k][j] =   ((1 + EPS) * copyCodeBook[i][j]);
			codebook[k+1][j] = ((1 - EPS) * copyCodeBook[i][j]);
		}
		k=k+2;
	}
	
	return;
}

void LBGAlgo()
{
	int round = 1;
	while(CBsize <= M)
	{
        printf("Codebook of size: %d\n", CBsize);

        k_means_algo();
        printf("Codebook generated by LBG - \n");

        //printing codebook
        if(CBsize == 0){
            return;
        }
        printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
        printf("in LBG algo\n");
		for(int i=0;i<CBsize;i++){
            for(int j=0;j<12;j++){
                printf("| %10.7lf ", codebook[i][j]);
            }
            printf("|\n");
        }
        printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
		
        round++;
        if(CBsize >= M) return;
		
		for(int i=0;i<M;i++)
		{
			for(int j=0;j<12;j++)
			{
				copyCodeBook[i][j] = codebook[i][j];										//copying the codebook
				codebook[i][j] = 0.0;														//clear codebook
			}
		}
		
        CBsize *=2;
        updateCodebook_LBG();
    }
}

//read codebook from "codebook.txt"
void read_codebook()
{
	FILE* pt_cb;
	if(change==0)pt_cb = fopen("codebook.txt", "r");
	else pt_cb = fopen("new_codebook.txt", "r");
	for(int i=0;i<M;i++) 
	{
		for(int j=0;j<p;j++)
		{
			long double x;
			fscanf(pt_cb, "%lf", &x);
			codebook[i][j] = x;
			// printf("%.10e ", codebook[i][j]);
			
		}
		// printf("\n");
		i++;
	}
	fclose(pt_cb);
}


//it's create livetesting file and create its observation sequence and finds its probability with all the digits
void live_testing()
{
	char input[100];
	char file[100];
	system("\"Recording_Module\\Recording_Module.exe\" 2 livetesting.wav livetesting.txt");
	sprintf(input,"livetesting.txt");

	FILE *read = fopen(input,"r");
	index_utter=0;
	if(read==NULL)
	{
		printf("can't read the file\n");
		return;
	}
	for(int k=0;k<=160;k++)ob[k]=0;

	create_observation_sequence(read,9,1);
	// printf(" T size / utterance =  %d ",index_utter);
	long double max_prob=DBL_MIN;
	int ind=0;
	FILE* pq;
	pq = fopen("alpha_in_testing_with_all_model.txt","w+");
	long double prob_array[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	for(int k=0;k<=9;k++) 
	{
		load_A_and_B(k);
		cal_alpha();
		
		for(int t=1;t<=index_utter;t++)
		{
			for(int w=1;w<=N;w++)
			{
				fprintf(pq,"%.32e ",alpha[t][w]);
			}
			fprintf(pq,"\n");
		}
		fprintf(pq,"\n\n");
		prob_array[k] = prob;
		if(max_prob<prob){
			max_prob = prob;
			ind =k;
		}
	}
	fclose(pq);
	for(int k=0;k<=9;k++)printf("digit =  %d , probability = %.32e\n",k,prob_array[k]);
	printf("  digit detect =  %d \n",ind);
	
}


//----------------------------------------------------main function ----------------------------------------------------------------------------------------------------
int _tmain(int argc, _TCHAR* argv[])
{
	printf(" want to create universe 1-yes and 0-No : ");
	scanf("%d",&flag);
	if(flag==1)
	{
		FILE* pt;
		pt = fopen("new_universe.txt","w+");
		for(int i=0;i<=9;i++)
		{
			for(int j=1;j<=25;j++)
			{
				if(change==0)strcpy(fileName,"../../dataset/234101051_E_");							// ../../dataset/234101051_E_0_1.txt
				else  strcpy(fileName,"../../my_digit/234101051_");									// ../../my_digit/0/234101051_0_E_1.txt
				sprintf(fileNo,"%d",i);
				strcat(fileName,fileNo);
				if(change==0)strcat(fileName,"_");
				else strcat(fileName,"_E_");
				sprintf(fileNo,"%d",j);
				strcat(fileName,fileNo);
				strcat(fileName,".txt");
				printf("%s\n", fileName);

				FILE *ptr;
				ptr = fopen(fileName,"r");

				if (ptr == NULL) {
					printf("Error while opening the file");
					exit(EXIT_FAILURE);
				}

				// Skipping the headers in the file

   				 if(change == 1)
				{
					char ignore[5000];
					for (int i = 0; i < 5; i++)fgets(ignore, sizeof(ignore),ptr);
				}
				
				long double normalized_data[40000];
				long double data[322],max=0;
				int count=0;

				while(!feof(ptr))
				{
					long double x;
					fscanf(ptr,"%lf",&x);
					if(abs(max) <abs(x))max=x;
					// if(count<500)dc_mean += x;												
					normalized_data[count]=x;												// storing the files data in normalized_data array
					count++;																//count will tell the no. of data's in the file
				}
				printf("no. of data = %d\n",count);

				for(int j=0;j<count;j++)												    //normalization					
				{
					normalized_data[j] = (normalized_data[j]*5000)/max;
				}

				long double wd;
				int file_index,temp=1,x=0;	
				
				int ind = 0;																
				for (int y = 0; y < count; y++) {
					if (x == 320)
						x = 0;
					wd = (0.54 - 0.46*cos((2*22/7*(long double)x)/319.0));
					x++;
					normalized_data[ind] = normalized_data[ind]*wd;
				}

				while (ind + 320 < count)
				{
						long double R[p+1];
						for(int w=0;w<=p;w++)
						{
							long double sum=0.0;
							for(int t=ind;t<ind+320-w;t++)
							{
								sum = sum+normalized_data[t]*normalized_data[t+w];
							}
							R[w] = sum;
						}
						ind += 80;
						durbins(R,pt);															//calling durbins function
					
				}
				//printf("data in this cis - %d\n",cnt_universe);
			}
			//printf("  \n data in this digit is  %d\n",cnt_universe);
		}
		fclose(pt);
		size_of_uni = cnt_universe;
	}	
	/*------------------------------ creating codebook from the universe.txt file using LBG -----------------------------------------*/
	if(flag==0)
	{
		read_universe();
		cnt_universe = size_of_uni;
	}
	printf("size of unniverse = %d  or  %d \n",cnt_universe,size_of_uni);
	

	
	if(flag==1)
	{
		//FindCentroidOfUniverse;
		for(int i=0;i<cnt_universe;i++)
		{ 
			for(int j=0;j<12;j++)
			{
				//printf("%lf  ",universe_matrix[i][j]);
				codebook[0][j]+=universe_matrix[i][j];											//finding avg of each column of the universe matrix
			}
		}
		for(int k = 0; k<12; k++)
		{
			codebook[0][k] /= (long double)cnt_universe;										//taking average value by dividing the total number of vectors of universe
			//printf("codebook[0][%d] =  %lf  ",k,codebook[0][k]);     
		}
		CBsize=1;
		LBGAlgo();

		FILE* cb_ptr;
		if(change==0)cb_ptr = fopen("codebook.txt","w+");
		else cb_ptr = fopen("new_codebook.txt","w+");
		printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			printf("FINAL\n");
			for(int i=0;i<CBsize;i++){
				for(int j=0;j<12;j++){
					//printf("| %.32e ", codebook[i][j]);
					fprintf(cb_ptr,"%.32e ",codebook[i][j]);
				}
				fprintf(cb_ptr,"\n");
				//printf("|\n");
			}
		printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
		fclose(cb_ptr);
		
	}
	else
	{
		read_codebook();
	}

	
	printf("Select type \n1.Manual testing \n2.live testing \n3.Train Model\nExit\n");
	scanf("%d",&type);
	switch(type)
	{
	case 1:testing_model();
		break;
	case 2:live_testing();
		break;
	case 3:training_model();
		break;
	case 4:return 0;
	}

	return 0;
}


void training_model()
{
	
	for(int i=0;i<=9;i++)
	{
		char avg_name_A[100],avg_num_A[100],avg_name_B[100],avg_num_B[100];
		FILE *avg_model_A,*avg_model_B;
		strcpy(avg_name_A,"../lamda_new/final_lamda/model_A_");						// ../lamda/final_lamda/model_A_0.txt
		strcpy(avg_name_B,"../lamda_new/final_lamda/model_B_");						// ../lamda/final_lamda/model_B_0.txt
		// strcpy(avg_name_A,"../my_lamda/final_lamda/model_A_");						// ../lamda/final_lamda/model_A_0.txt
		// strcpy(avg_name_B,"../my_lamda/final_lamda/model_B_");						// ../lamda/final_lamda/model_B_0.txt
		sprintf(avg_num_A,"%d",i);
		sprintf(avg_num_B,"%d",i);
		strcat(avg_name_A,avg_num_A);
		strcat(avg_name_B,avg_num_B);
		strcat(avg_name_A,".txt");
		strcat(avg_name_B,".txt");
		printf("%s\n", avg_name_A);
		avg_model_A = fopen(avg_name_A,"w+");
		avg_model_B = fopen(avg_name_B,"w+");

		for(int x=1;x<=N;x++){
			for(int y=1;y<=N;y++)avg_A[x][y]=0.0;
		}
		for(int x=1;x<=N;x++){
			for(int y=1;y<=M;y++)avg_B[x][y] = 0.0;
		}

		for(int j=1;j<=25;j++)
		{
			strcpy(fileName,"../../dataset/234101051_E_");								// ../../dataset/234101051_E_0_1.txt
			// strcpy(fileName,"../../my_digit/234101051_");							// ../../my_digit/0/234101051_0_E_1.txt
			sprintf(fileNo,"%d",i);
			strcat(fileName,fileNo);
			strcat(fileName,"_");
			// strcat(fileName,"_E_");
			sprintf(fileNo,"%d",j);
			strcat(fileName,fileNo);
			strcat(fileName,".txt");
			printf("%s\n", fileName);

			FILE *ptr;
			ptr = fopen(fileName,"r");

			if (ptr == NULL)
			{
				printf("Error while opening the file");
				exit(EXIT_FAILURE);
			}

			index_utter=0;
			create_observation_sequence(ptr,i,j);
			read_all_hmm_files();
		
			
			int itr=1;
			do
			{	prev_pstar = pstar;
				cal_alpha();
				cal_beta();
				cal_gamma();
				viterbi();
				calc_xi(); 
				reevaluate_model_parameters();
				if(pstar>prev_pstar)load_calculated_model();
				itr++;
			} while (pstar>prev_pstar && itr<100);

			printf("final pstar =  %.32e   \n",pstar);
			
			char fileModel_A[100],filenum_A[10];
			char fileModel_B[100],filenum_B[10];
			strcpy(fileModel_A,"../lamda_new/234101051_E_A_");								//../lamda/234101051_E_A_0_1.txt
			strcpy(fileModel_B,"../lamda_new/234101051_E_B_");								//../lamda/234101051_E_B_0_1.txt
			// strcpy(fileModel_A,"../my_lamda/234101051_E_A_");								//../lamda/234101051_E_A_0_1.txt
			// strcpy(fileModel_B,"../my_lamda/234101051_E_B_");								//../lamda/234101051_E_B_0_1.txt
			sprintf(filenum_A,"%d",i);
			sprintf(filenum_B,"%d",i);
			strcat(fileModel_A,filenum_A);
			strcat(fileModel_B,filenum_B);
			strcat(fileModel_A,"_");
			strcat(fileModel_B,"_");
			sprintf(filenum_A,"%d",j);
			sprintf(filenum_B,"%d",j);
			strcat(fileModel_A,filenum_A);
			strcat(fileModel_B,filenum_B);
			strcat(fileModel_A,".txt");
			strcat(fileModel_B,".txt");
			FILE *ptr_model_A,*ptr_model_B;
			ptr_model_A = fopen(fileModel_A,"w+");
			ptr_model_B = fopen(fileModel_B,"w+");

			store_model(ptr_model_A,ptr_model_B);
			fclose(ptr_model_A);
			fclose(ptr_model_B);
			
		}

		for(int x=1;x<=N;x++){
			for(int y=1;y<=N;y++){
				A[x][y] = avg_A[x][y]/25.0;
			}
		}

		for(int x=1;x<=N;x++){
			for(int y=1;y<=M;y++){
				B[x][y] = avg_B[x][y]/25.0;
			}
		}
		
		store_model(avg_model_A,avg_model_B);
		fclose(avg_model_A);
		fclose(avg_model_B);
	}

}

//this test the input file by creating observation sequence and finding its probability with all digits
void testing_model()
{
	for(int i=0;i<=9;i++)
	{
		
		long double prob_array[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		for(int j=21;j<=30;j++)
		{
			// for(int c=1;c<=N;c++)pi[c]=0.0;
			if(change==0)strcpy(fileName,"../../dataset/234101051_E_");					// ../../dataset/234101051_E_0_1.txt
			else strcpy(fileName,"../../my_digit/234101051_");							// ../../my_digit/0/234101051_0_E_1.txt
			sprintf(fileNo,"%d",i);
			strcat(fileName,fileNo);
			if(change==0)strcat(fileName,"_");
			else strcat(fileName,"_E_");
			sprintf(fileNo,"%d",j);
			strcat(fileName,fileNo);
			strcat(fileName,".txt");
			printf("-----------------------------------------------------------------------------------------------------------\n");
			printf("%s\n", fileName);

			FILE *ptr;
			ptr = fopen(fileName,"r");

			if(change == 1)
			{
				char ignore[5000];
				for (int i = 0; i < 5; i++)fgets(ignore, sizeof(ignore),ptr);
			}
			index_utter=0;

			for(int k=0;k<=160;k++)ob[k]=0;

			create_observation_sequence(ptr,i,j);										//create observation sequence and set utter_index
			// printf(" T size / utterance =  %d ",index_utter);
			long double max_prob=DBL_MIN;
			int ind=0;
			FILE* pq;
			pq = fopen("alpha_in_testing_with_all_model.txt","w+");
			for(int k=0;k<=9;k++) 
			{
				load_A_and_B(k);
				cal_alpha();
				
				for(int t=1;t<=index_utter;t++)
				{
					for(int w=1;w<=N;w++)
					{
						fprintf(pq,"%.32e ",alpha[t][w]);
					}
					fprintf(pq,"\n");
				}
				fprintf(pq,"\n\n");
				prob_array[k] = prob;
				if(max_prob<prob){
					max_prob = prob;
					ind =k;
				}
			}
			fclose(pq);
			for(int k=0;k<=9;k++)printf("digit =  %d , probability = %.32e\n",k,prob_array[k]);
			printf("  digit detect =  %d \n",ind);

		}
	}
}


//load the A and B matrix
void load_A_and_B(int k)
{
	FILE *fp_A,*fp_B;
	char name_A[100],num_A[10];
	char name_B[100],num_B[10];
	if(change==0)strcpy(name_A,"../lamda_new/final_lamda/model_A_");
	else strcpy(name_A,"../my_lamda/final_lamda/model_A_");
	sprintf(num_A,"%d",k);
	strcat(name_A,num_A);
	strcat(name_A,".txt");

	if(change==0)strcpy(name_B,"../lamda_new/final_lamda/model_B_");
	else strcpy(name_B,"../my_lamda/final_lamda/model_B_");
	sprintf(num_B,"%d",k);
	strcat(name_B,num_B);
	strcat(name_B,".txt");
	
	FILE* ch;
	ch = fopen("loading_A_and_B_values.txt","w+");
	fp_A = fopen(name_A,"r");
	fp_B = fopen(name_B,"r");
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++){
			fscanf(fp_A,"%lf",&A[i][j]);
			fprintf(ch,"%.32e ",A[i][j]);
		}
		fprintf(ch,"\n");
	}
	
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++){
			fscanf(fp_B,"%lf",&B[i][j]);
			fprintf(ch,"%.32e ",B[i][j]);
		}
		fprintf(ch,"\n");
	}
	fclose(ch);
	fclose(fp_A);
	fclose(fp_B);

}

//store the A and B matrix in TEXT file
void store_model(FILE *ptr_model_A,FILE *ptr_model_B)
{
	// printf("storing values\n");
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++){
			fprintf(ptr_model_A,"%.32e ",A[i][j]);
			avg_A[i][j] += A[i][j];
		}
		fprintf(ptr_model_A,"\n");
	}

	for(int i=1;i<=N;i++){
		for(int j=1;j<=M;j++){
			fprintf(ptr_model_B,"%.32e ",B[i][j]);
			avg_B[i][j] += B[i][j];
		}
		fprintf(ptr_model_B,"\n");
	}
}

//Read the A.txt  , B.txt and pi.txt files
void read_all_hmm_files()
{
	FILE *fpA,*fpB,*fpPI,*fpO;
	fpA = fopen("A.txt","r");
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			fscanf(fpA,"%lf",&A[i][j]);
		}
	}
	fclose(fpA);

	fpB = fopen("B.txt","r");
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++)
		{
			fscanf(fpB,"%lf",&B[i][j]);
		}
	}
	fclose(fpB);

	fpPI = fopen("pi.txt","r");
	for(int i=1;i<=N;i++)fscanf(fpPI,"%lf",&pi[i]);
	fclose(fpPI);
}


// create observation sequence and store in ob[] matrix and upadte the index_utter values
//	which stores the #obs seq. store in the array
void create_observation_sequence(FILE *ptr,int digitNo, int utterNo)
{
	FILE obs_seq[100],digi[100];
	strcpy(fileName,"../obs_sequence/seq_");									// ../../obs_sequence/seq_0_1.txt
	sprintf(fileNo,"%d",digitNo);
	strcat(fileName,fileNo);
	strcat(fileName,"_");
	sprintf(fileNo,"%d",utterNo);
	strcat(fileName,fileNo);
	strcat(fileName,".txt");

	FILE *store_cis,*store_obs;
	store_obs = fopen(fileName,"w+");
	store_cis = fopen("cis_of_single_utterance.txt","w+");
	long double normalized_data[40000];
	long double data[322],max=0,y;
	int count=0,t_count=0;

	while(!feof(ptr))
	{
		t_count++;
		if(type==2 && t_count<10000)continue;
		fscanf(ptr,"%lf",&y);
		if(abs(max) <abs(y))max=y;
		normalized_data[count]=y;												// storing the files data in normalized_data array
		count++;																//count will tell the no. of data's in the file
	}

	for(int j=0;j<count;j++)												    //normalization					
	{ 
		normalized_data[j] = (normalized_data[j]*5000)/max;
	}
	long double wd;
	int file_index,temp=1,x=0;													//temp will check 320 size
	

	index_utter=0;
	for(int j=0;j<=T;j++)
	{
		for(int s=0;s<12;s++)single_utter_cis[j][s]=0;
	}
	int ind = 0;																
	for (int y = 0; y < count; y++) {
		if (x == 320)
			x = 0;
		wd = (0.54 - 0.46*cos((2*22/7*(long double)x)/319.0));
		x++;
		normalized_data[ind] = normalized_data[ind]*wd;
	}

	while (ind + 320 < count)
	{
		long double R[p+1];
		for(int w=0;w<=p;w++)
		{
			long double sum=0.0;
			for(int t=ind;t<ind+320-w;t++)
			{
				sum = sum+normalized_data[t]*normalized_data[t+w];
			}
			R[w] = sum;
		}
		ind += 80;
		durbins(R,store_cis);																			//calling durbins function
		
	}
	

	fclose(store_cis);

	for(int i=0;i<index_utter;i++)													                     //looping over all vectors of univers 
	{
		long double MinDis = LDBL_MAX,distance;
		int index=1;
		long double current_utterance[12];
		long double dis_array[M];

		for(int a=0;a<M;a++)dis_array[a]=0;

		for(int k=0;k<12;k++)current_utterance[k] = single_utter_cis[i][k];
		for(int j=0;j<M;j++)
		{
			long double current_region[12];

			for(int y=0;y<12;y++)current_region[y] = codebook[j][y];

			dis_array[j] = tokhuraDistance(current_region,current_utterance);
			
		}
		for(int k=0;k<M;k++)
		{
			if(MinDis>dis_array[k])
			{
				MinDis = dis_array[k];
				index = k;
			}
		}
		ob[i+1] = index+1;  
		// printf("obs - ");                                                                     //observation sequence range 1-32;
		// printf("%d ",ob[i+1]);
		fprintf(store_obs,"%d ",ob[i+1]);
	}

	fclose(store_obs);
}



void calc_xi()
{
	long double sum=0;

	for(int t=1;t<= index_utter-1 ;t++)															//T-1
	{
		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=N;j++)
			{
				sum=0;
				for(int x=1;x<=N;x++)
				{
					for(int y=1;y<=N;y++)
					{
						sum+=alpha[t][x]*A[x][y]*B[y][ob[t+1]]*beta[t+1][y];
					}
				}
				xi[t][i][j] = alpha[t][i]*A[i][j]*B[j][ob[t+1]]*beta[t+1][j]/sum;
			}
		}
	}
}


//calculating gamma values using xita
void cal_gamma()
{
	long double sum=0;
	for(int t=1;t<=index_utter;t++)															//T
	{
		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=N;j++)
			{
				sum+=alpha[t][j]*beta[t][j];
			}
			gamma[t][i] = (alpha[t][i]*beta[t][i])/sum;
			sum=0;			
		}
		
	}
	return;
}


//Calculation of Beta variable
void cal_beta()
{
	//initialization_Beta();
	for(int i=1;i<=index_utter;i++)																		//T
	{
		for(int j=1;j<=N;j++)beta[i][j]=0;
	}

	for(int i=1;i<=N;i++)
		beta[index_utter][i] = 1.0;


	//induction_Beta();
	for(int t=index_utter-1;t>0;t--)																	//T-1
	{
		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=N;j++)
			{
				beta[t][i]+= A[i][j]*B[j][ob[t+1]]*beta[t+1][j];
			}
		}
	}
	return;
}

//Calculation of alpha variable to find the solution of problem number 1.
void cal_alpha()
{   
	
	//initialization_alpha();
	for(int i=1;i<=N;i++)
		alpha[1][i] = pi[i]*B[i][ob[1]];
	
	
	//induction_alpha();
	for(int t=1;t<=index_utter-1;t++)
	{
		long double sum=0;
		for(int j=1;j<=N;j++)
		{
			for(int i=1;i<=N;i++)
			{
				sum+= alpha[t][i]*A[i][j];
			}
			alpha[t+1][j] = sum*B[j][ob[t+1]];
			sum = 0;
		}
	}

	//termination_alpha();
	prob=0;
	for(int i=1;i<=N;i++)
	{
		prob+=alpha[index_utter][i];
	}
	//  printf("\nmy problality is = %.32elf\n",prob);

}


//finding best state sequence for an observation sequence
void viterbi()
{
	//initialization
	for(int i=1;i<=N;i++)
	{
		Delta[1][i] = pi[i]*B[i][ob[1]];
		Psi[i][1] = 0;
	}

	//recursion
	for(int t=2;t<=index_utter;t++)																		//T
	{
		for(int j=1;j<=N;j++)
		{
			long double max=0,index=0,ti;

			for(int i=1;i<=N;i++)
			{
				ti = Delta[t-1][i]*A[i][j];
				if(max< ti)
				{
					max= Delta[t-1][i]*A[i][j];
					index = i;
				}
			}
			Delta[t][j] = max*B[j][ob[t]];					//delta is the highest probability along a single path at time t
			Psi[t][j] = index;
		}
	}


	//termination
	long double max = 0;
    for(int i=1; i<=N; i++){
        if(Delta[index_utter][i] > max) {														//T
			max = Delta[index_utter][i];
			Q[index_utter] = i;
		}

        pstar = max;
    }
	
    for(int t = index_utter-1; t>0; t--){														//T-1
        Q[t] = Psi[t+1][Q[t+1]];
    }
}


//loading the model parameters with new calculated values
void load_calculated_model(){
	int i, j;

	//updating Abar
	for(i=1;i<=N;i++)
	{
		long double sum=0.0,max=0;
		int max_index=0;
		for(j=1;j<=N;j++)
		{
			sum += Abar[i][j];
			if(max<Abar[i][j])
			{
				max = Abar[i][j];
				max_index = j;
			}
		}
		// printf("sum of A[%d] = %.32e value going to add = %.32e  \n",i,sum,(1.0 - sum));
		if(sum!=1.0)
		{
			Abar[i][max_index] += (1.0 - sum);
		}
	}
	
	//updating Bbar
	for(int i=1;i<=N;i++)
	{
		long double sum=0.0,max=0;
		int max_index=0;
		for(j=1;j<=M;j++)
		{

			if(Bbar[i][j]<threshold)Bbar[i][j]=threshold;					///------------------
			sum += Bbar[i][j];
			if(max<Bbar[i][j])
			{
				max = Bbar[i][j];
				max_index = j;
			}
		}
		// printf("sum of B[%d] = %.32e  value going to add = %.32e \n",i,sum,(1.0 - sum));
		if(sum!=1.0)
		{
			Bbar[i][max_index] += (1.0 - sum);
		}
	}

	for(i=1;i<=N;i++){
		pi[i]=Pibar[i];
	}
	
	for(i=1;i<=N;i++){
		for(j=1;j<=N;j++){
			A[i][j]= Abar[i][j];
		}
	}
	for(i=1;i<=N;i++){
		for(j=1;j<=M;j++){
			B[i][j] = Bbar[i][j];
		}
	}
}

void reevaluate_model_parameters(){
	int i, j, k, t;
	long double sum1=0 , sum2 =0;

	//Re-evaluating Pi
	for(i=1;i<=N;i++){
		Pibar[i] = gamma[1][i];
	}
	
	for(i = 1; i<=N; i++){
		for(j = 1; j <= N; j++){
			long double t1 = 0, t2 = 0;
			for(t = 1; t <= T-1; t++){
				t1 += xi[t][i][j];
				t2 += gamma[t][i];
			}
			Abar[i][j] = t1/t2;
		}
	}
	//Re-evaluating B

	for(j=1;j<=N;j++){
		
		for(k=1;k<=M;k++){
			sum1 =0 , sum2 =0;
			for(t=1;t<=T;t++){
				sum1 = sum1 + gamma[t][j];
				if(ob[t]==k){
				
					sum2 = sum2 + gamma[t][j];				
				}
			}
			Bbar[j][k] = sum2/sum1;
			
		}
		
	}
}
