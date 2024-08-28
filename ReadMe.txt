========================================================================
    CONSOLE APPLICATION : 234101051_DigitRecognition Project Overview
========================================================================

AppWizard has created this 234101051_DigitRecognition application for you.

This file contains a summary of what you will find in each of the files that
make up your 234101051_DigitRecognition application.


234101051_DigitRecognition.vcxproj
    This is the main project file for VC++ projects generated using an Application Wizard.
    It contains information about the version of Visual C++ that generated the file, and
    information about the platforms, configurations, and project features selected with the
    Application Wizard.

234101051_DigitRecognition.vcxproj.filters
    This is the filters file for VC++ projects generated using an Application Wizard. 
    It contains information about the association between the files in your project 
    and the filters. This association is used in the IDE to show grouping of files with
    similar extensions under a specific node (for e.g. ".cpp" files are associated with the
    "Source Files" filter).

234101051_DigitRecognition.cpp
    This is the main application source file.

/////////////////////////////////////////////////////////////////////////////
Other standard files:

StdAfx.h, StdAfx.cpp
    These files are used to build a precompiled header (PCH) file
    named 234101051_DigitRecognition.pch and a precompiled types file named StdAfx.obj.

/////////////////////////////////////////////////////////////////////////////
Other notes:

AppWizard uses "TODO:" comments to indicate parts of the source code you
should add to or customize.

/////////////////////////////////////////////////////////////////////////////


The whole code is the collection of the modular functions
=> There are major 7 functions which are used in the main function and those are listed below:
    -> readuniverse()
        - The function reads the the universe from the universe.txt file
    
    -> LBG()
        - the LBG algorithm it generates the codebook of size 32 using the universe file generated at above function.

    -> training_model()
        - Using the training file (i.e. first 25 utterances of the all digits) it trains the HMM Models and generates the final models for every digits
    
    -> testing_model()
        - It tests all the test files (i.e. 20 to 30 utterances of all the digits) and shows the accuracy along with the winning probability of all the digits
    
    -> liveTesting()
        - This function records the digit live and recognises the spoken digit
    

    -> readCodebook()
        - this reads the final codebook in the codebook array to generate the observation sequences for all the files

    ->create_obesrvation_sequence()
	  - this create observation sequence and store in ob[] array

=> the file structure is as follows:
    -> "dataset" folder is having the dataset used for training and testing
    -> "lamda_new" stores the  train model
    -> "obs_seq" store the observation sequence 
