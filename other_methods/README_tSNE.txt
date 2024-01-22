************************************************
* README t-SNE (other methods)
************************************************

*******************  Author  *******************
Philippe Flores (flores.philipe'at'gmail.com)
https://github.com/philippeflores/fcm_ctflowhd
************************************************

****************  Disclaimer #1  ***************
This document is a tutorial to use the t-SNE method on raw FCM data. To use 
t-SNE on CTFlowHD results, please refer to the CTFlowHD folder.
************************************************

*****************  REFERENCES  *****************
VAN DER MAATEN, Laurens and HINTON, Geoffrey,
Visualizing data using t-SNE.
Journal of machine learning research, 2008.
************************************************

*******************  Step #1  ******************
The t-SNE method takes pre-processed files as an input. Therefore, if you 
want to process a FCS file (i.e. 'file.fcs'), you must perform the 
pre-processing steps. To do this, please refer to the folder 
'./fcm_ctflowhd/pre_processingR/'.
************************************************

*******************  Step #2  ******************
To analyze a pre-processed FCS dataset (i.e. 'file_CCT.fcs'), this file 
must be inside the folder './fcm_ctflowhd/other_methods/data/'. This folder 
is not tracked in the git repository for memory concerns. Because of that, 
it may not exist for your first use of t-SNE. In this case, you may create 
the folder './fcm_ctflowhd/other_methods/data/'.

If the folder './fcm_ctflowhd/other_methods/data' already exists, put the 
file 'file_CCT.fcs' in the folder './fcm_ctflowhd/other_methods/data'.
************************************************

*******************  Step #3  ******************
Open MATLAB. Change your working directory to 
'./fcm_ctflowhd/other_methods'. Open the script file named 'fcm_tSNE.m'.
************************************************

*******************  Step #4  ****************** >>> Lines 1 to 7
Set up the MATLAB workspace. Run the first section without any changes. It 
will set up the workspace by clearing it, close any open windows and clear 
up the command window. Finally, it will add to the working directory the 
functions needed to run this script.
************************************************

*******************  Step #5  ****************** >>> Lines 8 to 14
Select the file to be analyzed. If you want to analyze the file 
'file_CCT.fcs', you can either type the following options:
    - './data/file_CCT.fcs'
    - 'data/file_CCT.fcs'
    - 'file_CCT.fcs'
    - './data/file_CCT'
    - 'data/file_CCT'
    - 'file_CCT'
If you select a file that is not in the './fcm_ctflowhd/other_methods/data'
directory, the function will ask you to type a valid filename (if this is 
possible).
************************************************

*******************  Step #6  ****************** >>> Lines 15 to 22
Select the variables to build the t-SNE map. At line 16, this is possible 
to manually type the variables that will be used for the t-SNE map. To help
you in your choice, the variable names are stored in the variable 
'strVariables'.

If you choose:
    >>> indVar = [];
only the fluorescence variables are going to be selected.

At line 21, there is the possibility to change the labels for the 
variables. To do so, uncomment the line 21 and run the section again. 
Follow the instructions to rename each label.
************************************************

*******************  Step #7  ****************** >>> Lines 23 to 29
This section runs the t-SNE algorithm. Before running this, the 
hyperparameter 'perplexity' must be set up. This parameter can be chosen 
between 2 and N where N represents the number of rows of X. The t-SNE 
algorithm run slower for high perplexities values. However, for higher 
values of perplexity, cells are more separated.
************************************************

*******************  Step #8  ****************** >>> Lines 30 to the end
This section plots one t-SNE map for each marker in 'indVar'. 

For large 
datasets, it is possible to plot 1 cell every 'stepCloud' cells. To change 
this optional parameter to 30 for example, change line 35 for this:
    >>> plotTsne(Y,indVar,X,strLabel,'stepCloud',30)

You can also change the parameter 'strScreen' with the same method. You can
set this parameter with the following options:
    >>> ...,'strScreen','demiTop', ...) : half top of the screen (default),
    >>> ...,'strScreen','demiBot', ...) : half bottom of the screen,
    >>> ...,'strScreen','demiL', ...) : half left of the screen,
    >>> ...,'strScreen','demiR', ...) : half right of the screen,
    >>> ...,'strScreen','full', ...) : full screen.

You can also change the parameter 'xlimMan' with the same method. You can
set this parameter with any array of size 1x2 [a b] where a<b. The colormap
limits for every plot will be set with [a,b]. Default is [-0.5 4.5] and is 
suited for fluorescence values. If you want colomap limits picked 
automaticly, choose 'auto' for this parameter.

You can also change the parameter 'colMax' with the same method. You can
set this parameter with any number strictly positive. It will define the 
maximum number of plots along one column.
************************************************

*******  Example given in the manuscript  ******
In the manuscript, figure 5.17 features an example of t-SNE applied to a 
4-fluorescence controlled dataset. Feel free to ask me by e-mail for this 
dataset.
************************************************
