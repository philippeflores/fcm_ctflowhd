************************************************
* README GMM (other methods)
************************************************

*******************  Author  *******************
Philippe Flores (flores.philipe'at'gmail.com)
https://github.com/philippeflores/fcm_ctflowhd
************************************************

*****************  REFERENCES  *****************
McLachlan, G., and D. Peel,
Finite Mixture Models,
John Wiley & Sons, New York, 2000.
************************************************

*******************  Step #1  ******************
The fcm_GMM method takes pre-processed files as an input. Therefore, if you 
want to process a FCS file (i.e. 'file.fcs'), you must perform the 
pre-processing steps. To do this, please refer to the folder 
'./fcm_ctflowhd/pre_processingR/'.
************************************************

*******************  Step #2  ******************
To analyze a pre-processed FCS dataset (i.e. 'file_CCT.fcs'), this file 
must be inside the folder './fcm_ctflowhd/other_methods/data/'. This folder 
is not tracked in the git repository for memory concerns. Because of that, 
it may not exist for your first use of GMM. In this case, you may create 
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
Select the variables to build the GMM. At line 16, this is possible 
to manually type the variables. To help you in your choice, the variable 
names are stored in the variable 'strVariables'.

If you choose:
    >>> indVar = [];
only the fluorescence variables are going to be selected.

At line 21, there is the possibility to change the labels for the 
variables. To do so, uncomment the line 21 and run the section again. 
Follow the instructions to rename each label.
************************************************

*******************  Step #7  ****************** >>> Lines 23 to 28
This section permits to define the variable discrete intervals for the 
variables of  the analysis. There are several options that can be entered 
into the  function 'supportDist'. The goal of this function is to find the 
interval [a,b] for each variable in 'indVar'. First, the principle of 
'supportDist'  is to restrain a large interval [a0,b0] into a fitted 
interval [a,b]. Hence, a0 must be defined such that a0<a and b0 such that
b<b0.

For fluorescence variables, the values a0 and b0 are called 
'minFluo' and 'maxFluo'. Those values are set by default to [-5,10]. To 
change those values to 'newMin' and 'newMax', enter those parameters into 
the function:
    >>> newMin = -1; newMax = 6;
    >>> t = supportDist(X,indVar,100,'minFluo',newMin,'maxFluo',newMax);

For scatter variables, the values a0 and b0 are called 
'minScat' and 'maxScat' and are set by default to [0,1e+6]. To 
change those values, enter those parameters into the function:
    >>> newMin = 10000; newMax = 200000;
    >>> t = supportDist(X,indVar,100,'minScat',newMin,'maxScat',newMax);

It is possible to select how variable space is going to be chosen by 
choosing the parameter 'strMethod'. This parameter is set to 'none' but it 
can be change into 'fluo', 'label'  or 'man'.
    - 'none' (default):
For this option, all variables are considered fluorescence variables. This 
option is not suited if 'indVar' contains scatter variables. For all 
variables, the interval is set to [-0.5, 4.5].
    - 'fluo':
For this option, all variables are considered fluorescence variables. This 
option is not suited if 'indVar' contains scatter variables. For all 
variables, the interval is fitted from [minFluo, maxFluo]. 
    - 'label':
This option not suited if 'indVar' contains both scatter and fluorescence
variables. For fluorescence variables, the interval is fitted from 
[minFluo, maxFluo]. For scatter variables, the interval is fitted from 
[minScat, maxScat]. This function needs the 'strLabel' variable to detect 
if a variable is fluorescence or scatter. See below for more information.
    - 'man':
For this option, the intervals [a,b] can be chosen manually for each 
variable. This function needs the 'edges' variable as an input. See below 
for more information.

The option 'strLabel' is mandatory for the 'strMethod' option 'label'. This
variable should already exist at this point with the name 'strLabel'. 
Hence, it can be added as an input of supportDist:
    >>> t = supportDist(X,indVar,100,'strMethod','label','strLabel',strLabel);

The option 'edges' is mandatory for the 'strMethod' option 'man'. This
variable must be a variable of type cell of size 1xM where M is the number 
of variables in 'indVar'. The m-th element 'edges{m}' of edges must be a 
1x2 vector containing the values [a,b] of the corresponding m-th variable 
of 'indVar'. After creating such variable 'edges', it can be added as an 
input of supportDist:
    >>> edges = cell(1,M);
    >>> edges{1} = [...]; ... 
    >>> t = supportDist(X,indVar,100,'strMethod','man','edges',edges);
************************************************

*******************  Step #8  ****************** >>> Lines 29 to 35
This section builds the GMM. Before running this, the 
hyperparameter 'R' must be set up. It corresponds to the number of gaussian
mixtures for the model.
************************************************

*******************  Step #9  ****************** >>> Lines 36 to the end
This section plots the results of the  GMM. 

For large  datasets, it is possible to plot 1 cell every 'stepCloud' cells.
To change this optional parameter to 30 for example (default is 10), change 
for this:
    >>> stepCloud = 30;
    >>> plotGMMcyto(g,X,indVar,t,strLabel,'stepCloud',stepCloud)

You can also change the parameter 'strScreen' with the same method. You can
set this parameter with the following options:
    >>> ...,'strScreen','demiTop', ...) : half top of the screen (default),
    >>> ...,'strScreen','demiBot', ...) : half bottom of the screen,
    >>> ...,'strScreen','demiL', ...) : half left of the screen,
    >>> ...,'strScreen','demiR', ...) : half right of the screen,
    >>> ...,'strScreen','full', ...) : full screen.

By default, the function 'plot_GMM' plots a MxM grid of plots. On the 
diagonal, 1D-histograms are plotted. The global histogram and the 
histograms of each GMM component. On the rest of the grid, bivariate plots 
are featuring the cells along with the GMM component contours. If M is 
large, it can be usefull to plot a subset of these plots. To do so, one can
change the parameter 'indPlot' of the 'plot_GMM' function. 'indPlot' is a 
cell variable that contains the list of plots. For example, if one want to 
plot the 1D-histograms of variables 1 and 2 and the bivariate plots of 
variables 2 and 3, this must be typed:
    >>> indPlot = {1,2,[2,3]};
    >>> plotGMMcyto(g,X,indVar,t,strLabel,'indPlot',indPlot)

You can also change the parameter 'strXlim'. This parameter permits to 
choose on which intervals the variables are plotted. It is set by default 
to 'support' but can be changed to 'fluo' or 'auto':
    - 'support' (default):
If this option is chosen, variables are plotted with limits [a,b] defined 
by the function 'supportDist' (see Step 7).
    - 'fluo':
If this option is chosen, variables are plotted with limits [-0.5,4.5], it 
is not suited if scatter variables are plotted.
    - 'auto':
If this option is chosen, variables are plotted with the limits [a,b] 
defined automaticaly by MATLAB.

You can also change the parameter 'colMax' with the same method. You can
set this parameter with any number strictly positive. It will define the 
maximum number of plots along one column. This parameter will only affect 
figures where 'indPlot' is active.
************************************************

*******  Example given in the manuscript  ******
In the manuscript, figure 2.3 features an example of a GMM for a synthetic 
dataset.
************************************************
