************************************************
* README K-means (other methods)
************************************************

*******************  Author  *******************
Philippe Flores (flores.philipe'at'gmail.com)
https://github.com/philippeflores/fcm_ctflowhd
************************************************

****************  Disclaimer #1  ***************
This document is a tutorial to use the K-means method on raw FCM data. 
To use K-means on CTFlowHD results, please refer to the CTFlowHD folder.
************************************************

*****************  REFERENCES  *****************
James MacQueen and others.
Some methods for classification and analysis of multivariate observations.
Proceedings of the fifth Berkeley symposium on mathematical statistics and 
probability, 1967.
************************************************

*******************  Step #1  ******************
The K-means method takes pre-processed files as an input. Therefore, if you 
want to process a FCS file (i.e. 'file.fcs'), you must perform the 
pre-processing steps. To do this, please refer to the folder 
'./fcm_ctflowhd/pre_processingR/'.
************************************************

*******************  Step #2  ******************
To analyze a pre-processed FCS dataset (i.e. 'file_CCT.fcs'), this file 
must be inside the folder './fcm_ctflowhd/other_methods/data/'. This folder 
is not tracked in the git repository for memory concerns. Because of that, 
it may not exist for your first use of K-means. In this case, you may 
create the folder './fcm_ctflowhd/other_methods/data/'.

If the folder './fcm_ctflowhd/other_methods/data' already exists, put the 
file 'file_CCT.fcs' in the folder './fcm_ctflowhd/other_methods/data'.
************************************************

*******************  Step #3  ******************
Open MATLAB. Change your working directory to 
'./fcm_ctflowhd/other_methods'. Open the script file named 'fcm_Kmeans.m'.
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
Select the variables to apply the K-means method. At line 16, this is 
possible to manually type the variables that will be used for the method. 
To help you in your choice, the variable names are stored in the variable 
'strVariables'.

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
This section runs the K-means algorithm. Before running this, the 
hyperparameter 'K' must be set up. This parameter refers to the number of 
clusters.
************************************************

*******************  Step #9  ****************** Lines 36 to the end
This section plots the K-means clusters with a marginal visualization. 
The function 'plot_Kmeans' plots every 1D and 2D marginals possible from 
the set of variables in 'indVar'. This results in a MxM grid of plots. 
On the m-th diagonal of the grid, the 1D marginals of each K-means cluster 
is plotted as well as the empirical estimation of the corresponding 1D 
marginal. For other plots of indices (i,j), a 2D scatter plot of the cells 
is represented. After this, K-means clusters are plotted with contours 
of the 2D marginal distribution.

It is possible to plot only a portion of those lower-order marginals. This 
can be performed by entering 'indMarg' as an input of this function:
    >>> 'all' (default):
This option permits to plot every 1D and 2D marginals.
    >>> a variable of type cell:
If one inputs 'indMarg' as a variable of type cell containing the indices 
of the requested marginals, the function will only plots those marginals.
First, an 'indMarg' element can be a scalar corresponding to an index m
between 1 and M. This will results in the apparation of the m-th 1D 
marginal. Finally, an 'indMarg' element can be a 1x2 vector of indices 
(i,j) between 1 and M. This will results in the plot of the 2D marginal 
for variables i and j.
Example:
    >>> indMarg = {1,2,3,[1,2],[3,1],[2,3]};
    >>> plot_Kmeans(X,indVar,t,idx,C,strLabel,'indMarg',indMarg)
This will result in a figure that plots the 1D marginals of indices 1 2 and 
3 as well as 2D marginals (1,2) (3,1) and (2,3).

It is possible to separate 1D marginal plots from 2D marginal plots. To do 
so, you can enter the option 'boolSeparation' with the following:
    >>> indMarg = {1,2,3,[1,2],[3,1],[2,3]};
    >>> plot_Kmeans(X,indVar,t,idx,C,strLabel,'indMarg',indMarg, ...
            'boolSeparation',1)
This default choice is separating the two types of plots.
    >>> indMarg = {1,2,3,[1,2],[3,1],[2,3]};
    >>> plot_Kmeans(X,indVar,t,idx,C,strLabel,'indMarg',indMarg, ...
            'boolSeparation',0)
This default choice is not separating the two types of plots. However, 1D 
marginals will always be plotted before 2D marginals. This means that it is
not possible to plot alternatively a 1D marginal then a 2D and after coming 
back to a 1D marginal.

You can change the parameter 'strScreen' that permits to change the 
position and size of the figure. You can set this parameter with the 
following options:
    >>> ... ,'strScreen','halfTop', ...) : half top of the screen,
    >>> ... ,'strScreen','halfBot', ...) : half bottom of the screen,
    >>> ... ,'strScreen','halfL', ...) : half left of the screen,
    >>> ... ,'strScreen','halfR', ...) : half right of the screen,
    >>> ... ,'strScreen','full', ...) : full screen (default).

You can also change the parameter 'colMax' with the same method. You can
set this parameter with any number strictly positive. It will define the 
maximum number of plots along one column. This option will be ignored if 
all possible lower-order marginals are plotted.
    >>> indMarg = ...;
    >>> plot_Kmeans(X,indVar,t,idx,C,strLabel,'indMarg',indMarg,'colMax',2)

For large datasets, it is possible to plot 1 cell every 'stepCloud' cells. 
You can change this optional parameter to 30 (default is 50) for example:
    >>> plot_Kmeans(X,indVar,t,idx,C,strLabel,'stepCloud',30)
************************************************

*******  Example given in the manuscript  ******
In the thesis manuscript of the author, figure 5.15 features an example of 
one 2D marginal plot of Kmeans results for a flow cytometry controlled 
dataset.
************************************************
