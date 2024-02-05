************************************************
* README CTFlowHD_plot_Kmeans
************************************************

*******************  Author  *******************
Philippe Flores (flores.philipe@gmail.com)
https://github.com/philippeflores/fcm_ctflowhd
************************************************

*******************  Step #1  ******************
The script 'CTFlowHD_plot_Kmeans' is a visualization script. Therefore, it
requires that the PCTF3D method was applied prior to this script. Please 
refer to the script 'CTFlowHD_PCTF3D' if this is not the case.
************************************************

*******************  Step #2  ****************** Lines 1 to 9
As a first step, PCTF3D results must be loaded on your workspace. There are
two possible ways to do this.

First, you can load a PCTF3D result saved file from the directory 
'./fcm_CTFlowHD/save/'. In order to save PCTF3D results, please refer to 
the script 'CTFlowHD_PCTF3D' and its instructions 'README_PCTF3D'.

Secondly, you can run the script 'CTFlowHD_plot_Kmeans' right after running 
the script 'CTFlowHD_PCTF3D'. To do so, you have to skip the first two 
sections of the script 'CTFlowHD_plot_Kmeans' (Lines 1 to 9) as it clears 
your workspace.
************************************************

*******************  Step #3  ****************** Lines 10 to the end
This section plots results of PCTF3D with a K-means visualization. The 
K-means clusters are obtained from the rank-one terms of the CPD obtained 
by PCTF3D. The hyperparameter 'K' must be set up. It refers to the number 
of clusters.

To perform K-means on rank-one terms, the centroid of each component must 
be computed. To do so, it is possible to change the option 'strCentroid' 
of the function 'plot_NBMkMeans'. The two options for this parameter are 
'esp' and 'max':
    - 'esp' (default):
For this option, the centroid of a factor y{m}(:,r) is defined as the 
expected value of the 1D marginal y{m}(:,r).
    - 'max' :
For this option, the centroid of a factor y{m}(:,r) is defined as the 
centroid of the bin that has maximum probability y{m}(:,r).

You can change the parameter 'strScreen' that permits to change the 
position and size of the figure. You can set this parameter with the 
following options:
    >>> ...,'strScreen','halfTop', ...) : half top of the screen,
    >>> ...,'strScreen','halfBot', ...) : half bottom of the screen,
    >>> ...,'strScreen','halfL', ...) : half left of the screen,
    >>> ...,'strScreen','halfR', ...) : half right of the screen,
    >>> ...,'strScreen','full', ...) : full screen (default).

By default, the function 'plot_NBMkMeans' plots all possible plots of order
1 and 2. However, it is possible to plot only a portion of lower-order 
plots. This can be performed by entering *indMarg* as an input of this 
function:
    >>> 'all' (default):
This option permits to plot every 1D and 2D plots.
    >>> a variable of type cell:
If one inputs *indMarg* as a variable of type cell containing the indices 
of the requested plots, the function will only plots those. First, an 
*indMarg* element can be a scalar corresponding to an index m between 1 and 
*M*. This will results in the plot of the m-th 1D marginal. Finally, 
an *indMarg* element can be a 1x2 vector of indices (j,k) between 1 and 
*M*. This will results in the plot of the 2D plot of centroids for 
variables j and k. Example:
    >>> indMarg = {1,2,3,[1,2],[3,1],[2,3]};
    >>> plot_NBMkMeans(y,lambda,K,indVar,t,strLabel,'indMarg',indMarg)
This will result in a figure that plots the 1D marginals of indices 1 2 and 
3 as well as 2D plots for variables in (1,2) (3,1) and (2,3).

It is possible to separate 1D marginal plots from 2D plots. To do so, you 
can enter the option 'boolSeparation' with the following:
    >>> plot_NBMkMeans(y,lambda,K,indVar,t,strLabel,'indMarg',indMarg, ...
            'boolSeparation',0)
The default choice is separating the two types of plots.
    >>> plot_NBMkMeans(y,lambda,K,indVar,t,strLabel,'indMarg',indMarg, ...
            'boolSeparation',1)
This default choice is not separating the two types of plots. However, 1D 
marginals will always be plotted before 2D plots. This means that it is
not possible to plot alternatively a 1D marginal then a 2D plot and after 
coming back to a 1D marginal.

You can also change the parameter 'colMax' with the same method. You can
set this parameter with any strictly positive number. It will define the 
maximum number of plots along one column. This option will be ignored if 
all possible plots are considered.
    >>> plot_NBMkMeans(y,lambda,K,indVar,t,strLabel,'indMarg',indMarg, ...
            'colMax',3)
************************************************

*******  Example given in the manuscript  ******
In the thesis manuscript of the author, figure 5.15 features an example of 
one K-means plot for a flow cytometry controlled dataset.
************************************************
