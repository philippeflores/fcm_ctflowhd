************************************************
* README CTFlowHD_plot_Marg
************************************************

*******************  Author  *******************
Philippe Flores (flores.philipe@gmail.com)
https://github.com/philippeflores/fcm_ctflowhd
************************************************

*******************  Step #1  ******************
The script 'CTFlowHD_plot_Marg' is a visualization script. Therefore, it
requires that the PCTF3D method was applied prior to this script. Please 
refer to the script 'CTFlowHD_PCTF3D' if this is not the case.
************************************************

*******************  Step #2  ****************** Lines 1 to 9
As a first step, PCTF3D results must be loaded on your workspace WITH the 
observation matrix *X*. There are two possible ways to do this.

First, you can load a PCTF3D result saved file from the directory 
'./fcm_CTFlowHD/save/'. In order to save PCTF3D results, please refer to 
the script 'CTFlowHD_PCTF3D' and its instructions 'README_PCTF3D'. Do not 
forget that this script needs the observation matrix X to be saved. To save 
the variable X in a PCTF3D results saved file, you have to mention it as an 
option of the function 'savePCTF3D':
    >>> savePCTF3D("data")

Secondly, you can run the script 'CTFlowHD_plot_Marg' right after running 
the script 'CTFlowHD_PCTF3D'. To do so, you have to skip the first two 
sections of the script 'CTFlowHD_plot_Dendro' (Lines 1 to 9) as it clears 
your workspace.
************************************************

*******************  Step #4  ****************** Lines 10 to the 19
This section plots results of PCTF3D with a dendrogram visualization. It is
required to build a dendrogram to visualize afterwards the marginal of rank
one terms. You can close or keep this figure for the rest of the script.

After building the dendrogram between rank-1 components, the dendrogram is 
cut thanks to a threshold set by end users. Components are then clustered 
if together their normalized dendrogram distance is below the threshold.
You can add options to this function to change this visualization.

You can change the parameter 'strScreen' that permits to change the 
position and size of the figure. You can set this parameter with the 
following options:
    >>> ...,'strScreen','halfTop', ...) : half top of the screen,
    >>> ...,'strScreen','halfBot', ...) : half bottom of the screen,
    >>> ...,'strScreen','halfL', ...) : half left of the screen,
    >>> ...,'strScreen','halfR', ...) : half right of the screen,
    >>> ...,'strScreen','full', ...) : full screen (default).

You can also change the parameter 'strLanguage' that permits to change the 
language of the figure. You can set this parameter with the 
following options:
    >>> ...,'strLanguage','FR', ...) : for french,
    >>> ...,'strLanguage','FR', ...) : for english (default).

It is possible to choose how a distance between two rank-1 terms is 
computed. We propose in this framework 3 possible metrics:
    - 'max':
Each rank-one term is represented by a M vectors of size Ix1. First, the 
'max' metric finds the centroid of the bin which has maximum probability. 
This ouputs a Mx1 vector for each component of the decomposition. By 
stacking these vectors in a RxM matrix, it is possible to perform linkage 
by computing the Euclidian distance between vectors of maximum probability.
    - 'esp' (default):
Each rank-one term is represented by a M vectors of size Ix1. First, the 
'esp' metric finds the expected value of each 1D factor. This ouputs a Mx1 
vector for each component of the decomposition. By stacking these vectors 
in a RxM matrix, it is possible to perform linkage by computing the 
Euclidian distance between vectors of maximum probability.
    - 'corr':
For this metric, the distance between two rank-one terms is directly 
computed before the linkage. This distance is defined as the product of the
*M* correlations between the 1D vectors of probability. To enhance 
visualization, the distance computed for this metric are normalized. It 
does not affect the clustering method as we perform a hierarchical 
clustering afterwards.

It is possible to change the linkage method that permits to build the 
dendrogram. The possibilities for this parameter are the same as the 
'linkage' function developped by MATLAB. Hence, we give the section of the
documentation on that matter in the following:
'single'    --- nearest distance
'complete'  --- furthest distance
'average'   --- unweighted average distance (UPGMA) (also known as group 
                average)
'weighted'  --- weighted average distance (WPGMA)
'centroid'  --- unweighted center of mass distance (UPGMC)
'median'    --- weighted center of mass distance (WPGMC)
'ward'      --- inner squared distance (min variance algorithm)
In our flow cytometry framework, the 'complete' option is chosen as 
default.

Finally, it is possible to regroup rank-1 terms in the visualization by 
changing the option 'boolRegroup':
    >>> plot_Dendro(y,lambda,t,indVar,strLabel,threshDendro, ...
            'boolRegroup',1)
This option is 0 by default, meaning that rank 1 terms are plotted 
separately by default. However, rank-1 terms that have a distance below the
value *threshDendro* will have the same color. If 'boolRegroup' is set to 
1, clustered rank-one terms will be plotted as 1 component defined by the 
weighted sum of 1D factors.
************************************************

*******************  Step #3  ****************** Lines 20 to the end
This section plots the dendrogram clusters with a marginal visualization. 
The function 'plot_Marg' plots every 1D and 2D marginals possible from the 
set of variables in *indVar*. This results in a MxM grid of plots. On the 
m-th diagonal of the grid, the 1D marginals of each dendrogram cluster is 
plotted as well as the empirical estimation of the corresponding 1D 
marginal. For other plots of indices (j,k), a 2D scatter plot of the cells 
is represented. After this, dendrogram clusters are plotted with contours 
of the 2D marginal distribution. Contours are colored with the same color 
used in the dendrogram visualization (see Step 2).

It is possible to plot only a portion of those lower-order marginals. This 
can be performed by entering *indMarg* as an input of this function:
    >>> 'all' (default):
This option permits to plot every 1D and 2D marginals.
    >>> a variable of type cell:
If one inputs *indMarg* as a variable of type cell containing the indices 
of the requested marginals, the function will only plots those marginals.
First, an *indMarg* element can be a scalar corresponding to an index m
between 1 and *M*. This will results in the apparation of the m-th 1D 
marginal. Finally, an *indMarg* element can be a 1x2 vector of indices 
(j,k) between 1 and *M*. This will results in the plot of the 2D marginal 
for variables j and k.
Example:
    >>> indMarg = {1,2,3,[1,2],[3,1],[2,3]};
    >>> plot_Marg(y,lambda,compGroup,X,t,strLabel,indVar,'indMarg',indMarg)
This will result in a figure that plots the 1D marginals of indices 1 2 and 
3 as well as 2D marginals (1,2) (3,1) and (2,3).

It is possible to separate 1D marginal plots from 2D marginal plots. To do 
so, you can enter the option 'boolSeparation' with the following:
    >>> plot_Marg(y,lambda,compGroup,X,t,strLabel,indVar, ...
            'boolSeparation',1)
This default choice is separating the two types of plots.
    >>> plot_Marg(y,lambda,compGroup,X,t,strLabel,indVar, ...
            'boolSeparation',0)
This default choice is not separating the two types of plots. However, 1D 
marginals will always be plotted before 2D marginals. This means that it is
not possible to plot alternatively a 1D marginal then a 2D and after coming 
back to a 1D marginal.

Th dendrogram function 'plot_Dendro' returns an object *compGroup* that 
stores the dendrogram clusters. It is a 1xnCluster cell variable that 
contains the information of the *nCluster* groups of rank-one terms. By 
default, the function 'plot_Marg' is plotting every *compGroup* clusters 
but it is possible to plot a subset of clusters by entering the variable 
*indGroup* as an input of this function:
    >>> plot_Marg(y,lambda,compGroup,X,t,strLabel,indVar, ...
            'indGroup',[1,3])
This will only plot the first and third clusters.

You can change the parameter 'strScreen' that permits to change the 
position and size of the figure. You can set this parameter with the 
following options:
    >>> ...,'strScreen','halfTop', ...) : half top of the screen,
    >>> ...,'strScreen','halfBot', ...) : half bottom of the screen,
    >>> ...,'strScreen','halfL', ...) : half left of the screen,
    >>> ...,'strScreen','halfR', ...) : half right of the screen,
    >>> ...,'strScreen','full', ...) : full screen (default).

You can also change the parameter 'colMax' with the same method. You can
set this parameter with any strictly positive number. It will define the 
maximum number of plots along one column. This option will be ignored if 
all possible lower-order marginals are plotted.
    >>> indMarg = ...;
    >>> plot_Marg(y,lambda,compGroup,X,t,strLabel,indVar, ...
            'indMarg',indMarg,'colMax',2)

For large datasets, it is possible to plot 1 cell every *stepCloud* cells. 
You can change this optional parameter to 30 (default is 50) for example:
    >>> plot_Marg(y,lambda,compGroup,X,t,strLabel,indVar, ...
            'stepCloud',30)
************************************************

*******  Example given in the manuscript  ******
In the thesis manuscript of the author, figure 5.13 features an example of 
one 2D marginal plot for a flow cytometry controlled dataset.
************************************************
