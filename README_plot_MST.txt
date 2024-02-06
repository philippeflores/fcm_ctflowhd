************************************************
* README CTFlowHD_plot_MST
************************************************

*******************  Author  *******************
Philippe Flores (flores.philipe@gmail.com)
https://github.com/philippeflores/fcm_ctflowhd
************************************************

*******************  Step #1  ******************
The script 'CTFlowHD_plot_MST' is a visualization script. Therefore, it
requires that the PCTF3D method was applied prior to this script. Please 
refer to the script 'CTFlowHD_PCTF3D' if this is not the case.
************************************************

*******************  Step #2  ****************** Lines 1 to 9
As a first step, PCTF3D results must be loaded on your workspace. There are
two possible ways to do this.

First, you can load a PCTF3D result saved file from the directory 
'./fcm_CTFlowHD/save/'. In order to save PCTF3D results, please refer to 
the script 'CTFlowHD_PCTF3D' and its instructions 'README_PCTF3D'.

Secondly, you can run the script 'CTFlowHD_plot_MST' right after running 
the script 'CTFlowHD_PCTF3D'. To do so, you have to skip the first two 
sections of the script 'CTFlowHD_plot_MST' (Lines 1 to 9) as it clears 
your workspace.
************************************************

*******************  Step #3  ****************** Lines 10 to the 19
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

*******************  Step #4  ****************** Lines 20 to the end
This section plots the dendrogram clusters with a MST visualization. 
The function 'plot_MST' builds a Minimum Spanning Tree (MST) to visualize
how components are distributed. First, centroids for each component are 
computed which permit to create a centroid matrix of size RxM. This matrix 
enables the building of the MST. First, the MST is plotted with the 
colorization of the dendrogram previously plotted. Finally, for each marker
a colorized MST is plotted with mean value of fluorescence. 

To perform MST on rank-one terms, the centroid of each component must 
be computed. To do so, it is possible to change the option 'strCentroid' 
of the function 'plot_MST'. The two options for this parameter are 
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

You can also change the parameter 'colMax' with the same method. You can
set this parameter with any positive number. It will define the maximum 
number of plots along one column.
    >>> Y = plot_MST(y,lambda,compGroup,t,strLabel,indVar,'colMax',0);
************************************************

*******  Example given in the manuscript  ******
In the thesis manuscript of the author, figure 5.17 features an example of 
a MST visualization for a flow cytometry controlled dataset. This method of
visualization is used in the SPADE method used for flow cytometry data 
analysis.
************************************************
