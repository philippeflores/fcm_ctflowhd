************************************************
* README CTFlowHD_Visualizations
************************************************

*******************  Author  *******************
Philippe Flores (flores.philipe@gmail.com)
https://github.com/philippeflores/fcm_ctflowhd
************************************************

*******************  Step #1  ******************
The script 'CTFlowHD_Visualizations' is a visualization script. Therefore, 
it requires that the PCTF3D method was applied prior to this script. Please 
refer to the script 'CTFlowHD_PCTF3D' if this is not the case.
************************************************

*******************  Step #2  ****************** Lines 1 to 9
As a first step, PCTF3D results must be loaded on your workspace WITH the 
observation matrix *X*. There are two possible ways to do this.

First, you can load a PCTF3D result saved file from the directory 
'./fcm_CTFlowHD/save/'. In order to save PCTF3D results, please refer to 
the script 'CTFlowHD_PCTF3D' and its instructions 'README_PCTF3D'. Do not 
forget that some visualization methods used in this script needs the 
observation matrix *X* to be saved. To save the variable *X* in a PCTF3D 
results saved file, you have to mention it as an option of the function 
'savePCTF3D':
    >>> savePCTF3D("data")

Secondly, you can run the script 'CTFlowHD_Visualizations' right after 
running the script 'CTFlowHD_PCTF3D'. To do so, you have to skip the first 
two sections of the script 'CTFlowHD_Visualizations' (Lines 1 to 9) as it 
clears your workspace.
************************************************

*******************  Step #3  ****************** Lines 10 to the end
This script is a collection of all visualization tools present in this 
package. Each section of code after line 10 permits to plot a different 
visualization. Each visualization tool is thoroughly explained in its 
correponding README file attached to this git repository. In summary, here 
is the list of all the visualization tools contained in this script:
    - NBM raw visualization without clustering *README_plot_Black*
This visualization plots raw factors of the NBM without any clustering 
methods.
    - NBM visualization with dendrogram clustering *README_plot_Dendro*
This visualization orders and colors rank-one terms that are similar thanks
to a dendrogram tree build based on hierarchical clustering.
    - Dendrogram clustering with marginal visualizations *README_plot_Marg*
Based on a dendrogram clustering, lower-order marginals for each dendrogram
cluster are plotted.
    - Dendrogram clustering with t-SNE maps *README_plot_tSNE*
t-SNE maps of rank-one terms is built and colorized with marker
fluorescence.
    - Dendrogram clustering with MST visualizations *README_plot_MST*
A Minimum Spanning Tree (MST) of rank-one terms is built and colorized with 
marker fluorescence.
    - K-means clusters *README_plot_Kmeans*
The Kmeans clustering method is applied to rank-one terms and clustering
results are plotted for lower-order marginals.
************************************************
