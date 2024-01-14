************************************************
* README SPADE (other methods)
************************************************

*******************  Author  *******************
Philippe Flores (flores.philipe'at'gmail.com)
https://github.com/philippeflores/fcm_ctflowhd
************************************************

****************  Disclaimer #1  ***************
This document is a tutorial to use the SPADE method on raw FCM data. To use
SPADE on CTFlowHD results, please refer to the CTFlowHD folder.
************************************************

*****************  REFERENCES  *****************
QIU, Peng and SIMONDS, Erin F. and BENDALL, Sean C. and GIBBS, Kenneth Jr 
D. and BRUGGNER, Robert V. and LINDERMAN, Michael D. and SACHS, Karen and 
NOLAN, Garry P. and PLEVRITIS, Sylvia K., 
Extracting a cellular hierarchy from high-dimensional cytometry data with 
SPADE.
Nature biotechnology, 2011
************************************************

*******************  Step #1  ******************
The SPADE method takes pre-processed files as an input. Therefore, if you 
want to process a FCS file (i.e. 'file.fcs'), you must perform the 
pre-processing steps. To do this, please refer to the folder 
'./fcm_ctflowhd/pre_processingR/'.
************************************************

*******************  Step #2  ******************
To analyze a pre-processed FCS dataset (i.e. 'file_CCT.fcs'), this file 
must be inside the folder './fcm_ctflowhd/other_methods/SPADE/data/'. This 
folder is not tracked in the git repository for memory concerns. Because of 
that, it may not exist for your first use of SPADE. In this case, you may 
create the folder './fcm_ctflowhd/other_methods/SPADE/data/'.

If the folder './fcm_ctflowhd/other_methods/SPADE/data' already exists, put
the file 'file_CCT.fcs' in the folder 
'./fcm_ctflowhd/other_methods/SPADE/data'.
************************************************

*******************  Step #3  ******************
Open MATLAB. Change your working directory to 
'./fcm_ctflowhd/other_methods/SPADE'. Type in the command window the 
following command:
	>>> SPADE
This should open the SPADE window.
************************************************

*******************  Step #4  ******************
Select the working directory. On top of the SPADE window, select the 
'Browse' button to select the working directory. A window opens. Navigate 
in order to chose as working directory the folder 
'./fcm_ctflowhd/other_methods/SPADE/data'.

At this point, if a FCS file is present in 
'./fcm_ctflowhd/other_methods/SPADE/data', a window that mention this file 
should appear. Close this window. You can open this window again by 
clicking on the 'View File List' button of the SPADE window.
************************************************

*******************  Step #5  ******************
Change SPADE parameters. To change parameters of SPADE, click on the button 
'View / update algorithm parameters'.

SPADE features a memory buffer for parameters. From one use to another, it 
is possible that parameters are already entered. If this the case, verify 
that those settings are good.

If no parameters are entered, parameters should be chosen:
	1) 'markers used to build SPADE tree':
This parameter permits to select which FCM variables are chosen to build 
the SPADE tree. The list on the left is for variables not used while the 
right list shows variables used in the SPADE tree build. Use the arrows to 
select or remove a variable.
	2) 'Compensation option':
Because the file was pre-processed before SPADE application, select 
'ignore compensation'.
	3) 'Arcsinh Transformation Option':
Because the file was pre-processed before SPADE application, select 
'No transformation'.
	4) 'Local density calculation parameters':
Default parameters are fine. You can tweak them if necessary. Please check 
the SPADE article on that matter.
	5) 'Files used to build SPADE tree':
Same as for variables. SPADE offers the possibility to analyse several 
datasets all at once. Select the desired pool of files. On the bottom of 
this sub-window, the maximum allowable cells in pooled downsampled data 
can be chosen.
	6) 'Outlier and Target Densities (OD and TD)':
Default parameters are fine. You can tweak them if necessary. Please check 
the SPADE article on that matter.
	7) 'Clustering parameter':
It is possible to chose between an agglomerative algorithm and K-means for 
the clustering. Those clusters will be the nodes of the visualisation tree. 
Their number can be chosen too.

When all parameters are chosen, you can click on the 'Close' button at the 
bottom right of this window.
************************************************

*******************  Step #6  ******************
You are now back to the SPADE window. You must now 'Run SPADE Analysis'. To 
do so, click on the button 'One click for all'. SPADE is running and output 
messages are shown in the command window of MATLAB.
************************************************

*******************  Step #7  ******************
Once SPADE has run, visualisations are available by clicking on the 
'View resulting SPADE tree' button at the bottom of the SPADE window.

Some of the window features are presented in the following.
	1) SPADE tree
The SPADE tree is visualized at the top right of the visualisation window. 
Nodes can be selected for further analysis.
	2) Bivariate plots
At the bottom left of the window, bivariate scatter plots and a contour 
plots of cells are shown. You can chose the variables for the plots. When 
nodes are selected and added to annotation, you can update plots to show 
those specific cells on those bivariate plots.
	3) Annotation panel
At the top right of the window, it is possible to manage annotations and 
selected nodes. 
	4) Color options
At the bottom right of the window, it is possible to manage the 
colorization of the SPADE tree. To do this, check the box 'overlay 
information by coloring nodes'. You can then select the FCM variable so 
that each node is colorised regarding the mean value for this variable.	
    5) Export results
It is possible to export results with the bottom right panel. For example, 
is is possible to export colorised figures, export annotation cell 
proportion, clustering results, etc.
************************************************

*******  Example given in the manuscript  ******
In the manuscript, figure 5.18 features an example of SPADE applied to a
4-fluorescence controlled dataset. Feel free to ask me by e-mail for this 
dataset.
************************************************
