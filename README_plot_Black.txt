************************************************
* README CTFlowHD_plot_Black
************************************************

*******************  Author  *******************
Philippe Flores (flores.philipe@gmail.com)
https://github.com/philippeflores/fcm_ctflowhd
************************************************

*******************  Step #1  ******************
The script 'CTFlowHD_plot_Black' is a visualization script. Therefore, it
requires that the PCTF3D method was applied prior to this script. Please 
refer to the script 'CTFlowHD_PCTF3D' if this is not the case.
************************************************

*******************  Step #2  ****************** Lines 1 to 9
As a first step, PCTF3D results must be loaded on your workspace. There are
two possible ways to do this.

First, you can load a PCTF3D result saved file from the directory 
'./fcm_CTFlowHD/save/'. In order to save PCTF3D results, please refer to 
the script 'CTFlowHD_PCTF3D' and its instructions 'README_PCTF3D'.

Secondly, you can run the script 'CTFlowHD_plot_Black' right after running 
the script 'CTFlowHD_PCTF3D'. To do so, you have to skip the first two 
sections of the script 'CTFlowHD_plot_Black' (Lines 1 to 9) as it clears 
your workspace.
************************************************

*******************  Step #3  ****************** Lines 10 to the end
This section plots results of PCTF3D without any clustering method. Because
of this, all components are plotted in black. 

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
************************************************

*******  Example given in the manuscript  ******
In the manuscript, figure 5.10 features an example of a plot without 
clustering method for a controlled dataset.
************************************************
