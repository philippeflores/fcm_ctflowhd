************************************************
* README CTFlowHD_PCTF3D
************************************************

*******************  Author  *******************
Philippe Flores (flores.philipe'at'gmail.com)
https://github.com/philippeflores/fcm_ctflowhd
************************************************

*****************  REFERENCES  *****************
PCTF3D:
Philippe Flores, Guillaume Harlé, Anne-Béatrice Notarantonio, Konstantin 
Usevich, Maud d'Aveni, Stéphanie Grandemange, Marie-Thérèse Rubio, and 
David Brie.
Coupled tensor factorization for flow cytometry data analysis.
2022 IEEE 32nd International Workshop on Machine Learning for Signal 
Processing (MLSP), IEEE, 2022.

N-way toolbox:
Clauss Andersson, and Rasmus Bro.
The N-way toolbox for MATLAB.
Chemometrics and intelligent laboratory systems, Elsevier, 2000.
************************************************

****************  Disclaimer #1  ***************
The script 'CTFlowHD_PCTF3D' contains a group of functions that require the
installation of the package 'N-way' that is mentioned in the references.
************************************************

****************  Disclaimer #2  ***************
This script does not contain any visualization features. Clustering and 
visualization steps are included in other scripts.
************************************************

*******************  Step #1  ******************
The 'CTFlowHD_PCTF3D' script takes pre-processed files as an input. 
Therefore, if you want to process a FCS file (i.e. 'file.fcs'), you must 
perform the pre-processing steps. To do this, please refer to the folder 
'./fcm_ctflowhd/pre_processingR/'.
************************************************

*******************  Step #2  ******************
To analyze a pre-processed FCS dataset (i.e. 'file_CCT.fcs'), this file 
must be inside the folder './fcm_ctflowhd/data/'. This folder is not 
tracked in the git repository for memory concerns. Because of that, it may 
not exist for your first use of CTFlowHD_PCTF3D. In this case, you may 
create the folder './fcm_ctflowhd/data/'.

If the folder './fcm_ctflowhd/data' already exists, put the file 
'file_CCT.fcs' in the folder './fcm_ctflowhd/data'.
************************************************

*******************  Step #3  ******************
Open MATLAB. Change your working directory to 
'./fcm_ctflowhd/'. Open the script file named 'CTFlowHD_PCTF3D.m'.
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
If you select a file that is not in the './fcm_ctflowhd/data'
directory, the function will ask you to type a valid filename (if this is 
possible).
************************************************

*******************  Step #6  ****************** >>> Lines 15 to 22
Select the variables to run PCTF3D. At line 16, this is possible to 
manually type the variables. To help you in your choice, the variable names
are stored in the variable 'strVariables'.

If you choose:
    >>> indVar = [];
only the fluorescence variables are going to be selected.

At line 21, there is the possibility to change variable labels. To do so, 
uncomment the line 21 and run the section again. Follow the instructions to
rename each label.
************************************************

*******************  Step #7  ****************** >>> Lines 23 to 29
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

*******************  Step #8  ****************** >>> Lines 30 to 37
In this section, the subset of marginals to be considered in the coupling
of PCTF3D is set. To do so, you can choose both the strategy and the number
of triplets. After the subset of triplets 'calT' is defined, the resulted
3D marginals are computed.

To change the 'strStrategy', you have the following options:
    - 'full' :
All triplets are considered. For this option, the number of triplets is 
nchoosek(M,3) where M is the number of variables in 'indVar'.
    >>> strStrategy = 'full';
    >>> [calT,T] = createTriplets(M,indVar,strStrategy);
    - '+1' :
M triplets are considered. The first triplet is [1 2 3], the second is
[2 3 4], then [3 4 5], and so on until [M-2 M-1 M], [M-1 M 1] and 
[M 1 2]. Between two triplets, the operation +1 is applied.
    >>> strStrategy = '+1';
    >>> [calT,T] = createTriplets(M,indVar,strStrategy);
    - '+2' :
\lfloor M/2 \rfloor triplets are considered. The first triplet is [1 2 3], 
the second is [3 4 5], then [5 6 7], and so on until [M-2 M-1 M] if M is 
odd and [M-3 M-2 M-1], [M-2 M-1 M] if M is even. Between two triplets, the 
operation +2 is applied.
    >>> strStrategy = '+2';
    >>> [calT,T] = createTriplets(M,indVar,strStrategy);
    - 'rng' :
This strategy is random. It requires that the user also enters a number of 
triplets as an input of this function (see Example at the end). The number 
of triplets must be lower or equal to nchoosek(M,3) and greater or equal to
\lfloor M/2 \rfloor.
    >>> strStrategy = 'rng';
    >>> T = 20;
    >>> [calT,T] = createTriplets(M,indVar,strStrategy,'Tin',T);
    - 'bal' :
This strategy is balanced. This means that the resulting coupling is chosen
such that each variable in 'indVar' appears the same number of times. This 
number is defines by round(3*T/M) where T is the number of triplets and M
is the number of variables in 'indVar'. This strategy requires that the 
user also enters a number of triplets as an input of this function (see 
Example at the end). The number of triplets must be lower or equal to 
nchoosek(M,3) and greater or equal to \lfloor M/2 \rfloor.
    >>> strStrategy = 'bal';
    >>> T = 20;
    >>> [calT,T] = createTriplets(M,indVar,strStrategy,'Tin',T);
    - 'man' :
If you want to enter manually the subset of triplets, it is possible by 
using this option. Therefore, you will need to enter the parameter 'calTin'
as an input of this function. 'calTin' must be a 1xT cell variable. The 
t-th element 'calTin{t}' must be a 1x3 vector containing the triplet 
variables.
    >>> strStrategy = 'man';
    >>> calT = {[1,2,3],[1,2,4],[1,2,5],[1,3,4],[1,3,5],[1,4,5]};
    >>> [calT,T] = createTriplets(M,indVar,strStrategy,'calTin',calT);
    - 'A/B' :
This option will select randomly triplets with a proportion of A triplets 
out of B triplets. The number of triplets will be equal to 
round(nchoosek(M,3)*A/B).
    >>> strStrategy = '3/4';
    >>> [calT,T] = createTriplets(M,indVar,strStrategy);
For this example, if there are M=6 variables, T=15 triplets will be chosen
randomly. 
    - 'A/Bbal' :
Idem for a balanced selection of round(nchoosek(M,3)*A/B) triplets.
    >>> strStrategy = '1/4bal';
    >>> [calT,T] = createTriplets(M,indVar,strStrategy);
For this example, if there are M=6 variables, T=5 triplets will be chosen
in a balanced way.
    - '0.F' :
This option will select randomly triplets with a proportion of 0.F triplets
randomly hence round(nchoosek(M,3)*0.F).
    >>> strStrategy = '0.4';
    >>> [calT,T] = createTriplets(M,indVar,strStrategy);
For this example, if there are M=6 variables, T=8 triplets will be chosen
randomly. 
    - '0.Fbal' :
Idem for a balanced selection of round(nchoosek(M,3)*0.F) triplets.
    >>> strStrategy = '0.6bal';
    >>> [calT,T] = createTriplets(M,indVar,strStrategy);
For this example, if there are M=6 variables, T=12 triplets will be chosen
in a balanced way.
************************************************

*******************  Step #9  ****************** >>> Lines 38 to 48
This section builds the NBM of X thanks to the PCTF3D algorithm. The only 
mandatory argmument of this algorithm is the number of NBM components 'R'.
Moreover, it is possible to change the following parameters:
    - 'y0' (default --> random) :
It is possible to change the initialization of PCTF3D. 'y0' is the variable
that contains the M factors matrices where M is the number of variables in
'indVar'. Therefore, 'y0' must be a 1xM cell and each y0{m} must be a IxR
matrix.
    >>> y0 = ...;
    >>> [y,lambda] = PCTF3D(dataMarg,R,'y0',y0);
    - 'lambda0' (default --> random) :
It is possible to change the initialization of PCTF3D. 'lambda0' is the 
variable that contains the weight of each NBM component. Therefore, 
'lambda0' must be a Rx1 vector.
    >>> lambda0 = ...;
    >>> [y,lambda] = PCTF3D(dataMarg,R,'lambda0',lambda0);
    - 'T1' (default --> 1000) :
T1 permits to change the maximal number of outer iterations of ADMM.
    >>> T1 = 2000;
    >>> [y,lambda] = PCTF3D(dataMarg,R,'T1',T1);
    - 'T2' (default --> 20) :
T2 permits to change the maximal number of inner iterations of ADMM.
    >>> T2 = 50;
    >>> [y,lambda] = PCTF3D(dataMarg,R,'T2',T2);
    - 'eps':
eps permits to change the stopping criterion of ADMM.
    >>> eps = 10^(-7);
    >>> [y,lambda] = PCTF3D(dataMarg,R,'eps',eps);
    - 'tolNewY' (default --> 10^(-7)) :
tolNewY permits to change the stopping criterion of PCTF3D. if the 
difference between two iterations of PCTF3D is lower than tolNewY, PCTF3D
stops.
    >>> tolNewY = 10^(-6);
    >>> [y,lambda] = PCTF3D(dataMarg,R,'tolNewY',tolNewY);
    - 'boolPlot' (default --> 1) :
If 'boolPlot' is set to 0, no information will be plotted to the user.
    >>> [y,lambda] = PCTF3D(dataMarg,R,'boolPlot',0);
    - 'threshLambda' (default --> -1) :
It may be possible that resulting weights are exactly equal to 0 or very
small. To delete them as they do not contribute to the sum of NBM 
component, it is possible to change the parameter 'threshLambda'.
    >>> threshLambda = 1/(20*R);
    >>> [y,lambda] = PCTF3D(dataMarg,R,'threshLambda',threshLambda);
    - 'computeError' (default --> 0) :
This option permits to choose if or not the estimation errors are computed.
    >>> computeError = 1;
    >>> [y,lambda] = PCTF3D(dataMarg,R,'computeError',computeError);
************************************************

*******************  Step #10  ****************** >>> Lines 49 to the end
This section permits to save PCTF3D results. Saved files will appear in 
the folder './fcm_ctflowhd/save/'. This folder is not tracked in the git 
repository for memory concerns. Because of that, it may not exist for your 
first use of CTFlowHD_PCTF3D. In this case, you may create the folder 
'./fcm_ctflowhd/save/'. After this, you can run the function 'savePCTF3D'.

Without arguments, this function will save all the workspace of this script
except 'dataMarg' and 'X' which are too big for storage. The name of the 
saved file will be chosen as follows (with * between name of variables in 
the workspace): *strFile*_I*I*R*R*T*T*.mat. For the example of 'file_CCT' 
with 'I'=25, 'R'=10 and 'T'=4 the filename will be 'file_CCT_I25R10T4.mat'.

If you want to enter manually the output filename, you can enter it as a 
parameter of this function (with or without the prefix './save/'):
    >>> savePCTF3D('saveFile.mat')

If you want to save either marginal data or the data matrix X inside the 
saved file, it is possible by entering the following input in the function.
To save the observation matrix 'X':
    >>> savePCTF3D('data') 
To save the marginals 'dataMarg':
    >>> savePCTF3D('marg')
These inputs can be combined together. They can also be combined with the 
choice of the name of the saved file.
************************************************
