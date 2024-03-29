%% Set up the workspace
clear 
close all
clc

addpath("./include")

%% Load filename

filename = "file_CCT";

[X,fcsHdr,strFile] = loadFilename(filename);
strVariables = {fcsHdr.par.marker}';

%% Selection of variables

indVar = [];

[indVar,M,strLabel] = loadIndVar(X,indVar,fcsHdr);

strLabel = changeLabel(strLabel);

%% Definition of variable spaces

strMethod = 'fluo';
I = 25;

t = supportDist(X,indVar,I,'strMethod',strMethod);

%% Choice of subset of triplets and computation of 3D marginals

strStrategy = '1/2';
% T = 20;

[calTfull,Tfull] = createTriplets(M,indVar,'full');
tic, dataMargFull = computeMarginals(X,indVar,t,I,calTfull); timeMarg = toc;

if strcmpi(strStrategy,'full') && exist('calTfull','var')
    calT = calTfull;
else
    [calT,T] = createTriplets(M,indVar,strStrategy);
    tic, dataMarg = computeMarginals(X,indVar,t,I,calT); timeMarg = toc;
end

%% PCTF3D optimization algorithm

R = 20;

threshLambda = 1/(20*R);
T1 = 1000;
T2 = 20;

tic, [y,lambda,cost,error3,error1,flag,R] = PCTF3D(dataMarg,R, ...
    'threshLambda',threshLambda,'dataMargFull',dataMargFull,'computeError',0); timePCTF3D = toc; timePCTF3D = timePCTF3D+timeMarg;
fprintf("\nComputation time PCTF3D : %fs\n",timePCTF3D)

%% Saving results

savePCTF3D("data")
