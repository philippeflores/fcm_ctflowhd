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

strStrategy = 'full';
% T = 20;

[calT,T] = createTriplets(M,indVar,strStrategy);
dataMarg = computeMarginals(X,indVar,t,I,calT);

%% PCTF3D

R = 50;

threshLambda = 1/(20*R);
T1 = 2000;
T2 = 50;

[y,lambda,erreur3,erreur1,flag,R] = PCTF3D(dataMarg,R,'threshLambda',threshLambda);

%% Saving results

savePCTF3D
