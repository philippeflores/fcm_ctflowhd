%% Set up the workspace
clear 
close all
clc

addpath("../include")

%% Load filename

filename = "file_CCT";

[X,fcsHdr,strFile] = loadFilename(filename);
strVariables = {fcsHdr.par.marker}';

%% Selection of variables

indVar = [];

[indVar,M,strLabel] = loadIndVar(X,indVar,fcsHdr);

% strLabel = changeLabel(strLabel);

%% Definition of variable spaces

strMethod = 'fluo';

t = supportDist(X,indVar,100,'strMethod',strMethod);

%% Building the Gaussian Mixture Model

R = 4;		   % Nombre de composantes (Plus il y en a plus le temps de calcul est long)

tic, [g,y,lambda] = buildGMM(X,indVar,R,t); time_GMM = toc;
fprintf("Computation time for Gaussian Mixture Model : %fs\n",time_GMM)     

%% Plot the GMM with single cell bivariate plots

stepCloud = 1000;

close all
plot_GMM(g,X,indVar,t,strLabel,'stepCloud',stepCloud)


