%% Set up the workspace
clear 
close all
clc

addpath("../include")

%% Load filename

filename = "file_CCT";

[X,fcsHdr,strFile] = loadFilename(filename);
strVariables = {fcsHdr.par.marker}';

%% Selection of variables to build the K-means map

indVar = [];

[indVar,M,strLabel] = loadIndVar(X,indVar,fcsHdr);

% strLabel = changeLabel(strLabel);

%% Definition of variable spaces

strMethod = 'fluo';

t = supportDist(X,indVar,100,'strMethod',strMethod);

%% Building the K-means model

K = 3;

tic, [idx,C] = kmeans(X(:,indVar),K); time_Kmeans = toc;
fprintf("Computation time for K-means : %fs\n",time_Kmeans)     

%% Plot the Kmeans clusters

stepCloud = 10;

close all
plot_Kmeans(X,indVar,t,idx,C,strLabel,'stepCloud',stepCloud)
