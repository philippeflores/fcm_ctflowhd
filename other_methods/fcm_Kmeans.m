%% Set up the workspace
clear 
close all
clc

addpath("../include")

%% Load filename

filename = "file_CCT";

[X,fcsHdr,strFile] = loadFilename(filename);
strVariables = {fcsHdr.par.marker}';

%% Selection of variables to plot the K-means clusters

indVar = [];

[indVar,M,strLabel] = loadIndVar(X,indVar,fcsHdr);

% strLabel = changeLabel(strLabel);

%% Definition of variable spaces

strMethod = 'fluo';

t = supportDist(X,indVar,100,'strMethod',strMethod);

%% Define the K-means clusters

K = 3;

tic, [idx,C] = kmeans(X(:,indVar),K); time_Kmeans = toc;
fprintf("Computation time for K-means : %fs\n",time_Kmeans)     

%% Plot the K-means clusters

stepCloud = 10;

close all
plot_Kmeans(X,indVar,t,idx,C,strLabel,'stepCloud',stepCloud)
