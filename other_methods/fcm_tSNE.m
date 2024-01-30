%% Set up the workspace
clear 
close all
clc

addpath("../include")

%% Load filename

filename = "file_CCT";

[X,fcsHdr,strFile] = loadFilename(filename);
strVariables = {fcsHdr.par.marker}';

%% Selection of variables to build the t-SNE map

indVar = [];

[indVar,M,strLabel] = loadIndVar(X,indVar,fcsHdr);

% strLabel = changeLabel(strLabel);

%% Building the t-SNE map

perplexity = 200;

tic, Y = tsne(X(:,indVar),'Perplexity',perplexity); time_tSNE = toc;
fprintf("Computation time for tSNE : %fs\n",time_tSNE)     

%% Plot the t-SNE maps colorized with marker fluorescence

stepCloud = 1;

close all
plot_tSNE(Y,indVar,X,strLabel,'stepCloud',stepCloud)
