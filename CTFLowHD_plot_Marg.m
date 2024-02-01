clear
close all

addpath("./include")

%% Load PCTF3D results

load("save/file_CCT_I25R45T4.mat");

%% Create dendrogram visualization
close all

threshDendro = 0.5;
strDist = 'max';
strLink = 'single';

[y,lambda,compGroup,matDist,matLinkage,dendro] = plot_Dendro(y,lambda,t, ...
    indVar,strLabel,threshDendro,'strDist',strDist,'strLink',strLink,'strScreen','halfL');

%% Plot dendrogram clusters with marginal visualization

stepCloud = 10;

plot_Marg(y,lambda,compGroup,X,t,strLabel,indVar, ...
    'stepCloud',stepCloud,'strScreen','halfR','indMarg',{1,2,3,[2,3],[3,4],[1,3]})