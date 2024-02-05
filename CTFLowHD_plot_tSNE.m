clear
close all

addpath("./include")

%% Load PCTF3D results

load("save/file_CCT_I25R47T4.mat");

%% Create dendrogram visualization
close all

threshDendro = 0.5;
strDist = 'max';
strLink = 'single';

[y,lambda,compGroup,matDist,matLinkage,dendro] = plot_Dendro(y,lambda,t, ...
    indVar,strLabel,threshDendro,'strDist',strDist,'strLink',strLink,'strScreen','halfL');

%% Plot dendrogram clusters with t-SNE visualization

Y = plot_NBMtSNE(y,lambda,compGroup,t,strLabel,indVar,'colMax',0);