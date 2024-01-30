clear
close all

addpath("./include")

%% Load PCTF3D results

load("save/file_CCT_I25R47T4.mat");

%% Plot PCTF3D results 
close all

threshDendro = 0.5;
strDist = 'max';
strLink = 'single';

[y,lambda,compGroup,matDist,matLinkage,dendro] = plot_Dendro(y,lambda,t, ...
    indVar,strLabel,threshDendro,'strDist',strDist,'strLink',strLink);