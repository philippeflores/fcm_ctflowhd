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

%% Plot dendrogram clusters with Minimum Spanning Tree (MST) visualization

plot_MST(y,lambda,strLabel,compGroup,t,indVar,'strCentroid','max','colMax',0)