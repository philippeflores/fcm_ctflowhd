clear
close all

addpath("./include")

%% Load PCTF3D results

load("save/file_CCT_I25R45T4.mat");

%% NBM raw visualization without clustering
close all

[y,lambda] = plot_NBMblack(y,lambda,t,indVar,strLabel);

%% NBM visualization with dendrogram clustering 

threshDendro = 0.5;
strDist = 'max';
strLink = 'single';

[y,lambda,compGroup,matDist,matLinkage,dendro] = plot_Dendro(y,lambda,t, ...
    indVar,strLabel,threshDendro,'strDist',strDist,'strLink',strLink);

%% Dendrogram clustering with marginal visualizations 

stepCloud = 10;

plot_Marg(y,lambda,compGroup,X,t,strLabel,indVar, ...
    'stepCloud',stepCloud,'strScreen','halfR','indMarg',{1,2,3,[2,3],[3,4],[1,3]})

%% Dendrogram clustering with t-SNE maps

Y = plot_NBMtSNE(y,lambda,compGroup,t,strLabel,indVar,'colMax',0);

%% Dendrogram clustering with Minimum Spanning Tree (MST) visualizations

plot_MST(y,lambda,strLabel,compGroup,t,indVar,'strCentroid','max','colMax',0)

%% K-means clusters

K = 3;
indMarg = {1,2,3,[1,2],[2,3],[3,1]};

plot_NBMkMeans(y,lambda,K,indVar,t,strLabel,'indMarg',indMarg,'colMax',3)
