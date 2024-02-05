clear
close all

addpath("./include")

%% Load PCTF3D results

load("save/file_CCT_I25R47T4.mat");

%% Plot dendrogram clusters with marginal visualization
close all

K = 3;

indMarg = {1,2,3,[1,2],[2,3],[3,1]};

plot_NBMkMeans(y,lambda,K,indVar,t,strLabel,'indMarg',indMarg,'colMax',3)