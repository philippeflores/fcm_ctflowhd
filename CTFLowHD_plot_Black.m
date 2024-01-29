clear
close all

addpath("./include")

%% Load PCTF3D results

load("save/file_CCT_I25R10T4.mat");

%% Plot PCTF3D results 

close all

[y,lambda] = plot_NBMblack(y,lambda,t,indVar,strLabel);