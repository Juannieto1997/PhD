%% p_MeanCompleteMatrix
% Code to generate the figure with the complete matrix view of the mean
% activity according to the trigger data. Allows to select one single file
% and saving location.
%% Preparation of the workspace 
clc, close all; clear; 
% Add functio folder to the path
addpath("C:\Users\Juan\Documents\CODES\FUNCTIONS")
%path fo save the figures
str_saveFolder = uigetdir("D:\Figures\");
% Define extention to look for 
str_ext = 'CutSegments.mat';
% Find all files
[str_file,str_path]=uigetfile("D:\Data\");
%% loading data 
fprintf('reading file %s \n',str_file)
load (fullfile(str_path, str_file)); 
fprintf('File sucessfuly loaded \n')

f_MeanCompleteMatriz (c_cut,st_header.ChannOrder,str_saveFolder,str_file(1:10),[-50 100])

