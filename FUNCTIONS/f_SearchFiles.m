function [c_Files] = f_SearchFiles (str_Folder,str_Ext)
%% f_SearchFiles
% Function by: Juan Nieto
% function to search for all the files in the given folder and subfolder
% witht the given extention, it doesnt considers lower and upper case.
% INPUTS: 
%   str_Folder= Principal folder containing all the subfolders and files to
%   search on 
%   str_Ext = Extention to validate from the different files in each folder
%   and subfolder. 
% Outputs: 
%   c_Files = cell containing all the files that share the givven extention

%% initialize variables 
c_Files = {};
%% Find all the files with the extention
% Get all files and directories
st_Sub = dir(str_Folder);
% Remove . & .. folders
st_Sub = st_Sub (3:length(st_Sub));
% Separate folders and files
st_Files = st_Sub(~[st_Sub.isdir]); % all files in the folder
st_Folders =st_Sub([st_Sub.isdir]); % all subfolders in the folder
%% add all the files ending in the given extention to the cell with all the names
c_Files = st_Sub(endsWith({st_Files.name},str_Ext));
%% find all files in the subfolders
for idx = 1:length(st_Folders)
    str_currFolder = fullfile(st_Folders(idx).folder,st_Folders(idx).name);
    %look for all the cells in the subfolders.
    [c_Files] = [c_Files;f_SearchFiles(str_currFolder,str_Ext)];
end 

