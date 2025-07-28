%% Initialize workspace 
%clear workspace
clc;close all; clear; 
%add functions folder to workspace 
addpath("C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS")
%% Initialize Variables
% Variables for files
% Select folder using explorer
[str_path]=uigetdir("D:");
% get all folders ending with "RawData.mat"
% Extetion to look for
str_ext = 'eFP.mat';
[st_Files] = f_SearchFiles (str_path,str_ext);
% Excel file with the experiment information 
str_exc = 'C:\Users\Juan\Desktop\FUS Success.xlsx'; 
t_data = readtable (str_exc);
s_group = 12;str_unitsG = '(s)'; % How will the data be groupled (same parameter for all 4 = driving voltage, 12 = duration)
s_treat = 4; str_nameT = 'Driving Voltage'; str_unitsT = '(uV)';% what is chanign in each case 
% Variables for plot
str_Var = 'Amplitude'; str_unitsV = '(mV)'; 
str_chann = 'C8';
v_DataBox = []; 
v_TreaBox = [];
%% Extract Data from excel 
c_recs = table2cell(t_data(:,"Rec"));
v_group = table2array(t_data(:,s_group));
v_treat =table2array(t_data(:,s_treat));

%% Remove not necesary files 
%Remove files not in the excel
c_names = {st_Files.name};
b_keep = contains(c_names,c_recs);
st_Files(~b_keep) = [];
% Femove files that are not in the folder
c_names = split({st_Files.name},str_ext);
c_noEx  = {c_names{1,:,1}};
c_names ={st_Files.name},str_ext;
b_keep = contains(c_recs,c_noEx);
c_recs(~b_keep)= [];
v_treat(~b_keep)= [];
v_group(~b_keep)= [];
 
v_UnGroup = unique(v_group); % Unique treatments 


clearvars t_data c_noEx %clear table for ease of writing the fing code :)

for idxG = 1:length(v_UnGroup)
    s_currGroup = v_UnGroup(idxG); 
    b_treat = v_group==s_currGroup;
    if sum(b_treat)<=1
        continue 
    end 
    c_CurrRec = c_recs(b_treat);
    v_DataBox = []; 
    v_TreaBox = [];
    c_DataBox = {};
    c_TreatBox ={};
    for idxR = 1:length(c_CurrRec)
        str_CurrRec = c_CurrRec{idxR};
        s_CurrTreat = v_treat(find(contains(c_recs,str_CurrRec)));
        s_FileIdx = find(contains(c_names,str_CurrRec));
        str_Filename = fullfile(st_Files(s_FileIdx).folder,st_Files(s_FileIdx).name);
        load(str_Filename);
        s_ChanIdx = find(contains({st_features.ChannName},str_chann));
        v_DataBox = [v_DataBox st_features(s_ChanIdx).(str_Var)];
        v_TreaBox(length(v_TreaBox)+1:length(v_DataBox)) = s_CurrTreat;
        c_DataBox(idxR) = {st_features(s_ChanIdx).(str_Var)'};
        c_TreatBox(idxR) ={num2str(s_CurrTreat)};
    end 
    figure ();
    
    violin(c_DataBox,'xlabel',c_TreatBox)
    title(strcat(num2str(s_currGroup),str_unitsG))
    ylabel(strcat(str_Var,str_unitsV))
    xlabel(strcat(str_nameT,str_unitsT))

end 
