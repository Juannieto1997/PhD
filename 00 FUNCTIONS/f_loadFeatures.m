function v_Features = f_loadFeatures(str_folder,c_files)
v_Features = [];
for idxfile = 1:length(c_files)
    str_currFile = c_files{idxfile}; 
    str_file = fullfile(str_folder,str_currFile); 
    load(str_file)
    v_Features = [v_Features v_eFPAmp];
end 