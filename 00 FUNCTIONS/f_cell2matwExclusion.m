function v_data = f_cell2matwExclusion (c_Data,c_chann,c_exclude)
b_remove = contains(c_chann,c_exclude);
c_Data = {c_Data{~b_remove}};
v_data = [];
for idxChann = 1:length(c_Data)
    v_data = [v_data c_Data{idxChann}];
end 

