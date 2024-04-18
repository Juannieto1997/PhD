function [c_init] = f_AppendCells (c_init,c_append)
%% f_AppendCells
% Function by: Juan Nieto
% Simple for cicle to append cells into a new cell. 
for idx = 1:length(c_append)
    c_init{end+1} = c_append{idx};
end 