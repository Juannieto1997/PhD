function f_HeatMapMEA12x12(v_Data,c_Chann,v_caxis,str_title,str_label,str_SaveFile)

%% Prepare Data
% organize the channels in alphanumerical order
[c_Chann, v_idx] = f_organizeChannels (c_Chann);
% Recognize unique characters in channels
c_Let = {};
c_Num = {};
for idx = 1:length(v_idx)
    str_curr = c_Chann{idx}; 
    str_let = str_curr(1); 
    str_num = str_curr(2:3);
    if  isempty(find(contains(c_Let,str_let)))
        c_Let{end+1} = str_let; 
    end 
    if  isempty(find(contains(c_Num,str_num)))
        c_Num{end+1} = str_num; 
    end 
end 

[c_Let,~] = sortrows(c_Let');
[c_Num,~] = sortrows(c_Num');
% Organize the matriz acording to the previous index 
v_Data = v_Data(v_idx,:);
s_size = size(v_Data);
% Initialize Matrix 
m_ElecData = zeros(length(c_Num),length(c_Let));
%% Organize data in matrix 
for idxNum = 1:length(c_Num)
    for idxLet = 1:length(c_Let)
        str_currChann = strcat(c_Let{idxLet},c_Num{idxNum}); % Channel to look for 
        idxChan = find(contains(c_channels,str_currChann),1); 
        % Channel not found skip 
        if isempty(idxChan)
            continue
        end 
        m_ElecData(idxNum,idxLet) = v_Data(idxChan);
    end
end
