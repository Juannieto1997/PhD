function [c_ChanOrd, v_idx] = f_organizeChannels (c_Channels)
%% f_organizeChannels
% Function to organize the channels in a growing order acording to the
% nomenclature in the MEA system 120 grid. (Letter-ChannelNo)
% INPUTS: 
%   c_Channels  - Disorganized channels in MEA format ["EL_"LetterChannelNo] 
% OUTPUTLS:
%   c_ChanOrd   - Channels in the final order 
%   v_idx       - Index of the organized channels to organize the matriz. 
%% Validate there are no spaces

v_spidx = find(contains(c_Channels,' '));
for s_idx = 1:length(v_spidx)
    temp = c_Channels{s_idx};
    c_temp = split(temp,' '); 
    c_Channels{s_idx} = c_temp{end};
end 

%% Set all channels in the same format
% Organize so the data has the same lenght and single digits start with 0 
for idx = 1:length(c_Channels)
    str_temp = c_Channels{idx}; 
    c_part = split(str_temp,'_'); 
    str_temp = c_part{end}; 
    if length(str_temp) == 2
        str_temp(end+1) = str_temp(end);
        str_temp(end-1) = '0';
    end 
c_Channels{idx} = str_temp;
end 

%% Orgaize the data
[c_ChanOrd,v_idx] = sortrows(c_Channels');
if v_idx == 1
    v_idx = 1:length(c_ChanOrd);
end 