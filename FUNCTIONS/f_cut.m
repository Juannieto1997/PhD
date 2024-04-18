function [v_indBase,v_indResp,m_cut] = f_cut (v_data,v_locs,v_cutIdx)
%% m_cut
% FUNCTION: Cut a single channel according to the locations given by
% parameter, Caculates individual base line and stimulus response 
% INPUT: 
%   * v_data = vector containing the desired data to cut
%   * v_locs = Locations to cut 
%   * v_cutIdx = index to cut before and after the stimulation 
% OUTPUT: 
%   * m_cut = final matrix containing all the cut data
%   * v_indBase = 
%   * v_indResp = 
%% Define variables
% initialize the final matriz
m_cut = zeros (length(v_locs),sum(v_cutIdx)+1); 
v_indBase = zeros (1,length(v_locs));
v_indResp = zeros (1,length(v_locs));
for idxStim = 1:length(v_locs)
    s_currStim = v_locs(idxStim); % Get current stiulus 
    % Index of the activity of the stimulus
    %v_idxStim = [s_currStim-v_cutIdx(1)-3:s_currStim-3,s_currStim+15:s_currStim+v_cutIdx(2)+14];
    v_idxStim=[s_currStim-v_cutIdx(1):s_currStim+v_cutIdx(2)];
    % Add the corresponding information to the final matrix.
    if any(v_idxStim > length (v_data)) + any(v_idxStim < 0)
        continue
    end 
    v_cut = v_data(int64(v_idxStim));
    m_cut(idxStim,:) = v_cut-mean(v_cut);
    % Calculate baseline 
    s_base = mean(v_cut(1:5));
    v_indBase(idxStim) = s_base; 
    % Calculate response 
    s_reponse = max(v_cut(1,35:end));
    v_indResp(idxStim) = s_reponse;
end 