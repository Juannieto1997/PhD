function f_rasterplot (c_Locs,s_SRate)

v_allspikes = [];
v_axes = [0 0];
figure();
v_axes(1) = subplot(3,1,[1,2]);
title('Rasterplot')
ylabel('Channels')
xlabel('Time (s)')
for idxchann = 1: length(c_Locs)
    c_Spikes = c_Locs{idxchann}/s_SRate; 
    m_Spikesx = repmat(c_Spikes,3,1);
    m_Spikesy = nan(size(m_Spikesx));

    if ~isempty(m_Spikesy)
        m_Spikesy(1,:) = idxchann -1;
        m_Spikesy(2,:) = idxchann;
    end 

    plot(m_Spikesx,m_Spikesy,'k')
    v_allspikes = [v_allspikes c_Spikes];
end 
v_axes(2) = subplot(3,1,3);
title('histogram')
ylabel('Peaks')
xlabel('Time (s)')
hist(v_allspikes)

linkaxes(v_axes,'x')