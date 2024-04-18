function v_Filter = v_FFTfilt (v_data,s_SampRate,v_Rfreq)
%% v_FFTfilt
% function to make a fast fourier filter in a given set of data (TO BE
% VALIDATED)
%   INPUTS: 
%       v_data = Data to be filtered 
%       s_SampRate = sampling frequency 
%       v_freq = vector containing the frequency to remove, [low hig]
%   OUTPUT:
%       v_Filter = Filtered data with an fft Filter
%% Define the required variables 
v_freq = ([0:length(v_data)-1]/(length(v_data)-1))*s_SampRate; %frequency vector 
% frequencies to remove from the FFT 
v_ref = s_SampRate-flip(v_Rfreq);
v_remove = logical((v_freq > v_Rfreq(1)) .* (v_freq < v_Rfreq(2)) + (v_freq > v_ref(1)) .* (v_freq < v_ref(2)));
% calculate the FFT
v_FFT = fft(v_data);
% Remove the selected data 
v_FFT(v_remove) = v_FFT(v_remove)*10*exp(-10);
% Calculate the inverse FFT
v_Filter = real(ifft(v_FFT));
