function  [m_MorseWT,v_TFTime,v_FreqAxis] = f_TimeFreq (v_Data,v_TimeLims,s_SampRate,v_Lims)

s_TimeInf = v_TimeLims(1);
s_TimeSup = v_TimeLims(2);

s_MinFreqHz = v_Lims(1);
s_MaxFreqHz = v_Lims(2);
% s_FreqSeg = 128;
% s_NumOfCycles = 3.2;

s_FreqSeg = 100;
s_NumOfCycles = 1.5;

s_Magnitudes = 1;
s_SquaredMag = 0;
s_MakeBandAve = 0;
s_Phases = 0;
s_TimeStep = [];
s_TimeSupIni = s_TimeSup;

if s_TimeInf == 0

    s_SampInf = 1;
    s_SampSup = s_TimeSup*s_SampRate;
    s_Cut = 0;

elseif s_TimeInf >= (numel(v_Data)-1)/s_SampRate

    s_SampSup = numel(v_Data);
    s_SampInf = s_SampSup-s_TimeSupIni*s_SampRate;
    s_Cut = 0;
else

    s_SampInf = (s_TimeInf*s_SampRate)-2000;                               % To remove border effect
    s_SampSup = (s_TimeSup*s_SampRate)+2000;                               % To remove border effect
    s_Cut = 1;
end

v_SegWav = v_Data(s_SampInf:s_SampSup);

[m_MorseWT,~,v_FreqAxis] = ...
    f_MorseAWTransformMatlab(...
    v_SegWav, ...
    s_SampRate, ...
    s_MinFreqHz, ...
    s_MaxFreqHz, ...
    s_FreqSeg, ...
    s_NumOfCycles, ...
    s_Magnitudes, ...
    s_SquaredMag, ...
    s_MakeBandAve, ...
    s_Phases, ...
    s_TimeStep);

if s_Cut == 1

    m_MorseWT(:,end-2000:end)=[];
    m_MorseWT(:,1:2000)=[];

end

[~,s_Samp] =size(m_MorseWT);
v_TFTime = linspace(s_TimeInf,s_TimeSup,s_Samp);

end

