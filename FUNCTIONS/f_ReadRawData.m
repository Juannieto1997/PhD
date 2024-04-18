function [m_Data,stru_Header,s_SampRate] = f_ReadRawData(str_FileName)
% f_ReadRawData 
% Function by: David Henao 
% allows to read .RAW files from multichannels system 
% export files with heather in RAW format
%   INPUTS:
%       str_FileName = Name of the file as string including .raw in the name 
%   OUTPUTS:
%       m_Data = nxm Matrix containing all the data from the raw file, where n = Channels m = Data points
%       stru_header = structure containg the file header 
%       s_SampRate = sample rate of the file

%str_FileNameDataRaw = strcat(str_FileName,'.raw');
str_FileNameDataRaw = str_FileName;
% Get header
fileIDHeader = fopen(str_FileNameDataRaw,'r');
tline = fgetl(fileIDHeader);
disp('Reading header...')

while 1

    if contains(tline,'Sample')
        s_SampRate = str2double(tline(15:end));
    elseif contains(tline,'ADC')
        s_ADCzero = str2double(tline(11:end));
    elseif contains(tline,'V/AD')
        s_EI = str2double(tline(5:strfind(tline,'V/AD')-2));
    elseif contains(tline,'Streams')
        str_Chn = tline;
        [cll_Chan,~] = strsplit(str_Chn,';');
        s_NumElec = numel(cll_Chan);        
    elseif strcmp(tline,'EOH') == 1
        break
    end

    tline = fgetl(fileIDHeader);
end

fclose(fileIDHeader);
disp('Done!')

% Modify Chan 1 name

cll_Chan{1} = cll_Chan{1}(end-5:end);

% Build structure variable with header info

stru_Header.s_ADCzero = s_ADCzero;
stru_Header.s_EI = s_EI;
stru_Header.s_NumElec = s_NumElec;
stru_Header.ChannOrder = cll_Chan;

fprintf('\n');      

%% Get data

disp('Reading data...')

fileID = fopen(str_FileNameDataRaw,'r');

str_TextRead = fgetl(fileID);
cont = 1;
while strcmp(str_TextRead,'EOH')~=1

    fseek(fileID,cont,'bof');
    str_TextRead = fgetl(fileID);
    cont = cont+1;
end
cont = cont+4;

fseek(fileID,cont,'bof');
v_Data = fread(fileID,'uint16')';
v_Data = s_EI*(v_Data-s_ADCzero);
fclose(fileID);

disp('Done!')
fprintf('\n');

% Reshape data
m_Data = reshape(v_Data,[s_NumElec,numel(v_Data)/s_NumElec]);

end