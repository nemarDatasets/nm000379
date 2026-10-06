function [h_files, reref_files] = ascii_to_h5(ascii_dir)
% ASCII_TO_H5 converts raw data to h5 for easier loading
%
% INPUTS: ascii_dir - string, specifying directory with NK m00 files
%
% OUTPUTS: h_files - structure containing monopolar h5 files
%          reref_files - structure containing bipolar h5 files

% Copyright (C) 2021  James Kragel

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

m00_files = dir([ascii_dir '/*.m00']);
m00_files = {m00_files(:).name};

for f = 1:length(m00_files)
    
    
    [header,data]=hdrload(m00_files{f});
    data=data';
    
    spaces=findstr(' ',header(2,:));
    spaces(length(spaces)+1)=length(header(2,:));
    
    ct=1;
    channels = cell(length(spaces)-1,1);
    for n=1:length(spaces)-1
        channels{n}=header(2,spaces(n):spaces(n+1));
        ct=ct+1;
    end
    
    % sometimes ref is REF 1 instead of REF1, so spacing is off
    channels(strcmp(channels, ' 1 ')) = [];
    
    %the header can have leading spaces:
    if length(channels)==size(data,1)+1
        channels(1)=[];
    end
    
    %get the sampling rate:
    srate=1000/str2double(header(1,strfind(header(1,:),'SamplingInterval[ms]')+21));
    if srate==Inf
        srate=1000/str2double(header(1,strfind(header(1,:),'SamplingInterval[ms]')+21:strfind(header(1,:),'SamplingInterval[ms]')+24));
    end
    
    % some processing of channel names
    for c = 1:length(channels)
        channels{c} = strrep(channels{c},'*','');
        channels{c} = strrep(channels{c},' ','');
    end
   
    h_files{f} = to_h5(m00_files{f}, data, channels, srate);
    reref_files{f} = ecog_reref(h_files{f});
end

end

function filename = to_h5(ascii_filename, data, channels, srate)

[~, filestart] = fileparts(ascii_filename);
channels = char(deblank(channels));
filename = [filestart '.h5'];

if exist(filename, 'file')
    return
end

% channel info
string_to_h5(filename, 'channels', channels);

% eeg data
h5create(filename, '/eeg', size(data));
h5write(filename, '/eeg', data);

% sampling rate
h5create(filename, '/srate', size(srate));
h5write(filename, '/srate', srate);

end

