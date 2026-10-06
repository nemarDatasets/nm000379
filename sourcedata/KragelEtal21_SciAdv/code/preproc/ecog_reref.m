function reref_filename = ecog_reref(h5_file)
% ECOG_REREF computes bipolar pairs from sEEG in h5 format
%
% INPUTS: h5file - string, path to monopolar h5 file
%
% OUTPUTS: reref_filename - string, path to output bipolar h5 file

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



% assumes data are contact-REFX format

channels = h5read(h5_file, '/channels');
data     = h5read(h5_file, '/eeg');
srate     = h5read(h5_file, '/srate');

[contact, cno, i] = find_groups(channels);
[pairs, pair_cnos, out_idx] = create_bipolar(contact, cno, i);
reref_data = apply_reref(data, channels, out_idx);

pairs = [pairs, 'SYNC'];

% remove scalp and reference only pairs
to_remove = startsWith(pairs, 'Ref') | startsWith(pairs, 's');
pairs = pairs(~to_remove);
reref_data = reref_data(~to_remove,:);

reref_filename = to_h5(h5_file, reref_data, pairs, srate);

end


function filename = to_h5(orig_filename, data, pairs, srate)

[~, filestart] = fileparts(orig_filename);
pairs = char(deblank(pairs));
filename = [filestart '_bp.h5'];

if exist(filename, 'file')
    return
end

% channel info
string_to_h5(filename, 'pairs', pairs);

% eeg data
h5create(filename, '/eeg', size(data));
h5write(filename, '/eeg', data);

% sampling rate
h5create(filename, '/srate', size(srate));
h5write(filename, '/srate', srate);

end

function reref_data = apply_reref(data, channels, pair_cnos)

% find the sync pulses
is_sync = false(size(channels));
for c = 1:length(channels)
   if startsWith(strtrim(channels{c}),'DC') || startsWith(strtrim(channels{c}),'*DC')
       is_sync(c) = true;
   end
end

if sum(is_sync) > 1 || sum(is_sync) == 0
    error('Could not find correct sync channel');
end

reref_data = data(pair_cnos(:,1),:) - data(pair_cnos(:,2),:);

reref_data(end+1,:) = data(is_sync,:);

end

function [pairs, pair_cnos, out_idx] = create_bipolar(contact, cno, i)

[groups, ~, ib] = unique(i);

pairs = {}; 
pair_cnos = [];
out_idx = [];
orig_idx = 1:length(cno);

for g = 1:length(groups)
   idx = find(i == groups(g));
   
   for c = 1:length(idx)
      dist = cno(idx(c)) - cno(idx);
      
      for d = 1:length(dist)
          if dist(d) == -1
              pairs = [pairs [contact{idx(c)} '-' contact{idx(d)}]];
              pair_cnos = cat(1, pair_cnos, [cno(idx(c)) cno(idx(d))]);
              out_idx = cat(1, out_idx, [orig_idx(idx(c)) orig_idx(idx(d))]);
              
          end
      end
   end
   
end

end


function [contact, contactno, i] = find_groups(channels)
% FIND_GROUPS takes common ref and identifies groups by names

for c = 1:length(channels)
    
    fno = find( (channels{c} - '0') < 10 & ...
        (channels{c} - '0') >= 0, ...
        1, ...
        'first');
    
    sno = findstr(channels{c},'-');
    
    groups{c} = channels{c}(1:fno-1);
    contact{c} = strtrim(channels{c}(1:sno-1));
    contactno(c) = str2double(channels{c}(fno:sno-1));
    
end

[groups, ~, i] = unique(groups);

end