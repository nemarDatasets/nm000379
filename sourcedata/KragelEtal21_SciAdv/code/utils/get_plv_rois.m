function [is_hc, is_vis, is_dan, chansel, pairs] = get_plv_rois(subj)
% GET_PLV_ROIS helper function to return indices of different ROIs for PLV
% and PSI analyses
%
% INPUTS: subj - string, subject id
%
% OUTPUTS: is_hc   - logical, specifying hippocampal contacts
%          is_vis  - logical, specifying contains in vis. network
%          is_dan  - logical, specifying contains in dorsal attn. network
%          chansel - double, index of channel in original pairs
%          pairs   - string array, updated list of pair names

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

eeg_dir = ['..\..\data\eeg\' subj '\'];
fname = [eeg_dir subj '_sl_block1_bp.h5'];
    
pairs = h5read(fname, '/pairs');

for p = 1:length(pairs)
    pairs{p} = upper(pairs{p});
    pairs{p} = strrep(pairs{p},'*','');
end

is_hc = return_hc(subj, pairs);
try
    [is_dan, is_vis, pairs] = return_yeo(subj, pairs);
catch
    is_dan = false(size(pairs));
    is_vis = false(size(pairs));
end

to_keep = find(is_hc | is_dan | is_vis);

is_hc = is_hc(to_keep);
is_dan = is_dan(to_keep);
is_vis = is_vis(to_keep);

chansel = find_connections(subj, pairs(to_keep));

pairs = pairs(to_keep);

end


function chansel = find_connections(subj, pairs)


nchan = length(pairs);
chanindx = tril(true(nchan), -1);
cmbindx1 = repmat((1:nchan)', [1 nchan]);
cmbindx2 = repmat((1:nchan),  [nchan 1]);

chansel = [cmbindx1(chanindx) cmbindx2(chanindx)];

is_hc = return_hc(subj, pairs);

% only find combinations that are between Hc and DAN/VIS, (not within Hc and
% within DAN/VIS)

between_region = false(1, size(chansel,1));
for c = 1:size(chansel,1)
    if sum(is_hc(chansel(c,:))) == 1
        between_region(c) = true;
    end
end

chansel = chansel(between_region,:);

end

