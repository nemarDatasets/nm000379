function [is_ltc, pairs] = return_ltc(subj, pairs)
% RETURN_LTC Returns logical specifying if bipolars pairs are within the 
%   lateral temporal cortex
%
% INPUTS:  subj - string, contains subject identifier
%
%          pairs - cell array, contains strings of bipolars pairs to evaluate
%
% OUTPUTS: is_ltc - logical, specifying if the pair is in the LTC
%
%          pairs - cell array, string containing bipolar pair names

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

if ~exist('pairs','var')
    res_dir = ['..\..\data\eeg\' subj filesep];
    files = dir([res_dir filesep '*_bp.h5']);
    pairs = h5read([res_dir filesep files(1).name], '/pairs');
end

basedir = '..\..\data\localization\';
[x, y, z, contacts] = read_sub_coords(basedir, subj);
is_ltc = find_ltc([x y z]);
ltc_contacts = contacts(is_ltc);

is_ltc = to_bipolar(pairs, ltc_contacts);


end

function is_region = to_bipolar(pairs, contacts)

% and convert to bipolar
is_region = false(size(pairs));
for i = 1:length(is_region)
    for c = 1:length(contacts)
        if startsWith(deblank(pairs{i}), [contacts{c} '-'])
            is_region(i) = true;
        elseif endsWith(deblank(pairs{i}), ['-' contacts{c}])
            is_region(i) = true;
        elseif startsWith(deblank(pairs{i}), [contacts{c}]) && endsWith(deblank(pairs{i}), [contacts{c}])
            is_region(i) = true;
        end
    end
end

end

function [x, y, z, channel] = read_sub_coords(basedir, sub)

tab = readtable([basedir '/' sub '/' sub '_mni.csv']);
x = tab.x;
y = tab.y;
z = tab.z;
channel = tab.channel;

end

function is_ltc = find_ltc(coords)

ref_vol = spm_vol('..\..\data\imaging\HarvardOxford-cort-maxprob-thr0-1mm.nii');
vox = [coords ones(size(coords,1),1)]*inv(ref_vol.mat)';
roi = spm_sample_vol(ref_vol,vox(:,1), vox(:,2), vox(:,3), 0);

is_ltc = false(size(roi));
for r = 1:length(roi)
    if (roi(r) < 16 && roi(r) >= 8)
        is_ltc(r,1) = true;
    end
end


end

