function [is_dan, is_vis, pairs] = return_yeo(subj, pairs)
% RETURN_YEO Returns logical specifying if bipolars pairs are within the DAN
%   or VN per Yeo et al 2011
%
% INPUTS:  subj - string, contains subject identifier
%
%          pairs - cell array, contains strings of bipolars pairs to evaluate
%
% OUTPUTS: is_dan - logical, specifying if the pair is in the DAN
%
%          is_vn - logical, specifying if the pair is in the VN
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
[is_dan, is_vis] = find_yeo([x y z]);

dan_contacts = contacts(is_dan);
vis_contacts = contacts(is_vis);

is_dan = to_bipolar(pairs, dan_contacts);
is_vis = to_bipolar(pairs, vis_contacts);

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

function [is_dan, is_vis] = find_yeo(coords)

% YEO
ref_vol = spm_vol('..\..\data\imaging\Yeo2011_7Networks.nii');
vox = [coords ones(size(coords,1),1)]*inv(ref_vol.mat)';
yeo_roi = spm_sample_vol(ref_vol,vox(:,1), vox(:,2), vox(:,3), 0);

is_dan = yeo_roi == 3;
is_vis = yeo_roi == 1;

end

