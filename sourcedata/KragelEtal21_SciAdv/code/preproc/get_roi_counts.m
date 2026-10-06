function counts = get_roi_counts(subject)
% GET_ROI_COUNTS - counts contacts in various regions of interest
%
% INPUTS: subjects - string, specifying subject id
%
% OUTPUTS: counts - structure containing number of contacts per ROI

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


mtl_flag = true;
basedir = '..\..\data\localization\';

[x, y, z, hc_idx, channel] = read_sub_coords(basedir, subject, mtl_flag);
is_hc = hc_idx;

% find rois
n_elec = length(x);
is_vis = nan(1, length(n_elec));
is_dan = is_vis;
is_both = is_vis;
roi = is_vis;

for e = 1:n_elec
    [is_vis(e), is_dan(e), is_ltc(e), roi(e)] = get_roi([x(e) y(e) z(e)]);
end

is_ltc(hc_idx==1) = false;
is_dan(hc_idx==1) = false;
is_vis(hc_idx==1) = false;

counts.hc = sum(is_hc);
counts.ltc = sum(is_ltc);
counts.dan = sum(is_dan);
counts.vis = sum(is_vis);

end

function [x, y, z, hc_flag, channel] = read_sub_coords(basedir, sub, mtl_flag)

tab = readtable([basedir '/' sub '/' sub '_mni.csv']);
x = tab.x;
y = tab.y;
z = tab.z;
channel = tab.channel;

hc_flag = [];
amy_flag = [];

if mtl_flag
    hc_flag = return_hc(sub, tab.channel);
end

end

function [is_vis, is_dan, is_ltc, roi] = get_roi(coords)

% harvard oxford
ref_vol = spm_vol('..\..\data\imaging\HarvardOxford-cort-maxprob-thr0-1mm.nii');
vox = [coords ones(size(coords,1),1)]*inv(ref_vol.mat)';
roi = spm_sample_vol(ref_vol,vox(:,1), vox(:,2), vox(:,3), 0);

% YEO
ref_vol = spm_vol('..\..\data\imaging\Yeo2011_7Networks.nii');
vox = [coords 1]*inv(ref_vol.mat)';
yeo_roi = spm_sample_vol(ref_vol,vox(1), vox(2), vox(3), 0);

is_dan = yeo_roi == 3;
is_vis = yeo_roi == 1;

for r = 1:length(roi)
    if roi(r) < 16 && roi(r) >= 8
        is_ltc(r) = true;
    else
        is_ltc(r) = false;
    end
end


end
