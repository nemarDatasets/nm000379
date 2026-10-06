function res = get_bosc(subj, base_dir, vis_flag, ltc_flag)
% GET_BOSC helper function to load bosc results structure
%
% INPUTS:  subj - string, subject identifier
%
%		   base_dir - string, path to study directory
%
%		   vis_flag - logical, flag specifying whether to load DAN\VN data
%
%		   ltc_flag - logical, flag specifying whether to load LTC data
%
% OUTPUTS: res - structure, bosc results for this subject/roi pair

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


if ~exist('base_dir', 'var')
    base_dir = '..\..\scratch_data\';
end

if ~exist('vis_flag','var')
    vis_flag = false;
end

if ~exist('ltc_flag','var')
    ltc_flag = false;
end

bosc_dir = [base_dir 'bosc\'];

if ~vis_flag && ~ltc_flag
    load([bosc_dir subj '_bosc.mat'],'res');
elseif vis_flag
    load([bosc_dir subj '_vis_bosc.mat'],'res');
elseif ltc_flag
    load([bosc_dir subj '_ltc_bosc.mat'],'res');    
end

