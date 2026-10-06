function fs = get_srate(subj, base_dir)
% GET_SRATE helper function to load sampling rate rather than reading header
% from m00 which is time consuming
%
% INPUTS:  subj - string, subject identified
%
%		   base_dir - string, path to study directory
%
% OUTPUTS: fs - double, sampling rate


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
    base_dir = '../../data/';
end

h5_dir = [base_dir 'eeg/' subj '/'];
a = dir([h5_dir '/*_bp.h5']);

fs = h5read([h5_dir filesep a(1).name], '/srate');

end

