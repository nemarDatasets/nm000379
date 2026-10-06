function data = update_trialinfo(data, fix_events)
% UPDATE_TRIALINFO appends revisit and subsequent memory information for fieldtrip data
%
% INPUTS: data - fieldtrip data structure containing time-frequency data
%				 to compare
%
%		  events - structure, containing trial level event information
%
%		  fix_events - structure, containing fixation level envent information
%
%
% OUTPUTS: stats - fieldtrip structure containing permutations stats

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

% reorganize information in fixation even structure to fixation level
% observations
fix_info = ev2fix(fix_events, true);

% find where fixation events are in data
fix_idx = data.trialinfo(:,1)==2; % fixations

% init subsequent memory, and revisit vs. other
data.trialinfo(:, 2:3) = nan;
data.trialinfo(fix_idx, 2) = fix_info.resp==1; 	     % subsequent hit vs. miss
data.trialinfo(fix_idx, 3) = fix_info.repeated == 0; % revisit (not initial viewing)
data.trialinfo(isnan(fix_info.repeated), 3) = 2;     % other fixation

end

