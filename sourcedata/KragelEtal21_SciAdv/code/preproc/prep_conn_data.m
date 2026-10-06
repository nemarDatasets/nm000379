function prep_conn_data()
% PREP_CONN_DATA - prepares timeseries data for analysis in connectivity analysis
%
% INPUTS: none - simple wrapper for epoch_data
%
% OUTPUTS: none - outputs saved to temporary dir specified in epoch_data.m

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

include_vis = true; % can specify whether you want to include VN contacts

subjects = {'S1', 'S2', 'S4', 'S5', 'S6'};

for s = 1:length(subjects)

	% get events to analyze
	events = get_events(subjects{s});

	% saves processed data to folder specified in epoch_data
	epoch_data(subjects{s}, events, include_vis);
    
end
