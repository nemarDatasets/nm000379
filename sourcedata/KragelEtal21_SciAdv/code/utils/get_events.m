function events = get_events(subj, update_flag, base_dir)
% GET_EVENTS helper function to load events for a subject
%
% INPUTS:  subj - string, subject identifier
%
%          update_flag - logical, flag specifying whether to generate new events
%
%          base_dir - string, path to study directory
%
% OUTPUTS: events - structure, contains event information for trials and eye movements

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

if ~exist('update_flag','var')
    update_flag = false;
end

if ~exist('base_dir', 'var')
    base_dir = '..\..\scratch_data\';
end

ev_dir = [base_dir 'events\'];
load([ev_dir subj '_events.mat'],'events');

if update_flag

    ub = unique([events.block]);
    
    for b = 1:length(ub)
        if b == 1
            n_events = append_submem(events([events.block]==ub(b)));
        else
            n_events = [n_events append_submem(events([events.block]==ub(b)))];
        end
    end
    events = n_events; clear n_events

    save([ev_dir subj '_events.mat'],'events');

end

