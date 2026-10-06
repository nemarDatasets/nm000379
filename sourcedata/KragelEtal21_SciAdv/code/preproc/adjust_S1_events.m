function [events, adjust] = adjust_S1_events()
% ADJUST_S1_EVENTS adjusts events due to different timing of sync pulses at retrieval
%
% INPUTS: none
%
% OUTPUTS: events - structure, containing behavioral and eye events for S1
%          adjust - double, contains offset to adjust timing of events

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

events = get_events('S1');

ub = unique([events.block]);

adjust  = nan(length(ub), 48);
for b = 1:length(ub)
    if b == 1
        [n_events, adjust] = adjust_rec_onsets(events([events.block]==ub(b)));
    else
        [ev, adjust(b,:)] = adjust_rec_onsets(events([events.block]==ub(b)));
        n_events = [n_events ev];
        
    end
end
events = n_events;

end

function [events, adjust] = adjust_rec_onsets(events)

stim_idx = [events.stim_idx];
ids = unique(stim_idx);

for i = 1:length(ids)
    
    % and for recognition trials
    rec_trial_idx = strcmp({events.phase},'recog') & ...
        [events.stim_idx] == ids(i) & ...
        strcmp({events.type},'trial');
    
    if ~any(rec_trial_idx)
        continue
    end
    
    adjust(i) = (events(rec_trial_idx).end_time - events(rec_trial_idx).start_time) - ...
        (events(rec_trial_idx).rt + 800);
    
    events(rec_trial_idx).start_time = events(rec_trial_idx).start_time + adjust(i);
    
    rec_fix_idx = strcmp({events.phase},'recog') & ...
        [events.stim_idx] == ids(i) & ...
        strcmp({events.type},'fixation');
    
    if any(rec_fix_idx)
        events(rec_fix_idx).start_time = events(rec_fix_idx).start_time - adjust(i);
    end
    
    rec_sac_idx = strcmp({events.phase},'recog') & ...
        [events.stim_idx] == ids(i) & ...
        strcmp({events.type},'saccade');
    
    if any(rec_sac_idx)
        events(rec_sac_idx).start_time = events(rec_sac_idx).start_time - adjust(i);
    end
    
    rec_blink_idx = strcmp({events.phase},'recog') & ...
        [events.stim_idx] == ids(i) & ...
        strcmp({events.type},'blink');
    
    if any(rec_blink_idx)
        events(rec_blink_idx).start_time = events(rec_blink_idx).start_time - adjust(i);
    end
end


end



