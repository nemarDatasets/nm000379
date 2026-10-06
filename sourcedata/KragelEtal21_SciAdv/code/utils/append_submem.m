function events = append_submem(events)
% APPEND_SUBMEM for all events for the study, appends information to
% encoding events (saccades, fixations, and trials) about the behavioral
% outcomes at retrieval
%
% INPUTS:  events - structure, contains event information
%
% OUTPUTS: events - structure, contains updated event information with 
%                   additional fields

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

stim_idx = [events.stim_idx];
ids = unique(stim_idx);

for i = 1:length(ids)

    % identify the relevant events
    enc_idx = find(strcmp({events.phase},'encode') & ...
        [events.stim_idx] == ids(i));
    
    rec_idx = find(strcmp({events.phase},'recog') & ...
        [events.stim_idx] == ids(i) & ...
        strcmp({events.type},'fixation'));
    
    enc_trial_idx = enc_idx(ismember(enc_idx, ...
        find(strcmp({events.type},'trial'))));
    
    enc_fix_idx = enc_idx(ismember(enc_idx, ...
        find(strcmp({events.type},'fixation'))));
    
    if ~isempty(enc_idx) && ~isempty(enc_fix_idx) && ~isempty(enc_trial_idx)
        
        enc_fix_onsets = events(enc_fix_idx).start_time - ...
            events(enc_trial_idx).start_time; % these are in ms
        
        events(enc_fix_idx).precue = enc_fix_onsets < 800;
        
    end
    
    % and for recognition trials
    rec_trial_idx = find(strcmp({events.phase},'recog') & ...
        [events.stim_idx] == ids(i) & ...
        strcmp({events.type},'trial'));
    
    rec_fix_idx = rec_idx(ismember(rec_idx, ...
        find(strcmp({events.type},'fixation'))));

    rec_fix_onsets = events(rec_fix_idx).start_time - ...
        events(rec_trial_idx).start_time; % these are in ms
    
    events(rec_fix_idx).precue = rec_fix_onsets < 800;
    
    if ~isempty(enc_idx) && ~isempty(enc_fix_idx) % novel item
        % append temporal reinstatement/recall order information
        
        res = spotlight_cvp(events([enc_fix_idx rec_idx]));
        
        n_enc_fix = length(events(enc_fix_idx).mean_x);
        events(enc_fix_idx).reins = res.enc_reins(1:n_enc_fix);
        events(enc_fix_idx).lag_one = res.lag_one(1:n_enc_fix);
        events(enc_fix_idx).lag_minusone = res.lag_minusone(1:n_enc_fix);
        events(enc_fix_idx).reins = res.enc_reins(1:n_enc_fix);
        
        n_fix = length(events(rec_idx).mean_x);
        events(rec_idx).recall = res.recalls(1:n_fix);
        events(rec_idx).lag = res.lags(1:n_fix);
        
        % append information about replay during encoding to encoding
        % events
        res = spotlight_enc_cvp(events(enc_fix_idx));
        events(enc_fix_idx).enc_seq = res.recalls(1:n_enc_fix);
        events(enc_fix_idx).enc_repeated = res.enc_reins(1:n_enc_fix);
        events(enc_fix_idx).enc_replay = res.lag_one(1:n_enc_fix);
        events(enc_fix_idx).enc_rev_replay = res.lag_minusone(1:n_enc_fix);
    end
    
    n_fix = length(events(rec_idx).mean_x);
    resp  = events(rec_idx).resp;
    rt  = events(rec_idx).rt;
    
    for e = 1:length(enc_idx)
       events(enc_idx(e)).rt = rt;
       events(enc_idx(e)).resp = resp;
       events(enc_idx(e)).n_fix = n_fix;
    end
    
    % update n_fix on all rec events
    rec_idx = find(strcmp({events.phase},'recog') & ...
        [events.stim_idx] == ids(i));
    
    for e = 1:length(rec_idx)
        events(rec_idx(e)).n_fix = n_fix;
        if isempty(enc_idx)
            events(rec_idx(e)).recall = nan(1,n_fix);
            events(rec_idx(e)).lag = nan(1,n_fix);
        end
    end
    
end

end
