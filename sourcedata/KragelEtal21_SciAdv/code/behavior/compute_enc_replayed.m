function [n_revisited, n_replayed] = compute_enc_replayed(events)
% COMPUTE_ENC_REPLAYED for a subject (events), computes the number
% of fixations that were revisitations or replays within the study
% period
%
% INPUTS:  events - structure, contains event information
%
% OUTPUTS: n_revisited - number of revisitation fixations
%          
%          n_replayed - number of replayed fixations

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

% filter to only encoding fixation events
enc_fix_events = events(strcmp('fixation', {events.type}) & ...
    strcmp('encode',{events.phase}));

n_revisited = nan(1, length(enc_fix_events));
n_replayed = nan(1, length(enc_fix_events));

for i = 1:length(enc_fix_events)
    
    % compute percent of enc_repeated == 0 (i.e., second time viewing the
    % same location), that are replayed at encoding (enc_replay==1)

    seq = enc_fix_events(i).enc_seq;
    enc_repeated = enc_fix_events(i).enc_repeated;

    r_idx = find(enc_repeated==0);
    
    is_repeated = true(size(r_idx));
    is_replayed = false(size(r_idx));
    for j = 1:length(r_idx)
        
        this_transition = seq(j:j+1);
        
        if diff(this_transition)==0
           is_repeated(j) = false; % it is repeated, but we want distinct positions
           continue
        end
        
        transitions = strfind(seq, this_transition); % trick with strfind
        is_replayed(j) = length(transitions)>1; % sequence occurs more than once
        
    end
    n_revisited(i) = sum(is_repeated); % number of non initial
    n_replayed(i) = sum(is_replayed);
    
end
