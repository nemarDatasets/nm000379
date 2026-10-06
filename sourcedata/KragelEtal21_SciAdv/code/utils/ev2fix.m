function fix_info = ev2fix(fix_events, add_enc_lag, ex_fields)
% EV2FIX extracts trial level information for each fixation event
%
% INPUTS:  fix_events - string, specifies subject id to load data 
%
%          add_enc_lag - structure, bosc output structure from run_subject_bosc.m
%
%          ex_fields - structure, fields to include in output fixation structure
%
% OUTPUTS: fix_info - structure, contains vector format information for fields of interest

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

if ~exist('add_enc_lag', 'var')
    add_enc_lag = false;
end    

if ~exist('ex_fields', 'var')
    ex_fields = {'block', 'trial', 'old', 'salience', 'resp', 'rt'};
end

fix_info = struct();

for f = 1:length(ex_fields)
    fix_info.(ex_fields{f}) = build_fix_data(fix_events, ex_fields{f});
end

% for encoding, also taking reins, lag_one
fix_info.('viewed') = append_fixation_level(fix_events, 'reins');
fix_info.('lag_one') = append_fixation_level(fix_events, 'lag_one');
fix_info.('lag_minusone') = append_fixation_level(fix_events, 'lag_minusone');
fix_info.('precue') = append_fixation_level(fix_events, 'precue');

fix_info.('seq') = append_fixation_level(fix_events, 'enc_seq');
fix_info.('repeated') = append_fixation_level(fix_events, 'enc_repeated');
fix_info.('replay') = append_fixation_level(fix_events, 'enc_replay');
fix_info.('rev_replay') = append_fixation_level(fix_events, 'enc_rev_replay');

% for retrieval, also taking recall, and lag
fix_info.('recall') = append_fixation_level(fix_events, 'recall');
fix_info.('lag') = append_fixation_level(fix_events, 'lag');
fix_info.('lagshift') = append_lagshift(fix_events);

% for retrieval, ad repeat from recall sequence
fix_info.('rec_repeated') = append_repeat(fix_events);

% add fixation number
fix_info.('fix_no') = append_fixno(fix_events, false);
fix_info.('reverse_fix_no') = append_fixno(fix_events, true);

% add encoding lag
if add_enc_lag
    fix_info.('lag') = append_enclag(fix_events);    
end

% whether previous fixations was revisitation/repeated
fix_info.('prev_repeated') = append_prev_fixation_level(fix_events, 'enc_repeated');

end

function out = append_enclag(fix_events)

out = [];
for e = 1:length(fix_events)
    n_fix = size(fix_events(e).start_time,2);
    
    if isempty(fix_events(e).enc_seq)
        out = [out nan(1, n_fix)];
    else
        out = [out diff(fix_events(e).enc_seq) nan];
    end
end

end

function out = append_fixno(fix_events, rev_flag)

out = [];
for e = 1:length(fix_events)
    n_fix = size(fix_events(e).start_time,2);
    if rev_flag % reverse order, from end of trial
        out = [out fliplr(1:n_fix)];
    else
        out = [out 1:n_fix];
    end
end

end

function out = append_lagshift(fix_events)

out = [];

for e = 1:length(fix_events)
   n_fix = size(fix_events(e).start_time,2);
   
    if isempty(fix_events(e).lag)
        out = [out nan(1, n_fix)];
    else
        out = [out nan fix_events(e).lag(1:end-1)];
    end
end

end

function out = append_repeat(fix_events)

out = [];

for e = 1:length(fix_events)
   n_fix = size(fix_events(e).start_time,2);
   
    if isempty(fix_events(e).lag)
        out = [out nan(1, n_fix)];
    else
        out = [out nan diff(fix_events(e).recall)];
    end
end

end

function out = append_fixation_level(fix_events, fn)

out = [];
for e = 1:length(fix_events)
    n_fix = size(fix_events(e).start_time,2);
    
    if isempty(fix_events(e).(fn))
        out = [out nan(1, n_fix)];
    else
        out = [out fix_events(e).(fn)];
    end
end


end


function out = append_prev_fixation_level(fix_events, fn)

out = [];
for e = 1:length(fix_events)
    n_fix = size(fix_events(e).start_time,2);
    
    if isempty(fix_events(e).(fn))
        out = [out nan(1, n_fix)];
    else
        out = [out nan fix_events(e).(fn)(1:end-1)];
    end
end


end

function out = build_fix_data(fix_events, fn)

out = [];
for e = 1:length(fix_events)
    n_fix = size(fix_events(e).start_time,2);
    out = [out fix_events(e).(fn)*ones(1, n_fix)];
end


end

