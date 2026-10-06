function events = create_events(edf_dir, data_dir, subjid)
% CREATE_EVENTS - takes the eyelink data format and behav output and
% constructs events for behavioral/neural analysis
%
% INPUTS: edf_dir  - string, specifying directory with eyelink data
%         data_dir - string, specifying directory with behavioral data
%         subjid   - string, specifying subject id 
%
% OUTPUTS: events - structure containing all events necessary for analysis

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


% some constants for the study - trials per phase
n_enc = 24;
n_rec = 48;

% read in the relevant information from the data_dir
split_files = false;
try
    enc_info = dlmread([data_dir subjid '_enc_array.txt'],'\t');
    rec_info = dlmread([data_dir subjid '_rec_array.txt'],'\t');
    stim_info = load([data_dir subjid '_stim_array.txt']);
catch
    split_files = true; % for first pt
end

behav_output = dlmread([data_dir subjid '_output_array.txt'],'\t');
% read the recognition responses
recog_resp = behav_output(behav_output(:,5)~=0,5);
recog_rt = behav_output(behav_output(:,5)~=0,2);
recog_block = behav_output(behav_output(:,5)~=0,6);

if mod(length(recog_resp),48) ~=0
    warning('Check recognition responses in behavioral output.');
end

% for each 'trial' in the file, get the continuous data of interest

% take the output files and construct the information we need

n_trials_completed = size(behav_output,1);
n_blocks_completed = n_trials_completed/(n_enc+n_rec);

for b = 1:n_blocks_completed
       
    if split_files
        if b < 4
            enc_info = dlmread([data_dir subjid '_enc_array_s1.txt'],'\t');
            rec_info = dlmread([data_dir subjid '_rec_array_s1.txt'],'\t');
            stim_info = load([data_dir subjid '_stim_array_s1.txt']);
        else
            enc_info = dlmread([data_dir subjid '_enc_array_s2.txt'],'\t');
            rec_info = dlmread([data_dir subjid '_rec_array_s2.txt'],'\t');
            stim_info = load([data_dir subjid '_stim_array_s2.txt']);
        end
    end
    
    % encoding information for this block
    edf_struct = edfmex([edf_dir '/' subjid '_s' num2str(b) '.edf']);
    
    enc_duration = 3000;
    start_times = [edf_struct.RECORDINGS([edf_struct.RECORDINGS.state]==1).time];
    stop_times = [edf_struct.RECORDINGS([edf_struct.RECORDINGS.state]==0).time];
    durations = stop_times - start_times;
    delay = durations - enc_duration; 
    %delay = zeros(size(delay)); % include fixation period

    % parse the edf struct for actual events
    [trial, fix, saccade, blink] = get_eye_events(edf_struct, b, delay, 'encode');
    
    % recognition information for this block
    edf_struct = edfmex([edf_dir '/' subjid '_r' num2str(b) '.edf']);
    
    rec_duration = recog_rt(recog_block==b);
    start_times = [edf_struct.RECORDINGS([edf_struct.RECORDINGS.state]==1).time];
    stop_times = [edf_struct.RECORDINGS([edf_struct.RECORDINGS.state]==0).time];
    durations = stop_times - start_times;
    durations = double(durations');
    delay = durations - rec_duration;
    
    % parse the edf struct for actual events
    [tt, tf, ts, tb] = get_eye_events(edf_struct, b, delay, 'recog');
    
    fix = [fix tf];
    saccade = [saccade ts];
    blink = [blink tb];
    trial = [trial tt];
   
%   enc and rec_block used to find where to look in output
    enc_block = nan(1,n_enc*b);
    cnt=1;
    for bb = 1:b
        enc_block(cnt:cnt+n_enc-1) = bb;
        cnt=cnt+n_enc;
    end

    rec_block = nan(1,n_rec*b);
    cnt=1;
    for bb = 1:b
        rec_block(cnt:cnt+n_rec-1) = bb;
        cnt=cnt+n_rec;
    end  
    
    if b == 1
        events = [];
    end
    
    events = merge_events(events, ...
        trial, ...
        fix, ...
        blink, ...
        saccade, ...
        enc_info, ...
        enc_block, ...
        rec_info, ...
        rec_block, ...
        recog_resp, ...
        recog_rt, ...
        stim_info);
    
    
end

end

function [trial, fix, saccade, blink] = get_eye_events(edf_struct, block, ...
    delay, phase)
% GET_EYE_EVENTS returns eye events that occur during recordings in
% edf_struct

% events occuring before 'delay' in the trial are excluded. This is to
% potentially exclude eye events that occur during ISI or pre-cue periods

% delay should be a vector of times per trial

% filter the recording logs to start and stop events
start_times = [edf_struct.RECORDINGS([edf_struct.RECORDINGS.state]==1).time];
stop_times  = [edf_struct.RECORDINGS([edf_struct.RECORDINGS.state]==0).time];

% compute three eye event structures containing information per event

[trial, fix, saccade, blink] = init_eye_structs(length(start_times), block, phase);

for t = 1:length(start_times)
    
    % fixations - start or end during the trial
    fix_idx = (([edf_struct.FEVENT.entime] >= start_times(t) + delay(t) & ...
        [edf_struct.FEVENT.entime] < stop_times(t)) | ...
        ([edf_struct.FEVENT.sttime] >= start_times(t)  + delay(t) & ...
        [edf_struct.FEVENT.sttime] < stop_times(t))) & ...
        strcmp({edf_struct.FEVENT.codestring}, 'ENDFIX');
    
    fix(t) = populate_fields(fix(t), edf_struct, fix_idx);
    
    % saccades - start or end during the trial
    sac_idx = (([edf_struct.FEVENT.entime] >= start_times(t) + delay(t) & ...
        [edf_struct.FEVENT.entime] < stop_times(t)) | ...
        ([edf_struct.FEVENT.sttime] >= start_times(t) + delay(t) & ...
        [edf_struct.FEVENT.sttime] < stop_times(t))) & ...
        strcmp({edf_struct.FEVENT.codestring}, 'ENDSACC');
    
    % remove all ENDSACC produce by blinks
    sblink_idx = false(size(sac_idx));
    si = find(sac_idx);
    for i = 1:length(si)
       if contains(edf_struct.FEVENT(si(i)-1).codestring, 'BLINK')
           sblink_idx(si(i)) = true;
       end
    end
    sac_idx(sblink_idx)=false;
    
    saccade(t) = populate_fields(saccade(t), edf_struct, sac_idx);
    
    % blinks - start or end during the trial
    blink_idx = (([edf_struct.FEVENT.entime] >= start_times(t) + delay(t) & ...
        [edf_struct.FEVENT.entime] < stop_times(t)) | ...
        ([edf_struct.FEVENT.sttime] >= start_times(t) + delay(t) & ...
        [edf_struct.FEVENT.sttime] < stop_times(t))) & ...
        strcmp({edf_struct.FEVENT.codestring}, 'ENDBLINK');
    
    blink(t) = populate_fields(blink(t), edf_struct, blink_idx);
    
    % trial - start or end during the trial
    trial_idx = (([edf_struct.FEVENT.entime] >= start_times(t)-100 & ...
        [edf_struct.FEVENT.entime] < stop_times(t)) | ...
        ([edf_struct.FEVENT.sttime] >= start_times(t)-100 & ...
        [edf_struct.FEVENT.sttime] < stop_times(t))) & ...
        strcmp({edf_struct.FEVENT.codestring}, 'MESSAGEEVENT') & ...
        strcmp({edf_struct.FEVENT.message}, phase);
    
    trial(t) = populate_fields(trial(t), edf_struct, trial_idx);
    
end

end

function [trial, fix, saccade, blink] = init_eye_structs(n_events, block, phase)
% initialize output structures

for i = 1:n_events
    
    fix(i) = struct('phase', phase, ...
        'block', block, ...
        'trial', i, ...
        'start_time', nan, ...
        'end_time', nan, ....
        'start_x', nan, ...
        'start_y', nan, ...
        'end_x', nan, ...
        'end_y', nan, ...
        'mean_x', nan, ...
        'mean_y', nan, ...
        'sal_s', nan, ...
        'eye', 1, ...
        'n_fix', nan);
end

saccade = fix;
blink = fix;
trial = fix;

end

function ev_struct = populate_fields(ev_struct, edf_struct, idx)

ev_struct.start_time = [edf_struct.FEVENT(idx).sttime];
ev_struct.end_time = [edf_struct.FEVENT(idx).entime];
ev_struct.start_x = [edf_struct.FEVENT(idx).gstx];
ev_struct.start_y = [edf_struct.FEVENT(idx).gsty];
ev_struct.end_x = [edf_struct.FEVENT(idx).genx];
ev_struct.end_y = [edf_struct.FEVENT(idx).geny];
ev_struct.mean_x = [edf_struct.FEVENT(idx).gavx];
ev_struct.mean_y = [edf_struct.FEVENT(idx).gavy];

end

function events = merge_events(prev_events, trial, fix, blink, saccade, ...
    enc_info, enc_block, rec_info, rec_block, recog_resp, recog_rt, stim_info)

[trial.type] = deal('trial');
[fix.type] = deal('fixation');
[blink.type] = deal('blink');
[saccade.type] = deal('saccade');

events = [trial fix blink saccade];

r_idx = false(1,length(events));
for e = 1:length(events)
    if isempty(events(e).start_time)
        r_idx(e) = true;
    end
end

events(r_idx) = [];

% sort by first start_time
st = nan(1,length(events));
for e = 1:length(events)
    st(e) = events(e).start_time(1);
end

[~, i] = sort(st); % reorder by the start of each event
events = events(i);

ub = unique([events.block]);

for b = 1:length(ub)
    
    % update encoding events
    enc_out_idx = find(enc_block == ub(b));
    
    % this only works for full data, not per block
    enc_out_trial = 1:length(enc_out_idx); 
    
    enc_trial_idx = find([events.block]==ub(b) ...
        & strcmp({events.phase},'encode'));
    
    for e = 1:length(enc_trial_idx)
        
       idx = enc_out_idx(enc_out_trial == events(enc_trial_idx(e)).trial);
       
       events(enc_trial_idx(e)).stim_idx  = stim_info(enc_info( idx, 2), ...
                                                      enc_info( idx, 1));
       
       events(enc_trial_idx(e)).condition = enc_info( idx, 2);
       events(enc_trial_idx(e)).old       = enc_info( idx, 3);
       
       % check if the item has been presented before
       if ~isempty(prev_events)
           events(enc_trial_idx(e)).old = events(enc_trial_idx(e)).old | ismember(events(enc_trial_idx(e)).stim_idx, ...
               [prev_events.stim_idx]);
       end
       
       events(enc_trial_idx(e)).salience  = enc_info( idx, 4);
       events(enc_trial_idx(e)).resp      = nan;
       events(enc_trial_idx(e)).rt        = nan;
       
    end

    % update recognition events
    rec_out_idx = find(rec_block == ub(b));
    rec_out_trial = 1:length(rec_out_idx);
    
    rec_trial_idx = find([events.block]==ub(b) ...
        & strcmp({events.phase},'recog'));
    
    for e = 1:length(rec_trial_idx)
        
       idx = rec_out_idx(rec_out_trial == events(rec_trial_idx(e)).trial);
       
       events(rec_trial_idx(e)).stim_idx  = stim_info(rec_info( idx, 2), ...
                                                      rec_info( idx, 1));
       
       events(rec_trial_idx(e)).condition = rec_info( idx, 2);
       events(rec_trial_idx(e)).old       = rec_info( idx, 3);
       
       % check if the item has been presented before
       if ~isempty(prev_events)
           events(rec_trial_idx(e)).old = rec_info( idx, 3) | ismember(events(rec_trial_idx(e)).stim_idx,...
               [prev_events.stim_idx]);
       end
       events(rec_trial_idx(e)).salience  = rec_info( idx, 4);
       events(rec_trial_idx(e)).resp      = recog_resp( idx );
       events(rec_trial_idx(e)).rt        = recog_rt( idx );
       
    end
    
end

events = [events prev_events];

end