function data = epoch_data(subj, events, include_vis)
% EPOCH_DATA takes raw h5 bipolar pairs and runs some basic preprocessing
% and redefining into individual trials/fixations etc
%
% INPUTS: subj - string, subject id to process data
%         
%         events - structure, contains behavioral and eye events for this subj
%
% OUTPUTS: data - fieldtrip structure containing epoched data

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

% if we want to mask points near IEDs
filter_ied = false;

% some parameters - could be changed to inputs
fs = get_srate(subj);

% includes 2 sec (or more) of buffer
pre = 3000 * fs/1000; % 3000 ms
post = 5000 * fs/1000; % 5000 ms

blocks = unique([events.block]);
n_block = length(unique(blocks));

% from all events, separate trials, saccades, and fixations

% only look at encoding events
events = events(strcmp({events.phase},'encode'));

[events, fix_events, sac_events] = separate_events(events);

%%
eeg_dir = ['..\..\data\eeg\' subj '\'];

[enc_onsets, rec_onsets, rec_adjust] = get_onsets(subj);

for b = 1:n_block
    block_idx = find([events.block]==blocks(b));
    fname = [eeg_dir subj '_sl_block' num2str(blocks(b)) '_bp.h5'];
    
    % find hc for filter ieds
    pairs = h5read(fname, '/pairs');
    
    for p = 1:length(pairs)
        pairs{p} = upper(pairs{p});
        pairs{p} = strrep(pairs{p},'*','');
    end
    
    is_hc = return_hc(subj, pairs);
    try
        [is_dan, is_vis] = return_yeo(subj, pairs);
    catch
        is_dan = false(size(pairs));
        is_vis = false(size(pairs));
    end
    
    dat = h5read(fname, '/eeg');

    if include_vis
        chan_idx = is_hc | is_dan | is_vis;
    else
        chan_idx = is_hc | is_dan;
    end
       
    trl_trial = [];
    trl_sac   = [];
    trl_fix   = [];
    
    block_ev = events(block_idx);
    
    for be = 1:length(block_ev)
        if strcmp(block_ev(be).phase, 'encode')
            if block_ev(be).trial <= length(enc_onsets(b,:))
                ref_idx = enc_onsets(b,block_ev(be).trial);
                tri_idx = ref_idx;
            end
        elseif strcmp(block_ev(be).phase, 'recog')
            if block_ev(be).trial <= length(rec_onsets(b,:))
                ref_idx = rec_onsets(b,block_ev(be).trial);
                ref_adj = double(rec_adjust(b,block_ev(be).trial)/1000*fs);
                tri_idx = ref_idx + ref_adj;
            end
        end
        
        % find matching saccade/fix events
        
        sac_idx = [sac_events.block] == block_ev(be).block & ...
            [sac_events.trial] == block_ev(be).trial & ...
            strcmp({sac_events.phase}, block_ev(be).phase);
        
        if any(sac_idx)
            sac_rel_samples = double(sac_events(sac_idx).start_time - ...
                block_ev(be).start_time)/1000*fs;
        else
            sac_rel_samples = [];
        end
        
        fix_idx = [fix_events.block] == block_ev(be).block & ...
            [fix_events.trial] == block_ev(be).trial & ...
            strcmp({fix_events.phase}, block_ev(be).phase);
        
        if any(fix_idx)
            fix_rel_samples = double(fix_events(fix_idx).start_time - ...
                block_ev(be).start_time)/1000*fs;
        else
            fix_rel_samples = [];
        end
        
        if tri_idx+post < size(dat,2) % recording may have restarted
            % trial
            this_idx = (tri_idx - pre):(tri_idx + post);
            
            % exclude any timepoints within 1 s of IED
            
            trl_trial(block_idx(be),:) = [tri_idx-pre tri_idx+post pre];
            
            % filter ied
            if filter_ied
                ied_idx = ismember(this_idx, find(is_ied));
                epoched_trials(block_idx(be), :, ied_idx) = nan;
            end
            
            % saccade
            for s = 1:length(sac_rel_samples)
                this_idx = (ref_idx + sac_rel_samples(s) - pre):(ref_idx + sac_rel_samples(s) + post);
                
                %                     epoched_saccades(sac_events(sac_idx).idx(s), :, :) = dat(c, this_idx);
                trl_sac(sac_events(sac_idx).idx(s),:) = [ref_idx + sac_rel_samples(s) - pre ...
                    ref_idx + sac_rel_samples(s) + post ...
                    pre];
                
                % filter ied
                if filter_ied
                    ied_idx = ismember(this_idx, find(is_ied));
                    epoched_saccades(sac_events(sac_idx).idx(s), :, ied_idx) = nan;
                end
            end
            
            % fixation
            for f = 1:length(fix_rel_samples)
                this_idx = (ref_idx + fix_rel_samples(f) - pre):(ref_idx + fix_rel_samples(f) + post);
                
                %                     epoched_fixations(fix_events(fix_idx).idx(f), :, :) = dat(c, this_idx);
                trl_fix(fix_events(fix_idx).idx(f),:) = [ref_idx + fix_rel_samples(f) - pre ...
                    ref_idx + fix_rel_samples(f) + post ...
                    pre];
                
                % filter ied
                if filter_ied
                    ied_idx = ismember(this_idx, find(is_ied));
                    epoched_fixations(fix_events(fix_idx).idx(f), :, ied_idx) = nan;
                end
            end
            
        end
    end
    
    % init continuous data
    data{b}.label = pairs(chan_idx); % convert pairs to lower and remove all *
    data{b}.fsample = fs;
    data{b}.trial = {dat(chan_idx,:)};
    data{b}.time = {linspace(0,(size(dat,2)-1)/fs,size(dat,2))};
    
    % resample to 500 Hz
    cfg = [];
    cfg.resamplefs = 500;
    data{b} = ft_resampledata(cfg, data{b});
    
    % remove line noise and harmonics
    cfg = [];
    cfg.dftfreq       = [60 120 180];
    data{b} = ft_preprocessing(cfg, data{b});  
    
    % and split into trials
    
    trl_trial(all(trl_trial==0,2),:) = [];
    trl_sac(all(trl_sac==0,2),:) = [];
    trl_fix(all(trl_fix==0,2),:) = [];

    cfg = [];
    cfg.trl = cat(1, [trl_trial 1*ones(size(trl_trial,1),1)], ...
        [trl_fix 2*ones(size(trl_fix,1),1)]);%, ...
    
    if isempty(cfg.trl) % no trial data for this block
        data{b} = [];
        continue
    end
    
    cfg.trl(:,1:3) = round(cfg.trl(:,1:3) * (1/fs) * ...
        500);
    
    
    data{b} = ft_redefinetrial(cfg, data{b});
    
    trial_info(b).trl = cfg.trl;  

end

% and append data across all blocks
cfg = [];
cfg.keepsampleinfo = 'no';
if strcmp(subj, 'S6')
    data = ft_appenddata(cfg, data{[1 2 3 5]});
else
    data = ft_appenddata(cfg, data{:});
end
% ec125 remove block 4
if strcmp(subj, 'S6')
    events([events.block]==4) = [];
    fix_events([fix_events.block]==4) = [];
end

if ~exist('..\..\scratch_data\connectivity\', 'dir')
    mkdir('..\..\scratch_data\connectivity\');
end

if include_vis
    save(['..\..\scratch_data\connectivity\' subj '_data_noar_wvis.mat'], ...
          'data', ...
          'events', ...
          'fix_events', ...
          'trial_info', ...
          '-v7.3');
else
     save(['..\..\scratch_data\connectivity\' subj '_data_noar.mat'], ...
           'data', ...
           'events', ...
           'fix_events', ...
           'trial_info', ...
           '-v7.3');   
end

end

function [events, fix_events, sac_events] = separate_events(events)

% remove saccade and fixation events to unique structures
sac_idx = strcmp({events.type},'saccade');
fix_idx = strcmp({events.type},'fixation');

sac_events = events(sac_idx);
fix_events = events(fix_idx);

% update events to include saccade and fixation indices
sac_cnt = 0;
for se = 1:length(sac_events)
    if ~isempty(sac_events(se).start_time)
        sac_events(se).idx = (1:length(sac_events(se).start_time)) + sac_cnt;
        sac_cnt = sac_cnt + length(sac_events(se).start_time);
    end
end

fix_cnt = 0;
for fe = 1:length(fix_events)
    if ~isempty(fix_events(fe).start_time)
        fix_events(fe).idx = (1:length(fix_events(fe).start_time)) + fix_cnt;
        fix_cnt = fix_cnt + length(fix_events(fe).start_time);
    end
end

events(sac_idx | fix_idx) = [];
events(strcmp({events.type},'blink')) = [];

end
