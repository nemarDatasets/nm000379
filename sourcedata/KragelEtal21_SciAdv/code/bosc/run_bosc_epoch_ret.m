function [rec_res, events] = run_bosc_epoch_ret(subj, vis_flag, ltc_flag, ...
    theta_res)
%RUN_BOSC_EPOCH_RET epochs continuous bosc data into trial, fixation, and saccade events
%
% INPUTS:  subj - string, specifies subject id to load data 
%
%          vis_flag  - logical, flag whether to analyze contacts in DAN/VN
%
%          ltc_flag - logical, flag whether to analyze contacts in lateral temporal cortex
%
%          theta_res - structure, theta results structure, containing information about
%                      oscillatory peaks for this subjects, from bosc_find_theta.m
%
% OUTPUTS: rec_res - structure, contains epoched p_episode data for each event type
%
%          events - structure, contains event information

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

if ~exist('vis_flag','var')
    vis_flag = false;
end

if ~exist('ltc_flag','var')
    ltc_flag = false;
end

base_dir = '..\..\scratch_data\';
if ~vis_flag && ~ltc_flag
    res_file = [base_dir 'bosc' filesep subj '_bosc_res_ud_ret.mat'];
elseif vis_flag
    res_file = [base_dir 'bosc' filesep subj '_bosc_res_vis_ret.mat'];
elseif ltc_flag
    res_file = [base_dir 'bosc' filesep subj '_bosc_res_ltc_ret.mat'];
end

if exist(res_file, 'file') % load
    load(res_file, 'rec_res');
    events = get_events(subj, false, base_dir); % events for later stats
else % compute and save
    % get everything we need for this subject
    events = get_events(subj, false, base_dir);
    res = get_bosc(subj, base_dir, vis_flag, ltc_flag);
    
    base_dir = '..\..\';
    fs = get_srate(subj, base_dir);
    
    rec_events = events(strcmp({events.phase},'recog'));
    
    % find theta peaks
    if ~exist('theta_res','var')
        theta_res = bosc_find_theta(subj, res, vis_flag, ltc_flag);
    end
    
    rec_res = bosc_retrieval_stats(subj, ...
        res, ...
        rec_events, ...
        fs, ...
        theta_res, ...
        true, ... % filter ied
        vis_flag, ...
        ltc_flag);
    
    save(res_file, 'rec_res', '-v7.3');

end

end

function bosc_res = bosc_retrieval_stats(subj, res, events, fs, theta_res, ...
    filter_ied, vis_flag, ltc_flag)
% BOSC_RETRIEVAL_STATS epochs p_episode output for retrieval events
%
% INPUTS:  subj - string, specifies subject id to load data 
%
%          res - structure, bosc output structure from run_subject_bosc.m
%
%          events - structure, retrieval events
%
%          fs - double, sampling rate
%
%          theta_res - structure, theta results structure, containing information about
%                      oscillatory peaks for this subjects, from bosc_find_theta.m
%         
%          filter_ied - logical, flag specifying whether to exclude time periods with IEDs
%
%          vis_flag  - logical, flag whether to analyze contacts in DAN/VN
%
%          ltc_flag - logical, flag whether to analyze contacts in lateral temporal cortex
%
% OUTPUTS: bosc_res - structure, contains epoched saccade, fixation, and trial level data

if ~exist('filter_ied', 'var')
    filter_ied = true;
end

if ~exist('vis_flag', 'var')
    vis_flag = false;
end

if ~exist('ltc_flag', 'var')
    ltc_flag = false;
end

pre = 3000 * fs/1000;
post = 750 * fs/1000;

resp_lock = true;
%
n_block = size(res,1);
n_channel = size(res,2);

blocks = unique([events.block]);

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

eeg_dir = ['..\..\data\eeg\' subj '\'];

p_ep_low = nan(length(events), n_channel, length(-pre:post));
p_ep_high = nan(length(events), n_channel, length(-pre:post));

p_ep_sac_low = nan(length([sac_events.start_time]), n_channel, length(-pre:post));
p_ep_sac_high = nan(length([sac_events.start_time]), n_channel, length(-pre:post));

p_ep_fix_low = nan(length([fix_events.start_time]), n_channel,  length(-pre:post));
p_ep_fix_high = nan(length([fix_events.start_time]), n_channel,  length(-pre:post));

[enc_onsets, rec_onsets, rec_adjust] = get_onsets(subj);

for b = 1:n_block
    block_idx = find([events.block]==blocks(b));
    fname = [eeg_dir subj '_sl_block' num2str(blocks(b)) '_bp.h5'];
    
    pairs = h5read(fname, '/pairs');
    % find hc for filter ieds
    if ~vis_flag && ~ltc_flag
        
        is_hc = return_hc(subj, pairs);
        hc_idx = find(is_hc);
    elseif vis_flag
        pairs = strrep(pairs,'*','');
        pairs = strrep(pairs,' ','');
        [is_dan, is_vis, is_pfc, pairs] = return_yeo(subj, pairs);
        
        dan_idx = find(is_dan);
        vis_idx = find(is_vis);
        
        elec_idx = [dan_idx; vis_idx];
        hc_idx = unique(elec_idx);
    elseif ltc_flag
        pairs = strrep(pairs,'*','');
        pairs = strrep(pairs,' ','');
        [is_ltc, pairs] = return_ltc(subj, pairs);
        
        hc_idx = find(is_ltc);
    end
    
    for c = 1:n_channel
        
        potential_low = find(theta_res.chans==c & theta_res.freqs < 4);
        low = theta_res.peaks(potential_low);
        low_theta_idx = ismember(res(b,c).freqs, theta_res.freqs(potential_low(low==max(low))));
        
        potential_high = find(theta_res.chans==c & theta_res.freqs > 4 & theta_res.freqs < 10);
        high = theta_res.peaks(potential_high);
        high_theta_idx = ismember(res(b,c).freqs, theta_res.freqs(potential_high(high==max(high))));
        
        is_osc = res(b,c).is_osc;
        
        ied_times = round(res(b,c).ied_out.pos(res(b,c).ied_out.chan==hc_idx(c))*fs);
        is_ied = false(1, size(is_osc,2));
        is_ied(ied_times) = true;
        
        tc = ones(1, fs);
        is_ied = conv(is_ied, tc, 'same');
        is_ied(is_ied>0)=1;
        
        block_ev = events(block_idx);
        
        for be = 1:length(block_ev)
            if strcmp(block_ev(be).phase, 'encode')
                if block_ev(be).trial <= length(enc_onsets(b,:))
                    ref_idx = enc_onsets(b, block_ev(be).trial);
                    tri_idx = ref_idx;
                end
            elseif strcmp(block_ev(be).phase, 'recog')
                if block_ev(be).trial <= length(rec_onsets(b,:))
                    
                    ref_idx = rec_onsets(b, block_ev(be).trial);
                    ref_adj = double(rec_adjust(b,block_ev(be).trial)/1000*fs);
                    tri_idx = ref_idx + ref_adj;
                    
                    if resp_lock
                        % adjust this relative to the response of the trial,
                        % not the onset of the trial
                        tri_idx = tri_idx + round(block_ev(be).rt/1000*fs);
                        
                    end
                    
                end
            end
            
            % find matching saccade/fix events
            
            sac_idx = [sac_events.block] == block_ev(be).block & ...
                [sac_events.trial] == block_ev(be).trial & ...
                strcmp({sac_events.phase}, block_ev(be).phase);
            
            if any(sac_idx)
                if resp_lock
                    sac_rel_samples = double(sac_events(sac_idx).start_time - ...
                        block_ev(be).start_time)/1000*fs - round(block_ev(be).rt/1000*fs);
                else
                    sac_rel_samples = double(sac_events(sac_idx).start_time - ...
                        block_ev(be).start_time)/1000*fs;
                end
            else
                sac_rel_samples = [];
            end
            
            fix_idx = [fix_events.block] == block_ev(be).block & ...
                [fix_events.trial] == block_ev(be).trial & ...
                strcmp({fix_events.phase}, block_ev(be).phase);
            
            if any(fix_idx)
                if resp_lock
                    % adjust relative timing for response locking
                    fix_rel_samples = double(fix_events(fix_idx).start_time - ...
                        block_ev(be).start_time)/1000*fs - ...
                        round(block_ev(be).rt/1000*fs);
                else
                    fix_rel_samples = double(fix_events(fix_idx).start_time - ...
                        block_ev(be).start_time)/1000*fs;
                end
            else
                fix_rel_samples = [];
            end
            
            if tri_idx+post < size(is_osc,2) % recording may have restarted
                % trial
                this_idx = (tri_idx - pre):(tri_idx + post);
                
                % exclude any timepoints within 1 s of IED
                if any(low_theta_idx)
                    p_ep_low(block_idx(be), c, :) = is_osc(low_theta_idx, this_idx);
                end
                
                if any(high_theta_idx)
                    p_ep_high(block_idx(be), c, :) = is_osc(high_theta_idx, this_idx);
                end
                
                % filter ied
                if filter_ied
                    ied_idx = ismember(this_idx, find(is_ied));
                    p_ep_low(block_idx(be), c, ied_idx) = nan;
                    p_ep_high(block_idx(be), c, ied_idx) = nan;
                end
                
                % saccade
                for s = 1:length(sac_rel_samples)
                    this_idx = (ref_idx + sac_rel_samples(s) - pre):(ref_idx + sac_rel_samples(s) + post);
                    if max(this_idx) < size(is_osc,2)
                        if any(low_theta_idx)
                            p_ep_sac_low(sac_events(sac_idx).idx(s), c, :) = is_osc(low_theta_idx, this_idx);
                        end
                        if any(high_theta_idx)
                            p_ep_sac_high(sac_events(sac_idx).idx(s), c, :) = is_osc(high_theta_idx, this_idx);
                        end
                        
                        % filter ied
                        if filter_ied
                            ied_idx = ismember(this_idx, find(is_ied));
                            p_ep_sac_low(sac_events(sac_idx).idx(s), c, ied_idx) = nan;
                            p_ep_sac_high(sac_events(sac_idx).idx(s), c, ied_idx) = nan;
                        end
                        
                    end
                end
                
                % fixation
                for f = 1:length(fix_rel_samples)
                    this_idx = (ref_idx + fix_rel_samples(f) - pre):(ref_idx + fix_rel_samples(f) + post);
                    if max(this_idx) < size(is_osc,2)
                        if any(low_theta_idx)
                            p_ep_fix_low(fix_events(fix_idx).idx(f), c, :) = is_osc(low_theta_idx, this_idx);
                        end
                        if any(high_theta_idx)
                            p_ep_fix_high(fix_events(fix_idx).idx(f), c, :) = is_osc(high_theta_idx, this_idx);
                        end
                        
                        % filter ied
                        if filter_ied
                            ied_idx = ismember(this_idx, find(is_ied));
                            p_ep_fix_low(fix_events(fix_idx).idx(f), c, ied_idx) = nan;
                            p_ep_fix_high(fix_events(fix_idx).idx(f), c, ied_idx) = nan;
                        end
                    end
                end
                
            end
        end
        
    end
end

bosc_res.p_ep_low = p_ep_low;
bosc_res.p_ep_sac_low = p_ep_sac_low;
bosc_res.p_ep_fix_low = p_ep_fix_low;

bosc_res.p_ep_high = p_ep_high;
bosc_res.p_ep_sac_high = p_ep_sac_high;
bosc_res.p_ep_fix_high = p_ep_fix_high;

end

function [enc_onsets, rec_onsets] = get_trial_onsets(fname, bp_flag)

if bp_flag
    pairs = deblank(h5read(fname, '/pairs'));
    idx = strcmp(pairs, 'SYNC');
else
    channels = deblank(h5read(fname, '/channels'));
    idx = strcmp(channels, 'SYNC');
end

% read sync data
sync_data = h5read(fname, '/eeg', [find(idx) 1], [1 Inf]);

% find event starts
[~,locs]=findpeaks(sync_data,'MinPeakProminence',2);
onsets=locs-1;  %take point 1 samples before peak

% both of these are pre-cues (fixation, which start 750ms before
% the image appears)
enc_idx = 3:2:49; % 51 end encoding, 52,53 start rec
rec_idx = 54:1:101; % 102 end rec

enc_idx(enc_idx>length(onsets))=[]; % if run stops early
rec_idx(rec_idx>length(onsets))=[]; % if run stops early

enc_onsets = onsets(enc_idx);
rec_onsets = onsets(rec_idx);

end
