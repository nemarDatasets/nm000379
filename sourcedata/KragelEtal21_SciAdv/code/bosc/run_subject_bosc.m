function res = run_subject_bosc(subj)
% RUN_SUBJECT_BOSC used to detect theta oscillations for this subject
%
% INPUTS:  subj - string, specifies subject id to load data 
%
% OUTPUTS: res - structure, contains oscillation information for this subject
%                which is also written to file

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

% params for outputs
out_dir = '..\..\scratch_data\bosc\';

if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

% hc
ltc_flag = false; vis_flag = false;
res = run_bosc(subj, ltc_flag, vis_flag);
res_file = [out_dir subj '_bosc.mat'];
save(res_file, 'res', '-v7.3');

% ltc
ltc_flag = true; vis_flag = false;
run_bosc(subj, ltc_flag, vis_flag);
res_file = [out_dir subj '_ltc_bosc.mat'];
save(res_file, 'res', '-v7.3');

% vis/dan
ltc_flag = false; vis_flag = true;
run_bosc(subj, ltc_flag, vis_flag);
res_file = [out_dir subj '_vis_bosc.mat'];
save(res_file, 'res', '-v7.3');

end

function res = run_bosc(subj, ltc_flag, vis_flag)
%
% INPUTS:  subj - string, specifies subject id to load data 
%
% OUTPUTS: res - structure, contains oscillation information for this subject
%                which is also written to file

res_dir = ['..\..\data\eeg\' subj filesep];
files = dir([res_dir filesep '*_bp.h5']);

if strcmp(subj ,'S6')
    % exclude block4 as it isn't analyzed
    for f = 1:length(files)
        if contains(files(f).name,'block4')
            r_idx(f) = true;
        end
    end
    files = files(~r_idx); 
end

for f = 1:length(files) % each file is a block of data
    
    % load the raw data for this subject/block
    pairs = h5read([res_dir filesep files(f).name], '/pairs');
    eeg = h5read([res_dir filesep files(f).name], '/eeg');
    fs = h5read([res_dir filesep files(f).name], '/srate');
    
    % find ied from all channels
    [out, discharges]=spike_detector_hilbert_v23(eeg',fs);

    if ~vis_flag && ~ltc_flag
        % find hc contacts in this block
        is_hc = return_hc(subj, pairs);
        elec_idx = find(is_hc);
        out_idx = 1:length(elec_idx);
    elseif ltc_flag
        pairs = strrep(pairs,'*','');
        pairs = strrep(pairs,' ','');
        [is_ltc, pairs] = return_ltc(subj, pairs);
        
        ltc_idx = find(is_ltc);    
        elec_idx = unique(ltc_idx);
        out_idx = 1:length(elec_idx);
    elseif vis_flag
        pairs = strrep(pairs,'*','');
        pairs = strrep(pairs,' ','');
        [is_dan, is_vis, pairs] = return_yeo(subj, pairs);
        
        dan_idx = find(is_dan);
        vis_idx = find(is_vis);
        
        elec_idx = [dan_idx; vis_idx];
        elec_idx = unique(elec_idx);
        
        out_idx = 1:length(elec_idx);
    end   
    
    freqs = logspace(log10(1), log10(40), 50);
    
    for p = elec_idx'
    
        ied_times = round(out.pos(out.chan==p)*fs);
        is_ied = false(size(eeg(p,:))); 
        is_ied(ied_times) = true;
        
        tc = ones(1, fs);
        is_ied = conv(is_ied, tc, 'same');
        is_ied(is_ied>0)=1;
        % convolve with square waveform to exclude timepoints
        % within 1 sec of IED
        
        [B, T, F] = BOSC_tf(eeg(p,:), ... % signal
            freqs, ... %freqs
            fs, ... %sameplrate
            6); % number of wavelet cycles
        
        % padding for edge effects
        edge=ceil(6*fs/min(F));

        exclude_mask = is_ied;
        exclude_mask(1:edge) = true;
        exclude_mask(end-edge:end) = true;
        
        ps = nanmean(B(:, ~exclude_mask),2);
        ps = log10(ps');
        [ap_params, ap_ps] = robust_ap_fit(freqs, ps);
        ap_ps = 10.^ap_ps;
        
        [powthresh,durthresh] = BOSC_thresholds(fs, ...
            0.95, ...
            3, ...
            F, ...
            ap_ps);
        
        is_osc = nan(size(B));
        for fr = 1:size(B,1)
            is_osc(fr,:) = BOSC_detect(B(fr,:), ...
                powthresh(fr), ...
                durthresh(fr), ...
                fs);
        end
        
        res(f, out_idx(elec_idx == p)).is_osc = is_osc;
        res(f, out_idx(elec_idx == p)).freqs  = freqs;
        res(f, out_idx(elec_idx == p)).powthresh  = powthresh;
        res(f, out_idx(elec_idx == p)).durthresh  = durthresh;
        res(f, out_idx(elec_idx == p)).ap_ps  = ap_ps;
        res(f, out_idx(elec_idx == p)).ps  = ps;
        res(f, out_idx(elec_idx == p)).ied_out = out;
    end
    
end

end

