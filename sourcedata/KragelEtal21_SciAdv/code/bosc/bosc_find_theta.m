function theta_res = bosc_find_theta(subj, res, vis_flag, ltc_flag)
% BOSC_FIND_THETA returns information at theta peaks for a given subject
%
% INPUTS:  subj - string, specifies subject id to load data 
%
%          res - structure, bosc output structure from run_subject_bosc.m
%
%          vis_flag  - logical, flag whether to analyze contacts in DAN/VN
%
%          ltc_flag - logical, flag whether to analyze contacts in lateral temporal cortex
%
% OUTPUTS: theta_res - structure, contains information about theta oscillations

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

server_dir = ['..\..\data\eeg\' subj '\'];

if ~exist('vis_flag','var')
    vis_flag = false;
end

if ~exist('ltc_flag','var')
    ltc_flag = false;
end


n_block = size(res, 1);
blocks = 1:n_block;

n_channel = size(res, 2);

p_ep = nan(n_block, n_channel, length(res(1,1).freqs));

for b = 1:n_block
    
    fname = [server_dir subj '_sl_block' num2str(blocks(b)) '_bp.h5'];
    
    pairs = h5read(fname, '/pairs');
    fs = h5read(fname, '/srate');
    
    if ~vis_flag && ~ltc_flag
        is_hc = return_hc(subj, pairs);
        hc_idx = find(is_hc);
    elseif vis_flag
        pairs = strrep(pairs, '*', '');
        pairs = strrep(pairs, ' ', '');
        [is_dan, is_vis, pairs] = return_yeo(subj, pairs);
        
        dan_idx = find(is_dan);
        vis_idx = find(is_vis);
        
        hc_idx = [dan_idx; vis_idx];
        hc_idx = unique(hc_idx);
    elseif ltc_flag
        pairs = strrep(pairs,'*','');
        pairs = strrep(pairs,' ','');
        [is_ltc, pairs] = return_ltc(subj, pairs);
        
        hc_idx = find(is_ltc);
    end
    
    for c = 1:n_channel
        
        is_osc = res(b,c).is_osc;
        
        if isempty(is_osc)
            continue
        end
        
        ied_times = round(res(b,c).ied_out.pos(res(b,c).ied_out.chan==hc_idx(c))*fs);
        is_ied = false(1, size(is_osc,2));
        is_ied(ied_times) = true;
        
        tc = ones(1, fs);
        is_ied = conv(is_ied, tc, 'same');
        is_ied(is_ied>0)=1;
        
        p_ep(b,c,:) = nanmean(is_osc(:,~is_ied),2);
        
    end
end

for_peak = squeeze(nanmean(p_ep))';

peaks = []; freqs = []; chans = [];
for c = 1:size(for_peak,2)
   [pks, locs] = findpeaks(for_peak(:,c), ...
                           res(1,1).freqs, ...
                           'MinPeakProminence', .02);
                       
    peaks = [peaks; pks(:)];
    freqs = [freqs; locs(:)];
    chans = [chans; c*ones(size(pks(:)))];
    
end

theta_res.p_ep = p_ep;
theta_res.peaks = peaks;
theta_res.freqs = freqs;
theta_res.chans = chans;

