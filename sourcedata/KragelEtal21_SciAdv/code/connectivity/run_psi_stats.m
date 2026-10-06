function run_psi_stats()
% RUN_PSI_STATS computes phase-slope index and permutation stats at the
% subject level prior to group analysis
%
% INPUTS:  no inputs  - simple wrapper to run analysis over subjects  
%
% OUTPUTS: no outputs - saves psi data to an intermediate data directory

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

subjects = {'S1', 'S2', 'S4', 'S5', 'S6'};

for s = 1:length(subjects)

    res = run_subj_psi(subjects{s});

    if ~exist('..\..\scratch_data\connectivity\', 'dir')
        mkdir('..\..\scratch_data\connectivity\') 
    end

    save(['..\..\scratch_data\connectivity\' subjects{s} ' _psi_wvis_base.mat'], 'res');    
    
end

end

function res = run_subj_psi(subj)
% INPUTS:  subj  - string, subject id

% OUTPUTS: res - results structure with psi and permutation stats at the subject level

% load time-frequency representation to compute psi
load(['..\..\scratch_data\connectivity\' subj '_data_noar_wvis.mat'], ...
    'data', 'events', 'fix_events');

% updata data to include relevant behavioral info
data = update_trialinfo(data, fix_events);

%% use fieldtrip to compute fourier for csd

% select correct combinations
chansel = find_connections(subj, data);

% config for power and csd
cfg        = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
fboilim = round([1 10] .* 8.1920) + 1;
fboi    = fboilim(1):1:fboilim(2);
foi = (fboi-1) ./ 8.1920;
cfg.foi = foi(1:2:75);
cfg.tapsmofrq = 2;
cfg.t_ftimwin      = 4 ./ cfg.foi; % 4 cycles

cfg.toi = data.time{1}(data.time{1} >= 5.25 & data.time{1} <= 6.75);
cfg.output = 'powandcsd';
cfg.pad    = 'nextpow2';
cfg.keeptrials = 'yes';
cfg.precision = 'single'; %to save memory

% run these one channel combination at a time
channelcmb = data.label(chansel);

for i = 1:size(channelcmb,1)
    
    cfg.channelcmb = channelcmb(i,:);
    
    % new data struct with just these channels, needed to construct full
    % representation with fieldtrip code that doesn't excede memory limits
    tmp_cfg = [];
    tmp_cfg.channel = cfg.channelcmb;
    data_tmp = ft_selectdata(tmp_cfg, data);
    
    cfg.trials         = data.trialinfo(:,1) == 2 & (data.trialinfo(:,3) == 1); % revisit fixations
    
    % select only the channels we need.
    
    revisit_freq       = ft_freqanalysis(cfg, data_tmp);
    revisit_freq = ft_checkdata(revisit_freq, 'cmbrepresentation', 'full');
    revisit_freq.crsspctrm = single(revisit_freq.crsspctrm);
    
    cfg.trials = data.trialinfo(:,1) == 2 & (data.trialinfo(:,3) == 2); % other fixations
    other_freq = ft_freqanalysis(cfg, data_tmp);
    other_freq = ft_checkdata(other_freq, 'cmbrepresentation', 'full');
    other_freq.crsspctrm = single(other_freq.crsspctrm);

    % compute psi and permutation stats for difference between conditions (revisit, other)
    [res(i).psi_diff, res(i).psi_1, res(i).psi_2, res(i).null_mean, res(i).null_std, res(i).psi_z1, res(i).psi_z2]  = perm_diff_psi(revisit_freq, other_freq, [2 10]);
    res(i).label = revisit_freq.label;
end


end

function chansel = find_connections(subj, fix_freq)
% helper function to identify/find relevant connections

nchan = length(fix_freq.label);
chanindx = tril(true(nchan), -1);
cmbindx1 = repmat((1:nchan)', [1 nchan]);
cmbindx2 = repmat((1:nchan),  [nchan 1]);

chansel = [cmbindx1(chanindx) cmbindx2(chanindx)];

is_hc = return_hc(subj, fix_freq.label);

% only find combinations that are between Hc and DAN, (not within Hc and
% within DAN)

between_region = false(1, size(chansel,1));
for c = 1:size(chansel,1)
    if sum(is_hc(chansel(c,:))) == 1
        between_region(c) = true;
    end
end

chansel = chansel(between_region,:);

end

