function [stats, dat] = run_power_analysis()
% RUN_POWER_ANALYSIS compute power for revisit and other fixations and
% compares at the group level
%
% INPUTS: no inputs - this function will generate power results in 
%                      supplementary material
%
% OUTPUTS: stats - fieldtrip structure containing permutations stats
%
%          dat   - fieldtrip structure containing power in time-frequency
%                  format around fixation events

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

subjects = {'S1', 'S2','S3','S4', 'S5', 'S6'};

for s = 1:length(subjects)
    [pow_revisit, pow_other] = get_subj_pow(subjects{s});
    
    % log10
    pow_revisit.powspctrm = log10(pow_revisit.powspctrm);
    pow_other.powspctrm = log10(pow_other.powspctrm);
    
    % avg over event and chan
    revisit(1,s,:,:) = nanmean(nanmean(pow_revisit.powspctrm));
    other(1,s,:,:) = nanmean(nanmean(pow_other.powspctrm));
end

dat.dimord = 'chan_rpt_freq_time';
dat.label  = {'xx'};
dat.time   = -.75:.002:.75;
dat.powspctrm    = cat(2,revisit,other);
dat.freq = pow_revisit.freq;

c1_idx = false(1, size(dat.powspctrm,2));
c1_idx(1:6) = true;

% run permutation stats
stats = run_power_stats(dat, c1_idx, ~c1_idx);

% scale back
stats.powspctrm = nanmean(10.^dat.powspctrm(:, c1_idx, :, :),2) - nanmean(10.^dat.powspctrm(:, ~c1_idx, :, :),2);
stats.powspctrm = reshape(squeeze(stats.powspctrm), 1, size(stats.powspctrm,3), size(stats.powspctrm,4));


% and generate a figure

fig('FontSize', 14, ...
    'Font', 'Arial', ...
    'border','on', ...
    'units', 'centimeters', ...
    'width', 12, ...
    'height', 8);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.maskparameter = 'mask';
cfg.maskalpha  = .3;
cfg.zlim = [-5 5];
cfg.colormap = nawhimar;

cfg = ft_singleplotTFR(cfg, stats);

set(gca,'XTickLabel',[-500 0 500]);

xlabel('Time to Fixation (ms)');
ylabel('Frequency (Hz)');

hold on
hr = plot([0 0], get(gca,'YLim'), '--');
hr.Color = 'k';
hr.LineWidth = 2;
title('');

hc = colorbar;
ylabel(hc,'\DeltaPower (\muV^2)')

end

function [pow_revisit, pow_other] = get_subj_pow(subj)


load(['..\..\scratch_data\connectivity\' subj '_data_noar_wvis.mat'], ...
     'data', 'events', 'fix_events');

% updata data to include relevant behavioral info
data = update_trialinfo(data, fix_events);

cfg        = [];

% just hc channels
cfg.channel = data.label(startsWith(data.label, 'A') | ...
                startsWith(data.label, 'B') | ...
                startsWith(data.label, 'C') | ...
                startsWith(data.label, 'D'));

cfg.method = 'mtmconvol';
cfg.taper = 'hanning';

% frequencies to include
cfg.foi = logspace(log10(1), log10(200), 50);
cfg.tapsmofrq = 2;
cfg.t_ftimwin = 2 ./ cfg.foi; % 2 cycles

% peri-fixation window
cfg.toi = data.time{1}(data.time{1} >= 5.25 & data.time{1} <= 6.75);
cfg.output = 'pow';
cfg.pad    = 'nextpow2';
cfg.keeptrials = 'yes';
cfg.trials = data.trialinfo(:,1) == 2 & (data.trialinfo(:,3) == 1 | data.trialinfo(:,3) == 2); % fixations

% compute power
pow       = ft_freqanalysis(cfg, data);

% select revisitation fixations
cfg = [];
cfg.trials = pow.trialinfo(:,3)==1;
pow_revisit = ft_selectdata(cfg, pow);

% select other fixations
cfg = [];
cfg.trials = pow.trialinfo(:,3)==2;
pow_other = ft_selectdata(cfg, pow);


end
