function [stats, stats_tavg] = plot_group_plv_byroi()
% PLOT_GROUP_PLV_BYROI plots group synchrony results compute with
% the phase-locking value
%
% INPUTS:  no inputs, loads relevant data from temporary data directories
%
% OUTPUTS: stats - cell array of fieldtrip structures, contains permutation
%                  results for each network
%
%          stats_tavg - cell array of fieldtrip structures, contains permutation
%                  results for each network average over peri-fixation timepoints

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

subjects = {'S1','S2','S4','S5','S6'};
roi_names = {'VN', 'DAN'};

colors = [255 0 0; ... vis
          255 191 13]/255; % dan

fboilim = round([1 10] .* 8.1920) + 1;
fboi    = fboilim(1):1:fboilim(2);
foi = (fboi-1) ./ 8.1920;
foi = foi(1:2:75);

n_freq = length(foi);
n_tp = 751;

plv = nan(length(subjects), length(roi_names), n_freq, n_tp);
plv_z = plv;

for s = 1:length(subjects)
    
   load(['..\..\scratch_data\connectivity\' subjects{s} '_plv_base_wvis.mat'], 'res');
   
   p = 2*normcdf(-abs(res.plv_z));
   
   [is_hc, is_vis, is_dan, chansel] = get_plv_rois(subjects{s});
   
   roi_sel = [];
   roi_sel(1,:) = is_vis & ~is_hc;
   roi_sel(2,:) = is_dan & ~is_hc;
   
   for r = 1:length(roi_names)
       plv(s,r,:,:) = nanmean(res.plv(roi_sel(r,res.chansel(:,1))==1,:,:,:),1);
       plv_z(s,r,:,:) = nanmean(res.plv_z(roi_sel(r,res.chansel(:,1))==1,:,:,:),1);
   end
  
    
end

fig('FontSize', 14, ...
    'Font', 'Arial', ...
    'border','on', ...
    'units', 'centimeters', ...
    'width', 30, ...
    'height', 6);


for r = 1:(length(roi_names))
     
plv = reshape(plv_z(:,r,:,:), 1, size(plv_z,1), size(plv_z,3), size(plv_z,4));
plv_null = zeros(size(plv));

dat.dimord = 'chan_rpt_freq_time';
dat.label  = {'xx'};
dat.time   = -.75:.002:.75;
dat.plv    = cat(2,plv,plv_null);

fboilim = round([1 10] .* 8.1920) + 1;
fboi    = fboilim(1):1:fboilim(2);
foi = (fboi-1) ./ 8.1920;
dat.freq = foi(1:2:75);

c1_idx = false(1, size(dat.plv,2));
c1_idx(1:5) = true;


stats{r} = run_group_plv_stats(dat, c1_idx, ~c1_idx);

stats_tavg{r} = run_group_plv_stats(dat, c1_idx, ~c1_idx, false, true);

% comput effect size and simple t for text:

c1 = squeeze(nanmean(nanmean(dat.plv(:, c1_idx,stats_tavg{r}.mask,:),4),3));
c2 = squeeze(nanmean(nanmean(dat.plv(:, ~c1_idx,stats_tavg{r}.mask,:),4),3));

mes_stats(r) = mes(c1', 0, 'g1');

stats{r}.plv = reshape(nanmean(plv,2), size(stats{r}.prob));

subplot(1,length(roi_names)+1,r);

cfg = [];
cfg.parameter = 'plv';
cfg.maskparameter = 'mask';
cfg.maskstyle      = 'opacity';
cfg.maskalpha  = .2;
cfg.zlim = [-5 5];
cfg.colormap = parula;
    
ft_singleplotTFR(cfg, stats{r});

set(gca,'XTickLabel',[-500 0 500]);

xlabel('Time to Fixation (ms)');
ylabel('Frequency (Hz)');

hold on
hr = plot([0 0], get(gca,'YLim'), '--');
hr.Color = 'k';
hr.LineWidth = 2;

ht = title(roi_names{r});
ht.Color = colors(r,:);

hc = colorbar;
ylabel(hc,'Synchrony (Z)')


end

% and average over time

subplot(1,length(roi_names)+1,length(roi_names)+1);
hold on

y = nanmean(plv_z,4);
x = stats_tavg{1}.freq;
x = repmat(x, length(roi_names), 1);
y_m = squeeze(nanmean(y));
e   = squeeze(nanstd(y)./sqrt(sum(~isnan(y(:,:,1,1)))));

e_fbl = e'; e_fbl = reshape(e_fbl, size(e_fbl,1), 1, size(e_fbl,2));
[hl1, hp1] = boundedline(x', y_m', e_fbl);



for r = 1:length(hl1)
    hl1(r).Color = colors(r,:);
    hl1(r).LineWidth = 2;
    hp1(r).FaceColor = colors(r,:); %[0 255 0]/255;
    hp1(r).FaceAlpha = .4;
end


set(gca, 'XLim', [1 10], ...
         'YLim', [0 7], ...
         'LineWidth', 2);
     
hleg = legend(hl1,roi_names);
hleg.Location = 'best';
hleg.Box = 'off';

yval = [7 7.5];
for r = 1:2
l = bwlabeln(stats_tavg{r}.mask);

ul = unique(l);
for i = 1:max(ul)
    hl(i) = plot(stats_tavg{r}.freq(l==ul(i)), yval(r)*ones(size(stats_tavg{r}.freq(l==ul(i)))), '-');
    hl(i).Color = colors(r,:);
    hl(i).LineWidth = 2;
end
end

set(gca,'YLim',[0 8])

xlabel('Synchrony (Z)');
ylabel('Frequency (Hz)');


end

