function plot_psi_stats()
% PLOT_PSI_STATS plots phase slope index group results
%
% INPUTS:  no inputs, loads relevant data from temporary data directories
%
% OUTPUTS: no outputs, runs stats and generates figures

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

for s = 1:length(subjects)
    res(s) = gather_psi_stats(subjects{s});
end

% init some arrays
all_conn = nan(5,751);
all_dan = all_conn;
all_vis = all_conn;

rev_conn = nan(5,751);
rev_dan = rev_conn;
rev_vis = rev_conn;

other_conn = nan(5,751);
other_dan = other_conn;
other_vis = other_conn;

t_idx = 76:676; % time window before fixation when hc -> vn/dan is believed to occur

for s = 1:length(subjects)
    
    % collect rois
    [rois, is_hc, is_vis, is_dan, chansel, pairs] = get_plv_rois(subjects{s});
    
    dan_mask = is_dan(chansel(:,1))' & ~is_hc(chansel(:,1))';
    vis_mask = is_vis(chansel(:,1))' & ~is_hc(chansel(:,1))';
    
    [all_conn(s,:), rev_conn(s,:), other_conn(s,:)] = gen_null(res(s).psi_revisit, ...
        res(s).psi_other, ...
        res(s).psi_diff, ...
        t_idx);
    
    [all_vis(s,:), rev_vis(s,:), other_vis(s,:)] = gen_null(res(s).psi_revisit(vis_mask,:), ...
        res(s).psi_other(vis_mask,:), ...
        res(s).psi_diff(vis_mask,:), ...
        t_idx);
    
    [all_dan(s,:), rev_dan(s,:), other_dan(s,:)] = gen_null(res(s).psi_revisit(dan_mask,:), ...
        res(s).psi_other(dan_mask,:), ...
        res(s).psi_diff(dan_mask,:), ...
        t_idx);

end

% stats and effect sizes for reporting using mes toolbox

% revisit vs. other
stats_vis = group_perm_stats(rev_vis, other_vis, t_idx);
stats_dan = group_perm_stats(rev_dan, other_dan, t_idx);

idx = get_peak(stats_vis.stat(stats_vis.mask));
m_idx = find(stats_vis.mask);

c1 = nanmean(rev_vis(:, t_idx(m_idx(idx))),2);
c2 = nanmean(other_vis(:, t_idx(m_idx(idx))),2);
mes_stats_dv = mes(c1, c2, 'hedgesg', 'isdep', true); % dv, difference visual

% stats for other vs. psi of zero
stats_other_vis = group_perm_stats(other_vis, zeros(size(all_vis)), t_idx);
stats_other_dan = group_perm_stats(other_dan, zeros(size(all_vis)), t_idx);

idx = get_peak(-stats_other_vis.stat(stats_other_vis.mask));
m_idx = find(stats_other_vis.mask);
c1 = nanmean(other_vis(:, t_idx(m_idx(idx))),2);

idx = get_peak(stats_other_dan.stat(stats_other_dan.mask));
m_idx = find(stats_other_dan.mask);
c2 = nanmean(other_dan(:, t_idx(m_idx(idx))),2);
mes_stats_ov = mes(c1, 0, 'g1');
mes_stats_od = mes(c2, 0, 'g1');

% stats for revisit vs. psi of zero
stats_rev_vis = group_perm_stats(rev_vis, zeros(size(all_vis)), t_idx);
stats_rev_dan = group_perm_stats(rev_dan, zeros(size(all_vis)), t_idx);

idx = get_peak(stats_rev_dan.stat(stats_rev_dan.mask));
m_idx = find(stats_rev_dan.mask);
c2 = nanmean(rev_dan(:, t_idx(m_idx(idx))),2);

mes_stats_rv = mes(c1, 0, 'g1');
mes_stats_rd = mes(c2, 0, 'g1');

% and generate the figure

colors = [255 0 0; ... %vis
          255 191 13]/255; % dan

rev_color = [0 1 0];
other_color = [.3 .3 .3];
      
figure;

% vis
subplot(1,3,1);
hold on;

labels = {'Revisit'};
ylim = [-.5 .6];
[hlr_vis, hp_vis] = plot_single_region(rev_vis(:,t_idx), ...
    stats_rev_vis, ...
    ylim, ...
    rev_color, ...
    labels);
set(gca,'XLim',[-600 50]);

labels = {'Other'};
ylim = [-.6 .6];
[hlo_vis, hp_vis] = plot_single_region(other_vis(:,t_idx), ...
    stats_other_vis, ...
    ylim, ...
    other_color, ...
    labels);
set(gca,'XLim',[-600 50]);

hlo_vis.LineStyle = '--';

hleg = legend([hlr_vis, hlo_vis], {'Revisit','Other'});
hleg.Location = 'best';

ylabel('PSI (Z)');

% dan
subplot(1,3,2);

labels = {'Revisit'};
ylim = [-.5 .6];
[hlr_dan, hp_dan] = plot_single_region(rev_dan(:,t_idx), ...
    stats_rev_dan, ...
    ylim, ...
    rev_color, ...
    labels);
set(gca,'XLim',[-600 50]);

labels = {'Other'};
ylim = [-.6 .6];
[hlo_dan, hp_dan] = plot_single_region(other_dan(:,t_idx), ...
    stats_other_dan, ...
    ylim, ...
    other_color, ...
    labels);
set(gca,'XLim',[-600 50]);

hlo_dan.LineStyle = '--';
hleg = legend([hlr_dan, hlo_dan], {'Revisit','Other'});
hleg.Location = 'southwest';

ylabel('PSI (Z)');

subplot(1,3,3);
hold on;

subplot(1,3,3);

labels = {'VN'};
ylim = [-.5 .6];
[hl_vis, hp_vis] = plot_single_region(rev_vis(:,t_idx)-other_vis(:,t_idx), ...
    stats_vis, ...
    ylim, ...
    colors(1,:), ...
    labels);
set(gca,'XLim',[-600 50]);

hleg = legend(hl_vis, 'Visual');
hleg.Location = 'best';

labels = {'DAN'};
ylim = [-.6 .6];
[hl_dan, hp_dan] = plot_single_region(rev_dan(:,t_idx)-other_dan(:,t_idx), ...
    stats_dan, ...
    ylim, ...
    colors(2,:), ...
    labels);
set(gca,'XLim',[-600 50]);

hleg = legend([hl_dan, hl_vis], {'DAN', 'VN'});
hleg.Location = 'best';

set(gcf,'Position',[ 1175 588 737 267]);


end

function [hl1, hp1] = plot_single_region(dat, stats, ylim, colors, labels)

x = stats.time*1000;

y1 = nanmean(dat);
e1 =  nanstd(dat)./sqrt(sum(~isnan(dat(:,1))));

[hl1, hp1] = boundedline(x', y1', e1');
hold on;
hl1.Color = colors(1,:);
hl1.LineWidth = 2;
hp1.FaceColor = colors(1,:);
hp1.FaceAlpha = .2;

set(gca, 'XLim', [min(x) max(x)], ...
    'YLim', ylim, ...
    'LineWidth', 2, ...
    'XTick', -500:250:500);

xlabel('Time (ms)');
ylabel('\DeltaPSI (Z)');

hl3 = plot([0 0], ylim, '--');
hl3.Color = 'k'; hl3.LineWidth = 2;
uistack(hl3,'bottom')

if any(stats.mask) % fdr corrected
    
    [l, n] = bwlabeln(stats.mask);
    
    for i = 1:n
        xsig = x(l==i);
        
        yvals = y1+e1;
        yvals = yvals(l==i);
        yvals = max(yvals)+ max(yvals)/10;
        
        hsig = plot(xsig, min(ylim)*ones(size(xsig)), 'r-'); %yvals*ones(size(xsig))+.025
        hsig.LineWidth = 4;
    end
    
end

if any(stats.prob < .05 & ~stats.mask) % uncorrected
    
    [l, n] = bwlabeln(stats.prob < .05 & ~stats.mask);
    
    for i = 1:n
        xsig = x(l==i);
        
        yvals = y1+e1;
        yvals = yvals(l==i);
        yvals = max(yvals)+ max(yvals)/10;
        
        hsig = plot(xsig,  min(ylim)*ones(size(xsig)), 'y-'); %yvals*ones(size(xsig))+.025
        hsig.LineWidth = 4;
    end
    
end

hr = refline(0,0);
hr.Color = 'k';
uistack(hr, 'bottom');

hleg = legend(hl1, labels);
hleg.Box = 'off'; hleg.Location = 'southwest';

end

function stats = group_perm_stats(c1, c2, tidx)

if ~exist('tidx','var')
    tidx = 1:size(c1,2);
%     tidx = 1:751; % -.75 to .10 sec
end

dat.dimord = 'chan_rpt_time'; % chan is really subject but doesnt matt
dat.label  = {'xx'};
dat.time   = -.75:.002:.75;
dat.time   = dat.time(tidx);
dat.avg    = cat(1,c1,c2);
dat.avg    = reshape(dat.avg, 1, size(dat.avg,1), size(dat.avg,2));
dat.avg    = dat.avg(:,:,tidx);

c1_idx = false(1, size(dat.avg,2));
c1_idx(1:size(c1,1)) = true;

stats = run_group_psi_stats(dat, c1_idx, ~c1_idx);

end

function res = gather_psi_stats(subj)

[psi_revisit, psi_other, psi_diff, anat, label] = get_psi_by_network(subj);

res.label = label;

% some book-keeping
res.psi_diff = psi_diff;
res.psi_revisit = psi_revisit;
res.psi_other   = psi_other;

res.dan_con = anat.dan_con;
res.vis_con = anat.vis_con;

res.rev_vis = nanmean(psi_revisit(anat.vis_con,:),1);
res.oth_vis = nanmean(psi_other(anat.vis_con,:),1);

res.rev_dan = nanmean(psi_revisit(anat.dan_con,:),1);
res.oth_dan = nanmean(psi_other(anat.dan_con,:),1);

end

function [z_psi_diff, rev_sig, oth_sig] = gen_null(psi_rev, psi_oth, psi_diff, ...
    t_idx)
% this function adjusts the psi_diff value for selection based on showing a
% positive (or negative) PSI value over a given time period, and generates
% a new average (across channels meeting the criteria) measure of psi_diff

p_rev = 2*normcdf(-abs(psi_rev));
p_oth = 2*normcdf(-abs(psi_oth));

% old method, based on psi_diff
null_diff = nan(1000, size(psi_rev,2));
for n = 1:1000
    

    tmp_rev = psi_rev; % random sign flip
    tmp_oth = psi_oth;
    to_swap = round(rand(1,size(psi_rev,1))); % to exchange
    tmp_rev(to_swap==1, :) = psi_oth(to_swap==1, :);
    tmp_oth(to_swap==1, :) = psi_rev(to_swap==1, :);

    p_rev_tmp = p_rev;
    p_oth_tmp = p_oth;
    p_rev_tmp(to_swap==1, :) = p_oth(to_swap==1, :);
    p_oth_tmp(to_swap==1, :) = p_rev(to_swap==1, :);

	sig_idx = any(p_rev_tmp(:,t_idx) < .05 | p_oth_tmp(:,t_idx) < .05, 2);

    null_diff(n,:) = nanmean(tmp_rev(sig_idx==1,:) - tmp_oth(sig_idx==1,:));

end

sig_idx = any(p_rev(:,t_idx) < .05 | p_oth(:,t_idx) < .05, 2);
rev_sig_idx = any(p_rev(:,t_idx) < .05, 2);
oth_sig_idx = any(p_oth(:,t_idx) < .05, 2);
 
m_psi_diff = nanmean(psi_rev(sig_idx==1,:)-psi_oth(sig_idx==1,:),1);
z_psi_diff = (m_psi_diff - nanmean(null_diff))./nanstd(null_diff);

obs = nanmean(m_psi_diff(t_idx));
null = nanmean(null_diff(:,t_idx),2);

rev_sig = nanmean(psi_rev(sig_idx==1,:),1);
oth_sig = nanmean(psi_oth(sig_idx==1,:),1);

end

function [hl1, hp1, hl2, hp2] = plot_res(dat1, dat2, stats, ylim, colors, labels)
% subfunction to help with plotting

x = stats.time*1000;

y1 = nanmean(dat1);
e1 =  nanstd(dat1)./sqrt(size(dat1,1));

y2 = nanmean(dat2);
e2 = nanstd(dat2)./sqrt(size(dat2,1));

[hl1, hp1] = boundedline(x', y1', e1');

hl1.Color = colors(1,:);
hl1.LineWidth = 2;
hp1.FaceColor = colors(1,:);
hp1.FaceAlpha = .2;

hold on;
[hl2, hp2] = boundedline(x', y2', e2');

hl2.Color = colors(2,:);
hl2.LineWidth = 2;
hp2.FaceColor = colors(2,:);
hp2.FaceAlpha = .2;

set(gca, 'XLim', [min(x) max(x)], ...
    'YLim', ylim, ...
    'LineWidth', 2, ...
    'XTick', -500:250:500);

xlabel('Time (ms)');
ylabel('Phase Slope Index');

hl3 = plot([0 0], ylim, '--');
hl3.Color = 'k'; hl3.LineWidth = 2;
uistack(hl3,'bottom')

if any(stats.mask) % fdr corrected
    
    [l, n] = bwlabeln(stats.mask);
    
    for i = 1:n
        xsig = x(l==i);
        
        yvals = max(cat(1, y1+e1, y2+e2));
        yvals = yvals(l==i);
        yvals = max(yvals)+ max(yvals)/10;
        
        hsig = plot(xsig, yvals*ones(size(xsig))+.025, 'r-');
        hsig.LineWidth = 4;
    end
    
end

hleg = legend([hl1 hl2], labels);
hleg.Box = 'off'; hleg.Location = 'southeast';

end

function idx = get_peak(vec)

p = prctile(vec, 100);
idx = vec==p;

end
