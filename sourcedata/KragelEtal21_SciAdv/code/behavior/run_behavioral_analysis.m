function res = run_behavioral_analysis()
% RUN_BEHAVIORAL_ANALYSIS runs behavioral analysis on all Ss and runs group stats
%
% INPUTS:  no inputs, runs group behavioral analysis
%
% OUTPUTS: res - structure, contains behavioral results
%
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

subjects =  {'S1','S2','S3','S4','S5','S6'};

for s = 1:length(subjects)
       
    % load events
    events = get_events(subjects{s});
    
    % memory stats
    res(s) = compute_recog_performance(events);
    
    % cvp
    close_only = false;
    exclude_precue = true;
    cvp_res = spotlight_cvp(events, close_only, exclude_precue);
    
    cvp_out{s} = cvp_res.cvp;
    
    fix_enc = events(strcmp({events.type}, 'fixation') & ...
        strcmp({events.phase},'encode'));
    
    sal = [fix_enc.salience]==1;
    rec = [fix_enc.resp]==1;
    
    cvp_hit(s,:) = nanmean(cvp_res.cvp(rec==1,:));
    cvp_miss(s,:) = nanmean(cvp_res.cvp(rec==0,:));
    
    first_rec = [];
    for i = 1:size(cvp_res.recalls,1)
        idx = find(cvp_res.recalls(i,:)>0,1,'first');
        if ~isempty(idx)
            first_rec(i) = cvp_res.recalls(i,idx);
        else
            first_rec(i) = nan;
        end
    end
    
    for fr = 1:5
        cvp_start(s,fr,:) = nanmean(cvp_res.cvp(first_rec==fr,:));
    end
    
    cvp(s,1,:) = nanmean(cvp_res.cvp(sal==1 & rec==1,:));
    cvp(s,2,:) = nanmean(cvp_res.cvp(sal==1 & rec==0,:));
    cvp(s,3,:) = nanmean(cvp_res.cvp(sal==0 & rec==1,:));
    cvp(s,4,:) = nanmean(cvp_res.cvp(sal==0 & rec==0,:));
    
end

% compute some things reported in the results section
d_prime_mean = nanmean([res.all_d_prime]);
d_prime_sem  = nanstd([res.all_d_prime])/sqrt(length([res.all_d_prime]));
stats_dp = mes([res.all_d_prime]', 0, 'g1');

% overall hit and cr rates
tmp = cat(3, res.all_C);

hit = squeeze(tmp(1,1,:));
hit_m = nanmean(hit);
hit_se = nanstd(hit)./sqrt(length(hit));
stats_hit = mes(hit, 0.5, 'g1');

cr = squeeze(tmp(2,2,:));
cr_m = nanmean(cr);
cr_se = nanstd(cr)./sqrt(length(cr));
stats_cr = mes(cr, 0.5, 'g1');

% plot hit, cr 
fig('FontSize', 12, ...
    'Font', 'Arial', ...
    'border', 'on');

set(gcf,'Position', [ 17.9652   13.8498    8.0000   12.0000]);
subplot(2,2,3);
ax = gca;
[he, hs] = acc_plot(ax, 100*hit, 100*cr, {'Old','New'}, [50 100]);
xlabel('Condition');
ylabel('Accuracy (%)');
uistack(hs,'top')

[hl1, ht1] = add_siglines(ax, hit, .5*ones(size(hit)), [1 1], 90);
ht1.Color = 'r';
[hl1, ht1] = add_siglines(ax, cr, .5*ones(size(hit)), [2 2], 82);
ht1.Color = 'r';

% number of fixations at retrieval by novelty

n_fix_old = nan(1, length(res));
n_fix_new = nan(1, length(res));
for s = 1:length(res)
    n_fix_old(s) = nanmean([res(s).n_fix.n_old_high ...
        res(s).n_fix.n_old_low]);
    n_fix_new(s) = nanmean([res(s).n_fix.n_new_high ...
        res(s).n_fix.n_new_low]);
end

nf_old_mean = nanmean(n_fix_old);
nf_old_sem = nanstd(n_fix_old)/sqrt(length(n_fix_old));

nf_new_mean = nanmean(n_fix_new);
nf_new_sem = nanstd(n_fix_new)/sqrt(length(n_fix_new));

stats_nfix = mes(n_fix_new', n_fix_old', 'hedgesg','isDep',1);

% number of fixations at recognition accuracy
n_fix_hit = nan(1, length(res));
n_fix_miss = nan(1, length(res));
n_fix_cr = nan(1, length(res));
n_fix_fa = nan(1, length(res));

for s = 1:length(res)
    n_fix_hit(s) = nanmean([res(s).n_fix.n_old_high(res(s).n_fix.resp_old_high==1) ...
        res(s).n_fix.n_old_low(res(s).n_fix.resp_old_low==1)]);
    n_fix_miss(s) = nanmean([res(s).n_fix.n_old_high(res(s).n_fix.resp_old_high~=1) ...
        res(s).n_fix.n_old_low(res(s).n_fix.resp_old_low~=1)]);
 
    n_fix_cr(s) = nanmean([res(s).n_fix.n_new_high(res(s).n_fix.resp_new_high~=1) ...
        res(s).n_fix.n_new_low(res(s).n_fix.resp_new_low~=1)]);
    n_fix_fa(s) = nanmean([res(s).n_fix.n_new_high(res(s).n_fix.resp_new_high==1) ...
        res(s).n_fix.n_new_low(res(s).n_fix.resp_new_low==1)]); 
end

stats_nfix_recog = mes(n_fix_hit', n_fix_miss', 'hedgesg','isDep',1);

nf_hit_mean = nanmean(n_fix_hit);
nf_hit_sem = nanstd(n_fix_hit)/sqrt(length(n_fix_hit));

nf_miss_mean = nanmean(n_fix_miss);
nf_miss_sem = nanstd(n_fix_miss)/sqrt(length(n_fix_miss));

stats_nfix_lure = mes(n_fix_cr', n_fix_fa', 'hedgesg','isDep',1);

nf_cr_mean = nanmean(n_fix_cr);
nf_cr_sem = nanstd(n_fix_cr)/sqrt(length(n_fix_cr));

nf_fa_mean = nanmean(n_fix_fa);
nf_fa_sem = nanstd(n_fix_fa)/sqrt(length(n_fix_fa));

% rt
rt_hit = nan(1, length(res));
rt_miss = nan(1, length(res));
rt_cr = nan(1, length(res));
rt_fa = nan(1, length(res));

for s = 1:length(res)
    rt_hit(s) = nanmean([res(s).rt.dist_old_high(res(s).n_fix.resp_old_high==1) ...
        res(s).rt.dist_old_low(res(s).n_fix.resp_old_low==1)]);
    rt_miss(s) = nanmean([res(s).rt.dist_old_high(res(s).n_fix.resp_old_high~=1) ...
        res(s).rt.dist_old_low(res(s).n_fix.resp_old_low~=1)]);
 
    rt_cr(s) = nanmean([res(s).rt.dist_new_high(res(s).n_fix.resp_new_high~=1) ...
        res(s).rt.dist_new_low(res(s).n_fix.resp_new_low~=1)]);
    rt_fa(s) = nanmean([res(s).rt.dist_new_high(res(s).n_fix.resp_new_high==1) ...
        res(s).rt.dist_new_low(res(s).n_fix.resp_new_low==1)]);
end

stats_rt_recog = mes(rt_hit', rt_miss', 'hedgesg','isDep',1);

rt_hit_mean = nanmean(rt_hit);
rt_hit_sem = nanstd(rt_hit)/sqrt(length(rt_hit));

rt_miss_mean = nanmean(rt_miss);
rt_miss_sem = nanstd(rt_miss)/sqrt(length(rt_miss));

stats_rt_lure = mes(rt_cr', rt_fa', 'hedgesg','isDep',1);

rt_cr_mean = nanmean(rt_cr);
rt_cr_sem = nanstd(rt_cr)/sqrt(length(rt_cr));

rt_fa_mean = nanmean(rt_fa);
rt_fa_sem = nanstd(rt_fa)/sqrt(length(rt_fa));

stats_rt_hitfa = mes(rt_hit', rt_fa', 'hedgesg','isDep',1);
stats_rt_hitcr = mes(rt_hit', rt_cr', 'hedgesg','isDep',1);

% plot rt

fig('FontSize', 12, ...
    'Font', 'Arial', ...
    'border', 'on');

set(gcf,'Position', [ 17.9652   13.8498    8.0000   12.0000]);
subplot(2,1,2);
[he, hs] = errorplot(gca, rt_hit', ...
                          rt_miss', ...
                          rt_cr', ...
                          rt_fa', ...
                          {'Hit','Miss','CR','FA'}, ...
                          [0 18], ...
                          [17 13 15 13]);

xlabel('Trial Type');
ylabel('Response Time (sec)');

% plot n_fix
subplot(2,1,1);

[he, hs] = errorplot(gca, n_fix_hit', ...
                          n_fix_miss', ...
                          n_fix_cr', ...
                          n_fix_fa', ...
                          {'Hit','Miss','CR','FA'}, ...
                          [0 40], ...
                          [39 33 36 33]);

xlabel('Trial Type');
ylabel('Fixation Count (N)');

cvp_all = nan(6,10);
for s = 1:6
    cvp_all(s,:) = nanmean(cvp_out{s});
end

fig('FontSize', 12, ...
    'Font', 'Arial', ...
    'border', 'on', ...
    'units', 'centimeters', ...
    'width', 8, ...
    'height', 12);

cvp_all(:,5) = nan;
cvp_hit(:,5) = nan;
cvp_miss(:,5) = nan;

ax = subplot(211);
[hl, hp] = plot_single_cvp(cvp_all(:,1:9), ax);
set(gca,'YLim',[0 .3], ...
        'XTick',[-4 -3 -2 -1 1 2 3 4], ...
        'XTickLabel', {'-4','-3','-2','-1','+1','+2','+3','+4'});

ax = subplot(212);cla
[hl, hp] = plot_cvp_diff(cvp_hit, cvp_miss, ax);
set(gca,'YLim',[0 1], ...
        'XLim',[0 6]);

[stats_all, summ_all] = cvp_anova(cvp_all(:,1:9));

lags = [-1, 1];
for l = 1:length(lags)
    [stats_bylag(l), summ_bylag{l}] = cvp_condition_anova(cvp, lags(l));
end

% this concludes all behavioral measures, other than revisitations, which
% are handled in a different script

% information for supplemental table

for s = 1:length(res)
    
    nfix_enc = [res(s).n_enc_fix.n_hit_high ...
        res(s).n_enc_fix.n_hit_low ...
        res(s).n_enc_fix.n_miss_high ...
        res(s).n_enc_fix.n_miss_low];
    
    enc_n(s) = nanmean(nfix_enc);
    enc_rate(s) = enc_n(s)/3;
    
    nfix_rec = [res(s).n_fix.n_old_high ...
               res(s).n_fix.n_old_low ...
               res(s).n_fix.n_new_high ...
               res(s).n_fix.n_new_low];
           
           rt_rec = [res(s).rt.dist_old_high ...
               res(s).rt.dist_old_low ...
               res(s).rt.dist_new_high ...
               res(s).rt.dist_new_low];
           
           rec_n(s) = nanmean(nfix_rec);
           rec_rate(s) = nanmean(nfix_rec./rt_rec);
    
end


end

function [he, hs] = acc_plot(ax, cond1, cond2, label, ylim)

color = [.7 .7 .7];

x = [1 2];
for s = 1:size(cond1,1)
    y = [cond1(s) cond2(s)];% cond3(s) cond4(s)];
    hold on;
    hs(s) = plot(ax, x+.1*rand(size(x))-.05, y, 'o');
    hs(s).Color = color;
    hs(s).MarkerFaceColor = color;
    hs(s).MarkerEdgeColor = color;
    hs(s).MarkerSize = 2;  
end

color = 'k';

x = [1 2];
y = [nanmean(cond1), nanmean(cond2)];
e = [nanstd(cond1)./sqrt(sum(~isnan(cond1))), ...
    nanstd(cond2)./sqrt(sum(~isnan(cond2)))];

he = errorbar(ax, x, y, e);
he.LineStyle = 'none';
he.Marker = 'o';
he.Color = color; he.LineWidth = 2;
he.MarkerFaceColor = color;
he.MarkerEdgeColor = color;
he.MarkerSize = 4;

set(gca,'XLim', [.5 2.5], ...
        'XTick', [1 2], ...
        'YLim', ylim, ...
        'box', 'off', ...
        'LineWidth', 2, ...
        'XTickLabel', label);


end


function [he, hs] = errorplot(ax, cond1, cond2, cond3, cond4, label, ylim, ...
    yscale)

color = [.7 .7 .7];
hs = add_subjects(ax, cond1, cond2, cond3, cond4, color);
color = 'k';
he = add_errorbars(ax, cond1, cond2, cond3, cond4, color, label);
set(gca, 'YLim', ylim);

[hl, ht] = add_siglines(ax, cond1, cond4, [1 4], yscale(1));
[hl, ht] = add_siglines(ax, cond1, cond2, [1 2], yscale(2));
[hl, ht] = add_siglines(ax, cond1, cond3, [1 3], yscale(3));
[hl, ht] = add_siglines(ax, cond3, cond4, [3 4], yscale(4));

end

function [hl1, ht1] = add_siglines(ax, cond1, cond2, x, yscale)

hold on;

[p, h, stats] = signrank(cond1, cond2, 'method', 'approximate');

if p < 0.001
    txt = '***';
elseif p < .005
    txt = '**';
elseif p < .05
    txt = '*';
elseif p < .1
    txt = '~';
else
    hl1 = []; ht1 = [];
end

if ~exist('hl1','var')
    yd = diff(get(gca,'YTick')); yd = yd(1);
    ym = [cond1 cond2]; ym = max(ym(:));
    hl1 = plot(ax, x, [yscale yscale], 'k-');
    hl1.LineWidth = 2;
    ht1 = text(ax, mean(x), yscale+.2, txt, 'FontSize', 14);
    ht1.HorizontalAlignment = 'center';
end

end

function he = add_errorbars(ax, cond1, cond2, cond3, cond4, color, label)

x = [1 2 3 4];
y = [nanmean(cond1), nanmean(cond2), nanmean(cond3), nanmean(cond4)];
e = [nanstd(cond1)./sqrt(sum(~isnan(cond1))), ...
    nanstd(cond2)./sqrt(sum(~isnan(cond2))), ...
    nanstd(cond3)./sqrt(sum(~isnan(cond3))), ...
    nanstd(cond4)./sqrt(sum(~isnan(cond4)))];

he = errorbar(ax, x, y, e);
he.LineStyle = 'none';
he.Marker = 'o';
he.Color = color; he.LineWidth = 2;
he.MarkerFaceColor = color;
he.MarkerEdgeColor = color;
he.MarkerSize = 4;

set(gca,'XLim', [.5 4.5], ...
        'XTick', [1 2 3 4], ...
        'YLim', [0 0.6], ...
        'box', 'off', ...
        'LineWidth', 2, ...
        'XTickLabel', label);

end

function hs = add_subjects(ax, cond1, cond2, cond3, cond4, color)

x = [1 2];
for s = 1:size(cond1,1)
    y = [cond1(s) cond2(s)];% cond3(s) cond4(s)];
    hold on;
    hs(s) = plot(ax, x, y, 'o-');
    hs(s).Color = color;
    hs(s).MarkerFaceColor = color;
    hs(s).MarkerEdgeColor = color;
    hs(s).MarkerSize = 2;  
end

x = [3 4];
for s = 1:size(cond1,1)
    y = [cond3(s) cond4(s)];% cond3(s) cond4(s)];
    hold on;
    hs(s) = plot(ax, x, y, 'o-');
    hs(s).Color = color;
    hs(s).MarkerFaceColor = color;
    hs(s).MarkerEdgeColor = color;
    hs(s).MarkerSize = 2;  
end

end
