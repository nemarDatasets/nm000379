function stats = run_revisitation_controls()
% RUN_REVISITATION_CONTROLS evaluates whether revisitations are correlated
% with distance, salience, or temporal lag 
%
% INPUTS:  no inputs, runs group behavioral analysis
%
% OUTPUTS: stats - structure, contains statistics for effects of salience
%                  and lag on revisitation behavior
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

lags_all = []; dist_all = []; sal_all = []; n_replayed_all = []; n_repeated_all = [];
for s = 1:length(subjects)
       
    % load events
    events = get_events(subjects{s});
    
    
    [res, lags, sal, dist] = compute_dist_by_lag(events);
    dist = dist/40; % pix2deg
    
    % for reviewer question, how often are revisited locations also
    % replayed at encoding?
    [n_repeated, n_replayed] = compute_enc_replayed(events);
    
    [lag_p(s), sal_p(s)] =  compute_revisit_ranks(events);
    
    p_replayed_encoding(s) = sum(n_replayed)/sum(n_repeated);
    
    n_replayed_all = [n_replayed_all n_replayed];
    n_repeated_all = [n_repeated_all n_repeated];
    
    % dist confound
    dist_neg_far(s) = nanmean(dist(lags<-1));
    dist_neg_1(s) = nanmean(dist(lags==-1));
    
    dist_pos_far(s) = nanmean(dist(lags>1));
    dist_pos_1(s) = nanmean(dist(lags==1));
    
    dist_abs_far(s) = nanmean(dist(abs(lags)>1));
    dist_abs_1(s) = nanmean(dist(abs(lags)==1));
    
    % salience confound
    sal_neg_far(s) = nanmean(sal(lags<-1));
    sal_neg_1(s) = nanmean(sal(lags==-1));
    
    sal_pos_far(s) = nanmean(sal(lags>1));
    sal_pos_1(s) = nanmean(sal(lags==1));
    
    sal_abs_far(s) = nanmean(sal(abs(lags)>1));
    sal_abs_1(s) = nanmean(sal(abs(lags)==1));
    
    % within subject comparison
    lags_all = [lags_all lags];
    dist_all = [dist_all dist];
    sal_all = [sal_all sal];
    
    % distance
    % backward
    [h, p, c, st] = ttest2(dist(lags == -1), dist(lags < -1));
    
    dist_t_back(s) = st.tstat;
    dist_p_back(s) = p;
    dist_df_back(s) = st.df;
    
    % forward
    [h, p, c, st] = ttest2(dist(lags == 1), dist(lags > 1));
    
    dist_t_forw(s) = st.tstat;
    dist_p_forw(s) = p;
    dist_df_forw(s) = st.df;
    
    % all
    [h, p, c, st] = ttest2(dist(abs(lags) == 1), dist(abs(lags) > 1));
    
    dist_t_all(s) = st.tstat;
    dist_p_all(s) = p;
    dist_df_all(s) = st.df;
    
    % salience
    
    % backward
    [h, p, c, st] = ttest2(sal(lags == -1), sal(lags < -1));
    
    sal_t_back(s) = st.tstat;
    sal_p_back(s) = p;
    sal_df_back(s) = st.df;
    
    % forward
    [h, p, c, st] = ttest2(sal(lags == 1), sal(lags > 1));
    
    sal_t_forw(s) = st.tstat;
    sal_p_forw(s) = p;
    sal_df_forw(s) = st.df;
    
    % all
    [h, p, c, st] = ttest2(sal(abs(lags) == 1), sal(abs(lags) > 1));
    
    sal_t_all(s) = st.tstat;
    sal_p_all(s) = p;
    sal_df_all(s) = st.df;
    
    % gather data for between subject comparison
    dbl(s,:) = res.dbl_y/40;
    sbl(s,:) = res.sbl_y;
    n(s,:) = res.n;
    
end

% test whether revisitations are influenced by distance and/or salience

stats.revisit_sal = mes(sal_p', .5, 'g1');
stats.revisit_lag = mes(lag_p', .5, 'g1');

% figure here
figure;
set(gcf,'Position', [241   602   286   305]);

x = 1;
y = nanmean(sal_p);
e = nanstd(sal_p)/sqrt(length(sal_p));
y_all = sal_p;

he = errorbar(x, y, e);
he.Color = 'k';
he.LineWidth = 2;

hold on
hp = plot(x' + .05*randn(size(y_all')), y_all', 'o');
for p = 1:length(hp)
    hp(p).MarkerFaceColor = [.7 .7 .7];
    hp(p).MarkerEdgeColor = [.7 .7 .7];
    hp(p).Color = [.7 .7 .7];
    hp(p).MarkerSize = 2;
end

uistack(hp, 'bottom')

x = 2;
y = nanmean(lag_p);
e = nanstd(lag_p)/sqrt(length(lag_p));
y_all = lag_p;

he = errorbar(x, y, e);
he.Color = 'k';
he.LineWidth = 2;

hold on
hp = plot(x' + .05*randn(size(y_all')), y_all', 'o');
for p = 1:length(hp)
    hp(p).MarkerFaceColor = [.7 .7 .7];
    hp(p).MarkerEdgeColor = [.7 .7 .7];
    hp(p).Color = [.7 .7 .7];
    hp(p).MarkerSize = 2;
end

uistack(hp, 'bottom')

% siglines

% sal
[h, p, c, st] = ttest(sal_p, .5);

if p < .001
    txt = '***';
elseif p < .01
    txt = '**';
elseif p < 0.05
    txt = '*';
elseif p < .1
    txt = '~';
else
    hl = []; ht = [];
end

if p < .1
yd = diff(get(gca,'YTick')); yd = yd(1);
ym = sal_p; ym = max(ym(:));

ht = text(gca, 1, ym+.3*yd, txt, 'FontSize', 14);
ht.FontWeight = 'bold';
ht.HorizontalAlignment = 'center';
end


% lag
[h, p, c, st] = ttest(lag_p, .5);

if p < .001
    txt = '***';
elseif p < .01
    txt = '**';
elseif p < 0.05
    txt = '*';
elseif p < .1
    txt = '~';
else
    hl = []; ht = [];
end

if p < .1
yd = diff(get(gca,'YTick')); yd = yd(1);
ym = lag_p; ym = max(ym(:));

ht = text(gca, 2, ym+.3*yd, txt, 'FontSize', 14);
ht.FontWeight = 'bold';
ht.HorizontalAlignment = 'center';
end


hr = plot(get(gca,'XLim'), [.5 .5], 'k--');
uistack(he,'top');
he.Color = 'k';
he.LineWidth = 2;

set(gca, 'YLim', [0.3 .7], ...
         'Xlim', [0.5 2.5], ...
         'box', 'off', ...
         'linewidth', 2, ...
         'XTick', [1 2], ...
         'XTickLabel', {'Salience','Lag'}, ...
         'FontSize', 12, ...
         'FontName', 'Arial');

xlabel('Factor');
ylabel('Percentile Rank');


% percent of fixations to distinct locations less than 2 dva apart
numer = sum(dist_all  < 2 & lags_all ~=0);
denom = sum(lags_all ~=0 & ~isnan(lags_all));

% numer/denom = 0.1831

% percent of fixations to distinct locations less than 1.25 dva apart
numer = sum(dist_all  < 1.25 & lags_all ~=0);
denom = sum(lags_all ~=0 & ~isnan(lags_all));

% numer/denom = 0.018

% percent of revisitation that lead to replay at encoding
nansum(n_replayed_all)/nansum(n_repeated_all);

% 0.0811

stats.dist_neg=mes(dist_neg_far', dist_neg_1', 'hedgesg', 'isDep', true);
stats.dist_pos=mes(dist_pos_far', dist_pos_1', 'hedgesg', 'isDep', true);
stats.dist=mes(dist_abs_far', dist_abs_1', 'hedgesg', 'isDep', true);

stats.sal_neg=mes(sal_neg_far', sal_neg_1', 'hedgesg', 'isDep', true);
stats.sal_pos=mes(sal_pos_far', sal_pos_1', 'hedgesg', 'isDep', true);
stats.sal=mes(sal_abs_far', sal_abs_1', 'hedgesg', 'isDep', true);

y = nanmean(dbl);
e = nanstd(dbl)./sqrt(sum(~isnan(dbl)));

figure;
set(gcf,'Position', [241   602   286   305]);

% subplot(3,1,1);
he = errorbar(-4:4, y, e);

set(gca, 'YLim', [0 6], ...
         'Xlim', [-5 5], ...
         'box', 'off', ...
         'linewidth', 2, ...
         'XTick', -4:4, ...
         'FontSize', 12, ...
         'FontName', 'Arial');

hold on;

hr = plot(get(gca,'XLim'), [2 2], 'k--');
uistack(he,'top');
he.Color = 'k';
he.LineWidth = 2;

xlabel('Lag');
ylabel('Distance (DVA)');

figure;
set(gcf,'Position', [241   602   286   305]);

y = nanmean(10.^sbl);
e = nanstd(10.^sbl)./sqrt(sum(~isnan(10.^sbl)));

he = errorbar(-4:4, y, e);

set(gca, 'Xlim', [-5 5], ...
         'Ylim', [1e-13 1e-10], ...
         'box', 'off', ...
         'linewidth', 2, ...
         'XTick', -4:4, ...
         'FontSize', 12, ...
         'FontName', 'Arial', ...
         'yscale', 'log');

hold on;
he.Color = 'k';
he.LineWidth = 2;

xlabel('Lag');
ylabel('Salience (Probability of fixation)');
