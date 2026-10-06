function stats = run_revisitation_memory()
% RUN_REVISITATION_MEMORY computes effects of revisitation on subsequent
% fixation-sequence replay and recognition performance
%
% INPUTS:  no inputs, loads event structures to run analyses
%
% OUTPUTS: stats - structure, contains stats for subsequent recognition,
%                  as well as forward and reverse replay  

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

subjects = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6'};


for s = 1:length(subjects)
    
    events = get_events(subjects{s}, false); % based on events
    fix_idx = strcmp({events.type},'fixation');
    enc_idx = strcmp({events.phase}, 'encode');
    enc_fix_events = events(fix_idx & enc_idx);
    
    fix_info = ev2fix(enc_fix_events, true);
    
    %     p_revisit(s) = nanmean(fix_info.repeated == 0); % for reporting
    %     p_revisit_init(s) = nanmean(fix_info.repeated == 1); % for reporting
    %     p_revisit_all(s) = nanmean(~isnan(fix_info.repeated)); % for reporting
    
    p_revisit(s) = nanmean(get_trial_mean(fix_info, ...
        fix_info.old == 1, ...
        'repeated', ...
        false, ...
        false));
        
    p_revisit_hit(s) = nanmean(get_trial_mean(fix_info, ...
        fix_info.resp == 1, ...
        'repeated', ...
        false, ...
        false));
    
    p_revisit_miss(s) = nanmean(get_trial_mean(fix_info, ...
        fix_info.resp ~= 1, ...
        'repeated', ...
        false, ...
        false));
    
    rev_replayed_revisit(s) = nanmean(get_trial_mean(fix_info, ...
        ~isnan(fix_info.repeated), ...
        'lag_minusone', ...
        true, ...
        false));
    
    rev_replayed_other(s) = nanmean(get_trial_mean(fix_info, ...
        isnan(fix_info.repeated), ...
        'lag_minusone', ...
        true, ...
        false));
    
    replayed_revisit(s) = nanmean(get_trial_mean(fix_info, ...
        ~isnan(fix_info.repeated), ...
        'lag_one', ...
        true, ...
        false));
    
    replayed_other(s) = nanmean(get_trial_mean(fix_info, ...
        isnan(fix_info.repeated), ...
        'lag_one', ...
        true, ...
        false));
    
end

% overall proportion that are revisited
prop_revisit = p_revisit;
prop_revisit_m = nanmean(prop_revisit);
prop_revisit_sem = nanstd(prop_revisit)/sqrt(length(prop_revisit));

% revisitation sme
revisit_hit_m = nanmean(p_revisit_hit);
revisit_hit_sem = nanstd(p_revisit_hit)/sqrt(length(p_revisit_hit));

revisit_miss_m = nanmean(p_revisit_miss);
revisit_miss_sem = nanstd(p_revisit_miss)/sqrt(length(p_revisit_miss));

stats.revisit_sme = mes(p_revisit_hit', ...
    p_revisit_miss', ...
    'hedgesg', ...
    'isDep', 1);

% revisitation subs. replay
revisit_replay_m = nanmean(replayed_revisit);
revisit_replay_sem = nanstd(replayed_revisit)/sqrt(length(replayed_revisit));

other_replay_m = nanmean(replayed_other);
other_replay_sem = nanstd(replayed_other)/sqrt(length(replayed_other));

stats.revisit_replay = mes(replayed_revisit', ...
    replayed_other', ...
    'hedgesg', ...
    'isDep', 1);

% revisitation subs. rev replay
revisit_rev_replay_m = nanmean(rev_replayed_revisit);
revisit_rev_replay_sem = nanstd(rev_replayed_revisit)/sqrt(length(rev_replayed_revisit));

other_rev_replay_m = nanmean(rev_replayed_other);
other_rev_replay_sem = nanstd(rev_replayed_other)/sqrt(length(rev_replayed_other));

stats.revisit_rev_replay = mes(rev_replayed_revisit', ...
    rev_replayed_other', ...
    'hedgesg', ...
    'isDep', 1);


fig('FontSize', 12, ...
'Font', 'Arial', ...
'border', 'on');

set(gcf,'Position', [39.8727   11.0596    6.6410   10.5040]);

subplot(4,1,1);
[he, hs] = errorplot(gca, replayed_revisit', ...
    replayed_other', ...
    {'Revisitation','Other'}, ...
    [0 0.15]);

xlabel('Location Type');
ylabel('P(Replay)');

subplot(4,1,2);
[he, hs] = errorplot(gca, rev_replayed_revisit', ...
    rev_replayed_other', ...
    {'Revisitation','Other'}, ...
    [0 0.15]);

xlabel('Location Type');
ylabel('P(Replay)');


end

function m = get_trial_mean(fixations, idx, field_name, binary_flag, sum_flag)

% get the unique trials for the passed in fixation information
[ut, ~, trial] = unique([fixations.trial; fixations.block]', 'rows');

m = nan(size(ut,1),1);
for t = 1:size(ut,1) % iterate over trials
    
    % get the fixation numbers that meet the condition
    this_trial_idx = trial==t;
    denom = length(unique(fixations.seq(this_trial_idx))); %
    
    this_idx = find(this_trial_idx & idx');
    
    if isempty(this_idx) && ~binary_flag
        continue
    end
   
    tmp_seq = fixations.seq(this_idx);
    [~,ia] = unique(tmp_seq);
    this_idx = this_idx(ia);
   
    % find all fixations to the same positions
    this_m = nan(length(this_idx), 1);
    for i = 1:length(this_idx)
        
        this_seq = fixations.seq(this_idx(i));
        
        % same location as this one
        to_eval = this_trial_idx & fixations.seq' == this_seq;
        
        if binary_flag % values ==1
            if sum_flag
                this_m(i) = any(fixations.(field_name)(to_eval));
            else
                this_m(i) = any(fixations.(field_name)(to_eval));
            end
        else
            if sum_flag
                this_m(i) = any(~isnan(fixations.(field_name)(to_eval)));
            else
                this_m(i) = any(~isnan(fixations.(field_name)(to_eval))); %
            end
        end
    end

    if sum_flag
        m(t) = nansum(this_m)/length(this_m); % nanmean
    else
        m(t) = nansum(this_m)/denom;
    end
end


end

function [he, hs] = errorplot(ax, cond1, cond2, label, ylim)


color = 'k';
he = add_errorbars(ax, cond1, cond2, color, label);
set(gca, 'YLim', ylim);

color = [.7 .7 .7];
hs = add_subjects(ax, cond1, cond2, color);

[hl, ht] = add_siglines(ax, cond1, cond2);


end

function [hl, ht] = add_siglines(ax, cond1, cond2)

hold on;

[p, h, stats] = signrank(cond1, cond2, 'method', 'approximate');

if p < 0.05
    txt = '*';
elseif p < .1
    txt = '~';
else
    hl = []; ht = [];
    return
end

yd = diff(get(gca,'YTick')); yd = yd(1);
ym = [cond1 cond2]; ym = max(ym(:));
hl = plot(ax, [1 2], [ym+.5*yd ym+.5*yd], 'k-');
hl.LineWidth = 2;
ht = text(ax, 1.5, ym+.6*yd, txt, 'FontSize', 14);
ht.HorizontalAlignment = 'center';
end

function he = add_errorbars(ax, cond1, cond2, color, label)

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
    'YLim', [0 0.6], ...
    'box', 'off', ...
    'LineWidth', 2, ...
    'XTickLabel', label);

end

function hs = add_subjects(ax, cond1, cond2, color)

x = [1 2];
for s = 1:size(cond1,1)
    y = [cond1(s) cond2(s)];
    hold on;
    hs(s) = plot(ax, x, y, 'o-');
    hs(s).Color = color;
    hs(s).MarkerFaceColor = color;
    hs(s).MarkerEdgeColor = color;
    hs(s).MarkerSize = 2;
    
end


end

