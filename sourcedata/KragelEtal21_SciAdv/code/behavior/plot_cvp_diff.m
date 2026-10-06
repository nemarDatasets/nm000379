function [hl, hp] = plot_cvp_diff(cvp_hit, cvp_miss, ax)
% PLOT_CVP_DIFF plots comparison of CVP for hits and misses
%
% INPUTS:  cvp_hit - double, array (subject x lag) containing cvp for hits
%
%          cvp_miss - double, array (subject x lag) containing cvp for misses
%
%          ax - axis handle, specifies which axis to plot on (optional)
%
% OUTPUTS: hl - line handle, for modifying after plotting
%
%          hp - patch handle, for modifying after plotting

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

if ~exist('ax', 'var')
    
    fig('FontSize', 12, ...
        'Font', 'Arial', ...
        'border', 'on', ...
        'units', 'centimeters', ...
        'width', 12);
    
    ax = gca;
    
end

x = [1 2 4 5];

y = [cvp_hit(:,4) cvp_miss(:,4) ...
    cvp_hit(:,6) cvp_miss(:,6)];

y_m = [nanmean(cvp_hit(:,4)) nanmean(cvp_miss(:,4)) ...
    nanmean(cvp_hit(:,6)) nanmean(cvp_miss(:,6))];

e = [nanstd(cvp_hit(:,4))./sqrt(size(cvp_hit,1)) nanstd(cvp_miss(:,4))./sqrt(size(cvp_hit,1)) ...
    nanstd(cvp_hit(:,6))./sqrt(size(cvp_hit,1)) nanstd(cvp_miss(:,6))./sqrt(size(cvp_hit,1))];

% group err
for i = 1:4
    hl(i) = errorbar(ax, x(i), y_m(i), e(i), 'o');
    hold on
    hl(i).LineWidth = 2;
    if ismember(i, [1 3])
        hl(i).Color = 'r';
        hl(i).MarkerFaceColor = 'r';
        hl(i).MarkerEdgeColor = 'r';
         hl(i).MarkerSize = 4;
    else
        hl(i).Color = 'b';
        hl(i).MarkerFaceColor = 'b';
        hl(i).MarkerEdgeColor = 'b';
        hl(i).MarkerSize = 4;
    end
end
    
% group lines
idx_sets = [1 2; 3 4];
for i = 1:size(idx_sets,1)
    
    hl2(i) = plot(ax, x(idx_sets(i,:)), ...
        y_m(idx_sets(i,:)), '-');
    hold on
    hl2(i).LineWidth = 2;
    hl2(i).Color = 'k';
    uistack(hl2, 'bottom');
end

% subject points
idx_sets = [1 2; 3 4];
for i = 1:size(idx_sets,1)
    hp = plot(x(idx_sets(i,:))' + .3*rand(size(y(:, idx_sets(i,:))))', ...
        y(:,idx_sets(i,:))',...
        'o-');
    for p = 1:length(hp)
        hp(p).MarkerFaceColor = [.7 .7 .7];
        hp(p).MarkerEdgeColor = [.7 .7 .7];
        hp(p).MarkerSize = 2;
        hp(p).Color = [.7 .7 .7];
    end
    uistack(hp, 'bottom');
    
    % add siglines
    
    [h, ht] = add_siglines(ax, ...
        y(:,idx_sets(i,1)), ...
        y(:,idx_sets(i,2)), ...
        x(idx_sets(i,1)), ...
        x(idx_sets(i,2)), ...
        .5);
    
    if ~isempty(ht)
        ht.HorizontalAlignment = 'center';
        ht.FontWeight = 'bold';
    end
end

set(ax, 'XLim', [0 9], ...
    'YLim', [ 0 .4], ...
    'LineWidth', 2, ...
    'XTick', [1.5 4.5], ...
    'XTickLabel', {'-1','+1'},...
    'box','off');

xlabel('Lag');
ylabel('CVP');

hleg = legend(hl(1:2), {'Hit', 'Miss'});
hleg.Location = 'NorthWest';
hleg.LineWidth = 1;
hleg.Box = 'off';

end

function [hl, ht] = add_siglines(ax, cond1, cond2, x1, x2, d)

hold on;

p = signrank(cond1, cond2);

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
hl = plot(ax, [x1 x2], [ym+d*yd ym+d*yd], 'k-');
hl.LineWidth = 2;
ht = text(ax, mean([x1 x2]), ym+(d+d/5)*yd, txt, 'FontSize', 14);

end