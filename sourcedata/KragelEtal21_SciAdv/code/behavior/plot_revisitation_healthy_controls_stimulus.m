function plot_revisitation_healthy_controls_stimulus(dat, study, init_idx)
% PLOT_REVISITATION_HEALTHY_CONTROLS_STIMULUS generates plots showing relation between
% revisitation and subsequent fixations sequence reinstatement
%
% INPUTS:  data - structure, organized fixation sequences reinstatement for 
%                 each fixations type (see run_revisitation_controls.m)
%
% INPUTS:  data - structure, organized fixation sequences reinstatement for 
%                 each fixations type (see run_revisitation_controls.m)
%
%          study - double, specifies which dataset the subjects in dat belong to
%
%          init_idx - double, which index (lag) to examine
%
% OUTPUTS: none, generates figure

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

fig('FontSize', 16, ...
    'Font', 'Arial', ...
    'border','off', ...
    'units', 'centimeters', ...
    'width', 18, ...
    'height', 14);

set(gcf,'Position',[17.9652   16.5365   10.0012    9.3133]);

subplot(2,1,1);

ax = gca;

x = [6 7];
y = dat.null_reins(study==1,:);
[he1, hs1, hp, stats(1)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);

for s = 1:length(hs1)
    hs1(s).Color = [255 153 0]/255;
    hs1(s).MarkerFaceColor = [255 153 0]/255;
    hs1(s).MarkerEdgeColor = [255 153 0]/255;
end

uistack(hp,'top');
uistack(he1,'top');


x = [7.5 8.5];
y = dat.null_reins(study==2,:);
[he1, hs1, hp, stats(2)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);

for s = 1:length(hs1)
    hs1(s).Color = [255 153 0]/255;
    hs1(s).MarkerFaceColor = [255 153 0]/255;
    hs1(s).MarkerEdgeColor = [255 153 0]/255;
end

uistack(hp,'top');
uistack(he1,'top');

x = [9 10];
y = dat.null_reins(study==3,:);
[he1, hs1, hp, stats(3)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);

for s = 1:length(hs1)
    hs1(s).Color = [255 153 0]/255;
    hs1(s).MarkerFaceColor = [255 153 0]/255;
    hs1(s).MarkerEdgeColor = [255 153 0]/255;
end

uistack(hp,'top');
uistack(he1,'top');


x = [1 2];
y = dat.null_rev_reins(study==1,:);
[he1, hs1, hp, stats(4)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);
uistack(hp,'top');
uistack(he1,'top');

for s = 1:length(hs1)
    hs1(s).Color = [153 51 204]/255;
    hs1(s).MarkerFaceColor = [153 51 204]/255;
    hs1(s).MarkerEdgeColor = [153 51 204]/255;
end

x = [2.5 3.5];
y = dat.null_rev_reins(study==2,:);
[he1, hs1, hp, stats(5)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);
uistack(hp,'top');
uistack(he1,'top');

for s = 1:length(hs1)
    hs1(s).Color = [153 51 204]/255;
    hs1(s).MarkerFaceColor = [153 51 204]/255;
    hs1(s).MarkerEdgeColor = [153 51 204]/255;
end

x = [4 5];
y = dat.null_rev_reins(study==3,:);
[he1, hs1, hp, stats(6)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);
uistack(hp,'top');
uistack(he1,'top');

for s = 1:length(hs1)
    hs1(s).Color = [153 51 204]/255;
    hs1(s).MarkerFaceColor = [153 51 204]/255;
    hs1(s).MarkerEdgeColor = [153 51 204]/255;
end

set(gca,'XLim', [0.5 10.5], ...
'YLim', [0 .6], ...
'LineWidth', 2, ...
'box', 'off', ...
'YTick', [0:.2:.6], ...
'XTick', [1.5 3 4.5 6.5 8 9.5], ...
'XTickLabel', {'D1','D2','D3','D1','D2','D3'});
ylabel('P(Reinstatement)')
xlabel('Dataset');

set(gcf,'Position',[12.8852   16.2454   21.0344    9.3133]);

hleg = legend(hp, {'Revisitation','Other'});
hleg.Location = 'northwest';
hleg.Box = 'off';

subplot(2,1,2);

ax = gca;

x = [6 7];
y = dat.m_reins(study==1,:) - dat.null_reins(study==1,:);
[he1, hs1, hp, stats(1)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);

for s = 1:length(hs1)
    hs1(s).Color = [255 153 0]/255;
    hs1(s).MarkerFaceColor = [255 153 0]/255;
    hs1(s).MarkerEdgeColor = [255 153 0]/255;
end

uistack(hp,'top');
uistack(he1,'top');


x = [7.5 8.5];
y = dat.m_reins(study==2,:) - dat.null_reins(study==2,:);
[he1, hs1, hp, stats(2)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);

for s = 1:length(hs1)
    hs1(s).Color = [255 153 0]/255;
    hs1(s).MarkerFaceColor = [255 153 0]/255;
    hs1(s).MarkerEdgeColor = [255 153 0]/255;
end

uistack(hp,'top');
uistack(he1,'top');

x = [9 10];
y = dat.m_reins(study==3,:) - dat.null_reins(study==3,:);
[he1, hs1, hp, stats(3)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);

for s = 1:length(hs1)
    hs1(s).Color = [255 153 0]/255;
    hs1(s).MarkerFaceColor = [255 153 0]/255;
    hs1(s).MarkerEdgeColor = [255 153 0]/255;
end

uistack(hp,'top');
uistack(he1,'top');


x = [1 2];
y = dat.m_rev_reins(study==1,:) - dat.null_rev_reins(study==1,:);
[he1, hs1, hp, stats(4)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);
uistack(hp,'top');
uistack(he1,'top');

for s = 1:length(hs1)
    hs1(s).Color = [153 51 204]/255;
    hs1(s).MarkerFaceColor = [153 51 204]/255;
    hs1(s).MarkerEdgeColor = [153 51 204]/255;
end

x = [2.5 3.5];
y = dat.m_rev_reins(study==2,:) - dat.null_rev_reins(study==2,:);
[he1, hs1, hp, stats(5)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);
uistack(hp,'top');
uistack(he1,'top');

for s = 1:length(hs1)
    hs1(s).Color = [153 51 204]/255;
    hs1(s).MarkerFaceColor = [153 51 204]/255;
    hs1(s).MarkerEdgeColor = [153 51 204]/255;
end

x = [4 5];
y = dat.m_rev_reins(study==3,:) - dat.null_rev_reins(study==3,:);
[he1, hs1, hp, stats(6)] = errorplot(ax, x, y(:, init_idx), y(:,2), {'Revisitation','Other'}, [0 1]);
uistack(hp,'top');
uistack(he1,'top');

for s = 1:length(hs1)
    hs1(s).Color = [153 51 204]/255;
    hs1(s).MarkerFaceColor = [153 51 204]/255;
    hs1(s).MarkerEdgeColor = [153 51 204]/255;
end

set(gca,'XLim', [0.5 10.5], ...
'YLim', [-.2 .8], ...
'LineWidth', 2, ...
'box', 'off', ...
'YTick', [0:.2:.8], ...
'XTick', [1.5 3 4.5 6.5 8 9.5], ...
'XTickLabel', {'D1','D2','D3','D1','D2','D3'});
ylabel('\Delta P(Reinstatement)')
xlabel('Dataset');

end

function [he, hs, hp, stats] = errorplot(ax, x, cond1, cond2, label, ylim)


color = 'k';
[he, hp] = add_errorbars(ax, x, cond1, cond2, color, label);
set(gca, 'YLim', ylim);

color = [.7 .7 .7];
hs = add_subjects(ax, x, cond1, cond2, color);

[hl, ht, stats] = add_siglines(ax, x, cond1, cond2);


end

function [hl, ht, stats] = add_siglines(ax, x, cond1, cond2)

hold on;

[h, p, c, st] = ttest(cond1, cond2);

stats = mes(cond1, cond2, 'hedgesg', 'isdep', true);

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
    return
end

yd = diff(get(gca,'YTick')); yd = yd(1);
ym = [cond1 cond2]; ym = max(ym(:));
hl = plot(ax, x, [ym+.25*yd ym+.25*yd], 'k-');
hl.LineWidth = 2;

ht = text(ax, mean(x), ym+.3*yd, txt, 'FontSize', 14);
ht.FontWeight = 'bold';
ht.HorizontalAlignment = 'center';

end

function [he, hp] = add_errorbars(ax, x, cond1, cond2, color, label)

y = [nanmean(cond1), nanmean(cond2)];
e = [nanstd(cond1)./sqrt(sum(~isnan(cond1))), ...
    nanstd(cond2)./sqrt(sum(~isnan(cond2)))];

he = errorbar(ax, x, y, e);
he.LineStyle = '-';
he.Marker = 'none';
he.Color = color; he.LineWidth = 2;
he.MarkerFaceColor = color;
he.MarkerEdgeColor = color;
he.MarkerSize = 4;
hold on;
hp(1) = plot(x(1), y(1), 'o');
hp(1).MarkerFaceColor = [0 1 0];
hp(1).MarkerEdgeColor = [0 1 0];
hp(1).MarkerSize = 10;

hp(2) = plot(x(2), y(2), 'o');
hp(2).MarkerFaceColor = [.8 .8 .8];
hp(2).MarkerEdgeColor = [.8 .8 .8];
hp(2).MarkerSize = 10;

end

function hs = add_subjects(ax, x, cond1, cond2, color)

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
