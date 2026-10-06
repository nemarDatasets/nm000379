function [hl1, hp1, hl2, hp2] = plot_ft_res(dat1, dat2, stats, ylim, ...
    colors, labels, ax)
% PLOT_FT_RES plots permutation results for bosc timeseries data from fieldtrip
%
% INPUTS:  dat1 - structure, fieldtrip data structure for condition 1
%
%          dat2 - structure, fieldtrip data structure for condition 2
%
%          stats - structure, fieldtrip stats structure comparing dat1 to dat2
%
%          ylim  - double, vector containing y limits for plotting
%
%          colors - double, 2x3 array containing RGB values for each condition
%
%          labels - cell array, contains string labels for each condition
%
%          ax - optional axis handle to plot onto, if missing makes a new figure
%
% OUTPUTS: hl1 - handle for condition 1 line
%
%          hp1 - handle for condition 1 patch
%
%          hl2 - handle for condition 2 line
%
%          hp2 - handle for condition 2 patch

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

if ~exist('colors','var')
    colors(1,:) = [0 255 0]/255;
    colors(2,:) = [102 102 102]/255;
end

if ~exist('labels','var')
    labels = {'Revisitation','Other'};
end

if ~exist('ax', 'var')
fig('FontSize', 14, ...
'Font', 'Arial', ...
'border','on', ...
'units', 'centimeters', ...
'width', 12, ...
'height', 8);
ax = gca;
end
axes(ax);

x = stats.time*1000;
full_time = -750:1:750;

y1 = nanmean(dat1(:, ismember(full_time, x)));
e1 =  nanstd(dat1)./sqrt(size(dat1,1));
e1 =  e1(:, ismember(full_time, x));

y2 = nanmean(dat2(:, ismember(full_time, x)));
e2 = nanstd(dat2)./sqrt(size(dat2,1));
e2 = e2(:, ismember(full_time, x));

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
         'LineWidth', 2);
     
xlabel('Time (ms)');
ylabel('Theta Prevalence (Z)');

hl3 = plot([0 0], ylim, '--');
hl3.Color = 'k'; hl3.LineWidth = 2;
uistack(hl3,'bottom')

if any(stats.mask)
    
    [l, n] = bwlabeln(stats.mask);
    
    for i = 1:n
        xsig = x(l==i);
        hsig = plot(xsig, min(get(gca,'YLim'))*ones(size(xsig))+.025, 'r-');
        hsig.LineWidth = 4;
    end
    
end

hleg = legend([hl1 hl2], labels);
hleg.Box = 'off'; hleg.Location = 'best';

end

