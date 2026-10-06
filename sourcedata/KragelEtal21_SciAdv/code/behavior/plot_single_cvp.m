function [hl, hp] = plot_single_cvp(cvp, ax)
% PLOT_SINGLE_CVP plots comparison of CVP for hits and misses
%
% INPUTS:  cvp - double, array (subject x lag) containing cvp data
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

x = -4:4;

y = cvp(:,1:9);
y_m = squeeze(nanmean(cvp(:,1:9),1));
e = squeeze(nanstd(cvp(:,1:9)))/sqrt(size(cvp,1));

hl = errorbar(ax, x, y_m, e); hold on
hl.LineWidth = 2;
hl.Color = 'k';

hp = plot(x' + .5*rand(size(y')) - .25, y', 'o-');
for p = 1:length(hp)
    hp(p).MarkerFaceColor = [.7 .7 .7];
    hp(p).MarkerEdgeColor = [.7 .7 .7];
    hp(p).Color = [.7 .7 .7];
    hp(p).MarkerSize = 2;
end

uistack(hp, 'bottom');

set(ax, 'XLim', [-5 5], ...
    'YLim', [ 0 .4], ...
    'LineWidth', 2, ...
    'XTick', -4:2:4, ...
    'box','off');

xlabel('Lag');
ylabel('CVP');

end

