function plot_supp_osc()
% PLOT_SUPP_OSC plots hippocampal theta oscillations and 1\f fits for
% supplementary figure
%
% INPUTS:  none, this function loads all data from temporary files
%
% OUTPUTS: none, this function generates a figure for display

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
base_dir = '..\..\scratch_data\';

colors = [ 0 0 0; ...
           255 153 51 ; ...
           51 102 204; ...
           255 0 0 ; ...
           102 0 204; ...
           51 255 102]/255;

figure;ecnt=0;

for s = 1:length(subjects)

    res = get_bosc(subjects{s}, base_dir);
    freqs = res(1).freqs;

    ap_ps = []; ps = [];

    for e = 1:size(res,2) % electrodes

        ecnt=ecnt+1;
        subplot(4,8,ecnt);
        
        tmp_ap_ps = nan(size(res,1),50);
        tmp_ps = nan(size(res,1),50);
        for b = 1:size(res,1)
            tmp_ap_ps(b,:) = log10(res(b, e).ap_ps);
            tmp_ps(b,:) = res(b, e).ps;
        end
        
        ap_ps(e,:) = nanmean(tmp_ap_ps);
        ps(e,:) = nanmean(tmp_ps);
        
        h(1) = plot(freqs, ap_ps(e,:), '--');
        hold on;
        
        h(2) = plot(freqs, ps(e,:), '-');
        
        h(1).LineWidth = 2;
        h(1).Color = colors(s,:);
        h(2).LineWidth = 2;
        h(2).Color = colors(s,:);
        
        set(gca, 'Box', 'off', ...
                 'FontSize', 14, ...
                 'XTick', [0 10 20 30], ...
                 'XLim', [1 30]);
    end

end

end

function [pe, ps] = get_pe(theta_res, freqs)

ps = nan(length(theta_res), length(freqs));
pe = [];

for s = 1:length(theta_res)
    
    if isnan(theta_res(s).p_ep)
        continue
    end
    
    y = squeeze(nanmean(theta_res(s).p_ep));
    pe = cat(1, pe, y);
    
    is_lf = theta_res(s).freqs > 1 & theta_res(s).freqs < 10;
    y_idx = unique(theta_res(s).chans(is_lf));
    
    all_idx = unique(theta_res(s).chans);
    
    
    if isempty(all_idx)
        continue
    end
    
    ps(s,:) = nanmean(y(all_idx,:));
    
end

end
