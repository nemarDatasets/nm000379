function plot_osc_rois()
% PLOT_OSC_ROIS plots p_episode as a function of frequency for each ROI
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

% finds theta peaks, etc. for LTC, DAN, VIS, HC

subjects = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6'};
freqs = logspace(log10(1), log10(40), 50);

% ltc, hipp, and visual rois in separate mat files
load('..\..\scratch_data\bosc\ltc_theta_res.mat');
ltc_theta = append_subs(theta_res, subjects);

load('..\..\scratch_data\bosc\vis_theta_res.mat');
vis_theta = append_subs(theta_res, subjects([1 2 4 5 6]));

load('..\..\scratch_data\bosc\theta_res.mat');
hc_theta = append_subs(theta_res, subjects);

clear theta_res

p_e = nan(length(subjects), length(freqs));
p_all = [];

for s = 1:length(subjects)
    
    subj = subjects{s};

    pairs = get_pairs(subj); % get all bipolars
    
    % get pairs in ltc mat
    [is_ltc, ltc_chans] = get_ltc_pairs(subj, pairs);
    
    % get pairs in vis mat
    [out_is_vis, all_vis_chans, dan_sub, vis_sub] = get_vis_pairs(subj, pairs);
    
    dan_chans = all_vis_chans(dan_sub);
    vis_chans = all_vis_chans(vis_sub);
    
    hc_chans = pairs(return_hc(subj, pairs));
    
    
    [is_hc, is_dan, is_vis] = get_roi_labels(subj, ltc_chans);
    
    % return ltc
    roi_idx = true(size(ltc_chans));
    roi_idx(is_hc | is_dan | is_vis) = false;
    
    s_idx = strcmp({ltc_theta.subject}, subj);
    if any(s_idx)
        [ltc(s).p_ep, ltc(s).p_ep_m] = get_p_ep(ltc_theta(s_idx), roi_idx);
    end
    
    % return dan/vis
    [is_hc, is_dan, is_vis] = get_roi_labels(subj, all_vis_chans);
    
    % vis
    roi_idx = false(size(all_vis_chans));
    roi_idx(dan_sub) = true;
    roi_idx(is_hc) = false;
    
    s_idx = strcmp({vis_theta.subject}, subj);
    if any(s_idx)
        [vis(s).p_ep, vis(s).p_ep_m] = get_p_ep(vis_theta(s_idx), roi_idx);
    end
    
    % dan
    roi_idx = false(size(all_vis_chans));
    roi_idx(vis_sub) = true;
    roi_idx(is_hc) = false;
    
    s_idx = strcmp({vis_theta.subject}, subj);
    if any(s_idx)
        [dan(s).p_ep, dan(s).p_ep_m] = get_p_ep(vis_theta(s_idx), roi_idx);
    end
    
    % hipp
    roi_idx = true(size(hc_chans));
    s_idx = strcmp({hc_theta.subject}, subj);
    if any(s_idx)
        [hc(s).p_ep, hc(s).p_ep_m] = get_p_ep(hc_theta(s_idx), roi_idx);
    end
    
end

vis_pep = cat(1, vis.p_ep_m);
dan_pep = cat(1, dan.p_ep_m);
ltc_pep = cat(1, ltc.p_ep_m);
hc_pep  = cat(1, hc.p_ep_m);

% first plot, hc, ltc

colors = [0 123 255; ... hc
          255 13 255; ... ltc
            255 0 0; ... vis
          255 191 13]/255; % ltc
            
x = freqs;

y = [nanmean(hc_pep); nanmean(ltc_pep); nanmean(vis_pep); nanmean(dan_pep)];
e = [nanstd(hc_pep)./sqrt(size(hc_pep,1)); ...
    nanstd(ltc_pep)./sqrt(size(ltc_pep,1)); ...
    nanstd(vis_pep)./sqrt(size(vis_pep,1)); ...
    nanstd(dan_pep)./sqrt(size(dan_pep,1))];

fig('FontSize', 16, ...
    'Font', 'Arial', ...
    'border','off', ...
    'units', 'centimeters', ...
    'width', 8, ...
    'height', 6);

for r = 1:2
hold on;
[hl(r), hp(r)] = boundedline(x, y(r,:), e(r, :, :));
hl(r).Color = colors(r,:);
hl(r).LineWidth = 2;
hp(r).FaceColor = colors(r,:);
hp(r).FaceAlpha = .4;
end

set(gca,'LineWidth',2, ...
        'XLim',[1 40]);

xlabel('Frequency (Hz)');
ylabel('P_{episode}');

hleg = legend(hl, {'Hipp','LTC'});
hleg.Box = 'off';
hleg.Location = 'best';

% second plot, vis, dan

fig('FontSize', 16, ...
    'Font', 'Arial', ...
    'border','off', ...
    'units', 'centimeters', ...
    'width', 8, ...
    'height', 6);

for r = 3:4
hold on;
[hl(r), hp(r)] = boundedline(x, y(r,:), e(r, :, :));
hl(r).Color = colors(r,:);
hl(r).LineWidth = 2;
hp(r).FaceColor = colors(r,:);
hp(r).FaceAlpha = .4;
end

set(gca,'LineWidth',2, ...
        'XLim',[1 20]);

xlabel('Frequency (Hz)');
ylabel('P_{episode}');

hleg = legend(hl(3:4), {'VN','DAN'});
hleg.Box = 'off';
hleg.Location = 'best';

end

function [pep, pep_m] = get_p_ep(theta_res, roi_idx)
% returns p_ep for all channels in theta_res, and the mean over roi indices
    pep = squeeze(nanmean(theta_res.p_ep)); % avg over blocks
    
    is_lf = theta_res.freqs > 1 & theta_res.freqs < 10;
    y_idx = false(size(roi_idx));
    y_idx(theta_res.chans(is_lf)) = true;
    
    all_idx = false(size(roi_idx));
    all_idx(theta_res.chans) = true;
    
    pep_m = nanmean(pep(all_idx & y_idx & roi_idx,:),1);
    
end

function [is_hc, is_dan, is_vis] = get_roi_labels(subj, pairs)

[is_dan, is_vis] = return_yeo(subj, pairs);
is_hc = return_hc(subj, pairs);

is_dan(is_hc) = false;
is_vis(is_hc) = false;

end

function [is_ltc, ltc_chans] = get_ltc_pairs(subj, pairs)

pairs = strrep(pairs,'*','');
pairs = strrep(pairs,' ','');
[is_ltc, pairs] = return_ltc(subj, pairs);
ltc_idx = find(is_ltc);
elec_idx = unique(ltc_idx);
ltc_chans = pairs(elec_idx);

end

function [out_is_vis, vis_chans, dan_sub, vis_sub] = get_vis_pairs(subj, pairs)

pairs = strrep(pairs,'*','');
pairs = strrep(pairs,' ','');

[is_dan, is_vis, pairs] = return_yeo(subj, pairs);

dan_idx = find(is_dan);
vis_idx = find(is_vis);

elec_idx = [dan_idx; vis_idx];
elec_idx = unique(elec_idx);

out_is_vis = false(size(is_vis));
out_is_vis(elec_idx) = true;

vis_chans = pairs(out_is_vis);
dan_sub = is_dan(out_is_vis);
vis_sub = is_vis(out_is_vis);

end

function pairs = get_pairs(subj)

server_dir = ['..\..\data\eeg\' subj '\'];
fname = [server_dir subj '_sl_block1_bp.h5'];
pairs = h5read(fname, '/pairs');

end

function st = append_subs(st, subjects)

for s = 1:length(subjects)
    st(s).subject = subjects{s};
end

end
