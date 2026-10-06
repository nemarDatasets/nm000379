function res = run_subj_plv_baseline(subj)
% RUN_SUBJ_PLV_BASELINE computes phase-locking value connectivity metric at the
% subject level and saves to an intermediate file
%
% INPUTS:  subj - string, subject identifier to run analysis on
%
% OUTPUTS: res - results structure with plv and permutation data

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

load(['..\..\scratch_data\connectivity\' subj '_data_noar_wvis.mat'], ...
    'data', 'events', 'fix_events');

% updata data to include relevant behavioral info
data = update_trialinfo(data, fix_events);

% config for fourier representation 
cfg        = [];
cfg.method = 'mtmconvol';
cfg.taper  = 'hanning';

fboilim = round([1 10] .* 8.1920) + 1;
fboi    = fboilim(1):1:fboilim(2);
foi     = (fboi-1) ./ 8.1920;
cfg.foi = foi(1:2:75);

cfg.tapsmofrq = 2;
cfg.t_ftimwin = 4 ./ cfg.foi; % 4 cycles

cfg.toi        = data.time{1}(data.time{1} >= 5.25 & data.time{1} <= 6.75);
cfg.output     = 'fourier';
cfg.pad        = 'nextpow2';
cfg.channelcmb = 'all';
cfg.trials     = data.trialinfo(:,1) == 2 & (data.trialinfo(:,3) == 1 | data.trialinfo(:,3) == 2); % fixations
cfg.channel    = data.label(is_hc | is_dan);
fix_freq       = ft_freqanalysis(cfg, data);

% data no longer needed
clear data

chansel = find_connections(subj, fix_freq);

rois = get_plv_rois(subj);

angles = single(angle(fix_freq.fourierspctrm));

% free some memory
clear fix_freq

aplv = nan(size(chansel,1), size(angles, 3), size(angles, 4));
aplv_z = nan(size(chansel,1), size(angles, 3), size(angles, 4));
anull_m = nan(size(chansel,1), size(angles, 3), size(angles, 4));
anull_s = nan(size(chansel,1), size(angles, 3), size(angles, 4));

for c = 1:size(chansel,1)
    
    % compute plv
    plv = abs(sum(exp(1i* (angles(:, chansel(c,1), :, :) - angles(:, chansel(c,2), :, :)))))/ ...
        size(angles,1);
    plv = squeeze(plv);
    
    null_plv = nan(size(plv,1), size(plv,2), 1000);
    all_null_angles = squeeze(angles(:, [chansel(c,1) chansel(c,2)], :, :));

    parfor t = 1:size(angles,4)
        null_angles = squeeze(all_null_angles(:,:,:,t));
        null_angles = repmat(null_angles, [1 1 1 1000]);
        
        for p = 1:1000
            shuff = randperm(size(null_angles,1));
            null_angles(:,1,:,p) = null_angles(shuff, 1, :, p);
        end
        
        null_plv(:,t,:) = abs(sum(exp(1i* (null_angles(:,1,:,:) - null_angles(:, 2, :, :)))))/ ...
            size(null_angles,1);
    end
    
    null_m = nanmean(null_plv,3);
    null_s = nanstd(null_plv,[],3);
    
    plv_z = (plv - null_m)./null_s;
    
    aplv(c,:,:,:) = plv;
    aplv_z(c,:,:,:) = plv_z;
    anull_m(c,:,:,:) = null_m;
    anull_s(c,:,:,:) = null_s;
    
end

res.plv = aplv;
res.plv_z = aplv_z;
res.null_m = anull_m;
res.null_s = anull_s;

res.rois = rois;
res.chansel = chansel;

end


function chansel = find_connections(subj, fix_freq)


nchan = length(fix_freq.label);
chanindx = tril(true(nchan), -1);
cmbindx1 = repmat((1:nchan)', [1 nchan]);
cmbindx2 = repmat((1:nchan),  [nchan 1]);

chansel = [cmbindx1(chanindx) cmbindx2(chanindx)];

is_hc = return_hc(subj, fix_freq.label);

% only find combinations that are between Hc and DAN/VIS, (not within Hc and
% within DAN/VIS)

between_region = false(1, size(chansel,1));
for c = 1:size(chansel,1)
    if sum(is_hc(chansel(c,:))) == 1
        between_region(c) = true;
    end
end

chansel = chansel(between_region,:);

end
