function [psi_diff, psi_1, psi_2, null_mean, null_std, psiz_1, psiz_2] = perm_diff_psi(c1_data, c2_data, freqbins)
% RUN_PSI_STATS computes phase-slope index and permutation stats and the
% subject level prior to group analysis
%
% INPUTS:  c1_data  - fieldtrip power and csd representation for condition 1
%
%          c2_data  - fieldtrip power and csd representation for condition 2
%
%          freqbins - double vector, frequency bins to comput PSI over
%
% OUTPUTS: psi_diff - double, standardized (z-scored) difference
%
%          psi_1 - double, observed psi for condition 1
%
%          psi_2 - double, observed psi for condition 2
%
%          null_mean - double, mean of null distribution used to
%          standardize difference
%
%          null_std - double, standard dev of null distribution used to
%          standardize difference
%
%          psiz_1 - double, standardized psi for condition 1
%
%          psiz_2 - double, standardized psi for condition 2

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

% This function was modified from the original data2psiX.m by Mike X Cohen.

n_permutes = 1000;

% data dimensions
[ntrial_1, nchan, nchan, nhz, npnts] = size(c1_data.crsspctrm);
[ntrial_2, nchan, nchan, nhz, npnts] = size(c2_data.crsspctrm);

full_cs = cat(1, c1_data.crsspctrm, c2_data.crsspctrm);

cs_1 = c1_data.crsspctrm; % can easily check this by comparing keeptrials output from ft_freqanalysis
cs_2 = c2_data.crsspctrm;

% define frequencies
hz  = c1_data.freq;

null_mean = zeros(size(freqbins,1),npnts);
null_std = zeros(size(freqbins,1),npnts);

%% compute PSI

% convert full_cs to gpuArray for faster processing
full_cs = gpuArray(full_cs);

for freqbini=1:size(freqbins,1) % this is currently always 1
    
    % find freq indices of requested frequency bands
    freqidx = dsearchn(hz',freqbins(freqbini,1)) : dsearchn(hz',freqbins(freqbini,2));
    nfidx   = length(freqidx);
    
    if nfidx<4
        warning('There are fewer than four frequency bins. Consider using longer time windows or wider frequency bands.')
    end
    
    [psi_1, psi_2, psiz_1, psiz_2] = return_psi(squeeze(mean(cs_1)), squeeze(mean(cs_2)), nchan, nfidx, npnts, freqidx, freqbini, true);
    psi_diff = psi_1 - psi_2;
    
    
    %% permutation testing
    
    nulldist = nan(size(psi_diff,2), 1000);
    
    for permi=1:n_permutes
        
        perm_idx = randperm(size(full_cs,1));
        null_c1 = squeeze(mean(full_cs(perm_idx(1:ntrial_1), :, :, :, :)));
        null_c2 = squeeze(mean(full_cs(perm_idx(ntrial_1+1:end), :, :, :, :)));
        
        [null_1, null_2] = return_psi(gather(null_c1), ...
            gather(null_c2), ...
            nchan, ...
            nfidx, ...
            npnts, ...
            freqidx, ...
            freqbini, ...
            false);
        
        nulldist(:,permi) = null_1 - null_2;
        
    end
    
    null_mean(freqbini,:) = mean(nulldist,2);
    null_std(freqbini,:) = std(nulldist,[],2);
    
    psi_diff(freqbini,:) =  (psi_diff(freqbini,:) - null_mean(freqbini,:)) ...
        ./ null_std(freqbini,:);
   
    
end

end

function [psi_1, psi_2, psiz_1, psiz_2] = return_psi(cs_1, cs_2, nchan, ...
    nfidx, npnts, freqidx, freqbini, run_perm)

% temporary phase-frequency matrix from this frequency band
psiz_1 = nan;
psiz_2 = nan;

pp_1 = zeros(nfidx,npnts);
pp_2 = pp_1;

for fi=1:nfidx
    for t = 1:npnts
        
        d1 = sqrt(cs_1(1,1,freqidx(fi),t)*cs_1(2,2,freqidx(fi),t));
        d2 = sqrt(cs_2(1,1,freqidx(fi),t)*cs_2(2,2,freqidx(fi),t));
        
        pp_1(fi,t) = cs_1(1,2,freqidx(fi),t)/d1;
        pp_2(fi,t) = cs_2(1,2,freqidx(fi),t)/d2;

    end
end

% average phase slope for each frequency band, for each condition of
% interest

psi_1 = nan(1,npnts);
psi_2 = psi_1;

for t = 1:npnts
    psi_1(freqbini,t) = sum(imag(conj(pp_1(1:end-1,t)).*pp_1(2:end,t)));
    psi_2(freqbini,t) = sum(imag(conj(pp_2(1:end-1,t)).*pp_2(2:end,t)));
end

if run_perm
    
    n_permutes = 1000;
    nulldist_1 = nan(1,npnts, n_permutes);
    nulldist_2 = nan(1,npnts, n_permutes);
    for t = 1:npnts
        for i=1:n_permutes
            nulldist_1(:,t,i) = sum(imag(conj(pp_1(randperm(nfidx),t)).*pp_1(randperm(nfidx),t)));
            nulldist_2(:,t,i) = sum(imag(conj(pp_2(randperm(nfidx),t)).*pp_2(randperm(nfidx),t)));
        end
    end
    
    n1_mean = mean(nulldist_1,3);
    n1_std = std(nulldist_1,[],3);
    psiz_1 = (psi_1-n1_mean) ./ n1_std;

    n2_mean = mean(nulldist_2,3);
    n2_std = std(nulldist_2,[],3);
    psiz_2 = (psi_2-n2_mean) ./ n2_std;
 
end

end