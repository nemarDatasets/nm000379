function revisit_res = run_figrim_revisitation_memory(res, null_res)
% RUN_FIGRIM_REVISITATION_MEMORY organizes figrim results for statistical comparison
%
% INPUTS:  res - structure, figrim results structure with fixation-sequence
%                replay information
%
%          null_res - structure, figrim results structure with fixation-sequence
%                replay information based on similar viewing across subjects
%
%
% OUTPUTS: revisit_res - structure, data structure containing average fixation-
%                sequence replay information broken down by condition

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

% initialize arrays
hit_m_viewed = nan(size(res.n_fix_hit,1), size(res.n_fix_hit,2), 2);
hit_m_reins = nan(size(res.n_fix_hit,1), size(res.n_fix_hit,2), 2);
hit_m_rev_reins = nan(size(res.n_fix_hit,1), size(res.n_fix_hit,2), 2);
hit_n_viewed = nan(size(res.n_fix_hit,1), size(res.n_fix_hit,2), 2);

miss_m_viewed = nan(size(res.n_fix_hit,1), size(res.n_fix_hit,2), 2);
miss_m_reins = nan(size(res.n_fix_hit,1), size(res.n_fix_hit,2), 2);
miss_m_rev_reins = nan(size(res.n_fix_hit,1), size(res.n_fix_hit,2), 2);
miss_n_viewed = nan(size(res.n_fix_hit,1), size(res.n_fix_hit,2), 2);

for s = 1:size(res.n_fix_hit,1) % subject
    for f = 1:size(res.n_fix_hit,2) % image/file
        
        % start with initial viewing only, just organize what we need
        % to compute
        
        n_fix_hit = res.n_fix_hit(s,f);
        viewed_hit = res.sub_view_hit(s,f,1:n_fix_hit);
        reinstated_hit = res.sub_reins_hit(s,f,1:n_fix_hit);
        rev_reinstated_hit = res.sub_rev_reins_hit(s,f,1:n_fix_hit);
        replay_hit = res.enc_replay_hit(s,f,1:n_fix_hit);
        repeated_hit = res.enc_repeated_hit(s,f,1:n_fix_hit);
        enc_ord_hit = res.enc_ord_hit(s,f,1:n_fix_hit);
        
        n_fix_miss = res.n_fix_miss(s,f);
        viewed_miss = res.sub_view_miss(s,f,1:n_fix_miss);
        reinstated_miss = res.sub_reins_miss(s,f,1:n_fix_miss);
        rev_reinstated_miss = res.sub_rev_reins_miss(s,f,1:n_fix_miss);
        replay_miss = res.enc_replay_miss(s,f,1:n_fix_miss);
        repeated_miss = res.enc_repeated_miss(s,f,1:n_fix_miss);
        enc_ord_miss = res.enc_ord_miss(s,f,1:n_fix_miss);
        
        if n_fix_hit > 0
        [hit_m_viewed(s,f,:), hit_m_reins(s,f,:), hit_m_rev_reins(s,f,:), ...
            hit_n_viewed(s,f,:)] = compute_res(viewed_hit, ...
            repeated_hit, ...
            reinstated_hit, ...
            rev_reinstated_hit, ...
            replay_hit, ...
            enc_ord_hit);
        end
        
        if n_fix_miss > 0
        [miss_m_viewed(s,f,:), miss_m_reins(s,f,:), miss_m_rev_reins(s,f,:), ...
            miss_n_viewed(s,f,:)] = compute_res(viewed_miss, ...
            repeated_miss, ...
            reinstated_miss, ...
            rev_reinstated_miss, ...
            replay_miss, ...
            enc_ord_miss);
        end
        
    end
end

% average over images
hit_m_viewed = squeeze(nanmean(hit_m_viewed,2));
hit_m_reins = squeeze(nanmean(hit_m_reins,2));
hit_m_rev_reins = squeeze(nanmean(hit_m_rev_reins,2));

miss_m_viewed = squeeze(nanmean(miss_m_viewed,2));
miss_m_reins = squeeze(nanmean(miss_m_reins,2));
miss_m_rev_reins = squeeze(nanmean(miss_m_rev_reins,2));

revisit_res.m_viewed = nanmean(cat(3, hit_m_viewed, miss_m_viewed),3);
revisit_res.m_reins = nanmean(cat(3, hit_m_reins, miss_m_reins),3);
revisit_res.m_rev_reins = nanmean(cat(3, hit_m_rev_reins, miss_m_rev_reins),3);

% and stimulus-driven null res
revisit_res.null_viewed = nan(size(null_res.n_fix,1), size(null_res.n_fix,2), 2);
revisit_res.null_reins = nan(size(null_res.n_fix,1), size(null_res.n_fix,2), 2);
revisit_res.null_rev_reins =  nan(size(null_res.n_fix,1), size(null_res.n_fix,2), 2);

for s = 1:size(null_res.n_fix,1)% subject
    
    for i = 1:size(null_res.n_fix,2) % image

        n_other = length(null_res.n_fix(s,i).dat);
        
        if n_other == 0 % didn't view this image
            continue
        end
        
        tmp_viewed = nan(length(null_res.n_fix(s,i).dat), 2);
        tmp_reins = nan(length(null_res.n_fix(s,i).dat), 2);
        tmp_rev_reins = nan(length(null_res.n_fix(s,i).dat), 2);
        tmp_n_viewed = nan(length(null_res.n_fix(s,i).dat), 2);
        
        for d = 1:length(null_res.n_fix(s,i).dat)
            n_fix = null_res.n_fix(s,i).dat(d);
            
            viewed = null_res.sub_view(s,i).dat(d,1:n_fix);
            reinstated =  null_res.sub_reins(s,i).dat(d,1:n_fix);
            rev_reinstated = null_res.sub_rev_reins(s,i).dat(d,1:n_fix);
            replayed = null_res.enc_replay(s,i).dat(d,1:n_fix);
            repeated = null_res.enc_repeated(s,i).dat(d,1:n_fix);
            enc_ord = null_res.enc_seq(s,i).dat(d,1:n_fix);
            
            [tmp_viewed(d,:), tmp_reins(d,:), tmp_rev_reins(d,:),
            tmp_n_viewed(d,:)] = compute_res(viewed, ...
                                             repeated, ...
                                             reinstated, ...
                                             rev_reinstated, ...
                                             replayed, ...
                                             enc_ord);
        
        end
        
        revisit_res.null_viewed(s,i,:) = nanmean(tmp_viewed);
        revisit_res.null_reins(s,i,:) = nanmean(tmp_reins);
        revisit_res.null_rev_reins(s,i,:) = nanmean(tmp_rev_reins);
        
    end

end

% and average over images
revisit_res.null_viewed =  squeeze(nanmean(revisit_res.null_viewed,2));
revisit_res.null_reins =  squeeze(nanmean(revisit_res.null_reins,2));
revisit_res.null_rev_reins =  squeeze(nanmean(revisit_res.null_rev_reins,2));

end

function [m_viewed, m_reins, m_rev_reins, n_viewed] = compute_res(viewed, ...
    repeated, reinstated, rev_reinstated, replay, enc_seq)

m_viewed = nan(1,2);
m_reins = nan(1,2);
m_rev_reins = nan(1,2);
n_viewed = nan(1,2);

[n_viewed(1), m_viewed(1)] = get_trial_stats(enc_seq, repeated==1, viewed, false);
[n_viewed(2), m_viewed(2)] = get_trial_stats(enc_seq, isnan(repeated), viewed, false);

[~, m_reins(1)] = get_trial_stats(enc_seq, repeated==1, reinstated, false);
[~, m_reins(2)] = get_trial_stats(enc_seq, isnan(repeated), reinstated, false);

[~, m_rev_reins(1)] = get_trial_stats(enc_seq, repeated==1, rev_reinstated, false);
[~, m_rev_reins(2)] = get_trial_stats(enc_seq, isnan(repeated), rev_reinstated, false);

end

function [n, m] = get_trial_stats(enc_seq, idx, dat, binary_flag)

% only once per position
unique_seq = unique(enc_seq(idx));

% find all fixations to the same positions
this_m = nan(length(unique_seq), 1);
for i = 1:length(unique_seq)
    
    this_seq = unique_seq(i);
    
    % same location as this one
    to_eval = enc_seq == this_seq;
    
    if ~binary_flag
        this_m(i) = any(dat(to_eval));
    else
        this_m(i) = any(dat(to_eval)==1);
    end
end

m = nanmean(this_m); % average response (e.g., replay at retrieval) for those locations
n = length(this_m); % number of locations that were revisited

end
