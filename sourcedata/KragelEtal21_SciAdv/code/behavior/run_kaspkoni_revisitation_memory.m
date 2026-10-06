function revisit_res = run_kaspkoni_revisitation_memory(res, null_res, dataset)
% RUN_KASPKONI_REVISITATION_MEMORY organizes results for statistical comparison
%
% INPUTS:  res - structure, figrim results structure with fixation-sequence
%                replay information
%
%          null_res - structure, figrim results structure with fixation-sequence
%                replay information based on similar viewing across subjects
%
%          dataset - string, specifying the name of the dataset
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


for s = 1:size(res.(dataset).n_fix,1)
    for f = 1:size(res.(dataset).n_fix,2)
        for i = 1:size(res.(dataset).n_fix,3)
            
            % observed
            n_fix = res.(dataset).n_fix(s,f,i);
            
            viewed = res.(dataset).sub_view(s,f,i,1:n_fix);
            reinstated = res.(dataset).sub_reins(s,f,i,1:n_fix);
            rev_reinstated = res.(dataset).sub_rev_reins(s,f,i,1:n_fix);            
            
            enc_seq = res.(dataset).enc_seq(s,f,i,1:n_fix);
            replay = res.(dataset).enc_replay(s,f,i,1:n_fix);
            repeated = res.(dataset).enc_repeated(s,f,i,1:n_fix);
            rev_replay = res.(dataset).enc_rev_replay(s,f,i,1:n_fix); % not looking at this
            
                        
            m_viewed(s,f,i,1) = get_trial_stats(enc_seq, repeated==1, viewed, false);
            m_viewed(s,f,i,2) = get_trial_stats(enc_seq, isnan(repeated), viewed, false);
            
            m_reins(s,f,i,1) = get_trial_stats(enc_seq, repeated==1, reinstated, false);
            m_reins(s,f,i,2) = get_trial_stats(enc_seq, isnan(repeated), reinstated, false);
            
            m_rev_reins(s,f,i,1) = get_trial_stats(enc_seq, repeated==1, rev_reinstated, false);
            m_rev_reins(s,f,i,2) = get_trial_stats(enc_seq, isnan(repeated), rev_reinstated, false);
            
            % other subjects
            n_fix = res.(dataset).n_fix(s,f,i);
            
            tmp_viewed = nan(n_fix, 2);
            tmp_reins = nan(n_fix, 2);
            tmp_rev_reins = nan(n_fix, 2);
            
            for np = 1:size(null_res.(dataset).enc_seq(s,f,i).dat,1)

            n_fix = null_res.(dataset).n_fix(s,f,i).dat(np);
            viewed = null_res.(dataset).sub_view(s,f,i).dat(np, 1:n_fix);
            reinstated = null_res.(dataset).sub_reins(s,f,i).dat(np,1:n_fix);
            rev_reinstated = null_res.(dataset).sub_rev_reins(s,f,i).dat(np, 1:n_fix);
                
            enc_seq = null_res.(dataset).enc_seq(s,f,i).dat(np, 1:n_fix);
            replay = null_res.(dataset).enc_replay(s,f,i).dat(np, 1:n_fix);
            repeated = null_res.(dataset).enc_repeated(s,f,i).dat(np, 1:n_fix);
            
            tmp_viewed(np, 1) = get_trial_stats(enc_seq, repeated==1, viewed, false);
            tmp_viewed(np, 2) = get_trial_stats(enc_seq, isnan(repeated), viewed, false);
            
            tmp_reins(np, 1) = get_trial_stats(enc_seq, repeated==1, reinstated, false);
            tmp_reins(np, 2) = get_trial_stats(enc_seq, isnan(repeated), reinstated, false);
            
            tmp_rev_reins(np, 1) = get_trial_stats(enc_seq, repeated==1, rev_reinstated, false);
            tmp_rev_reins(np, 2) = get_trial_stats(enc_seq, isnan(repeated), rev_reinstated, false);
            
            end
            
            null_viewed(s,f,i,:) = nanmean(tmp_viewed);            
            null_reins(s,f,i,:) = nanmean(tmp_reins);            
            null_rev_reins(s,f,i,:) = nanmean(tmp_rev_reins);            

            diff_viewed(s,f,i,:) = squeeze(m_viewed(s,f,i,:))' - nanmean(tmp_viewed);            
            diff_reins(s,f,i,:) = squeeze(m_reins(s,f,i,:))' - nanmean(tmp_reins);            
            diff_rev_reins(s,f,i,:) = squeeze(m_rev_reins(s,f,i,:))' - nanmean(tmp_rev_reins);
            
        end
    end
end

m_viewed = squeeze(nanmean(nanmean(m_viewed,2),3));
m_reins = squeeze(nanmean(nanmean(m_reins,2),3));
m_rev_reins = squeeze(nanmean(nanmean(m_rev_reins,2),3));

revisit_res.m_viewed = m_viewed;
revisit_res.m_reins = m_reins;
revisit_res.m_rev_reins = m_rev_reins;

% null/stim
null_viewed = squeeze(nanmean(nanmean(null_viewed,2),3));
null_reins = squeeze(nanmean(nanmean(null_reins,2),3));
null_rev_reins = squeeze(nanmean(nanmean(null_rev_reins,2),3));

revisit_res.null_reins = null_reins;
revisit_res.null_rev_reins = null_rev_reins;

% difference
diff_viewed = squeeze(nanmean(nanmean(diff_viewed,2),3));
diff_reins = squeeze(nanmean(nanmean(diff_reins,2),3));
diff_rev_reins = squeeze(nanmean(nanmean(diff_rev_reins,2),3));

end

function [m, n] = get_trial_stats(enc_seq, idx, dat, binary_flag)

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
