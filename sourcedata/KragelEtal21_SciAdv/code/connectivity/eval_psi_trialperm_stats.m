function eval_psi_trialperm_stats()
% EVAL_PSI_TRIALPERM_STATS reruns main psi analysis using permutation at the fixation
% level to determine significance
%
% INPUTS:  no inputs - loads relevant data from temporary data directories
%
% OUTPUTS: no outputs - runs stats and generates figures

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

subjects = {'S1','S2','S4','S5','S6'};

for s = 1:length(subjects)
    res(s) = gather_psi_stats(subjects{s});
end

% init some structures
all_conn = nan(5,751);
all_dan = all_conn;
all_vis = all_conn;

rev_conn = nan(5,751);
rev_dan = rev_conn;
rev_vis = rev_conn;

other_conn = nan(5,751);
other_dan = other_conn;
other_vis = other_conn;

t_idx = 76:676;

for s = 1:length(subjects)
       
    % correct rois
    [rois, is_hc, is_vis, is_dan, chansel, pairs] = get_plv_rois(subjects{s});
    
    dan_mask = is_dan(chansel(:,1))' & ~is_hc(chansel(:,1))';
    vis_mask = is_vis(chansel(:,1))' & ~is_hc(chansel(:,1))';
    
    
    [vis_psi_diff(s,:), vis_psi_rev(s,:), vis_psi_oth(s,:), ...
        vis_null_psi_diff(s,:,:), vis_null_psi_rev(s,:,:), vis_null_psi_oth(s,:,:)] = gen_null(res(s), vis_mask);

    [dan_psi_diff(s,:), dan_psi_rev(s,:), dan_psi_oth(s,:), ...
        dan_null_psi_diff(s,:,:), dan_null_psi_rev(s,:,:), dan_null_psi_oth(s,:,:)] = gen_null(res(s), dan_mask);
    
  
end


vis_idx =  t_idx(round(median(find(stats_vis.mask))));

[~,~,~,st] = ttest(vis_psi_diff(:, vis_idx));
vis_diff = nanmean(st.tstat);

[~,~,~,st] = ttest(vis_null_psi_diff(:, vis_idx,:));
vis_null_diff = st.tstat;

n_pos = sum(vis_diff <= squeeze(vis_null_diff)');
p_pos = (n_pos+1)/1001;
n_neg = sum(vis_diff >= squeeze(vis_null_diff)');
p_neg = (n_neg+1)/1001;
p_vis = min(p_pos, p_neg);

% dan revisit vs. all trials

dan_idx = t_idx(round(median(find(stats_rev_dan.mask))));

[~,~,~,st] = ttest(dan_psi_rev(:, dan_idx));
dan_rev = nanmean(st.tstat);

[~,~,~,st] = ttest(dan_null_psi_rev(:, dan_idx, :));
dan_null_rev = st.tstat;

n_pos = sum(dan_rev <= squeeze(dan_null_rev)');
p_pos = (n_pos+1)/1001;
n_neg = sum(dan_rev >= squeeze(dan_null_rev)');
p_neg = (n_neg+1)/1001;

p_dan_rev = min(p_pos, p_neg);

% and dan, other vs. all

dan_idx = t_idx(round(median(find(stats_other_dan.mask))));

[~,~,~,st] = ttest(dan_psi_oth(:, dan_idx));
dan_oth = nanmean(st.tstat);

[~,~,~,st] = ttest(dan_null_psi_oth(:, dan_idx, :));
dan_null_oth = st.tstat;

n_pos = sum(dan_oth <= squeeze(dan_null_oth)');
p_pos = (n_pos+1)/1001;
n_neg = sum(dan_oth >= squeeze(dan_null_oth)');
p_neg = (n_neg+1)/1001;

p_dan_oth = min(p_pos, p_neg);

% vis, other vs. all
 
vis_idx = t_idx(round(median(find(stats_other_vis.mask))));

[~,~,~,st] = ttest(vis_psi_oth(:, vis_idx));
vis_oth = nanmean(st.tstat);

[~,~,~,st] = ttest(vis_null_psi_oth(:, vis_idx, :));
vis_null_oth = st.tstat;

n_pos = sum(vis_oth <= squeeze(vis_null_oth)');
p_pos = (n_pos+1)/1001;
n_neg = sum(vis_oth >= squeeze(vis_null_oth)');
p_neg = (n_neg+1)/1001;

p_vis_oth = min(p_pos, p_neg);

stats_vis = group_perm_stats(vis_psi_rev, vis_psi_oth, t_idx);
stats_dan = group_perm_stats(rev_dan, other_dan, t_idx);

% mes stats for vis and dan

idx = get_peak(stats_vis.stat(stats_vis.mask));
m_idx = find(stats_vis.mask);

c1 = nanmean(rev_vis(:, t_idx(m_idx(idx))),2);
c2 = nanmean(other_vis(:, t_idx(m_idx(idx))),2);
mes_stats_dv = mes(c1, c2, 'hedgesg', 'isdep', true);

% looking averaged across conditions
stats_other_vis = group_perm_stats(other_vis, zeros(size(all_vis)), t_idx);
stats_other_dan = group_perm_stats(other_dan, zeros(size(all_vis)), t_idx);

idx = get_peak(-stats_other_vis.stat(stats_other_vis.mask));
m_idx = find(stats_other_vis.mask);
c1 = nanmean(other_vis(:, t_idx(m_idx(idx))),2);

idx = get_peak(stats_other_dan.stat(stats_other_dan.mask));
m_idx = find(stats_other_dan.mask);
c2 = nanmean(other_dan(:, t_idx(m_idx(idx))),2);

mes_stats_ov = mes(c1, 0, 'g1');
mes_stats_od = mes(c2, 0, 'g1');

stats_rev_vis = group_perm_stats(rev_vis, zeros(size(all_vis)), t_idx);
stats_rev_dan = group_perm_stats(rev_dan, zeros(size(all_vis)), t_idx);

idx = get_peak(stats_rev_dan.stat(stats_rev_dan.mask));
m_idx = find(stats_rev_dan.mask);
c2 = nanmean(rev_dan(:, t_idx(m_idx(idx))),2);

% mes_stats_rv = mes(c1, 0, 'g1'); - nonsig
mes_stats_rd = mes(c2, 0, 'g1');

end

function stats = group_perm_stats(c1, c2, tidx)

if ~exist('tidx','var')
    tidx = 1:size(c1,2);
%     tidx = 1:751; % -.75 to .10 sec
end

dat.dimord = 'chan_rpt_time'; % chan is really subject but doesnt matt
dat.label  = {'xx'};
dat.time   = -.75:.002:.75;
dat.time   = dat.time(tidx);
dat.avg    = cat(1,c1,c2);
dat.avg    = reshape(dat.avg, 1, size(dat.avg,1), size(dat.avg,2));
dat.avg    = dat.avg(:,:,tidx);

c1_idx = false(1, size(dat.avg,2));
c1_idx(1:size(c1,1)) = true;

stats = run_group_psi_stats(dat, c1_idx, ~c1_idx);

end

function res = gather_psi_stats(subj)

[psi_revisit, psi_other, psi_revisit_null, psi_other_null, anat, label] = get_psi_trialperm_by_network(subj);

res.label = label;

% some book-keeping
res.psi_diff = psi_revisit - psi_other;
res.psi_revisit = psi_revisit;
res.psi_other   = psi_other;

res.psi_revisit_null = psi_revisit_null;
res.psi_other_null = psi_other_null;

res.dan_con = anat.dan_con;
res.vis_con = anat.vis_con;

res.rev_vis = nanmean(psi_revisit(anat.vis_con,:),1);
res.oth_vis = nanmean(psi_other(anat.vis_con,:),1);

res.rev_dan = nanmean(psi_revisit(anat.dan_con,:),1);
res.oth_dan = nanmean(psi_other(anat.dan_con,:),1);

res.rev_null_vis = nanmean(psi_revisit_null(anat.vis_con,:,:),1);
res.oth_null_vis = nanmean(psi_other_null(anat.vis_con,:,:),1);

res.rev_null_dan = nanmean(psi_revisit_null(anat.dan_con,:,:),1);
res.oth_null_dan = nanmean(psi_other_null(anat.dan_con,:,:),1);

end

function [m_psi_diff, m_psi_rev, m_psi_oth, null_psi_diff, null_psi_rev, null_psi_oth] = gen_null(res, c_idx)

psi_diff = res(1).psi_revisit(c_idx,:) - res(1).psi_other(c_idx,:);

% diff
null_m = squeeze(mean(res(1).psi_revisit_null(c_idx,:,:) - res(1).psi_other_null(c_idx,:,:),3));
null_sd = squeeze(std(res(1).psi_revisit_null(c_idx,:,:) - res(1).psi_other_null(c_idx,:,:),[],3));
psi_diff_z = (psi_diff-null_m)./null_sd;
psi_diff_null_z = (res(1).psi_revisit_null(c_idx,:,:) - res(1).psi_other_null(c_idx,:,:) - null_m)./null_sd;

% rev
null_m = squeeze(mean(res(1).psi_revisit_null(c_idx,:,:), 3));
null_sd = squeeze(std(res(1).psi_revisit_null(c_idx,:,:),[],3));
psi_rev_z = (res(1).psi_revisit(c_idx,:)-null_m)./null_sd;
psi_rev_null_z = (res(1).psi_revisit_null(c_idx,:,:)-null_m)./null_sd;

% other
null_m = squeeze(mean(res(1).psi_other_null(c_idx,:,:), 3));
null_sd = squeeze(std(res(1).psi_other_null(c_idx,:,:),[],3));
psi_oth_z = (res(1).psi_other(c_idx,:)-null_m)./null_sd;
psi_oth_null_z = (res(1).psi_other_null(c_idx,:,:)-null_m)./null_sd;

% only looks at channels with some directional effect
p_rev = 2*normcdf(-abs(psi_rev_z));
p_oth = 2*normcdf(-abs(psi_oth_z));

p_rev_null = 2*normcdf(-abs(psi_rev_null_z));
p_oth_null = 2*normcdf(-abs(psi_oth_null_z));

t_idx = 76:676;

% observed - difference
sig_idx = any(p_rev(:,t_idx) < .05 | p_oth(:,t_idx) < .05, 2);
m_psi_diff = nanmean(psi_diff_z(sig_idx==1,:),1);
m_psi_rev = nanmean(psi_rev_z(sig_idx==1,:),1);
m_psi_oth = nanmean(psi_oth_z(sig_idx==1,:),1);


% and null
null_psi_diff = nan(751, 1000);
null_psi_rev = nan(751, 1000);
null_psi_oth = nan(751, 1000);

for i =  1:size(psi_rev_null_z,3)
    null_psi_diff(:,i) = nanmean(psi_diff_null_z(sig_idx==1,:,i),1);
    null_psi_rev(:,i) = nanmean(psi_rev_null_z(sig_idx==1,:,i),1);
    null_psi_oth(:,i) = nanmean(psi_oth_null_z(sig_idx==1,:,i),1);
end

end

function idx = get_peak(vec)

p = prctile(vec, 100);
idx = vec==p;

end
