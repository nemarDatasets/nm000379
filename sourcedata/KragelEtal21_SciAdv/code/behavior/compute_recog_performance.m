function  stats = compute_recog_performance(events)
% COMPUTE_RECOG_PERFORMANCE computes sensitivity (d') and rt as a function
% of condition
%
% INPUTS:  events - structure, contains event information
%
% OUTPUTS: stats - structure, results structure containing recognition
%                  performance metrics for this subject

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

% compute d-prime for low and high salience cues
recog_trials = events(strcmp({events.type},'trial') & ...
    strcmp({events.phase},'recog'));

old = [recog_trials.old];
sal = [recog_trials.salience];
resp = [recog_trials.resp];
resp(resp==4)=0; resp(resp==2)=nan;

[stats.C(1,:,:), stats.d_prime(1), stats.bias(1)] = recog_stats(old(sal==0), resp(sal==0));
[stats.C(2,:,:), stats.d_prime(2), stats.bias(2)] = recog_stats(old(sal==1), resp(sal==1));

[stats.all_C, stats.all_d_prime, stats.all_bias] = recog_stats(old(sal==0), resp(sal==0));

rt = [recog_trials.rt]/1000;

% ranksum
stats.rt = run_ranksum(rt, old, sal);

stats.rt.dist_old_high = rt(old==1 & sal==1);
stats.rt.dist_new_high = rt(old==0 & sal==1);
stats.rt.dist_old_low = rt(old==1 & sal==0);
stats.rt.dist_new_low = rt(old==0 & sal==0);

% compute the number of fixations for low and high salience cues
recog_fixations = events(strcmp({events.type},'fixation') & ...
    strcmp({events.phase},'recog'));

old = [recog_fixations.old];
sal = [recog_fixations.salience];
resp = [recog_fixations.resp];

n_fix = nan(1,length(recog_fixations));
for t = 1:length(recog_fixations)
    n_fix(t) = length(recog_fixations(t).start_time);
end

stats.n_fix = run_ranksum(n_fix, old, sal);

stats.n_fix.n_old_high = n_fix(old==1 & sal==1);
stats.n_fix.n_new_high = n_fix(old==0 & sal==1);
stats.n_fix.n_old_low = n_fix(old==1 & sal==0);
stats.n_fix.n_new_low = n_fix(old==0 & sal==0);

stats.n_fix.resp_old_high = resp(old==1 & sal==1);
stats.n_fix.resp_new_high = resp(old==0 & sal==1);
stats.n_fix.resp_old_low = resp(old==1 & sal==0);
stats.n_fix.resp_new_low = resp(old==0 & sal==0);

% run analysis of n_fix at encoding
enc_fixations = events(strcmp({events.type},'fixation') & ...
    strcmp({events.phase},'encode'));

hit = [enc_fixations.resp] == 1;
sal = [enc_fixations.salience];

n_enc_fix = nan(1,length(enc_fixations));
for t = 1:length(enc_fixations)
    n_enc_fix(t) = length(enc_fixations(t).start_time);
end

% below is for SME, doesn't show much
stats.n_enc_fix = run_ranksum(n_enc_fix, hit, sal);

stats.n_enc_fix.n_hit_high = n_enc_fix(hit==1 & sal==1);
stats.n_enc_fix.n_miss_high = n_enc_fix(hit==0 & sal==1);
stats.n_enc_fix.n_hit_low = n_enc_fix(hit==1 & sal==0);
stats.n_enc_fix.n_miss_low = n_enc_fix(hit==0 & sal==0);

end

function [C, d_prime, bias] = recog_stats(old, resp)
% RECOG_STATS
% Stanislaw and Todorov 1999 Behav Res Meth

C = confusionmat(old, resp);
trial_cnt = sum(C,2);
C(1,:) = C(1,:)./trial_cnt(1);
C(2,:) = C(2,:)./trial_cnt(2);

for x = 1:2
    for y = 1:2
        if C(x,y) == 1
            C(x,y) = 1-.5/trial_cnt(x);
        elseif C(x,y) == 0
            C(x,y) = .5/trial_cnt(x);
        end
    end
end

pHit = C(2,2);
pFA = C(1,2);

zHit = -sqrt(2)*erfcinv(2*pHit);
zFA  = -sqrt(2)*erfcinv(2*pFA);
d_prime = zHit - zFA;
bias = -.5*(zHit + zFA);

end

function res = run_ranksum(stat, old, sal)

[p, ~, st] = ranksum(stat(old==1), stat(old==0));

res.p_oldvnew = p;
res.z_oldvnew = st.zval;

[p, ~, st] = ranksum(stat(sal==1), stat(sal==0));
res.p_highvlow = p;
res.z_highvlow = st.zval;

[p, ~, st] = ranksum(stat(old==1 & sal==1), stat(old==0 & sal==1), 'method', 'approximate');
res.p_high_oldvnew = p;
res.z_high_oldvnew = st.zval;

[p, ~, st] = ranksum(stat(old==1 & sal==0), stat(old==0 & sal==0), 'method', 'approximate');
res.p_low_oldvnew = p;
res.z_low_oldvnew = st.zval;

[p, ~, st] = ranksum(stat(old==1 & sal==1), stat(old==1 & sal==0), 'method', 'approximate');
res.p_old_highvlow = p;
res.z_old_highvlow = st.zval;

[p, ~, st] = ranksum(stat(old==0 & sal==1), stat(old==0 & sal==0), 'method', 'approximate');
res.p_new_highvlow = p;
res.z_new_highvlow = st.zval;


end
