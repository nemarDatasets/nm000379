function res = spotlight_cvp(events, close_only, exclude_precue)
% SPOTLIGHT CVP computes conditional viewing probabilities for
% reinstatement analysis
%
% INPUTS:  events - structure, events for a subject
%
%          close_only - logical, whether to exclude novel locations
%                       when considering transitions at test
%
%          exclude_precue - logical, whether to exclude fixations in
%                           the period before the scene appears 
%
% OUTPUTS: res - structure, contains results about reinstatement analysis

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

if ~exist('exclude_precue','var')
    exclude_precue = false;
end

if ~exist('close_only','var')
    close_only = false;
end

% from the events structure, grab encoding and retrieval fixations
fix_enc = events(strcmp({events.type}, 'fixation') & ...
    strcmp({events.phase},'encode'));

fix_rec = events(strcmp({events.type}, 'fixation') & ...
    strcmp({events.phase},'recog'));

% set some initial parameters
skip_only = false;
exclude_first = false;
thr = 80;

% maximum number of lags to evaluate in the CVP
max_n = 5;

% initialize for memory allocation
cvp = nan(length(fix_enc), ...
    2*max_n);

for i = 1:length(fix_enc) % image
    
    % where is this image in the recall sequence
    r_idx = [fix_rec.stim_idx] == fix_enc(i).stim_idx & ...
        [fix_rec.block] == fix_enc(i).block;
    
    if exclude_precue
        enc_precue = fix_enc(i).precue;
        rec_precue = fix_rec(r_idx).precue;
    else
        enc_precue = false(size(fix_enc(i).precue));
        rec_precue = false(size(fix_rec(r_idx).precue));
    end
        
        
    rec_mask = true(size(rec_precue));
    
    enc_seq.fix = [fix_enc(i).mean_x(~enc_precue); fix_enc(i).mean_y(~enc_precue)];
    rec_seq.fix = [fix_rec(r_idx).mean_x(~rec_precue & rec_mask); fix_rec(r_idx).mean_y(~rec_precue & rec_mask)];
    
    n_fix_enc(i) = size(enc_seq.fix,2);
    n_fix_rec(r_idx) = size(rec_seq.fix,2);
    
    p_rep(i,:) = compute_prob_repeat(enc_seq, rec_seq, thr);

    
    [cvp(i,:), enc_reins(i,:), actual(i,:), possible(i,:), lag_one(i,:), lag_minusone(i,:), recalls(i,:), lags(i,:)] = compute_cvp(enc_seq, ...
        rec_seq, ...
        thr, ...
        max_n, ...
        skip_only, ...
        exclude_first, ...
        close_only);    
    
end

cvp_h = cvp([fix_enc.resp]==1,:);
cvp_m = cvp([fix_enc.resp]~=1,:);

res.cvp = cvp;
res.cvp_h = cvp_h;
res.cvp_m = cvp_m;

% output information to update encoding events with 'subsequent
% reinstatement'
res.enc_reins = enc_reins;
res.lag_one   = lag_one;
res.lag_minusone = lag_minusone;
res.recalls   = recalls;
res.lags      = lags;


end

function p_rep = compute_prob_repeat(s1, s2, thr)
% compare proportion of repeats for novel vs. repeated items

% p_rep is the probability of repeatedly viewing old vs. new locations

if isempty(s1) || isempty(s2) || size(s2.fix,2)==1
    p_rep = [nan nan];
    return
end

distmat = dist(s2.fix);
rmat = distmat < thr;
repeat = [nan; diag(rmat,-1)];

% compute 'serial position' of recalls
distmat = dist([s1.fix s2.fix]);
recmat = distmat < thr;
tmp = recmat(size(s1.fix,2)+1:end, 1:size(s1.fix,2));

tmp_dist = distmat(size(s1.fix,2)+1:end, 1:size(s1.fix,2));
recalls = nan(1, size(tmp_dist,1));

for e = 1:size(tmp_dist,1)
    this_dist = tmp_dist(e,:);
    idx = find(this_dist==min(this_dist) & tmp(e,:));
    if ~isempty(idx)
        recalls(e) = idx;
    else
        recalls(e) = 0; % new
    end
end

p_rep(1) = nanmean(repeat(recalls~=0));% old
p_rep(2) = nanmean(repeat(recalls==0));% new

end

function [cvp, enc_reins, actual, possible, lag_one, lag_minusone, recalls_out, lags_out] = compute_cvp(s1, s2, thr, ...
    max_n, skip_only, exclude_first, close_only)

if isempty(s1) || isempty(s1.fix) || isempty(s2)
    cvp = nan(size(-max_n+1:1:max_n-1));
    cvp(end+1) = nan;
    actual = cvp;
    possible = cvp;
    enc_reins = nan(1,200);
    lag_one = nan(1,200);
    lag_minusone = nan(1,200);
    recalls_out = nan(1,200);
    lags_out = nan(1,200);
    return
end

distmat = dist([s1.fix s2.fix]);
recmat = distmat < thr;
tmp = recmat(size(s1.fix,2)+1:end, 1:size(s1.fix,2));

% compute 'serial position' of recalls
tmp_dist = distmat(size(s1.fix,2)+1:end, 1:size(s1.fix,2));
recalls = nan(1,size(tmp_dist,1));
for e = 1:size(tmp_dist,1)
    this_dist = tmp_dist(e,:);
    idx = find(this_dist==min(this_dist) & tmp(e,:));
    if ~isempty(idx)
        recalls(e) = idx;
    else
        recalls(e) = 0; % new
    end
end

recalls(recalls==0)=nan;

mask = true(size(recalls));
lags = [diff(recalls,[],2) nan];

if skip_only
    
    % first jump forward
    [~,j] = find(lags > 1, 1, 'first');
    mask = false(size(mask));
    
    % include skipped item and the next item to compute transitions
    mask(j+1) = true;
    mask(j+2) = true;
    
end


% output some encoding information
reinstated = unique(recalls(~isnan(lags)));
enc_reins = nan(1,200);
if ~isempty(reinstated)
    enc_reins(reinstated) = 1;
    enc_reins(enc_reins~=1) = 0;
end

lag_one = nan(1,200);
lag_minusone = nan(1,200);
f_reinstated = unique(recalls(lags==1));
r_reinstated = unique(recalls(lags==-1));

if ~isempty(f_reinstated)
    lag_one(f_reinstated) = 1;
    lag_one(lag_one~=1) = 0;
end

if ~isempty(r_reinstated)
    lag_minusone(r_reinstated) = 1;
    lag_minusone(lag_minusone~=1) = 0;
end


if close_only
    idx = ~isnan(recalls);
    recalls = recalls(idx); 
    lags = [diff(recalls,[],2) nan];
end

% output some retrieval information
recalls_out = zeros(1,200);
if ~close_only
    recalls_out(1:length(recalls)) = recalls;
else
    recalls_out(idx) = recalls;
end

lags_out = nan(1,200);
if ~close_only
        lags_out(1:length(lags)) = lags;
else
    lags_out(idx) = lags;
end

to_mask_pres = true(size(s1.fix(1,:)));
from_mask_pres =  true(size(s1.fix(1,:)));

from_mask_rec = ~isnan(recalls);
to_mask_rec = ~isnan(recalls);

params.to_mask_pres = to_mask_pres;
params.from_mask_pres = from_mask_pres;

[actual, possible] = ...
    conditional_transitions(recalls, ...
                            from_mask_rec, ...
                            to_mask_rec, ...
                            @lag, ...
                            @possible_transitions, ...
                            1, params);
                        
possible = [possible{:}];
                        

all_possible_transitions = -max_n + 1 : max_n - 1;
n_actual = collect(actual, all_possible_transitions);
n_possible = collect(possible, all_possible_transitions);

% compute prob of novel
all_poss = unique(possible);
n_act = collect(actual, all_poss);
n_poss = collect(possible, all_poss);
prob_new = 1-sum(n_act)/sum(n_poss);

cvp = n_actual ./ n_possible;
cvp = [cvp prob_new];
actual = [n_actual nan];
possible = [n_possible nan];

end

function d = lag(sp1, sp2, params)
  d = sp2 - sp1;
end