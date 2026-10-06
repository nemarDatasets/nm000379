function [res, p_enc_rep2, p_rec_rep2] = spotlight_enc_cvp(events, to_repeat, exclude_precue)
% SPOTLIGHT CVP computes conditional viewing probabilities for
% reinstatement within the study period - finds revisitation fixations
% as well
%
% INPUTS:  events - structure, events for a subject
%
%          to_repeat - logical, whether to only consider fixations to
%                      revisited locations
%
%          exclude_precue - logical, whether to exclude fixations in
%                           the period before the scene appears 
%
% OUTPUTS: res - structure, contains results about reinstatement analysis
%
%          p_enc_rep2 - double, proportion of revisitations at study
%
%          p_rec_rep2 - double, proportion of revisitations at test

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

if ~exist('to_repeat','var')
    to_repeat = false;
end

% from the events structure, grab encoding fixations

fix_enc = events(strcmp({events.type}, 'fixation') & ...
    strcmp({events.phase},'encode'));


fix_rec = events(strcmp({events.type}, 'fixation') & ...
    strcmp({events.phase},'recog'));

% set some initial parameters
thr = 120;
max_n = 5;

% initialize for memory allocation
cvp = nan(length(fix_enc), ...
    2*max_n);

for i = 1:length(fix_enc) % image
    
    % where is this image in the recall sequence
    r_idx = [fix_rec.stim_idx] == fix_enc(i).stim_idx & ...
        [fix_rec.block] == fix_enc(i).block;
    
    if exclude_precue
        precue = fix_enc(i).precue;
        rec_precue = fix_rec(r_idx).precue;
    else
        precue = false(size(fix_enc(i).precue));
        rec_precue = false(size(fix_rec(r_idx).precue));
    end
    
    rec_mask = true(size(rec_precue));
    enc_seq.fix = [fix_enc(i).mean_x(~precue); fix_enc(i).mean_y(~precue)];
    rec_seq.fix = [fix_rec(r_idx).mean_x(~rec_precue & rec_mask); fix_rec(r_idx).mean_y(~rec_precue & rec_mask)];
    
    [cvp(i,:), enc_reins(i,:), actual(i,:), possible(i,:), lag_one(i,:), lag_minusone(i,:), recalls(i,:), lags(i,:)] = compute_cvp(enc_seq, ...
        enc_seq, ...
        thr, ...
        max_n, ...
        to_repeat);
    
    [~, rec_reins(i,:), ~, ~, ~, ~, rs(i,:), rec_lags(i,:)] = compute_cvp(rec_seq, ...
        rec_seq, ...
        thr, ...
        max_n, ...
        to_repeat);
    
    
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

enc_rep = []; enc_lags = [];
for i = 1:size(enc_reins,1)
    n_fix = sum(recalls(i,:)~=0);
       
    enc_rep = [enc_rep enc_reins(i,1:n_fix)];
    enc_lags = [enc_lags lags(i,1:n_fix)];
end

rec_rep = []; r_lags= [];
for i = 1:size(rec_reins,1)
    n_fix = sum(rs(i,:)~=0);
    r_lags = [r_lags [nan rec_lags(i,2:n_fix-1)]];
    rec_rep = [rec_rep  [nan rec_reins(i,2:n_fix-1)]];
end

p_enc_rep = sum(enc_rep==0)/length(enc_rep);
p_rec_rep = sum(rec_rep==0)/length(rec_rep);

p_enc_rep2 = sum(enc_rep(enc_lags~=0)==0)/length(enc_rep(enc_lags~=0));
p_rec_rep2 = sum(~(r_lags==1 | r_lags==-1 | r_lags==0) & rec_rep==0)/ sum(~(r_lags==1 | r_lags==-1 | r_lags==0));

end

function [cvp, enc_reins, actual, possible, lag_one, lag_minusone, recalls_out, lags_out] = compute_cvp(s1, s2, thr, ...
    max_n, to_repeat)

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
novel_mask = nan(1, size(tmp_dist,1));

scnt = 1;
for e = 1:size(tmp_dist,1)
    this_dist = tmp_dist(e,:); this_dist(e) = nan;
    idx = find(this_dist==min(this_dist) & tmp(e,:));
    if ~isempty(idx)
        if (e > idx)
            recalls(e) = recalls(idx);
            novel_mask(e) = false;
        else
            recalls(e) = scnt;
            novel_mask(e) = true;
            scnt=scnt+1;
        end
    else
        recalls(e) = scnt; % new
        novel_mask(e) = true;
        scnt=scnt+1;
    end
end

recalls(recalls==0)=nan;

lags = [diff(recalls,[],2) nan];

% output some encoding information
reinstated = unique(recalls(~novel_mask));
enc_reins = nan(1,200);
if ~isempty(reinstated)
    for e = 1:length(reinstated)
        enc_reins(recalls == reinstated(e) & novel_mask) = 1; % later viewed in this sequence
        enc_reins(recalls == reinstated(e) & ~novel_mask) = 0;
    end
end

lag_one = nan(1,200);
to_mask = [nan novel_mask(1:end-1)];
reinstated = unique(recalls(~novel_mask & lags ==1 & to_mask == 0));
if ~isempty(reinstated)
    for e = 1:length(reinstated)
        lag_one(recalls == reinstated(e) & novel_mask) = 1;
        lag_one(recalls == reinstated(e) & ~novel_mask) = 0;
    end
end

lag_minusone = nan(1,200);
to_mask = [nan novel_mask(1:end-1)];
reinstated = unique(recalls(~novel_mask & lags ==-1 & to_mask == 0));
if ~isempty(reinstated)
    for e = 1:length(reinstated)
        lag_minusone(recalls == reinstated(e) & novel_mask) = 1;
        lag_minusone(recalls == reinstated(e) & ~novel_mask) = 0;
    end
end

% output some retrieval information
recalls_out = zeros(1,200);
recalls_out(1:length(recalls)) = recalls;


lags_out = nan(1,200);
lags_out(1:length(lags)) = lags;

to_mask_pres = true(size(1:scnt));
from_mask_pres =  true(size(1:scnt));

if to_repeat
    from_mask_rec = novel_mask == 1; %& lags~= 0
    to_mask_rec = novel_mask == 0;
else % from repeated items to repeated
    from_mask_rec = novel_mask == 0; %& lags~= 0;
    to_mask_rec = novel_mask == 0; %true(size(novel_mask));
end

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