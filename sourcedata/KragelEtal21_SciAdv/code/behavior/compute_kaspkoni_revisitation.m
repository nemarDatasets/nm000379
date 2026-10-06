function [res, null_res] = compute_kaspkoni_revisitation()
% COMPUTE_KASPKONI_REVISITATION computes fixation-sequence replay for revisitation
% and other fixations
%
% INPUTS:  none, loads data from disk and runs analysis
%
% OUTPUTS: res - structure, contains behavioral results
%
%          null_res - structure, contains results for stimulus-driven viewing
%                     behaviors obtained by permuting fixation sequences 
%                     across subjects viewing the same scene

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

run_perm = false; % for running stimulus-driven as well

% dataset 1

fixmat = get_fixmat('etdb_v1.0.hdf5', 'Memory I');
[fix_all, t_fix, subs, filenos, iter] = convert_data(fixmat);
subs = unique(fixmat.SUBJECTINDEX);

for f = 1:length(subs)
    filter_struct.SUBJECTINDEX = subs(f);
    
    [fix_seq, idx] = prep_seq(fixmat, filter_struct);
    
    iter = fixmat.iteration(idx);
    category  = fixmat.category(idx);
    filenumber = fixmat.filenumber(idx);
    
    
    fcnt = 1;
    for s = 1:length(fix_seq)
        nf = length(fix_seq(s).dur);
        % include repetition
        fix_seq(s).iter = iter(fcnt:fcnt+nf-1);
        fix_seq(s).iter = fix_seq(s).iter';
        
        % include stimulus category
        fix_seq(s).category = category(fcnt:fcnt+nf-1);
        fix_seq(s).category = fix_seq(s).category';
        
        % and fileno
        fix_seq(s).fileno = filenumber(fcnt:fcnt+nf-1);
        fix_seq(s).fileno = fix_seq(s).fileno';
        
        fcnt = fcnt+nf;
    end
    
    iter = nan(1,length(fix_seq));
    category = nan(1,length(fix_seq));
    fileno = nan(1,length(fix_seq));
    for s = 1:length(fix_seq)
        iter(s) = mode(fix_seq(s).iter);
        category(s) = mode(fix_seq(s).category);
        fileno(s) = mode(fix_seq(s).fileno);
    end
    
    % compute cvp, based off of previous viewing
    
    [id, ~, ib] = unique([category; fileno]','rows');
    
    thr = 84;
    max_n = 10;
    
    combs = nchoosek(1:5, 2);
    iter_lags = diff(combs, [], 2);
    
    for i = 1:size(combs,1)
        for stim = 1:size(id,1)
            
            idx_1 = category == id(stim,1) & ...
                fileno  == id(stim,2) & ...
                iter == combs(i,1);
            
            % adding permutation (subject) nulls for stimulus-driven viewing
            
            null_idx_1 = find([fix_all.SUBJECTINDEX] ~= subs(f) & ...
                [fix_all.category] == id(stim,1) & ...
                [fix_all.filenumber] == id(stim,2) & ...
                [fix_all.iteration] == combs(i,1));
            
            idx_2 = category == id(stim,1) & ...
                fileno  == id(stim,2) & ...
                iter == combs(i,2);
        
            if ~any(idx_1) || ~any(idx_2)
                continue
            end
            
            % use idx_1 and idx_2 to examine whether fixations are
            % reinstated or replayed during the second viewing
            
            [res.mem1.cvp(f,stim,i,:), enc_reins, ~, ~, lag_one, lag_minusone] = compute_cvp(fix_seq(idx_1), ...
                fix_seq(idx_2), ...
                thr, ...
                max_n, ...
                false, ...
                false);           
            
            res.mem1.n_fix(f,stim,i) = size(fix_seq(idx_1).fix,2);
            res.mem1.sub_view(f,stim,i,:) = enc_reins;
            res.mem1.sub_reins(f,stim,i,:) = lag_one;
            res.mem1.sub_rev_reins(f,stim,i,:) = lag_minusone;
   
            % reinstatement/replay events within the encoding period
            [~, enc_reins, ~, ~, lag_one, lag_minusone, ...
                recalls_out, ~] = compute_cvp(fix_seq(idx_1), ...
                fix_seq(idx_1), ... fix_seq(idx_2), ...
                thr, ...
                max_n, ...
                false, ...
                true);
            
            res.mem1.enc_seq(f,stim,i,:)  = recalls_out;
            res.mem1.enc_repeated(f,stim,i,:)  = enc_reins;
            res.mem1.enc_replay(f,stim,i,:)  = lag_one;
            res.mem1.enc_rev_replay(f,stim,i,:)  = lag_minusone;
            
            if run_perm

                for n = 1:length(null_idx_1)
                    [~, enc_reins, ~, ~, lag_one, lag_minusone, ...
                        recalls_out, ~] = compute_cvp(fix_all(null_idx_1(n)), ...
                        fix_all(null_idx_1(n)), ... %fix_seq(idx_1), ...
                        thr, ...
                        max_n, ...
                        false, ...
                        true);
                    
                    null_res.mem1.enc_seq(f,stim,i).dat(n,:) = recalls_out;
                    null_res.mem1.enc_repeated(f,stim,i).dat(n,:) = enc_reins;
                    null_res.mem1.enc_replay(f,stim,i).dat(n,:)  = lag_one;
                    null_res.mem1.enc_rev_replay(f,stim,i).dat(n,:)  = lag_minusone;
                    
                    [null_res.mem1.cvp(f,stim,i).dat(n,:), enc_reins, ~, ~, lag_one, lag_minusone, ...
                     recalls_out, ~] = compute_cvp(fix_all(null_idx_1(n)), ...
                        fix_seq(idx_2), ...
                        thr, ...
                        max_n, ...
                        false, ...
                        false);
                    
                    null_res.mem1.n_fix(f,stim,i).dat(n,:) = size(fix_all(null_idx_1(n)).fix,2);
                    null_res.mem1.sub_view(f,stim,i).dat(n,:) = enc_reins;
                    null_res.mem1.sub_reins(f,stim,i).dat(n,:) = lag_one;
                    null_res.mem1.sub_rev_reins(f,stim,i).dat(n,:) = lag_minusone;
                    
                end
                
                
                null_res.mem1.category(f,stim,i) = id(stim,1);
                null_res.mem1.fileno(f,stim,i) = id(stim,2);
                
            end

            % recall to recall
            
            % reinstatement/replay events within the retrieval period
            [~, within_rec_reins, ~, ~, within_rec_lag_one, within_rec_lag_minusone, ...
                within_rec_recalls_out, within_rec_lags_out] = compute_cvp(fix_seq(idx_2), ...
                fix_seq(idx_2), ... fix_seq(idx_2), ...
                thr, ...
                max_n, ...
                false, ...
                true);
            
            res.mem1.rec_seq(f,stim,i,:)  = within_rec_recalls_out;
            res.mem1.rec_repeated(f,stim,i,:)  = within_rec_reins;
            res.mem1.rec_replay(f,stim,i,:)  = within_rec_lag_one;
            res.mem1.rec_rev_replay(f,stim,i,:)  = within_rec_lag_minusone;
            
            res.mem1.category(f,stim,i) = id(stim,1);
            res.mem1.fileno(f,stim,i) = id(stim,2);
            
        end
    end
end

res.mem1.iter_lags = iter_lags;
res.mem1.combs = combs;

% dataset 2

fixmat = get_fixmat('etdb_v1.0.hdf5', 'Memory II');
[fix_all, t_fix, subs, filenos, iter] = convert_data(fixmat);

subs = unique(fixmat.SUBJECTINDEX);

for f = 1:length(subs)
    filter_struct.SUBJECTINDEX = subs(f);
    
    [fix_seq, idx] = prep_seq(fixmat, filter_struct);
    
    iter = [fixmat.iteration(idx)];
    category  = [fixmat.category(idx)];
    filenumber = [fixmat.filenumber(idx)];
    condition = [fixmat.condition(idx)];
    
    fcnt = 1;
    for s = 1:length(fix_seq)
        nf = length(fix_seq(s).dur);
        % include repetition
        fix_seq(s).iter = iter(fcnt:fcnt+nf-1);
        fix_seq(s).iter = fix_seq(s).iter';
        
        % include stimulus category
        fix_seq(s).category = category(fcnt:fcnt+nf-1);
        fix_seq(s).category = fix_seq(s).category';
        
        % and fileno
        fix_seq(s).fileno = filenumber(fcnt:fcnt+nf-1);
        fix_seq(s).fileno = fix_seq(s).fileno';
        
        % and fileno
        fix_seq(s).condition = condition(fcnt:fcnt+nf-1);
        fix_seq(s).condition = fix_seq(s).condition';
        
        fcnt = fcnt+nf;
    end
    
    iter = nan(1,length(fix_seq));
    category = nan(1,length(fix_seq));
    fileno = nan(1,length(fix_seq));
    condition = nan(1,length(fix_seq));
    
    for s = 1:length(fix_seq)
        iter(s) = mode(fix_seq(s).iter);
        category(s) = mode(fix_seq(s).category);
        fileno(s) = mode(fix_seq(s).fileno);
        condition(s) = mode(fix_seq(s).condition);
    end
    
    % compute cvp, based off of previous viewing
    
    [id, ~, ib] = unique([category; fileno; condition]','rows');
    
    thr = 104;
    max_n = 10;
    
    combs = nchoosek(1:5, 2);
    iter_lags = diff(combs, [], 2);
    
    for i = 1:size(combs,1)
        for stim = 1:length(id)
            
            idx_1 = category == id(stim,1) & ...
                fileno  == id(stim,2) & ...
                iter == combs(i,1);
            
            idx_2 = category == id(stim,1) & ...
                fileno  == id(stim,2) & ...
                iter == combs(i,2);
            
            if ~any(idx_1) || ~any(idx_2)
                continue
            end
            
             % adding permutation (subject) nulls
            
            null_idx_1 = find([fix_all.SUBJECTINDEX] ~= subs(f) & ...
                [fix_all.category] == id(stim,1) & ...
                [fix_all.filenumber] == id(stim,2) & ...
                [fix_all.iteration] == combs(i,1));
            
            % use idx_1 and idx_2 to examine whether fixations are
            % reinstated or replayed during the second viewing
            
            [~, enc_reins, ~, ~, lag_one, lag_minusone] = compute_cvp(fix_seq(idx_1), ...
                fix_seq(idx_2), ...
                thr, ...
                max_n, ...
                false, ...
                false);

            res.mem2.n_fix(f,stim,i) = size(fix_seq(idx_1).fix,2);
            res.mem2.sub_view(f,stim,i,:) = enc_reins;
            res.mem2.sub_reins(f,stim,i,:) = lag_one;
            res.mem2.sub_rev_reins(f,stim,i,:) = lag_minusone;
   
            % reinstatement/replay events within the encoding period
            [~, enc_reins, ~, ~, lag_one, lag_minusone, ...
                recalls_out, ~] = compute_cvp(fix_seq(idx_1), ...
                fix_seq(idx_1), ... fix_seq(idx_2), ...
                thr, ...
                max_n, ...
                false, ...
                true);
            
            res.mem2.category(f,stim,i) = id(stim,1);
            res.mem2.fileno(f,stim,i) = id(stim,2);
            res.mem2.condition(f,stim,i) = id(stim,3);
            
            null_res.mem2.category(f,stim,i) = id(stim,1);
            null_res.mem2.fileno(f,stim,i) = id(stim,2);
            null_res.mem2.condition(f,stim,i) = id(stim,3);
            
            res.mem2.enc_seq(f,stim,i,:)  = recalls_out;
            res.mem2.enc_repeated(f,stim,i,:)  = enc_reins;
            res.mem2.enc_replay(f,stim,i,:)  = lag_one;
            res.mem2.enc_rev_replay(f,stim,i,:)  = lag_minusone;
            
            if run_perm

                for n = 1:length(null_idx_1)
                    [~, enc_reins, ~, ~, lag_one, lag_minusone, ...
                     recalls_out, ~] = compute_cvp(fix_all(null_idx_1(n)), ...
                                                fix_all(null_idx_1(n)), ... %fix_seq(idx_1), ...
                                                thr, ...
                                                max_n, ...
                                                false, ...
                                                true);
                            
                    null_res.mem2.enc_seq(f,stim,i).dat(n,:) = recalls_out;
                    null_res.mem2.enc_repeated(f,stim,i).dat(n,:) = enc_reins;
                    null_res.mem2.enc_replay(f,stim,i).dat(n,:)  = lag_one;
                    null_res.mem2.enc_rev_replay(f,stim,i).dat(n,:)  = lag_minusone;
                    
                    % also reins
                     [null_res.mem2.cvp(f,stim,i).dat(n,:), enc_reins, ~, ~, lag_one, lag_minusone, ...
                     recalls_out, ~] = compute_cvp(fix_all(null_idx_1(n)), ...
                        fix_seq(idx_2), ...
                        thr, ...
                        max_n, ...
                        false, ...
                        false);
                    
                    null_res.mem2.n_fix(f,stim,i).dat(n,:) = size(fix_all(null_idx_1(n)).fix,2);
                    null_res.mem2.sub_view(f,stim,i).dat(n,:) = enc_reins;
                    null_res.mem2.sub_reins(f,stim,i).dat(n,:) = lag_one;
                    null_res.mem2.sub_rev_reins(f,stim,i).dat(n,:) = lag_minusone;
                    
                end

            end
                      
            % reinstatement/replay events within the retrieval period
            [~, within_rec_reins, ~, ~, within_rec_lag_one, within_rec_lag_minusone, ...
                within_rec_recalls_out, within_rec_lags_out] = compute_cvp(fix_seq(idx_2), ...
                fix_seq(idx_2), ... fix_seq(idx_2), ...
                thr, ...
                max_n, ...
                false, ...
                true);
            
            res.mem2.rec_seq(f,stim,i,:)  = within_rec_recalls_out;
            res.mem2.rec_repeated(f,stim,i,:)  = within_rec_reins;
            res.mem2.rec_replay(f,stim,i,:)  = within_rec_lag_one;
            res.mem2.rec_rev_replay(f,stim,i,:)  = within_rec_lag_minusone;
            
        end
    end
end

res.mem2.iter_lags = iter_lags;
res.mem2.combs = combs;

end

function [fix_seq, idx] = prep_seq(fixmat, filter_struct)
% PREP SEQ - prepares fixation and duration information

% extract the trials from the overall fixation structure
idx = false(length(filter_struct), size(fixmat.x,1));
fns = fieldnames(filter_struct);
for f = 1:length(fns)
    idx(f,:) = ismember(fixmat.(fns{f}), filter_struct.(fns{f}));
end
if size(idx,1) ~= 1
    idx = all(idx);
end

trials = unique(fixmat.trial(idx));

for s = 1:length(trials)
    t_idx = fixmat.trial == trials(s) & idx';
    
    fix_seq(s).fix = [fixmat.x(t_idx), fixmat.y(t_idx)]';
    fix_seq(s).dur = [fixmat.end(t_idx) - fixmat.start(t_idx)]';
    
    
    for f = 1:length(fns)
        fix_seq(s).(fns{f}) = filter_struct.(fns{f});
    end
    
end


end

function [cvp, enc_reins, actual, possible, lag_one, lag_minusone, recalls_out, lags_out] = compute_cvp(s1, s2, thr, ...
    max_n, to_repeat, enc2enc_flag)

if isempty(s1) || isempty(s1.fix) || isempty(s2)
    cvp = nan(size(-max_n+1:1:max_n-1));
    cvp(end+1) = nan;
    actual = cvp;
    possible = cvp;
    enc_reins = nan(1,200);
    lag_one = nan(1,200);
    lag_minusone = nan(1,200);
    return
end

distmat = dist([s1.fix s2.fix]);
recmat = distmat < thr;
tmp = recmat(size(s1.fix,2)+1:end, 1:size(s1.fix,2));

% compute 'serial position' of recalls
tmp_dist = distmat(size(s1.fix,2)+1:end, 1:size(s1.fix,2));
recalls = nan(1,size(tmp_dist,1));
novel_mask = nan(1, size(tmp_dist,1));

if enc2enc_flag
    scnt = 1;
    for e = 1:size(tmp_dist,1)
        this_dist = tmp_dist(e,:);
        if e <= size(tmp_dist,2) % added for null condition where lengths may not match
        this_dist(e) = nan;
        end
        idx = find(this_dist==min(this_dist) & tmp(e,:), 1, 'first');
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
else
    for e = 1:size(tmp_dist,1)
        this_dist = tmp_dist(e,:);
        idx = find(this_dist==min(this_dist) & tmp(e,:), 1, 'first');
        if ~isempty(idx)
            recalls(e) = idx;
        else
            recalls(e) = 0; % new
        end
    end
    scnt = size(s1.fix,2);
end

recalls(recalls==0)=nan;

lags = [diff(recalls,[],2) nan];

% output some encoding information

if enc2enc_flag
    reinstated = unique(recalls(~novel_mask));
    enc_reins = nan(1,200);
    if ~isempty(reinstated)
        for e = 1:length(reinstated)
            enc_reins(recalls == reinstated(e) & novel_mask) = 1; % later viewed in this sequence
            enc_reins(recalls == reinstated(e) & ~novel_mask) = 0;
        end
    end
else
    reinstated = unique(recalls(~isnan(lags)));
    enc_reins = nan(1,200);
    if ~isempty(reinstated)
        enc_reins(reinstated) = 1;
        enc_reins(enc_reins~=1) = 0;
    end
end

if enc2enc_flag
    to_mask = [nan novel_mask(1:end-1)];
    
    % lag one
    reinstated = unique(recalls(~novel_mask & lags ==1 & to_mask == 0));
    lag_one = nan(1,200);
    if ~isempty(reinstated)
        for e = 1:length(reinstated)
            lag_one(recalls == reinstated(e) & novel_mask) = 1;
            lag_one(recalls == reinstated(e) & ~novel_mask) = 0;
        end
    end
    
    % lag minus one
    reinstated = unique(recalls(~novel_mask & lags ==-1 & to_mask == 0));
    lag_minusone = nan(1,200);
    if ~isempty(reinstated)
        for e = 1:length(reinstated)
            lag_minusone(recalls == reinstated(e) & novel_mask) = 1;
            lag_minusone(recalls == reinstated(e) & ~novel_mask) = 0;
        end
    end    
else
    % lag_one
    lag_one = nan(1,200);
    reinstated = unique(recalls(lags==1));
    if ~isempty(reinstated)
        lag_one(reinstated) = 1;
        lag_one(lag_one~=1) = 0;
    end
    
    % lag_minusone
    lag_minusone = nan(1,200);
    reinstated = unique(recalls(lags==-1));
    if ~isempty(reinstated)
        lag_minusone(reinstated) = 1;
        lag_minusone(lag_minusone~=1) = 0;
    end
end

% output some retrieval information
recalls_out = zeros(1,200);
recalls_out(1:length(recalls)) = recalls;

lags_out = nan(1,200);
lags_out(1:length(lags)) = lags;

to_mask_pres = true(size(1:scnt));
from_mask_pres =  true(size(1:scnt));

to_mask_rec = ~isnan(recalls);
from_mask_rec = ~isnan(recalls);

if enc2enc_flag
    if to_repeat
        from_mask_rec = novel_mask == 1; %& lags~= 0
        to_mask_rec = novel_mask == 0;
    else % from repeated items to repeated
        from_mask_rec = novel_mask == 0; %& lags~= 0;
        to_mask_rec = novel_mask == 0; %true(size(novel_mask));
    end
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

function value_counts = collect(data_matrix, values)

% sanity!
if ~exist('data_matrix', 'var')
  error('You must pass a data_matrix')
elseif ~exist('values', 'var')
  error('You must pass a vector of values to count')
elseif ~isvector(values)
  error('values must be a vector')
end

value_counts = zeros(1, length(values));
for i = 1:length(values)
  value_counts(i) = count(data_matrix, values(i));
end

end

function value_count = count(data_matrix, value)


% sanity checks:
if length(value) > 1
  error('value must be a scalar')
end

% very simple (perhaps almost naive) implementation:
value_count = nnz(data_matrix == value);

end

function d = lag(sp1, sp2, params)
  d = sp2 - sp1;
end


function [fix_all, t_fix, subs, filenos, iter] = convert_data(fixmat)
% creates a fixation structure to work with CVP code

subs = unique(fixmat.SUBJECTINDEX);
filenos = unique(fixmat.filenumber);
iter = unique(fixmat.iteration);

% add condition for mem1
if ~isfield(fixmat,'condition')
    fixmat.condition = ones(size(fixmat.fix));
end

[uni_fix, j] = unique([fixmat.SUBJECTINDEX fixmat.filenumber fixmat.iteration fixmat.category fixmat.condition],'rows');

all_subs = uni_fix(:,1);
all_file = uni_fix(:,2);
all_iter = uni_fix(:,3);
all_cat = uni_fix(:,4);
all_cond = uni_fix(:,5);


for idx = 1:length(j)
    filter_struct.SUBJECTINDEX = all_subs(idx);
    filter_struct.filenumber = all_file(idx);
    filter_struct.iteration = all_iter(idx);
    filter_struct.category = all_cat(idx);
    filter_struct.condition = all_cond(idx);
    
    fix_all(idx) = prep_seq(fixmat, filter_struct);
end

filter_struct = [];
for f = 1:length(subs)
    filter_struct.SUBJECTINDEX = subs(f);
    [fix_seq, idx] = prep_seq(fixmat, filter_struct);
    
    iter = fixmat.iteration(idx);
    category  = fixmat.category(idx);
    filenumber = fixmat.filenumber(idx);
    condition = fixmat.condition(idx);
    
    fcnt = 1;
    for s = 1:length(fix_seq)
        nf = length(fix_seq(s).dur);
        % include repetition
        t_fix(f).fix_seq(s).iter = iter(fcnt:fcnt+nf-1);
        t_fix(f).fix_seq(s).iter = t_fix(f).fix_seq(s).iter';
        
        % include stimulus category
        t_fix(f).fix_seq(s).category = category(fcnt:fcnt+nf-1);
        t_fix(f).fix_seq(s).category = t_fix(f).fix_seq(s).category';
        
        % and fileno
        t_fix(f).fix_seq(s).fileno = filenumber(fcnt:fcnt+nf-1);
        t_fix(f).fix_seq(s).fileno = t_fix(f).fix_seq(s).fileno';
        
        % and condition, if appropriate
        t_fix(f).fix_seq(s).condition = condition(fcnt:fcnt+nf-1);
        t_fix(f).fix_seq(s).condition = t_fix(f).fix_seq(s).condition';
        
        fcnt = fcnt+nf;
    end
    
    iter = nan(1, length(fix_seq));
    category = nan(1, length(fix_seq));
    fileno = nan(1, length(fix_seq));
    condition = nan(1, length(fix_seq));
    for s = 1:length(fix_seq)
        iter(s) = mode(t_fix(f).fix_seq(s).iter);
        category(s) = mode(t_fix(f).fix_seq(s).category);
        fileno(s) = mode(t_fix(f).fix_seq(s).fileno);
        condition(s) = mode(t_fix(f).fix_seq(s).condition);
    end
    
    t_fix(f).fix_seq = fix_seq;
    t_fix(f).iter = iter;
    t_fix(f).category = category;
    t_fix(f).fileno = fileno;
    t_fix(f).condition = condition;
    
end

end