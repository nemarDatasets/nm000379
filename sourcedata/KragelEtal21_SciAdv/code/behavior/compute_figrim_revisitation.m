function [res, null_res] = compute_figrim_revisitation()
% COMPUTE_FIGRIM_REVISITATION computes fixation-sequence replay for revisitation
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

% if you want to also compute stimulus-driven measures (slow)
run_perm = false; 

obj = load('allImages_release.mat');
allImages = obj.allImages; clear obj

% some parameters for computing cvp curves
thr = 66;
max_n = 5;

% initialize arrays
cvp = nan(length(allImages(1).userdata), ...
    length(allImages), ...
    2*max_n);
cvp_hit = cvp;
cvp_miss = cvp;

actual = cvp; actual2 = actual; actual3 = actual;
possible = cvp; possible2 = possible; possible3 = possible;

actual_hit = cvp; actual2_hit = actual; actual3_hit = actual;
possible_hit = cvp; possible2_hit = possible; possible3_hit = possible;

actual_miss = cvp; actual2_miss = actual; actual3_miss = actual;
possible_miss = cvp; possible2_miss = possible; possible3_miss = possible;

cvp2 = cvp;
cvp2_hit = cvp;
cvp2_miss = cvp;

cvp3 = cvp;
cvp3_hit = cvp;
cvp3_miss = cvp;

n_fix = nan(length(allImages(1).userdata), length(allImages));
n_fix2 = n_fix; n_fix3 = n_fix;

p_rep2 = nan(length(allImages(1).userdata), length(allImages), 2);
p_rep3 = p_rep2;


enc_reins = nan(length(allImages(1).userdata), length(allImages), 200);
enc_reins_hit = enc_reins; enc_reins_miss = enc_reins;
enc_reins2 = enc_reins; enc_reins3 = enc_reins;

fr = nan(length(allImages(1).userdata), length(allImages), 5);

for s = 1:length(allImages(1).userdata) % subject
    
    for i = 1:length(allImages) % image
        
        if ~isempty(allImages(i).userdata(s).fixations)
            
            %% gather observed sequences for this image/subject
            enc_seq.fix = allImages(i).userdata(s).fixations.enc';
            rec_seq.fix = allImages(i).userdata(s).fixations.rec';
            rec2_seq.fix = allImages(i).userdata(s).fixations.rec2';
            
            n_fix(s,i) = size(enc_seq.fix,2);
            n_fix2(s,i) = size(rec_seq.fix,2);
            n_fix3(s,i) = size(rec2_seq.fix,2);
            
            %% gather all other subject's fixations to this image
            n_subs = length(allImages(i).userdata);
            diff_subs = setdiff(1:n_subs, s);
            
            % remove subjects without fixations for this image
            bs = false(size(diff_subs));
            for d = 1:length(diff_subs)
                if isempty(allImages(i).userdata(diff_subs(d)).fixations)
                    bs(d) = true;
                end
            end
            diff_subs = diff_subs(~bs);
            
            % get fix for each of these
            clear null_enc_seq null_rec_seq null_rec2_seq            
            for d = 1:length(diff_subs)
                null_enc_seq(d).fix = allImages(i).userdata(diff_subs(d)).fixations.enc';
                null_rec_seq(d).fix = allImages(i).userdata(diff_subs(d)).fixations.rec';
                null_rec2_seq(d).fix = allImages(i).userdata(diff_subs(d)).fixations.rec2';
            end
            
            
            %% final recognition stats
            
            fr(s,i,:) = final_recog(enc_seq, rec_seq, rec2_seq, thr);
            
            %% info about subs. reinstatement and forward/rev. replay - actual data
            [~, er, ~, ~, lo, lmo, ro] = compute_cvp(enc_seq, ...
                rec_seq, ...
                thr, ...
                max_n, ...
                false, ...
                false);
            
            [~, er2, ~, ~, lo2, lmo2, ro2] = compute_cvp(enc_seq, ...
                rec2_seq, ...
                thr, ...
                max_n, ...
                false, ...
                false);
            
            [~, er3, ~, ~, lo3, lmo3, ro3] = compute_cvp(rec_seq, ...
                rec2_seq, ...
                thr, ...
                max_n, ...
                false, ...
                false);
            
            % enc2enc
            n_fix(s,i) = size(enc_seq.fix,2);
            sub_view(s,i,:) = er;
            sub_reins(s,i,:) = lo;
            sub_rev_reins(s,i,:) = lmo;
            
            % rec2rec
            n_fix2(s,i) = size(rec_seq.fix,2);
            sub_view2(s,i,:) = er2;
            sub_reins2(s,i,:) = lo2;
            sub_rev_reins2(s,i,:) = lmo2;
            
            % rec2rec2
            n_fix3(s,i) = size(rec2_seq.fix,2);
            sub_view3(s,i,:) = er3;
            sub_reins3(s,i,:) = lo3;
            sub_rev_reins3(s,i,:) = lmo3;
            
            %% reinstatement/replay events within the encoding period - actual data
            
            [cvp(s,i,:), enc_reins(s,i,:), actual(s,i,:), possible(s,i,:), lag_one, lag_minusone, ...
                recalls_out] = compute_cvp(enc_seq, enc_seq, thr, max_n, false, true);
            [cvp2(s,i,:), enc_reins2(s,i,:), actual2(s,i,:), possible2(s,i,:), lag_one2, lag_minusone2, ...
                recalls_out2] = compute_cvp(rec_seq, rec_seq, thr, max_n, false, true);
            [cvp3(s,i,:), enc_reins3(s,i,:), actual3(s,i,:), possible3(s,i,:), lag_one3, lag_minusone3, ...
                recalls_out3] = compute_cvp(rec2_seq, rec2_seq, thr, max_n, false, true);
            
            enc_ord(s,i,:)  = recalls_out;
            enc_replay(s,i,:)  = lag_one;
            enc_rev_replay(s,i,:)  = lag_minusone;
            
            enc_ord2(s,i,:)  = recalls_out2;
            enc_replay2(s,i,:)  = lag_one2;
            enc_rev_replay2(s,i,:)  = lag_minusone2;
            
            enc_ord3(s,i,:)  = recalls_out3;
            enc_replay3(s,i,:)  = lag_one3;
            enc_rev_replay3(s,i,:)  = lag_minusone3;
            
            if allImages(i).userdata(s).SDT(2) == 1
                
                % enc2enc cvp information
                cvp_hit(s,i,:) = cvp(s,i,:);
                actual_hit(s,i,:) = actual(s,i,:);
                possible_hit(s,i,:) = possible(s,i,:);
                
                n_fix_hit(s,i) = size(enc_seq.fix,2);
                sub_view_hit(s,i,:) = sub_view(s,i,:);
                sub_reins_hit(s,i,:) = sub_reins(s,i,:);
                sub_rev_reins_hit(s,i,:) = sub_rev_reins(s,i,:);
                
                enc_ord_hit(s,i,:)  = enc_ord(s,i,:);
                enc_repeated_hit(s,i,:)  = enc_reins(s,i,:);
                enc_replay_hit(s,i,:)  = enc_replay(s,i,:);
                enc_rev_replay_hit(s,i,:)  = enc_rev_replay(s,i,:);
                
            elseif allImages(i).userdata(s).SDT(2) == 3
                cvp_miss(s,i,:) = cvp(s,i,:);
                actual_miss(s,i,:) = actual(s,i,:);
                possible_miss(s,i,:) = possible(s,i,:);
                
                n_fix_miss(s,i) = size(enc_seq.fix,2);
                sub_view_miss(s,i,:) = sub_view(s,i,:);
                sub_reins_miss(s,i,:) = sub_reins(s,i,:);
                sub_rev_reins_miss(s,i,:) = sub_rev_reins(s,i,:);
                
                enc_ord_miss(s,i,:)  = enc_ord(s,i,:);
                enc_repeated_miss(s,i,:)  = enc_reins(s,i,:);
                enc_replay_miss(s,i,:)  = enc_replay(s,i,:);
                enc_rev_replay_miss(s,i,:)  = enc_rev_replay(s,i,:);
            end
            
            if allImages(i).userdata(s).SDT(3) == 1
                cvp2_hit(s,i,:) = cvp2(s,i,:);
                actual2_hit(s,i,:) = actual2(s,i,:);
                possible2_hit(s,i,:) = possible2(s,i,:);
                
                n_fix2_hit(s,i) = size(rec_seq.fix,2);
                sub_view2_hit(s,i,:) = sub_view2(s,i,:);
                sub_reins2_hit(s,i,:) = sub_reins2(s,i,:);
                sub_rev_reins2_hit(s,i,:) = sub_rev_reins2(s,i,:);
                
                enc_ord2_hit(s,i,:)  = enc_ord2(s,i,:);
                enc_repeated2_hit(s,i,:)  = enc_reins2(s,i,:);
                enc_replay2_hit(s,i,:)  = enc_replay2(s,i,:);
                enc_rev_replay2_hit(s,i,:)  = enc_rev_replay2(s,i,:);
                
            elseif allImages(i).userdata(s).SDT(3) == 3
                cvp2_miss(s,i,:) = cvp2(s,i,:);
                actual2_miss(s,i,:) = actual2(s,i,:);
                possible2_miss(s,i,:) = possible2(s,i,:);
                
                n_fix2_miss(s,i) = size(rec_seq.fix,2);
                sub_view2_miss(s,i,:) = sub_view2(s,i,:);
                sub_reins2_miss(s,i,:) = sub_reins2(s,i,:);
                sub_rev_reins2_miss(s,i,:) = sub_rev_reins2(s,i,:);
                
                enc_ord2_miss(s,i,:)  = enc_ord2(s,i,:);
                enc_repeated2_miss(s,i,:)  = enc_reins2(s,i,:);
                enc_replay2_miss(s,i,:)  = enc_replay2(s,i,:);
                enc_rev_replay2_miss(s,i,:)  = enc_rev_replay2(s,i,:);
            end
            
            if allImages(i).userdata(s).SDT(3) == 1
                cvp3_hit(s,i,:) = cvp3(s,i,:);
                actual3_hit(s,i,:) = actual3(s,i,:);
                possible3_hit(s,i,:) = possible3(s,i,:);
                
                n_fix3_hit(s,i) = size(rec2_seq.fix,2);
                sub_view3_hit(s,i,:) = sub_view3(s,i,:);
                sub_reins3_hit(s,i,:) = sub_reins3(s,i,:);
                sub_rev_reins3_hit(s,i,:) = sub_rev_reins3(s,i,:);
                
                enc_ord3_hit(s,i,:)  = enc_ord3(s,i,:);
                enc_repeated3_hit(s,i,:)  = enc_reins3(s,i,:);
                enc_replay3_hit(s,i,:)  = enc_replay3(s,i,:);
                enc_rev_replay3_hit(s,i,:)  = enc_rev_replay3(s,i,:);
                
            elseif allImages(i).userdata(s).SDT(3) == 3
                cvp3_miss(s,i,:) = cvp3(s,i,:);
                actual3_miss(s,i,:) = actual3(s,i,:);
                possible3_miss(s,i,:) = possible3(s,i,:);
                
                n_fix3_miss(s,i) = size(rec2_seq.fix,2);
                sub_view3_miss(s,i,:) = sub_view3(s,i,:);
                sub_reins3_miss(s,i,:) = sub_reins3(s,i,:);
                sub_rev_reins3_miss(s,i,:) = sub_rev_reins3(s,i,:);
                
                enc_ord3_miss(s,i,:)  = enc_ord3(s,i,:);
                enc_repeated3_miss(s,i,:)  = enc_reins3(s,i,:);
                enc_replay3_miss(s,i,:)  = enc_replay3(s,i,:);
                enc_rev_replay3_miss(s,i,:)  = enc_rev_replay3(s,i,:);
            end
            
            %% and stimulus-driven results, from using other subjects data 
            %% at encoding, but the 'current' subject at 'retrieval'
            
            if run_perm

                for d = 1:length(diff_subs)
                    
                    %% info about subs. reinstatement and forward/rev. replay - actual data
                    
                    [~, ner, ~, ~, nlo, nlmo] = compute_cvp(null_enc_seq(d), ...
                        rec_seq, ...
                        thr, ...
                        max_n, ...
                        false, ...
                        false);
                    
                    [~, ner2, ~, ~, nlo2, nlmo2] = compute_cvp(null_enc_seq(d), ...
                        rec2_seq, ...
                        thr, ...
                        max_n, ...
                        false, ...
                        false);
                    
                    [~, ner3, ~, ~, nlo3, nlmo3] = compute_cvp(null_rec_seq(d), ...
                        rec2_seq, ...
                        thr, ...
                        max_n, ...
                        false, ...
                        false);
                    
                    % enc2rec
                    null_res.n_fix(s,i).dat(d) = size(null_enc_seq(d).fix,2);
                    null_res.sub_view(s,i).dat(d,:) = ner;
                    null_res.sub_reins(s,i).dat(d,:) = nlo;
                    null_res.sub_rev_reins(s,i).dat(d,:) = nlmo;
                     
                    % enc2rec
                    null_res.n_fix2(s,i).dat(d) = size(null_enc_seq(d).fix,2);
                    null_res.sub_view2(s,i).dat(d,:) = ner2;
                    null_res.sub_reins2(s,i).dat(d,:) = nlo2;
                    null_res.sub_rev_reins2(s,i).dat(d,:) = nlmo2;
                     
                    % rec2rec2
                    null_res.n_fix3(s,i).dat(d) = size(null_rec_seq(d).fix,2);
                    null_res.sub_view3(s,i,:).dat(d,:) = ner3;
                    null_res.sub_reins3(s,i,:).dat(d,:) = nlo3;
                    null_res.sub_rev_reins3(s,i,:).dat(d,:) = nlmo3;
                    
                    %% and then for the null/stimulus-driven analysis, compute revisitation
                    
                    % reinstatement/replay events within the encoding periods - null data
                               
                    % first to second viewing
                    [~, nenc_reins, ~, ~, nlag_one, nlag_minusone, ...
                        nrecalls_out] = compute_cvp(null_enc_seq(d), ...
                                                   null_enc_seq(d), ...
                                                   thr, ...
                                                   max_n, ...
                                                   false, ...
                                                   true);
                    
                    null_res.n_fix(s,i).dat(d) = size(null_enc_seq(d).fix,2);
                    null_res.hit(s,i).dat(d) = allImages(i).userdata(s).SDT(2) == 1;
                    null_res.hit2(s,i).dat(d) = allImages(i).userdata(s).SDT(3) == 1;
                    
                    null_res.enc_seq(s,i).dat(d,:) = nrecalls_out;
                    null_res.enc_repeated(s,i).dat(d,:) = nenc_reins;
                    null_res.enc_replay(s,i).dat(d,:)  = nlag_one;
                    null_res.enc_rev_replay(s,i).dat(d,:)  = nlag_minusone;
                    
                    % second viewing
                    [~, nenc_reins, ~, ~, nlag_one, nlag_minusone, ...
                        nrecalls_out] = compute_cvp(null_rec_seq(d), ...
                        null_rec_seq(d), ...
                        thr, ...
                        max_n, ...
                        false, ...
                        true);
                    
                    null_res.n_fix_2(s,i).dat(d) = size(null_rec_seq(d).fix,2);
                    null_res.hit_2(s,i).dat(d) = allImages(i).userdata(s).SDT(2) == 1;
                    null_res.hit2_2(s,i).dat(d) = allImages(i).userdata(s).SDT(3) == 1;
                    
                    null_res.enc_seq_2(s,i).dat(d,:) = nrecalls_out;
                    null_res.enc_repeated_2(s,i).dat(d,:) = nenc_reins;
                    null_res.enc_replay_2(s,i).dat(d,:)  = nlag_one;
                    null_res.enc_rev_replay_2(s,i).dat(d,:)  = nlag_minusone;
                    
                end

            end
            
        end
    end
end


res.n_fix = n_fix;
res.n_fix2 = n_fix2;
res.n_fix3 = n_fix3;

% prob repeating location
res.p_rep2 = p_rep2;
res.p_rep3 = p_rep3;

% overall stats

res.cvp_h = squeeze(nanmean(cvp_hit,2));
res.cvp_m = squeeze(nanmean(cvp_miss,2));

res.cvp2_h = squeeze(nanmean(cvp2_hit,2));
res.cvp2_m = squeeze(nanmean(cvp2_miss,2));

res.cvp3_h = squeeze(nanmean(cvp3_hit,2));
res.cvp3_m = squeeze(nanmean(cvp3_miss,2));


% --hits--

% 1
res.n_fix_hit = n_fix_hit;
res.sub_view_hit = sub_view_hit;
res.sub_reins_hit = sub_reins_hit;
res.sub_rev_reins_hit = sub_rev_reins_hit;

res.enc_ord_hit = enc_ord_hit;
res.enc_repeated_hit = enc_repeated_hit;
res.enc_replay_hit = enc_replay_hit;
res.enc_rev_replay_hit = enc_rev_replay_hit;

% 2

res.n_fix2_hit = n_fix2_hit;
res.sub_view2_hit = sub_view2_hit;
res.sub_reins2_hit = sub_reins2_hit;
res.sub_rev_reins2_hit = sub_rev_reins2_hit;

res.enc_ord2_hit = enc_ord2_hit;
res.enc_repeated2_hit = enc_repeated2_hit;
res.enc_replay2_hit = enc_replay2_hit;
res.enc_rev_replay2_hit = enc_rev_replay2_hit;

% 3
res.n_fix3_hit = n_fix3_hit;
res.sub_view3_hit = sub_view3_hit;
res.sub_reins3_hit = sub_reins3_hit;
res.sub_rev_reins3_hit = sub_rev_reins3_hit;

res.enc_ord3_hit = enc_ord3_hit;
res.enc_repeated3_hit = enc_repeated3_hit;
res.enc_replay3_hit = enc_replay3_hit;
res.enc_rev_replay3_hit = enc_rev_replay3_hit;

% --misses--

% 1
res.n_fix_miss = n_fix_miss;
res.sub_view_miss = sub_view_miss;
res.sub_reins_miss = sub_reins_miss;
res.sub_rev_reins_miss = sub_rev_reins_miss;

res.enc_ord_miss = enc_ord_miss;
res.enc_repeated_miss = enc_repeated_miss;
res.enc_replay_miss = enc_replay_miss;
res.enc_rev_replay_miss = enc_rev_replay_miss;

% 2

res.n_fix2_miss = n_fix2_miss;
res.sub_view2_miss = sub_view2_miss;
res.sub_reins2_miss = sub_reins2_miss;
res.sub_rev_reins2_miss = sub_rev_reins2_miss;

res.enc_ord2_miss = enc_ord2_miss;
res.enc_repeated2_miss = enc_repeated2_miss;
res.enc_replay2_miss = enc_replay2_miss;
res.enc_rev_replay2_miss = enc_rev_replay2_miss;

% 3
res.n_fix3_miss = n_fix3_miss;
res.sub_view3_miss = sub_view3_miss;
res.sub_reins3_miss = sub_reins3_miss;
res.sub_rev_reins3_miss = sub_rev_reins3_miss;

res.enc_ord3_miss = enc_ord3_miss;
res.enc_repeated3_miss = enc_repeated3_miss;
res.enc_replay3_miss = enc_replay3_miss;
res.enc_rev_replay3_miss = enc_rev_replay3_miss;

% orig code

res.cvp_hit_all = squeeze(nanmean(cat(3, res.cvp_h, ...
    res.cvp2_h, ...
    res.cvp3_h),3));

res.cvp_miss_all = squeeze(nanmean(cat(3, res.cvp_m, ...
    res.cvp2_m, ...
    res.cvp3_m),3));

res.cvp = cvp;
res.cvp2 = cvp2;
res.cvp3 = cvp3;

res.enc_reins = enc_reins;
res.enc_reins_hit = enc_reins_hit;
res.enc_reins_miss = enc_reins_miss;

res.actual = actual;
res.actual_hit = actual_hit;
res.actual_miss = actual_miss;
res.possible = possible;
res.possible_hit = possible_hit;
res.possible_miss = possible_miss;

res.actual2 = actual2;
res.actual2_hit = actual2_hit;
res.actual2_miss = actual2_miss;
res.possible2 = possible2;
res.possible2_hit = possible2_hit;
res.possible2_miss = possible2_miss;

res.actual3 = actual3;
res.actual3_hit = actual3_hit;
res.actual3_miss = actual3_miss;
res.possible3 = possible3;
res.possible3_hit = possible3_hit;
res.possible3_miss = possible3_miss;

end

function out = final_recog(s1, s2, s3, thr)
% computes proportion of fixations at r2:
% out contains

% 1 - proportions from e1
% 2 - proportions from r1
% 3 - proportions from only e1
% 4 - proportions from only r1
% 5 - proportions from e1 and r1

if isempty(s1) || isempty(s1.fix) || ...
        isempty(s2) || isempty(s2.fix) || ...
        isempty(s3) || isempty(s3.fix)
    out = nan(1,5);
    return
end

distmat = dist([s1.fix s2.fix s3.fix]);
recmat = distmat < thr;

s1_idx = 1:size(s1.fix,2);
s2_idx = (size(s1.fix,2)+1):(size(s1.fix,2)+size(s2.fix,2));

tmp = recmat(size(s1.fix,2)+size(s2.fix,2)+1:end, ...
    1:size(s1.fix,2)+size(s2.fix,2));


% compute 'serial position' of recalls
tmp_dist = distmat(size(s1.fix,2)+size(s2.fix,2)+1:end, ...
    1:size(s1.fix,2)+size(s2.fix,2));

recalls = nan(2,size(tmp_dist,1));
for e = 1:size(tmp_dist,1)
    this_dist = tmp_dist(e,:);
    idx = find(this_dist(s1_idx)==min(this_dist(s1_idx)) & tmp(e,s1_idx));
    if ~isempty(idx)
        recalls(1,e) = idx;
    else
        recalls(1,e) = 0; % new
    end
    
    idx = find(this_dist(s2_idx)==min(this_dist(s2_idx)) & tmp(e,s2_idx));
    if ~isempty(idx)
        recalls(2,e) = idx;
    else
        recalls(2,e) = 0; % new
    end
end


% compute 'repetition number' of recalls

iter_rec = false(size(tmp,1),2);

for e = 1:size(tmp,1)
    this_rec = tmp(e,:);
    idx = find(this_rec);
    if any(ismember(idx, s1_idx))
        iter_rec(e,1) = true;
    end
    
    if any(ismember(idx, s2_idx))
        iter_rec(e,2) = true;
    end
end

out = [nanmean(iter_rec,1) ...
    nanmean(iter_rec(:,1)==1 & iter_rec(:,2)==0) ...
    nanmean(iter_rec(:,1)==0 & iter_rec(:,2)==1) ...
    nanmean(all(iter_rec,2))];

end

function p_rep = compute_prob_repeat(s1, s2, thr)
% added 8/6/2019 to compare proportion of repeats for novel vs. repeated
% items

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
    max_n, to_repeat, enc2enc_flag)

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

if enc2enc_flag
    scnt = 1;
    for e = 1:size(tmp_dist,1)
        this_dist = tmp_dist(e,:);
        this_dist(e) = nan;

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
    from_mask_rec = novel_mask == 1;
    to_mask_rec = novel_mask == 0;
else % from repeated items to repeated
    from_mask_rec = novel_mask == 0;
    to_mask_rec = novel_mask == 0;
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


function [transits_array] = possible_transitions(serial_position, ...
    prior_recalls, ...
    ~, params)


to_mask_pres = params.to_mask_pres;
from_mask_pres = params.from_mask_pres;

% transitions from an intrusion, a repeat, or a masked out position
% are invalid, so there are no possible transitions from them
if serial_position < 1 || ~from_mask_pres(serial_position)
    transits_array = [];
    return
end

% calculate all possible lags
transits_array = find(to_mask_pres) - serial_position;

end

function [trans_row, cond_cell] = conditional_transitions(data_row, ...
    from_mask, to_mask, transit_func, condition, step, params)


% sanity checks
if ~islogical(from_mask) || ~islogical(to_mask)
    error('Masks must be logical variables');
end
if step < 1
    error('Non-positive steps are not supported');
end
if ~isa(condition, 'function_handle')
    error('condition must be a function handle');
end

% initialization: trans_row should have as many elements as
% data_row - step, as should cond_cell, since that is the maximum
% number of transitions and conditions we can calculate
row_length = size(data_row, 2);
trans_row = NaN(1, row_length - step);
cond_cell = cell(1, row_length - step);

for i = 1:row_length - step
    if from_mask(i) && to_mask(i + step)
        % this transition is not masked out
        from_pt = data_row(i);
        to_pt = data_row(i + step);
        
        % calculate the current transition and condition and append
        % them to trans_row and cond_cell
        transition = transit_func(from_pt, to_pt, params);
        poss_trans = condition(from_pt, data_row(1:i-1), transition, params);
        
        if ~isempty(transition) && any(transition == poss_trans)
            % the transition meets the condition; include it in trans_row
            trans_row(i) = transition;
        end
        
        % whether or not the transition meets the condition, we update
        % cond_cell and priors
        cond_cell{i} = poss_trans;
    end
end

end

function d = lag(sp1, sp2, params)
d = sp2 - sp1;
end