function  [auc_enc, auc_rec, auc_center_enc, auc_center_rec] = compute_model_performance(subjid)
% COMPUTE_MODEL_PERFORMANCE computers classification performance (fixated/non-fixated)
% for DGII and centerbias models
%
% INPUTS:  subjid - string, identifier of subject to analyze
%
% OUTPUTS: auc_enc - double, prediction performance (auc) for DGII at study
%
%          auc_rec - double, prediction performance (auc) for DGII at test
%
%          auc_center_enc - double, prediction performance (auc) for centerbias
%                           model at study
%
%          auc_center_rec - double, prediction performance (auc) for centerbias
%                           model at test

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

% path that contains dgii model outputs in hdf5 format
model_dir = '..\..\dgii_model\';

events = get_events(subjid);

fix_e = events(strcmp({events.type}, 'fixation'));     
fix_r = fix_e(strcmp({fix_e.phase},'recog'));
fix_e = fix_e(strcmp({fix_e.phase},'encode'));

res.n_fix = [fix_e.n_fix];
res.rt = [fix_e.rt];
res.sal = [fix_e.salience];
res.resp = [fix_e.resp];
 
res.r_n_fix = [fix_r.n_fix];
res.r_rt = [fix_r.rt];
res.r_sal = [fix_r.salience];
res.r_resp = [fix_r.resp];

% encoding
auc_enc = nan(1, length(fix_e)); 
for n = 1:length(fix_e)
   
   this_model = [model_dir sprintf('%012d',fix_e(n).stim_idx) '.h5'];
   this_model = h5read(this_model, '/dgii_lpd');
       
   this_xy = [fix_e(n).start_x(~fix_e(n).precue); fix_e(n).start_y(~fix_e(n).precue)];
   
  if isempty(this_xy)
    auc_enc(n) = nan;
    outmap_enc(n,:,:) = nan;
    continue
  end
   
  this_xy = round(this_xy);
  this_xy(this_xy<=0)=1;
  this_xy(1,this_xy(1,:)>1920)=1920;
  this_xy(2,this_xy(2,:)>1080)=1080;
   
  [this_full_fixmap, this_fixmap] = gaze2binary(this_xy);
   
  auc_enc(1,n) = AUC_Judd(this_model', this_fixmap, 1, 0);

  outmap(n,:,:) = imgaussfilt(this_fixmap, 40); % for null center model

end

% centerbias for encoding

for n = 1:length(fix_e)
   
  this_xy = [fix_e(n).start_x(~fix_e(n).precue); fix_e(n).start_y(~fix_e(n).precue)];
   
   
  if isempty(this_xy)
     auc_enc(n) = nan;
     outmap_enc(n,:,:) = nan;
     continue
  end
   
  this_xy = round(this_xy);
  this_xy(this_xy<=0)=1;
  this_xy(1,this_xy(1,:)>1920)=1920;
  this_xy(2,this_xy(2,:)>1080)=1080;
   
  [this_full_fixmap, this_fixmap] = gaze2binary(this_xy);
   
  this_center = squeeze(nanmean(outmap(~ismember(1:length(fix_e), n),:,:)));
  this_center = this_center./sum(this_center(:)); % normalize
   
  auc_center_enc(1,n) = AUC_Judd(this_center, this_fixmap, 1, 0);
   
end


% retrieval 
for n = 1:length(fix_r)
   
   this_model = [model_dir sprintf('%012d',fix_r(n).stim_idx) '.h5'];
   this_model = h5read(this_model, '/dgii_lpd');
   
   this_xy = [fix_r(n).start_x(~fix_r(n).precue); fix_r(n).start_y(~fix_r(n).precue)];
   
   if isempty(this_xy)
       auc_rec(n) = nan;
       outmap_rec(n,:,:) = nan;
       continue
   end
   
   this_xy = round(this_xy);
   this_xy(this_xy<=0)=1;
   this_xy(1,this_xy(1,:)>1920)=1920;
   this_xy(2,this_xy(2,:)>1080)=1080;
   
   [this_full_fixmap, this_fixmap] = gaze2binary(this_xy);
   
   auc_rec(n) = AUC_Judd(this_model', this_fixmap, 1, 0);
   outmap_rec(n,:,:) = imgaussfilt(this_fixmap, 40); % for null center model

end

% centerbias for retrieval

for n = 1:length(fix_r)
   
   this_xy = [fix_r(n).start_x(~fix_r(n).precue); fix_r(n).start_y(~fix_r(n).precue)];
   
  if isempty(this_xy)
    auc_center_rec(n) = nan;
    continue
  end
   
  this_xy = round(this_xy);
  this_xy(this_xy<=0)=1;
  this_xy(1,this_xy(1,:)>1920)=1920;
  this_xy(2,this_xy(2,:)>1080)=1080;
   
  [this_full_fixmap, this_fixmap] = gaze2binary(this_xy);
   
  this_center = squeeze(nanmean(outmap_rec(~ismember(1:length(fix_r), n),:,:)));
  this_center = this_center./sum(this_center(:)); % normalize
   
   auc_center_rec(1,n) = AUC_Judd(this_center, this_fixmap, 1, 0);
   
end

end

function [full_fixmap, fixmap] = gaze2binary(xy)

% full screen
fixmap = zeros(1080,1920);

for i = 1:size(xy,2)
fixmap(xy(2,i),xy(1,i)) = 1;
end

% and the image
full_fixmap = fixmap;
fixmap = fixmap(61:end-60,321:end-320);

end

