function events = append_salience(events)
% APPEND_SALIENCE updates event information to specify the
% relative (to the rest of the image) salience for each
% saccade or fixation that occurs
%
% INPUTS: events - structure, contains behavioral and eye events
%
% OUTPUTS: events - structure, updated with visual salience per dgii

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


model_dir = '..\..\dgii_model\';

fix_idx = find(strcmp({events.type},'fixation'));
for n = 1:length(fix_idx)
   
   this_model = [model_dir sprintf('%012d',events(fix_idx(n)).stim_idx) '.h5'];
   this_model = h5read(this_model, '/dgii_lpd');
   this_model = this_model';
   this_model = imresize(this_model, [960, 1280]);
   
    %%% Below is only necessary for display
    % this_image = [stim_dir sprintf('%012d',fix_e(n).stim_idx) '.jpg'];
    % this_image = imread(this_image);
    % 
    % model_full = nan(1080,1920);
    % model_full(61:end-60,171:end-170) = imresize(this_model',2);
    %    
    % image_full = uint8(175*ones(1080,1920,3));
    % image_full(61:end-60,171:end-170,:) = imresize(this_image,2);
   
   this_xy = [events(fix_idx(n)).start_x; events(fix_idx(n)).start_y];
   this_xy = round(this_xy);
   this_xy(this_xy<=0)=1;
   
   sal_s = nan(1, size(this_xy,2));
   for e = 1:size(this_xy,2)
       if this_xy(1,e)>1920 || this_xy(2,e)>1080 % outside of image
           continue
       end
       
       [~, fixmap] = gaze2binary(this_xy(:,e));
       if ~any(fixmap==1) % outside of image
           continue
       end
       sal_s(e) = this_model(fixmap==1); % should only be one value
       
   end
   events(fix_idx(n)).sal_s = sal_s;
end

% add for saccades
sac_idx = find(strcmp({events.type},'saccade'));

for n = 1:length(sac_idx)
   
   this_model = [model_dir sprintf('%012d',events(sac_idx(n)).stim_idx) '.h5'];
   this_model = h5read(this_model, '/dgii_lpd');
   this_model = this_model';
   this_model = imresize(this_model, [960, 1280]);
   
    %%% Below is only necessary for display
    % this_image = [stim_dir sprintf('%012d',fix_e(n).stim_idx) '.jpg'];
    % this_image = imread(this_image);
    % 
    % model_full = nan(1080,1920);
    % model_full(61:end-60,171:end-170) = imresize(this_model',2);
    %    
    % image_full = uint8(175*ones(1080,1920,3));
    % image_full(61:end-60,171:end-170,:) = imresize(this_image,2);
   
    % endpoint of the saccade
   this_xy = [events(sac_idx(n)).end_x; events(sac_idx(n)).end_y];
   this_xy = round(this_xy);
   this_xy(this_xy<=0)=1;
   
   sal_s = nan(1, size(this_xy,2));
   for e = 1:size(this_xy,2)
       
       if this_xy(1,e)>1920 || this_xy(2,e)>1080 % outside of image
           continue
       end
       
       [~, fixmap] = gaze2binary(this_xy(:,e));
       if ~any(fixmap==1) % outside of image
           continue
       end
       sal_s(e) = this_model(fixmap==1); % should only be one value
       
   end
   events(sac_idx(n)).sal_s = sal_s;
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


