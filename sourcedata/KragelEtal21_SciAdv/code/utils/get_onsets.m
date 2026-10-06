function [enc_onsets, rec_onsets, rec_adjust] = get_onsets(subj, base_dir, eeg_dir, bp_flag)
% GET_ONSETS this function gets onsets in samples for encoding and
% recognition trials
%
% INPUTS:  subj - string, subject identifier
%
%          base_dir - string, path to scratch data directory
%
%          eeg_dir - string, path to eeg directory
%
%          bp_flag - logical, flag specifying whether to look at bipolar or monopolar contacts
%
% OUTPUTS: enc_onsets - double, array of onsets for each encoding trial
%
%          rec_onsets - double, array of onsets for each retrieval trial
%
%          rec_adjust - double, array of adjustments to make for retrieval trials for 
%                       different experiment versions (S1)

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

if ~exist('base_dir', 'var')
    base_dir = '..\..\scratch_data\';
end

if ~exist('eeg_dir', 'var')
    eeg_dir = '..\..\data\eeg\';
end

if ~exist('bp_flag', 'var')
    bp_flag = true;
end

onset_file = [base_dir 'onsets\' subj '_onsets.mat'];

if exist(onset_file, 'file')
    load(onset_file, 'enc_onsets', 'rec_onsets', 'rec_adjust');
else
    % compute and save
    tmp_dir = [eeg_dir filesep subj filesep];
    if bp_flag
        files = dir([tmp_dir subj '*_bp.h5']);
    end
    
    n_blocks = length(files);
    enc_onsets = nan(8,24); rec_onsets = nan(8,48);
    for i = 1:8
        for  n = 1:n_blocks
            if contains(files(n).name, ['block' num2str(i)])
                fname = [files(n).folder filesep files(n).name];
                [a, b] = get_trial_onsets(fname, bp_flag);
                enc_onsets(i,1:length(a)) = a;
                rec_onsets(i,1:length(b)) = b;
            end
        end
    end
    
    % and adjustment for recognition trial onsets for S1
    if isequal(subj, 'S1')
        [~, rec_adjust] = adjust_S1_events;
    else
        rec_adjust = zeros(size(rec_onsets));
    end
    
    % missing encoding data for half of block 4, nan onsets
    if isequal(subj, 'S6')
        enc_onsets(4,:) = nan;
        rec_onsets(4,:) = nan;
    end
    
    save(onset_file, 'enc_onsets', 'rec_onsets', 'rec_adjust');
end

end

function [enc_onsets, rec_onsets] = get_trial_onsets(fname, bp_flag)
% GET_TRIAL_ONSETS, finds encoding and retrieval onsets from an h5 file
% with iEEG data

if bp_flag
    pairs = deblank(h5read(fname, '/pairs'));
    idx = strcmp(pairs, 'SYNC');
else
    channels = deblank(h5read(fname, '/channels'));
    idx = strcmp(channels, 'SYNC');
end

% read sync data
sync_data = h5read(fname, '/eeg', [find(idx) 1], [1 Inf]);

% find event starts
[~,locs]=findpeaks(sync_data,'MinPeakProminence',2);
onsets=locs-1;  %take point 1 samples before peak

% both of these are pre-cues (fixation, which start 750ms before
% the image appears)
enc_idx = 3:2:49; % 51 end encoding, 52,53 start rec
rec_idx = 54:1:101; % 102 end rec

enc_idx(enc_idx>length(onsets))=[]; % if run stops early
rec_idx(rec_idx>length(onsets))=[]; % if run stops early

enc_onsets = onsets(enc_idx);
rec_onsets = onsets(rec_idx);

end