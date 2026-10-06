function [psi_revisit, psi_other, psi_diff, anat, label] = get_psi_by_network(subj)
% GET_PSI_BY_NETWORK loads psi data and organizes by connection
%
% INPUTS:  subj - string, specifies subject id to load data 
%
% OUTPUTS: psi_revisit - double, average psi_z for each connection around revisitation fixations
%
%          psi_other - double, average psi_z for each connection around other fixations
%
%          psi_diff - double, z-scored difference in psi between revisit and other fixations
%
%          anat - structure specifying anatomical information for channel pairs
%
%          label - cell array of strings, labels for channel pairs

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

load( ['..\..\scratch_data\connectivity\' subj '_psi_wvis_base.mat'], 'res');

label = cat(2, res.label)';

n_conn = length(res);

anat.dan_con = false(1,n_conn);
anat.vis_con = false(1, n_conn);

psi_revisit = nan(n_conn, 751);
psi_other = nan(n_conn, 751);
psi_diff = nan(n_conn, 751);


for i = 1:n_conn
   
    % find out of the connection is to DAN or VIS network
    
    is_hc = return_hc(subj, res(i).label);
    [is_dan, is_vis, pairs] = return_yeo(subj, res(i).label);
    
    is_dan(is_hc) = false;
    is_vis(is_hc) = false;
    
    if is_hc(1) % to make sure direction of psi is consistent across all pairs
        sign_flip = false;
    else
        sign_flip = true;
    end
    
    % assign which type of pair this is, for later summaries
    if any(is_dan)
        anat.dan_con(i) = true;
    end
    
    if any(is_vis)
        anat.vis_con(i) = true;
    end
    
    if ~sign_flip
        psi_revisit(i,:) = res(i).psi_z1;
        psi_other(i,:) = res(i).psi_z2;
        psi_diff(i,:) = res(i).psi_diff;
    else
        psi_revisit(i,:) = -res(i).psi_z1;
        psi_other(i,:) = -res(i).psi_z2;
        psi_diff(i,:) = -res(i).psi_diff;
    end

end

end
