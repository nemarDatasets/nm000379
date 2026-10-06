function [lag_p, sal_p] = compute_revisit_ranks(events)
% COMPUTE_REVISIT_RANKS computes the rank order of revisitation
% fixations based on the temporal lag or salience of previously
% viewed locations
%
% INPUTS:  events - structure, contains fixations information
%
% OUTPUTS: lag_p - double, average rank of revisitations according to lag
%          
%          sal_p - double, average rank of revisitations according to salience

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

% filter to only encoding fixation events
enc_fix_events = events(strcmp('fixation', {events.type}) & ...
    strcmp('encode',{events.phase}));

for i = 1:length(enc_fix_events)
    % compute percent of enc_repeated == 0 (i.e., second time viewing the
    % same location), that are replayed at encoding (enc_replay==1)

    seq = enc_fix_events(i).enc_seq;
    enc_repeated = enc_fix_events(i).enc_repeated;
    sal = enc_fix_events(i).sal_s;
    precue = enc_fix_events(i).precue;
        
    repeated = find(enc_repeated==0); % when they returned to previous location
    
    sal_rank = nan(1, length(repeated));
    lag_rank = nan(1, length(repeated));
    for s = 1:length(repeated)
        
        if ~isnan(sal(repeated(s))) 
        sal_rank(s) = percentile_rank(sal(repeated(s)),sal(1:repeated(s)));
        end
                
        u_locs = unique(seq(1:repeated(s)));
        poss_dist = nan(1, length(u_locs)); % lags from current item to all others
        for u = 1:length(u_locs)
            
            if seq(repeated(s)) == u_locs(u)
                poss_dist(u) = repeated(s) - find(seq(1:(repeated(s)-1)) == u_locs(u), 1, 'last');
            else
                poss_dist(u) = repeated(s) - find(seq(1:repeated(s)) == u_locs(u), 1, 'last');
            end
        end
        
        act_dist = poss_dist(u_locs == seq(repeated(s)));
        lag_rank(s) = percentile_rank(-act_dist, -poss_dist);
    end
    
    if ~isempty(lag_rank)
        lag_p(i) = nanmean(lag_rank);
         sal_p(i) = nanmean(sal_rank);
    else
        lag_p(i) = nan;
        sal_p(i) = nan;
    end  
   
    
end

% average over trials
lag_p = nanmean(lag_p);
sal_p = nanmean(sal_p);
