function [res, lags, sal, dist] = compute_dist_by_lag(events)
% COPMUTE_DIST_BY_LAG for a subject (events), computes distance by lag
%
% INPUTS:  events - structure, contains event information
%
% OUTPUTS: res - structure, contains subject medians and standard errors
%                for all measures
%
%          lags - double, transition lags for all test period fixations
%
%          sal - double, salience for all test period fixations
%
%          dist - double, distance between test period fixations

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


% filter to only retrieval fixation events
rec_fix_events = events(strcmp('fixation', {events.type}) & ...
    strcmp('recog',{events.phase}));

% loop over events and append distance from previous saccade

close_flag = false; % close paths between 'recalled' locations
for e = 1:length(rec_fix_events)
    rec_fix_events(e).dist = nan;
    rec_fix_events(e).sal_s = [nan rec_fix_events(e).sal_s(1:end-1)]; % shift by one for salience of next item
    rec_fix_events(e) = append_distance(rec_fix_events(e), close_flag);
end

lags = [rec_fix_events.lag];
sal  = [rec_fix_events.sal_s];
dist = [rec_fix_events.dist];
precue = [rec_fix_events.precue];

tr = precue==1; % exclude fixations in the pre-cue period
ul = -4:4;

res.dbl_y = nan(1, length(ul));
res.dbl_e = nan(1, length(ul));
res.n = nan(1, length(ul));

for i = 1:length(ul)
   res.dbl_y(i) = nanmedian(dist(lags==ul(i))); 
   res.dbl_e(i) = nanstd(dist(lags==ul(i)))./sqrt(sum(lags==ul(i)));
   
   res.sbl_y(i) = nanmedian(sal(lags==ul(i)));
   res.sbl_e(i) = nanstd(sal(lags==ul(i)))./sqrt(sum(lags==ul(i)));
   
   res.n(i) = sum(lags==ul(i));
end

end

function event = append_distance(event, close_flag)

if ~isequal(length(event.precue), length(event.lag))
    event.precue = [];
end

dist = nan(size(event.lag));
for f  = 1:length(event.lag)
    if ~isnan(event.lag(f))
        
        if ~close_flag
            dist(f) = sqrt( (event.mean_x(f) - event.mean_x(f+1))^2 + ...
                (event.mean_y(f) - event.mean_y(f+1))^2 );
        else
            tmp = find(find(~isnan(event.lag))>f,1,'first');
             dist(f) = sqrt( (event.mean_x(f) - event.mean_x(f+tmp))^2 + ...
                (event.mean_y(f) - event.mean_y(f+tmp))^2 );
        end

    end
end

event.dist = dist;

end