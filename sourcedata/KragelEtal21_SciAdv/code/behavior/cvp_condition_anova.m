function [stats, summ] = cvp_condition_anova(cvp, pos)
% CVP_CONDITION_ANOVA runs a repeated measures anova with Direction and Lag
% as factors
%
% INPUTS:  cvp - double, conditional viewing probabilites from lags 
%                -4 to 4
%
%           pos - double, lags to include in analysis 
%
% OUTPUTS: stats - structure, from mes toolbox reporting summary 
%                  statistics and effect sizes
%
%          summ  - cell array, table describing anova results

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

Y = [];
X = [];

ref_lags = -4:4;
lag_idx = find(ref_lags == pos);

for s = 1:size(cvp,1)
    for cond = 1:size(cvp,2)
        cvp_tmp = cvp(s,cond,:);
        
        if ~isnan(cvp_tmp(lag_idx))
            Y = cat(1, Y, cvp_tmp(lag_idx));
            
            if cond == 1
                X = cat(1, X, [1 1]);
            elseif cond == 2
                X = cat(1, X, [1 0]);
            elseif cond == 3
                X = cat(1, X, [0 1]);
            else
                X = cat(1, X, [0 0]);
            end
            
        end
    end
end

[X, i] = sortrows(X);
Y = Y(i);

[stats, summ] = mes2way(Y, X, 'partialeta2', ...
    'isDep',[1 1], ...
    'fName',{'Salience', 'Memory'}, ...
    'nBoot', 1000);

end

