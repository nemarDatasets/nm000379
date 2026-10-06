function [stats, summ] = cvp_anova(cvp)
% CVP_ANOVA runs a repeated measures anova with Direction and Lag
% as factors
%
% INPUTS:  cvp - double, conditional viewing probabilites from lags 
%                -4 to 4
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
for s = 1:size(cvp,1)
    
    cvp_tmp = cvp(s,:);
        
    if ~any(isnan(cvp_tmp([1 2 3 4 6 7 8 9])))
        for l = [1 2 3 4 6 7 8 9]
            
            Y = cat(1, Y, cvp_tmp(l));
            X = cat(1, X, [sign(ref_lags(l)) abs(ref_lags(l))]);
        end
    end
end

[X, i] = sortrows(X);
Y = Y(i);


[stats, summ] = mes2way(Y, X, 'partialeta2', ...
    'isDep',[1 1], ...
    'fName',{'Direction', 'Lag'}, ...
    'nBoot', 1000);

end

