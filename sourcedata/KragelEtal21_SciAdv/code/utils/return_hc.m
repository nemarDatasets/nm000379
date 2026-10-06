function is_hc = return_hc(subjid, label)
% RETURN_HC looks up whether label (pairs or mono) are in hippocampus
%
% INPUTS:  label - string, contains contact names
%
% OUTPUTS: is_hc - logical, specifying where/if the label is in hippocampus
%

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

if isequal(subjid, 'S1')
    contacts = {'B1','B2','B3','C1','D1','AL''1'};
elseif isequal(subjid, 'S2')
    contacts = {'*B1','*B2', ...
                '*C1','*C2', ...
                '*D1','*D2', ...
                'B1','B2', ...
                'C1','C2', ...
                'D1','D2'};
elseif isequal(subjid, 'S3')
    contacts = {'B1','B2','C1','C2','D1','D2'};
elseif isequal(subjid, 'S4')    
    contacts = {'C1','C2','D1','D2'};
elseif isequal(subjid, 'S5')
    contacts = {'B1','B2','D1','D2'};
elseif isequal(subjid, 'S6')
    contacts = {'B1','B2','B3','C1','C2'};
end

is_hc = false(size(label));
for i = 1:length(is_hc)
    for c = 1:length(contacts)
        if startsWith(deblank(label{i}), [contacts{c} '-'])
            is_hc(i) = true;
        elseif endsWith(deblank(label{i}), ['-' contacts{c}])
            is_hc(i) = true;
        elseif startsWith(deblank(label{i}), [contacts{c}]) && endsWith(deblank(label{i}), [contacts{c}])
            is_hc(i) = true;
        end
    end
end
    

