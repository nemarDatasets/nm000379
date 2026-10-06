function is_contact = return_contact(label, contacts)
% RETURN_CONTACT looks up where label (pairs or mono) are in the contacts
%
% INPUTS:  label - string, pattern looking for within contacts
%
%          contacts - cell array with strings for contacts IDs
%
% OUTPUTS: is_contact - logical, specifying where/if the label is in contacts
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

is_contact = false(size(label));
for i = 1:length(is_contact)
    for c = 1:length(contacts)
        if startsWith(deblank(label{i}), [contacts{c} '-'])
            is_contact(i) = true;
        elseif endsWith(deblank(label{i}), ['-' contacts{c}])
            is_contact(i) = true;
        elseif startsWith(deblank(label{i}), [contacts{c}]) && endsWith(deblank(label{i}), [contacts{c}])
            is_contact(i) = true;
        elseif contains(contacts{c}, deblank(label{i}))
            is_contact(i) = true;
        end
    end
end
    

