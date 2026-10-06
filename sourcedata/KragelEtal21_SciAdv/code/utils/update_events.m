function update_events()
% UPDATE_EVENTS helper function to update events for all subjects and
% write to disk
%
% INPUTS: none, loads events from disk and updates them
%
% OUTPUTS: none, loads events from disk and updates them

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

subjects = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6'};
for s = 1:length(subjects)
    get_events(subjects{s}, true);
end

end

