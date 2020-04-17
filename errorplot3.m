%% Error plot3
% Plots semi-transparent error in current figure

%  Copyright 2019 Philip Jean-Richard-dit-Bressel, UNSW Sydney

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function errorplot3(LCI, UCI, xlims, col, alp)

if ~isempty(xlims)
   xl = linspace(xlims(1),xlims(2),length(LCI));
else
   xl = 1:length(LCI);
end

xu = flip(xl);
x = [xl xu];

y = [LCI flip(UCI)];
patch(x,y, col, 'linestyle', 'none','FaceAlpha', alp);