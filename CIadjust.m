function [adjLCI,adjUCI] = CIadjust(LCI,UCI,est,n,adj_type)
%% Adjusts CI according to adj_type
   % Type 1 = extend CI from mean by sqrt((n)/(n-1))
      % needs all inputs
   % Type 2 = expand CI by sqrt((n)/(n-1))
      % doesn't need mean

%%Inputs:
% LCI = lower confidence limit vector
% UCI = upper confidence limit vector
% est = parameter estimate  vector (eg. mean, only necessary for Type 1)
% n = experimental n
% adj_type = adjustment type (see above)

%%Output:
% adjLCI = adjusted LCI
% adjUCI = adjusted UCI

%% By Philip Jean-Richard-dit-Bressel, UNSW Sydney, 2020
% Feel free to use with citation: Jean-Richard-dit-Bressel et al. (2020). https://doi.org/10.3389/fnmol.2020.00014

%% GNU
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

%%
if adj_type == 1
   %% Extend CI from mean
   CI_fix = sqrt((n)/(n-1));
   fprintf(['CI extended from mean by ' num2str(CI_fix*100) 'pc\n']);
   adjUCI = (UCI-est).*CI_fix + est;
   adjLCI = est - (est-LCI)./CI_fix; 
   
elseif adj_type == 2
   %% Expand CI
   CI_fix = sqrt((n)/(n-1));
   fprintf(['CI expanded by ' num2str(CI_fix*100) 'pc\n']);
   CIchange = ((UCI - LCI).*CI_fix - (UCI - LCI))/2; %/2
   adjUCI = UCI+CIchange;
   adjLCI = LCI-CIchange;
end
