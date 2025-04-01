%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ForceBasedInternalForce
% Copyright (C) 2023  Ahmet Alper Parker
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TR] = TransMatrix (x1,y1,z1,x2,y2,z2,alfa)

beta= atan2(z2-z1,sqrt((x2-x1)^2+(y2-y1)^2));
gamma=atan2(y2-y1,x2-x1);



TTR=[ cos(beta)*cos(gamma)   cos(beta)*sin(gamma)   sin(beta);
      sin(gamma)*sin(alfa)-cos(gamma)*sin(beta)*cos(alfa)  -sin(alfa)*cos(gamma)-cos(alfa)*sin(gamma)*sin(beta)  cos(alfa)*cos(beta);
      sin(gamma)*cos(alfa)+cos(gamma)*sin(beta)*sin(alfa)   -cos(alfa)*cos(gamma)+sin(gamma)*sin(beta)*sin(alfa)  -sin(alfa)*cos(beta)];



TR=[TTR zeros(3,3) zeros(3,3) zeros(3,3);
    zeros(3,3) TTR zeros(3,3) zeros(3,3);
    zeros(3,3) zeros(3,3) TTR zeros(3,3);
    zeros(3,3) zeros(3,3) zeros(3,3) TTR];
