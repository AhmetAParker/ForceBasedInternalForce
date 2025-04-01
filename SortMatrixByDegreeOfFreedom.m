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

function [K,IX] = SortMatrixByDegreeOfFreedom(S,KnownDispVec)

% This function Sorts a Matrix S to a given Known Displacement Vector
% as KnownDispVec
% Example SortByDegreeOfFreedom([1:5;6:10;11:15;16:20;21:25],[3 5])


sofS=length(S);

KnownDispVector=zeros(sofS,1);

for aa=KnownDispVec
    KnownDispVector(aa)=1;
end


[SA,IX] = sort(KnownDispVector);

sortedS1 = S(IX,:);

sortedS2 = sortedS1(:,IX);

K=sortedS2;
end