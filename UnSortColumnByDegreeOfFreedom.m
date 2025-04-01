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

function [Sold] = UnSortColumnByDegreeOfFreedom(sortA,ind)

% This function Unsorts a sorted Matrix by
% [K,IX]=SortByDegreeOfFreedom(S,KnownDispVec) with
% Sold=UnsorByDegreeOfFreedom(K,IX);

unsorted = 1:length(sortA);

newInd(ind) = unsorted;

sortedS1 = sortA(newInd,:);


Sold=sortedS1;


end