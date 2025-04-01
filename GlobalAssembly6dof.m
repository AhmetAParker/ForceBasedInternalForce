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

function [ S ] = GlobalAssembly6dof( S,k,nodeids )

dim=length(nodeids);
 
m=zeros(1,6*dim);
n=zeros(1,6*dim);

e=0;
f=0;

m=[6*nodeids(1)-5 6*nodeids(1)-4 6*nodeids(1)-3 6*nodeids(1)-2 6*nodeids(1)-1 6*nodeids(1)];
n=[6*nodeids(1)-5 6*nodeids(1)-4 6*nodeids(1)-3 6*nodeids(1)-2 6*nodeids(1)-1 6*nodeids(1)];


for aa=2:dim;

m=[m 6*nodeids(aa)-5 6*nodeids(aa)-4 6*nodeids(aa)-3 6*nodeids(aa)-2 6*nodeids(aa)-1 6*nodeids(aa)];
n=[n 6*nodeids(aa)-5 6*nodeids(aa)-4 6*nodeids(aa)-3 6*nodeids(aa)-2 6*nodeids(aa)-1 6*nodeids(aa)];


end


    
    for kk=m
e=e+1;
for lll=1:length(n)
    for ll=n(lll)
    f=f+1;
    S(kk,ll)=S(kk,ll)+k(e,f);
    end
end
f=0;
end

end

