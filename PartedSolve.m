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

function [Displacements,Reactions] = PartedSolve(K,F,Restraints,D21)

LK=length(K);
LR=length(Restraints);
LF=length(F);


[Knew1,~]=SortMatrixByDegreeOfFreedom(K,Restraints);


K11=Knew1(1:LK-LR,1:LK-LR);
K12=Knew1(1:LK-LR,LK-LR+1:LK);
K21=Knew1(LK-LR+1:LK,1:LK-LR);
K22=Knew1(LK-LR+1:LK,LK-LR+1:LK);


[Fnew1,IXF1]=SortColumnByDegreeOfFreedom(F,Restraints);


F11=Fnew1(1:LF-LR,1);
% F21=Fnew1(LF-LR+1:LF,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D11=K11^-1*(F11-K12*D21);
% % Solution of Equations with Incomplete Cholesky Factorization with
% % preconditioned conjugate gradient method
% K11sparse=sparse(K11);
% icholK11=ichol(K11sparse);
% D11=pcg(K11,(F11-K12*D21),1e-8,1000,icholK11,icholK11');
% %%%%%%%%%%%%%%%%%%%%%%5%%%%%

F11new=zeros(length(F11),1);
F21new=(K21*D11+K22*D21);

Displacements=UnSortColumnByDegreeOfFreedom([D11;D21],IXF1);
Reactions=UnSortColumnByDegreeOfFreedom([F11new;F21new],IXF1);


end

