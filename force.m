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


clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force based frame element with no element iterations
% Only Structural iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this example contains only on e element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-d coordinateds of node 1 and 2
% node 1 is (x1,y1,z1) node 2 is (x2,y2,z2)
x1=0;
y1=0;
z1=0;
x2=6;
y2=0;
z2=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha001 is the orientation angle with respect to the default orientation
alpha001=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G is shear modulus, I1 is the respective inertia to find section rigidity
% G*I1 (in OpenSees GJ)
G=1;
I1=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a minimum amount of strain (trial strain) to get stress from fibers
% elasticity (the strain should be small for linear elastic respnse
strainlin=0.0000001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for this example
% number of integration points for element is 6
numintpoint=6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for Newton-Raphson iterations
% number of increment is 10
% number of iterations is 25
% incremental-iterative formulation
numinc=10;
numiter=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section description
% A is area of fibers
% there are four fibers with 0.01 area
A=[0.01;0.01;0.01;0.01];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p is structural deformation vector, here it is initialized with zeros for
% increment calculations
p=zeros(12,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross section axes y and z shown, this is a right-handed axis, you can
% determine x acording to y and z
%         y
          %
          %          
      %%%%%%%
      %  %  %
% z %%%%%%  %   
      %  %  %
      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fiber layout (their center according to y and z axes
yy1=-0.05;
zz1=+0.05;
yy2=-0.05;
zz2=-0.05;
yy3=+0.05;
zz3=-0.05;
yy4=+0.05;
zz4=+0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I is the cross section definition matrix
I=[-yy1 zz1 1;
   -yy2 zz2 1;
   -yy3 zz3 1;
   -yy4 zz4 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lengthI=length(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L is length of frame 
L=((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T is the matrix that helps finding the 10 joint forces of element (except
% the torsional forces (2 force in each joint) these will total 12 element
% force, 6 in 1st joint, 6 in 2nd joint)

T=[0 1/L 0 0 1 0 -1/L 0 0 0;
   0 1/L 0 0 0 0 -1/L 0 0 1;
   0 0 -1/L 1 0 0 0 1/L 0 0;
   0 0 -1/L 0 0 0 0 1/L 1 0;
  -1 0 0 0 0 1 0 0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TransM is the transformation matrix for local to global or global to
% local transformations
TransM=TransMatrix(x1,y1,z1,x2,y2,z2,alpha001);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J is Jacobian of element
J=L/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detJ=det(J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fele is Flexibility matrix of element, here it is initialized with zeros
% for future calculations
FEle=zeros(5,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lgwt is Gauss-Legendre quadrature function
% by Greg von Winckel (2023). Legendre-Gauss Quadrature Weights and Nodes 
% (https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes),
% MATLAB Central File Exchange. Retrieved March 3, 2023.
% -1 and +1 are isoparametric element's end locations
% intpoint is integration points
% weight is weigths of integration point locations
[intpoint,weight]=lgwt(numintpoint,-1,+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first, element's linear elastic response, here we integrate according to
% integration points and weights below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ai=1:numintpoint
r=intpoint(ai,1);
w=weight(ai,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b is force interpolation function
b(:,:,ai)=[1/2*(r-1) 1/2*(r+1) 0 0 0;
   0 0 1/2*(r-1) 1/2*(r+1) 0;
   0 0 0 0 1];    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for trial strain defined above, we calculate fiber's elasticity modulus,
% in here this should be linear elastic elasticity modulus
for fb=1:lengthI
E(fb,ai)=(210000000000*0.05*strainlin+(210000000000-210000000000*0.05)*strainlin/(1+(strainlin/(355000000/210000000000))^20)^(1/20))/strainlin;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k is section stiffness
k(:,:,ai)=transpose(I)*(diag(E(:,ai).*A))*I;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEle is element flexibility matric
% FEle is calculated accoding to elastic response here since
% we are in elastic range (initial calculation for Newton-Raphson iteration
FEle=FEle+transpose(b(:,:,ai))*k(:,:,ai)^-1*b(:,:,ai)*detJ*w;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KEle is Stiffness of element (inverse of Flexibility)
KEle=FEle^-1;
% Here our Stiffness is 5by5 it will be 10 by 10 after transformation with
% T matrix,
% then we will add torsional values to stiffness to make it 12 by 12
KEle51=KEle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KE is 10 by 10 stiffness matrix
KE=transpose(T)*KEle*T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Knew is 12 by 12 matrix for stiffness (initialization
KNew=zeros(12,12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here at the end, we obtain 12 by 12 stiffness matrix of element
% first we insert KE to KNew
KNew(1:10,1:10)=KE;
% then torsional values to stiffness matrix
KNew(11,11)=+G*I1/L;
KNew(11,12)=-G*I1/L;
KNew(12,11)=-G*I1/L;
KNew(12,12)=+G*I1/L;
% then we arrange stiffness matrix to appropriate degree of freedoms
KStiff=KNew(:,[1 2 3 11 4 5 6 7 8 12 9 10]);
KStiff=KStiff([1 2 3 11 4 5 6 7 8 12 9 10],:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KStiff is the local stiffness matrix (12 by 12)
% TransM is transformation matrix
% Here we transform local stiffness matrix from local to global
KEle=transpose(TransM)*KStiff*TransM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KEle1 is to represent first element in calculations, there may be other
% elements, however here we have only one element
KEle1=KEle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KStr is Structure's stiffness matrix, here we initialized it
% since we have only one element, it is 12 by 12
KStr=zeros(12,12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GlobalAssembly6dof is a functions to assemble global element stiffness
% matrices to the global structural stiffness matrix
KStr=GlobalAssembly6dof(KStr,KEle1,[1 2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P is load vector, it is the external forces applied to the relative
% degrees of freedoms
P=[0;0;0;0;0;0;18500000;1000;1000;1000;1000;1000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PartedSolve is a solver for known displacements
% you may also get reactions from function (not in this example)
% the known displacements for degrees of freedom 1 2 3 4 5 6 are zero to
% show fixed supports
deltap=PartedSolve(KStr,P/numinc,transpose([1 2 3 4 5 6]),transpose([0 0 0 0 0 0]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total structural displacement is optained from difference displacement
% which is deltap
p=p+deltap;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton-Raphson iterations with load increments
% Static,load control analysis
% No element iterations, only structural iterations
% Element stiffness is found in one step
% other elements can be assembled to the structure if you have
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:numinc %number of newton raphson increment
for jj=1:numiter %number of newton raphson iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P is structural force, p is structural displacement
% Q is element forve, q is element displacement
% D is section force (not used here, d is section displacement
% EE is fiber force (stress), e is fiber displacement (strain)
% E is Elasticity modulus of fiber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we have one element, so element displacement is eaqual to structural displacement    
q=p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KEle=KEle51;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q is transformed from global to local
q=TransM*q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qTorsion is torsional values of global displacement of element
% we store them for future use
qTorsion=q([4 10],:);

% we get 10 by 10 displacements
q=q([1 2 3 5 6 7 8 9 11 12],:);

% we get 5 by 5 displacements of local element displacements 
q=T*q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of element force vector
FEle=zeros(5,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of internal resisting force vector
ri=zeros(5,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integtation point and weigth calculations for element
[intpoint,weight]=lgwt(numintpoint,-1,+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integration for element
for ai=1:numintpoint
r=intpoint(ai,1);
w=weight(ai,1);    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b is force interpolation matrix
b(:,:,ai)=[1/2*(r-1) 1/2*(r+1) 0 0 0;
   0 0 1/2*(r-1) 1/2*(r+1) 0;
   0 0 0 0 1];   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d is section deformation
% we found it from Kele*q which is element force
% element force is interpolated by b (b*element force
% we do these two steps by (b*KEle*q)
% then we founf section deformation by multiplying
% inverse of k with section force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is like F=k*x so, x = k^-1*F
% F is force in the above line, k is stiffness, x is displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d(:,ai)=k(:,:,ai)^-1*(b(:,:,ai)*KEle*q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e is fibers' displacement (strain of the fibers
% we find it using I and section displacement d
e(:,ai)=I*d(:,ai);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we used Menegetto-Pinto material model for fibers (first cycle of
% response, we simplified the material model for this example
% we used strains of the fibers (e) to find the elasticities of fibers
% (they are probably nonlinear now this is of cource according to the strain (please consult the material model))
for fb=1:lengthI
E(fb,ai)=(210000000000*0.05*e(fb,ai)+(210000000000-210000000000*0.05)*e(fb,ai)/(1+(e(fb,ai)/(355000000/210000000000))^20)^(1/20))/e(fb,ai);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section stiffness is calculated according to the new elasticity moduluses
% they can be nonlinear now (according to the model of course)
k(:,:,ai)=transpose(I)*(diag(E(:,ai).*A))*I;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ri internal reaction force is calculated here by displacement of section
% times force interpolation matrix
ri=ri+transpose(b(:,:,ai))*d(:,ai)*detJ*w; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Element Flexibility matrix is calculated  by b and k
% we used detJ because of the isoparametric conversion
FEle=FEle+transpose(b(:,:,ai))*k(:,:,ai)^-1*b(:,:,ai)*detJ*w;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% again transformations etc., please look above comments
KEle=FEle^-1;
KEle51=KEle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KE=transpose(T)*KEle*T;
KNew=zeros(12,12);
KNew(1:10,1:10)=KE;
KNew(11,11)=+G*I1/L;
KNew(11,12)=-G*I1/L;
KNew(12,11)=-G*I1/L;
KNew(12,12)=+G*I1/L;
KStiff=KNew(:,[1 2 3 11 4 5 6 7 8 12 9 10]);
KStiff=KStiff([1 2 3 11 4 5 6 7 8 12 9 10],:);
KEle=transpose(TransM)*KStiff*TransM; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KEle1=KEle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local to global transformation etc.
ri=transpose(T)*ri;
ri(11,1)=0;
ri(12,1)=0;
ri=ri([1 2 3 11 4 5 6 7 8 12 9 10],:);
ri=transpose(TransM)*ri;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly of Structural Stiffness
KStr=zeros(12,12);
KStr=GlobalAssembly6dof(KStr,KEle1,[1 2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External force is calculated from incrementation
PExt=(P*ii/numinc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PInt is the resisting force
PInt=ri;
% PInt 1 to 6 is zero because of fixed 1st joint
PInt(1:6,1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the unbalanced force is calculated
PUnb=PExt-PInt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% structural displacement is calculated from Stiffness and Unbalanced force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% it is like x=k^-1*F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=PartedSolve(KStr,PUnb,transpose([1 2 3 4 5 6]),transpose([0 0 0 0 0 0]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U is displcament vector of the element
U(:,ii)=p;
% ForceF is the applied force vector of the element
ForceF(:,ii)=P*ii/numinc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot 2nd joints x force-displacement graph
plot (U(7,:),ForceF(7,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

