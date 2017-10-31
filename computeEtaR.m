function etaR = computeEtaR(Sol,fsT,k)

%computeEtaR: computes residual-based error estimator 
%
%Usage:
%
%etaR = computeEtaR(Sol,f,k)
%
%Comments:
%    Sol contains the mesh data and the solution data. f is the vector of the 
%    residual of each element. k is the diffusion coefficient. For the Allen-Cahn
%    equation, k = solver.epsilon^2, where solver.epsilon is the interfacial thickness
%    parameter.
%
%    The function returns the column vector etaR where etaR(J) is the
%    squared error indicator associated with the j-th element. These values
%    may be used to mark triangles for refinement. In particular, the 
%    value of the residual error estimator is given by sqrt(sum(etaR)).
%    
%Remark:
%
%    This program is a supplement to the paper 
%    >> Efficient Implementation of Adaptive P1-FEM in Matlab <<
%    by S. Funken, D. Praetorius, and P. Wissgott. The reader should 
%    consult that paper for more information.   
%
%Authors:
% 
%    S. Funken, D. Praetorius, P. Wissgott  10-07-08
%
%Modified by:
%
%    Vaibhav Joshi, Rajeev K. Jaiman
%
%

elements = Sol.elem ;
x = Sol.phiAlpha ;
coordinates = Sol.node ;
dirichlet2edges = Sol.Dirichlet ;
neumann = Sol.Neumann ;
edge = [Sol.elem(:,[1,2]); Sol.elem(:,[1,3]); Sol.elem(:,[2,3])];
edge2nodes = unique(sort(edge,2),'rows');

element2edges = [sort(Sol.elem(:,[1 2]),2) sort(Sol.elem(:,[1 3]),2) sort(Sol.elem(:,[2 3]),2)] ;
[~,indx1] = ismember(element2edges(:,1:2),edge2nodes,'rows');
[~,indx2] = ismember(element2edges(:,3:4),edge2nodes,'rows');
[~,indx3] = ismember(element2edges(:,5:6),edge2nodes,'rows');
clear element2edges
element2edges = [indx1 indx3 indx2] ;

%*** First vertex of elements and corresponding edge vectors
c1  = coordinates(elements(:,1),:);
d21 = coordinates(elements(:,2),:) - c1;
d31 = coordinates(elements(:,3),:) - c1;
%*** Vector of element volumes 2*|T|
J = [d21(:,1), d31(:,1), ...
     d21(:,2), d31(:,2)] ;
volume = ( J(:,1).*J(:,4) -J(:,2).*J(:,3) ) ;
Nx = [-1 ; 1 ; 0 ];
Ny = [-1 ; 0 ; 1 ];
nen = 3 ;
DNDx = ((J(:,4))*Nx'+ (-J(:,3))*Ny')./repmat(volume,1,nen);
DNDy = ((-J(:,2))*Nx'+(J(:,1))*Ny')./repmat(volume,1,nen);

dPhidx = DNDx(:,1).*x(elements(:,1)) + DNDx(:,2).*x(elements(:,2)) + DNDx(:,3).*x(elements(:,3)) ;
dPhidy = DNDy(:,1).*x(elements(:,1)) + DNDy(:,2).*x(elements(:,2)) + DNDy(:,3).*x(elements(:,3)) ;
% %*** Compute curl(uh) = (-duh/dy, duh/dx)
curl = [-dPhidy dPhidx] ;

%*** Compute edge terms hE*(duh/dn) for uh
dudn21 = sum(d21.*curl,2);
dudn13 = -sum(d31.*curl,2);
dudn32 = -(dudn13+dudn21);
etaR = accumarray(element2edges(:),[dudn21;dudn32;dudn13],[size(edge2nodes,1) 1]);

%*** Incorporate Neumann data
%%% Since the Neumann boundary condition on the order parameter phi is zero, this is not needed !!

%*** Incorporate Dirichlet data
etaR(dirichlet2edges) = 0;
etaR = k.*etaR ;
%*** Assemble edge contributions of indicators
etaR = sum(etaR(element2edges).^2,2);
%*** Add volume residual to indicators
% fsT = feval(f,(c1+(d21+d31)/3));
etaR = etaR + (0.5*volume.*fsT).^2 ;
