%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %	
%	NONLINEAR ADAPTIVE VARIATIONAL PARTITIONED (NAVP) SOLVER TO MODEL     %
%	 TWO-PHASE FLOWS USING PHASE-FIELD METHOD EMPLOYING CONSERVATIVE      %
%	                    ALLEN-CAHN EQUATION                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	The  current  matlab  solver  solves  the  coupled  incompressible    %
%	Navier-Stokes  equations  in   the   Petrov-Galerkin   variational    %
%	framework  with the phase-field  modeling  employing  conservative    %
%	Allen-Cahn equation to  model  the  two-phase  flows for a general    %
%	unstructured mesh using  an adaptive  algorithm. The details about    %
%	the formulation can be found in the following research articles:      %
%                                                                         %
%	(1) Vaibhav Joshi, Rajeev K. Jaiman. "A  positivity preserving and    %
%	    conservative  variational  scheme  for phase-field modeling of    %
%	    two-phase  flows."  Journal  of  Computational  Physics (Under    %
%	    review). https://arxiv.org/abs/1710.09831                         %
%                                                                         %
%	(2) Vaibhav Joshi,  Rajeev  K. Jaiman. "An  adaptive  unstructured    %
%	    procedure via positivity preserving variational method for the    %
%	    conservative  Allen-Cahn  phase-field  model."  Journal of Co-    %
%	    mputational Physics (Under review).                               %
%	    https://arxiv.org/abs/1710.10406                                  %
%                                                                         %
%	Moreover,  the  adaptive  algorithm  adopted  in  the  code is the    %
%	newest  vertex  bisection method  with coarsening of Chen et. al.,    %
%	the details of which can be found in:                                 %
%                                                                         %
%	(3) L. Chen, C. Zhang.  "AFEM@MATLAB: A MATLAB package of adaptive    %
%	    finite  element  methods."  Technical  report.  University  of 	  %
%	    Maryland.                                                         %
%                                                                         %
%	(4) L. Chen, C. Zhang. "iFEM: An innovative finite  element method    %
%	    package in MATLAB." Technical report. University of California 	  %
%	    at Irvine.                                                        %
%                                                                         %
%	(5) S. Funken, D. Praetorius, P. Wissgott. " Efficient implementa-	  %
%	    tion of adaptive P1-FEM in Matlab." Computational  Methods  in	  %
%	    Applied Mathematics. Vol. 11 (2011), No. 4, pp. 460-490.		  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 	
%	Some important comments about the code:                               %
%                                                                         %
%	(1) The code consists of the following files:				          %
%	    - multiPhaseFem_NAVP.m : The main function consis-                %
%	      ting of the complete variational incompressible Navier-Stokes   %
%	      and Allen-Cahn solvers including the time loop, the non-linear  %
%	      iteration loop and the convergence criteria.			          %
%	    - computeEtaR.m : This code has been modified from reference (5). %
%	      It evaluates the error indicator based on the residual error	  %
%	      estimates of the Allen-Cahn equation.				              %
%	    - bisection.m : Modified from reference (3), this refined  the	  %
%	      mesh using the newest vertex bisection algorithm.			      %
%	    - coarsening.m : Modified from reference (3), this coarsens the	  %
%	      mesh.                                                           %
%                                                                         %
%	    As a supplement to the solver files, we also provide a  simple	  %
%	    setup of a sloshing tank problem discussed in references (1,2)	  %
%	    consisting of data required by the solver as input:  		      %
%	    - coord.mat : This file contains the coordinates  of  all  the	  %
%	      nodes in the domain of the problem.                             %
%	    - conn.mat : This file contains the connectivities of the elements%
%	      contained in the mesh of the problem.                           %
%	    - BCLeft.mat   ---     These files contain the nodal Ids of the	  %
%	    - BCRight.mat    | 	   nodes contained in the boundaries (Left,	  %
%	    - BCTop.mat      |---> Right, Top, Bottom) such that each entry	  %
%	    - BCBottom.mat   | 	   is an edge of the boundary connecting two  %
%	    - BCAll.mat    ---	   nodes. 'All' contains all the edges.       %
%                                                                         %
%	(2) The code is vectorized for efficient computations.			      %
%	(3) It takes coordinates, connectivity and nodal file data  of the	  %
%	    boundaries as input (An example of the sloshing  tank  problem	  %
%	    is  given  in  this  code).  These  data  files  can be easily 	  %
%	    generated using the renowned meshing tool Gmsh.			          %
%	(4) For trying a different multiphase problem,  some  parts of the 	  %
%	    code need to be changed according  to  the boundary conditions 	  %
%	    and domain of the problem.  This has been signaled in the code 	  %
%	    (where changes are to be made).                                   %
%   (5) Depending on the solver.outFreq which indicates the  frequency    %
%       at which the data will be dumped, *.plt  files  are  generated    %
%       which can be visualized in tecplot.                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;
tic
wrkDir = './' ;
problemString = 'st' ;
elemType = '3Tri' ;
problemType = '2D' ;
if (strcmp(problemType,'3D'))
    Flag3D = 1;
else
    Flag3D = 0;
end

%% Get coordinate, connectivity, nodal file data
crdstr = strcat(wrkDir,'coord.mat') ;
load(crdstr);
cnnstr = strcat(wrkDir,'conn.mat') ;
load(cnnstr);
BC1str = strcat(wrkDir,'BCLeft.mat') ;
BC2str = strcat(wrkDir,'BCRight.mat') ;
BC3str = strcat(wrkDir,'BCTop.mat') ;
BC4str = strcat(wrkDir,'BCBottom.mat') ;
BC5str = strcat(wrkDir,'BCAll.mat') ;
load(BC1str);
load(BC2str);
load(BC3str);
load(BC4str);
load(BC5str);
ndof = size(coord,1) ;

%% Define solver specific details

% Nonlinear iteration data:
solver.nLIterMin = 1 ;
solver.nLIterMax = 10 ;
solver.dt = 0.01 ;
solver.maxSteps = 2000 ;
solver.rhoinfty = 1.0 ;  % Generalized alpha damping parameter
solver.nLTol1 = 5e-4 ;   % Tolerance for error in Navier-Stokes and Allen-Cahn equations
solver.nLTol2 = 1e-3 ;   % Tolerance for error in residual error estimates
solver.epsilon = 0.01 ;  % Interfacial thickness parameter for the phase-field model

% Fluid properties:
fluid.dens1 = 1000 ;  
fluid.dens2 = 1 ;
fluid.visc1 = 1.0 ;
fluid.visc2 = 0.01 ;
fluid.gravFrc = [0, -1, 0];

% Initial boundary conditions:
fluid.vel0 = [0, 0, 0];
fluid.pres0 = 0.0 ;
fluid.phi0 = 0.0 ;

% Output data generation details:
solver.outFreq = 1 ;

% Boundary conditions:
% Dirichlet: 1
% Neumann: 2
% Symmetry/slip: 3
% No boundary (2D): 0
% variable type for velocity: 1 for ux, 2 for uy, 3 for uz

% Change the following appropriately for satisfying boundary conditions
bc.left.type = 3;
bc.right.type = 3 ;
bc.top.type = 3 ;
bc.bottom.type = 3 ;
bc.side_1.type = 0 ;
bc.side_2.type = 0 ;

bc.left.nNodes = size(unique(BCLeft),1);
bc.left.nodes = unique(BCLeft);
bc.left.nodes(bc.left.nodes==0) = [];
if (Flag3D == 1)
bc.left.value = [0.0 0.0 0.0] ;
bc.left.var = [1 2 3] ;
else
bc.left.value = [0.0 0.0] ;
bc.left.var = [1 2] ;
end

bc.right.nNodes = size(unique(BCRight),1);
bc.right.nodes = unique(BCRight);
bc.right.nodes(bc.right.nodes==0) = [];
if (Flag3D == 1)
bc.right.value = [0.0 0.0 0.0] ;
bc.right.var = [1 2 3] ;
else
bc.right.value = [0.0 0.0] ;
bc.right.var = [1 2] ;
end

bc.top.nNodes = size(unique(BCTop),1);
bc.top.nodes = unique(BCTop);
bc.top.nodes(bc.top.nodes==0) = [];
if (Flag3D == 1)
bc.top.value = [0.0 0.0 0.0] ;
bc.top.var = [1 2 3] ;
else
bc.top.value = [0.0 0.0] ;
bc.top.var = [1 2] ;
end

bc.bottom.nNodes = size(unique(BCBottom),1);
bc.bottom.nodes = unique(BCBottom);
bc.bottom.nodes(bc.bottom.nodes==0) = [];
if (Flag3D == 1)
bc.bottom.value = [0.0 0.0 0.0] ;
bc.bottom.var = [1 2 3] ;
else
bc.bottom.value = [0.0 0.0] ;
bc.bottom.var = [1 2] ;
end

dirichlet = [] ;
neumann = [BCLeft; BCRight; BCTop; BCBottom]; % For order parameter phi (for solving Allen-Cahn equation)
neumann = unique(sort(neumann,2),'rows');

%% Define the quadrature data for the solver
if strcmp(elemType,'3Tri')
    gP = ...
   [1/6,  1/6
    2/3,  1/6
    1/6,  2/3] ;
    gW = ...
   [1/3,  1/3,  1/3] ;
 
    N(:,1) = (1.-gP(:,1)-gP(:,2)) ;
    N(:,2) = (gP(:,1)) ;
    N(:,3) = (gP(:,2)) ;
    
    Nx(:,1) = -ones(3,1) ;
    Nx(:,2) =  ones(3,1) ;
    Nx(:,3) =  zeros(3,1) ; 
    Ny(:,1) = -ones(3,1) ;
    Ny(:,2) =  zeros(3,1) ;
    Ny(:,3) =  ones(3,1) ;    
elseif strcmp(elemType,'4Quad')
    gP = ...
   [-5.7735026918962584E-01, -5.7735026918962584E-01
     5.7735026918962584E-01, -5.7735026918962584E-01
    -5.7735026918962584E-01,  5.7735026918962584E-01
     5.7735026918962584E-01,  5.7735026918962584E-01] ;
    gW = [1, 1, 1, 1 ] ;
    
    N(:,1) = 0.25.*(1-gP(:,1)).*(1-gP(:,2)) ;
    N(:,2) = 0.25.*(1+gP(:,1)).*(1-gP(:,2)) ;
    N(:,3) = 0.25.*(1+gP(:,1)).*(1+gP(:,2)) ;
    N(:,4) = 0.25.*(1-gP(:,1)).*(1+gP(:,2)) ;
    
    Nx(:,1) = -0.25.*(1-gP(:,2)) ;
    Nx(:,2) =  0.25.*(1-gP(:,2)) ;
    Nx(:,3) =  0.25.*(1+gP(:,2)) ;
    Nx(:,4) = -0.25.*(1+gP(:,2)) ;
    Ny(:,1) = -0.25.*(1-gP(:,1)) ;
    Ny(:,2) = -0.25.*(1+gP(:,1)) ;
    Ny(:,3) =  0.25.*(1+gP(:,1)) ;
    Ny(:,4) =  0.25.*(1-gP(:,1)) ;
elseif strcmp(elemType,'6Prism')
    gP = ...
   [3.3333333333333331E-01,  3.3333333333333331E-01, -5.7735026918962584E-01
    1.3333333333333333E+00,  3.3333333333333331E-01, -5.7735026918962584E-01
    3.3333333333333331E-01,  1.3333333333333333E+00, -5.7735026918962584E-01
    3.3333333333333331E-01,  3.3333333333333331E-01,  5.7735026918962584E-01
    1.3333333333333333E+00,  3.3333333333333331E-01,  5.7735026918962584E-01
    3.3333333333333331E-01,  1.3333333333333333E+00,  5.7735026918962584E-01] ;
    gW = ...
   [6.6666666666666667E-01,  6.6666666666666667E-01,  6.6666666666666667E-01,...
    6.6666666666666667E-01,  6.6666666666666667E-01,  6.6666666666666667E-01] ;

    N(:,1) = 0.25*(2-gP(:,1)-gP(:,2)).*(1-gP(:,3)) ;
    N(:,2) = 0.25*(gP(:,1)).*(1-gP(:,3)) ;
    N(:,3) = 0.25*(gP(:,2)).*(1-gP(:,3)) ;
    N(:,4) = 0.25*(2-gP(:,1)-gP(:,2)).*(1+gP(:,3)) ;
    N(:,5) = 0.25*(gP(:,1)).*(1+gP(:,3)) ;
    N(:,6) = 0.25*(gP(:,2)).*(1+gP(:,3)) ;
    
    Nx(:,1) = -0.25.*(1-gP(:,3)) ;
    Nx(:,2) =  0.25.*(1-gP(:,3)) ;
    Nx(:,3) =  zeros(6,1) ;
    Nx(:,4) = -0.25.*(1+gP(:,3)) ;
    Nx(:,5) =  0.25.*(1+gP(:,3)) ;
    Nx(:,6) =  zeros(6,1) ;    
    Ny(:,1) = -0.25.*(1-gP(:,3)) ;
    Ny(:,2) =  zeros(6,1) ;
    Ny(:,3) =  0.25.*(1-gP(:,3)) ;
    Ny(:,4) = -0.25.*(1+gP(:,3)) ;
    Ny(:,5) =  zeros(6,1) ;
    Ny(:,6) =  0.25.*(1+gP(:,3)) ;    
    Nz(:,1) = -0.25.*(2-gP(:,1)-gP(:,2)) ;
    Nz(:,2) = -0.25.*(gP(:,1)) ;
    Nz(:,3) = -0.25.*(gP(:,2)) ;
    Nz(:,4) =  0.25.*(2-gP(:,1)-gP(:,2)) ;
    Nz(:,5) =  0.25.*(gP(:,1)) ;
    Nz(:,6) =  0.25.*(gP(:,2)) ;
elseif strcmp(elemType,'8Hex')
    gP = ...
   [-5.7735026918962584E-01, -5.7735026918962584E-01, -5.7735026918962584E-01
     5.7735026918962584E-01, -5.7735026918962584E-01, -5.7735026918962584E-01
    -5.7735026918962584E-01,  5.7735026918962584E-01, -5.7735026918962584E-01
     5.7735026918962584E-01,  5.7735026918962584E-01, -5.7735026918962584E-01
    -5.7735026918962584E-01, -5.7735026918962584E-01,  5.7735026918962584E-01
     5.7735026918962584E-01, -5.7735026918962584E-01,  5.7735026918962584E-01
    -5.7735026918962584E-01,  5.7735026918962584E-01,  5.7735026918962584E-01
     5.7735026918962584E-01,  5.7735026918962584E-01,  5.7735026918962584E-01] ;
    gW = ...
   [1, 1, 1, 1, 1, 1, 1, 1 ] ;

    N(:,1) = 0.125.*(1-gP(:,1)).*(1-gP(:,2)).*(1-gP(:,3)) ;
    N(:,2) = 0.125.*(1+gP(:,1)).*(1-gP(:,2)).*(1-gP(:,3)) ;
    N(:,3) = 0.125.*(1+gP(:,1)).*(1+gP(:,2)).*(1-gP(:,3)) ;
    N(:,4) = 0.125.*(1-gP(:,1)).*(1+gP(:,2)).*(1-gP(:,3)) ;
    N(:,5) = 0.125.*(1-gP(:,1)).*(1-gP(:,2)).*(1+gP(:,3)) ;
    N(:,6) = 0.125.*(1+gP(:,1)).*(1-gP(:,2)).*(1+gP(:,3)) ;
    N(:,7) = 0.125.*(1+gP(:,1)).*(1+gP(:,2)).*(1+gP(:,3)) ;
    N(:,8) = 0.125.*(1-gP(:,1)).*(1+gP(:,2)).*(1+gP(:,3)) ;

    Nx(:,1) = -0.125.*(1-gP(:,2)).*(1-gP(:,3)) ;
    Nx(:,2) =  0.125.*(1-gP(:,2)).*(1-gP(:,3)) ;
    Nx(:,3) =  0.125.*(1+gP(:,2)).*(1-gP(:,3)) ;
    Nx(:,4) = -0.125.*(1+gP(:,2)).*(1-gP(:,3)) ;
    Nx(:,5) = -0.125.*(1-gP(:,2)).*(1+gP(:,3)) ;
    Nx(:,6) =  0.125.*(1-gP(:,2)).*(1+gP(:,3)) ;
    Nx(:,7) =  0.125.*(1+gP(:,2)).*(1+gP(:,3)) ;
    Nx(:,8) = -0.125.*(1+gP(:,2)).*(1+gP(:,3)) ; 
    Ny(:,1) = -0.125.*(1-gP(:,1)).*(1-gP(:,3)) ;
    Ny(:,2) = -0.125.*(1+gP(:,1)).*(1-gP(:,3)) ;
    Ny(:,3) =  0.125.*(1+gP(:,1)).*(1-gP(:,3)) ;
    Ny(:,4) =  0.125.*(1-gP(:,1)).*(1-gP(:,3)) ;
    Ny(:,5) = -0.125.*(1-gP(:,1)).*(1+gP(:,3)) ;
    Ny(:,6) = -0.125.*(1+gP(:,1)).*(1+gP(:,3)) ;
    Ny(:,7) =  0.125.*(1+gP(:,1)).*(1+gP(:,3)) ;
    Ny(:,8) =  0.125.*(1-gP(:,1)).*(1+gP(:,3)) ; 
    Nz(:,1) = -0.125.*(1-gP(:,1)).*(1-gP(:,2)) ;
    Nz(:,2) = -0.125.*(1+gP(:,1)).*(1-gP(:,2)) ;
    Nz(:,3) = -0.125.*(1+gP(:,1)).*(1+gP(:,2)) ;
    Nz(:,4) = -0.125.*(1-gP(:,1)).*(1+gP(:,2)) ;
    Nz(:,5) =  0.125.*(1-gP(:,1)).*(1-gP(:,2)) ;
    Nz(:,6) =  0.125.*(1+gP(:,1)).*(1-gP(:,2)) ;
    Nz(:,7) =  0.125.*(1+gP(:,1)).*(1+gP(:,2)) ;
    Nz(:,8) =  0.125.*(1-gP(:,1)).*(1+gP(:,2)) ;
else
    warning('elemType not found!!');
end

Nx = Nx' ;
Ny = Ny' ;
if (Flag3D == 1)
Nz = Nz' ;
end
nQuad = length(gW) ;

%% Derived quantities

% Gen-alpha parameters
pmc.alphaM = 0.5*(3-solver.rhoinfty)/(1+solver.rhoinfty) ;
pmc.alpha = 1/(1+solver.rhoinfty) ;
pmc.gamma = 0.5 + pmc.alphaM - pmc.alpha ;

nodeId = coord(:,1) ;
crd = coord(:,2:4) ;
cnn = conn(:,1:end) ;
nen = size(cnn,2) ;
nElem = size(cnn,1) ;
ndof = size(crd,1);

if (Flag3D == 1)
Sol.u = zeros(ndof,3,1);
Sol.uDot = zeros(ndof,3,1);
Sol.u(:,1,:) = fluid.vel0(1) ;
Sol.u(:,2,:) = fluid.vel0(2) ;
Sol.u(:,3,:) = fluid.vel0(3) ;
Sol.uAlpha = zeros(ndof,3,1) ;
Sol.uDotAlpha = zeros(ndof,3,1) ;
Sol.uPrev = Sol.u ;
Sol.uDotPrev = Sol.uDot ;
else
Sol.u = zeros(ndof,2,1);
Sol.uDot = zeros(ndof,2,1);
Sol.u(:,1,:) = fluid.vel0(1) ;
Sol.u(:,2,:) = fluid.vel0(2) ;
Sol.uAlpha = zeros(ndof,2,1) ;
Sol.uDotAlpha = zeros(ndof,2,1) ;
Sol.uPrev = Sol.u ;
Sol.uDotPrev = Sol.uDot ;
end

Sol.p = fluid.pres0*ones(ndof,1);
% Initial profile for phi (Change this for different initial profile for phi)
Sol.phi = -tanh( (crd(:,2) - (1.01 + 0.1*sin((crd(:,1) - 0.5)*pi)))/(sqrt(2)*solver.epsilon));
Sol.phiDot = zeros(ndof,1) ;

type = uint8(ones(size(crd,1),1));
crd(:,3) = [];
Sol.node = crd ;
Sol.elem = cnn ;
Sol.type = type ;
Sol.Dirichlet = dirichlet ;
Sol.Neumann = neumann ;

Sol.phiAlpha = zeros(ndof,1) ;
Sol.phiDotAlpha = zeros(ndof,1) ;
Sol.pPrev = Sol.p ;
Sol.phiPrev = Sol.phi ;
Sol.phiDotPrev = Sol.phiDot ;

% Refinement of the background mesh to get initial refinement
for i=1:8
    for i=1:nen
        xxf(:,i) = crd(cnn(:,i),1);
        yyf(:,i) = crd(cnn(:,i),2);
        if (Flag3D == 1)
        zzf(:,i) = crd(cnn(:,i),3);
        end
        src(:,i) = Sol.phi(cnn(:,i),1) ;
    end
    for p = 1:nQuad  
      if (Flag3D ~= 1)
            J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                 yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
            if size(J,2)==1
                J = J';
            end
            volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );   
      else
            J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                 yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                 zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
            if size(J,2)==1
                 J = J';
            end
            volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                      J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                      J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
      end
      volume = abs(volume);
      if (Flag3D ~= 1)
                DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
                DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
      else
                DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
                      (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
                      (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
                DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
                      (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
                      (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
                DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
                      (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
                      (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);
      end

      localGradPhiX  = sum(repmat(DNDx(p,:),nElem,1).*src,2);
      localGradPhiY  = sum(repmat(DNDy(p,:),nElem,1).*src,2);
    
      volume = abs(volume);

      gradphi =  sqrt( localGradPhiX.^2 + localGradPhiY.^2 ) ;
      etaRPhi(:,p) = abs(gradphi) ;
    end
    etaRPhi = sum(etaRPhi,2) ;

    theta = 0.61 ;
    theta_c = 0.3 ;
    
    [Sol,etaRPhi] = bisection(Sol,etaRPhi.^2,theta);
    crd = [Sol.node zeros(size(Sol.node,1),1)] ;
    cnn = Sol.elem ;
    ndof = size(crd,1) ;
    nElem = size(cnn,1) ;
    clear xxf yyf zzf src localSrc etaRPhi
    % Initial profile for phi (Change this for different initial profile for phi)
    Sol.phi = -tanh( (crd(:,2) - (1.01 + 0.1*sin((crd(:,1) - 0.5)*pi)))/(sqrt(2)*solver.epsilon));
end
% Initial profile for phi (Change this for different initial profile for phi)
Sol.phi = -tanh( (crd(:,2) - (1.01 + 0.1*sin((crd(:,1) - 0.5)*pi)))/(sqrt(2)*solver.epsilon));
Sol.phiDot = zeros(ndof,1) ;
Sol.phiPrev = Sol.phi ;
Sol.phiDotPrev = Sol.phiDot ;

%% Time loop starts here
for timeStep = 1:solver.maxSteps

    fprintf('Time step:%d\n',timeStep);

    % Predict the solution
    Sol.u = Sol.uPrev ;
    Sol.uDot = (pmc.gamma - 1)/pmc.gamma * Sol.uDotPrev ;
    Sol.p = Sol.pPrev ;
    Sol.phi = Sol.phiPrev ;
    Sol.phiDot = (pmc.gamma - 1)/pmc.gamma * Sol.phiDotPrev ;
    
    % Nonlinear iterations start here
    for nLiter = 1:solver.nLIterMax
        
        iif = zeros(nen^2*nElem,1); jjf = zeros(nen^2*nElem,1);
        index = 0;
        for i = 1:nen
            for j = 1:nen
                iif(index+1:index+nElem) = double(cnn(:,i)); 
                jjf(index+1:index+nElem) = double(cnn(:,j));  
                index = index + nElem;
            end
        end
        
	% Define the boundary nodes (Change this according to the problem)
        zerotypenodes = find(Sol.type == 0);
        bc.left.nodes = find(Sol.node(:,1)< 0.0+1e-8) ;
        bc.left.nodes = setdiff(bc.left.nodes, zerotypenodes)' ;
        bc.right.nodes = find(Sol.node(:,1)> 1.0-1e-8) ;
        bc.right.nodes = setdiff(bc.right.nodes, zerotypenodes)' ;
        bc.top.nodes = find(Sol.node(:,2)> 1.5-1e-8) ;
        bc.top.nodes = setdiff(bc.top.nodes, zerotypenodes)' ;
        bc.bottom.nodes = find(Sol.node(:,2)< 0.0+1e-8) ;
        bc.bottom.nodes = setdiff(bc.bottom.nodes, zerotypenodes)' ;

        % Satisfy Dirichlet boundary condition for velocity and pressure (Change this according to the problem)
        if (bc.left.type == 1)
            dirichletNodes = bc.left.nodes' ;
            if (Flag3D == 1)
            Sol.u(dirichletNodes,bc.left.var(1),1) = bc.left.value(1).*ones(size(bc.left.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.left.var(2),1) = bc.left.value(2).*ones(size(bc.left.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.left.var(3),1) = bc.left.value(3).*ones(size(bc.left.nodes',1),1) ;
            else
            Sol.u(dirichletNodes,bc.left.var(1),1) = bc.left.value(1).*ones(size(bc.left.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.left.var(2),1) = bc.left.value(2).*ones(size(bc.left.nodes',1),1) ;
            end
        elseif (bc.left.type == 3)
            dirichletNodes = bc.left.nodes' ;
            Sol.u(dirichletNodes,1,1) = zeros(size(bc.left.nodes',1),1) ;
        end
        if (bc.right.type == 1)
            dirichletNodes = bc.right.nodes' ;
            if (Flag3D == 1)
            Sol.u(dirichletNodes,bc.right.var(1),1) = bc.right.value(1).*ones(size(bc.right.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.right.var(2),1) = bc.right.value(2).*ones(size(bc.right.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.right.var(3),1) = bc.right.value(3).*ones(size(bc.right.nodes',1),1) ;
            else
            Sol.u(dirichletNodes,bc.right.var(1),1) = bc.right.value(1).*ones(size(bc.right.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.right.var(2),1) = bc.right.value(2).*ones(size(bc.right.nodes',1),1) ;
            end
        elseif (bc.right.type == 3)
            dirichletNodes = bc.right.nodes' ;
            Sol.u(dirichletNodes,1,1) = zeros(size(bc.right.nodes',1),1) ;
        end
        if (bc.top.type == 1)
            dirichletNodes = bc.top.nodes' ;
            if (Flag3D == 1)
            Sol.u(dirichletNodes,bc.top.var(1),1) = bc.top.value(1).*ones(size(bc.top.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.top.var(2),1) = bc.top.value(2).*ones(size(bc.top.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.top.var(3),1) = bc.top.value(3).*ones(size(bc.top.nodes',1),1) ;
            else
            Sol.u(dirichletNodes,bc.top.var(1),1) = bc.top.value(1).*ones(size(bc.top.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.top.var(2),1) = bc.top.value(2).*ones(size(bc.top.nodes',1),1) ;
            end
        elseif (bc.top.type == 3)
            dirichletNodes = bc.top.nodes' ;
            Sol.u(dirichletNodes,2,1) = zeros(size(bc.top.nodes',1),1) ;
        end
        if (bc.bottom.type == 1)
            dirichletNodes = bc.bottom.nodes' ;
            if (Flag3D == 1)
            Sol.u(dirichletNodes,bc.bottom.var(1),1) = bc.bottom.value(1).*ones(size(bc.bottom.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.bottom.var(2),1) = bc.bottom.value(2).*ones(size(bc.bottom.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.bottom.var(3),1) = bc.bottom.value(3).*ones(size(bc.bottom.nodes',1),1) ;
            else
            Sol.u(dirichletNodes,bc.bottom.var(1),1) = bc.bottom.value(1).*ones(size(bc.bottom.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.bottom.var(2),1) = bc.bottom.value(2).*ones(size(bc.bottom.nodes',1),1) ;
            end
        elseif (bc.bottom.type == 3)
            dirichletNodes = bc.bottom.nodes' ;
            Sol.u(dirichletNodes,2,1) = zeros(size(bc.bottom.nodes',1),1) ;
        end
        if (bc.side_1.type == 1)
            dirichletNodes = bc.side_1.nodes' ;
            Sol.u(dirichletNodes,bc.side_1.var,1) = bc.side_1.value.*ones(size(bc.side_1.nodes,1),1) ;
        elseif (bc.side_1.type == 3)
            dirichletNodes = bc.side_1.nodes' ;
            Sol.u(dirichletNodes,3,1) = zeros(size(bc.side_1.nodes,1),1) ;
        end
        if (bc.side_2.type == 1)
            dirichletNodes = bc.side_2.nodes' ;
            Sol.u(dirichletNodes,bc.side_2.var,1) = bc.side_2.value.*ones(size(bc.side_2.nodes,1),1) ;
        elseif (bc.side_2.type == 3)
            dirichletNodes = bc.side_2.nodes' ;
            Sol.u(dirichletNodes,3,1) = zeros(size(bc.side_2.nodes,1),1) ;
        end
        
        % 2D dirichlet condition in 3D mesh 
        if (Flag3D == 1)
        Sol.u(:,3,1) = zeros(size(crd,1),1) ;
        end
        Sol.p(bc.top.nodes',1) = 0.0 ;
        
        
        % Interpolate for alpha values for Gen-alpha
        Sol.uAlpha = Sol.uPrev + pmc.alpha.*(Sol.u - Sol.uPrev) ;
        Sol.uDotAlpha = Sol.uDotPrev + pmc.alphaM.*(Sol.uDot - Sol.uDotPrev) ;
        Sol.phiAlpha = Sol.phiPrev + pmc.alpha.*(Sol.phi - Sol.phiPrev) ;
        Sol.phiDotAlpha = Sol.phiDotPrev + pmc.alphaM.*(Sol.phiDot - Sol.phiDotPrev) ;
        
        %% Navier-Stokes equations
        xxf = zeros(size(cnn));
        yyf = zeros(size(cnn));
        zzf = zeros(size(cnn));
        ux = zeros(size(cnn));
        uy = zeros(size(cnn));
        if (Flag3D == 1)
        uz = zeros(size(cnn));
        end

        for i=1:nen
            xxf(:,i) = crd(cnn(:,i),1);
            yyf(:,i) = crd(cnn(:,i),2);
            if (Flag3D == 1)
            zzf(:,i) = crd(cnn(:,i),3);
            end
            ux(:,i) =  Sol.uAlpha(cnn(:,i),1,1) ;
            uxDot(:,i) = Sol.uDotAlpha(cnn(:,i),1,1) ;
            uy(:,i) =  Sol.uAlpha(cnn(:,i),2,1) ;
            uyDot(:,i) = Sol.uDotAlpha(cnn(:,i),2,1) ;
            if (Flag3D == 1)
            uz(:,i) =  Sol.uAlpha(cnn(:,i),3,1) ;
            uzDot(:,i) = Sol.uDotAlpha(cnn(:,i),3,1) ;
            end
            phi(:,i) = Sol.phiAlpha(cnn(:,i),1) ;
            pres(:,i) = Sol.p(cnn(:,i),1) ;
        end
        
        % Galerkin terms for Navier-Stokes equations
        sA = [] ;
        sA1 = zeros(nen^2*nElem,nQuad); 
        sA2 = zeros(nen^2*nElem,nQuad);
        sA3 = zeros(nen^2*nElem,nQuad);

        sA4 = zeros(nen^2*nElem,nQuad);
        sA5 = zeros(nen^2*nElem,nQuad);
        sA6 = zeros(nen^2*nElem,nQuad);
        
        sA7 = zeros(nen^2*nElem,nQuad); 
        sA8 = zeros(nen^2*nElem,nQuad);
        sA9 = zeros(nen^2*nElem,nQuad);

        sA10 = zeros(nen^2*nElem,nQuad);
        sA11 = zeros(nen^2*nElem,nQuad);
        sA12 = zeros(nen^2*nElem,nQuad);
        
        sA13 = zeros(nen^2*nElem,nQuad);
        
        sA14 = zeros(nen^2*nElem,nQuad);
        sA15 = zeros(nen^2*nElem,nQuad);
        sA16 = zeros(nen^2*nElem,nQuad);
        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
            end
            
            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            
            ind1 = find(localPhi > 1.0) ;
            ind2 = find(localPhi < -1.0) ;
            ind3 = find(localPhi <= 1.0 & localPhi >= -1.0) ;
            dens = zeros(size(localPhi,1),1) ;
            visc = zeros(size(localPhi,1),1) ;
            
            dens(ind1) = fluid.dens1 ;
            visc(ind1) = fluid.visc1 ;
            dens(ind2) = fluid.dens2 ;
            visc(ind2) = fluid.visc2 ;
            dens(ind3) = fluid.dens1*0.5*(1+localPhi(ind3)) + fluid.dens2*0.5*(1-localPhi(ind3)) ;
            visc(ind3) = fluid.visc1*0.5*(1+localPhi(ind3)) + fluid.visc2*0.5*(1-localPhi(ind3)) ;
            clear ind1 ind2 ind3
            
            if (Flag3D ~= 1)
                DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
                DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);   
            else
                DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
                      (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
                      (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
                DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
                      (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
                      (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
                DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
                      (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
                      (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);
                  
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);
                locUZ  = sum(repmat(N(p,:),nElem,1).*uz,2);
            end
            
            index = 0;
            for i = 1:nen
                for j = 1:nen
                    % Galerkin inertia term
                    Mij = gW(p)*(N(p,i)*N(p,j));
                    Mij = Mij.*dens ;
                    Mij = Mij.*volume;
                    sA(index+1:index+nElem,p) = Mij;
                    
                    % Galerkin convection term
                    Aij_1 = gW(p)*N(p,i)*(locUX.*DNDx(:,j));
                    Aij_2 = gW(p)*N(p,i)*(locUY.*DNDy(:,j));
                    Aij_1 = Aij_1.*volume.*dens;
                    Aij_2 = Aij_2.*volume.*dens;
                    sA1(index+1:index+nElem,p) = Aij_1;
                    sA2(index+1:index+nElem,p) = Aij_2;
                    
                    % Galerkin viscous diffusion term
                    Aij_4 = gW(p)*(DNDx(:,i).*DNDx(:,j));
                    Aij_5 = gW(p)*(DNDy(:,i).*DNDy(:,j));
                    Aij_7 = gW(p)*(DNDx(:,i).*DNDy(:,j));
                    Aij_4 = Aij_4.*volume.*visc;
                    Aij_5 = Aij_5.*volume.*visc;
                    Aij_7 = Aij_7.*volume.*visc;
                    sA4(index+1:index+nElem,p) = Aij_4;
                    sA5(index+1:index+nElem,p) = Aij_5;
                    sA7(index+1:index+nElem,p) = Aij_7;
                    
                    % Galerkin pressure term (gradP)
                    Aij_10 = gW(p)*(DNDx(:,i).*N(p,j));
                    Aij_11 = gW(p)*(DNDy(:,i).*N(p,j));
                    Aij_10 = Aij_10.*volume;
                    Aij_11 = Aij_11.*volume; 
                    sA10(index+1:index+nElem,p) = Aij_10;
                    sA11(index+1:index+nElem,p) = Aij_11;
                    
                    % Galerkin source term (gravity force)
                    Aij_13 = gW(p)*(N(p,i)*N(p,j));
                    Aij_13 = Aij_13.*dens.*volume ;
                    sA13(index+1:index+nElem,p) = Aij_13;
                    
                    % Galerkin continuity term
                    Aij_14 = gW(p)*(N(p,i)).*(DNDx(:,j)) ;
                    Aij_15 = gW(p)*(N(p,i)).*(DNDy(:,j)) ;
                    Aij_14 = Aij_14.*volume ;
                    Aij_15 = Aij_15.*volume ;
                    sA14(index+1:index+nElem,p) = Aij_14;
                    sA15(index+1:index+nElem,p) = Aij_15;
                         
                    if (Flag3D == 1)
                    Aij_3 = gW(p)*(locUZ.*DNDz(:,j)*N(p,i));
                    Aij_3 = Aij_3.*volume.*dens;
                    sA3(index+1:index+nElem,p) = Aij_3;
                    Aij_6 = gW(p)*(DNDz(:,i).*DNDz(:,j));
                    Aij_8 = gW(p)*(DNDy(:,i).*DNDz(:,j));
                    Aij_9 = gW(p)*(DNDx(:,i).*DNDz(:,j));
                    Aij_6 = Aij_6.*volume.*visc;
                    Aij_8 = Aij_8.*volume.*visc;
                    Aij_9 = Aij_9.*volume.*visc;
                    sA6(index+1:index+nElem,p) = Aij_6;
                    sA8(index+1:index+nElem,p) = Aij_8;
                    sA9(index+1:index+nElem,p) = Aij_9;
                    Aij_12 = gW(p)*(DNDz(:,i).*N(p,j));
                    Aij_12 = Aij_12.*volume;
                    sA12(index+1:index+nElem,p) = Aij_12;
                    Aij_16 = gW(p)*(N(p,i)).*(DNDz(:,j)) ;
                    Aij_16 = Aij_16.*volume ;
                    sA16(index+1:index+nElem,p) = Aij_16;
                    end

                    index = index + nElem;
                end
            end
        end
        sA = sum(sA,2);
        sA1 = sum(sA1,2);
        sA2 = sum(sA2,2);
        sA4 = sum(sA4,2);
        sA5 = sum(sA5,2);
        sA7 = sum(sA7,2);
        sA10 = sum(sA10,2);
        sA11 = sum(sA11,2);
        sA13 = sum(sA13,2);
        sA14 = sum(sA14,2);
        sA15 = sum(sA15,2);


        if (Flag3D == 1)
        sA3 = sum(sA3,2);
        sA6 = sum(sA6,2);
        sA8 = sum(sA8,2);
        sA9 = sum(sA9,2);
        sA12 = sum(sA12,2);
        sA16 = sum(sA16,2);
        end
        
        % Assemble the matrix
        Mf = sparse(iif,jjf,sA,ndof,ndof);  
        ZeroF = sparse(ndof,ndof);
        
        A1 = sparse(iif,jjf,sA1,ndof,ndof);
        A2 = sparse(iif,jjf,sA2,ndof,ndof);
        A4 = sparse(iif,jjf,sA4,ndof,ndof);
        A5 = sparse(iif,jjf,sA5,ndof,ndof);
        A7 = sparse(iif,jjf,sA7,ndof,ndof);
        A10 = sparse(iif,jjf,sA10,ndof,ndof);
        A11 = sparse(iif,jjf,sA11,ndof,ndof);
        A13 = sparse(iif,jjf,sA13,ndof,ndof);
        A14 = sparse(iif,jjf,sA14,ndof,ndof);
        A15 = sparse(iif,jjf,sA15,ndof,ndof);
        
        if (Flag3D == 1)
        A3 = sparse(iif,jjf,sA3,ndof,ndof);
        A6 = sparse(iif,jjf,sA6,ndof,ndof);
        A8 = sparse(iif,jjf,sA8,ndof,ndof);
        A9 = sparse(iif,jjf,sA9,ndof,ndof);
        A12 = sparse(iif,jjf,sA12,ndof,ndof);
        A16 = sparse(iif,jjf,sA16,ndof,ndof);
        end
        
        if (Flag3D == 1)
        Mf = [Mf ZeroF ZeroF ZeroF;...
              ZeroF Mf ZeroF ZeroF;...
              ZeroF ZeroF Mf ZeroF;...
              ZeroF ZeroF ZeroF ZeroF];
        Conv = [A1+A2+A3 ZeroF ZeroF ZeroF;...
                ZeroF A1+A2+A3 ZeroF ZeroF;...
                ZeroF ZeroF A1+A2+A3 ZeroF;...
                ZeroF ZeroF ZeroF ZeroF];
        Kf = [2*A4+A5+A6 A7' A9' ZeroF;...
              A7 A4+2*A5+A6 A8' ZeroF;...
              A9 A8 A4+A5+2*A6 ZeroF;...
              ZeroF ZeroF ZeroF ZeroF];
        Fp = [ZeroF ZeroF ZeroF A10;...
              ZeroF ZeroF ZeroF A11;...
              ZeroF ZeroF ZeroF A12;...
              ZeroF ZeroF ZeroF ZeroF];
        GMassCons = [ZeroF ZeroF ZeroF ZeroF;...
              ZeroF ZeroF ZeroF ZeroF;...
              ZeroF ZeroF ZeroF ZeroF;...
              A14 A15 A16 ZeroF];
        Src = [A13.*fluid.gravFrc(1) ZeroF ZeroF ZeroF;...  
               ZeroF A13.*fluid.gravFrc(2) ZeroF ZeroF;...
               ZeroF ZeroF A13.*fluid.gravFrc(3) ZeroF;...
               ZeroF ZeroF ZeroF ZeroF];
        else
        Mf = [Mf ZeroF ZeroF;...
              ZeroF Mf ZeroF;...
              ZeroF ZeroF ZeroF];
        Conv = [A1+A2 ZeroF ZeroF;...
                ZeroF A1+A2 ZeroF;...
                ZeroF ZeroF ZeroF];
        Kf = [2*A4+A5 A7' ZeroF;...
              A7 A4+2*A5  ZeroF;...
              ZeroF ZeroF ZeroF];
        Fp = [ZeroF ZeroF A10;...
              ZeroF ZeroF A11;...
              ZeroF ZeroF ZeroF];
        GMassCons = [ZeroF ZeroF ZeroF;...
              ZeroF ZeroF ZeroF;...
              A14 A15 ZeroF];
        Src = [A13.*fluid.gravFrc(1) ZeroF ZeroF;...  
               ZeroF A13.*fluid.gravFrc(2) ZeroF;...;...
               ZeroF ZeroF ZeroF];
        end
        Mflow = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Mf;

        MatlabFlow = Mflow + Conv  + Kf - Fp + GMassCons;
        
        clear sA Mij sA1 sA2 sA3 Aij_1 Aij_2 Aij_3 sA4 sA5 sA6 sA7 sA8 sA9
        clear sA10 sA11 sA12 sA13 sA14 sA15 sA16
        clear Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12 Aij_13 Aij_14 Aij_15 Aij_16
        clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16
        
        uAlpha = Sol.uAlpha(:,:,1) ;
        uAlpha = [uAlpha(:); zeros(ndof,1)];
        uDotAlpha = Sol.uDotAlpha(:,:,1) ;
        uDotAlpha = [uDotAlpha(:); zeros(ndof,1)];
        if (size(Sol.p,1) == 1 & size(Sol.p,2)>1)
            Sol.p = [Sol.p]';
        end
    
        if (Flag3D == 1)
        p1 = [zeros(ndof*3,1);Sol.p];
        gravVec = [ones(3*ndof,1); zeros(ndof,1)];
        else
        p1 = [zeros(ndof*2,1);Sol.p];
        gravVec = [ones(2*ndof,1); zeros(ndof,1)];
        end
   
        RHS = -(Mf * uDotAlpha(:)) - (Conv + Kf + GMassCons)* uAlpha(:) ;
        RHS = RHS + (Fp * p1) + (Src)*gravVec ;
        
        % Petrov-Galerkin stabilization terms for Navier-Stokes equations
        % Part 1: (tauM/rho)(rho*u.gradN).(rho*dU/dt) + 
        %         (tauM/rho)(rho*u.gradN).(rho*u.gradU) +
        %         (tauM/rho)(gradN).(rho*dU/dt)
        %
        sA1 = zeros(nen^2*nElem,nQuad); 
        sA2 = zeros(nen^2*nElem,nQuad);
        sA3 = zeros(nen^2*nElem,nQuad);

        sA4 = zeros(nen^2*nElem,nQuad);
        sA5 = zeros(nen^2*nElem,nQuad);
        sA6 = zeros(nen^2*nElem,nQuad);      
        sA7 = zeros(nen^2*nElem,nQuad); 
        sA8 = zeros(nen^2*nElem,nQuad);
        sA9 = zeros(nen^2*nElem,nQuad);
        sA10 = zeros(nen^2*nElem,nQuad);
        sA11 = zeros(nen^2*nElem,nQuad);
        sA12 = zeros(nen^2*nElem,nQuad);
        
        sA13 = zeros(nen^2*nElem,nQuad); 
        sA14 = zeros(nen^2*nElem,nQuad);
        sA15 = zeros(nen^2*nElem,nQuad);
        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
                xsInv(:,4) = J(:,1).*J(:,4) - J(:,2).*J(:,3) ;
                xsInv(:,4) = 1./xsInv(:,4) ;
                xsInv(:,1) = xsInv(:,4).* J(:,4) ;
                xsInv(:,2) = xsInv(:,4).* (-J(:,2)) ;
                xsInv(:,3) = xsInv(:,4).* (-J(:,3)) ;
                xsInv(:,4) = xsInv(:,4).* J(:,1) ;
                
                if (strcmp(elemType,'3Tri'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                elseif (strcmp(elemType,'4Quad'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                end
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                      
                xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
                xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
                xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
                xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
                xsInv(:,9) = 1./xsInv(:,9) ;
                xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
                xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
                xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
                xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
                xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
                xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
                xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
                xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
                xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
                if (strcmp(elemType,'4Tet'))
                    c1 = 1.259921049894873E0 ;
                    c2 = 6.299605249474365D-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * (xsInv(:,4) + xsInv(:,7)) ;
                    a2 = c1 * xsInv(:,4) + c2 * (xsInv(:,1) + xsInv(:,7)) ;
                    a3 = c1 * xsInv(:,7) + c2 * (xsInv(:,1) + xsInv(:,4)) ;
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,2) + c2 * (xsInv(:,5) + xsInv(:,8)) ;
                    a2 = c1 * xsInv(:,5) + c2 * (xsInv(:,2) + xsInv(:,8)) ;
                    a3 = c1 * xsInv(:,8) + c2 * (xsInv(:,2) + xsInv(:,5)) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,3) + c2 * (xsInv(:,6) + xsInv(:,9)) ;
                    a2 = c1 * xsInv(:,6) + c2 * (xsInv(:,3) + xsInv(:,9)) ;
                    a3 = c1 * xsInv(:,9) + c2 * (xsInv(:,3) + xsInv(:,6)) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* a3;
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
                elseif (strcmp(elemType,'6Prism'))
                    c1 = 1.154700538379252E0 ;
                    c2 = 5.773502691896259E-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
                    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
                    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
                    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
                    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
                    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
                elseif (strcmp(elemType,'8Hex'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,4).*xsInv(:,4) + xsInv(:,7).*xsInv(:,7);
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,5).*xsInv(:,5) + xsInv(:,8).*xsInv(:,8);
                    gijDown(:,3) = xsInv(:,3).*xsInv(:,3) + xsInv(:,6).*xsInv(:,6) + xsInv(:,9).*xsInv(:,9);
                    gijDown(:,4) = xsInv(:,1).*xsInv(:,2) + xsInv(:,4).*xsInv(:,5) + xsInv(:,7).*xsInv(:,8);
                    gijDown(:,5) = xsInv(:,2).*xsInv(:,3) + xsInv(:,5).*xsInv(:,6) + xsInv(:,8).*xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1).*xsInv(:,3) + xsInv(:,4).*xsInv(:,6) + xsInv(:,7).*xsInv(:,9);
                end
            end
            
            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            
            ind1 = find(localPhi > 1.0) ;
            ind2 = find(localPhi < -1.0) ;
            ind3 = find(localPhi <= 1.0 & localPhi >= -1.0) ;
            dens = zeros(size(localPhi,1),1) ;
            visc = zeros(size(localPhi,1),1) ;
            
            dens(ind1) = fluid.dens1 ;
            visc(ind1) = fluid.visc1 ;
            dens(ind2) = fluid.dens2 ;
            visc(ind2) = fluid.visc2 ;
            dens(ind3) = fluid.dens1*0.5*(1+localPhi(ind3)) + fluid.dens2*0.5*(1-localPhi(ind3)) ;
            visc(ind3) = fluid.visc1*0.5*(1+localPhi(ind3)) + fluid.visc2*0.5*(1-localPhi(ind3)) ;
            clear ind1 ind2 ind3
            
            if (Flag3D ~= 1)
                DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
                DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2); 
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,3) .* locUY ;
                VG2	= gijDown(:,3) .* locUX + gijDown(:,2) .* locUY ;

                VGV	= VG1 .* locUX + VG2 .* locUY ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + 2 * ( gijDown(:,3).^2 );
            else
                DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
                      (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
                      (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
                DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
                      (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
                      (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
                DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
                      (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
                      (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);
                locUZ  = sum(repmat(N(p,:),nElem,1).*uz,2);
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
                VG2	= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
                VG3	= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

                VGV	= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
                    2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
            end
            tauM = (2/solver.dt)^2 + ((visc./dens).^2).* GG * 36 + abs(VGV);    
            tauM = tauM.^(0.5);
            tauM = 1./tauM;
            
            index = 0;
            for i = 1:nen
                for j = 1:nen
                    % 1st term
                    Aij_1 = gW(p)*(locUX.*DNDx(:,i)*N(p,j));
                    Aij_2 = gW(p)*(locUY.*DNDy(:,i)*N(p,j));
                    Aij_1 = Aij_1.*volume.*tauM.*dens;
                    Aij_2 = Aij_2.*volume.*tauM.*dens;
                    sA1(index+1:index+nElem,p) = Aij_1;
                    sA2(index+1:index+nElem,p) = Aij_2;
                    
                    % 2nd term
                    Aij_4 = gW(p)*locUX.*locUX.*(DNDx(:,i).*DNDx(:,j));
                    Aij_5 = gW(p)*locUX.*locUY.*(DNDx(:,i).*DNDy(:,j));
                    Aij_7 = gW(p)*locUY.*locUX.*(DNDy(:,i).*DNDx(:,j));
                    Aij_8 = gW(p)*locUY.*locUY.*(DNDy(:,i).*DNDy(:,j)); 
                    Aij_4 = Aij_4.*volume.*tauM.*dens;
                    Aij_5 = Aij_5.*volume.*tauM.*dens;
                    Aij_7 = Aij_7.*volume.*tauM.*dens;
                    Aij_8 = Aij_8.*volume.*tauM.*dens;
                    sA4(index+1:index+nElem,p) = Aij_4;
                    sA5(index+1:index+nElem,p) = Aij_5;
                    sA7(index+1:index+nElem,p) = Aij_7;
                    sA8(index+1:index+nElem,p) = Aij_8;

                    % 3rd term
                    Aij_13 = gW(p)*(DNDx(:,i).*N(p,j));
                    Aij_14 = gW(p)*(DNDy(:,i).*N(p,j));
                    Aij_13 = Aij_13.*volume.*tauM;
                    Aij_14 = Aij_14.*volume.*tauM; 
                    sA13(index+1:index+nElem,p) = Aij_13;
                    sA14(index+1:index+nElem,p) = Aij_14;
      
                    if (Flag3D == 1)
                    Aij_3 = gW(p)*(locUZ.*DNDz(:,i)*N(p,j));
                    Aij_3 = Aij_3.*volume.*tauM.*dens;
                    sA3(index+1:index+nElem,p) = Aij_3;
                    Aij_6 = gW(p)*locUX.*locUZ.*(DNDx(:,i).*DNDz(:,j));
                    Aij_9 = gW(p)*locUY.*locUZ.*(DNDy(:,i).*DNDz(:,j));
                    Aij_10 = gW(p)*locUZ.*locUX.*(DNDz(:,i).*DNDx(:,j));
                    Aij_11 = gW(p)*locUZ.*locUY.*(DNDz(:,i).*DNDy(:,j));
                    Aij_12 = gW(p)*locUZ.*locUZ.*(DNDz(:,i).*DNDz(:,j));
                    Aij_6 = Aij_6.*volume.*tauM.*dens;
                    Aij_9 = Aij_9.*volume.*tauM.*dens;
                    Aij_10 = Aij_10.*volume.*tauM.*dens;
                    Aij_11 = Aij_11.*volume.*tauM.*dens;
                    Aij_12 = Aij_12.*volume.*tauM.*dens;
                    sA6(index+1:index+nElem,p) = Aij_6;
                    sA9(index+1:index+nElem,p) = Aij_9;
                    sA10(index+1:index+nElem,p) = Aij_10;
                    sA11(index+1:index+nElem,p) = Aij_11;
                    sA12(index+1:index+nElem,p) = Aij_12;
                    Aij_15 = gW(p)*(DNDz(:,i).*N(p,j));
                    Aij_15 = Aij_15.*volume.*tauM;
                    sA15(index+1:index+nElem,p) = Aij_15;
                    end

                    index = index + nElem;
                end
            end
        end
        sA1 = sum(sA1,2);
        sA2 = sum(sA2,2);       
        sA4 = sum(sA4,2);
        sA5 = sum(sA5,2);
        sA7 = sum(sA7,2);
        sA8 = sum(sA8,2);
        sA13 = sum(sA13,2);
        sA14 = sum(sA14,2);

        if (Flag3D == 1)
        sA3 = sum(sA3,2);
        sA6 = sum(sA6,2);
        sA9 = sum(sA9,2);
        sA10 = sum(sA10,2);
        sA11 = sum(sA11,2);
        sA12 = sum(sA12,2);
        sA15 = sum(sA15,2);
        end
        
        % Assemble the matrix        
        A1 = sparse(iif,jjf,sA1,ndof,ndof);
        A2 = sparse(iif,jjf,sA2,ndof,ndof);
        A4 = sparse(iif,jjf,sA4,ndof,ndof);
        A5 = sparse(iif,jjf,sA5,ndof,ndof);
        A7 = sparse(iif,jjf,sA7,ndof,ndof);
        A8 = sparse(iif,jjf,sA8,ndof,ndof);
        A13 = sparse(iif,jjf,sA13,ndof,ndof);
        A14 = sparse(iif,jjf,sA14,ndof,ndof);
        
        if (Flag3D == 1)
        A3 = sparse(iif,jjf,sA3,ndof,ndof);
        A6 = sparse(iif,jjf,sA6,ndof,ndof);
        A9 = sparse(iif,jjf,sA9,ndof,ndof);
        A10 = sparse(iif,jjf,sA10,ndof,ndof);
        A11 = sparse(iif,jjf,sA11,ndof,ndof);
        A12 = sparse(iif,jjf,sA12,ndof,ndof);
        A15 = sparse(iif,jjf,sA15,ndof,ndof);
        end
        
        if (Flag3D == 1)
        stabK1 = [A1+A2+A3 ZeroF ZeroF ZeroF;...
                  ZeroF A1+A2+A3 ZeroF ZeroF;...
                  ZeroF ZeroF A1+A2+A3 ZeroF;...
                  ZeroF ZeroF ZeroF ZeroF];
        stabK2 = [A4+A5+A6+A7+A8+A9+A10+A11+A12 ZeroF ZeroF ZeroF;...
                  ZeroF A4+A5+A6+A7+A8+A9+A10+A11+A12 ZeroF ZeroF;...
                  ZeroF ZeroF A4+A5+A6+A7+A8+A9+A10+A11+A12 ZeroF;...
                  ZeroF ZeroF ZeroF ZeroF];
        stabG1 = [ZeroF ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF ZeroF ZeroF;...
                  A13 A14 A15 ZeroF];
        else
        stabK1 = [A1+A2 ZeroF ZeroF;...
                  ZeroF A1+A2 ZeroF;...
                  ZeroF ZeroF ZeroF];
        stabK2 = [A4+A5+A7+A8 ZeroF ZeroF;...
                  ZeroF A4+A5+A7+A8 ZeroF;...
                  ZeroF ZeroF ZeroF];
        stabG1 = [ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF ZeroF;...
                  A13 A14 ZeroF];
        end
        stabK1Flow = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*stabK1;
        stabG1Flow = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*stabG1;

        MatlabFlow = MatlabFlow + (stabK1Flow + stabK2 + stabG1Flow);
        
        clear sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 sA13 sA14 sA15
        clear Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
        clear Aij_13 Aij_14 Aij_15
        clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 xsInv
        
        RHS = RHS -(stabK1 + stabG1) * uDotAlpha(:) - (stabK2) * uAlpha(:) ;

        % Petrov-Galerkin stabilization terms for Navier-Stokes equations
        % Part 2: (grad.N).tauC.(rho*(grad.U)) + 
        %         (tauM/rho)(rho*u.gradN).(gradP) +
        %         (tauM/rho)(-rho*gravFrc)
        sA1 = zeros(nen^2*nElem,nQuad); 
        sA2 = zeros(nen^2*nElem,nQuad);
        sA3 = zeros(nen^2*nElem,nQuad);
        sA4 = zeros(nen^2*nElem,nQuad);
        sA5 = zeros(nen^2*nElem,nQuad);
        sA6 = zeros(nen^2*nElem,nQuad); 
        
        sA7 = zeros(nen^2*nElem,nQuad); 
        sA8 = zeros(nen^2*nElem,nQuad);
        sA9 = zeros(nen^2*nElem,nQuad);
        sA10 = zeros(nen^2*nElem,nQuad);
        sA11 = zeros(nen^2*nElem,nQuad);
        sA12 = zeros(nen^2*nElem,nQuad);
        sA13 = zeros(nen^2*nElem,nQuad); 
        sA14 = zeros(nen^2*nElem,nQuad);
        sA15 = zeros(nen^2*nElem,nQuad);
        
        sA16 = zeros(nen^2*nElem,nQuad); 
        sA17 = zeros(nen^2*nElem,nQuad);
        sA18 = zeros(nen^2*nElem,nQuad);
        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
                xsInv(:,4) = J(:,1).*J(:,4) - J(:,2).*J(:,3) ;
                xsInv(:,4) = 1./xsInv(:,4) ;
                xsInv(:,1) = xsInv(:,4).* J(:,4) ;
                xsInv(:,2) = xsInv(:,4).* (-J(:,2)) ;
                xsInv(:,3) = xsInv(:,4).* (-J(:,3)) ;
                xsInv(:,4) = xsInv(:,4).* J(:,1) ;
                
                if (strcmp(elemType,'3Tri'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                elseif (strcmp(elemType,'4Quad'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                end
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                      
                xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
                xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
                xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
                xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
                xsInv(:,9) = 1./xsInv(:,9) ;
                xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
                xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
                xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
                xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
                xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
                xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
                xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
                xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
                xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
                if (strcmp(elemType,'4Tet'))
                    c1 = 1.259921049894873E0 ;
                    c2 = 6.299605249474365D-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * (xsInv(:,4) + xsInv(:,7)) ;
                    a2 = c1 * xsInv(:,4) + c2 * (xsInv(:,1) + xsInv(:,7)) ;
                    a3 = c1 * xsInv(:,7) + c2 * (xsInv(:,1) + xsInv(:,4)) ;
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,2) + c2 * (xsInv(:,5) + xsInv(:,8)) ;
                    a2 = c1 * xsInv(:,5) + c2 * (xsInv(:,2) + xsInv(:,8)) ;
                    a3 = c1 * xsInv(:,8) + c2 * (xsInv(:,2) + xsInv(:,5)) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,3) + c2 * (xsInv(:,6) + xsInv(:,9)) ;
                    a2 = c1 * xsInv(:,6) + c2 * (xsInv(:,3) + xsInv(:,9)) ;
                    a3 = c1 * xsInv(:,9) + c2 * (xsInv(:,3) + xsInv(:,6)) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* a3;
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
                elseif (strcmp(elemType,'6Prism'))
                    c1 = 1.154700538379252E0 ;
                    c2 = 5.773502691896259E-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
                    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
                    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
                    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
                    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
                    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
                elseif (strcmp(elemType,'8Hex'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,4).*xsInv(:,4) + xsInv(:,7).*xsInv(:,7);
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,5).*xsInv(:,5) + xsInv(:,8).*xsInv(:,8);
                    gijDown(:,3) = xsInv(:,3).*xsInv(:,3) + xsInv(:,6).*xsInv(:,6) + xsInv(:,9).*xsInv(:,9);
                    gijDown(:,4) = xsInv(:,1).*xsInv(:,2) + xsInv(:,4).*xsInv(:,5) + xsInv(:,7).*xsInv(:,8);
                    gijDown(:,5) = xsInv(:,2).*xsInv(:,3) + xsInv(:,5).*xsInv(:,6) + xsInv(:,8).*xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1).*xsInv(:,3) + xsInv(:,4).*xsInv(:,6) + xsInv(:,7).*xsInv(:,9);
                end
            end
            
            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            
            ind1 = find(localPhi > 1.0) ;
            ind2 = find(localPhi < -1.0) ;
            ind3 = find(localPhi <= 1.0 & localPhi >= -1.0) ;
            dens = zeros(size(localPhi,1),1) ;
            visc = zeros(size(localPhi,1),1) ;
            
            dens(ind1) = fluid.dens1 ;
            visc(ind1) = fluid.visc1 ;
            dens(ind2) = fluid.dens2 ;
            visc(ind2) = fluid.visc2 ;
            dens(ind3) = fluid.dens1*0.5*(1+localPhi(ind3)) + fluid.dens2*0.5*(1-localPhi(ind3)) ;
            visc(ind3) = fluid.visc1*0.5*(1+localPhi(ind3)) + fluid.visc2*0.5*(1-localPhi(ind3)) ;
            clear ind1 ind2 ind3
            
            if (Flag3D ~= 1)
                DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
                DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2); 
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,3) .* locUY ;
                VG2	= gijDown(:,3) .* locUX + gijDown(:,2) .* locUY ;

                VGV	= VG1 .* locUX + VG2 .* locUY ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + 2 * ( gijDown(:,3).^2 );
            else
                DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
                      (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
                      (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
                DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
                      (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
                      (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
                DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
                      (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
                      (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);
                locUZ  = sum(repmat(N(p,:),nElem,1).*uz,2);
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
                VG2	= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
                VG3	= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

                VGV	= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
                    2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
            end
            tauM = (2/solver.dt)^2 + ((visc./dens).^2).* GG * 36 + abs(VGV);    
            tauM = tauM.^(0.5);
            tauC = tauM ./((gijDown(:,1) + gijDown(:,2) + gijDown(:,3)));
            tauM = 1./tauM;
            
            index = 0;
            for i = 1:nen
                for j = 1:nen
                    % 1st term
                    Aij_1 = gW(p)*(DNDx(:,i).*DNDx(:,j));
                    Aij_2 = gW(p)*(DNDy(:,i).*DNDy(:,j)); 
                    Aij_4 = gW(p)*(DNDx(:,i).*DNDy(:,j));
                    Aij_1 = Aij_1.*volume.*tauC.*dens;
                    Aij_2 = Aij_2.*volume.*tauC.*dens;
                    Aij_4 = Aij_4.*volume.*tauC.*dens;
                    sA1(index+1:index+nElem,p) = Aij_1;
                    sA2(index+1:index+nElem,p) = Aij_2;
                    sA4(index+1:index+nElem,p) = Aij_4;
                    
                    % 2nd term
                    Aij_7 = gW(p)*locUX.*(DNDx(:,i).*DNDx(:,j));
                    Aij_8 = gW(p)*locUX.*(DNDx(:,i).*DNDy(:,j));
                    Aij_10 = gW(p)*locUY.*(DNDy(:,i).*DNDx(:,j));
                    Aij_11 = gW(p)*locUY.*(DNDy(:,i).*DNDy(:,j));
                    Aij_7 = Aij_7.*volume.*tauM;
                    Aij_8 = Aij_8.*volume.*tauM;
                    Aij_10 = Aij_10.*volume.*tauM;
                    Aij_11 = Aij_11.*volume.*tauM;
                    sA7(index+1:index+nElem,p) = Aij_7;
                    sA8(index+1:index+nElem,p) = Aij_8;
                    sA10(index+1:index+nElem,p) = Aij_10;
                    sA11(index+1:index+nElem,p) = Aij_11;
                    
                    % 3rd term
                    Aij_16 = gW(p)*(locUX.*DNDx(:,i))*N(p,j) ;
                    Aij_17 = gW(p)*(locUY.*DNDy(:,i))*N(p,j) ;        
                    Aij_16 = Aij_16.*volume.*tauM.*dens ;
                    Aij_17 = Aij_17.*volume.*tauM.*dens ;      
                    sA16(index+1:index+nElem,p) = Aij_16;
                    sA17(index+1:index+nElem,p) = Aij_17;     

                    if (Flag3D == 1)
                    Aij_3 = gW(p)*(DNDz(:,i).*DNDz(:,j));
                    Aij_5 = gW(p)*(DNDy(:,i).*DNDz(:,j));
                    Aij_6 = gW(p)*(DNDz(:,i).*DNDx(:,j));
                    Aij_3 = Aij_3.*volume.*tauC.*dens;
                    Aij_5 = Aij_5.*volume.*tauC.*dens;
                    Aij_6 = Aij_6.*volume.*tauC.*dens;
                    sA3(index+1:index+nElem,p) = Aij_3;
                    sA5(index+1:index+nElem,p) = Aij_5;
                    sA6(index+1:index+nElem,p) = Aij_6;
                    Aij_9 = gW(p)*locUX.*(DNDx(:,i).*DNDz(:,j));
                    Aij_12 = gW(p)*locUY.*(DNDy(:,i).*DNDz(:,j));
                    Aij_13 = gW(p)*locUZ.*(DNDz(:,i).*DNDx(:,j));
                    Aij_14 = gW(p)*locUZ.*(DNDz(:,i).*DNDy(:,j));
                    Aij_15 = gW(p)*locUZ.*(DNDz(:,i).*DNDz(:,j));
                    Aij_9 = Aij_9.*volume.*tauM;
                    Aij_12 = Aij_12.*volume.*tauM;
                    Aij_13 = Aij_13.*volume.*tauM;
                    Aij_14 = Aij_14.*volume.*tauM;
                    Aij_15 = Aij_15.*volume.*tauM;
                    sA9(index+1:index+nElem,p) = Aij_9;
                    sA12(index+1:index+nElem,p) = Aij_12;
                    sA13(index+1:index+nElem,p) = Aij_13;
                    sA14(index+1:index+nElem,p) = Aij_14;
                    sA15(index+1:index+nElem,p) = Aij_15;
                    Aij_18 = gW(p)*(locUZ.*DNDz(:,i))*N(p,j) ;
                    Aij_18 = Aij_18.*volume.*tauM.*dens ;
                    sA18(index+1:index+nElem,p) = Aij_18;
                    end

                    index = index + nElem;
                end
            end
        end
        sA1 = sum(sA1,2);
        sA2 = sum(sA2,2);       
        sA4 = sum(sA4,2);
        sA7 = sum(sA7,2);
        sA8 = sum(sA8,2);
        sA10 = sum(sA10,2);
        sA11 = sum(sA11,2);
        sA16 = sum(sA16,2);
        sA17 = sum(sA17,2);   

        if (Flag3D == 1)
        sA3 = sum(sA3,2);
        sA5 = sum(sA5,2);
        sA6 = sum(sA6,2);
        sA9 = sum(sA9,2);
        sA12 = sum(sA12,2);
        sA13 = sum(sA13,2);
        sA14 = sum(sA14,2);
        sA15 = sum(sA15,2);
        sA18 = sum(sA18,2);
        end
        
        % Assemble the matrix        
        A1 = sparse(iif,jjf,sA1,ndof,ndof);
        A2 = sparse(iif,jjf,sA2,ndof,ndof);
        A4 = sparse(iif,jjf,sA4,ndof,ndof);
        A7 = sparse(iif,jjf,sA7,ndof,ndof);
        A8 = sparse(iif,jjf,sA8,ndof,ndof);
        A10 = sparse(iif,jjf,sA10,ndof,ndof);
        A11 = sparse(iif,jjf,sA11,ndof,ndof);
        A16 = sparse(iif,jjf,sA16,ndof,ndof);
        A17 = sparse(iif,jjf,sA17,ndof,ndof);      
        
        if (Flag3D == 1)
        A3 = sparse(iif,jjf,sA3,ndof,ndof);
        A5 = sparse(iif,jjf,sA5,ndof,ndof);
        A6 = sparse(iif,jjf,sA6,ndof,ndof);
        A9 = sparse(iif,jjf,sA9,ndof,ndof);
        A12 = sparse(iif,jjf,sA12,ndof,ndof);
        A13 = sparse(iif,jjf,sA13,ndof,ndof);
        A14 = sparse(iif,jjf,sA14,ndof,ndof);
        A15 = sparse(iif,jjf,sA15,ndof,ndof);
        A18 = sparse(iif,jjf,sA18,ndof,ndof);
        end
        
        if (Flag3D == 1)
        stabC1 = [A1 A4 A6' ZeroF;...
                  A4' A2 A5 ZeroF;...
                  A6 A5' A3 ZeroF;...
                  ZeroF ZeroF ZeroF ZeroF];
        stabF1 = [ZeroF ZeroF ZeroF A7+A10+A13;...
                  ZeroF ZeroF ZeroF A8+A11+A14;...
                  ZeroF ZeroF ZeroF A9+A12+A15;...
                  ZeroF ZeroF ZeroF ZeroF];
        Src3 = [(A16+A17+A18).*fluid.gravFrc(1) ZeroF ZeroF ZeroF;...
                ZeroF (A16+A17+A18).*fluid.gravFrc(2) ZeroF ZeroF;...
                ZeroF ZeroF (A16+A17+A18).*fluid.gravFrc(3) ZeroF;...
                ZeroF ZeroF ZeroF ZeroF];
        else
        stabC1 = [A1 A4 ZeroF;...
                  A4' A2 ZeroF;...
                  ZeroF ZeroF ZeroF];
        stabF1 = [ZeroF ZeroF A7+A10;...
                  ZeroF ZeroF A8+A11;...
                  ZeroF ZeroF ZeroF];
        Src3 = [(A16+A17).*fluid.gravFrc(1) ZeroF ZeroF;...
                ZeroF (A16+A17).*fluid.gravFrc(2) ZeroF;...
                ZeroF ZeroF ZeroF];
        end

        MatlabFlow = MatlabFlow + stabC1 + (stabF1) ;
        
        clear sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 sA13 sA14 sA15 sA16 sA17 sA18
        clear Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
        clear Aij_13 Aij_14 Aij_15 Aij_16 Aij_17 Aij_18
        clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16 A17 A18 xsInv
        
        RHS = RHS - stabC1 * uAlpha(:) ;
        RHS = RHS - stabF1 * p1 + Src3 *gravVec ;
        
        % Petrov-Galerkin stabilization terms for Navier-Stokes equations
        % Part 3: (tauM/rho)(gradN).(rho*u.(gradU)) + 
        %         (tauM/rho)(gradN).(gradP) +
        %         (tauM/rho)(gradN).(-rho.gravVec)
        sA1 = zeros(nen^2*nElem,nQuad); 
        sA2 = zeros(nen^2*nElem,nQuad);
        sA3 = zeros(nen^2*nElem,nQuad);
        sA4 = zeros(nen^2*nElem,nQuad);
        sA5 = zeros(nen^2*nElem,nQuad);
        sA6 = zeros(nen^2*nElem,nQuad); 
        sA7 = zeros(nen^2*nElem,nQuad); 
        sA8 = zeros(nen^2*nElem,nQuad);
        sA9 = zeros(nen^2*nElem,nQuad);
        
        sA10 = zeros(nen^2*nElem,nQuad);
        sA11 = zeros(nen^2*nElem,nQuad);
        sA12 = zeros(nen^2*nElem,nQuad);
        
        sA13 = zeros(nen^2*nElem,nQuad);
        sA14 = zeros(nen^2*nElem,nQuad);
        sA15 = zeros(nen^2*nElem,nQuad);
        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
                xsInv(:,4) = J(:,1).*J(:,4) - J(:,2).*J(:,3) ;
                xsInv(:,4) = 1./xsInv(:,4) ;
                xsInv(:,1) = xsInv(:,4).* J(:,4) ;
                xsInv(:,2) = xsInv(:,4).* (-J(:,2)) ;
                xsInv(:,3) = xsInv(:,4).* (-J(:,3)) ;
                xsInv(:,4) = xsInv(:,4).* J(:,1) ;
                
                if (strcmp(elemType,'3Tri'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                elseif (strcmp(elemType,'4Quad'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                end
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                      
                xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
                xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
                xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
                xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
                xsInv(:,9) = 1./xsInv(:,9) ;
                xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
                xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
                xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
                xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
                xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
                xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
                xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
                xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
                xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
                if (strcmp(elemType,'4Tet'))
                    c1 = 1.259921049894873E0 ;
                    c2 = 6.299605249474365D-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * (xsInv(:,4) + xsInv(:,7)) ;
                    a2 = c1 * xsInv(:,4) + c2 * (xsInv(:,1) + xsInv(:,7)) ;
                    a3 = c1 * xsInv(:,7) + c2 * (xsInv(:,1) + xsInv(:,4)) ;
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,2) + c2 * (xsInv(:,5) + xsInv(:,8)) ;
                    a2 = c1 * xsInv(:,5) + c2 * (xsInv(:,2) + xsInv(:,8)) ;
                    a3 = c1 * xsInv(:,8) + c2 * (xsInv(:,2) + xsInv(:,5)) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,3) + c2 * (xsInv(:,6) + xsInv(:,9)) ;
                    a2 = c1 * xsInv(:,6) + c2 * (xsInv(:,3) + xsInv(:,9)) ;
                    a3 = c1 * xsInv(:,9) + c2 * (xsInv(:,3) + xsInv(:,6)) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* a3;
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
                elseif (strcmp(elemType,'6Prism'))
                    c1 = 1.154700538379252E0 ;
                    c2 = 5.773502691896259E-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
                    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
                    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
                    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
                    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
                    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
                elseif (strcmp(elemType,'8Hex'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,4).*xsInv(:,4) + xsInv(:,7).*xsInv(:,7);
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,5).*xsInv(:,5) + xsInv(:,8).*xsInv(:,8);
                    gijDown(:,3) = xsInv(:,3).*xsInv(:,3) + xsInv(:,6).*xsInv(:,6) + xsInv(:,9).*xsInv(:,9);
                    gijDown(:,4) = xsInv(:,1).*xsInv(:,2) + xsInv(:,4).*xsInv(:,5) + xsInv(:,7).*xsInv(:,8);
                    gijDown(:,5) = xsInv(:,2).*xsInv(:,3) + xsInv(:,5).*xsInv(:,6) + xsInv(:,8).*xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1).*xsInv(:,3) + xsInv(:,4).*xsInv(:,6) + xsInv(:,7).*xsInv(:,9);
                end
            end
            
            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            
            ind1 = find(localPhi > 1.0) ;
            ind2 = find(localPhi < -1.0) ;
            ind3 = find(localPhi <= 1.0 & localPhi >= -1.0) ;
            dens = zeros(size(localPhi,1),1) ;
            visc = zeros(size(localPhi,1),1) ;
            
            dens(ind1) = fluid.dens1 ;
            visc(ind1) = fluid.visc1 ;
            dens(ind2) = fluid.dens2 ;
            visc(ind2) = fluid.visc2 ;
            dens(ind3) = fluid.dens1*0.5*(1+localPhi(ind3)) + fluid.dens2*0.5*(1-localPhi(ind3)) ;
            visc(ind3) = fluid.visc1*0.5*(1+localPhi(ind3)) + fluid.visc2*0.5*(1-localPhi(ind3)) ;
            clear ind1 ind2 ind3
            
            if (Flag3D ~= 1)
                DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
                DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2); 
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,3) .* locUY ;
                VG2	= gijDown(:,3) .* locUX + gijDown(:,2) .* locUY ;

                VGV	= VG1 .* locUX + VG2 .* locUY ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + 2 * ( gijDown(:,3).^2 );
            else
                DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
                      (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
                      (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
                DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
                      (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
                      (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
                DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
                      (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
                      (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);
                locUZ  = sum(repmat(N(p,:),nElem,1).*uz,2);
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
                VG2	= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
                VG3	= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

                VGV	= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
                    2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
            end
            tauM = (2/solver.dt)^2 + ((visc./dens).^2).* GG * 36 + abs(VGV);    
            tauM = tauM.^(0.5);
            tauM = 1./tauM;
            
            index = 0;
            for i = 1:nen
                for j = 1:nen
                    % 1st term
                    Aij_1 = gW(p)*locUX.*(DNDx(:,i).*DNDx(:,j));
                    Aij_2 = gW(p)*locUY.*(DNDx(:,i).*DNDy(:,j));
                    Aij_4 = gW(p)*locUX.*(DNDy(:,i).*DNDx(:,j));
                    Aij_5 = gW(p)*locUY.*(DNDy(:,i).*DNDy(:,j));
                    Aij_1 = Aij_1.*volume.*tauM;
                    Aij_2 = Aij_2.*volume.*tauM;
                    Aij_4 = Aij_4.*volume.*tauM;
                    Aij_5 = Aij_5.*volume.*tauM;
                    sA1(index+1:index+nElem,p) = Aij_1;
                    sA2(index+1:index+nElem,p) = Aij_2;
                    sA4(index+1:index+nElem,p) = Aij_4;
                    sA5(index+1:index+nElem,p) = Aij_5;

                    % 2nd term
                    Aij_10 = gW(p)*(DNDx(:,i).*DNDx(:,j));
                    Aij_11 = gW(p)*(DNDy(:,i).*DNDy(:,j));
                    Aij_10 = Aij_10.*volume.*tauM.*(1./dens);
                    Aij_11 = Aij_11.*volume.*tauM.*(1./dens);
                    sA10(index+1:index+nElem,p) = Aij_10;
                    sA11(index+1:index+nElem,p) = Aij_11;
                    
                    % 3rd term
                    Aij_13 = gW(p)*(DNDx(:,i)*N(p,j));
                    Aij_14 = gW(p)*(DNDy(:,i)*N(p,j));
                    Aij_13 = Aij_13.*volume.*tauM ;
                    Aij_14 = Aij_14.*volume.*tauM ;
                    sA13(index+1:index+nElem,p) = Aij_13;
                    sA14(index+1:index+nElem,p) = Aij_14;
                    
                    if (Flag3D == 1)
                    Aij_3 = gW(p)*locUZ.*(DNDx(:,i).*DNDz(:,j));
                    Aij_6 = gW(p)*locUZ.*(DNDy(:,i).*DNDz(:,j));
                    Aij_7 = gW(p)*locUX.*(DNDz(:,i).*DNDx(:,j));
                    Aij_8 = gW(p)*locUY.*(DNDz(:,i).*DNDy(:,j));
                    Aij_9 = gW(p)*locUZ.*(DNDz(:,i).*DNDz(:,j));
                    Aij_3 = Aij_3.*volume.*tauM;
                    Aij_6 = Aij_6.*volume.*tauM;
                    Aij_7 = Aij_7.*volume.*tauM;
                    Aij_8 = Aij_8.*volume.*tauM;
                    Aij_9 = Aij_9.*volume.*tauM;
                    sA3(index+1:index+nElem,p) = Aij_3;
                    sA6(index+1:index+nElem,p) = Aij_6;
                    sA7(index+1:index+nElem,p) = Aij_7;
                    sA8(index+1:index+nElem,p) = Aij_8;
                    sA9(index+1:index+nElem,p) = Aij_9;
                    Aij_12 = gW(p)*(DNDz(:,i).*DNDz(:,j));
                    Aij_12 = Aij_12.*volume.*tauM.*(1./dens);
                    sA12(index+1:index+nElem,p) = Aij_12;
                    Aij_15 = gW(p)*(DNDz(:,i)*N(p,j));
                    Aij_15 = Aij_15.*volume.*tauM ;
                    sA15(index+1:index+nElem,p) = Aij_15;
                    end

                    index = index + nElem;
                end
            end
        end
        sA1 = sum(sA1,2);
        sA2 = sum(sA2,2);       
        sA4 = sum(sA4,2);
        sA5 = sum(sA5,2);
        sA10 = sum(sA10,2);
        sA11 = sum(sA11,2);
        sA13 = sum(sA13,2);
        sA14 = sum(sA14,2);

        if (Flag3D == 1)
        sA3 = sum(sA3,2);
        sA6 = sum(sA6,2);
        sA7 = sum(sA7,2);
        sA8 = sum(sA8,2);
        sA9 = sum(sA9,2);
        sA12 = sum(sA12,2);
        sA15 = sum(sA15,2);
        end
        
        % Assemble the matrix        
        A1 = sparse(iif,jjf,sA1,ndof,ndof);
        A2 = sparse(iif,jjf,sA2,ndof,ndof);
        A4 = sparse(iif,jjf,sA4,ndof,ndof);
        A5 = sparse(iif,jjf,sA5,ndof,ndof);
        A10 = sparse(iif,jjf,sA10,ndof,ndof);
        A11 = sparse(iif,jjf,sA11,ndof,ndof);
        A13 = sparse(iif,jjf,sA13,ndof,ndof);
        A14 = sparse(iif,jjf,sA14,ndof,ndof);
        
        if (Flag3D == 1)
        A3 = sparse(iif,jjf,sA3,ndof,ndof);
        A6 = sparse(iif,jjf,sA6,ndof,ndof);
        A7 = sparse(iif,jjf,sA7,ndof,ndof);
        A8 = sparse(iif,jjf,sA8,ndof,ndof);
        A9 = sparse(iif,jjf,sA9,ndof,ndof);
        A12 = sparse(iif,jjf,sA12,ndof,ndof);
        A15 = sparse(iif,jjf,sA15,ndof,ndof);
        end
        
        if (Flag3D == 1)
        stabG2 = [ZeroF ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF ZeroF ZeroF;...
                  A1+A2+A3 A4+A5+A6 A7+A8+A9 ZeroF];
        stabC1 = [ZeroF ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF ZeroF A10+A11+A12];
        Src4 = [ZeroF ZeroF ZeroF ZeroF;...
                ZeroF ZeroF ZeroF ZeroF;...
                ZeroF ZeroF ZeroF ZeroF;...
                A13.*fluid.gravFrc(1) A14.*fluid.gravFrc(2) A15.*fluid.gravFrc(3) ZeroF] ;
        else
        stabG2 = [ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF ZeroF;...
                  A1+A2 A4+A5 ZeroF];
        stabC1 = [ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF ZeroF;...
                  ZeroF ZeroF A10+A11];
        Src4 = [ZeroF ZeroF ZeroF;...
                ZeroF ZeroF ZeroF;...
                A13.*fluid.gravFrc(1) A14.*fluid.gravFrc(2) ZeroF] ;
        end

        MatlabFlow = MatlabFlow + (stabG2 + stabC1);
        
        clear sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 
        clear Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
        clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12
        
        RHS = RHS - stabG2 * uAlpha(:) ;
        RHS = RHS - stabC1 * p1 + Src4 * gravVec;
      
        % Discontinuity capturing stabilization terms for Navier-Stokes equations
        % Part 4: 
        %
        sA1 = zeros(nen^2*nElem,nQuad); 
        sA2 = zeros(nen^2*nElem,nQuad);
        sA3 = zeros(nen^2*nElem,nQuad);
        sA4 = zeros(nen^2*nElem,nQuad);
        sA5 = zeros(nen^2*nElem,nQuad);
        sA6 = zeros(nen^2*nElem,nQuad); 
        sA7 = zeros(nen^2*nElem,nQuad); 
        sA8 = zeros(nen^2*nElem,nQuad);
        sA9 = zeros(nen^2*nElem,nQuad);
        
        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
                xsInv(:,4) = J(:,1).*J(:,4) - J(:,2).*J(:,3) ;
                xsInv(:,4) = 1./xsInv(:,4) ;
                xsInv(:,1) = xsInv(:,4).* J(:,4) ;
                xsInv(:,2) = xsInv(:,4).* (-J(:,2)) ;
                xsInv(:,3) = xsInv(:,4).* (-J(:,3)) ;
                xsInv(:,4) = xsInv(:,4).* J(:,1) ;
                
                if (strcmp(elemType,'3Tri'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;

                    gijUp(:,1) = J(:,1).*J(:,1) + J(:,2).*J(:,2) ;
                    gijUp(:,2) = J(:,3).*J(:,3) + J(:,4).*J(:,4) ;
                    gijUp(:,3) = J(:,1).*J(:,3) + J(:,2).*J(:,4) ;
                elseif (strcmp(elemType,'4Quad'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                    gijUp(:,1) = J(:,1).*J(:,1) + J(:,2).*J(:,2) ;
                    gijUp(:,2) = J(:,3).*J(:,3) + J(:,4).*J(:,4) ;
                    gijUp(:,3) = J(:,1).*J(:,3) + J(:,2).*J(:,4) ;
                end
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                      
                xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
                xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
                xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
                xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
                xsInv(:,9) = 1./xsInv(:,9) ;
                xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
                xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
                xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
                xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
                xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
                xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
                xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
                xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
                xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
                if (strcmp(elemType,'4Tet'))
                    c1 = 1.259921049894873E0 ;
                    c2 = 6.299605249474365D-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * (xsInv(:,4) + xsInv(:,7)) ;
                    a2 = c1 * xsInv(:,4) + c2 * (xsInv(:,1) + xsInv(:,7)) ;
                    a3 = c1 * xsInv(:,7) + c2 * (xsInv(:,1) + xsInv(:,4)) ;
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,2) + c2 * (xsInv(:,5) + xsInv(:,8)) ;
                    a2 = c1 * xsInv(:,5) + c2 * (xsInv(:,2) + xsInv(:,8)) ;
                    a3 = c1 * xsInv(:,8) + c2 * (xsInv(:,2) + xsInv(:,5)) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,3) + c2 * (xsInv(:,6) + xsInv(:,9)) ;
                    a2 = c1 * xsInv(:,6) + c2 * (xsInv(:,3) + xsInv(:,9)) ;
                    a3 = c1 * xsInv(:,9) + c2 * (xsInv(:,3) + xsInv(:,6)) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* a3;
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
                elseif (strcmp(elemType,'6Prism'))
                    c1 = 1.154700538379252E0 ;
                    c2 = 5.773502691896259E-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
                    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
                    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
                    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
                    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
                    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
                elseif (strcmp(elemType,'8Hex'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,4).*xsInv(:,4) + xsInv(:,7).*xsInv(:,7);
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,5).*xsInv(:,5) + xsInv(:,8).*xsInv(:,8);
                    gijDown(:,3) = xsInv(:,3).*xsInv(:,3) + xsInv(:,6).*xsInv(:,6) + xsInv(:,9).*xsInv(:,9);
                    gijDown(:,4) = xsInv(:,1).*xsInv(:,2) + xsInv(:,4).*xsInv(:,5) + xsInv(:,7).*xsInv(:,8);
                    gijDown(:,5) = xsInv(:,2).*xsInv(:,3) + xsInv(:,5).*xsInv(:,6) + xsInv(:,8).*xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1).*xsInv(:,3) + xsInv(:,4).*xsInv(:,6) + xsInv(:,7).*xsInv(:,9);
                    
                    gijUp(:,1) = J(:,1).*J(:,1) + J(:,2).*J(:,2) + J(:,3).*J(:,3);
                    gijUp(:,2) = J(:,4).*J(:,4) + J(:,5).*J(:,5) + J(:,6).*J(:,6);
                    gijUp(:,3) = J(:,7).*J(:,7) + J(:,8).*J(:,8) + J(:,9).*J(:,9);
                    gijUp(:,4) = J(:,1).*J(:,4) + J(:,2).*J(:,5) + J(:,3).*J(:,6);
                    gijUp(:,5) = J(:,4).*J(:,7) + J(:,5).*J(:,8) + J(:,6).*J(:,9);
                    gijUp(:,6) = J(:,1).*J(:,7) + J(:,2).*J(:,8) + J(:,3).*J(:,9);
                end
            end
            
            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            
            ind1 = find(localPhi > 1.0) ;
            ind2 = find(localPhi < -1.0) ;
            ind3 = find(localPhi <= 1.0 & localPhi >= -1.0) ;
            dens = zeros(size(localPhi,1),1) ;
            visc = zeros(size(localPhi,1),1) ;
            
            dens(ind1) = fluid.dens1 ;
            visc(ind1) = fluid.visc1 ;
            dens(ind2) = fluid.dens2 ;
            visc(ind2) = fluid.visc2 ;
            dens(ind3) = fluid.dens1*0.5*(1+localPhi(ind3)) + fluid.dens2*0.5*(1-localPhi(ind3)) ;
            visc(ind3) = fluid.visc1*0.5*(1+localPhi(ind3)) + fluid.visc2*0.5*(1-localPhi(ind3)) ;
            clear ind1 ind2 ind3
            
            if (Flag3D ~= 1)
                DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
                DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2); 
                locUXDot  = sum(repmat(N(p,:),nElem,1).*uxDot,2);
                locUYDot  = sum(repmat(N(p,:),nElem,1).*uyDot,2); 
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,3) .* locUY ;
                VG2	= gijDown(:,3) .* locUX + gijDown(:,2) .* locUY ;

                VGV	= VG1 .* locUX + VG2 .* locUY ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + 2 * ( gijDown(:,3).^2 );
            else
                DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
                      (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
                      (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
                DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
                      (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
                      (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
                DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
                      (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
                      (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);
                locUZ  = sum(repmat(N(p,:),nElem,1).*uz,2);
                
                locUXDot  = sum(repmat(N(p,:),nElem,1).*uxDot,2);
                locUYDot  = sum(repmat(N(p,:),nElem,1).*uyDot,2);
                locUZDot  = sum(repmat(N(p,:),nElem,1).*uzDot,2);
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
                VG2	= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
                VG3	= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

                VGV	= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
                    2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
            end
            tauM = (2/solver.dt)^2 + ((visc./dens).^2).* GG * 36 + abs(VGV);    
            tauM = tauM.^(0.5);
            tauM = 1./tauM;

            index = 0;
            for i = 1:nen
                for j = 1:nen
                    
                    locGradPX  = sum(DNDx.*pres,2);
                    locGradPY  = sum(DNDy.*pres,2);

                    locGradU(:,1)  = sum(DNDx.*ux,2);
                    locGradU(:,2)  = sum(DNDx.*uy,2);
                    locGradU(:,4)  = sum(DNDy.*ux,2);
                    locGradU(:,5)  = sum(DNDy.*uy,2);
                    if (Flag3D==1)
                    locGradPZ  = sum(DNDz.*pres,2);
                    locGradU(:,3)  = sum(DNDx.*uz,2);
                    locGradU(:,6)  = sum(DNDy.*uz,2);
                    locGradU(:,7)  = sum(DNDz.*ux,2);
                    locGradU(:,8)  = sum(DNDz.*uy,2);
                    locGradU(:,9)  = sum(DNDz.*uz,2);
                    ResU(:,1) = dens.*locUXDot + dens.*(locUX.*locGradU(:,1) + locUY.*locGradU(:,4) + ...
                     locUZ.*locGradU(:,7)) + locGradPX - dens.*fluid.gravFrc(1) ;
                    ResU(:,2) = dens.*locUYDot + dens.*(locUX.*locGradU(:,2) + locUY.*locGradU(:,5) + ...
                     locUZ.*locGradU(:,8)) + locGradPY - dens.*fluid.gravFrc(2) ;
                    ResU(:,3) = dens.*locUZDot + dens.*(locUX.*locGradU(:,3) + locUY.*locGradU(:,6) + ...
                     locUZ.*locGradU(:,9)) + locGradPZ - dens.*fluid.gravFrc(3) ;
                    else
                    ResU(:,1) = dens.*locUXDot + dens.*(locUX.*locGradU(:,1) + locUY.*locGradU(:,4)) + ...
                      locGradPX - dens.*fluid.gravFrc(1) ;
                    ResU(:,2) = dens.*locUYDot + dens.*(locUX.*locGradU(:,2) + locUY.*locGradU(:,5)) + ...
                      locGradPY - dens.*fluid.gravFrc(2) ; 
                    end
            
                    if (Flag3D==1)
                    gGradV(:,1) = gijUp(:,1).*locGradU(:,1) + gijUp(:,4).*locGradU(:,4) + gijUp(:,6).*locGradU(:,7) ;
                    gGradV(:,2) = gijUp(:,4).*locGradU(:,1) + gijUp(:,2).*locGradU(:,4) + gijUp(:,5).*locGradU(:,7) ;
                    gGradV(:,3) = gijUp(:,6).*locGradU(:,1) + gijUp(:,5).*locGradU(:,4) + gijUp(:,3).*locGradU(:,7) ;
                    gGradV(:,4) = gijUp(:,1).*locGradU(:,2) + gijUp(:,4).*locGradU(:,5) + gijUp(:,6).*locGradU(:,8) ;
                    gGradV(:,5) = gijUp(:,4).*locGradU(:,2) + gijUp(:,2).*locGradU(:,5) + gijUp(:,5).*locGradU(:,8) ;
                    gGradV(:,6) = gijUp(:,6).*locGradU(:,2) + gijUp(:,5).*locGradU(:,5) + gijUp(:,3).*locGradU(:,8) ;
                    gGradV(:,7) = gijUp(:,1).*locGradU(:,3) + gijUp(:,4).*locGradU(:,6) + gijUp(:,6).*locGradU(:,9) ;
                    gGradV(:,8) = gijUp(:,4).*locGradU(:,3) + gijUp(:,2).*locGradU(:,6) + gijUp(:,5).*locGradU(:,9) ;
                    gGradV(:,9) = gijUp(:,6).*locGradU(:,3) + gijUp(:,5).*locGradU(:,6) + gijUp(:,3).*locGradU(:,9) ;
            
                    dcFct = locGradU(:,1).*gGradV(:,1) + locGradU(:,4).*gGradV(:,2) ...
                        + locGradU(:,7).*gGradV(:,3) + locGradU(:,2).*gGradV(:,4) ...
                        + locGradU(:,5).*gGradV(:,5) + locGradU(:,8).*gGradV(:,6) ...
                        + locGradU(:,3).*gGradV(:,7) + locGradU(:,6).*gGradV(:,8) ...
                        + locGradU(:,9).*gGradV(:,9) +eps;
                    else
                    gGradV(:,1) = gijUp(:,1).*locGradU(:,1) + gijUp(:,3).*locGradU(:,4)  ;
                    gGradV(:,2) = gijUp(:,3).*locGradU(:,1) + gijUp(:,2).*locGradU(:,4)  ;
                    gGradV(:,4) = gijUp(:,1).*locGradU(:,2) + gijUp(:,3).*locGradU(:,5)  ;
                    gGradV(:,5) = gijUp(:,3).*locGradU(:,2) + gijUp(:,2).*locGradU(:,5)  ;
            
                    dcFct = locGradU(:,1).*gGradV(:,1) + locGradU(:,4).*gGradV(:,2) ...
                        + locGradU(:,2).*gGradV(:,4) + locGradU(:,5).*gGradV(:,5) + eps ;  
                    end
              
                    dcFct(dcFct>eps) = 1./dcFct(dcFct>eps) ;
                    dcFct(dcFct<=eps) = 0.0 ;
                    if (Flag3D==1)
                    dcFct = dcFct.*( ResU(:,1).*ResU(:,1) + ResU(:,2).*ResU(:,2) + ResU(:,3).*ResU(:,3) ) ;
                    else
                    dcFct = dcFct.*( ResU(:,1).*ResU(:,1) + ResU(:,2).*ResU(:,2) ) ;
                    end
                    dcQuad = tauM.*dcFct./(dens.^2) ;
                    dcConst = 1./(3 * tauM) ;
                    dcFct = min( dcConst, 2.*dcQuad ) ;
                    
                    
                    if (Flag3D == 1)
                    Aij_1 = gW(p)*(DNDx(:,i).*gijUp(:,1).*DNDx(:,j));
                    Aij_2 = gW(p)*(DNDx(:,i).*gijUp(:,4).*DNDy(:,j));
                    Aij_4 = gW(p)*(DNDy(:,i).*gijUp(:,4).*DNDx(:,j));
                    Aij_5 = gW(p)*(DNDy(:,i).*gijUp(:,2).*DNDy(:,j));
                    Aij_1 = Aij_1.*volume.*dens.*dcFct;
                    Aij_2 = Aij_2.*volume.*dens.*dcFct;
                    Aij_4 = Aij_4.*volume.*dens.*dcFct;
                    Aij_5 = Aij_5.*volume.*dens.*dcFct;
                    sA1(index+1:index+nElem,p) = Aij_1;
                    sA2(index+1:index+nElem,p) = Aij_2;
                    sA4(index+1:index+nElem,p) = Aij_4;
                    sA5(index+1:index+nElem,p) = Aij_5;
                    Aij_3 = gW(p)*(DNDx(:,i).*gijUp(:,6).*DNDz(:,j));
                    Aij_6 = gW(p)*(DNDy(:,i).*gijUp(:,5).*DNDz(:,j));
                    Aij_7 = gW(p)*(DNDz(:,i).*gijUp(:,6).*DNDx(:,j));
                    Aij_8 = gW(p)*(DNDz(:,i).*gijUp(:,5).*DNDy(:,j));
                    Aij_9 = gW(p)*(DNDz(:,i).*gijUp(:,3).*DNDz(:,j));
                    Aij_3 = Aij_3.*volume.*dens.*dcFct;
                    Aij_6 = Aij_6.*volume.*dens.*dcFct;
                    Aij_7 = Aij_7.*volume.*dens.*dcFct;
                    Aij_8 = Aij_8.*volume.*dens.*dcFct;
                    Aij_9 = Aij_9.*volume.*dens.*dcFct;
                    sA3(index+1:index+nElem,p) = Aij_3;
                    sA6(index+1:index+nElem,p) = Aij_6;
                    sA7(index+1:index+nElem,p) = Aij_7;
                    sA8(index+1:index+nElem,p) = Aij_8;
                    sA9(index+1:index+nElem,p) = Aij_9;
                    else
                    Aij_1 = gW(p)*(DNDx(:,i).*gijUp(:,1).*DNDx(:,j));
                    Aij_2 = gW(p)*(DNDx(:,i).*gijUp(:,3).*DNDy(:,j));
                    Aij_4 = gW(p)*(DNDy(:,i).*gijUp(:,3).*DNDx(:,j));
                    Aij_5 = gW(p)*(DNDy(:,i).*gijUp(:,2).*DNDy(:,j));
                    Aij_1 = Aij_1.*volume.*dens.*dcFct;
                    Aij_2 = Aij_2.*volume.*dens.*dcFct;
                    Aij_4 = Aij_4.*volume.*dens.*dcFct;
                    Aij_5 = Aij_5.*volume.*dens.*dcFct;
                    sA1(index+1:index+nElem,p) = Aij_1;
                    sA2(index+1:index+nElem,p) = Aij_2;
                    sA4(index+1:index+nElem,p) = Aij_4;
                    sA5(index+1:index+nElem,p) = Aij_5;    
                    end

                    index = index + nElem;
                end
            end
        end
        sA1 = sum(sA1,2);
        sA2 = sum(sA2,2);       
        sA4 = sum(sA4,2);
        sA5 = sum(sA5,2);

        if (Flag3D == 1)
        sA3 = sum(sA3,2);
        sA6 = sum(sA6,2);
        sA7 = sum(sA7,2);
        sA8 = sum(sA8,2);
        sA9 = sum(sA9,2);
        end
        
        % Assemble the matrix        
        A1 = sparse(iif,jjf,sA1,ndof,ndof);
        A2 = sparse(iif,jjf,sA2,ndof,ndof);
        A4 = sparse(iif,jjf,sA4,ndof,ndof);
        A5 = sparse(iif,jjf,sA5,ndof,ndof);
        
        if (Flag3D == 1)
        A3 = sparse(iif,jjf,sA3,ndof,ndof);
        A6 = sparse(iif,jjf,sA6,ndof,ndof);
        A7 = sparse(iif,jjf,sA7,ndof,ndof);
        A8 = sparse(iif,jjf,sA8,ndof,ndof);
        A9 = sparse(iif,jjf,sA9,ndof,ndof);
        end
        
        if (Flag3D == 1)
        stabDC = [A1+A2+A3+A4+A5+A6+A7+A8+A9 ZeroF ZeroF ZeroF;...
                  ZeroF A1+A2+A3+A4+A5+A6+A7+A8+A9 ZeroF ZeroF;...
                  ZeroF ZeroF A1+A2+A3+A4+A5+A6+A7+A8+A9 ZeroF;...
                  ZeroF ZeroF ZeroF ZeroF];
        else
        stabDC = [A1+A2+A4+A5 ZeroF ZeroF;...
                  ZeroF A1+A2+A4+A5 ZeroF;...
                  ZeroF ZeroF ZeroF];
        end

        MatlabFlow = MatlabFlow + (stabDC);
        
        clear sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 
        clear Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
        clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12
        
        RHS = RHS - stabDC *uAlpha(:) ;

        % Get free nodes where solution is to be solved (Change this according to the problem)
        zerotypenodes = find(Sol.type == 0);
	% Free nodes for y-velocity 
        freeNodes1 = unique([bc.top.nodes';bc.bottom.nodes' ]);
        freeNodes1 = setdiff(1:size(crd,1),[freeNodes1; zerotypenodes]);
	% Free nodes for x-velocity
        freeNodes2 = unique([bc.left.nodes';bc.right.nodes' ]);
        freeNodes2 = setdiff(1:size(crd,1),[freeNodes2; zerotypenodes]);
	% Free nodes for pressure
        freeNodes3 = setdiff(1:size(crd,1),[bc.top.nodes'; zerotypenodes]) ;
	% Free nodal matrix with [ x-vel; y-vel; z-vel; pressure ] format (Note here that z-vel is null for 2D problem)
        if (Flag3D == 1)
        freeNodes = [freeNodes2';freeNodes1' + size(crd,1); freeNodes3' + 3*size(crd,1)];
        else
        freeNodes = [freeNodes2';freeNodes1' + size(crd,1); freeNodes3' + 2*size(crd,1)];
        end

        
        result = Sol.uAlpha(:,:,1);
        result = result(:);
        result = [result;Sol.p];
        resultDot = Sol.uDotAlpha(:,:,1);
        resultDot = resultDot(:);
        resultDot = [resultDot;Sol.p];

        % Solve the linear system
        linSol1 = MatlabFlow(freeNodes,freeNodes)\RHS(freeNodes);
        result(freeNodes) = result(freeNodes) + linSol1;
        resultDot(freeNodes) = resultDot(freeNodes) + (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*linSol1 ;
        if (Flag3D == 1)
        Sol.uAlpha(:,:,1) = reshape(result(1:3*ndof),[],3);
        Sol.uDotAlpha(:,:,1) = reshape(resultDot(1:3*ndof),[],3);
        Sol.p = result((3*ndof+1):(4*ndof));
        else
        Sol.uAlpha(:,:,1) = reshape(result(1:2*ndof),[],2);
        Sol.uDotAlpha(:,:,1) = reshape(resultDot(1:2*ndof),[],2);
        Sol.p = result((2*ndof+1):(3*ndof));    
        end
         
        
        % Update the solution
        Sol.u = Sol.uPrev + (1/pmc.alpha)*( Sol.uAlpha - Sol.uPrev ) ;
        Sol.uDot = Sol.uDotPrev + (1/pmc.alphaM)*( Sol.uDotAlpha - Sol.uDotPrev ) ;
        
        normIndicator1 =  norm(linSol1)/norm(result(freeNodes)) ;
        fprintf('NS: %e, ', normIndicator1);
        clear freeNodes1 freeNodes2 freeNodes3
        clear result resultDot
        
        %% Allen-Cahn equation
        % Interpolate for alpha values for velocity
        Sol.uAlpha = Sol.uPrev + pmc.alpha.*(Sol.u - Sol.uPrev) ;

        xxf = zeros(size(cnn));
        yyf = zeros(size(cnn));
        zzf = zeros(size(cnn));
        ux = zeros(size(cnn));
        uy = zeros(size(cnn));
        if (Flag3D == 1)
        uz = zeros(size(cnn));
        end

        for i=1:nen
            xxf(:,i) = crd(cnn(:,i),1);
            yyf(:,i) = crd(cnn(:,i),2);
            if (Flag3D == 1)
            zzf(:,i) = crd(cnn(:,i),3);
            end
            ux(:,i) =  Sol.uAlpha(cnn(:,i),1,1) ;
            uy(:,i) =  Sol.uAlpha(cnn(:,i),2,1) ;
            if (Flag3D == 1)
            uz(:,i) =  Sol.uAlpha(cnn(:,i),3,1) ;
            end
            phi(:,i) = Sol.phiAlpha(cnn(:,i),1) ;
            phiPrev(:,i) = Sol.phiPrev(cnn(:,i),1) ;
            phiDot(:,i) = Sol.phiDotAlpha(cnn(:,i),1) ;
        end
        
        % Lagrange multiplier term and mass
        sB1 = [];
        sB2 = [];
        sB3 = [];
        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
            end
            
            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            localPhiPrev = sum(repmat(N(p,:),nElem,1).*phiPrev,2);
            
            timeA = pmc.alpha ;
            F1 = 0.25*( (1./timeA^3)*localPhi.^3 - (3./timeA^3 - 4./timeA^2)*((localPhi).^2).*localPhi + ...
                      +(3./timeA^3 - 8./timeA^2 + 6./timeA)*localPhi.*localPhi.^2 - (2./timeA).*localPhi + ...
                      +(-1./timeA^3 + 4./timeA^2 - 6./timeA + 4.).*localPhiPrev.^3 + (2./timeA - 4.).*localPhiPrev ) ;  
            F2 = 0.5*( (1./3.)*( (1./timeA^2).*(localPhi).^2 + (-2./timeA^2 + 3./timeA).*localPhi.*localPhiPrev + ...
                      + (1./timeA^2 - 3./timeA + 3.).*(localPhiPrev).^2 ) - 1 ) ;
            Bij1 = gW(p)*F1.*volume;
            Bij2 = gW(p)*F2.*volume;
            Bij3 = gW(p)*localPhiPrev.*volume ;
            sB1(1:nElem,p) = Bij1 ;
            sB2(1:nElem,p) = Bij2 ;
            sB3(1:nElem,p) = Bij3 ;
        end
        sB1 = sum(sum(sB1,2));
        sB2 = sum(sum(sB2,2));
        sB3 = sum(sum(sB3,2));
        lagPar = sB1/(sB2) ;
        if ( (abs(lagPar) < eps) || sB2==0.0 )
            lagPar = 0.0 ;
        end
        MassPhi = sB3 ;
        clear sB1 sB2 sB3 Bij1 Bij2 Bij3
        
        % Galerkin terms for Allen-Cahn equation
        sA = [] ;
        sA1 = zeros(nen^2*nElem,nQuad); 
        sA2 = zeros(nen^2*nElem,nQuad);
        sA3 = zeros(nen^2*nElem,nQuad);

        sA4 = zeros(nen^2*nElem,nQuad);
        sA5 = zeros(nen^2*nElem,nQuad);
        sA6 = zeros(nen^2*nElem,nQuad);
        
        sA7 = zeros(nen^2*nElem,nQuad); 
        sA8 = zeros(nen^2*nElem,nQuad);
        sA9 = zeros(nen^2*nElem,nQuad);

        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
            end
            
            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            localPhiPrev = sum(repmat(N(p,:),nElem,1).*phiPrev,2);
            
            if (Flag3D ~= 1)
                DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
                DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);   
            else
                DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
                      (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
                      (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
                DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
                      (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
                      (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
                DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
                      (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
                      (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);
                locUZ  = sum(repmat(N(p,:),nElem,1).*uz,2);
            end
            
            index = 0;
            for i = 1:nen
                for j = 1:nen
                    % Galerkin inertia term
                    Mij = gW(p)*(N(p,i)*N(p,j));
                    Mij = Mij.*volume;
                    sA(index+1:index+nElem,p) = Mij;
                    
                    % Galerkin convection term
                    Aij_1 = gW(p)*(locUX.*DNDx(:,j)*N(p,i));
                    Aij_2 = gW(p)*(locUY.*DNDy(:,j)*N(p,i));
                    Aij_1 = Aij_1.*volume;
                    Aij_2 = Aij_2.*volume;
                    sA1(index+1:index+nElem,p) = Aij_1;
                    sA2(index+1:index+nElem,p) = Aij_2;
                    
                    % Galerkin diffusion term
                    Aij_4 = gW(p)*(DNDx(:,i).*DNDx(:,j));
                    Aij_5 = gW(p)*(DNDy(:,i).*DNDy(:,j));
                    Aij_4 = Aij_4.*volume.*solver.epsilon^2;
                    Aij_5 = Aij_5.*volume.*solver.epsilon^2;
                    sA4(index+1:index+nElem,p) = Aij_4;
                    sA5(index+1:index+nElem,p) = Aij_5;
                    
                    % Galerkin reaction term 
                    timeA = pmc.alpha ;
                    reacJac = 0.25*( (3/timeA^3)*localPhi.^2 - 2*(3/timeA^3 - 4/timeA^2)*(localPhi).*localPhiPrev + ...
                                    +(3/timeA^3 - 8/timeA^2 + 6/timeA)*localPhiPrev.^2 - 2/timeA ) ...
                              - lagPar*0.5*( (2/(3*timeA^2)).*localPhi + (1/3)*(-2/timeA^2 + 3/timeA).*localPhiPrev );                    
                    Aij_7 = gW(p)*(reacJac).*(N(p,i)*N(p,j));
                    Aij_7 = Aij_7.*volume;
                    sA7(index+1:index+nElem,p) = Aij_7;
                    
                    % Galerkin source term matrices 
                    reac = 0.25*( (1/timeA^3)*localPhi.^2 - (3/timeA^3 - 4/timeA^2)*(localPhi).*localPhiPrev + ...
                            +(3/timeA^3 - 8/timeA^2 + 6/timeA)*localPhiPrev.^2 - 2/timeA ) ...
                            - lagPar*0.5*( (1/(3*timeA^2)).*localPhi + (1/3)*(-2/timeA^2 + 3/timeA).*localPhiPrev );                   
                    Aij_8 = gW(p)*(reac).*(N(p,i)*N(p,j));
                    Aij_8 = Aij_8.*volume;
                    sA8(index+1:index+nElem,p) = Aij_8;
                    
                    source = -0.25*( (-1/timeA^3 + 4/timeA^2 - 6/timeA + 4).*localPhiPrev.^3 + (2/timeA - 4).*localPhiPrev) + ...
                            + lagPar*0.5*( (1/3)*(1/timeA^2 - 3/timeA + 3).*(localPhiPrev.^2) - 1 ) ;                   
                    Aij_9 = gW(p)*(source).*(N(p,i)*N(p,j));
                    Aij_9 = Aij_9.*volume;
                    sA9(index+1:index+nElem,p) = Aij_9;
 
                    if (Flag3D == 1)
                    Aij_3 = gW(p)*(locUZ.*DNDz(:,j)*N(p,i));
                    Aij_3 = Aij_3.*volume;
                    sA3(index+1:index+nElem,p) = Aij_3;
                    Aij_6 = gW(p)*(DNDz(:,i).*DNDz(:,j));
                    Aij_6 = Aij_6.*volume.*solver.epsilon^2;
                    sA6(index+1:index+nElem,p) = Aij_6;
                    end

                    index = index + nElem;
                end
            end
        end
        sA = sum(sA,2);
        sA1 = sum(sA1,2);
        sA2 = sum(sA2,2);
        sA4 = sum(sA4,2);
        sA5 = sum(sA5,2);
        sA7 = sum(sA7,2);
        sA8 = sum(sA8,2);
        sA9 = sum(sA9,2);

        if (Flag3D == 1)
        sA3 = sum(sA3,2);
        sA6 = sum(sA6,2);
        end
        
        % Assemble the matrix
        Mf = sparse(iif,jjf,sA,ndof,ndof);
        
        A1 = sparse(iif,jjf,sA1,ndof,ndof);
        A2 = sparse(iif,jjf,sA2,ndof,ndof);
        A4 = sparse(iif,jjf,sA4,ndof,ndof);
        A5 = sparse(iif,jjf,sA5,ndof,ndof);
        A7 = sparse(iif,jjf,sA7,ndof,ndof);
        A8 = sparse(iif,jjf,sA8,ndof,ndof);
        A9 = sparse(iif,jjf,sA9,ndof,ndof);
        
        if (Flag3D == 1)
        A3 = sparse(iif,jjf,sA3,ndof,ndof);
        A6 = sparse(iif,jjf,sA6,ndof,ndof);
        end
        
        if (Flag3D == 1)
        Mf = [Mf];
        Conv = [A1+A2+A3];
        Kf = [A4+A5+A6];
        Rf = [A7];
        else
        Mf = [Mf];
        Conv = [A1+A2];
        Kf = [A4+A5];
        Rf = [A7];
        end
        MAC = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Mf;

        MatlabAC = MAC + Conv + Kf + Rf;

        phiAlpha = Sol.phiAlpha(:,1) ;
        phiDotAlpha = Sol.phiDotAlpha(:,1);
        
        RHS_AC = -(Mf * phiDotAlpha(:)) - (Conv + Kf) * phiAlpha(:) ;
        RHS_AC = RHS_AC - (A8) * phiAlpha(:) + (A9) * ones(ndof,1) ;
        
        clear sA Mij sA1 sA2 sA3 Aij_1 Aij_2 Aij_3 sA4 sA5 sA6 sA7 sA8 sA9
        clear sA10 sA11 sA12
        clear Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
        clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12
        
        % Petrov-Galerkin stabilization terms for Allen-Cahn equation
        % (tau)(u.gradN).dphi/dt + 
        % (tau)(u.gradN).(u.grad phi) 
        % (tau)(u.gradN).(F'(phi)) + 
        %-(tau)(u.gradN).beta(t).sqrt(F(phi))
        %
        sA1 = zeros(nen^2*nElem,nQuad); 
        sA2 = zeros(nen^2*nElem,nQuad);
        sA3 = zeros(nen^2*nElem,nQuad);
        sA4 = zeros(nen^2*nElem,nQuad);
        sA5 = zeros(nen^2*nElem,nQuad);
        sA6 = zeros(nen^2*nElem,nQuad); 
        sA7 = zeros(nen^2*nElem,nQuad); 
        sA8 = zeros(nen^2*nElem,nQuad);
        sA9 = zeros(nen^2*nElem,nQuad);
        
        sA10 = zeros(nen^2*nElem,nQuad);
        sA11 = zeros(nen^2*nElem,nQuad);
        sA12 = zeros(nen^2*nElem,nQuad);
        
        sA13 = zeros(nen^2*nElem,nQuad);
        sA14 = zeros(nen^2*nElem,nQuad);
        sA15 = zeros(nen^2*nElem,nQuad);
        sA16 = zeros(nen^2*nElem,nQuad); 
        sA17 = zeros(nen^2*nElem,nQuad); 
        sA18 = zeros(nen^2*nElem,nQuad);
        sA19 = zeros(nen^2*nElem,nQuad);
        sA20 = zeros(nen^2*nElem,nQuad);
        sA21 = zeros(nen^2*nElem,nQuad);
        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
                xsInv(:,4) = J(:,1).*J(:,4) - J(:,2).*J(:,3) ;
                xsInv(:,4) = 1./xsInv(:,4) ;
                xsInv(:,1) = xsInv(:,4).* J(:,4) ;
                xsInv(:,2) = xsInv(:,4).* (-J(:,2)) ;
                xsInv(:,3) = xsInv(:,4).* (-J(:,3)) ;
                xsInv(:,4) = xsInv(:,4).* J(:,1) ;
                
                if (strcmp(elemType,'3Tri'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                elseif (strcmp(elemType,'4Quad'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                end
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                      
                xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
                xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
                xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
                xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
                xsInv(:,9) = 1./xsInv(:,9) ;
                xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
                xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
                xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
                xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
                xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
                xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
                xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
                xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
                xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
                if (strcmp(elemType,'4Tet'))
                    c1 = 1.259921049894873E0 ;
                    c2 = 6.299605249474365D-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * (xsInv(:,4) + xsInv(:,7)) ;
                    a2 = c1 * xsInv(:,4) + c2 * (xsInv(:,1) + xsInv(:,7)) ;
                    a3 = c1 * xsInv(:,7) + c2 * (xsInv(:,1) + xsInv(:,4)) ;
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,2) + c2 * (xsInv(:,5) + xsInv(:,8)) ;
                    a2 = c1 * xsInv(:,5) + c2 * (xsInv(:,2) + xsInv(:,8)) ;
                    a3 = c1 * xsInv(:,8) + c2 * (xsInv(:,2) + xsInv(:,5)) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,3) + c2 * (xsInv(:,6) + xsInv(:,9)) ;
                    a2 = c1 * xsInv(:,6) + c2 * (xsInv(:,3) + xsInv(:,9)) ;
                    a3 = c1 * xsInv(:,9) + c2 * (xsInv(:,3) + xsInv(:,6)) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* a3;
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
                elseif (strcmp(elemType,'6Prism'))
                    c1 = 1.154700538379252E0 ;
                    c2 = 5.773502691896259E-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
                    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
                    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
                    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
                    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
                    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
                elseif (strcmp(elemType,'8Hex'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,4).*xsInv(:,4) + xsInv(:,7).*xsInv(:,7);
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,5).*xsInv(:,5) + xsInv(:,8).*xsInv(:,8);
                    gijDown(:,3) = xsInv(:,3).*xsInv(:,3) + xsInv(:,6).*xsInv(:,6) + xsInv(:,9).*xsInv(:,9);
                    gijDown(:,4) = xsInv(:,1).*xsInv(:,2) + xsInv(:,4).*xsInv(:,5) + xsInv(:,7).*xsInv(:,8);
                    gijDown(:,5) = xsInv(:,2).*xsInv(:,3) + xsInv(:,5).*xsInv(:,6) + xsInv(:,8).*xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1).*xsInv(:,3) + xsInv(:,4).*xsInv(:,6) + xsInv(:,7).*xsInv(:,9);
                end
            end
            
            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            localPhiPrev = sum(repmat(N(p,:),nElem,1).*phiPrev,2);   
               
            if (Flag3D ~= 1)
                DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
                DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2); 
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,3) .* locUY ;
                VG2	= gijDown(:,3) .* locUX + gijDown(:,2) .* locUY ;

                VGV	= VG1 .* locUX + VG2 .* locUY ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + 2 * ( gijDown(:,3).^2 );
            else
                DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
                      (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
                      (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
                DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
                      (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
                      (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
                DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
                      (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
                      (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);
                locUZ  = sum(repmat(N(p,:),nElem,1).*uz,2);
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
                VG2	= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
                VG3	= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

                VGV	= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
                    2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
            end
            timeA = pmc.alpha ;
            reac = 0.25*( (1/timeA^3)*localPhi.^2 - (3/timeA^3 - 4/timeA^2)*(localPhi).*localPhiPrev + ...
                            +(3/timeA^3 - 8/timeA^2 + 6/timeA)*localPhiPrev.^2 - 2/timeA ) ...
                            - lagPar*0.5*( (1/(3*timeA^2)).*localPhi + (1/3)*(-2/timeA^2 + 3/timeA).*localPhiPrev ); 
            tauP = (2/solver.dt)^2 + (solver.epsilon^4).* GG * 9 + abs(VGV) + reac.^2;    
            tauP = tauP.^(0.5);
            tauP = 1./tauP;
            
            index = 0;
            for i = 1:nen
                for j = 1:nen
                    % 1st term
                    Aij_1 = gW(p)*(locUX.*DNDx(:,i)*N(p,j));
                    Aij_2 = gW(p)*(locUY.*DNDy(:,i)*N(p,j));
                    Aij_1 = Aij_1.*volume.*tauP ;
                    Aij_2 = Aij_2.*volume.*tauP ;
                    sA1(index+1:index+nElem,p) = Aij_1;
                    sA2(index+1:index+nElem,p) = Aij_2;
                    
                    % 2nd term
                    Aij_4 = gW(p)*locUX.*locUX.*(DNDx(:,i).*DNDx(:,j));
                    Aij_5 = gW(p)*locUX.*locUY.*(DNDx(:,i).*DNDy(:,j));
                    Aij_7 = gW(p)*locUY.*locUX.*(DNDy(:,i).*DNDx(:,j));
                    Aij_8 = gW(p)*locUY.*locUY.*(DNDy(:,i).*DNDy(:,j)); 
                    Aij_4 = Aij_4.*volume.*tauP ;
                    Aij_5 = Aij_5.*volume.*tauP ;
                    Aij_7 = Aij_7.*volume.*tauP ;
                    Aij_8 = Aij_8.*volume.*tauP ;
                    sA4(index+1:index+nElem,p) = Aij_4;
                    sA5(index+1:index+nElem,p) = Aij_5;
                    sA7(index+1:index+nElem,p) = Aij_7;
                    sA8(index+1:index+nElem,p) = Aij_8;
                    
                    % 3rd term
                    timeA = pmc.alpha ;
                    reacJac = 0.25*( (3/timeA^3)*localPhi.^2 - 2*(3/timeA^3 - 4/timeA^2)*(localPhi).*localPhiPrev + ...
                                    +(3/timeA^3 - 8/timeA^2 + 6/timeA)*localPhiPrev.^2 - 2/timeA ) ...
                              - lagPar*0.5*( (2/(3*timeA^2)).*localPhi + (1/3)*(-2/timeA^2 + 3/timeA).*localPhiPrev );                      
                    Aij_13 = gW(p)*(reacJac).*(locUX.*DNDx(:,i)*N(p,j));
                    Aij_14 = gW(p)*(reacJac).*(locUY.*DNDy(:,i)*N(p,j));
                    Aij_13 = Aij_13.*volume.*tauP;
                    Aij_14 = Aij_14.*volume.*tauP;
                    sA13(index+1:index+nElem,p) = Aij_13;
                    sA14(index+1:index+nElem,p) = Aij_14;
                    
                    % Galerkin source term matrices 
                    reac = 0.25*( (1/timeA^3)*localPhi.^2 - (3/timeA^3 - 4/timeA^2)*(localPhi).*localPhiPrev + ...
                            +(3/timeA^3 - 8/timeA^2 + 6/timeA)*localPhiPrev.^2 - 2/timeA ) ...
                            - lagPar*0.5*( (1/(3*timeA^2)).*localPhi + (1/3)*(-2/timeA^2 + 3/timeA).*localPhiPrev );                  
                    Aij_16 = gW(p)*(reac).*(locUX.*DNDx(:,i)*N(p,j));
                    Aij_17 = gW(p)*(reac).*(locUY.*DNDy(:,i)*N(p,j));
                    Aij_16 = Aij_16.*volume.*tauP;
                    Aij_17 = Aij_17.*volume.*tauP;
                    sA16(index+1:index+nElem,p) = Aij_16;
                    sA17(index+1:index+nElem,p) = Aij_17;
                    
                    source = -0.25*( (-1/timeA^3 + 4/timeA^2 - 6/timeA + 4).*localPhiPrev.^3 + (2/timeA - 4).*localPhiPrev) + ...
                            + lagPar*0.5*( (1/3)*(1/timeA^2 - 3/timeA + 3).*(localPhiPrev.^2) - 1 ) ;                   
                    Aij_19 = gW(p)*(source).*(locUX.*DNDx(:,i)*N(p,j));
                    Aij_20 = gW(p)*(source).*(locUY.*DNDy(:,i)*N(p,j));
                    Aij_19 = Aij_19.*volume.*tauP;
                    Aij_20 = Aij_20.*volume.*tauP;
                    sA19(index+1:index+nElem,p) = Aij_19;
                    sA20(index+1:index+nElem,p) = Aij_20;
                    
                    if (Flag3D == 1)
                    Aij_3 = gW(p)*(locUZ.*DNDz(:,i)*N(p,j));
                    Aij_3 = Aij_3.*volume.*tauP ;
                    sA3(index+1:index+nElem,p) = Aij_3;
                    Aij_6 = gW(p)*locUX.*locUZ.*(DNDx(:,i).*DNDz(:,j));
                    Aij_9 = gW(p)*locUY.*locUZ.*(DNDy(:,i).*DNDz(:,j));
                    Aij_10 = gW(p)*locUZ.*locUX.*(DNDz(:,i).*DNDx(:,j));
                    Aij_11 = gW(p)*locUZ.*locUY.*(DNDz(:,i).*DNDy(:,j));
                    Aij_12 = gW(p)*locUZ.*locUZ.*(DNDz(:,i).*DNDz(:,j));
                    Aij_6 = Aij_6.*volume.*tauP ;
                    Aij_9 = Aij_9.*volume.*tauP ;
                    Aij_10 = Aij_10.*volume.*tauP ;
                    Aij_11 = Aij_11.*volume.*tauP ;
                    Aij_12 = Aij_12.*volume.*tauP ;
                    sA6(index+1:index+nElem,p) = Aij_6;
                    sA9(index+1:index+nElem,p) = Aij_9;
                    sA10(index+1:index+nElem,p) = Aij_10;
                    sA11(index+1:index+nElem,p) = Aij_11;
                    sA12(index+1:index+nElem,p) = Aij_12;
                    Aij_15 = gW(p)*(reacJac).*(locUZ.*DNDz(:,i)*N(p,j));
                    Aij_15 = Aij_15.*volume.*tauP;
                    sA15(index+1:index+nElem,p) = Aij_15;
                    Aij_18 = gW(p)*(reac).*(locUZ.*DNDz(:,i)*N(p,j));
                    Aij_18 = Aij_18.*volume.*tauP;
                    sA18(index+1:index+nElem,p) = Aij_18;
                    Aij_21 = gW(p)*(source).*(locUZ.*DNDz(:,i)*N(p,j));
                    Aij_21 = Aij_21.*volume.*tauP;
                    sA21(index+1:index+nElem,p) = Aij_21;
                    end

                    index = index + nElem;
                end
            end
        end
        sA1 = sum(sA1,2);
        sA2 = sum(sA2,2);       
        sA4 = sum(sA4,2);
        sA5 = sum(sA5,2);
        sA7 = sum(sA7,2);
        sA8 = sum(sA8,2);
        sA13 = sum(sA13,2);
        sA14 = sum(sA14,2);       
        sA16 = sum(sA16,2);
        sA17 = sum(sA17,2);
        sA19 = sum(sA19,2);
        sA20 = sum(sA20,2);

        if (Flag3D == 1)
        sA3 = sum(sA3,2);
        sA6 = sum(sA6,2);
        sA9 = sum(sA9,2);
        sA10 = sum(sA10,2);
        sA11 = sum(sA11,2);
        sA12 = sum(sA12,2);
        sA15 = sum(sA15,2);
        sA18 = sum(sA18,2);
        sA21 = sum(sA21,2);
        end
        
        % Assemble the matrix        
        A1 = sparse(iif,jjf,sA1,ndof,ndof);
        A2 = sparse(iif,jjf,sA2,ndof,ndof);
        A4 = sparse(iif,jjf,sA4,ndof,ndof);
        A5 = sparse(iif,jjf,sA5,ndof,ndof);
        A7 = sparse(iif,jjf,sA7,ndof,ndof);
        A8 = sparse(iif,jjf,sA8,ndof,ndof);
        A13 = sparse(iif,jjf,sA13,ndof,ndof);
        A14 = sparse(iif,jjf,sA14,ndof,ndof);
        A16 = sparse(iif,jjf,sA16,ndof,ndof);
        A17 = sparse(iif,jjf,sA17,ndof,ndof);
        A19 = sparse(iif,jjf,sA19,ndof,ndof);
        A20 = sparse(iif,jjf,sA20,ndof,ndof);
        
        if (Flag3D == 1)
        A3 = sparse(iif,jjf,sA3,ndof,ndof);
        A6 = sparse(iif,jjf,sA6,ndof,ndof);
        A9 = sparse(iif,jjf,sA9,ndof,ndof);
        A10 = sparse(iif,jjf,sA10,ndof,ndof);
        A11 = sparse(iif,jjf,sA11,ndof,ndof);
        A12 = sparse(iif,jjf,sA12,ndof,ndof);
        A15 = sparse(iif,jjf,sA15,ndof,ndof);
        A18 = sparse(iif,jjf,sA18,ndof,ndof);
        A21 = sparse(iif,jjf,sA21,ndof,ndof);
        end
        
        if (Flag3D == 1)
        stabAC1 = [A1+A2+A3];
        stabAC2 = [A4+A5+A6+A7+A8+A9+A10+A11+A12];
        stabAC3 = [A13+A14+A15];
        stabAC4 = [A16+A17+A18];
        stabAC5 = [A19+A20+A21];
        else
        stabAC1 = [A1+A2];
        stabAC2 = [A4+A5+A7+A8];
        stabAC3 = [A13+A14];
        stabAC4 = [A16+A17];
        stabAC5 = [A19+A20];
        end
        stabAC1Mat = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*stabAC1;

        MatlabAC = MatlabAC + stabAC1Mat + stabAC2 + stabAC3;
        
        clear sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 sA13 sA14 sA15
        clear sA16 sA17 sA18 sA19 sA20 sA21
        clear Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
        clear Aij_13 Aij_14 Aij_15 Aij_16 Aij_17 Aij_18 Aij_19 Aij_20 Aij_21
        clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16 A17 A18 A19 A20 A21
        
        RHS_AC = RHS_AC -(stabAC1) * phiDotAlpha(:) - (stabAC2 + stabAC4) * phiAlpha(:)  ;
        RHS_AC = RHS_AC +(stabAC5) * ones(ndof,1) ;
        
        % PPV stabilization terms for Allen-Cahn equation
        % 
        %
        sA1 = zeros(nen^2*nElem,nQuad); 
        sA2 = zeros(nen^2*nElem,nQuad);
        sA3 = zeros(nen^2*nElem,nQuad);
        sA4 = zeros(nen^2*nElem,nQuad);
        sA5 = zeros(nen^2*nElem,nQuad);
        sA6 = zeros(nen^2*nElem,nQuad); 
        sA7 = zeros(nen^2*nElem,nQuad); 
        sA8 = zeros(nen^2*nElem,nQuad);
        sA9 = zeros(nen^2*nElem,nQuad);
        
        sA10 = zeros(nen^2*nElem,nQuad);
        sA11 = zeros(nen^2*nElem,nQuad);
        sA12 = zeros(nen^2*nElem,nQuad);
        
        sA13 = zeros(nen^2*nElem,nQuad);
        sA14 = zeros(nen^2*nElem,nQuad);
        sA15 = zeros(nen^2*nElem,nQuad);
        sA16 = zeros(nen^2*nElem,nQuad); 
        sA17 = zeros(nen^2*nElem,nQuad); 
        sA18 = zeros(nen^2*nElem,nQuad);
        sA19 = zeros(nen^2*nElem,nQuad);
        sA20 = zeros(nen^2*nElem,nQuad);
        sA21 = zeros(nen^2*nElem,nQuad);
        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
                xsInv(:,4) = J(:,1).*J(:,4) - J(:,2).*J(:,3) ;
                xsInv(:,4) = 1./xsInv(:,4) ;
                xsInv(:,1) = xsInv(:,4).* J(:,4) ;
                xsInv(:,2) = xsInv(:,4).* (-J(:,2)) ;
                xsInv(:,3) = xsInv(:,4).* (-J(:,3)) ;
                xsInv(:,4) = xsInv(:,4).* J(:,1) ;
                
                if (strcmp(elemType,'3Tri'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                elseif (strcmp(elemType,'4Quad'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                end
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                      
                xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
                xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
                xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
                xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
                xsInv(:,9) = 1./xsInv(:,9) ;
                xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
                xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
                xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
                xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
                xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
                xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
                xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
                xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
                xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
                if (strcmp(elemType,'4Tet'))
                    c1 = 1.259921049894873E0 ;
                    c2 = 6.299605249474365D-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * (xsInv(:,4) + xsInv(:,7)) ;
                    a2 = c1 * xsInv(:,4) + c2 * (xsInv(:,1) + xsInv(:,7)) ;
                    a3 = c1 * xsInv(:,7) + c2 * (xsInv(:,1) + xsInv(:,4)) ;
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,2) + c2 * (xsInv(:,5) + xsInv(:,8)) ;
                    a2 = c1 * xsInv(:,5) + c2 * (xsInv(:,2) + xsInv(:,8)) ;
                    a3 = c1 * xsInv(:,8) + c2 * (xsInv(:,2) + xsInv(:,5)) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,3) + c2 * (xsInv(:,6) + xsInv(:,9)) ;
                    a2 = c1 * xsInv(:,6) + c2 * (xsInv(:,3) + xsInv(:,9)) ;
                    a3 = c1 * xsInv(:,9) + c2 * (xsInv(:,3) + xsInv(:,6)) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* a3;
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
                elseif (strcmp(elemType,'6Prism'))
                    c1 = 1.154700538379252E0 ;
                    c2 = 5.773502691896259E-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
                    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
                    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
                    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
                    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
                    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
                elseif (strcmp(elemType,'8Hex'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,4).*xsInv(:,4) + xsInv(:,7).*xsInv(:,7);
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,5).*xsInv(:,5) + xsInv(:,8).*xsInv(:,8);
                    gijDown(:,3) = xsInv(:,3).*xsInv(:,3) + xsInv(:,6).*xsInv(:,6) + xsInv(:,9).*xsInv(:,9);
                    gijDown(:,4) = xsInv(:,1).*xsInv(:,2) + xsInv(:,4).*xsInv(:,5) + xsInv(:,7).*xsInv(:,8);
                    gijDown(:,5) = xsInv(:,2).*xsInv(:,3) + xsInv(:,5).*xsInv(:,6) + xsInv(:,8).*xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1).*xsInv(:,3) + xsInv(:,4).*xsInv(:,6) + xsInv(:,7).*xsInv(:,9);
                end
            end
            
            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            localPhiDot  = sum(repmat(N(p,:),nElem,1).*phiDot,2);
            localPhiPrev = sum(repmat(N(p,:),nElem,1).*phiPrev,2);   
               
            if (Flag3D ~= 1)
                DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
                DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2); 
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,3) .* locUY ;
                VG2	= gijDown(:,3) .* locUX + gijDown(:,2) .* locUY ;

                VGV	= VG1 .* locUX + VG2 .* locUY ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + 2 * ( gijDown(:,3).^2 );
            else
                DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
                      (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
                      (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
                DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
                      (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
                      (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
                DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
                      (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
                      (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);
                locUZ  = sum(repmat(N(p,:),nElem,1).*uz,2);
                
                VG1 = gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
                VG2	= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
                VG3	= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

                VGV	= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

                GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
                    2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
            end
            timeA = pmc.alpha ;
            reac = 0.25*( (1/timeA^3)*localPhi.^2 - (3/timeA^3 - 4/timeA^2)*(localPhi).*localPhiPrev + ...
                            +(3/timeA^3 - 8/timeA^2 + 6/timeA)*localPhiPrev.^2 - 2/timeA ) ...
                            - lagPar*0.5*( (1/(3*timeA^2)).*localPhi + (1/3)*(-2/timeA^2 + 3/timeA).*localPhiPrev ); 
            source = -0.25*( (-1/timeA^3 + 4/timeA^2 - 6/timeA + 4).*localPhiPrev.^3 + (2/timeA - 4).*localPhiPrev) + ...
                            + lagPar*0.5*( (1/3)*(1/timeA^2 - 3/timeA + 3).*(localPhiPrev.^2) - 1 ) ;   
            tauP = (2/solver.dt)^2 + (solver.epsilon^4).* GG * 9 + abs(VGV) + reac.^2;    
            tauP = tauP.^(0.5);
            tauP = 1./tauP;
            
            if (Flag3D==1)
            charLen = 2.0 * sqrt(locUX.^2 + locUY.^2 + locUZ.^2)./sqrt(VGV) ;
            else
            charLen = 2.0 * sqrt(locUX.^2 + locUY.^2)./sqrt(VGV) ;
            end
            
            localGradPhiX  = sum(DNDx.*phi,2);
            localGradPhiY  = sum(DNDy.*phi,2);
            if (Flag3D==1)
            localGradPhiZ  = sum(DNDz.*phi,2);
            Resphi = localPhiDot + locUX.*localGradPhiX + locUY.*localGradPhiY + ...
                     locUZ.*localGradPhiZ + reac.* localPhi - source ;
            absGradPhi = sqrt(localGradPhiX.^2 + localGradPhiY.^2 + localGradPhiZ.^2) ;
            modVel = sqrt(locUX.^2 + locUY.^2 + locUZ.^2) ;
            else
            Resphi = localPhiDot + locUX.*localGradPhiX + locUY.*localGradPhiY + ...
                     reac.* localPhi - source ;
            absGradPhi = sqrt(localGradPhiX.^2 + localGradPhiY.^2) ;
            modVel = sqrt(locUX.^2 + locUY.^2) ;
            end
            
            RPhibyGradPhi = abs(Resphi)./absGradPhi ;
            RPhibyGradPhi(absGradPhi==0) = 0 ;
            
            tildeReac = 1/(pmc.alpha*solver.dt) + reac ;
            
            chi = 2.0./ (abs(tildeReac).*charLen + 2.*modVel) ;
            ksAdd = max( abs(modVel - tauP.*modVel.*tildeReac).*charLen/2.0 ...
                    - (solver.epsilon^2 + tauP.*(modVel.^2)) + (tildeReac).*charLen.^2/6, 0.0) ;
            kcAdd = max( abs(modVel ).*charLen/2.0 ...
                    - (solver.epsilon^2) + (tildeReac).*charLen.^2/6, 0.0) ;    
            sdk = chi.*RPhibyGradPhi.*ksAdd ;
            cdk = chi.*RPhibyGradPhi.*kcAdd ;
            if (modVel ~= 0)
            uTens1 = (locUX.^2)./ modVel.^2 ;
            uTens2 = (locUY.^2)./ modVel.^2 ;
            uTens4 = (locUX.*locUY)./ modVel.^2 ;
            if (Flag3D == 1)
            uTens3 = (locUZ.^2)./ modVel.^2 ;
            uTens5 = (locUY.*locUZ)./ modVel.^2 ;
            uTens6 = (locUX.*locUZ)./ modVel.^2 ;
            end
            else
            uTens1 = 0.0.*locUX ;
            uTens2 = 0.0.*locUX ;
            uTens4 = 0.0.*locUX ;
            if (Flag3D == 1)
            uTens3 = 0.0.*locUX ;
            uTens5 = 0.0.*locUX ;
            uTens6 = 0.0.*locUX ;
            end
            end
            
            index = 0;
            for i = 1:nen
                for j = 1:nen
                    % Streamline terms
                    Aij_1 = gW(p)*sdk.*uTens1.*(DNDx(:,i).*DNDx(:,j));
                    Aij_2 = gW(p)*sdk.*uTens4.*(DNDx(:,i).*DNDy(:,j));
                    Aij_4 = gW(p)*sdk.*uTens4.*(DNDy(:,i).*DNDx(:,j));
                    Aij_5 = gW(p)*sdk.*uTens2.*(DNDy(:,i).*DNDy(:,j)); 
                    Aij_1 = Aij_1.*volume ;
                    Aij_2 = Aij_2.*volume ;                   
                    Aij_4 = Aij_4.*volume ;
                    Aij_5 = Aij_5.*volume ;                   
                    sA1(index+1:index+nElem,p) = Aij_1;
                    sA2(index+1:index+nElem,p) = Aij_2;                 
                    sA4(index+1:index+nElem,p) = Aij_4;
                    sA5(index+1:index+nElem,p) = Aij_5;

                    % Crosswind terms
                    Aij_10 = gW(p)*cdk.*(1 - uTens1).*(DNDx(:,i).*DNDx(:,j));
                    Aij_11 = gW(p)*cdk.*(-uTens4).*(DNDx(:,i).*DNDy(:,j));  
                    Aij_13 = gW(p)*cdk.*(-uTens4).*(DNDy(:,i).*DNDx(:,j));
                    Aij_14 = gW(p)*cdk.*(1 - uTens2).*(DNDy(:,i).*DNDy(:,j));      
                    Aij_10 = Aij_10.*volume ;
                    Aij_11 = Aij_11.*volume ;   
                    Aij_13 = Aij_13.*volume ;
                    Aij_14 = Aij_14.*volume ;  
                    sA10(index+1:index+nElem,p) = Aij_10;
                    sA11(index+1:index+nElem,p) = Aij_11;
                    sA13(index+1:index+nElem,p) = Aij_13;
                    sA14(index+1:index+nElem,p) = Aij_14;
 
                    if (Flag3D == 1)
                    Aij_3 = gW(p)*sdk.*uTens6.*(DNDx(:,i).*DNDz(:,j));
                    Aij_6 = gW(p)*sdk.*uTens5.*(DNDy(:,i).*DNDz(:,j));
                    Aij_7 = gW(p)*sdk.*uTens6.*(DNDz(:,i).*DNDx(:,j));
                    Aij_8 = gW(p)*sdk.*uTens5.*(DNDz(:,i).*DNDy(:,j)); 
                    Aij_9 = gW(p)*sdk.*uTens3.*(DNDz(:,i).*DNDz(:,j));
                    Aij_3 = Aij_3.*volume ;
                    Aij_6 = Aij_6.*volume ;
                    Aij_7 = Aij_7.*volume ;
                    Aij_8 = Aij_8.*volume ;
                    Aij_9 = Aij_9.*volume ;
                    sA3(index+1:index+nElem,p) = Aij_3;
                    sA6(index+1:index+nElem,p) = Aij_6;
                    sA7(index+1:index+nElem,p) = Aij_7;
                    sA8(index+1:index+nElem,p) = Aij_8;
                    sA9(index+1:index+nElem,p) = Aij_9;
                    Aij_12 = gW(p)*cdk.*(-uTens6).*(DNDx(:,i).*DNDz(:,j));
                    Aij_15 = gW(p)*cdk.*(-uTens5).*(DNDy(:,i).*DNDz(:,j));
                    Aij_16 = gW(p)*cdk.*(-uTens6).*(DNDz(:,i).*DNDx(:,j));
                    Aij_17 = gW(p)*cdk.*(-uTens5).*(DNDz(:,i).*DNDy(:,j)); 
                    Aij_18 = gW(p)*cdk.*(1 - uTens3).*(DNDz(:,i).*DNDz(:,j));
                    Aij_12 = Aij_12.*volume ;
                    Aij_15 = Aij_15.*volume ;
                    Aij_16 = Aij_16.*volume ;
                    Aij_17 = Aij_17.*volume ;
                    Aij_18 = Aij_18.*volume ;
                    sA12(index+1:index+nElem,p) = Aij_12;
                    sA15(index+1:index+nElem,p) = Aij_15;
                    sA16(index+1:index+nElem,p) = Aij_16;
                    sA17(index+1:index+nElem,p) = Aij_17;
                    sA18(index+1:index+nElem,p) = Aij_18;
                    end

                    index = index + nElem;
                end
            end
        end
        sA1 = sum(sA1,2);
        sA2 = sum(sA2,2);       
        sA4 = sum(sA4,2);
        sA5 = sum(sA5,2);
        sA10 = sum(sA10,2);
        sA11 = sum(sA11,2);
        sA13 = sum(sA13,2);
        sA14 = sum(sA14,2);

        if (Flag3D == 1)
        sA3 = sum(sA3,2);
        sA6 = sum(sA6,2);
        sA7 = sum(sA7,2);
        sA8 = sum(sA8,2);
        sA9 = sum(sA9,2);
        sA12 = sum(sA12,2);
        sA15 = sum(sA15,2);
        sA16 = sum(sA16,2);
        sA17 = sum(sA17,2);
        sA18 = sum(sA18,2);
        end
        
        % Assemble the matrix        
        A1 = sparse(iif,jjf,sA1,ndof,ndof);
        A2 = sparse(iif,jjf,sA2,ndof,ndof);
        A4 = sparse(iif,jjf,sA4,ndof,ndof);
        A5 = sparse(iif,jjf,sA5,ndof,ndof);
        A10 = sparse(iif,jjf,sA10,ndof,ndof);
        A11 = sparse(iif,jjf,sA11,ndof,ndof);
        A13 = sparse(iif,jjf,sA13,ndof,ndof);
        A14 = sparse(iif,jjf,sA14,ndof,ndof);
        
        if (Flag3D == 1)
        A3 = sparse(iif,jjf,sA3,ndof,ndof);
        A6 = sparse(iif,jjf,sA6,ndof,ndof);
        A7 = sparse(iif,jjf,sA7,ndof,ndof);
        A8 = sparse(iif,jjf,sA8,ndof,ndof);
        A9 = sparse(iif,jjf,sA9,ndof,ndof);
        A12 = sparse(iif,jjf,sA12,ndof,ndof);
        A15 = sparse(iif,jjf,sA15,ndof,ndof);
        A16 = sparse(iif,jjf,sA16,ndof,ndof);
        A17 = sparse(iif,jjf,sA17,ndof,ndof);
        A18 = sparse(iif,jjf,sA18,ndof,ndof);
        end
        
        if (Flag3D == 1)
        stabAC_PPV = [A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+...
                      A16+A17+A18];
        else
        stabAC_PPV = [A1+A2+A4+A5+A10+A11+A13+A14];
        end
        

        MatlabAC = MatlabAC + stabAC_PPV;
        RHS_AC = RHS_AC - (stabAC_PPV) * phiAlpha(:) ;
        clear xsInv gijDown phiPrev gijUp

	% Free nodes for Allen-Cahn equation (All the nodes are to be solved since Neumann condition is satisfied on all the boundaries)
        zerotypenodes = find(Sol.type == 0);
        freeNodes = setdiff(1:size(crd,1),zerotypenodes);
        freeNodes = [freeNodes'];

        result = Sol.phiAlpha(:,1);
        result = result(:);
        resultDot = Sol.phiDotAlpha(:,1);
        resultDot = resultDot(:);

	% Solve the linear system
        linSol2 = MatlabAC(freeNodes,freeNodes)\RHS_AC(freeNodes);
        result(freeNodes) = result(freeNodes) + linSol2;
        resultDot(freeNodes) = resultDot(freeNodes) + (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*linSol2 ;

        Sol.phiAlpha(:,:,1) = reshape(result(1:ndof),[],1);
        Sol.phiDotAlpha(:,:,1) = reshape(resultDot(1:ndof),[],1);
        
        
        % Update the solution
        Sol.phi = Sol.phiPrev + (1/pmc.alpha)*( Sol.phiAlpha - Sol.phiPrev ) ;
        Sol.phiDot = Sol.phiDotPrev + (1/pmc.alphaM)*( Sol.phiDotAlpha - Sol.phiDotPrev ) ;
        
        normIndicator2 = norm(linSol2)/norm(result(freeNodes)) ;
        fprintf('AC: %e, Mass: %e, ', normIndicator2,MassPhi);
        clear result resultDot

        clear linSol1 linSol2
        clear Sol.uNew Sol.uDotNew Sol.uNewPrev Sol.uDotNewPrev Sol.phiAlpha phi
        clear uxDot uyDot uzDot pres locGradU ResU gGradV
        
        %% Adaptive algorithm
        % Indicator for adaptivity
        
        Sol.uAlpha = Sol.uPrev + pmc.alpha.*(Sol.u - Sol.uPrev) ;
        Sol.uDotAlpha = Sol.uDotPrev + pmc.alphaM.*(Sol.uDot - Sol.uDotPrev) ;
        Sol.phiAlpha = Sol.phiPrev + pmc.alpha.*(Sol.phi - Sol.phiPrev) ;
        Sol.phiDotAlpha = Sol.phiDotPrev + pmc.alphaM.*(Sol.phiDot - Sol.phiDotPrev) ;
        for i=1:nen
            xxf(:,i) = crd(cnn(:,i),1);
            yyf(:,i) = crd(cnn(:,i),2);
            if (Flag3D == 1)
            zzf(:,i) = crd(cnn(:,i),3);
            end
            ux(:,i) =  Sol.uAlpha(cnn(:,i),1,1) ;
            uy(:,i) =  Sol.uAlpha(cnn(:,i),2,1) ;
            phi(:,i) = Sol.phiAlpha(cnn(:,i),1) ;
            phiDot(:,i) = Sol.phiDotAlpha(cnn(:,i),1) ;
            phiPrev(:,i) = Sol.phiPrev(cnn(:,i),1) ;
        end
        clear lagPar ;
        
        % Lagrange multiplier term
        sB1 = [];
        sB2 = [];
        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
            end
            
            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            localPhiPrev = sum(repmat(N(p,:),nElem,1).*phiPrev,2);
            
            timeA = pmc.alpha ;
            F1 = 0.25*( (1./timeA^3)*localPhi.^3 - (3./timeA^3 - 4./timeA^2)*((localPhi).^2).*localPhi + ...
                      +(3./timeA^3 - 8./timeA^2 + 6./timeA)*localPhi.*localPhi.^2 - (2./timeA).*localPhi + ...
                      +(-1./timeA^3 + 4./timeA^2 - 6./timeA + 4.).*localPhiPrev.^3 + (2./timeA - 4.).*localPhiPrev ) ;  
            F2 = 0.5*( (1./3.)*( (1./timeA^2).*(localPhi).^2 + (-2./timeA^2 + 3./timeA).*localPhi.*localPhiPrev + ...
                      + (1./timeA^2 - 3./timeA + 3.).*(localPhiPrev).^2 ) - 1 ) ;
            Bij1 = gW(p)*F1.*volume;
            Bij2 = gW(p)*F2.*volume;
            sB1(1:nElem,p) = Bij1;
            sB2(1:nElem,p) = Bij2;
        end
        sB1 = sum(sum(sB1,2));
        sB2 = sum(sum(sB2,2));
        lagPar = sB1/(sB2) ;
        if ( (abs(lagPar) < eps) || sB2==0.0 )
            lagPar = 0.0 ;
        end
        clear sB1 sB2 Bij1 Bij2 volume localPhi localPhiPrev
        
        
        for p = 1:nQuad  
            if (Flag3D ~= 1)
               J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
                    yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) ); 
                xsInv(:,4) = J(:,1).*J(:,4) - J(:,2).*J(:,3) ;
                xsInv(:,4) = 1./xsInv(:,4) ;
                xsInv(:,1) = xsInv(:,4).* J(:,4) ;
                xsInv(:,2) = xsInv(:,4).* (-J(:,2)) ;
                xsInv(:,3) = xsInv(:,4).* (-J(:,3)) ;
                xsInv(:,4) = xsInv(:,4).* J(:,1) ;
                
                if (strcmp(elemType,'3Tri'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                elseif (strcmp(elemType,'4Quad'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
                    gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
                end
            else
                J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
                     yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
                     zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
                if size(J,2)==1
                    J = J';
                end
                volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                          J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                          J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                      
                xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
                xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
                xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
                xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
                xsInv(:,9) = 1./xsInv(:,9) ;
                xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
                xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
                xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
                xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
                xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
                xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
                xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
                xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
                xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
                if (strcmp(elemType,'4Tet'))
                    c1 = 1.259921049894873E0 ;
                    c2 = 6.299605249474365D-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * (xsInv(:,4) + xsInv(:,7)) ;
                    a2 = c1 * xsInv(:,4) + c2 * (xsInv(:,1) + xsInv(:,7)) ;
                    a3 = c1 * xsInv(:,7) + c2 * (xsInv(:,1) + xsInv(:,4)) ;
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,2) + c2 * (xsInv(:,5) + xsInv(:,8)) ;
                    a2 = c1 * xsInv(:,5) + c2 * (xsInv(:,2) + xsInv(:,8)) ;
                    a3 = c1 * xsInv(:,8) + c2 * (xsInv(:,2) + xsInv(:,5)) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
    
                    a1 = c1 * xsInv(:,3) + c2 * (xsInv(:,6) + xsInv(:,9)) ;
                    a2 = c1 * xsInv(:,6) + c2 * (xsInv(:,3) + xsInv(:,9)) ;
                    a3 = c1 * xsInv(:,9) + c2 * (xsInv(:,3) + xsInv(:,6)) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* a3;
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
                elseif (strcmp(elemType,'6Prism'))
                    c1 = 1.154700538379252E0 ;
                    c2 = 5.773502691896259E-01;
    
                    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
                    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
                    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
                    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
                    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
                    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
                    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
                    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
                    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
                    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
                    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
                elseif (strcmp(elemType,'8Hex'))
                    gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,4).*xsInv(:,4) + xsInv(:,7).*xsInv(:,7);
                    gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,5).*xsInv(:,5) + xsInv(:,8).*xsInv(:,8);
                    gijDown(:,3) = xsInv(:,3).*xsInv(:,3) + xsInv(:,6).*xsInv(:,6) + xsInv(:,9).*xsInv(:,9);
                    gijDown(:,4) = xsInv(:,1).*xsInv(:,2) + xsInv(:,4).*xsInv(:,5) + xsInv(:,7).*xsInv(:,8);
                    gijDown(:,5) = xsInv(:,2).*xsInv(:,3) + xsInv(:,5).*xsInv(:,6) + xsInv(:,8).*xsInv(:,9);
                    gijDown(:,6) = xsInv(:,1).*xsInv(:,3) + xsInv(:,4).*xsInv(:,6) + xsInv(:,7).*xsInv(:,9);
                end
            end

            volume = abs(volume);
            localPhi  = sum(repmat(N(p,:),nElem,1).*phi,2);
            localPhiDot  = sum(repmat(N(p,:),nElem,1).*phiDot,2); 
            localPhiPrev = sum(repmat(N(p,:),nElem,1).*phiPrev,2); 
               
            if (Flag3D ~= 1)
                DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
                DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2); 

            else
                DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
                      (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
                      (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
                DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
                      (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
                      (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
                DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
                      (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
                      (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);
    
                locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
                locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);
                locUZ  = sum(repmat(N(p,:),nElem,1).*uz,2);
    
            end
            timeA = pmc.alpha ;
            reacJac = 0.25*( (3/timeA^3)*localPhi.^2 - 2*(3/timeA^3 - 4/timeA^2)*(localPhi).*localPhiPrev + ...
                                    +(3/timeA^3 - 8/timeA^2 + 6/timeA)*localPhiPrev.^2 - 2/timeA ) ...
                              - lagPar*0.5*( (2/(3*timeA^2)).*localPhi + (1/3)*(-2/timeA^2 + 3/timeA).*localPhiPrev );  
            reac = 0.25*( (1/timeA^3)*localPhi.^2 - (3/timeA^3 - 4/timeA^2)*(localPhi).*localPhiPrev + ...
                            +(3/timeA^3 - 8/timeA^2 + 6/timeA)*localPhiPrev.^2 - 2/timeA ) ...
                            - lagPar*0.5*( (1/(3*timeA^2)).*localPhi + (1/3)*(-2/timeA^2 + 3/timeA).*localPhiPrev ); 
            source = -0.25*( (-1/timeA^3 + 4/timeA^2 - 6/timeA + 4).*localPhiPrev.^3 + (2/timeA - 4).*localPhiPrev) + ...
                            + lagPar*0.5*( (1/3)*(1/timeA^2 - 3/timeA + 3).*(localPhiPrev.^2) - 1 ) ;     

            localGradPhiX  = sum(DNDx.*phi,2);
            localGradPhiY  = sum(DNDy.*phi,2);
            if (Flag3D==1)
            localGradPhiZ  = sum(DNDz.*phi,2);
            Resphi = localPhiDot + locUX.*localGradPhiX + locUY.*localGradPhiY + ...
                     locUZ.*localGradPhiZ + reac.* localPhi  - source ;
            absGradPhi = sqrt(localGradPhiX.^2 + localGradPhiY.^2 + localGradPhiZ.^2) ;
            modVel = sqrt(locUX.^2 + locUY.^2 + locUZ.^2) ;
            else
            Resphi = localPhiDot + locUX.*localGradPhiX + locUY.*localGradPhiY + ...
                     reac.* localPhi  - source ;
            absGradPhi = sqrt(localGradPhiX.^2 + localGradPhiY.^2) ;
            modVel = sqrt(locUX.^2 + locUY.^2) ;
            end
            etaRPhi(:,p) = Resphi ;
            xx(:,p) = localPhi - localPhiPrev ;
        end
        fsT = sum(etaRPhi,2) ;
        xx = sum(xx,2) ;

        theta = 0.5 ;
        theta_c = 0.05 ;

	% Compute the error indicator for each element
        etaR = computeEtaR(Sol,fsT,solver.epsilon^2) ;

        activeNode = find(Sol.type~=0); 
        nonActive = setdiff(size(crd,1),activeNode);
        fprintf('EtaR: %e, ',sqrt(sum(etaR))) ;
        fprintf('nElem: %d, ',size(Sol.elem,1)) ;
        fprintf('ndofs: %d\n',size(activeNode,1)) ;
        
	% Convergence criteria for Nonlinear Adaptive Variational Partitioned (NAVP) procedure
        flagCoarsen = 0 ;
        flagRefine = 0 ;
        flagBreak = 0 ;
        if ( sqrt(sum(etaR)) < solver.nLTol2  && normIndicator1 <= solver.nLTol1 && normIndicator2 <= solver.nLTol1 )
            flagRefine = 0 ;
            flagCoarsen = 1 ;
            flagBreak = 1 ;
        elseif ( sqrt(sum(etaR)) < solver.nLTol2 && normIndicator1 > solver.nLTol1 && normIndicator2 > solver.nLTol1 )
            flagRefine = 0 ;
            flagCoarsen = 0 ;
        elseif ( nLiter == solver.nLIterMax )
            flagRefine = 0 ;
            flagCoarsen = 1 ;
            flagBreak = 1 ;
        elseif ( sqrt(sum(etaR)) > solver.nLTol2 )
            flagRefine = 1 ;
            flagCoarsen = 0 ;
        end

        crd(:,3) = [];
        
	% Coarsen the grid
        if (flagCoarsen == 1)
        [Sol,etaR] = coarsening(Sol,sqrt(etaR),theta_c);
        etaR = etaR.^2 ;
        end
        
	% Refine the grid through bisection
        if (flagRefine == 1)
        [Sol,etaR] = bisection(Sol,sqrt(etaR),theta);
        etaR = etaR.^2 ;
        end
        
	% Update all coordinates and connectivity, dof, nElem, etc.
        crd = [Sol.node zeros(size(Sol.node,1),1)] ;
        clear cnn ;
        cnn = Sol.elem ;
        ndof = size(crd,1) ;
        nElem = size(cnn,1) ;
        clear xxf yyf zzf src localSrc etaRPhi RhsVec localRhs xx
        clear Sol.phi Sol.phiDot Sol.phiPrev Sol.phiDotPrev
        
        Sol.u = Sol.uPrev + (1/pmc.alpha)*( Sol.uAlpha - Sol.uPrev ) ;
        Sol.uDot = Sol.uDotPrev + (1/pmc.alphaM)*( Sol.uDotAlpha - Sol.uDotPrev ) ;
        Sol.phi = Sol.phiPrev + (1/pmc.alpha)*( Sol.phiAlpha - Sol.phiPrev ) ;
        Sol.phiDot = Sol.phiDotPrev + (1/pmc.alphaM)*( Sol.phiDotAlpha - Sol.phiDotPrev ) ;
        
	% Boundary nodes on the new grid
        zerotypenodes = find(Sol.type == 0);
        bc.left.nodes = find(Sol.node(:,1)< 0.0+1e-8) ;
        bc.left.nodes = setdiff(bc.left.nodes, zerotypenodes)' ;
        bc.right.nodes = find(Sol.node(:,1)> 1.0-1e-8) ;
        bc.right.nodes = setdiff(bc.right.nodes, zerotypenodes)' ;
        bc.top.nodes = find(Sol.node(:,2)> 1.5-1e-8) ;
        bc.top.nodes = setdiff(bc.top.nodes, zerotypenodes)' ;
        bc.bottom.nodes = find(Sol.node(:,2)< 0.0+1e-8) ;
        bc.bottom.nodes = setdiff(bc.bottom.nodes, zerotypenodes)' ;
        % Satisfy Dirichlet boundary condition on the new grid for velocity and pressure
        if (bc.left.type == 1)
            dirichletNodes = bc.left.nodes' ;
            if (Flag3D == 1)
            Sol.u(dirichletNodes,bc.left.var(1),1) = bc.left.value(1).*ones(size(bc.left.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.left.var(2),1) = bc.left.value(2).*ones(size(bc.left.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.left.var(3),1) = bc.left.value(3).*ones(size(bc.left.nodes',1),1) ;
            else
            Sol.u(dirichletNodes,bc.left.var(1),1) = bc.left.value(1).*ones(size(bc.left.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.left.var(2),1) = bc.left.value(2).*ones(size(bc.left.nodes',1),1) ;
            end
        elseif (bc.left.type == 3)
            dirichletNodes = bc.left.nodes' ;
            Sol.u(dirichletNodes,1,1) = zeros(size(bc.left.nodes',1),1) ;
        end
        if (bc.right.type == 1)
            dirichletNodes = bc.right.nodes' ;
            if (Flag3D == 1)
            Sol.u(dirichletNodes,bc.right.var(1),1) = bc.right.value(1).*ones(size(bc.right.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.right.var(2),1) = bc.right.value(2).*ones(size(bc.right.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.right.var(3),1) = bc.right.value(3).*ones(size(bc.right.nodes',1),1) ;
            else
            Sol.u(dirichletNodes,bc.right.var(1),1) = bc.right.value(1).*ones(size(bc.right.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.right.var(2),1) = bc.right.value(2).*ones(size(bc.right.nodes',1),1) ;
            end
        elseif (bc.right.type == 3)
            dirichletNodes = bc.right.nodes' ;
            Sol.u(dirichletNodes,1,1) = zeros(size(bc.right.nodes',1),1) ;
        end
        if (bc.top.type == 1)
            dirichletNodes = bc.top.nodes' ;
            if (Flag3D == 1)
            Sol.u(dirichletNodes,bc.top.var(1),1) = bc.top.value(1).*ones(size(bc.top.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.top.var(2),1) = bc.top.value(2).*ones(size(bc.top.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.top.var(3),1) = bc.top.value(3).*ones(size(bc.top.nodes',1),1) ;
            else
            Sol.u(dirichletNodes,bc.top.var(1),1) = bc.top.value(1).*ones(size(bc.top.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.top.var(2),1) = bc.top.value(2).*ones(size(bc.top.nodes',1),1) ;
            end
        elseif (bc.top.type == 3)
            dirichletNodes = bc.top.nodes' ;
            Sol.u(dirichletNodes,2,1) = zeros(size(bc.top.nodes',1),1) ;
        end
        if (bc.bottom.type == 1)
            dirichletNodes = bc.bottom.nodes' ;
            if (Flag3D == 1)
            Sol.u(dirichletNodes,bc.bottom.var(1),1) = bc.bottom.value(1).*ones(size(bc.bottom.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.bottom.var(2),1) = bc.bottom.value(2).*ones(size(bc.bottom.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.bottom.var(3),1) = bc.bottom.value(3).*ones(size(bc.bottom.nodes',1),1) ;
            else
            Sol.u(dirichletNodes,bc.bottom.var(1),1) = bc.bottom.value(1).*ones(size(bc.bottom.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.bottom.var(2),1) = bc.bottom.value(2).*ones(size(bc.bottom.nodes',1),1) ;
            end
        elseif (bc.bottom.type == 3)
            dirichletNodes = bc.bottom.nodes' ;
            Sol.u(dirichletNodes,2,1) = zeros(size(bc.bottom.nodes',1),1) ;
        end
        if (bc.side_1.type == 1)
            dirichletNodes = bc.side_1.nodes' ;
            Sol.u(dirichletNodes,bc.side_1.var,1) = bc.side_1.value.*ones(size(bc.side_1.nodes,1),1) ;
        elseif (bc.side_1.type == 3)
            dirichletNodes = bc.side_1.nodes' ;
            Sol.u(dirichletNodes,3,1) = zeros(size(bc.side_1.nodes,1),1) ;
        end
        if (bc.side_2.type == 1)
            dirichletNodes = bc.side_2.nodes' ;
            Sol.u(dirichletNodes,bc.side_2.var,1) = bc.side_2.value.*ones(size(bc.side_2.nodes,1),1) ;
        elseif (bc.side_2.type == 3)
            dirichletNodes = bc.side_2.nodes' ;
            Sol.u(dirichletNodes,3,1) = zeros(size(bc.side_2.nodes,1),1) ;
        end
        
        % 2D dirichlet condition in 3D mesh 
        if (Flag3D == 1)
        Sol.u(:,3,1) = zeros(size(crd,1),1) ;
        end
        Sol.p(bc.top.nodes,1) = 0.0 ;
        
        clear etaRPhi phi phiDot convVel phiPrev xsInv gijDown rhs indicator iif jjf
        
	% Break the code (Come out of nonlinear iteration loop)
        if (flagBreak == 1)
            break;
        end
        
    end
    fprintf('\n');

    % Copy current time solution to previous time solution
    Sol.uPrev = Sol.u ;
    Sol.uDotPrev = Sol.uDot ;
    Sol.pPrev = Sol.p ;
    Sol.phiPrev = Sol.phi ;
    Sol.phiDotPrev = Sol.phiDot ;

    % Output in *.plt format to visualize in tecplot 
    if (mod(timeStep,solver.outFreq)==0)
    crd1 = [[1:size(Sol.node,1)]' Sol.node];
    uu = Sol.u ;
    pp = Sol.p ;
    phiphi = Sol.phi ;
    elemelem = Sol.elem ;
    
    indxFree = find(Sol.type==0) ;
    crd1(indxFree,:) = [];
    uu(indxFree,:,:) = [];
    pp(indxFree,:) = [];
    phiphi(indxFree,:) = [];
    indx2 = find(crd1(:,1)~=0) ;
    map = [crd1(:,1) indx2];
    elemelem = elemelem(:) ;
    [~,ii] = ismember(elemelem,map(:,1));
    elemelem = reshape(ii,[],3);
    crd1(:,1) = [] ;
        
    data = [crd1 uu pp phiphi];
    clear uu pp phiphi ;
    filename = sprintf('%s/%s.%d.plt',wrkDir,problemString,timeStep);
    fileId = fopen(filename,'w');

    FileTitle = 'multiphase plot';

    fprintf(fileId,' TITLE = \"%s\"\n',FileTitle);

    if (strcmp(FileTitle,'multiphase plot') & Flag3D == 1)
        fprintf(fileId,' VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"p\", \"phi\"\n');
    elseif (strcmp(FileTitle,'multiphase plot') & Flag3D ~= 1)
        fprintf(fileId,' VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"p\", \"<greek>f</greek>\"\n');
    end
    
    if (strcmp(elemType,'3Tri'))
        fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',(timeStep-1)*solver.dt,size(crd1,1),size(elemelem,1));
    elseif (strcmp(elemType,'4Quad'))
        fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n',(timeStep-1)*solver.dt,size(crd1,1),size(elemelem,1));
    elseif (strcmp(elemType,'4Tet'))
        fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON\n',(timeStep-1)*solver.dt,size(crd1,1),size(elemelem,1));
    elseif (strcmp(elemType,'6Prism'))
        fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n',(timeStep-1)*solver.dt,size(crd1,1),size(elemelem,1));
    elseif (strcmp(elemType,'8Hex'))
        fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n',(timeStep-1)*solver.dt,size(crd1,1),size(elemelem,1));
    end
    
    if (Flag3D == 1)
        fprintf(fileId,'%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n',data');
    else
        fprintf(fileId,'%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n',data');
    end

    if (strcmp(elemType,'3Tri'))
        fprintf(fileId,'%d %d %d\n',elemelem');
    elseif (strcmp(elemType,'4Quad'))
        fprintf(fileId,'%d %d %d %d\n',elemelem');
    elseif (strcmp(elemType,'4Tet'))
        fprintf(fileId,'%d %d %d %d\n',elemelem');
    elseif (strcmp(elemType,'6Prism'))
        fprintf(fileId,'%d %d %d %d %d %d\n',elemelem');
    elseif (strcmp(elemType,'8Hex'))
        fprintf(fileId,'%d %d %d %d %d %d %d %d\n',elemelem');
    end
    clear crd1 elemelem ;
    fclose(fileId);
    end
    
    toc
end
toc
fprintf('\n');
