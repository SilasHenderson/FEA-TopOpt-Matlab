// Notes on Vladmir Uskov's Matlab Code

clear; clc; close all;
set(gcf, 'color', [.6, .8, .7], 'menubar', 'none');
set(gca, 'units', 'normal', 'position', [.02, 0, .96, 1]);
topo2Skip(100,30,0.4,1.8,100);

function topo2Skip(nelx,nely,volfrac,rs,nloop)
penal = 3;                            % Penalization factor 
eps0 = 1e-3;                          % Skip level
eta = 0.1;                            % Damping coefficient, initial
etamax = 1;                           % and max values
nu = 0.3;                             % Poisson's ratio
%% Isoparametric element stiffness matrix
a = [12 3 -6 3 0 -6 -3 -3];
b = [-4 3 -2 -9 4 2 -3 9];
k = (a+nu*b)/(24*(1-nu^2));
i1 = [1 2 3 8; 2 1 4 5; 3 4 1 7; 8 5 7 1];
i2 = [6 7 5 4; 7 6 8 3; 5 8 6 2; 4 3 2 6];
KE = k([i1 i2; i2 i1]);         
%% Prepare for assembly of the global stiffness matrix
nel = nelx*nely;                      % Number of elements
nnodes = (1+nelx)*(1+nely);           % Number of nodes
ndof = 2*nnodes;                      % Number of DOFs
nodenrs = reshape(1:nnodes,1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,[],1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nel,1);
iK = reshape(kron(edofMat,ones(8,1))',[],1);
jK = reshape(kron(edofMat,ones(1,8))',[],1);
%% Loads and supports for half MBB beam
F = sparse(2,1,-1,ndof,1);
fixeddof = union(1:2:2*(nely+1),ndof);
%% Initialization
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
lm = 0;
reta = (etamax/eta)^(1/nloop);        % Augment of damping
freedof = setdiff(1:ndof,fixeddof);   % Free DOFs
U = zeros(ndof,1);                    % Displacements
x = repmat(volfrac,nel,1);            % Densities
%% Iterations
for loop = 1:nloop
	%% Gaussian filtering of densities
	x = imgaussfilt(reshape(x,nely,nelx),rs);	x = x(:);
  %% Sensitivity analysis
	y = x.^penal;	                      % Young�s modulus
	y(y < eps0*volfrac) = 0;            % is 0 for skipped elements
	K = sparse(iK,jK,KE(:)*y');         % Global stiffness matrix
	d = diag(K);                        % Diagonal
	skip = find(d == 0);                % Skipped DOFs
  remain = setdiff(freedof,skip);     % Remaining DOFs
	K = K(remain,remain);
  U(remain) = (K+K')\F(remain)*2;     % Solve for remaining DOFs only
  ue = U(edofMat);
  ce = max(0,y.*sum(ue*KE.*ue,2));    % Element energy
  %% Optimality criteria update of densities
	x = x.*(ce.^eta);
	xNew = @(a) min(1,x*exp(a));
	delta = @(a) sum(xNew(a)) - volfrac*nel;
	lm = fzero(delta,lm);
	x = xNew(lm);
	eta = min(etamax,eta*reta);         % New value of damping coefficient
  %% Print info
  fprintf('%u\t C: %.6g\t Skip: %u%%\n', ...
    loop,sum(ce),round(100*numel(skip)/numel(freedof)));
  %% Plot densities
	imagesc(reshape(1-x,nely,nelx));
  colormap(hot); caxis([0 1]); axis equal; axis off; drawnow;
end
end
