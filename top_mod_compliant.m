%   - generate lower half of a 'gripper' mechanism

%   This is modified version of Sigmund's 99-line topopt script. 
%   It is much slower and probably contains errors.

%   Also, 'k0' matrix is derived w/ symbolics within script

%   Re-organized sections as follows:

%       Options: 
%           choose size, p, filter radius, etc.
%       
%       Optimization: 
%           run functions in loop
%       
%       Functions: 
%           - Assemble k0 matrix
%           - Finite element analysis to find U
%           - Find compliance and compliance gradient
%           - Filter compliance gradient
%           - Optimize element densities (OC)
%           - Plot element densities

clear;

%--------------------------- Options-----------------------------------

% Mesh
nelx = 60;
nely = 30;
volfrac = .2;
x = ones(nely,nelx);

% Optimization
p = 3;
filterradius = 1.5;
iterations = 100;

% Force (in=1, out=2)
F = sparse(2*(nely+1)*(nelx+1),2); 

% Force in (column 1)
F(1,1) = 1;            

% Force out (column 2)
F( (2*nely+1)*nelx+20, 2) = 1;

% Fixed degrees of freedom
fixed_dof  = union( ...
     2:2*(nely+1):2*(nely+1)*(nelx+1)*3/4, ...
     2*(nely+1)-1:2*(nely+1));

% Passive elements
x_passive = zeros(nely,nelx);
x_passive(1:nely/3, 2/3*nelx:nelx) = 1;

min_x = .001;
x(find(x_passive)) = min_x;
  
% Make 'k0' matrix

k0 = Derive_k0;    

%---------------------------- Optimization -------------------------------

% Run iterations

for iter = 1:iterations 
    
    % Find displacements
    U = FindDisplacement(nelx,nely,x,p,k0,F,fixed_dof);                    
    
    % Find compliance (c) and compliance 'derivative' (dc)
    [c,dc] = FindCompliance(nelx,nely,p,x,k0,U);                 
    
    % Filter compliance derivative
    dc = ConeFilter(nelx,nely,filterradius,x,dc);       
    
    % Update Areas
    x  = OC(nelx,nely,x,volfrac,dc,x_passive,min_x);        
    
    % Plot current structure
    plotx(x); 
    
end                                                         

% -------------------------- Functions -------------------------------

% Derive_k0: derive w/symbolics, return non-symbolic matrix

function k0 = Derive_k0
    
    % Poisson's ratio, elastic modulus
    nu = .3; 
    E  =  1;                 
    
    % Natural coordinates
    syms xi eta                        

    % Shape function vector 
    N1234 = 1/4*[ ...               
        (1-xi)*(1-eta);          
        (1+xi)*(1-eta);
        (1+xi)*(1+eta);
        (1-xi)*(1+eta)];  
    
    Nxi   = diff(N1234,xi);  
    Neta  = diff(N1234,eta);

    % Strain-displacement matrix
    B = [ ...
        Nxi(1),  0,       Nxi(2),  0,       Nxi(3),  0,       Nxi(4),  0    ;
        0,       Neta(1), 0,       Neta(2), 0,       Neta(3), 0,       Neta(4);
        Neta(1), Nxi(1),  Neta(2), Nxi(2),  Neta(3), Nxi(3),  Neta(4), Nxi(4)];
 
    % Elasticity matrix
    Cmod = 1/(1-nu^2).*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    J = 4;
    
    % k0 matrix
    k0 = int(int(B'*Cmod*B,xi,-1,1),eta,-1,1); 
    k0 = J*E*double(k0);          
end

% Find Displacement

function U = FindDisplacement(nelx,nely,x,penal,k0,F,fixeddofs)
    
    % Init K and U
    K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
    U = zeros( 2*(nely+1)*(nelx+1),2);
   
    % Assemble K global
    for elx = 1:nelx                
        for ely = 1:nely
            
            upleftnode  = (nely+1)*(elx-1)+ely; 
            uprightnode = (nely+1)* elx   +ely;
            
            el_dof = [ ...
                2*upleftnode-1;  2*upleftnode;    
                2*uprightnode-1; 2*uprightnode; 
                2*uprightnode+1; 2*uprightnode+2; 
                2*upleftnode+1;  2*upleftnode+2];
            
            K(el_dof,el_dof) = ...
                K(el_dof,el_dof) + x(ely,elx)^penal*k0;  
        end
    end
    
    % Add extra stiffness to in/out pts
    extra1 = find(F(:,1));
    extra2 = find(F(:,2));
    K(extra1,extra1) = K(extra1,extra1) + .1;
    K(extra2,extra2) = K(extra2,extra2) + .1;
    
    % Apply boundary conditions, solve
    alldofs     = 1:2*(nely+1)*(nelx+1);
    freedofs    = setdiff(alldofs,fixeddofs);
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);     
    U(fixeddofs,:)= 0;
end

% FindCompliance: find compliance and compliance gradient

function [c,dc] = FindCompliance(nelx,nely,p,x,k0,U)                                
    
    c  = 0.; 
    dc = zeros(nely,nelx);   
    
    % Loop across each element
    for ely = 1:nely                
        for elx = 1:nelx
            
            n1 = (nely+1)*(elx-1)+ely;    
            n2 = (nely+1)* elx   +ely;
            Ue1 = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
            Ue2 = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],2);
        
            % add element c's
            c = c + x(ely,elx)^p*Ue2'*k0*Ue1;   
            
            % assign dc/dx el entry to dc matrix
            dc(ely,elx) = -p*x(ely,elx)^(p-1)*Ue2'*k0*Ue1;
        end
    end
end

% Conefilter:  filter compliance gradient

function dcfiltered = ConeFilter(nelx,nely,radius,x,dc)
        
    % Average each elem's dc with weighted dc's from neighbor elements
    dcfiltered=zeros(nely,nelx);
    
    % Loop over each element 
    for i = 1:nelx                                     
        for j = 1:nely                              
            
            % get 'weight' (cone height)
            cone_height_sum=0.0;             
            for k = max(i-floor(radius),1):min(i+floor(radius),nelx)
                for l = max(j-floor(radius),1):min(j+floor(radius),nely) 
                    
                    coneheight = radius-sqrt((i-k)^2+(j-l)^2);      
                    cone_height_sum = cone_height_sum+max(0,coneheight);
                    
                    dcfiltered(j,i) = dcfiltered(j,i) + ...
                        max(0,coneheight)*x(l,k)*dc(l,k);
                end
            end
            
            % Average element dc and neighbors dc
            dcfiltered(j,i) = dcfiltered(j,i)/(x(j,i)*cone_height_sum);
        end
    end
end

% OC: Update Areas with Optimality Criteria Method

function [x_new]=OC(nelx,nely,x,vol_frac,dc,x_passive,x_min)  

    % Lambda min, max, and 'move' parameter
    L_min = 0; 
    L_max = 100000; 
    move = 0.1;
    
    % Update areas with bisection method
        
    while (L_max-L_min)/(L_max+L_min)>1e-4 && L_max+.0001>1e-40
        
        x_new(find(x_passive)) = x_min;
        
        lambda = (L_max + L_min)/2;
        
        x_new = max(0.001, ...
            max(x-move, ...
                min(1., ...
                    min(x+move,x.*(max(1e-10,-dc./lambda)).^0.3))));
    
        if sum(sum(x_new)) - vol_frac*nelx*nely > 0
            L_min = lambda;
        else
            L_max = lambda;
        end
        
        x_new(find(x_passive)) = x_min;
    end 
end

% PlotUpdate: draw structure with current element densities

function plotx(x)
    cla;
    colormap(gray); 
    imagesc(-x); 
    pause(.0001);
end
