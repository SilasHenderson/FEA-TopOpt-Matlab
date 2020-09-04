%   Compliant mechanism synthesis:
%   - generate lower half of a 'gripper' mechanism

%   This is a heavily modified version of Sigmund's 99-line
%   topopt script.  It is slower, and may contain errors.

%   Re-organized sections as follows:

%       (A) options: choose size, p, filter radius, etc.
%       (B) optimization: run functions (1),(2),(3),(4),(5) in 'for loop'
%       (C) functions: 
%           (0)     assemble k0 matrix
%           (1)     finite element analysis to find U
%           (2)     find compliance and compliance gradient
%           (3)     filter compliance gradient
%           (4)     optimize element densities (OC)
%           (5)     plot element densities

clear
%------------------------(A) options-------------------------------------%

% Mesh properties
nelx = 90;
nely = 60;
volfrac = .2;
x = rand*ones(nely,nelx);

% Optimization options
p = 3;
filterradius = 1.5;
iterations = 60;

% Force vectors
F = sparse(2*(nely+1)*(nelx+1),2); 

% Force in (column 1)
F(1,1) = 1 ;            

% Force out (column 2)
F((2*nely+1)*(nelx)+2*10,2) = 1;

% Fixed degrees of freedom
fixeddofs1 = 2:2*(nely+1):2*(nely+1)*(nelx+1)*3/4;
fixeddofs2 = 2*(nely+1)-1:2*(nely+1);
fixeddofs = union(fixeddofs1,fixeddofs2);

% Passive elements
passive = zeros(nely,nelx);
passive(1:nely/3,2/3*nelx:nelx) = 1;
minp = .001;
x(find(passive)) = minp;
  
%------------------------(B) optimization--------------------------------%

k0 = k0make;                                                 %(0) 

for iter = 1:iterations                                      % start loop  

    [U]    = fea(nelx,nely,x,p,k0,F,fixeddofs);              % (1) displacement        
    [c,dc] = compliance(nelx,nely,p,x,k0,U);                 % (2) compliance
    [dc]   = conefilter(nelx,nely,filterradius,x,dc);        % (3) filter
    [x]    = OC(nelx,nely,x,volfrac,dc,passive,minp);        % (4) update
    plotx(x);                                                % (5) plot
    iter;
end                                                          % end loop

% ----------------------(C) functions------------------------------------%

% (0a) 'k0make':  assemble k0 matrix

function k0 = k0make
    
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

% (1) 'fea':  finite element analysis for finding U

function [U] = fea(nelx,nely,x,penal,k0,F,fixeddofs)
    
    % Init K and U
    K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
    U = zeros(2*(nely+1)*(nelx+1),2);
   
    % Assemble K global
    for elx = 1:nelx                
        for ely = 1:nely
            
            upleftnode  = (nely+1)*(elx-1)+ely; 
            uprightnode = (nely+1)* elx   +ely;
            
            edof = [ ...
                2*upleftnode-1;  2*upleftnode;    
                2*uprightnode-1; 2*uprightnode; 
                2*uprightnode+1; 2*uprightnode+2; 
                2*upleftnode+1;  2*upleftnode+2];
            
            K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*k0;  
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

% (2) 'compliance' :  find compliance and compliance gradient

function [c,dc] = compliance(nelx,nely,p,x,k0,U)                                
    
    c = 0.; 
    dc = zeros(nely,nelx);   
    
    % Loop over each element
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

% (3) 'conefilter':  filter compliance gradient

function [dcfiltered]=conefilter(nelx,nely,radius,x,dc)
        
    % Average each element dc with weighted dcs from neighbor elements
    dcfiltered=zeros(nely,nelx);
    
    % Loop over each element 
    for i = 1:nelx                                     
        for j = 1:nely                              
            
            % Weight (cone height)
            coneheightsum=0.0;             
            for k = max(i-floor(radius),1):min(i+floor(radius),nelx)
                for l = max(j-floor(radius),1):min(j+floor(radius),nely) 
                    coneheight = radius-sqrt((i-k)^2+(j-l)^2);      
                    coneheightsum = coneheightsum+max(0,coneheight);
                    dcfiltered(j,i) = dcfiltered(j,i) + max(0,coneheight)*x(l,k)*dc(l,k);
                end
            end
            
            % Average element dc and neighbors dc
            dcfiltered(j,i) = dcfiltered(j,i)/(x(j,i)*coneheightsum);
        end
    end
end

% (4) 'optimize' : bisection method to find good lambda for dC = dC/lambda

function [xnew]=OC(nelx,nely,x,volfrac,dc,passive,minp)  

    lambdamin = 0; 
    lambdamax = 100000; 
    move = 0.1;
    
    while (lambdamax-lambdamin)/(lambdamax+lambdamin) > 1e-4 && ...
           lambdamax + .0001 > 1e-40
    
        xnew(find(passive)) = minp;
        lambda = 0.5*(lambdamax+lambdamin);
        xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*...
            (max(1e-10,-dc./lambda)).^0.3))));
    
        if sum(sum(xnew)) - volfrac*nelx*nely > 0
            lambdamin = lambda;
        else
            lambdamax = lambda;
        end
        
        xnew(find(passive)) = minp;
    end 
end

% (5) plotx :  plot element densities

function plotx(x)
    colormap(gray); 
    imagesc(-x); 
    pause(.0001);
end
