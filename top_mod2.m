%  User-friendly '----' line code
%     re-organized sections as follows:

%       (A) options: choose size, p, filter radius, etc.
%       (B) optimization: run functions (1),(2),(3),(4),(5) in 'for loop'
%       (C) functions:     
%                       (0)     assemble k0 matrix
%                       (1)     finite element analysis to find U
%                       (2)     find compliance and compliance gradient
%                       (3)     filter compliance gradient
%                       (4)     optimize element densities
%                       (5)     plot element densities

%------------------------(A) options-------------------------------------%

% mesh properties
    nelx = 60;
    nely = 40;
    volfrac = .1;
    x = rand*ones(nely,nelx);

% optimization options
    p = 3;
    filterradius = 1.5;
    iterations = 50;

% force vector
    F = sparse(2*(nely+1)*(nelx+1),1); 
    F(2*(nely+1)*(nelx+1)) = -1;
    F(2*(nely+1)*(nelx/2)) = -5;
    F(2*(nely+1)*nelx+2)   = -.1;

% fixed degrees of freedom
    fixeddofs = 1:nely;

% passive elements
    passive = zeros(nely,nelx);
%     passive(nely/4:nely*3/4,nelx/4:nelx*3/4) = 1;
    minp = .001;
%     x(find(passive)) = minp;

%------------------------(B) optimization--------------------------------%
k0 = k0make;                                                 %(0)                                                            %(0b)
for iter = 1:iterations                                      % start loop    
    [U]    = fea(nelx,nely,x,p,k0,F,fixeddofs);              %(1)         
    [c,dc] = compliance(nelx,nely,p,x,k0,U);                 %(2)
    [dc]   = conefilter(nelx,nely,filterradius,x,dc);        %(3)
    [x]    = OC(nelx,nely,x,volfrac,dc,passive,minp);        %(4)
    plotx(x);                                                %(5)
end                                                          % end loop

% ----------------------(C) functions------------------------------------%

% (0) assemble k0 matrix
function k0 = k0make
    nu = .3; E = 1;                     % poisson's ratio, elastic modulus
    
    syms xi eta                         % natural coordinates    

N1234 = 1/4*[(1-xi)*(1-eta);        % shape function vector     
             (1+xi)*(1-eta);
             (1+xi)*(1+eta);
             (1-xi)*(1+eta)];         
Nxi   = diff(N1234,xi);  Neta  = diff(N1234,eta);

                                    % strain-displacement matrix
B = [Nxi(1)   0     Nxi(2)    0    Nxi(3)  0      Nxi(4)    0    ;
     0       Neta(1)  0     Neta(2)  0     Neta(3)   0    Neta(4);
     Neta(1) Nxi(1) Neta(2) Nxi(2) Neta(3) Nxi(3) Neta(4) Nxi(4)];
 
                                    % elasticity matrix
Cmod = 1/(1-nu^2).*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
J = 4;
k0 = int(int(B'*Cmod*B,xi,-1,1),eta,-1,1); 
k0 = J*E*double(k0);                % k0 matrix --> (dk/dx)
end

% (1) finite element analysis for finding U
function [U]=fea(nelx,nely,x,p,k0,F,fixeddofs)
                                    % K,U set at zeros
K = sparse(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
U = zeros (2*(nely+1)*(nelx+1),1);
                                    % Assemble K global
for elx = 1:nelx                    % loop over elements, add k0els to K
  for ely = 1:nely
    upleftnode   =   (nely+1)*(elx-1) + ely; 
    uprightnode  =   (nely+1)*(elx)   + ely;
    lowrightnode =   (nely+1)*(elx)   + ely + 1;
    lowleftnode  =   (nely+1)*(elx-1) + ely + 1;
    edofs = [2*upleftnode-1;    2*upleftnode;    
             2*uprightnode-1;   2*uprightnode; 
             2*lowrightnode-1;  2*lowrightnode; 
             2*lowleftnode-1;   2*lowleftnode];
    K(edofs,edofs) = K(edofs,edofs) + x(ely,elx)^p*k0;
  end
end
                                    % apply boundary conditions, solve
alldofs       = 1:2*(nely+1)*(nelx+1);
freedofs      = setdiff(alldofs,fixeddofs);
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);     
U(fixeddofs,:)= 0;
end

% (2) find compliance and compliance gradient
function [c,dc] = compliance(nelx,nely,p,x,k0,U)                                
    c = 0.; dc = zeros(nely,nelx);          
    for ely = 1:nely                  % loop over each element,
        for elx = 1:nelx
            
            upleftnode   = (nely+1)*(elx-1) + ely; 
            uprightnode  = (nely+1)*(elx)   + ely;
            lowrightnode = (nely+1)*(elx)   + ely + 1;
            lowleftnode  = (nely+1)*(elx-1) + ely + 1;
            
            edofs = [ ...
                2*upleftnode-1;    2*upleftnode  ;    
                2*uprightnode-1;   2*uprightnode ; 
                2*lowrightnode-1;  2*lowrightnode; 
                2*lowleftnode-1;   2*lowleftnode];
            
            Ue = U(edofs);
                                    % add element c's
            c = c + x(ely,elx)^p*Ue'*k0*Ue;   
                                    % assign dc/dxel entry to dc matrix
            dc(ely,elx) = -p*x(ely,elx)^(p-1)*Ue'*k0*Ue;
        end
    end
end

% (3) filter compliance gradient
function [dcfiltered]=conefilter(nelx,nely,radius,x,dc)
                                    % average each element dc with
                                    % weighted dcs from neighbor elements
    dcfiltered=zeros(nely,nelx);
    for i = 1:nelx                      % loop over each element                
        for j = 1:nely                              
            coneheightsum=0.0;              % weight (cone height)
            for k = max(i-floor(radius),1):min(i+floor(radius),nelx)
                for l = max(j-floor(radius),1):min(j+floor(radius),nely) 
                    coneheight = radius-sqrt((i-k)^2+(j-l)^2);      
                    coneheightsum = coneheightsum+max(0,coneheight);
                    dcfiltered(j,i) = dcfiltered(j,i) + max(0,coneheight)*x(l,k)*dc(l,k);
                end
            end                             % average element dc and neighbors dc
            
            dcfiltered(j,i) = dcfiltered(j,i)/(x(j,i)*coneheightsum);
        end
    end
end

% (4) bisection method to find good lambda for dC = dC/lambda
function [xnew]=OC(nelx,nely,x,volfrac,dc,passive,minp)  

    lambdamin = 0; lambdamax = 100000; move = 0.2;
    while (lambdamax-lambdamin > 1e-4)
        xnew(find(passive)) = minp;
        lambda = 0.5*(lambdamax+lambdamin);
        xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lambda)))));
        if sum(sum(xnew)) - volfrac*nelx*nely > 0
            lambdamin = lambda;
        else
            lambdamax = lambda;
        end
        xnew(find(passive)) = minp;
    end 
end

% (5) plot element densities
function plotx(x)
    colormap(gray); 
    imagesc(-x); 
    axis equal; axis tight; axis off;
    pause(.001);
end
