% Derive k0 matrix from Sigmund 99-line TopOpt, with symbolics

% poisson's ratio, elastic modulus
nu = .3; 
E  = 1;                     

% natural coordinates  
syms xi eta                           

% shape function vector    
N1234 = 1/4*[ ...
    (1-xi)*(1-eta);       
    (1+xi)*(1-eta);
    (1+xi)*(1+eta);
    (1-xi)*(1+eta)];      
         
% shape function derivative
Nxi   = diff(N1234,xi);  
Neta  = diff(N1234,eta);

% strain-displacement matrix
B = [ ...
    Nxi(1),  0,       Nxi(2),  0,       Nxi(3),  0,       Nxi(4),  0;
    0,       Neta(1), 0,       Neta(2), 0,       Neta(3), 0,       Neta(4);
    Neta(1), Nxi(1),  Neta(2), Nxi(2),  Neta(3), Nxi(3),  Neta(4), Nxi(4)];
 
% elasticity matrix
Cmod = 1/(1-nu^2).*[...
    1,  nu, 0; 
    nu, 1,  0; 
    0,  0, (1-nu)/2];

% k0 matrix
J  = 4;
k0 = int(int(B'*Cmod*B,xi,-1,1),eta,-1,1); 
k0 = J*E*double(k0);     

% print
disp(k0);
