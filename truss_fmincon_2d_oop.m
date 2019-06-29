clear; clc; close all;  global truss;

% --------------------------- Build Truss -------------------------- %
truss.n = [0, 2; 2, 2; 4, 2;
           0, 0; 2, 0; 4, 0];
     
truss.el = [1, 2; 2, 3; 4, 5; 5, 6; 1, 5;
            2, 6; 2, 4; 5, 3; 2, 5; 3, 6];
        
truss.dof = 1:numel(truss.n);           truss.dof([1, 2, 7, 8]) = [];
truss.F = zeros(numel(truss.n), 1);     truss.F(12) = -5;     
truss.f = truss.F(truss.dof);           truss.E = 1000;
truss.C = 1000000;  

truss.maxVol = 10;    
truss.area = ones(1, numel(truss.el)/2)/numel(truss.el)/2*truss.maxVol;

for e = 1:numel(truss.el)/2  
    truss.line(e) = line(0, 0);
    n1 = truss.el(e, 1);                               
    n2 = truss.el(e, 2);   
    dx = truss.n(n2,1) - truss.n(n1,1);
    dy = truss.n(n2,2) - truss.n(n1,2);  
    truss.len(1,e) = norm([dx, dy]);
        
    k4 = [dx*dx, dx*dy;
          dx*dy, dy*dy]*truss.E/truss.len(e)^3;
    
    k0 = zeros(numel(truss.n), numel(truss.n));
    k0([2*n1-1, 2*n1, 2*n2-1, 2*n2], ...
       [2*n1-1, 2*n1, 2*n2-1, 2*n2]) = [k4,-k4;-k4,k4];
        truss.k0(:,:,e) = k0(truss.dof,truss.dof);
        truss.k0El(:,:,e) = [k4, -k4; -k4, k4];
end
 
lb = zeros(size(truss.area));
ub = 10*ones(size(truss.area));

xOpt = fmincon(@compliance, truss.area, truss.len, truss.maxVol, ...
        [], [], lb, ub);

truss.area = xOpt;  

trussPlot;
% --------------------- Compliance Function -------------------------- %
function C = compliance(x)
    
    global truss; truss.area = x;
    truss.k = zeros(truss.dofN, truss.dofN);
    for e = 1:numel(truss.area)
        truss.k = truss.k + truss.k0(:,:,e).*x(e);
    end
    
    truss.U = zeros(numel(truss.n), 1);
    truss.U(truss.dof) = truss.k\truss.f;
    
    C = truss.F'*truss.U;
    
    for e = 1:numel(truss.el)/2
        n1 = truss.el(e, 1);
        n2 = truss.el(e, 2);
        uEl  = truss.U([2*n1 - 1, 2*n1, 2*n2 - 1, 2*n2]);
        kEl  = truss.k0El(:, :, e)*truss.area(1, e);
        dC(1, e) = - uEl'*kEl*uEl;  
    end
end

function trussPlot
    global truss;
    nodeU = reshape(truss.n', numel(truss.n), 1) + truss.U;
    for e = 1:numel(truss.el)/2
    n1 = truss.el(e,1);
    n2 = truss.el(e,2);      
          
    set(truss.line(e), ...
        'XData', [nodeU(2*n1 - 1, 1), nodeU(2*n2 - 1, 1)], ...
        'YData', [nodeU(2*n1,     1), nodeU(2*n2,     1)], ...
        'linewidth', 10*truss.area(e)/max(truss.area));
    end
    drawnow; grid on;
    set(gca, 'units', 'normalized', 'position', [.05 .05 .9 .9]);
end
