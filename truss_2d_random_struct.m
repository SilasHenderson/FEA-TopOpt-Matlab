% Truss Opt Test (Random Search) Silas Henderson
clear; close all;

truss = trussMake();
% ---------------------------- Loop ------------------------------ %
tic;
while toc < 10
    truss = areaSearch(truss);
    truss = plotU(truss);
    pause(.01);
end

% ------------------------ Build Truss -------------------------- %
function truss = trussMake()
    truss.el   = [1, 2; 2, 3];      truss.n = [0, 0; 1, 1; 0, 1];   
    truss.area = 10*rand(1, 2);     truss.dof = 3:4;
    truss.E    = 100;               truss.F = [0, 0, 0, -20, 0, 0]';
    truss.C    = 100000000; 
    
    truss.K0El = zeros(numel(truss.n), numel(truss.n));
    for e = 1:numel(truss.el)/2   
        n1 = truss.el(e, 1);
        n2 = truss.el(e, 2);   
        dx = truss.n(n2,1) - truss.n(n1,1);
        dy = truss.n(n2,2) - truss.n(n1,1);  
        truss.len(e) = norm([dx, dy]);
        k4 = [dx*dx, dx*dy;
              dx*dy, dy*dy]*truss.E/truss.len(e)^3;
        truss.K0El([2*n1-1, 2*n1, 2*n2-1, 2*n2], ...
                   [2*n1-1, 2*n1, 2*n2-1, 2*n2], e) = [k4,-k4;-k4,k4];
    end
    
    truss.ax = axes('XLim', [-.2, 1.2], 'XGrid', 'on', ...
                    'YLim', [-.2, 1.2], 'YGrid', 'on');                
    for e = 1:numel(truss.el)/2
        truss.linesU(e) = line('color', 'red',  'parent', truss.ax);
    end
end
 
% ------------------------ Test Random Areas --------------------- %
function truss = areaSearch(truss)
    testArea = truss.area + [1, -1]*(1 - 2*rand)*.1;
    testArea(testArea < 0) = 0;
    FRed = truss.F(truss.dof);
    
    K = zeros(numel(truss.dof), numel(truss.dof));
    for e = 1:numel(truss.el)/2
        K = truss.K0El(truss.dof, truss.dof, e)*testArea(e) + K;
    end
    
    truss.U = zeros(numel(truss.n), 1);
    truss.U(truss.dof) = inv(K)*FRed;
    C = truss.F'*truss.U;
    if C < truss.C
        truss.C = C;
        truss.area = testArea;
    end
end

%----------------------------- Plot Truss ---------------------- %
function truss = plotU(truss)
nodeU = truss.n;
for n = 1:numel(truss.n)/2
     nodeU(n,1) = truss.n(n,1) + truss.U(2*n-1, 1);
     nodeU(n,2) = truss.n(n,2) + truss.U(2*n  , 1);
end

for e = 1:numel(truss.el)/2
    n1 = truss.el(e,1);
    n2 = truss.el(e,2);                       
    set(truss.linesU(e), ...
        'XData', [nodeU(n1,1), nodeU(n2,1)], ...
        'YData', [nodeU(n1,2), nodeU(n2,2)], ...
        'linewidth', truss.area(e));
end
end
