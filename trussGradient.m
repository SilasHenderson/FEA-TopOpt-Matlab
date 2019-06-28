clear; clc; close all;

% ----------------- Define Elements, Nodes, Forces, Dofs -------------- %

for n =  1:20, nodes(n,:) = [n,        0]; end
for n = 21:40, nodes(n,:) = [n - 20, -.5]; end

nodes(41, :) = [6,  3];
nodes(42, :) = [14, 3];

for e =  1:19
    elements(4*e - 3, :) = [     e,  e +   1];
    elements(4*e - 2, :) = [     e,  e +  21];
    elements(4*e - 1, :) = [20 + e,  e +   1];
    elements(4*e    , :) = [20 + e , e +  21];
end

elements(77,:) = [ 1, 41];   elements(78,:) = [41, 10];
elements(79,:) = [41, 42];   elements(80,:) = [10, 42];
elements(81,:) = [42, 20]; 

dofs   = 1:84;               dofs([1, 2, 39, 40, 41, 42, 79, 80]) = [];    
forces = zeros(84, 1);       forces(42:2:80,1) =  -1;
 
truss    = trussMake(elements, nodes, dofs, forces);
truss    = kSolve(truss);
truss.C0 = sum(truss.c);

% ----------------------------- Loop --------------------------------- %
tic;
while toc < 20
    truss.area0 = truss.area;
    truss.c0    = truss.c;
              
    truss.area  = truss.area + truss.c0./truss.area*rand;
    truss.area(truss.area < .01) = .001;     
    truss.area = truss.area/dot(truss.area, truss.len)*truss.maxVol; 
    
    truss = kSolve(truss);   

    if sum(truss.c) > sum(truss.c0)
        truss.area = truss.area0; 
    end
    
    uPlot(truss);
    truss.iter = truss.iter + 1;    
    pause(.0001);
end

%------------------------ Build Truss -------------------------- %
function truss = trussMake(elements, nodes, dofs, forces)
    truss.n    = nodes;              truss.E = 1000;
    truss.el   = elements;           truss.C = 1000000;  
    truss.dof  = dofs;               truss.F = forces; 
    truss.dofN = numel(dofs);        truss.f = truss.F(truss.dof);
          
    for e = 1:numel(truss.el)/2  

        n1 = truss.el(e, 1);                               
        n2 = truss.el(e, 2);   
        dx = truss.n(n2,1) - truss.n(n1,1);
        dy = truss.n(n2,2) - truss.n(n1,2);  
        truss.len(1,e) = norm([dx, dy]);
        
        k4 = [dx*dx, dx*dy;
              dx*dy, dy*dy]*truss.E/truss.len(e)^3;
   
        k0 = zeros(numel(nodes), numel(nodes));
        k0([2*n1-1, 2*n1, 2*n2-1, 2*n2], ...
           [2*n1-1, 2*n1, 2*n2-1, 2*n2]) = [k4,-k4;-k4,k4];
        truss.k0(:,:,e) = k0(truss.dof,truss.dof);
        truss.k0El(:,:,e) = [k4, -k4; -k4, k4];
    end
    
    truss.maxVol = 200;    
    truss.area   = rand(1, numel(elements)/2);
    truss.area   = truss.area*(truss.maxVol/dot(truss.len, truss.area));
  
    truss.iter = 0;
    truss.fig  = figure('color', [.2, .2, .2], 'menubar',   'none');
    truss.ax   = axes(   'XLim',     [-1, 21],   'XGrid',     'on', ...
                         'YLim',      [-5, 5],   'YGrid',     'on', ...
                       'parent',    truss.fig,   'units', 'normal', ...
                     'position', [0, 0, 1, .9]);     
            
    truss.info = annotation('textbox', 'string', 'bananas', ...
        'units', 'normal', 'position', [0, .9,  1, .1], 'color', 'white');
    for e = 1:numel(truss.el)/2
        truss.linesU(e) = line('color', [.3 .3 .6],  'parent', truss.ax);
    end
end

% ------------------------ Test Random Areas --------------------- %
function truss = kSolve(truss)
    truss.k = zeros(truss.dofN, truss.dofN);
    for e = 1:numel(truss.area)
        truss.k = truss.k + truss.k0(:,:,e).*truss.area(e);
    end
    
    truss.U = zeros(numel(truss.n), 1);
    truss.U(truss.dof) = truss.k\truss.f;
    
    for e = 1:numel(truss.el)/2
        n1 = truss.el(e, 1);
        n2 = truss.el(e, 2);
        uEl  = truss.U([2*n1 - 1, 2*n1, 2*n2 - 1, 2*n2]);
        kEl  = truss.k0El(:, :, e)*truss.area(1, e);
        truss.c(1, e) = uEl'*kEl*uEl;  
    end
end

 %----------------------------- Plot Truss ---------------------- %
function truss = uPlot(truss)
    nodeU = reshape(truss.n', numel(truss.n), 1) + truss.U;
    for e = 1:numel(truss.el)/2
        n1 = truss.el(e,1);
        n2 = truss.el(e,2);      
          
        set(truss.linesU(e), ...
            'XData', [nodeU(2*n1 - 1, 1), nodeU(2*n2 - 1, 1)], ...
            'YData', [nodeU(2*n1,     1), nodeU(2*n2,     1)], ...
            'linewidth', truss.area(e));
    end
    set(truss.info, 'string', ...
    { sprintf('compliance0: %6.3f iterations:%d', truss.C0, truss.iter), ...
      sprintf('complianceI: %6.3f', sum(truss.c))});
end