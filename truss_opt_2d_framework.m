% --- Truss Opt Visual Testing Framework  (Silas Henderson, IUPUI) ----- %

% This script outlines a minimum framework for visually testing
% optimization routines for 2d truss structures in Matlab.

% All persistent variables are stored in a 'truss' object. This provides
% a flexible way for functions to set and get variables, while eliminating
% risk associated with global variables (unbounded scope).

% Graphics objects are also stored in the 'truss' object.  
% The graphics primitives (figure, axes, line, annotation) are pre-defined.
% During animation, only the line's x-y data and the annotation's text content
% are updated (as opposed to deleting/creating a new plot every iteration).
% This dramatically reduces graphics overhead cost.

% Truss Object (Core Variables)

% -- truss.n:       node-coordinate matrix
% -- truss.el:      element-node matrix
% -- truss.dof:     degrees of freedom
% -- truss.area:    vector of element areas
% -- truss.E:       Young's modulus
% -- truss.F:       global force vector
% -- truss.maxVol:  maximum volume of truss.
% -- truss.iter:    number of iterations performed
% -- truss.K0El:    3-tensor of element K0 global matrices. For example,
%                   truss.K0El(:, :, 1) is element 1's global K0 matrix.
% -- truss.C:       compliance under loading, dot(u, f);

% Truss Object (Graphics Objects)

% -- truss.fig:    figure (graphics container)
% -- truss.ax:     axes   (graphics container)
% -- truss.info:   textBox annotation
% -- truss.linesU: Array of line objects (one line for each element)

% Functions:

% -- trussMake():  Assembles and returns a 'truss' object                  
% -- areaSearch(): Tests a random small change in bar areas.  
%                  If this change results in lower compliance, 
%                  the area persists.
% -- plotU():      Update truss structure representation 
%                  Update data print-out

% Program:         The upper figure window displays information:
%                  -- number of iterations,
%                  -- initial compliance
%                  -- current compliance
%                  The axes plots the deformed truss, with 
%                  linewidths proportional to areas of each bar.

% Instructions:    The parameters in 'setup' and 'loop' can be modified
%                  without modifying the functions.  As long as the
%                  element/node/dof/force assembly is valid, the animation
%                  will run.  Simulation time is set on line 72 (seconds).

% Optimization:    Immediate improvement can be made and visualized, 
%                  by implementing a better algorithm.
%                  For instance, try searching along area's gradient 
%                  w.r.t. compliance.

%-------------------------------- Setup ---------------------------- %
clear; clc; close all;

nodes = ...
        [1, 1; 2, 1;
   0, 0; 1, 0; 2, 0; 3, 0];

elements = [3, 1; 1, 2; 2, 6;
         1, 4; 1, 5; 2, 4; 2, 5;
      3, 4;       4, 5;       5, 6];
      
dofs     = [1, 2, 3, 4,        7,  8, 9, 10];
forces   = [0, 0, 0, 0, 0, 0,  0,-15, 0,-15, 0, 0]';
truss    = trussMake(elements, nodes, dofs, forces);

% ---------------------------- Loop ------------------------------ %
tic;
while toc < 5
    truss = areaSearch(truss);
    truss = plotU(truss);
    pause(.01);
end

% ------------------------ Build Truss -------------------------- %
function truss = trussMake(elements, nodes, dofs, forces)
    truss.el  = elements;       truss.n = nodes;
    truss.dof = dofs;           truss.F = forces;
    truss.E = 100;              truss.C = 100000000; 
    
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
    
    truss.maxVol = 100;          
    truss.area   = 10*rand(1, numel(elements)/2);     
    truss.area   = truss.area*(truss.maxVol/(truss.area*truss.len'));
    
    truss.iter = 0;
    truss.fig = figure('color', [.2, .2, .2], 'menubar', 'none');
    truss.ax = axes('XLim', [-.2, 3.2], 'XGrid', 'on', ...
                    'YLim', [-.5, 1.5], 'YGrid', 'on', ...
                  'parent',  truss.fig, 'units', 'normal', ...
                'position', [0, 0, 1, .9]);     
            
    truss.info = annotation('textbox', 'string', 'bananas', ...
        'units', 'normal', 'position', [0, .9,  1, .1], 'color', 'white');
    for e = 1:numel(truss.el)/2
        truss.linesU(e) = line('color', [.3 .3 .6],  'parent', truss.ax);
    end
end
 
% ------------------------ Test Random Areas --------------------- %
function truss = areaSearch(truss)
    testArea = truss.area + ...
       (2*rand(1, numel(truss.area)) - ones(1, numel(truss.area)))*rand^2;
    
    testArea(testArea < .02) = .02;
    testArea = truss.maxVol*testArea/dot(testArea, truss.len);
    
    FRed = truss.F(truss.dof);
    
    K = zeros(numel(truss.dof), numel(truss.dof));
    for e = 1:numel(truss.el)/2
        K = truss.K0El(truss.dof, truss.dof, e)*testArea(e) + K;
    end
    
    truss.U = zeros(numel(truss.n), 1);
    truss.U(truss.dof) = inv(K)*FRed;
    
    if truss.iter == 0
        truss.C0 = truss.F'*truss.U;
        truss.C  = truss.F'*truss.U;
    else
        C = truss.F'*truss.U;
        if C < truss.C
            truss.C = C;
            truss.area = testArea;
        end
    end
    truss.iter = truss.iter + 1;
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

set(truss.info, 'string', ...
    { sprintf('compliance0: %6.3f iterations:%d', truss.C0, truss.iter), ...
      sprintf('complianceI: %6.3f', truss.C)});
end
