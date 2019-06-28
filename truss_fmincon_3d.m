% min-truss-opt (3d) [Silas Henderson]

clear; clc; close all; 
global nodes elements dofs;
global F K0 lines iter; iter = 0;
% ------------------------- Parameters ------------------------------- %
nodes = [0, 0, 0; 0, 0, 1; 1, 0, 1; 1, 0, 0;
         0, 1, 0; 0, 1, 1; 1, 1, 1; 1, 1, 0;];
      
elements = [1, 5; 1, 6; 1, 8;
            2, 6; 2, 5; 2, 7;
            3, 7; 3, 6; 3, 8;
            4, 8; 4, 7; 4, 5;
            5, 6; 6, 7; 7, 8; 8, 5; 6, 8; 7, 5];
          
node_count    = numel(nodes)/3;
element_count = numel(elements)/2;

dofs = 13:24;

F = zeros(3*node_count, 1); F([end, end - 3, end - 6, end - 9]) = -3;

E          = 1000;
test_areas = rand(1, element_count);
area_min   = zeros(1, element_count);
area_max   = 10*ones(1, element_count);
volume_max = 10;

K0 = zeros(numel(nodes), numel(nodes), numel(test_areas));

fig = figure('color', [.2, .2, .2]);
ax  = axes(  'XGrid', 'on', 'XLim', [-.2, 1.2], ...
             'YGrid', 'on', 'YLim', [-.2, 1.2], ...  
             'ZGrid', 'on', 'ZLim', [-.5, 1.2]); view(140, 18);

line([nodes(1, 1), nodes(2, 1), nodes(3, 1), nodes(4, 1)], ...
     [nodes(1, 2), nodes(2, 2), nodes(3, 2), nodes(4, 2)], ...
     [nodes(1, 3), nodes(2, 3), nodes(3, 3), nodes(4, 3)], ...
     'linestyle', 'none', 'marker', 'o', 'color', 'black', ...
     'markersize', 10, 'markerfacecolor', [.2 .8 .2]); fArrow;
 
% -------------------------- Assembly --------------------------------- %
for el = 1:element_count
    lines(el) = line(0, 0, 'parent', ax);
    
    node_a = elements(el, 1);                    
    node_b = elements(el, 2);
    dx = nodes(node_b,1) - nodes(node_a,1);
    dy = nodes(node_b,2) - nodes(node_a,2);  
    dz = nodes(node_b,3) - nodes(node_a,3);
    
    lengths(el) = norm([dx, dy, dz]);
      
    element_dofs = [3*node_a - 2, 3*node_a - 1, 3*node_a, ...
                    3*node_b - 2, 3*node_b - 1, 3*node_b];
    
    k9 = [dx*dx, dx*dy, dx*dz;
          dy*dx, dy*dy, dy*dz;
          dz*dx, dz*dy, dz*dz]*E/lengths(el)^4;
    
    K0(:, :, el) = zeros(numel(nodes), numel(nodes)); 
    K0(element_dofs, element_dofs, el) = [k9,-k9;-k9,k9];
end

% ----------------------- FMinCon ---------------------------------- %
options = optimoptions('fmincon', 'outputfcn', @trussPlot);

area_optimal = fmincon(@compliance, test_areas, lengths, volume_max, ...
               [], [], area_min, area_max, [], options);
  
% -------------------------- arrow  ------------------------------- %
function fArrow
    global F nodes; X = []; Y = []; Z = [];
    for f = 1:numel(F)/3
        fx = F(3*f - 2);   fy = F(3*f - 1);   fz = F(3*f);
        fLen = sqrt(fx*fx + fy*fy + fz*fz);
        if fLen ~= 0
            nx = nodes(f, 1);   fx = fx/fLen/2;  X = [X, nx, nx + fx];
            ny = nodes(f, 2);   fy = fy/fLen/2;  Y = [Y, ny, ny + fy];
            nz = nodes(f, 3);   fz = fz/fLen/3;  Z = [Z, nz, nz + fz]; 
        
            for ang = 0:.1:2*pi
                X = [X, nx + fx, nx + fx + .04*cos(ang), nan];
                Y = [Y, ny + fy, ny + fy + .04*sin(ang), nan];
                Z = [Z, nz + fz, nz + fz + .1, nan];
            end
        end
    end
    line('XData', X, 'YData', Y, 'ZData', Z, ...
       'linewidth', 2, 'color', [.8, .2, .2]);
end          

% --------------------------- Compliance ------------------------------ %
function C = compliance(test_areas)
    global nodes elements dofs F K0 U;
    K = zeros(numel(nodes), numel(nodes));  
    for el = 1:numel(elements)/2
        K = K + test_areas(el)*K0(:, :, el);
    end
    U = zeros(numel(nodes), 1);
    U(dofs) = K(dofs, dofs)\F(dofs);  
    C = F'*U;
end

% ----------------------------- Plot ---------------------------------- %
function stop = trussPlot(test_areas, ~, ~)
    stop = false;  global nodes elements U lines;
       
    nodeU = reshape(nodes', numel(nodes), 1) + U;
    for el = 1:numel(elements)/2
        node_a = elements(el,1);
        node_b = elements(el,2);      
        
        set(lines(el), ...
        'XData', [nodeU(3*node_a - 2), nodeU(3*node_b - 2)], ...
        'YData', [nodeU(3*node_a - 1), nodeU(3*node_b - 1)], ...
        'ZData', [nodeU(3*node_a),     nodeU(3*node_b)    ], ...
        'linewidth', 5*test_areas(el)/max(test_areas));
    end; drawnow;
end