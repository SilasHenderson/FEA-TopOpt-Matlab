% Ten-Bar Truss (lite) [Silas Henderson]

clear; clc; close all; 
global nodes elements dofs;
global F K0 C_too info lines iter; iter = 0;
% ------------------------- Parameters ------------------------------- %
nodes   = [0, 2; 2, 2; 4, 2;
           0, 0; 2, 0; 4, 0];
      
elements = [1, 2; 2, 3; 4, 5; 5, 6; 1, 5;
            2, 6; 2, 4; 5, 3; 2, 5; 3, 6];

node_count    = numel(nodes)/2;
element_count = numel(elements)/2;

dofs = 1:2*node_count;          dofs([1, 2, 7, 8]) = [];
F = zeros(2*node_count, 1);     F(12) = -5;     

E          = 1000;
test_areas = ones(1, element_count);
area_min   = zeros(1, element_count);
area_max   = 10*ones(1, element_count);
volume_max = 10;

K0 = zeros(numel(nodes), numel(nodes), numel(test_areas));

fig  = figure('color', [.2, .2, .2]);
ax   = axes('XGrid', 'on', 'XLim', [-.5, 4.5],   'units',  'normalized', ...
            'YGrid', 'on', 'YLim', [-.5, 2.5],'position', [0, 0, 1, .9]);                
info = annotation('textbox',   'string',  ' ',   'units', 'normal', ...
                 'position',  [0, .9,  1, .1],   'color', 'white');
             
% -------------------------- Assembly --------------------------------- %
for el = 1:element_count  
    lines(el) = line(0, 0, 'parent', ax);
    
    node_a = elements(el, 1);                    
    node_b = elements(el, 2);
    dx = nodes(node_b,1) - nodes(node_a,1);
    dy = nodes(node_b,2) - nodes(node_a,2);  
    lengths(el) = norm([dx, dy]);
       
    element_dofs = [2*node_a - 1, 2*node_a, 2*node_b - 1, 2*node_b];
    
    k4 = [dx*dx, dx*dy;
          dx*dy, dy*dy]*E/lengths(el)^3;
    
    K0(:, :, el) = zeros(numel(nodes), numel(nodes)); 
    K0(element_dofs, element_dofs, el) = [k4,-k4;-k4,k4];
end

% ----------------------- FMinCon ---------------------------------- %
options = optimoptions('fmincon', 'outputfcn', @trussPlot);

area_optimal = fmincon(@compliance, test_areas, lengths, volume_max, ...
               [], [], area_min, area_max, [], options);
    
% --------------------------- Compliance ------------------------------ %
function C = compliance(test_areas)
    global nodes elements dofs F K0 U C_too iter; iter = iter + 1;
    K = zeros(numel(nodes), numel(nodes));  
    for el = 1:numel(elements)/2
        K = K + test_areas(el)*K0(:, :, el);
    end
    U = zeros(numel(nodes), 1);
    U(dofs) = K(dofs, dofs)\F(dofs);  
    C = F'*U; C_too = C;
end

% ----------------------------- Plot ---------------------------------- %
function stop = trussPlot(test_areas, ~, ~)
    stop = false;  global nodes elements U lines info C_too iter;
    
    set(info, 'String', ...
        sprintf('Compliance: %8.4f \nEvals: %8.4f', C_too, iter));
    
    nodeU = reshape(nodes', numel(nodes), 1) + U;
    for el = 1:numel(elements)/2
        node_a = elements(el,1);
        node_b = elements(el,2);      
        
        set(lines(el), ...
             'XData', [nodeU(2*node_a - 1), nodeU(2*node_b - 1)], ...
             'YData', [nodeU(2*node_a),     nodeU(2*node_b)    ], ...
             'linewidth', 10*test_areas(el));
    end; drawnow;
end