% truss-opt with fmincon (3d) [Silas Henderson]
clear; clc; close all; 
global nodes elements dofs test_areas element_count node_count;
global F K0 lines flines lengths iter volume_max; iter = 0;
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
F = zeros(3*node_count, 1); 
for fs = dofs
    F(fs) = 30*(1 - 2*rand);
end
E          = 1000;
test_areas = rand(1, element_count);
area_min   = zeros(1, element_count);
area_max   = 10*ones(1, element_count);
volume_max = 10;
K0 = zeros(numel(nodes), numel(nodes), numel(test_areas));
% ------------------------------ Gui ---------------------------------- %
fig = figure('color', [.2, .2, .2]);
%'units', 'normalized', ...
%             'position', [.2 .2 .6 .6]);
ax  = axes(  'XGrid', 'on', 'XLim', [-.4, 1.4], ...
             'YGrid', 'on', 'YLim', [-.4, 1.4], ...  
             'ZGrid', 'on', 'ZLim', [-.4, 1.4], ...
             'color', [.95 .97 1], ...
             'gridcolor', [0 0 0], ...
             'units', 'normalized', 'outerposition', ...
             [-.05, -.05, 1.1, 1]);
view(140, 18);  
         
line([nodes(1, 1), nodes(2, 1), nodes(3, 1), nodes(4, 1)], ...
     [nodes(1, 2), nodes(2, 2), nodes(3, 2), nodes(4, 2)], ...
     [nodes(1, 3), nodes(2, 3), nodes(3, 3), nodes(4, 3)], ...
     'linestyle', 'none', 'marker', 'o', 'color', 'black', ...
     'markersize', 10, 'markerfacecolor', [.4 .6 .4]);
 
flines = line('linewidth', 2, 'color', [.8, .2, .2]);
info = annotation('textbox', 'position', [0, .9, 1, .1], ...
    'string', ' random-truss-opt with fmincon', 'fontsize', 16, ...
    'backgroundcolor', [.3, .2, .1], 'color', [.8 .8 1]);      
          
newKey = uicontrol(      'style','pushbutton', ...
           'units',     'normal',  'position', [.81 .91 .18 .08], ...
          'string',        'new',  'callback',              @opt, ...
 'backgroundcolor', [.4, .5, .3],  'fontname',           'arial', ...
 'foregroundcolor',    [1 1 1],  'fontsize',                 14);   
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
compliance(test_areas);
trussPlot(test_areas);
% ----------------------- FMinCon ---------------------------------- %
function opt(~, ~)
    global element_count node_count lengths F;
    dofs = 13:24;
    F = zeros(3*node_count, 1); 
    for fs = dofs
        F(fs) = 30*(1 - 2*rand);
    end
    test_areas = rand(1, element_count);
    area_min   = zeros(1, element_count);
    area_max   = 10*ones(1, element_count);
    volume_max = 10;
    fmincon(@compliance, test_areas, lengths, ...
        volume_max, [], [], area_min, area_max, [], optimoptions( ...
        'fmincon', 'outputfcn', @trussPlot, 'display',   'iter'));
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
    stop = false;  global nodes elements U lines flines F;
       
    nodeU = reshape(nodes', numel(nodes), 1) + U;
    for el = 1:numel(elements)/2
        node_a = elements(el,1);
        node_b = elements(el,2);      
        
        ua = test_areas(el)/max(test_areas);
        
        if test_areas(el) > 0
            set(lines(el), ...
                'XData', [nodeU(3*node_a - 2), nodeU(3*node_b - 2)], ...
                'YData', [nodeU(3*node_a - 1), nodeU(3*node_b - 1)], ...
                'ZData', [nodeU(3*node_a),     nodeU(3*node_b)    ], ...
                'linewidth', ua*7, ...
                'color', [.9 .9 1] - ua*[.9 .9 1]);
        end
    end
    
    X = []; Y = []; Z = [];
    for f = 1:numel(F)/3
        
        fx = F(3*f - 2);   fy = F(3*f - 1);   fz = F(3*f);
        fLen = sqrt(fx*fx + fy*fy + fz*fz);
        
        if fLen ~= 0
            ux = nodes(f, 1) + U(3*f - 2);  fx = fx/80;
            uy = nodes(f, 2) + U(3*f - 1);  fy = fy/80;
            uz = nodes(f, 3) + U(3*f);      fz = fz/80;
            
            X = [X, ux, ux + fx, nan]; 
            Y = [Y, uy, uy + fy, nan]; 
            Z = [Z, uz, uz + fz, nan];
            
           for i = 1:200
                rad = cross([fx;fy;fz], ones(3, 1) - 2*rand(3, 1)); 
                rad = rad/20; 
                
                X = [X, ux + fx, ux + fx*.8 + rad(1), nan];
                Y = [Y, uy + fy, uy + fy*.8 + rad(2), nan];
                Z = [Z, uz + fz, uz + fz*.8 + rad(3), nan];
           end
        end 
    end
    set(flines, 'XData', X, 'YData', Y, 'ZData', Z);
    drawnow; pause(.05);
end