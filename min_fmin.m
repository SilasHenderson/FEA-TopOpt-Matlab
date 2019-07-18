% 3-bar truss (version 2)

clear; clc; close all;

% ------------------- Write Parameters --------------------- %
E = 30e6;  f = 10000;

% node_labels:   n1x n2y n2x n2y n3x n3y n3x n4y
% dof_labels:      1   2   3   4   5   6   7   8
node_coords   = [  0,  0;  0, 10; 10, 10; 10,  0];                      
node_forces   = [  0; -f;  0;  0;  0;  0;  0;  0]; 
dof           = [  1, 2];                 

% element_labels:    el1, el2, el3 
element_nodes     = [1,2; 1,3; 1,4]; 
element_lengths   = [nan, nan, nan];
element_areas_0   = [  2,   2,   2];
element_areas_min = [  0,   0,   0];
element_areas_max = [ 10,  10,  10];

% constraint:
volume_max = 10;

% k0 element matrices, stacked in a 'k0_el_tensor'
% ..............................                                
% ...........____el3__..........  k1, k2, k3 = ...
% ________ _|__el2__  |_________  k0 matrices for each element
% _ _    _|__el1__  | |    
%  |    |         | |_|   /           
% (1)   | k0_el   |_|   (3)      example: get k0_el with ...
% _|_   |_________|     /                         
%      /---(2)---/

k0_el_tensor = zeros(numel(node_coords), ...       % |k0..| <-- k0 matrices
                     numel(node_coords), ...       % |....|,
                     numel(element_areas));        % -------> [el1,el2,el3]

% ------------------ Assemble K0 Tensor --------------------- %           
for el = 1:length(element_nodes)           
    
    node_a = element_nodes(el, 1);        
    node_b = element_nodes(el, 2);      
    
    node_a_x = node_coords(node_a, 1);    
    node_a_y = node_coords(node_a, 2);
    node_b_x = node_coords(node_b, 1);    
    node_b_y = node_coords(node_b, 2);
    
    el_len_x = node_b_x - node_a_x;        
    el_len_y = node_b_y - node_a_y;        
        
    el_len = norm([el_len_x, el_len_y]);   
    
    c = el_len_x/el_len;                
    s = el_len_y/el_len;               
    
    k0_el_local = E/el_len*[ c*c,  c*s, -c*c, -c*s;        
                             c*s,  s*s, -c*s, -s*s;           
                            -c*c, -c*s,  c*c,  c*s;            
                            -c*s, -c*s,  c*s,  s*s];           
     
    element_dof = [2*node_a-1, 2*node_a, 2*node_b-1, 2*node_b]; 
    k0_el_tensor(element_dof, element_dof, el) = k0_el_local;         
end

options = optimoptions('fmincon', 'outputfcn', @trussPlot);
area_optimal = fmincon(@compliance, element_areas, element_lengths, ...
    volume_max, [], [], area_min, area_max, [], options);
    
% --------------------------- Compliance ------------------------------ %
function compliance = ComplianceFind(element_areas)

    global k0_el_tensor node_coordinates element_nodes forces dof
    K = zeros(numel(node_coordinates), numel(node_coordinates));  
    for el = 1:length(element_nodes)
        K = K + element_areas(el)*k0_el_tensor(:, :, el);
    end
    
    U = zeros(numel(nodes), 1);
    U(dofs) = K(dofs, dofs)\forces(dof);  
    compliance = F'*U;
end

