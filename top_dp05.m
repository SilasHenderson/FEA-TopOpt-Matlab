% dp_05 (Silas Henderson)

clear; clc; close all;  global canvas;
% ----------------------------- Options ------------------------------- %
nelx =   20;     volfrac       = .5;    
nely =   40;     iterations    = 50;  
p    =    3;     filter_radius =  2; 
minp = .001;     E             =  1;

fixed_count = 0;
for i = 1:2*(nely+1):2*(nely+1)*(nelx+1)
    fixed_count = fixed_count + 1;  dof_fixed(fixed_count) = i;     
    fixed_count = fixed_count + 1;  dof_fixed(fixed_count) = i + 1;
end

f_global_1 = sparse(2*(nely+1)*(nelx+1),1); f_global_1(2*(nely+1)-1) = 1;
f_global_2 = sparse(2*(nely+1)*(nelx+1),1); f_global_2(2*(nely+1)*nelx-1,1) = -1;

% ----------------------------- Setup -------------------------------- %
dof_global = 1:2*(nely+1)*(nelx+1);
dof_free   = setdiff(dof_global,dof_fixed);
x          = volfrac.*ones(nely*nelx,1);

[el_nodes,el_dof, node_coords] = find_el_nodes(nelx,nely);
set_canvas(node_coords, el_nodes);
k0 = k0make(E); 

% ----------------------------- Solve -------------------------------- %
for i = 1:50
    
    [c_1, dc_1,u_global_1] = solve(x,p,k0,f_global_1,dof_free,el_dof);
    [c_2, dc_2,u_global_2] = solve(x,p,k0,f_global_2,dof_free,el_dof);
    dc = (dc_1+dc_2)/2;
    x  = oc_update(x,volfrac,dc);  
    
    [c, dc,u_global] = solve(x,p,k0,f_global_1+f_global_2,dof_free,el_dof);
    set(canvas.el_patch, 'facevertexcdata', 1-[x,x,x]);
    plot_u(u_global, node_coords,c, x);
    drawnow;
end                                                 

% --------------------------- Functions-------------------------------%
function k0 = k0make(E)
    nu = .3; syms xi eta;                          
    N1234 = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);(1+xi)*(1+eta);(1-xi)*(1+eta)];         
    Nxi   = diff(N1234,xi);  Neta  = diff(N1234,eta);
    B = [Nxi(1) 0     Nxi(2)    0    Nxi(3)  0      Nxi(4)    0    ;
         0      Neta(1)  0     Neta(2)  0     Neta(3)   0    Neta(4);
        Neta(1) Nxi(1) Neta(2) Nxi(2) Neta(3) Nxi(3) Neta(4) Nxi(4)];
    Cmod = 1/(1-nu^2).*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];  J = 4;
    k0 = int(int(B'*Cmod*B,xi,-1,1),eta,-1,1); 
    k0 = J*E*double(k0);     
end

function [el_nodes, el_dof, node_coords] = find_el_nodes(nelx,nely)
    el_count = 0;
    for ex = 1:nelx 
        for ey = 1:nely
            el_count = el_count + 1;
            upleft   = (nely+1)*(ex-1)+ey+1; upright  = (nely+1)*ex+ey+1;
            lowleft  = (nely+1)*(ex-1)+ey;   lowright = (nely+1)*ex+ey;   
            el_nodes(el_count,:) = [lowleft, lowright, upright, upleft];
            el_dof(  el_count,:) = ...
                [2*lowleft-1, 2*lowleft, 2*lowright-1, 2*lowright, ...
                 2*upright-1, 2*upright, 2*upleft - 1, 2*upleft];
        end
    end
    
    node_count = 0;
    for nx = 0:nelx
        for ny = 0:nely
            node_count = node_count + 1;
            node_coords(node_count,:) = [nx, ny];
        end
    end
end

function [c, dc, u_global] = solve(x,p,k0,f_global,dof_free,el_dof)
    
    k_global = sparse(numel(f_global), numel(f_global));
    u_global = zeros (numel(f_global),1); 
    
    for el = 1:numel(x)
        edof = el_dof(el,:);
        k_global(edof,edof) = k_global(edof,edof) + x(el)^p*k0;
    end
    
    u_global(dof_free,:) = k_global(dof_free, dof_free)\f_global(dof_free,:);
    c = 0.;     
    
    for el = 1:numel(x)
        u_el   = u_global(el_dof(el,:));
        c      = c + x(el,1)^p*u_el'*k0*u_el; 
        dc(el) = -p*x(el,1)^(p-1)*u_el'*k0*u_el;
    end
end

function xnew = oc_update(x,volfrac,dc)  
    
    lambdamin = 0; 
    lambdamax = 100000; 
    move = 0.2;
    
    while (lambdamax-lambdamin > 1e-4)
        lambda = 0.5*(lambdamax+lambdamin);
        xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc'./lambda)))));
        if sum(xnew) - volfrac*numel(x) > 0, lambdamin = lambda;
        else, lambdamax = lambda;
        end
    end
end

% ---------------------- Drawing Functions ------------------------- %
function set_canvas(node_coords, el_nodes)
    global canvas;  el_num = numel(el_nodes)/4;
    x_max = max(node_coords(:,1)) + 1;
    y_max = max(node_coords(:,2)) + 1; xy_max = max(x_max, y_max);
    canvas.fig  = figure('menubar', 'none',  'color', [ .2, .3, .2]);
    canvas.axis = axes( ...
     'units',    'normal',  'TickLen', [0,0], ...
      'box',         'on', 'position', [0.01,.01,.98,.88], ...
     'XLim', [-1, xy_max],    'XGrid', 'on', 'XTick', 0: xy_max, ...
     'YLim', [-1, xy_max],    'YGrid', 'on', 'YTick', 0: xy_max);

    canvas.info = annotation('textbox',   'string',   'compliance', ...
                 'color',    [1, 1, 1],    'units',       'normal',...
       'backgroundcolor', [.2, .2, .2], 'position', [0, .9, 1, .1]);
  
    canvas.el_patch = patch('vertices', node_coords,     'faces', el_nodes, ...
                     'facevertexcdata', ones(el_num,1), ...
                     'facecolor',   'flat', 'edgecolor', 'none', ...
                   'facevertexalphadata', 1);
    canvas.node_dots = line('XData', [], 'YData',     [], ...
     'linestyle', 'none',  'marker','.', 'markersize', 5, 'color', [.4, .5, .8]);
end

function plot_u(u_global, node_coords, c, x)
     global canvas; 
      X = []; Y = [];
     for n = 1:length(node_coords)
         x = node_coords(n,1) + u_global(2*n-1, 1); 
         y = node_coords(n,2) + u_global(2*n, 1);
         X = [X,x];
         Y = [Y,y];
     end
     set(canvas.node_dots,'XData', X, 'YData', Y);    
     set(canvas.info, 'string', ... 
        {sprintf('compliance = %4.4f',c);
         sprintf('volume = %5.4f', sum(x))});  
end