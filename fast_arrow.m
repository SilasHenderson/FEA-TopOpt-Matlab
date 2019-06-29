clear; clc; close all;

figure('color', [.2, .2, .2]);
axes('XLim', [-2, 2], 'YLim', [-2, 2], 'ZLim', [-2, 2]);
view(140, 20); grid on;

tic;
while toc < 1
    arrow(0, 0, 0, 2 - 4*rand, 2 - 4*rand, 2 - 4*rand);
    drawnow;
end
   
function arrow(x0, y0, z0, x1, y1, z1)
    X = [x0, x1, nan]; Y = [y0, y1, nan]; Z = [z0, z1, nan];
    u = [x1 - x0; y1 - y0; z1 - z0];
    for i = 1:200
        uv = cross(u, ones(3, 1) - 2*rand(3, 1)); uv = uv/norm(uv)/6;    
        X = [X, x1, x1 - (x1 - x0)/7 + uv(1), nan];
        Y = [Y, y1, y1 - (y1 - y0)/7 + uv(2), nan];
        Z = [Z, z1, z1 - (z1 - z0)/7 + uv(3), nan];
    end
    line(X, Y, Z, 'linewidth', 2, 'color', [.8 .2 .2]);
end