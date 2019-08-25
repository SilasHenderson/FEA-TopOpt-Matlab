% --- Fast-Plot(3d) -- (Silas Henderson) ----------------------------- %

clear; clc; close all; global canvas;
canvas.fig = figure('color', [.8, .8, .8],     'menubar', 'none', ...
                     'name', 'fast-plot3', 'numbertitle',  'off', ...
              'keypressfcn',   @keyboard); 
canvas.el = .5;  canvas.az = 2.7;
                canvas.camDist = 20;

canvas.ax = axes( 'units','normal','position', [.01 .01 .98 .88], ...
                   'XLim', [-1, 1],   'XGrid', 'on', 'box', 'on', ...
                   'YLim', [-1, 1],   'YGrid', 'on', ...
                   'ZLim', [-1, 1],   'ZGrid', 'on', ...
   'TickLabelInterpreter', 'latex','TickLength', [0 0], ...
        'DataAspectRatio', [1 1 1],'projection', 'perspective', ...
     'PlotBoxAspectRatio', [1 1 1],'CameraViewAngle', 30);
         
canvas.line = line(0, 0, 'parent', canvas.ax, 'linewidth', 1);     

canvas.inputBox = uicontrol( ...
        'style',   'edit', 'fontsize', 16, 'string', 'sin(x*y + t)/(cos(t/2) + 2)', ...
        'units', 'normal', 'position', [0.21 0.9 0.58 0.08]);

canvas.tDisp = annotation('textbox',    'string',        't = 0', ...
                'units',   'normal',  'position', [.8 .9 .2 .08], ...
                'color', [.2 .3 .2], 'linestyle',         'none', ...
          'interpreter',     'none',  'fontsize',             14);   
    
canvas.topKey = uicontrol( 'style', 'pushbutton', ...
     'string',   'plot',   'callback',    @plotInput, ...
 'fontweight',   'bold',   'backgroundcolor', [.2 .3 .4], ...
   'fontname',  'arial',   'foregroundcolor', [1 1 1], ...
   'fontsize',       12,   'units', 'normal', ...
   'position', [.01 .9 .18 .08]);

camSet;

function camSet   
    global canvas;
    set(gca,'CameraPosition', ...
       canvas.camDist*[cos(canvas.el)*cos(canvas.az), ...
                       cos(canvas.el)*sin(canvas.az), ...
                       sin(canvas.az)], ...
               'CameraTarget', [0 0 0], ...
      'XLim', canvas.camDist*[-1, 1]/5, ...
      'YLim', canvas.camDist*[-1, 1]/5, ...
      'ZLim', canvas.camDist*[-1, 1]/5);      
end

function keyboard(~, event)
    global canvas;
    switch event.Key
        case 'leftarrow',  canvas.az = canvas.az + .1;
        case 'rightarrow', canvas.az = canvas.az - .1;
        case 'downarrow',  canvas.el = canvas.el + .1;
        case 'uparrow',    canvas.el = canvas.el - .1;   
        case 'shift',      canvas.camDist = canvas.camDist*3/2;
        case 'return',     canvas.camDist = canvas.camDist*2/3;
    end
    camSet;
    pause(.001);
end   
   
function plotInput(~, ~)
    global canvas;
    text = canvas.inputBox.String;
    fun  = str2func(strcat('@(x, y, t)', text)); 
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    dx   = (xlim(2) - xlim(1))/100;   
    dy   = (ylim(2) - ylim(1))/100;
    for t = 0:.05:10
        I = 1;
        for x = xlim(1):dx:xlim(2)
            for y = ylim(1):dy:ylim(2)
                X(I) = x;
                Y(I) = y;
                Z(I) = fun(x, y, t);  
                I = I + 1;
            end
            X(I) = NaN; Y(I) = NaN; Z(I) = NaN;  I = I + 1;
        end
        for y = ylim(1):dy:ylim(2)
            for x = xlim(1):dx:xlim(2)
                X(I) = x;
                Y(I) = y;
                Z(I) = fun(x, y, t);  
                I = I + 1;
            end
            X(I) = NaN; Y(I) = NaN; Z(I) = NaN;  I = I + 1;
        end
        set(canvas.line,  'XData', X, 'YData', Y, 'ZData', Z);
        set(canvas.tDisp,'string', sprintf('t = %4.2f', t));
        pause(.01);
    end
end