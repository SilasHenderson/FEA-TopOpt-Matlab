% --- Fast-Plot (2d) -- (Silas Henderson) ------------------------------ %

clear; clc; close all; global canvas;
canvas.fig = figure('color', [.2, .2, .2],     'menubar', 'none', ...
                     'name',  'fast-plot', 'numbertitle',  'off', ...
              'keypressfcn',   @keyboard); 

canvas.ax  = axes( 'units', 'normal', 'position', [.01 .01 .98 .88], ...
           'XAxisLocation', 'origin',     'XLim', [-1, 1], 'XGrid', 'on', ...
           'YAxisLocation', 'origin',     'YLim', [-1, 1], 'YGrid', 'on', ...
              'TickLength',    [0 0], 'TickLabelInterpreter', 'latex');
  
canvas.line = line(0, 0, 'parent', canvas.ax, 'linewidth', 3);     

canvas.inputBox = uicontrol(...
        'style',   'edit', 'fontsize', 16, ...
        'units', 'normal', 'position', [0.21 0.9 0.58 0.08]);

canvas.topKey = uicontrol( 'style', 'pushbutton', ...
     'string',   'plot',   'callback',    @plotInput, ...
 'fontweight',   'bold',   'backgroundcolor', [.2 .3 .4], ...
   'fontname',  'arial',   'foregroundcolor', [1 1 1], ...
   'fontsize',       12,   'units', 'normal', ...
   'position', [.01 .9 .18 .08]);

canvas.tDisp = annotation('textbox',    'string',        't = 0', ...
                  'units', 'normal',  'position', [.8 .9 .2 .08], ...
                  'color',  [1 1 1], 'edgecolor',     [.2 .2 .2], ...
            'interpreter',   'none',  'fontsize',            14);

function keyboard(~, event)
    xlim = get(gca, 'XLim');  dx = xlim(2) - xlim(1);
    ylim = get(gca, 'YLim');  dy = ylim(2) - ylim(1);
    switch event.Key
        case 'leftarrow',  set(gca, 'XLim', xlim - dx/7*[1 1]);
        case 'rightarrow', set(gca, 'XLim', xlim + dx/7*[1 1]);
        case 'downarrow',  set(gca, 'YLim', ylim - dy/7*[1 1]);
        case 'uparrow',    set(gca, 'YLim', ylim + dy/7*[1 1]);  
        case 'shift'
            set(gca, 'XLim', xlim + dx/4*[-1, 1]);
            set(gca, 'YLim', ylim + dx/4*[-1, 1]);
        case 'return'
            set(gca, 'XLim', xlim + dx/4*[ 1,-1]);
            set(gca, 'YLim', ylim + dx/4*[ 1,-1]);
    end
    pause(.001);
end   
   
function plotInput(~, ~)
    global canvas;
    text = canvas.inputBox.String;
    fun  = str2func(strcat('@(x, t)', text)); 
    xlim = get(gca, 'XLim');
    dx   = (xlim(2) - xlim(1))/100;   
    
    for t = 0:.05:10
        xI = 1;
        for x = xlim(1):dx:xlim(2)
            y(xI) = fun(x, t);  
            xI = xI + 1;
        end
        set(canvas.line,  'XData', xlim(1):dx:xlim(2), 'YData', y);
        set(canvas.tDisp,'string', sprintf('t = %4.2f', t));
        pause(.01);
    end
end