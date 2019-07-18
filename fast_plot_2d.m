% --- fast-plot (2d), Silas Henderson IUPUI ---
clear; clc; close all;
global canvas;
canvas.t  =   0;
canvas.T  =  10;
canvas.dt = .05;
canvas.optionsOn = 0;
canvas.fig = figure('color', [.2, .2, .2],     'menubar', 'none', ...
                     'name',  'fast-plot', 'numbertitle',  'off', ...
              'keypressfcn',    @keyboard); 
 
canvas.altAx   = axes('units', 'normal', 'XLim', [-.01, .98], ...
        'position', [.01 .01, .98, .88], 'YLim', [ .01, .99], ...
        'color', [.7 .8 .9], 'visible',   'off', 'TickLen', [0 0]);      
        
          
canvas.ax  = axes('units', 'normal', 'position', [.01 .01, .98, .88], ...
  'XAxisLocation', 'origin', 'XLim', [-1, 1],'XGrid', 'on', ...
  'YAxisLocation', 'origin', 'YLim', [-1, 1],'YGrid', 'on', ...
  'TickLength', [0 0], 'TickLabelInterpreter', 'latex');
 
canvas.line      = line(0, 0, 'parent', canvas.ax, 'linewidth', 3);     
canvas.inputBox  = uicontrol(...
       'style',   'edit', 'fontsize', 16, ...
        'units','normal', 'position', [0.25 0.9 0.58 0.08]);
for i = 1:6   
    canvas.optionsBox(i) = uicontrol('style', 'edit', 'fontsize', 16, ...
        'units', 'normal', 'position', [.72, i/10 - .04, .24, .08], ...
        'visible', 'off');
    
    canvas.lines(i) = line(0, 0, 'parent', canvas.ax, ...
                  'linewidth', 2, 'color', rand(1,3));     
end
infoText = {'Welcome to fast-plot!';  "";
                      "How to use:";  "";
    "Write an expression f(x,t) in the upper textbox";        "";
    "Or, several expressions f(x) in the boxes to the right"; "";                    
    "Press *Enter/Shift* to Zoom";
    "Press *Keyboard Arrows* to Pan"};
canvas.optionsBox(7) = annotation('textbox', 'units', 'normal', ...
    'position', [.05 .06 .6 .78], 'String', infoText, 'Color', [.3 .1 .1], ...
    'visible', 'off', 'fontsize', 14, 'backgroundcolor', [1 1 1], ...
    'facealpha', .7, 'interpreter', 'latex');
canvas.topKey = uicontrol(  'style',        'pushbutton', ...
      'units', 'normal', 'position',    [.01 .9 .13 .08], ...
     'string','options', 'callback',         @optionsKey, ...
 'fontweight',   'bold', 'backgroundcolor',   [.4 .4 .4], ...
   'fontname',  'arial', 'foregroundcolor',      [1 1 1], ...
  'fontsize',       12);
canvas.topKey2 = uicontrol( 'style',        'pushbutton', ...
      'units', 'normal', 'position',    [.15 .9 .09 .08], ...
     'string',   'plot', 'callback',          @plotInput, ...
 'fontweight',   'bold', 'backgroundcolor',   [.2 .3 .4], ...
   'fontname',  'arial', 'foregroundcolor',      [1 1 1], ...
  'fontsize',       12);
canvas.tDisp = annotation('textbox',    'string',            '', ...
                  'units', 'normal',  'position', [.85 .9 .15 .08], ...
                  'color',  [1 1 1], 'edgecolor',       [.2 .2 .2], ...
            'interpreter',   'none',  'fontsize',              14);
function optionsKey(~, ~)
    global canvas;
    if canvas.optionsOn == 0
        canvas.optionsOn = 1;
        set(canvas.ax, 'color', [.7 .8 .9]);
        for i = 1:7
            set(canvas.optionsBox(i), 'visible', 'on');
        end
    else
        canvas.optionsOn = 0;
        set(canvas.ax, 'color', [1 1 1]);
        for i = 1:7
            set(canvas.optionsBox(i), 'visible', 'off');
        end
        staticPlot;
    end
end
        
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
   
function staticPlot(~, ~)
    global canvas;
    xlim = get(gca, 'XLim');
    dx   = (xlim(2) - xlim(1))/100;   
    for i = 1:6
        xI   = 1;
        text = canvas.optionsBox(i).String;
        set(canvas.lines(i), 'XData', [], 'YData', []);
        if isempty(text) == 0
            fun  = str2func(strcat('@(x)', text));    
            for x = xlim(1):dx:xlim(2)
                y(xI) = fun(x);  
                xI = xI + 1;
            end
            set(canvas.lines(i), 'XData', xlim(1):dx:xlim(2), 'YData', y);
        end
    end
end
function plotInput(~, ~)
    global canvas;
    if canvas.optionsOn == 0
        text  = canvas.inputBox.String;
        xlim  = get(gca, 'XLim');
        dx    = (xlim(2) - xlim(1))/100;  
        xI    = 1;
        if isempty(text) == 0
            if any(text == 't')
                fun   = str2func(strcat('@(x, t)', text)); 
             
                T  = canvas.T;
                dt = canvas.dt;
                t  = canvas.t;
    
                while t < T
                    t = t + dt;
                    xI = 1;
                    for x = xlim(1):dx:xlim(2)
                        y(xI) = fun(x, t);  
                        xI = xI + 1;
                    end
                    set(canvas.line,   'XData', xlim(1):dx:xlim(2), 'YData', y);
                    set(canvas.tDisp, 'string', sprintf('t=%4.1f', t));
                    pause(.001);
                end
            else
                fun   = str2func(strcat('@(x)', text)); 
                for x = xlim(1):dx:xlim(2)  
                    y(xI) = fun(x);
                    xI = xI + 1;
                end
                set(canvas.line,   'XData', xlim(1):dx:xlim(2), 'YData', y);
                set(canvas.tDisp, 'string', '');
            end
        end
    else
        staticPlot;
    end
end