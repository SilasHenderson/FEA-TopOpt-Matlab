% M.L. Paint {Silas Henderson 2019}
clear; clc; close all;
% ------------------------- Setup ------------------------------------ %
global mouseClick;   mouseClick = 0;                        % Mouse
global mouseDown;    mouseDown  = 0;
Colors = [.8 .6 .6  1 .5  0  1  0  0;                       
          .6 .8 .6  1 .5  0  0  1  0;
          .6 .6 .8  1 .5  0  0  0  1]';
[brush,  color,   bSize, ...
 panel, labels, mlpLogo] = startCanvas(Colors);             
% ----------------------------- Loop ----------------------------------- %
global programOn; programOn = 1;
while programOn == 1
    figPos = get(gcf, 'Position');                          % Position
    cPos   = get(  0, 'PointerLocation');        
    x      = (cPos(1) - figPos(1))/(figPos(3));        
    y      = (cPos(2) - figPos(2))/(figPos(4)); 
    set(brush, 'MarkerSize',  bSize,  'XData', x,...
          'MarkerFaceColor',  color,  'YData', y);    
      
    if mouseClick == 1                                      % Mouse Click
        Line = animatedline('LineWidth', bSize, 'Color', color);
        uistack(brush, 'top');
        if x > .9 && y < .9                                 % -- Side Bar                             
            color = Colors(ceil(10*y),:);
            uistack(  panel, 'top');  uistack(labels, 'top');
            uistack(mlpLogo, 'top');  uistack( brush, 'top');
        end
        if y > .9                                           % -- Top Bar
            uistack(  panel, 'top');  uistack(labels, 'top');
            uistack(mlpLogo, 'top');  uistack( brush, 'top');
            if x < .1                                       % -- Save
                imageTensor = getfield(getframe(gcf),'cdata');
                imwrite(imageTensor,uiputfile({'*.jpg'}))
                close all; return;
            elseif x > .15 && x < .25                       % -- New
                [brush,  color,   bSize, ...
                 panel, labels, mlpLogo] = startCanvas(Colors);
                 pause(.1); continue;
            end
        end 
        mouseClick = 0;                                     
    end
    if mouseDown == 1                                       % Mouse Down
        if y > .9
            if x > .3 && x < .45
                bSize = bSize + 1;
            elseif x > .45 && x < .6
                bSize = max(bSize - 1, 4); 
            end
            set(Line, 'LineWidth', bSize);
        end
        addpoints(Line, x, y);                           
    end
    pause(.01);                                         
end
delete(gcf);
% --------------------------- Functions ------------------------------ %
function [brush,  color,   bSize, ...
          panel, labels, mlpLogo] = startCanvas(Colors)
clf;                                                         % Vars
bSize  = 8;     color  = [0,0,0];              
set(  0, 'units', 'normalized');
set(gcf, 'WindowButtonDownFcn',   @mouseClickFcn,   'Color',      [0 0 0], ...
           'WindowButtonUpFcn', @mouseClickUpFcn, 'MenuBar',       'None', ...
             'CloseRequestFcn',     @closeReqFcn, 'pointer',     'custom', ...
           'PointerShapeCData',       NaN(16,16),   'units', 'normalized');
set(gca, 'Units',  'Normal',    'XLim', [0 1], ...          
      'Position', [0 0 1 1],    'YLim', [0 1], ...     
 'PickableParts',    'None', 'TickLen', [0 0], 'box', 'on')
for i = 1:9                                                  % Side Bar
    panel(i) = boxPatch(.95, .1*i  -.05, ...
                        .05,        .05, Colors(i,:),'black');
end
panel(10) = boxPatch(   .5,  .95,   .5,  .05, .15*ones(1,3), zeros(1,3));   
panel(11) = boxPatch(  .45, .005,  .45, .005,  .4*ones(1,3), zeros(1,3));
panel(12) = boxPatch(  .45, .895,  .45, .005,  .4*ones(1,3), zeros(1,3)); 
panel(13) = boxPatch( .005,  .45, .005,  .45,  .4*ones(1,3), zeros(1,3)); 
panel(14) = boxPatch( .895,  .45, .005,  .45,  .4*ones(1,3), zeros(1,3));
labels(1) = boxPatch(.05, .97, .035, .01, [1 1 1], [1 1 1]); % -- Save
labels(2) = boxPatch(.05, .93, .035, .01, [1 1 1], [1 1 1]);
labels(3) = boxPatch( .2, .95, .035, .03, [1 1 1], [1 1 1]); % -- New
labels(4) = boxPatch(.35, .95,  .03, .01, [1 1 1], [1 1 1]); % -- Size -    
labels(5) = boxPatch(.35, .95,  .01, .03, [1 1 1], [1 1 1]); 
labels(6) = boxPatch( .5, .95,  .03, .01, [1 1 1], [1 1 1]); % -- Size +
mlpLogo = text(.6, .95,  'M.L. Paint',    'color', 'white', ...
           'fontUnits',  'normalized', 'fontsize',     .08);
brush = line('linestyle',    'none', ...                    % Brush
            'MarkerSize', bSize, 'marker', 'o', ...
       'MarkerEdgeColor', 'white');
end
function mouseClickFcn(~,~)
global mouseClick; mouseClick = 1;
global mouseDown;  mouseDown  = 1;
end
function mouseClickUpFcn(~,~)
global mouseDown;  mouseDown = 0;
end
   
function closeReqFcn(~,~)
global programOn; programOn = 0;
delete(gcf);
end
function box = boxPatch(x, y, w, h, fColor, eColor)
box = patch('Vertices', [x - w, x - w, x + w, x + w;
                         y - h, y + h, y + h, y - h]', ...
               'Faces', [    1,     2,     3,     4] , ...
           'FaceColor', fColor, 'EdgeColor', eColor);     
end                                    