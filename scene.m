% MiniScape {Silas Henderson}
% Keyboard: {[up,down],[left,right],    [s],    [1,2,3]} -->
%           {   [move],      [turn], [view], [new game]}
clear; close all; 
global boxNum; boxNum = 50;
newGame;
% --------------------------- Functions ------------------------------- %
function newGame
clf;
global boxNum;                                                 % Vars
global camPos;      camPos     = [1,1,1]; 
global camAng;      camAng     =       0;  
global viewSwitch;  viewSwitch =       1;
set(gcf,'KeyPressFcn',     @controls,   'Color', [0 0 0], ...  % Fig/Ax
              'Units',  'Normalized', 'toolbar',  'none', ...
           'Position', [.1 .1 .6 .8], 'menubar',  'none');                                  
set(gca, 'projection', 'Perspective', 'Visible',   'Off', ...
    'CameraViewAngle',            50, ...
 'PlotBoxAspectRatio',       [1 1 1], ... 
    'DataAspectRatio',       [1 1 1]); 
for i = 1:boxNum                                                % Boxes 
    h = 5*rand + 5*(rand^4);
    boxPatch([100*rand; 100*rand; h], 3*rand,3*rand, h) 
end
goalPos = [20 + 80*rand, 20 + 80*rand, h];                      % Goal Box
goalBox = boxPatch(goalPos', 3*rand, 3*rand, h);
patch('Faces', [1 2 3 4], ...                                   % Floor
   'Vertices', 100*[1 0 0; 1 1 0; 0 1 0; 0 0 0], ...
  'FaceColor', [.5 .5 .5] + .2*rand(1,3)); 
                                            
camDots = line('Marker', '.', 'Markersize', 20);                % Markers
tic;   
while norm(camPos(1,1:2) - goalPos(1,1:2)) > 2                  
    
    camTar = camPos + 10*[cosd(camAng), sind(camAng), 0];       % Target
    
    set(goalBox,'FaceColor',[.5 + .5*cos(toc)*cos(toc*1.1),...  % Color
                             .5 + .5*cos(toc)*sin(toc*1.1),...
                             .5 + .5*sin(toc*1.1)]);
    if viewSwitch == 1                                          % View: FP
        set(gca,'CameraPosition', camPos,...        
                  'CameraTarget', camTar)
    else
        set(gca,'CameraPosition', [50,50,100], ...              % View: Arial
                  'CameraTarget', [50 50   0])
        set(camDots,     'XData', [camPos(1) camTar(1)], ...
                         'YData', [camPos(2) camTar(2)], 'ZData', [1 1])
    end
    pause(.02);
end
annotation('textbox', 'fontsize', 120, 'string', 'win!')
end
function box = boxPatch(o,l,w,h)                                % Box Maker
vertices = repmat(o,1,8) + [-l -l -l -l  l  l  l  l  
                            -w -w  w  w -w -w  w  w
                            -h  h  h -h -h  h  h -h];
faces = [1 5 1 3 1 2; 
         2 6 2 4 4 3;
         3 7 6 8 8 7; 
         4 8 5 7 5 6]';  
box = patch('Vertices', vertices',     'Faces', faces,...
           'FaceColor', rand(1,3), 'FaceAlpha', .8);
end
function controls(~,event)                                      % Buttons
 global camAng camPos viewSwitch boxNum
 switch event.Key
     case 'leftarrow'                                       
        camAng = camAng + 1;
     case 'rightarrow'                                                                 
        camAng = camAng - 1;   
     case 'downarrow'                                       
        camPos = camPos - [cosd(camAng), sind(camAng), 0];
     case 'uparrow'                                         
        camPos = camPos + [cosd(camAng), sind(camAng), 0];
     case 's'                                               
         if viewSwitch == 1
             viewSwitch = 0;
         else
             viewSwitch = 1;
         end
     case '1'
         boxNum = 20;  newGame;
     case '2'
         boxNum = 50;  newGame;
     case '3'
         boxNum = 100; newGame;
     case '0'
         close all;
 end
end