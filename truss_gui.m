%% Initialize Globals    
clear; clc; close all;
global programOn;           programOn       = 1;            % Mouse/ program On
global mouseClick;          mouseClick      = 0;                      
global mouseDown;           mouseDown       = 0;
global mouseUp;             mouseUp         = 0;  
infoBoxSwitch = 0;
labelSwitch   = 0;
buildMode     = 'truss';
gridSize      = [20 20];
E             = 200;
trussCount    = 0;
forceCount    = 0;
fixedPtCount  = 0;
[fig,    ax, sideBar, sideColor,  topBar,  topColor, ...
        tLine, oLine,   uLine,     fLine,  fArrow] = Canvas(gridSize);
%% Loop
while programOn == 1
    
    [figPos, axPos, gridPos, sideKey, topKey, paperSwitch] = posFcn(fig,ax);
    
    if mouseClick == 1
        mouseDown  = 1;
        mouseClick = 0;
        
        if paperSwitch == 1
            anchorPos = gridPos;
        end
        
        if topKey > 0 && topKey < 5
            set(topBar(topKey), 'Color', [.6 .4 .4]);
        end
        
        switch topKey
            case 1; clearCanvas(tLine, uLine, fLine, oLine, fArrow); 
            case 2; [elNode, oNode, fNode, nodeCoord] = getENFO(tLine, fLine, oLine); 
                     if labelSwitch == 0
                        labelSwitch = 1;
                        label = labelMake(elNode,nodeCoord,ax);    
                     else
                        labelSwitch = 0;
                        delete(label);
                     end
            case 3;  if infoBoxSwitch == 0
                        infoBoxSwitch = 1;
                        [elNode, oNode, fNode, nodeCoord] = getENFO(tLine, fLine, oLine); 
                        [K0, F, gK0, gF, uLine] = kfMake(elNode, nodeCoord, oNode, fNode, E, uLine);
                        infoBox = infoBoxFcn(fig, ax, E, elNode, nodeCoord, gK0, K0, gF, fNode, oNode);
                     else
                        infoBoxSwitch = 0;
                        delete(infoBox);
                     end
           case 4;      [elNode, oNode, fNode, nodeCoord] = getENFO(tLine, fLine, oLine);
                        [K0, F, gK0, gF, uLine] = kfMake(elNode, nodeCoord, oNode, fNode, E, uLine);
        end
  
        if sideKey > 0 && sideKey < 10                              % Side Key:                  
            set(sideBar(sideKey), 'BackGroundColor', [.8 .8 .8]);
        end
        
        switch sideKey         
            case 9; buildMode =   'truss';                      % Side Key: Truss Mode
            case 8; buildMode =   'force';                      % Side Key: Force Mode;
            case 7; buildMode = 'fixedPt';                      % Side Key:  
            case 6 
                xNew = get(gca, 'XLim') + [1 -1];
                yNew = get(gca, 'YLim') + [1 -1];
                set(gca,'XTick', min(xNew):max(xNew), 'XLim', xNew, ...
                        'YTick', min(yNew):max(yNew), 'YLim', yNew);
            case 5
                xNew = get(gca, 'XLim') + [-1 1];
                yNew = get(gca, 'YLim') + [-1 1];
                set(gca,'XTick', min(xNew):max(xNew), 'XLim', xNew, ...
                        'YTick', min(yNew):max(yNew), 'YLim', yNew);
            case 4
                yNew = get(gca, 'YLim') + [1 1];
                set(gca, 'YTick', min(yNew):max(yNew), 'YLim', yNew);
            case 3
                yNew = get(gca, 'YLim') + [-1 -1];
                set(gca, 'YTick', min(yNew):max(yNew), 'YLim', yNew);
            case 2 
                xNew = get(gca, 'XLim') + [-1 -1];
                set(gca, 'XTick', min(xNew):max(xNew), 'XLim', xNew);    
            case 1
                xNew = get(gca, 'XLim') + [1 1];
                set(gca, 'XTick', min(xNew):max(xNew), 'XLim', xNew);
        end
    end
    
    % Mouse Down
    if mouseDown == 1 
        if paperSwitch == 1
            switch buildMode
                case 'truss'
                    set(tLine(trussCount + 1), ...
                        'XData', [anchorPos(1), gridPos(1)], ...
                        'YData', [anchorPos(2), gridPos(2)]);
                case 'force'  
                    set(fLine(forceCount + 1), ...
                        'XData', [anchorPos(1), gridPos(1)], ...
                        'YData', [anchorPos(2), gridPos(2)]);
                    fArrow = fArrowSet(fArrow, gridPos, anchorPos,  forceCount);
                case 'fixedPt' 
                    set(oLine(fixedPtCount + 1), ...
                        'XData', [anchorPos(1), gridPos(1)], ...
                        'YData', [anchorPos(2), gridPos(2)]);
            end
        end
    end
    
    % Mouse Up              
    if mouseUp == 1
        mouseUp   = 0;
        mouseDown = 0;
        barReset(topBar, topColor, sideBar, sideColor);
        if paperSwitch == 1
            anchorDist = norm(gridPos - anchorPos);
        	switch buildMode
            case 'truss'
                if anchorDist > .8
                    trussCount = trussCount + 1;
                else
                    set(tLine(trussCount + 1),   'XData', [], 'YData', []);
                end        
            case 'force'
                if anchorDist > .8
                    forceCount = forceCount + 1;
                else
                    set(fLine(forceCount),   'Xdata', [], 'YData', []);
                end
            case 'fixedPt'
                if anchorDist < .8
                    fixedPtCount = fixedPtCount + 1;
                else
                    set(oLine(fixedPtCount + 1), 'Xdata', [], 'YData', []);
                end
            end
        end
    end
    pause(.01);
end
%% Mouse Pos Fcn
function [figPos, axPos, gridPos, sideKey, topKey, paperSwitch] = posFcn(fig,ax)
cPos     =  get(  0, 'PointerLocation');
figSize  =  get(fig,        'Position');
axSize   =  get(ax,         'Position');
xGrid    = get(ax, 'XLim');
yGrid    = get(ax, 'YLim');
gridSize = [xGrid(2) - xGrid(1), yGrid(2) - yGrid(1)];
figPos   = (cPos   - figSize(1:2))./figSize(3:4);
axPos    = (figPos -  axSize(1:2))./axSize(3:4).*gridSize + [xGrid(1), yGrid(1)];
gridPos  = floor([.5 .5] + axPos);
if figPos(1) > .9 && figPos(2) < .9
    sideKey = floor(10*figPos(2) + 1);
else
    sideKey = 0;
end
if figPos(2) > .9 && figPos(2) < 1 
    topKey = floor(10*figPos(1)) - 5;
else
    topKey = 0;
end
if figPos(1) > 0  && figPos(1) < .9 && figPos(2) > 0 && figPos(2) < .9
    paperSwitch = 1;
else
    paperSwitch = 0;
end
end
%% Arrow Set Function
function fArrow = fArrowSet(fArrow, gridPos, anchorPos, forceCount)
xStar = (anchorPos(1) - gridPos(1));
yStar = (anchorPos(2) - gridPos(2));
L     = sqrt(xStar^2 + yStar^2);
theta = atan2(yStar,xStar);
arrowL = gridPos + L/8*[cos(theta + .4), sin(theta + .4)];            
arrowR = gridPos + L/8*[cos(theta - .4), sin(theta - .4)];
set(fArrow(forceCount + 1), 'XData',  [arrowL(1), gridPos(1), arrowR(1)], ...
                            'YData',  [arrowL(2), gridPos(2), arrowR(2)]); 
end
%% Make [nodeCoord, elNode, fNodes, oNodes] from Plot Data (getENFO) 
function [elNode, oNode, fNode, nodeCoord] = getENFO(tLine, fLine, oLine)    
elCoord = [];
elCount = 0;
for i = 1:20                                            % Get Element Coordinates
    x = get(tLine(i), 'XData');
    y = get(tLine(i), 'YData');
    if numel(x) == 2
        elCount = elCount + 1;
        elCoord(:,elCount) = [x(1); y(1); x(2); y(2)];
    end
end
nodeCoord   = [];
uniqueCount = 0;
allNode = [elCoord(1:2,:), elCoord(3:4,:)];              % Find Unique Nodes
for n1 = 1:numel(allNode)/2                                
    unique = 1; 
    n1Pos  = allNode(:, n1);
    for n0 = 1:numel(nodeCoord)/2                              
        n0Pos = nodeCoord(:,n0);
        if norm(nodeCoord(:,n0) - allNode(:,n1)) < .1 
            unique = 0;
        end
    end
    if unique == 1                                           
        uniqueCount = uniqueCount + 1;
        nodeCoord(:, uniqueCount) = n1Pos; 
    end
end
elNode = zeros(2,elCount);                                      % generate elNodes
for n = 1:numel(nodeCoord)/2
    for e = 1:elCount
        if norm(nodeCoord(:,n) - elCoord(1:2,e)) < .1
            elNode(1,e) = n;
        end
        if norm(nodeCoord(:,n) - elCoord(3:4, e)) < .1
            elNode(2,e) = n;
        end
    end
end
% Part 2: get oNodes
oNode = [];                                % start oNodes as empty
oCount = 0;
for o = 1:numel(oLine)                             % Loop through oLine
    ox = get(oLine(o), 'XData');                   % Get oLine x,y
    oy = get(oLine(o), 'YData');
    if numel(ox) > 0
    oPos = [ox(1); oy(1)];
        for n = 1:numel(nodeCoord)/2                % Loop through node coordinates
            if norm(oPos - nodeCoord(:,n)) < .1     % if oPos = nodePos, add fixed pt                    
                oCount = oCount + 1;
                oNode(:,oCount) = n;
            end
        end
    end
end
% Part 3: get fNodes [3xN:] row 1: [fNodeNum; fx; fy]
fCount = 0;
fNode = [];
for f = 1:numel(fLine)                             
    px = get(fLine(f), 'XData');                 
    py = get(fLine(f), 'YData');
    if numel(px) > 0
        pPos = [px(1);py(1)];
        fMag = [px(2) - px(1); py(2) - py(1)];
        for n = 1:numel(nodeCoord)/2                % Loop through node coordinates
            if norm(pPos - nodeCoord(:,n)) < .1     % if oPos = nodePos, add fixed pt                    
                fCount = fCount + 1;
                fNode(:,fCount) = [n;fMag];
            end
        end
    end
end
end
%% Assembly Function: get K0 and F from Node Data (KOF)  
function [K0, F, gK0, gF, uLine] = kfMake(elNode, nodeCoord, oNode, fNode, E, uLine)
gDof = numel(nodeCoord);
gK0  =  zeros(gDof,gDof);
for el = 1:numel(elNode)/2
    nodeANum = elNode(1,el);
    nodeBNum = elNode(2,el);
        
    nodeAPos = nodeCoord(:,nodeANum);
    nodeBPos = nodeCoord(:,nodeBNum);
    
    dx = nodeBPos(1) - nodeAPos(1);
    dy = nodeBPos(2) - nodeAPos(2);
    L    = sqrt(dx^2 + dy^2);
    ccEL = E*(dx^2)/L^3;
    csEL = E*(dx*dy)/L^3;
    ssEL = E*(dy^2)/L^3;
    k0 = [ccEL,  csEL, -ccEL, -csEL;
          csEL,  ssEL, -csEL, -ssEL;
         -ccEL, -csEL,  ccEL,  csEL;
         -csEL, -ssEL,  csEL,  ssEL];
    
    dofA = [2*nodeANum-1, 2*nodeANum];
    dofB = [2*nodeBNum-1, 2*nodeBNum];
    
    gK0(dofA, dofA) = k0(1:2,1:2) + gK0(dofA, dofA);
    gK0(dofA, dofB) = k0(1:2,3:4) + gK0(dofA, dofB);
    gK0(dofB, dofA) = k0(3:4,1:2) + gK0(dofB, dofA);
    gK0(dofB, dofB) = k0(3:4,3:4) + gK0(dofB, dofB);
end
 
oDof = [];                                      % get fixed dofs
for i = 1:numel(oNode)
    oDof = [oDof, 2*oNode(i) - 1];
    oDof = [oDof, 2*oNode(i)];
end
gF = zeros(gDof, 1);                            % get forces
for i = 1:numel(fNode)/3
     gF(2*fNode(1,i) - 1) = fNode(2,i);
     gF(2*fNode(1,i))     = fNode(3,i);
end
activeDof       = 1:gDof;
activeDof(oDof) = [];
F               = gF(activeDof);
K0              = gK0(activeDof, activeDof);
U = inv(K0)*(F);
gU = zeros(numel(nodeCoord), 1);
for i = 1:numel(activeDof)
    gU(activeDof(i)) = U(i);
end
uCoord = zeros(2,numel(nodeCoord)/2);
for i = 1:numel(gU)/2
    uCoord(1,i) = gU(2*i - 1);
    uCoord(2,i) = gU(2*i);
end
dispCoord = nodeCoord + uCoord;
for i = 1:numel(elNode)/2
    n1 = dispCoord(:, elNode(1,i));
    n2 = dispCoord(:, elNode(2,i));
    set(uLine(i), 'XData', [n1(1), n2(1)], 'YData', [n1(2), n2(2)]);
end
end
%% Info Box    
function infoBox = infoBoxFcn(fig, ax, E, elNode, nodeCoord, ...
                                        gK0, K0, gF, fNode, oNode)
dof       = numel(nodeCoord);
gridLimX  = get(ax, 'XLim');
gridLimY  = get(ax, 'YLim');
w = gridLimX(2) - gridLimX(1);
h = gridLimY(2) - gridLimY(1);
% Background
infoBox(1) = patch('Faces', [1 2 3 4], 'Vertices', ...
      [gridLimX(1),   gridLimX(2), gridLimX(2), gridLimX(1);
       gridLimY(1),   gridLimY(1), gridLimY(2), gridLimY(2)]', ...
       'FaceColor',    [.1 .1 .1], ...
       'EdgeColor',    [.4 .5 .6]);
% Data
for i = 2:10
    infoBox(i) = text(  .25, .95, ' ', ...
            'FontName', 'FixedWidth',    'Color',   [.8 .8 .8], ... 
           'FontUnits', 'Normalized', 'FontSize',          .05, ...
 'HorizontalAlignment',       'left',    'Units', 'normalized', ...
   'VerticalAlignment',        'top');
end
infoString = {sprintf('%12s %3.0f',         'E:',   string(E)), ...
              sprintf('%12s %3.0f',       'els:',   numel(elNode)/2), ...
              sprintf('%12s %3.0f',     'nodes:', numel(nodeCoord)/2), ...
              sprintf('%12s %3.0f',    'forces:', numel(fNode)/3)};       
set(infoBox(2), 'String', infoString, 'HorizontalAlignment', 'right');
 
colormap(gray)
set(infoBox(3), 'String',  'K0 (non-reduced)', ...
              'Position',  [.6 .95]);
infoBox(4) = imagesc('XData', [gridLimX(1) +  .6*w, gridLimX(1) + .9*w], ...
                     'YData', [gridLimY(1) + .55*h, gridLimY(1) + .85*h], ...
                     'CData', -abs(flipud(gK0)));  
                  
fString = sprintf(repmat('%3.0f', 1, numel(gF)), gF);
set(infoBox(5), 'String', strcat('F:',fString), 'Position', [.1 .7]);
elString = {                                       'elNode:'; 
     sprintf(repmat('%3.0f', 1, numel(elNode)/2), elNode(1,:));
     sprintf(repmat('%3.0f', 1, numel(elNode)/2), elNode(2,:))};
set(infoBox(6),    'String', elString, 'Position', [.1 .5]); 
nodeString = {                                       'nodeCoord:'; 
     sprintf(repmat('%4.0d', 1, numel(nodeCoord)/2), nodeCoord(1,:));
     sprintf(repmat('%4.0d', 1, numel(nodeCoord)/2), nodeCoord(2,:))};
set(infoBox(7),    'String', nodeString, 'Position', [.1 .3]); 
if numel(oNode) > 0
    oString = strcat('oNode:', sprintf(repmat('%4.0f', 1, numel(oNode)), oNode));
    set(infoBox(8),  'String', oString, 'Position', [.5 .4])
end
end
%% Node labels
function label = labelMake(elNode,nodeCoord,ax)
labelCount = 0;
for i = 1:numel(nodeCoord)/2
    labelCount = labelCount + 1;
    x = nodeCoord(1,i);
    y = nodeCoord(2,i);
    label(labelCount) = text('Position',       [x, y],     'String',  string(i), ...
                                'color',  [ .8 .8 .8],  'fontunits',   'normal', ...   
                      'backgroundcolor',   [.1 .1 .1], 'fontweight',     'bold', ...         
                            'edgeColor',   [ 0  0  0],   'fontsize',        .04, ...
                               'margin',           2);
end
for i = 1:numel(elNode)/2
    labelCount = labelCount + 1;
    n1    = nodeCoord(:, elNode(1,i));
    n2    = nodeCoord(:, elNode(2,i));
    midPt = (n1 + n2)/2;
    
    label(labelCount) = text(midPt(1), midPt(2), string(i), ...
               'color',  [.1 .1 .1],    'fontUnits', 'normal',  ...
     'backgroundcolor',  [.8 .8 .8],   'fontweight',   'bold',  ...
           'edgeColor',  [ 0  0  0],     'fontsize',     .04, ...
       'Margin', 2);    
end
x = get(ax, 'XLim');
y = get(ax, 'YLim');
label(labelCount + 1) = line( x,  [y(1) + 1, y(1) + 1], ...
                             'linewidth', 2, 'color', [.2 .3 .4]);
label(labelCount + 2) = line([x(1) + 1, x(1) + 1],  y, ...
                             'linewidth', 2, 'color', [.2 .3 .4]);
labelCount = labelCount + 2;
numx = 6;
intx = floor((max(x) - min(x))/numx);
for xText = min(x) + intx:intx:max(x)
    labelCount = labelCount + 1;
    label(labelCount) = text(xText, y(1) + 1, string(xText), ...
                                 'color', [.3 .1 .1], ...
                   'horizontalalignment', 'center', 'verticalalignment', 'top', ...
                   'fontunits', 'normalized', 'fontsize', .03, 'clipping', 'off');
end
numy = 6;
inty = floor((max(x) - min(x))/numx);
for yText = min(y) + inty:inty: max(y)
    labelCount = labelCount + 1;
    label(labelCount) = text(x(1) + .9, yText, string(yText), ...
                               'color', [.3 .1 .1], ...
                 'horizontalAlignment', 'Right', 'VerticalAlignment', 'Middle', ...
                 'fontunits', 'normalized', 'fontsize', .03, 'Clipping', 'on', ...
                 'parent', ax);
end
end
%% Mouse Functions
function mouseClickDownCall(~,~)                                    % Click Call
global mouseClick mouseDown
mouseClick = 1;
mouseDown  = 1;
end
function mouseClickUpCall(~,~)                                      % Up Call
global mouseUp mouseDown
mouseUp   = 1;
mouseDown = 0;
end
function closeReqCall(~,~)                                          % Close Call
global programOn
programOn = 0;
delete(gcf);
end
%% Top Bar/ Side Bar reset
function barReset(topBar, topColor, sideBar, sideColor)
for i = 1:7
  set(topBar(i),            'Color', topColor(i,:));
end
for i = 1:9
    set(sideBar(i), 'BackGroundColor', sideColor(i,:));
end
end
%% Clear Canvas
function clearCanvas(tLine, uLine, fLine, oLine, fArrow)
fExist = 0;
for i = 1:numel(fLine)
    if numel(get(fLine(i), 'Xdata')) > 1
        fExist = 1;
    end
end
if fExist == 1
   for i = 1:numel(fLine)
        set(uLine(i),    'XData', [], 'YData', []);
        set(fLine(i),    'XData', [], 'YData', []);
        set(fArrow(i),   'XData', [], 'YData', []);
   end
else
    for i = 1:numel(fLine)
        set(oLine(i), 'XData', [], 'YData', []);
        set(tLine(i), 'XData', [], 'YData', []);
    end
end
end
%% CanvasMake
function [fig,    ax, sideBar, sideColor,  topBar,  topColor, ...
        tLine, oLine,   uLine,     fLine,  fArrow] = Canvas(gridSize)
axisPosition = [.02 .02 .86  .89];                                  % -- ax pos
set(0,   'Units', 'Normalized');                                    % -- root
fig = figure('WindowButtonDownFcn', @mouseClickDownCall, ...            % -- fig
               'WindowButtonUpFcn',   @mouseClickUpCall, ...   
                 'CloseRequestFcn',       @closeReqCall, ...
                           'units',        'normalized', ...
                           'color',       [.14 .14 .14], ...
                         'menubar',              'none', ...
                        'position',      [.25 .15 .6 .7]); 
                
ax = axes('pos',    axisPosition,     'color',      [.8 .8 .8], ...   % -- ax
        'XGrid',            'on',     'YGrid',            'on', ...
         'XLim', [0 gridSize(1)],      'YLim', [0 gridSize(2)], ...
        'XTick',   0:gridSize(1),     'YTick',   0:gridSize(2), ...
    'gridcolor', [.45  .48  .45], 'gridalpha',              .6, ...
      'TickLen',           [0 0],    'parent',             fig);
startPts = [.6 .7 .8 .9  0  .2 .4 ];
lens     = [.1 .1 .1 .1 .2  .2 .2 ];
topStrings = {  'new',  'grid',  'data', 'soln', ...
      '  truss\Thetapt', '\delta',  '\delta '};
topColor = [.7 .6 .7 .6 .4 .7 .7;
            .7 .5 .7 .5 .5 .7 .7;
            .7 .4 .7 .4 .6 .7 .7]';  
  
for i = 1:numel(startPts)                                           % Top Bar
    topBar(i) = annotation( 'textbox', [startPts(i), .94, lens(i), .06], ... 
   'verticalAlignment',      'middle',  'fontUnits',           'normal', ...
     'BackGroundColor', [.01 .01 .01],   'fontsize',                .05, ...    
           'edgecolor', [.01 .01 .01],      'Color',     topColor(i,:), ...   
               'String', topStrings(i));                                             
end
annotation('Rectangle',  [0, .93, 1, .005], ...                % Divider
           'Facecolor',         [.4 .5 .6], ...
           'EdgeColor',        [.4 .5 .6]);
                                            
sideStrings ={'\rightarrow',  '\leftarrow', '\downarrow', '\uparrow', ...
                   '-', '+', '\partial',     '\it{f}',    '//'};
  
sideColor = [.6  .7   .6 ;  % green
              .6  .7   .6 ;  % green
              .64 .64  .64;  % gray
              .64 .64  .64;  % gray
              .68 .68  .6 ;  % tan 
              .68 .68  .6 ;  % tan 
              .6  .68  .68;  % blue
              .6  .68  .68;  % blue
              .6  .68  .68]; % blue
for e = 1:9                                                         % SideBar
sideBar(e) = annotation('TextBox', [.9, (.1*e -.082), .0855, .095], ...                         
                'BackGroundColor',   sideColor(e, :),      'Color',   [0 0 0], ...        
            'HorizontalAlignment',          'center',  'FontUnits',  'normal', ... 
              'VerticalAlignment',          'middle',   'FontSize',       .08, ...      
                         'string',    sideStrings(e), 'FontWeight',     'bold');
                         
end
set(sideBar(9), 'Position', [.9, .82, .085, .09], 'fontsize', .08);
for i = 1:20                                                            
    tLine(i) = line('color', [.1 .2 .1], 'linewidth',   6, ...         % Truss Lines
                    'XData',          0,  'userdata',   i, ...
                    'YData',          0);         
    fLine(i) = line('color', [.6 .2 .2], 'linewidth',   4,  ...        % Force Lines
                    'XData',         [], 'userdata',    i,  ...
                    'YData',         []);
    oLine(i) = line('color',    [0 0 0],   'linestyle',    'none', ... % FixedPt Markers
                    'XData',         [],    'userdata',         i, ...
                    'YData',         [],      'marker', 'diamond', ...
          'markerFaceColor', [.9 .9 .9],  'markerSize',        16, ...
                'linewidth',         2);   
    fArrow(i)= line('color', [.6 .2 .2],  'linewidth',    4, ...            % Force Arrows
                    'XData',         [],   'userdata',    i, ...
                    'YData',         []);   
    uLine(i) = line('color',  [.2 .2 1],  'lineWidth',   2, ...            % u Lines
                    'XData',         [], ...
                    'YData',         []);
end
end