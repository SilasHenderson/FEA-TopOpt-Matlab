% Asteroids(ish) (Silas Henderson)

clear; clc; close all;
fig = figure(  'WindowKeyPressFcn', @keyDownFcn,   'color', 'black', ...
             'WindowKeyReleaseFcn',   @keyUpFcn, 'menubar',  'none');
        
ax  = axes('units',  'normal', 'position', [.01 .01 .98 .88], ...
           'color',   'white',  'TickLen', [0 0], ...
           'XGrid',      'on',    'YGrid', 'on', ...
            'XLim', [-10, 10],     'YLim', [-10, 10]);

info = annotation('textbox', 'units', 'normal', 'position', [0 .9 1 .1], ...
         'string', 'Points:   0', 'color',  'white', 'fontsize', 18);       
        
ship  = line('linewidth', 5, 'color', [.3 .3 .3]); 

laser = line(nan,nan,'linestyle',  'none',      'color', [.4 .4 .8], ...
                        'marker',     '.', 'markersize',        20);          
rocks = line(nan,nan,'linestyle',  'none',      'color', [.3 .5 .2],...
                        'marker',     '.', 'markersize',        50);           
           
global laserOn shipVel shipPhiVel;     
           
shipPos  = [0, 0];  shipVel    = 0;    
shipPhi  = 0;       shipPhiVel = 0;
laserPos = [];      laserVel   = [];
rocksPos = [];      rocksVel   = [];

t50 = 0;   t200 = 0; 
t   = 0;   t0   = 0;   
tic;       t20  = 0;    dt = 0;

points = 0;

while toc < 100
    t0 = t;
    t  = toc;
    dt = t - t0;
    
    if abs(t50 - t) > .05
        t50 = t;
        if laserOn == 1
            laserPos = [laserPos; shipPos];
            laserVel = [laserVel; [cos(shipPhi)/10, sin(shipPhi)/10]];
        end
    end
    if abs(t200 - t) > .2
        t200 = t;
        rocksPos = [rocksPos; [ - 9*rand,    -9.9]];
        rocksVel = [rocksVel; [  rand/60, rand/30]];
    end
     
    shipPhi = shipPhi + shipPhiVel;
    shipPos = shipPos + shipVel*[cos(shipPhi), sin(shipPhi)];
       
    x = [ shipPos(1) + 1.5*cos(shipPhi), ...
          shipPos(1) + cos(shipPhi + 2*pi/3), ...
          shipPos(1) + cos(shipPhi + 4*pi/3), ...
          shipPos(1) + 1.5*cos(shipPhi) ];
        
   y = [  shipPos(2) + 1.5*sin(shipPhi), ...
          shipPos(2) + sin(shipPhi + 2*pi/3), ...
          shipPos(2) + sin(shipPhi + 4*pi/3), ...
          shipPos(2) + 1.5*sin(shipPhi) ];
    
    laserPos = laserPos + dt*100*laserVel; 
    rocksPos = rocksPos + dt*100*rocksVel;
    
    if numel(laserPos) > 0
        i = 1;
        while i < numel(laserPos)/2
            if abs(laserPos(i, 1)) > 10 || abs(laserPos(i, 2)) > 10
                laserPos(i,:) = [];
                laserVel(i,:) = [];
                points = points + 500;
                set(info, 'string', sprintf('Points %10d', points));
            else
                i = i + 1;
            end
        end
    end
    
    if numel(rocksPos) > 0
        j = 1;
        while j < numel(rocksPos)/2
            if abs(rocksPos(j, 1)) > 10 || abs(rocksPos(j, 2)) > 10
                   rocksPos(j,:) = [];
                   rocksVel(j,:) = [];
            else
                j = j + 1;
            end
        end
        
        j = 1;
        while j < numel(rocksPos)/2
            jj = 1;
            while jj < numel(laserPos)/2
                if norm(laserPos(jj,:) - rocksPos(j,:)) < 1
                    laserPos(jj,:) = [];
                    laserVel(jj,:) = [];
                    rocksPos(j,:)  = [];
                    rocksVel(j,:)  = [];
                    break;
                else
                    jj = jj + 1;
                end
            end
            j = j + 1;
        end
        
        for i = 1:numel(rocksPos)/2
            if norm(rocksPos(i,:) - shipPos) < 1
                color = .7*rand(1, 3);
                set(gcf, 'color', color);
                set(info, 'edgecolor', color);
                
                points = points - 5000;
                set(info, 'string', sprintf('Points %10d', points));
            end
        end
    end
    
    if t - t20 > 1/60
        set(ship, 'XData', x, 'YData', y); 
        if numel(rocksPos > 0)
            set(rocks, 'XData', rocksPos(:, 1)', ...
                       'YData', rocksPos(:, 2)');
        end
        if numel(laserPos > 0)
            set(laser, 'XData', laserPos(:, 1)', ...
                       'YData', laserPos(:, 2)');
        end
        drawnow;
    end
end

function keyDownFcn(~, event)
    global shipVel shipPhiVel laserOn;
    switch event.Key
        case 'uparrow',     shipVel    =  .05;
        case 'downarrow',   shipVel    = -.05; 
        case 'rightarrow',  shipPhiVel = -.05; 
        case 'leftarrow',   shipPhiVel =  .05; 
        case 'space',       laserOn    =    1;
    end
end

function keyUpFcn(~, event)
    global shipVel shipPhiVel laserOn;
    switch event.Key
        case 'uparrow',     shipVel    = 0; shipPhiVel = 0;
        case 'downarrow',   shipVel    = 0; shipPhiVel = 0;
        case 'rightarrow',  shipPhiVel = 0; shipVel = 0;
        case 'leftarrow',   shipPhiVel = 0; shipVel = 0;
        case 'space',       laserOn    = 0; shipVel = 0; shipPhiVel = 0;
    end
end