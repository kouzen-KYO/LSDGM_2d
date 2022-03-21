clc

makevideo = 0;                                      % Make video or not
startStep = 1;                                      % Step started to watch animation
skipStep = 1;                                      % Step skipped to watch animation
endStep = 550;                                  % Step endded to watch animation
velo = .1;

colorMode = zeros(numAll, 3);
colorMode(:) = 0;

hold on

if makevideo
    v = VideoWriter(['vx(',string,').avi'], 'Uncompressed AVI');
    v.FrameRate = 30;
    open(v);
end

for ii = startStep : skipStep : endStep
    startInd = (ii - 1) * numAll + 1;
    endInd = ii * (numAll);
    posNow = pos(startInd : endInd, 2); % COM positions in this step
    veloNow = velocity(startInd : endInd, 1);
    
    if ii == startStep
        f = scatter(veloNow, posNow);
        f.Marker = '.';
        f.SizeData = 50;
        f.MarkerEdgeAlpha = 0.9;
        f.CData = colorMode(1 : numAll, :);
        
        axis([-velo velo ymin ymax])
        xlabel('$V_x$', 'interpreter', 'latex', 'Fontsize', 16);
        ylabel('$H/d_s$', 'interpreter', 'latex', 'Fontsize', 16);
        box on
        grid on
        set(gca, 'fontname', 'times new roman', 'fontsize', 16)
    else
        f.XData = veloNow;
        f.YData = posNow;
        
        drawnow;
        
        if makevideo
            frame = getframe(gca);
            writeVideo(v, frame);
        end
    end
%     pause(.2);
end

if makevideo
    close(v);
end
