clc
close all

readOpt = 1; % 1: Compression, 2: Shear
switch readOpt
    case 1
        filePath = "C:\Users\kyoko\Desktop\LSDEM_2D\test\";
                pos = importdata([filePath + 'positions_caicos_2.dat']);
                rotation = importdata([filePath + 'rotations_caicos_2.dat']);
%         pos = importdata([filePath + 'positions_iso.dat']);
%         rotation = importdata([filePath + 'rotations_iso.dat']);
%         wallPos = importdata([filePath + 'wallPos_caicos_2.dat']);
    case 2
        filePath = "";
        pos = importdata([filePath + 'positions_caicos_3.dat']);
        rotation = importdata([filePath + 'rotations_caicos_3.dat']);
        velocity = importdata([filePath + 'velocities_caicos_3.dat']);
        
end
mkVideo = 0;

xmin = -20; xmax = 1100; ymin = -10; ymax = 500;
numAll = 150;
copyNum = numAll/parTypeNum;
radAll = repmat(parTypeRad, copyNum, 1);

relPos = parTypePos - repmat(parTypeCOM, 1, parSurPoint(1));
relPosAll = repmat(relPos, copyNum, 1);
relRotAll = repmat(parTypeRot, copyNum, 1);

temPosX = nan(numAll * (numSurPts + 1), 1);
temPosY = nan(numAll * (numSurPts + 1), 1);

VideoFrameRate = 20;
if mkVideo
    vid = VideoWriter('LSDEM(Irregular).avi', 'Uncompressed AVI');
    vid.FrameRate = VideoFrameRate;
    open(vid);
end

%%
for ii = 1 : 1 : 350
    disp(ii);
    startInd = (ii - 1) * numAll + 1;
    endInd = ii * (numAll);
    posNow = pos(startInd : endInd, :); % COM positions in this step
    rotNow = rotation(startInd : endInd) - relRotAll;
    %     veloNow = velocity(startInd : endInd, :);
    clf;
    for jj = 1 : numAll
        xx = relPosAll(jj, 1 : 2 :end) * cos(rotNow(jj))...
            - relPosAll(jj, 2 : 2 :end) * sin(rotNow(jj)) + posNow(jj, 1);
        yy = relPosAll(jj, 1 : 2 :end) * sin(rotNow(jj))...
            + relPosAll(jj, 2 : 2 :end) * cos(rotNow(jj)) + posNow(jj, 2);
        temPosX((jj - 1) * (numSurPts + 1) + 1 : jj * (numSurPts + 1) - 1) = xx';
        temPosY((jj - 1) * (numSurPts + 1) + 1 : jj * (numSurPts + 1) - 1) = yy';
    end
    
    hold on
    %     % Velocity part
    %     scatter(posNow(:, 1), posNow(:, 2), 50, veloNow(:, 1), 'filled');
    %     colormap(jet)
    %     colorbar; %caxis([-0.5 0.5]);
    
    % Particle part
    line(temPosX, temPosY, 'linewidth', 1, 'color', 'k')
    
    %     % Wall part
    %     line([0 1500], [0 0], 'linestyle', '-', 'color', 'r', 'linewidth', 1)
    %     line([0 1500], [wallPos(ii, 1) wallPos(ii, 1)], 'linestyle', '-', 'color', 'r', 'linewidth', 1)
    
    %     set(gcf, 'position', [300 150 900 650], 'color', 'w');
    set(gcf, 'color', 'w');
    axis equal;
    axis off;
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    if ii < 20
               % pause(.3);
    end
    drawnow;
    if mkVideo
        frame = getframe(gca);
        writeVideo(vid, frame);
    end
end

if mkVideo
    close(vid);
end