% Must complied after script(LSDEM2D_01_ReadInputFile)
clc
% close all

pos = importdata('positions_caicos001.dat');
rotation = importdata('rotations_caicos001.dat');
xmin = -10; xmax = 750; ymin = -10; ymax = 1400;
numAll = 750;
copyNum = numAll/parTypeNum;
radAll = repmat(parTypeRad, copyNum, 1);

relPos = parTypePos - repmat(parTypeCOM, 1, parSurPoint(1));
relPosAll = repmat(relPos, copyNum, 1);
relRotAll = repmat(parTypeRot, copyNum, 1);

temPosX = nan(numAll * (numSurPts + 1), 1);
temPosY = nan(numAll * (numSurPts + 1), 1);

% transX = -(545.07320305-634.37729014 - 5.5);
% transY = -(784.81266662-656.22997566 + 5.5) + 645.3;
transX = 0; transY = 0;

for ii = 1
    startInd = (ii - 1) * numAll + 1;
    endInd = ii * (numAll);
    posNow = pos(startInd : endInd, :) + repmat([transX transY], numAll, 1); % COM positions in this step
    rotNow = rotation(startInd : endInd) - relRotAll;

    for jj = 1 : numAll
        xx = relPosAll(jj, 1 : 2 :end) * cos(rotNow(jj))...
            - relPosAll(jj, 2 : 2 :end) * sin(rotNow(jj)) + posNow(jj, 1);
        yy = relPosAll(jj, 1 : 2 :end) * sin(rotNow(jj))...
            + relPosAll(jj, 2 : 2 :end) * cos(rotNow(jj)) + posNow(jj, 2);
        temPosX((jj - 1) * (numSurPts + 1) + 1 : jj * (numSurPts + 1) - 1) = xx';
        temPosY((jj - 1) * (numSurPts + 1) + 1 : jj * (numSurPts + 1) - 1) = yy';
    end
    
    line(temPosX, temPosY, 'linewidth', 1, 'color', 'r')
    set(gcf, 'position', [300 150 800 600], 'color', 'w');
    axis equal;
%     axis off;
    xlim([xmin xmax]);
    ylim([ymin ymax]);
end