clc
clear
close all

numSurPts = 151;
fid = fopen('caicos001.dat', 'r');
a = fgetl(fid);
parTypeNum = str2num(a);

% Initilizations for other information
parTypeMass = zeros(parTypeNum, 1);
parTypeCOM = zeros(parTypeNum, 2);
parTypeVel = zeros(parTypeNum, 2);
parTypeIner = zeros(parTypeNum, 1);
parTypeRot = zeros(parTypeNum, 1);
parLSCOM = zeros(parTypeNum, 2);
parSurPoint = zeros(parTypeNum, 1);
parTypePos = zeros(parTypeNum, 2*numSurPts); % Number of surface points x 2
parTypeRad = zeros(parTypeNum, 1);

parGridNum = zeros(parTypeNum, 2);
parPhiVal = cell(parTypeNum, 1);

for ii = 1 : parTypeNum
    % Every particle has 15 rows for information
    % 01 -- Mass
    a = fgetl(fid);
    parTypeMass(ii) = str2num(a);
    % 02 -- Actual COM !!!
    a = fgetl(fid);
    parTypeCOM(ii, :) = str2num(a);
    % 03 -- Velocity
    a = fgetl(fid);
    parTypeVel(ii, :) = str2num(a);
    % 04 -- Moment of inertia
    a = fgetl(fid);
    parTypeIner(ii) = str2num(a);
    % 05 -- Rotation
    a = fgetl(fid);
    parTypeRot(ii) = str2num(a);
    % 06 -- Omega
    a = fgetl(fid);
    a = str2num(a);
    % 07 -- COM for level-set!!!
    a = fgetl(fid);
    parLSCOM(ii, :) = str2num(a);
    % 08 -- Number of surface points
    a = fgetl(fid);
    parSurPoint(ii) = str2num(a);
    % 09 -- Positions of surface points
    a = fgetl(fid);
    a = str2num(a);
    parTypePos(ii, :) = a;
    % 10 -- Bounding box radius
    a = fgetl(fid);
    parTypeRad(ii, 1) = str2num(a);
    % 11 -- Grid number in x and y axes
    a = fgetl(fid);
    parGridNum(ii, :) = str2num(a);
    % 12 -- \Phi values in every nodes
    a = fgetl(fid);
    parPhiVal{ii, 1} = str2num(a);
    
    for jj = 1 : 3
        a = fgetl(fid);
        a = str2num(a);
    end
end

%% View positions of all particles
viewOpt = 3; % 1 => Total with one color; 2 => Total with different colors;
% 3 => Single with phi values; 4 => Single with contour; 5 => Single with black bgc
switch viewOpt
    case 1
        aaa = parTypePos';
        aaa = aaa(:);
        scatter(aaa(1 : 2 : end), aaa(2 : 2 : end));
    case 2
        hold on
        for ii = 1 : parTypeNum
            scatter(parTypePos(ii, 1 : 2 : end), parTypePos(ii, 2 : 2 : end));
        end
    case 3
        graInd = 3;
        xxx = (0 : parGridNum(graInd, 1) - 1)'; yyy = 0 : parGridNum(graInd, 2) - 1;
        xx = repmat(xxx, 1, parGridNum(graInd, 2)); yy = repmat(yyy, parGridNum(graInd, 1), 1);
        xx1 = xx(:); yy1 = yy(:);
        rotNow = -parTypeRot(graInd);
        relPos = parTypePos(graInd, :) - repmat(parTypeCOM(graInd, :), 1, numSurPts);
        
        temX = relPos(1 : 2 : end) * cos(rotNow)...
            - relPos(2 : 2 : end) * sin(rotNow) + parLSCOM(graInd, 1);
        temY = relPos(1 : 2 : end) * sin(rotNow)...
            + relPos(2 : 2 : end) * cos(rotNow) + parLSCOM(graInd, 2);
        %         temX = parTypePos(graInd, 1 : 2 : end); temY = parTypePos(graInd, 2 : 2 : end)
        
        relLS = parLSCOM(graInd, :) - parTypeCOM(graInd, :);
        %         temLSX =  relLS(1) * cos(rotNow) - relLS(2) * sin(rotNow) + parLSCOM(graInd, 1);
        %         temLSY = relLS(1) * sin(rotNow) + relLS(2) * cos(rotNow) + parLSCOM(graInd, 2);
        temLSX = parLSCOM(graInd, 1);
        temLSY = parLSCOM(graInd, 2);
        
        zz1 = parPhiVal{graInd, 1}(:);
%                 scatter(xx1(zz1 <= 0), yy1(zz1 <= 0), 80, zz1(zz1 <= 0), 'filled')
scatter(xx1(zz1 > 0), yy1(zz1 > 0), 80, zz1(zz1 > 0), 'filled')
%         scatter(xx1, yy1, 80, zz1, 'filled')
        colormap(jet);
        colorbar;
        hold on
        scatter(temLSX, temLSY, 200, 'kx', 'linewidth', 2);
        plot(temX, temY, 'k', 'linewidth', 2)
        %         caxis([-10 0])
        %         colormap(flipud(jet))
        colormap(jet)
    case 4
        graInd = 1;
        xxx = (0 : parGridNum(graInd, 1) - 1)'; yyy = 0 : parGridNum(graInd, 2) - 1;
        xx = repmat(xxx, 1, parGridNum(graInd, 2));
        yy = repmat(yyy, parGridNum(graInd, 1), 1);
        zz = reshape(parPhiVal{graInd, 1}(:), parGridNum(graInd, 1), parGridNum(graInd, 2));
        contour(xx, yy, zz)
        colormap(jet);
        colorbar;
        hold on
        scatter(parLSCOM(graInd, 1), parLSCOM(graInd, 2), 200, 'kx', 'linewidth', 2);
    case 5
        graInd = 21;
        rotNow = -parTypeRot(graInd);
        relPos = parTypePos(graInd, :) - repmat(parTypeCOM(graInd, :), 1, 50);
        
        temX = relPos(1 : 2 : end) * cos(rotNow)...
            - relPos(2 : 2 : end) * sin(rotNow) + parLSCOM(graInd, 1);
        temY = relPos(1 : 2 : end) * sin(rotNow)...
            + relPos(2 : 2 : end) * cos(rotNow) + parLSCOM(graInd, 2);
        
        pp = polyshape(temX, temY);
        % Case 1
                f = plot(pp, 'FaceColor', 'k', 'edgecolor', 'none', 'FaceAlpha', 1);
                set(gcf, 'color', 'w')
        
        % Case 2
%         f = plot(pp, 'FaceColor', 'w', 'edgecolor', 'none', 'FaceAlpha', 1);
%         set(gcf, 'color', 'k')
        
        axis off
    otherwise
        error('No such options.');
end
axis equal
box on
% print(gcf, '-dpng', '-r600', 'irregular02');