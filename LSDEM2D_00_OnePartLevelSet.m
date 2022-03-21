clear
clc
close all

% Some initial settings
realArea = 400; % The particle area we roughly want (not actual area but close to)
parType = 5; % 1: disk, 2: ellipse, 3: triangle, 4: rectangle, 5: pentagon, 6: hexagon, 7: Octogon, 99: other irregular shapes (no analytical LS-functions)
numSurPoint = 150; % Number of surface points to describe the particles

%% Read the particle image for making input file
img = imread('pentagon.png');
level = graythresh(img);
im = im2bw(img, level);     % Binary values
BW = ~im;                   % Reverse 0 and 1 to make particle region denotes as 1
lbl = bwlabel(BW, 4);       % Return a matrix recorded every connected region (equal to particle number)
dist_map = bwdist(~BW);     % Get the Euclidean distance map of this particle

% Locate the distance map for this particle
[xIndex, yIndex] = find(lbl == 1);
for ii = 1 : length(xIndex)
    cData(ii) = dist_map(xIndex(ii), yIndex(ii));
end

% Find the pixel points of this particle and scale the particle size to what we want
[a, b] = find(BW == 1);
virtualArea = length(a);
incFac = sqrt(virtualArea/realArea);    % Scale factor
xIndex = xIndex/incFac; yIndex = yIndex/incFac; cData = cData/incFac;
% Shift the particle to make sure it is enclosed in a grid from (0, 0) to (xmax, ymax)
shiftX = min(xIndex) - 1; shiftY = min(yIndex) - 1;
xIndex = xIndex - shiftX; yIndex = yIndex - shiftY;
xmax = ceil(max(xIndex)) + 1; ymax = ceil(max(yIndex)) + 1;

%% Detect surface points
I_grey = rgb2gray(img);     % Transfer raw image to grey image
BW1 = edge(I_grey, 'Canny'); % Find the edge points
[B1, L] = bwboundaries(BW1, 'noholes');
B1{1, 1} = B1{1, 1}/incFac;
B1{1, 1}(:, 1) = B1{1, 1}(:, 1) - shiftX;
B1{1, 1}(:, 2) = B1{1, 1}(:, 2) - shiftY;

totalNum = size(B1{1, 1}, 1);
skipStep = round(totalNum/numSurPoint);
surPtsX = B1{1, 1}(1 : skipStep : end, 1);
surPtsY = B1{1, 1}(1 : skipStep : end, 2);
if length(surPtsX) > numSurPoint
%     surPtsX(end) = []; surPtsY(end) = [];
numSurPoint = length(surPtsX);
end
surPts = zeros(1, numSurPoint * 2);
surPts(1 : 2 : end) = surPtsX;
surPts(2 : 2 : end) = surPtsY;

%% Geometric properties
density = 1;
pp = polyshape(B1{1, 1}(:, 1), B1{1, 1}(:, 2));
[a, b] = boundingbox(pp);
[bdxmin, bdxmax, bdymin, bdymax] = deal(a(1), a(2), b(1), b(2));
totalArea = area(pp);

%  01 -- Mass
mass = totalArea * density;

% 02 -- COM (REAL) and Moment of Inertia
cx = 0; cy = 0;

for ii = 1 : length(B1{1, 1}(:, 1))
    if ii == 1
        x0 = B1{1, 1}(end, 1); y0 = B1{1, 1}(end, 2);
        x1 = B1{1, 1}(1, 1); y1 = B1{1, 1}(1, 2);
    else
        x0 = B1{1, 1}(ii - 1, 1); y0 = B1{1, 1}(ii - 1, 2);
        x1 = B1{1, 1}(ii, 1); y1 = B1{1, 1}(ii, 2);
    end
    
    cx = cx + (x0 + x1) * (x0 * y1 - x1 * y0);
    cy = cy + (y0 + y1) * (x0 * y1 - x1 * y0);
end
cx = abs(cx/6/totalArea); cy = abs(cy/6/totalArea);

% 04 -- Moment of inertia
iX = 0; iY = 0;
for ii = 1 : length(B1{1, 1}(:, 1))
    if ii == length(B1{1, 1}(:, 1))
        x0 = B1{1, 1}(end, 1); y0 = B1{1, 1}(end, 2);
        x1 = B1{1, 1}(1, 1); y1 = B1{1, 1}(1, 2);
    else
        x0 = B1{1, 1}(ii, 1); y0 = B1{1, 1}(ii, 2);
        x1 = B1{1, 1}(ii + 1, 1); y1 = B1{1, 1}(ii + 1, 2);
    end
    iX = iX + (y0^2 + y0*y1 + y1^2) * (x0 * y1 - x1 * y0);
    iY = iY + (x0^2 + x0*x1 + x1^2) * (x0 * y1 - x1 * y0);
end

iX = abs(iX)/12; iY = abs(iY)/12;
inertia = (iX + iY - totalArea * (cx^2 + cy^2)) * density;

%% Level set for some special shapes
xVal = zeros((xmax + 1) * (ymax + 1), 1);
yVal = zeros((xmax + 1) * (ymax + 1), 1);
index = 0;
for ii = 1 : ymax + 1
    yPos = ii - 1;
    for jj = 1 : xmax + 1
        index = index + 1;
        xPos = jj - 1;
        xVal(index) = xPos; yVal(index) = yPos;
    end
end
width = max(B1{1, 1}(:, 1)) - min(B1{1, 1}(:, 1));
height = max(B1{1, 1}(:, 2)) - min(B1{1, 1}(:, 2));
switch parType
    case 1 % 01 -- Disk
        lsVal = LSFunc_01_Circle(width, height, cx, cy, xmax, ymax);
    case 2 % 02 -- Ellipse
        lsVal = LSFunc_02_Ellipse(width, height, cx, cy, xmax, ymax);
    case 3 % 03 -- Triangle
        lsVal = LSFunc_03_Triangle(width, height, cx, cy, xmax, ymax);
    case 4 % 04 -- Rectangle
        lsVal = LSFunc_04_Box(width, height, cx, cy, xmax, ymax);
    case 5 % 05 -- Regular Pentagon
        lsVal = LSFunc_05_Pentagon(width, height, cx, cy, xmax, ymax);
    case 6 % 06 -- Regular Hexagon
        lsVal = LSFunc_06_Hexagon(width, height, cx, cy, xmax, ymax);
    case 7 % 07 -- Regular Octogon
        lsVal = LSFunc_07_Octogon(width, height, cx, cy, xmax, ymax);
    case 99 % 99 -- Arbitrary shapes
        lsVal = nan(1, (xmax + 1) * (ymax + 1));
        
        % LS Values for points inside the particle
        index = 0;
        for ii = 1 : ymax + 1
            yPos = ii - 1;
            for jj = 1 : xmax + 1
                index = index + 1;
                xPos = jj - 1;
                
                lsVal(index) = mean(cData(xIndex > xPos*0.99 & xIndex < xPos*1.01 &...
                    yIndex > yPos*0.99 & yIndex < yPos*1.01 ));
            end
        end
        
        % LS Values for points outside the particle 
        index = 0;
        for ii = 1 : ymax + 1
            yPos = ii - 1;
            for jj = 1 : xmax + 1
                index = index + 1;
                if ~isnan(lsVal(index)); continue; end
                xPos = jj - 1;
                temMat = B1{1, 1} - [xPos yPos];
                dis = vecnorm(temMat, 2, 2);
                lsVal(index) = -(min(dis));
            end
        end
        
        lsVal = -lsVal;
    otherwise
        error('No such particle type');
end

%% Finally view the level setting values (by color) and surface (by black connected points)
colormap(jet);
hold on
scatter(xVal, yVal, 50, lsVal, 'filled');
plot(surPtsX, surPtsY, 'o-', 'linewidth', 1, 'color', 'k');
axis equal
colorbar
axis off
set(gca, 'color', 'w')

figure
colormap(jet)
plot(surPtsX, surPtsY, 'o-', 'linewidth', 1, 'color', 'k');
axis equal;hold on
temX = reshape(xVal, xmax + 1, ymax + 1);
temY = reshape(yVal, xmax + 1, ymax + 1);
temZ = reshape(lsVal, xmax + 1, ymax + 1);
% scatter(xVal(lsVal <= 0), yVal(lsVal <= 0), 50, lsVal(lsVal <= 0), 'filled');
contour(temX, temY, temZ, 'showtext', 'on')

% % Generate
% pgon1 = nsidedpoly(8, 'Center', [0 0], 'SideLength', 10);
% newVert = zeros(8, 2);
% newVert(:, 1) = pgon1.Vertices(:, 1) * cos(-pi/2) - pgon1.Vertices(:, 2) * sin(-pi/2);
% newVert(:, 2) = pgon1.Vertices(:, 1) * sin(-pi/2) + pgon1.Vertices(:, 2) * cos(-pi/2) ;
% pgon1.Vertices = newVert;
% pp = plot(pgon1);
% axis equal
% pp.FaceAlpha = 1;
% pp.FaceColor = 'k';
% pp.EdgeColor = 'none';
% axis off
% set(gcf,'color','w');
% print(gcf, '-dpng', '-r600', 'octogon');