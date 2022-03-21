clc
close all

savePath = "";
writeOpt = 2; % 1: Cover; 2: Not cover
position = [cx cy];
velocity = zeros(1, 2);
rotation = 0;
omega = zeros(1, 1);
lsCOM = [cx cy];
[C, boundRad] = MATLAB_minboundcircle(surPtsX, surPtsY);
boundRad = boundRad * 1.1; % Safety factor, depende on the error between real shape and approximated shape.
gridN = [xmax + 1 ymax + 1];

kn = 30000.000000; kt = 27000.000000; mu = 0.100000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================  WRITE .dat FILE FOR CALCULATION ======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = "caicos001.dat";
if writeOpt == 1
    fid = fopen(savePath + filename, 'w');
    fprintf(fid, '%d\n', 1);
elseif writeOpt == 2
    fid = fopen(savePath + filename, 'at');
    %fprintf(fid, '\n');
end

for ii = 1
    % Every particle has 15 rows for information
    % 01 -- Mass
    fprintf(fid, '%.6f\n', mass(ii));
    % 02 -- Actual COM !!!
    fprintf(fid, '%.6f %.6f\n', position(ii, 1), position(ii, 2));
    % 03 -- Velocity
    fprintf(fid, '%.6f %.6f\n', velocity(ii, 1), velocity(ii, 2));
    % 04 -- Moment of inertia
    fprintf(fid, '%.6f\n', inertia(ii));
    % 05 -- Rotation
    fprintf(fid, '%.6f\n', rotation(ii));
    % 06 -- Omega
    fprintf(fid, '%.6f\n', omega(ii));
    % 07 -- COM for level-set!!!
    fprintf(fid, '%.6f %.6f\n', lsCOM(ii, 1), lsCOM(ii, 2));
    % 08 -- Number of surface points
    fprintf(fid, '%d\n', numSurPoint(ii));
    % 09 -- Positions of surface points
    for jj = 1 : numSurPoint(ii)
        fprintf(fid, '%.6f %.6f ', surPts(ii, 2*jj - 1), surPts(ii, 2*jj));
    end
    fprintf(fid, '\n');
    % 10 -- Bounding box radius
    fprintf(fid, '%.6f\n', boundRad(ii));
    % 11 -- Grid number in x and y axes
    fprintf(fid, '%d %d\n', gridN(ii, 1), gridN(ii, 2));
    % 12 -- \Phi values in every nodes
    for jj = 1 : gridN(ii, 1) * gridN(ii, 2)
        fprintf(fid, '%.6f ', lsVal(jj));
    end
    fprintf(fid, '\n');
    % 13 -- Normal stiffness
    fprintf(fid, '%.6f\n', kn);
    % 14 -- Shear stiffness
    fprintf(fid, '%.6f\n', kt);
    % 15 -- Friction constant
    fprintf(fid, '%.6f\n', mu);
end
fclose(fid);

fprintf('Model initialization completed !!!!\n');
fprintf('Model initialization completed !!!!\n');