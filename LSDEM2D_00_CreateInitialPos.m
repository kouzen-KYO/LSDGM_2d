clc
clear

numX = 30; numY = 10;
numAll = numX * numY;
pos = zeros(numAll, 2);
rot = zeros(numAll, 1);
bd = [30 1000 30 300];
x = linspace(bd(1, 1), bd(1, 2), numX);
y = linspace(bd(1, 3), bd(1, 4), numY);
index = 0;
for ii = 1 : 1 : numY
    yp = y(ii);
    for jj = 1 : 1 : numX
        index = index + 1;
        xp = x(jj);
        pos(index, :) = [xp yp];
    end 
end
rowrank = randperm(size(pos, 1));
pos = pos(rowrank, :);

fid = fopen('positions_caicos001.dat', 'w');
for ii = 1 : numAll
    fprintf(fid, '%.6f %.6f\n', pos(ii, 1), pos(ii, 2));
end
fclose(fid);
save('rotations_caicos001.dat', '-ascii', 'rot');