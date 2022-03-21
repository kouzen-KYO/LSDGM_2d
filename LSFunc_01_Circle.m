function lsVal = LSFunc_01_Circle(width, height, cx, cy, gridX, gridY)

ra = width/2; rb = height/2;
r = (ra + rb)/2;
lsVal = zeros(1, (gridX + 1) * (gridY + 1));
xp = (0 : gridX) - cx; yp =  (0 : gridY) - cy;

index = 0;
for jj = 1 : length(yp)
    thisY = yp(jj);
    for ii = 1 : length(xp)
        thisX = (xp(ii));
        p = [thisX, thisY];
        index = index + 1;
        
        lsVal(index) = norm(p) - r;
    end

end