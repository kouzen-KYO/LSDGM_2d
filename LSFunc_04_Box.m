function lsVal = LSFunc_04_Box(width, height, cx, cy, gridX, gridY)

lsVal = zeros(1, (gridX + 1) * (gridY + 1));
xp = (0 : gridX) - cx; yp =  (0 : gridY) - cy;
index = 0;
for jj = 1 : length(yp)
    thisY = abs(yp(jj));
    for ii = 1 : length(xp)
        thisX = abs(xp(ii));
        index = index + 1;
        if thisX <= width/2 && thisY <= height/2
            val = -min(width/2 - thisX, height/2 - thisY);
        else
            val = sqrt(max(0, thisX - width/2)^2 + max(0, thisY - height/2)^2);
        end
        lsVal(index) = val;
    end
end
end