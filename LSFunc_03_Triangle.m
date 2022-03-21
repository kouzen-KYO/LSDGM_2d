function lsVal = LSFunc_03_Triangle(width, height, cx, cy, gridX, gridY)
k = sqrt(3);
r = width/2;
lsVal = zeros(1, (gridX + 1) * (gridY + 1));
xp = (0 : gridX) - cx; yp =  (0 : gridY) - cy;
index = 0;
for jj = 1 : length(yp)
    thisY = yp(jj) + r/k;
    for ii = 1 : length(xp)
        thisX = abs(xp(ii)) - r;
        p = [thisX, thisY];
        index = index + 1;
        if p(1) + k*p(2) > 0
            p = [p(1) - k*p(2), -k * p(1) - p(2)]/2.0;
        end
        if p(1) > 0
            aa = 0;
        elseif p(1) < -2*r
            aa = -2*r;
        else
            aa = p(1);
        end
        p(1) = p(1) - aa;
        val = -norm(p) * sign(p(2));
        lsVal(index) = val;
    end
end

end