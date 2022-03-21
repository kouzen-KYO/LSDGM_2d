function lsVal = LSFunc_05_Pentagon(width, height, cx, cy, gridX, gridY)
r = height/(1 + 1/cos(deg2rad(36)));
k = [0.809016994,0.587785252,0.726542528]; % pi/5 : cos, sin, tan
lsVal = zeros(1, (gridX + 1) * (gridY + 1));
xp = (0 : gridX) - cx; yp =  (0 : gridY) - cy;
index = 0;
for jj = 1 : length(yp)
    thisY = -yp(jj);
    for ii = 1 : length(xp)
        thisX = abs(xp(ii));
        index = index + 1;
        
        p = [thisX, thisY] - 2 * min(dot([-k(1), k(2)], [thisX, thisY]), 0.0)*[-k(1), k(2)];
        p = p - 2 * min(dot([k(1), k(2)], p), 0.0) * [k(1), k(2)];
        if p(1) > r * k(3)
            aa = r * k(3);
        elseif p(1) < -r * k(3)
            aa = -r * k(3);
        else
            aa = p(1);
        end
        p = p - [aa, r];
        val = norm(p) * sign(p(2));
        lsVal(index) = val;
    end
end

end