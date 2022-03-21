function lsVal = LSFunc_07_Octogon(width, height, cx, cy, gridX, gridY)
r = height/2;
% pi/8: cos, sin, tan.
k = [-0.9238795325,...  % sqrt(2 + sqrt(2))/2 
     0.3826834323,...   % sqrt(2 - sqrt(2))/2
     0.4142135623 ];    % sqrt(2) - 1
 
lsVal = zeros(1, (gridX + 1) * (gridY + 1));
xp = (0 : gridX) - cx; yp =  (0 : gridY) - cy;
index = 0;

for jj = 1 : length(yp)
    thisY = abs(yp(jj));
    for ii = 1 : length(xp)
        thisX = abs(xp(ii));
        index = index + 1;
        
        p = [thisX, thisY] - 2 * min(dot([k(1), k(2)], [thisX, thisY]), 0.0)*[k(1), k(2)]; 
        p = p - 2 * min(dot([-k(1), k(2)], p), 0.0)*[-k(1), k(2)];      
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