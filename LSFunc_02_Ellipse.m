function lsVal = LSFunc_02_Ellipse(width, height, cx, cy, gridX, gridY)
a = width; b = height;
ab = [a b];
lsVal = zeros(1, (gridX + 1) * (gridY + 1));
xp = (0 : gridX) - cx; yp =  (0 : gridY) - cy;

index = 0;
for jj = 1 : length(yp)
    thisY = yp(jj);
    for ii = 1 : length(xp)
        thisX = (xp(ii));
        p = [thisX, thisY];
        p = abs(p);
        index = index + 1;
 
        if p(1) > p(2)
            p = [p(2) p(1)]; ab = [ab(2) ab(1)];
        end
        
        l = ab(2) * ab(2) - ab(1) * ab(1);
        
        m = ab(1) * p(1)/l;
        n = ab(2) * p(2)/l;
        m2 = m * m;
        n2 = n * n;
        
        c = (m2 + n2 - 1.0)/3.0;
        c3 = c * c * c;
        
        d = c3 + m2 * n2;
        q = d  + m2 * n2;
        g = m  + m * n2;
        
        if d < 0.0
            h = acos(q/c3)/3.0;
            s = cos(h) + 2.0;
            t = sin(h) * sqrt(3.0);
            rx = sqrt( m2 - c*(s + t) );
            ry = sqrt( m2 - c*(s - t) );
            co = ry + sign(l)*rx + abs(g)/(rx*ry);
        else
            h = 2 * m * n * sqrt(d);
            s = sign(q + h) * abs(q + h)^(1/3);
            t = sign(q - h) * abs(q - h)^(1/3);
            rx = -(s+t) - c*4.0 + 2.0*m2;
            ry =  (s-t)*sqrt(3.0);
            rm = sqrt( rx*rx + ry*ry );
            co = ry/sqrt(rm-rx) + 2.0*g/rm;
        end
        co = (co - m)/2.0;
        
        si = sqrt(max(1.0 - co*co, 0.0));
        r = ab.* [co, si];
        
        lsVal(index) = norm(r - p) * sign(p(2) - r(2));
    end
end

end