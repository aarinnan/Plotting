function h = arrow( coor, sz)
%h = arrow( coor, sz)
% Draws an arrow through the coordinates given in 'coor' in the current
% axes
%
%INPUT:
% coor  At least 2 by 2 with starting and ending position given in columns.
%        Optionally it can be longer, and then the arrows goes through all
%        the given coordinates
% sz    Size of the arrow heads. Should be two numbers (one for the 
%        beginning, and one for the end). If no arrowhead is wanted, set
%        the value to 'NaN'
%
%OUTPUT:
% h     The handle of the arrow, making it possible to change it's looks

%AAR 200529

%Get the limits of the axis in order to know how the sizes of the arrows
%should be made
lim = axis;
rn = [ lim(2) - lim(1) lim(4) - lim(3)];

if size( coor, 1) == 3 %size( coor, 1) > 2
    if size( coor, 1) == 3
        coor_old = coor;
        temp = [ coor( 1, 1)^2 - coor( 2, 1)^2 + coor( 1, 2)^2 - coor( 2, 2)^2 2 * ( coor( 2, 1) - coor( 1, 1) ) 2 * ( coor( 2, 2) - coor( 1, 2) );
            coor( 1, 1)^2 - coor( 3, 1)^2 + coor( 1, 2)^2 - coor( 3, 2)^2 2 * ( coor( 3, 1) - coor( 1, 1)) 2 * ( coor( 3, 2) - coor( 1, 2) )];
        temp( :, 1) = -temp( :, 1);
        temp( 1, :) = temp( 1, :)/ temp( 1, 2);
        temp( 2, :) = temp( 2, :) - temp( 1, :) * temp( 2, 2);
        temp( 2, :) = temp( 2, :)/ temp( 2, 3);
        temp( 1, :) = temp( 1, :) - temp( 2, :) * temp( 1, 3);
        xm = temp( :, 1);
        r = sqrt( ( coor( 1, 1) - xm(1))^2 + ( coor( 1, 2) - xm(2))^2);
        if abs( coor( 3, 2) - coor( 1, 2)) > abs( coor( 3, 1) - coor( 1, 1) )
            yax = linspace( coor( 1, 2), coor( 3, 2), 100);
            new = r^2 - ( yax - xm(2)).^2;
            xax = xm(1) + sqrt( new);
        else
            xax = linspace( coor( 1, 1), coor( 3, 1), 100);
            new = r^2 - ( xax - xm(1)).^2;
            yax = xm(2) + sqrt( new);
        end
        if yax(1) ~= coor( 1, 2)
            yax = xm(2) - sqrt( new);
        end
        h(1) = plot( xax, yax, '-');
        temp = [xax/ rn(1); yax/ rn(2)];
        temp = temp - temp( :, 1);
        temp = sqrt( temp( 1, :).^2 + temp( 2, :).^2);
        st = find( temp > sz(1), 1);
        coor = [ xax(1) yax(1); xax( st) yax(st)]; 
        temp = [xax/ rn(1); yax/ rn(2)];
        temp = temp - temp( :, end);
        temp = sqrt( temp( 1, :).^2 + temp( 2, :).^2);
        st = find( temp < sz(1), 1);
        coor = [ coor; xax( st) yax( st); xax( end) yax( end)];
%         coor = [ vec( xax) vec( yax)];
    else
        error( 'It has not (yet) been implemented for stranger arrows')
    end
else
    h(1) = plot( coor( :, 1), coor( :, 2), '-');
end

%The arrowheads
%The first arrowhead
n = 2;
if ~isnan( sz(1))
    %Calculate the angle of the arrow
    %Make them into normalized units
    y = ( coor( 2, 2) - coor( 1, 2))/ rn(2);
    x = ( coor( 2, 1) - coor( 1, 1))/ rn(1);
    u = [x 0]'; v = [x y]';
    t = acos( ( u' * v)/ (sqrt( u'*u) * sqrt( v' * v)) );
    if isnan( t)
        if u' * v == 0
            if y > 0
                t = pi/2;
            elseif y < 0
                t = -pi/2;
            end
        end
    end
    if x < 0 && y > 0
        t = pi/ 2 + ( pi/2 - t);
    end
    q = cos( t + [-.25 .25]);
    w = sqrt( sz(1)^2./q.^2 - sz(1)^2) .* sign( q);
    cor = ( sz(1)./ sqrt( sz(1)^2 + w.^2) ) .* sign( q);
    if t == -pi/2
        cor = -cor;
    end
    h(n) = plot( coor( 1, 1) + [0 sz(1) * cor(1) * rn(1)], coor( 1, 2) + [0 cor(1) * w(1) * rn(2)], '-');
    n = n + 1;
    h(n) = plot( coor( 1, 1) + [0 sz(1) * cor(2) * rn(1)], coor( 1, 2) + [0 cor(2) * w(2) * rn(2)], '-');
    n = n + 1;
end

%The second arrowhead
if ~isnan( sz(2))
    y = ( coor( end, 2) - coor( end - 1, 2))/ rn(2);
    x = ( coor( end, 1) - coor( end - 1, 1))/ rn(1);
    u = [x 0]'; v = [x y]';
    t = acos( ( u' * v)/ (sqrt( u'*u) * sqrt( v' * v)) );
    if isnan( t)
        if u' * v == 0
            if y > 0
                t = pi/2;
            elseif y < 0
                t = -pi/2;
            end
        end
    end    
    if x < 0 && y > 0
        t = pi/ 2 + ( pi/2 - t);
    end
    q = cos( t + [-.25 .25]);
    w = sqrt( sz(2)^2./q.^2 - sz(2)^2) .* sign( q);
    cor = ( sz(2)./ sqrt( sz(2)^2 + w.^2)) .* sign( q);
    if t == -pi/2
        cor = -cor;
    end
    h(n) = plot( coor( end, 1) - [0 sz(2) * cor(1) * rn(1)], coor( end, 2) - [0 cor(1) * w(1) * rn(2)], '-');
    n = n + 1;
    h(n) = plot( coor( end, 1) - [0 sz(2) * cor(2) * rn(1)], coor( end, 2) - [0 cor(2) * w(2) * rn(2)], '-');
    n = n + 1;
end
    