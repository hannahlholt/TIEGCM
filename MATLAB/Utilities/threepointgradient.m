function [dydx] = threepointgradient(x,y,linear)
%If "linear" equals 1, then the scale height is calculated with the 
%exponential technique. If "linear" = 0 then dy/dx is returned with a three
%point gradient technique.

points = length(x);
dydx = zeros(1,points);

if linear == 0
for i=1:points
    if i==1 %First Point gradient technique
        coeff1 = (2*x(1)-x(2)-x(3))/((x(1)-...
            x(2))*(x(1)-x(3)));
        
        coeff2 = (2*x(1)-x(1)-x(3))/((x(2)-...
            x(1))*(x(2)-x(3)));
        
        coeff3 = (2*x(1)-x(1)-x(2))/((x(3)-...
            x(1))*(x(3)-x(2)));
        
        dydx(1) = (y(1)*coeff1+y(2)*coeff2+...
            y(3)*coeff3);
        
        
    elseif i==points %Last point gradient technique
        coeff1 = (2*x(i)-x(i-1)-x(i))/((x(i-2)-...
            x(i-1))*(x(i-2)-x(i)));
        
        coeff2 = (2*x(i)-x(i-2)-x(i))/((x(i-1)-...
            x(i-2))*(x(i-1)-x(i)));
        
        coeff3 = (2*x(i)-x(i-2)-x(i-1))/((x(i)-...
            x(i-2))*(x(i)-x(i-1)));
        
        dydx(i) = (y(i-2)*coeff1+y(i-1)*coeff2+...
            y(i)*coeff3);

        
    else %Middle Points gradient technique
        coeff1 = (2*x(i)-x(i)-x(i+1))/((x(i-1)-...
            x(i))*(x(i-1)-x(i+1)));
        
        coeff2 = (2*x(i)-x(i-1)-x(i+1))/((x(i)-...
            x(i-1))*(x(i)-x(i+1)));
        
        coeff3 = (2*x(i)-x(i-1)-x(i))/((x(i+1)-...
            x(i-1))*(x(i+1)-x(i)));
        
        dydx(i) = (y(i-1)*coeff1+y(i)*coeff2+...
            y(i+1)*coeff3);
    end
end
end

if linear==1
    % -----Scale Height-----
    lnn_bf = reallog(y);
    dydx = zeros(points,1);
    for i=1:points
        if i==1
           dydx(i)=-(lnn_bf(2)-lnn_bf(1))./...
                            (x(2)-x(1));
        elseif i==points
            dydx(i)=-(lnn_bf(i)-lnn_bf(i-1))./...
                            (x(i)-x(i-1));
        else
            dydx(i)=-(lnn_bf(i+1)-lnn_bf(i-1))./...
                            (x(i+1)-x(i-1));
        end
    end
    dydx=1./dydx;
end

end

