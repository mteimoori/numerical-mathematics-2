function [phi] = a11e02matrixelement(i,j,k,l,h,order)
    if(order == 1)
        phi = quad2d(dot(a11e02pyramid(i,j,h,order),a11e02pyramid(k,l,h,order)),0,1,0,1);
    elseif(order == 0)
        phi = sqrt(sum(abs(dot(a11e02pyramid(i,j,h,order),a11e02pyramid(k,l,h,order))).^2));
    end
end