function U = a07e02LaxFrie( t,xmin,xmax,v,f,h,k )

l=t/k;
lamda=k/h;
xh=(xmin-l:h:xmax+l+1)';

U = v(xh);


while l > 0

    for i=2:size(U)-1
       
          U(i) = 0.5 .* (U(i+1)+U(i-1)) - lamda/2 .* (f(U(i+1)) - f(U(i-1))); 
          U(1) = 0.5 .* (U(2)+0) - lamda/2 .* (f(U(i+1)) - 0);
          U(end) = 0.5 .* (0+U(i-1)) - lamda/2 .* (0 - f(U(end-1)));
          
    end
l = l-1;
end

end

