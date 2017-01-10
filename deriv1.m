function A = deriv1( A,C,i,j,k,l,I,J,n )

coeff = [1,-1]*(n-1)^2/4;
Im = I-1; Ip = I+1;
Jm = J-1; Jp = J+1;
if I==1
    Im = n;
elseif I==n
    Ip = 1;
end
if J==1
    Jm = n;
elseif J==n
    Jp = 1;
end
if l==1 
    xl = [Ip,Im];
    yl = [J,J];
else
    xl = [I,I];
    yl = [Jp,Jm];
end
if j==1 
    xj = [Ip,Im];
    yj = [J,J];
else
    xj = [I,I];
    yj = [Jp,Jm];
end
for s=1:2   
    A((i-1)*n^2+(I-1)*n+J,(k-1)*n^2+(xl(s)-1)*n+yl(s)) ...
    = A((i-1)*n^2+(I-1)*n+J,(k-1)*n^2+(xl(s)-1)*n+yl(s)) + C(xj(s),yj(s))*coeff(s);
end

end

