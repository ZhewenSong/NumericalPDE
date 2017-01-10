function A = deriv2( A,C,i,j,k,l,I,J,n )


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
if l==1 && j==1 % d^2/dx^2
    x = [Im,I,I,Ip];
    y = [J,J,J,J];
    coeff = [1,-1,-1,1]*(n-1)^2;
elseif l==2 && j==2 % d^2/dy^2
    x = [I,I,I,I];
    y = [Jm,J,J,Jp];
    coeff = [1,-1,-1,1]*(n-1)^2;
else % d^2/dxdy
    x = [Im,Im,Ip,Ip];
    y = [Jm,Jp,Jm,Jp];
    coeff = [1,-1,-1,1]*(n-1)^2/4;
end

for s=1:4   
    A((i-1)*n^2+(I-1)*n+J,(k-1)*n^2+(x(s)-1)*n+y(s)) ...
    = A((i-1)*n^2+(I-1)*n+J,(k-1)*n^2+(x(s)-1)*n+y(s)) + C(I,J)*coeff(s);
end


end

