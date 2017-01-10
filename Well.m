h=0.1;
x=(h:h:100)';
l = 0;
N=size(x,1);
V = @(r) 1/r; 
A=zeros(N,N);
%{
A(1,1) = 1/h^2+V(x(1));
A(1,2) = -1/2/h^2;
A(N,N) = 1/h^2+V(x(N));
A(N-1,N) = -1/2/h^2;
for i=2:N-1
    A(i,i-1)=-1/2/h^2;
    A(i,i)=1/h^2+V(x(i));
    A(i,i+1)=-1/2/h^2;
end
%}

for i=1:N
    if i>1 
        A(i,i-1) = -1/h^2 + 1/x(i)/h; 
    end
    if i<N 
        A(i,i+1) = -1/h^2 - 1/x(i)/h; 
    end
    A(i,i) = 2/h^2 - 1/x(i) + l*(l+1)/x(i)^2;
end
[eigfunc,val] = eig(A);
eigval = zeros(N,1);
for i=1:N
    eigval(i) = val(i,i);
end
[eigval,index] = sort(eigval);
eigfunc = eigfunc/sqrt(h);
plot(x,eigfunc(:,index(3)));