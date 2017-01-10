err=[];
h=1e-2; 
k=1e-3; % timestep
x=(h:h:1-h)';
J=size(x,1);
N=200;
U=(10.*sin(pi*x)); % Initial condition
lambda=k/h^2;

a=-lambda*ones(J,1);
b=2*(1+lambda)*ones(J,1);
c=-lambda*ones(J,1);
a(1)=0; c(J)=0; % Boundary condition

for n=1:N
   d=zeros(J,1);
   d(1)=lambda*U(2)+2*(1-lambda)*U(1);
   for j=2:J-1
       d(j)=lambda*U(j+1)+2*(1-lambda)*U(j)+lambda*U(j-1);
   end
   d(J)=2*(1-lambda)*U(J)+lambda*U(J-1);
   U=TridiagonalMatrixSolver(a,b,c,d);
   
   if ~mod(n,50)
       plot(x,U,'r-',x,10*exp(-pi^2*n*k).*sin(pi*x),'b--');hold on
   end

end
