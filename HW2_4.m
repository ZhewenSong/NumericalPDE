h=1e-2; % mesh
k=1e-4; % timestep
epsilon=0.05;
lambda=k/h^2;

figure
x=(h:h:1-h)';
y=(h:h:1-h)';
I=size(x,1);
J=size(y,1);
N=2000;
U=zeros(I,J);
for i=1:I 
    for j=1:J
        if 4*(x(i)-0.5)^2+16*(y(j)-0.5)^2<=1 % Initial condition
            U(i,j)=1;
        else
            U(i,j)=-1;
        end
    end
end

U0=-1; % Boundary condition

a=-lambda*ones(J,1);
b=2*(1+lambda)*ones(J,1);
c=-lambda*ones(J,1);
a(1)=0; c(J)=0; 
alpha=k/4/epsilon;

for n=1:N
    
for i=1:I  % Step 1: U^\star=C(U^n,k/2)
    for j=1:J
        u=0;
        delta=(alpha*u^3+(1-alpha)*u-alpha*U(i,j)*(1-U(i,j)^2)-U(i,j))...
            /(3*alpha*u^2+1-alpha);
        while abs(delta)>1e-6  % Newton's iterative method for root finding
            u=u-delta;
            delta=(alpha*u^3+(1-alpha)*u-alpha*U(i,j)*(1-U(i,j)^2)-U(i,j))...
                /(3*alpha*u^2+1-alpha);
        end
        U(i,j)=u;
    end
end

for i=1:I  % Step 2: U^\star\star=B(U^\star,k)
    d=zeros(J,1);
    d(1)=lambda*U(i,2)+2*(1-lambda)*U(i,1)+lambda*U0-(-lambda)*U0;
    for j=2:J-1
        d(j)=lambda*U(i,j+1)+2*(1-lambda)*U(i,j)+lambda*U(i,j-1);
    end
    d(J)=2*(1-lambda)*U(i,J)+lambda*U(i,J-1)+lambda*U0-(-lambda)*U0;
    U(i,:)=TridiagonalMatrixSolver(a,b,c,d);
end

for j=1:J % Step 3: U^\star\star\star=A(U^\star\star,k)
    d=zeros(I,1);
    d(1)=lambda*U(2,j)+2*(1-lambda)*U(1,j)+lambda*U0-(-lambda)*U0;
    for i=2:I-1
        d(i)=lambda*U(i+1,j)+2*(1-lambda)*U(i,j)+lambda*U(i-1,j);
    end
    d(I)=2*(1-lambda)*U(I,j)+lambda*U(I-1,j)+lambda*U0-(-lambda)*U0;
    U(:,j)=TridiagonalMatrixSolver(a,b,c,d);
end

for i=1:I % Step 4: U^{n+1}=C(U^\star\star\star,k/2)
    for j=1:J
        u=0;
        delta=(alpha*u^3+(1-alpha)*u-alpha*U(i,j)*(1-U(i,j)^2)-U(i,j))...
            /(3*alpha*u^2+1-alpha);
        while abs(delta)>1e-6
            u=u-delta;
            delta=(alpha*u^3+(1-alpha)*u-alpha*U(i,j)*(1-U(i,j)^2)-U(i,j))...
                /(3*alpha*u^2+1-alpha);
        end
        U(i,j)=u;
    end
end

U_whole=-1*ones(I+2,J+2); 
U_whole(2:I+1,2:J+1)=U;
%mesh((0:h:1),(0:h:1),U_whole);
%zlim([-1,1])
%shg;
%pause(k)

if ~mod(n,50)
    contour((0:h:1),(0:h:1),U_whole,[0,0])
    axis equal
    hold on
end

end
