Err=[];
for p=[5, 6, 7, 8, 9, 10]
    N=2^p;
    h=1/N;
    A=zeros(N-1,N-1); % Defining the 5-diagonal matrix A
    A(1,1:3)=[7 -4 1];
    A(N-1,N-3:N-1)=[1 -4 7];
    A(2,1:4)=[-4 6 -4 1];
    A(N-2,N-4:N-1)=[1 -4 6 -4];
    for i=3:N-3
    A(i,i-2:i+2)=[1 -4 6 -4 1];
    end

    x=[1:(N-1)]*h;
    % f evaluated at grid points
    F=exp(x).*(x.^4+14*x.^3+49*x.^2+32*x-12); 
    U=A\F'*h^4; % solving U by inverting matrix A
    g=(exp(x).*x.^2.*(1-x).^2)'; % analytical solution 
    figure
    plot(x,g,x,U);
    legend('Analytical','Numerical')
    
    Err=[Err max(abs(U-g))];
    %Err=[Err max(abs(A*U-A*g))];
end
figure
plot(log(2.^[-5, -6, -7, -8, -9, -10]),log(Err),'o')
xlabel('log(h)')
ylabel('log(e_h)')
