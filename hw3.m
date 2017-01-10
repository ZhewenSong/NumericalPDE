figure
k = 1e-9;
H = [1e-1,1e-2,1e-3,1e-4];
err=[];
for h = H
x = (h:h:1-h)';
N = size(x,1);
B = zeros(N);
B(1,1) = 1 - 2*k/h^2;
B(1,2) = k/h^2;
B(N,N-1) = k/h^2;
B(N,N) = 1 - 2*k/h^2;
for i=2:N-1
    B(i,i) = 1 - 2*k/h^2;
    B(i,i-1) = k/h^2;
    B(i,i+1) = k/h^2;
end
u = 0.5*sin(6*pi*x) + sin(pi*x);

T = 10;
U = zeros(N,T);
U_exact = zeros(N,T);
for t=1:T
    u = B*u;
    U(:,t) = u;
    U_exact(:,t) = 0.5*exp(-36*pi^2*k*t).*sin(6*pi*x) ...
        + exp(-pi^2*k*t).*sin(pi*x);
    %plot([0;x;1],[0;u;0],[0;x;1],[0;U_exact(:,t);0])
    %pause(0.1);
end

err = [err max(max(abs(U_exact-U)))];
end
scatter(log10(H),log10(err))