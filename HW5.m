T=1;
epsilon=0.0008;
N=512;

h=T/N;
x=h:h:1;
V=10;

figure
U=exp(-25*(x-0.5).^2).*exp(-sqrt(-1)/5/epsilon.*log(exp(5*(x-0.5))+exp(-5*(x-0.5))));
t=0.54;
Uhat=zeros(1,N);
K=-N/2+1:N/2;
for k=1:N
    Uhat(k) = h*FFT(U,K(k)).*exp(-sqrt(-1).*K(k).^2*epsilon/2*t);
end
U1=zeros(1,N);
for j = 1:N
    U1(j) = IFFT([Uhat(N/2+1:N),Uhat(1:N/2)],j)/T;
end

U=U1.*exp(-sqrt(-1)*V/epsilon*t/2);

plot(x,abs(U).^2)
