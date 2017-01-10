T=1;
N=16;
h=T/N;
x=h:h:T;
K=-N/2+1:N/2;


%u=exp(-25*(x-0.5).^2).*exp(-sqrt(-1)/5/epsilon.*log(exp(5*(x-0.5))+exp(-5*(x-0.5))));
u=sin(x*2*pi/T);
uhat=[];
for k=1:N
    uhat(k)=0;
    for j=1:N
        uhat(k)=uhat(k)+h*exp(-sqrt(-1)*2*pi/T*K(k)*x(j))*u(j);
    end
end
%{

uhat=[];
for k=1:N
    uhat(k)=FFT(u,k)/N;
end

%}
unew=[];
for j=1:N
    unew(j)=0;
    for k=1:N
        unew(j)=unew(j)+exp(sqrt(-1)*2*pi/T*K(k)*x(j))*uhat(k);
    end
    %unew(j)=unew(j)+(exp(-sqrt(-1)*K(N)*x(j))+exp(sqrt(-1)*K(N)*x(j)))*uhat2(N)/T/2;
end

%{
unew=[];
for j=1:N
    unew(j)=IFFT(uhat,j);%+uhat(N)*(exp(-sqrt(-1)*N/2*x(j))-exp(sqrt(-1)*N/2*x(j)))/2;
    
end
%}
plot(x,real(u),x,real(unew))
