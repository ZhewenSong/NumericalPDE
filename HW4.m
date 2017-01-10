clear;clc;
h=0.02;
k=h*0.1;
x=(-1:h:1)';
J=size(x,1);
U1=zeros(1,J);
U2=zeros(1,J);
%%%% Initial Conditions  %%%%
for j=1:J
    if x(j)>=0 
        U1(j)=0.2; U2(j)=0; 
    else
        U1(j)=1.0; U2(j)=0;
    end
end

LLF_scheme(U1,U2,x,J,k,h);

k=h*0.1;
Roe_scheme(U1,U2,x,J,k,h);

title('\Delta x=0.002; \Delta t_{LLF}=0.0002, \Delta t_{Roe}=0.0001')
leg=legend('Density (LLF)','Velocity (LLF)','Density (Roe)','Velocity (Roe)','Location','northwest');
set(leg,'FontSize',13)