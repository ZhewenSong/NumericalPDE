
figure
K=[1e-2,1e-3,1e-4]; % timesteps
% Backward Euler
subplot(3,2,1)
err=[];
for k=K
t=0:k:1;
y=zeros(size(t));
y(1)=1;
for n=1:size(t,2)-1
    y(n+1)=(100*k*t(n+1)^3+3*k*t(n+1)^2+y(n))/(1+100*k);
end
plot(t,log10(abs(y-(t.^3+exp(-100.*t))))); hold on
err=[err log10(abs(y(n)-(t(n)^3+exp(-100*t(n)))))];
end


title('Backward Euler')
legend('k=1e-2','k=1e-3','k=1e-4')
xlabel('t')
ylabel('lg|Error|')
subplot(3,2,2)
plot(log10(K),err,'o')
xlabel('lg(k)')
ylabel('lg|Error|')

% Trapezoidal
subplot(3,2,3)
err=[];
for k=K
t=0:k:1;
y=zeros(size(t));
y(1)=1;
for n=1:size(t,2)-1
    y(n+1)=(y(n)+k/2*(100*t(n+1)^3+3*t(n+1)^2 ...
    +100*t(n)^3-100*y(n)+3*t(n)^2))/(1+50*k);
end
plot(t,log10(abs(y-(t.^3+exp(-100.*t))))); hold on
err=[err log10(abs(y(n)-(t(n)^3+exp(-100*t(n)))))];
end
legend('k=1e-2','k=1e-3','k=1e-4')
title('Trapezoidal')
xlabel('t')
ylabel('lg|Error|')
subplot(3,2,4)
plot(log10(K),err,'o')
xlabel('lg(k)')
ylabel('lg|Error|')
% RK4
subplot(3,2,5)
err=[];
for k=K
t=0:k:1;
y=zeros(size(t));
y(1)=1;
for n=1:size(t,2)-1
    Y1=y(n);
    Y2=y(n)+k/2*(100*(t(n)^3-Y1)+3*t(n)^2);
    Y3=y(n)+k/2*(100*((t(n)+k/2)^3-Y2)+3*(t(n)+k/2)^2);
    Y4=y(n)+k*(100*((t(n)+k/2)^3-Y3)+3*(t(n)+k/2)^2);
    y(n+1)=y(n)+k/6*((100*(t(n)^3-Y1)+3*t(n)^2) ...
    +2*(100*((t(n)+k/2)^3-Y2)+3*(t(n)+k/2)^2) ...
    +2*(100*((t(n)+k/2)^3-Y3)+3*(t(n)+k/2)^2) ...
    +(100*((t(n)+k)^3-Y4)+3*(t(n)+k)^2));
end
plot(t,log10(abs(y-(t.^3+exp(-100.*t))))); hold on
err=[err log10(abs(y(n)-(t(n)^3+exp(-100*t(n)))))];
end
legend('k=1e-2','k=1e-3','k=1e-4')
title('RK4')
xlabel('t')
ylabel('lg|Error|')
subplot(3,2,6)
plot(log10(K),err,'o')
xlabel('lg(k)')
ylabel('lg|Error|')


