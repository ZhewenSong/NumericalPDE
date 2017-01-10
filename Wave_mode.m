h=1e-1; % mesh
k=5e-2; % timestep

v=0;
x=(-3:h:3)'; I=size(x,1);
y=(-3:h:3)'; J=size(y,1);
N=500;
U=zeros(I,J);
U1=zeros(I,J);
U2=zeros(I,J);
for i=1:I
    for j=1:J
        if x(i)^4+y(j)^4>2^4 %x(i)^2+y(j)^2>2^2
            U1(i,j)=0;
            U2(i,j)=0;
        else
            U1(i,j)=exp(-pi*(x(i)^2+y(j)^2));
            U2(i,j)=k*v+U1(i,j);
        end
    end
end

[X,Y]=meshgrid(x);
mesh(X,Y,U)
zlim([-1,1])
xlim([-2,2])
ylim([-2,2])
axis tight
set(gca,'nextplot','replacechildren','visible','off')
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,20) = 0;

for n=1:N
    for i=2:I-1
        for j=2:J-1
            if abs(x(i))+abs(y(j))>2 %x(i)^2+y(j)^2>2^2
                U(i,j)=0;
            else
                U(i,j)=k^2/h^2*(U1(i+1,j)+U1(i-1,j)+U1(i,j+1)+U1(i,j-1)-4*U1(i,j))+2*U1(i,j)-U2(i,j);
            end
                        
                %U(i,j,t)=mu*(U(i+1,j,t-1)+U(i-1,j,t-1)+U(i,j+1,t-1)+U(i,j-1,t-1)-4*U(i,j,t-1))+U(i,j,t-1);
        end
    end
    U2=U1;
    U1=U;

    mesh(X,Y,U);
    zlim([-1,1])
    xlim([-2,2])
    ylim([-2,2])
    f = getframe;
    im(:,:,1,n) = rgb2ind(f.cdata,map,'nodither');

end
imwrite(im,map,'DancingPeaks.gif','DelayTime',0,'LoopCount',inf)
