k=1e-5;
h=1e-2;
x=(0:h:1)';
y=(-10:h:10)';
I=size(x,1);
J=size(y,1);
gamma=1.0;
alpha=2.0;
N=200;
phi=zeros(I,J);
for i=1:I
    for j=1:J
        phi(i,j)=rand()*2-1;
    end
end
[X,Y]=meshgrid(x,y);
%f = getframe;
%[im,map] = rgb2ind(f.cdata,256,'nodither');
%im(1,1,1,20) = 0;
for n=1:N
    for i=0:I-1
        for j=0:J-1
            dx23=(phi(mod(i+1,I)+1,j+1)^3-2*phi(i+1,j+1)^3+phi(mod(i-1,I)+1,j+1)^3)/h^2;
            dy23=(phi(i+1,mod(j+1,J)+1)^3-2*phi(i+1,j+1)^3+phi(i+1,mod(j-1,J)+1)^3)/h^2;
            dx2=(phi(mod(i+1,I)+1,j+1)-2*phi(i+1,j+1)+phi(mod(i-1,I)+1,j+1))/h^2;
            dy2=(phi(i+1,mod(j+1,J)+1)-2*phi(i+1,j+1)+phi(i+1,mod(j-1,J)+1))/h^2;
            dx4=(phi(mod(i+2,I)+1,j+1)-4*phi(mod(i+1,I)+1,j+1)+6*phi(i+1,j+1)-4*phi(mod(i-1,I)+1,j+1)+phi(mod(i-2,I)+1,j+1))/h^4;
            dy4=(phi(i+1,mod(j+2,J)+1)-4*phi(i+1,mod(j+1,J)+1)+6*phi(i+1,j+1)-4*phi(i+1,mod(j-1,J)+1)+phi(i+1,mod(j-2,J)+1))/h^4;
            phi(i+1,j+1)=phi(i+1,j+1)+k*(alpha*(dx23+dy23)-(dx2+dy2)-gamma*(dx4+dy4));
        end
    end
    mesh(X,Y,phi)
    axis equal
    axis off
    view(0,90)
    zlim([0,1])
    xlim([-10,10])
    ylim([-10,10])
    %f = getframe;
    %im(:,:,1,n) = rgb2ind(f.cdata,map,'nodither');
end
%imwrite(im,map,'PhaseSeparation.gif','DelayTime',0,'LoopCount',inf)

