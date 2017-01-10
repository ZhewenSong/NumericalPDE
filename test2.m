clear;
A=[  2     1     0
     1     3    -1
     0    -1     1];
U_l=[1 2 3]'; 
U_r=[3,2,1]';
[R,L]=eig(A);
for i=1:3
    lambda(i)=L(i,i);
end

dt=0.01;
dx=0.01;
t=0:dt:2;
x=-10:dx:10;
U=zeros(3,size(x,2));
U2=zeros(3,size(x,2));

for n=1:size(t,2)
    for j=1:size(x,2)
        U(:,j)=0;
        I=[];
        for i=1:3
            if lambda(i)*t(n)<x(j)
                I=[I,i];
            end
        end
        if size(I,1)==0
            U(:,j)=U_l;
        else
            I=max(I);
            for i=1:I
                U(:,j)=U(:,j)+(R(:,i)'*U_r)*R(:,i);
            end
            for i=I+1:3
                U(:,j)=U(:,j)+(R(:,i)'*U_l)*R(:,i);
            end
        end
    end
    
    
    for j=1:size(x,2)
        if x(j)/t(n)<lambda(1)
            U2(:,j)=U_l;
        elseif x(j)/t(n)<lambda(2)
            U2(:,j)=U_l+(R(:,1)'*(U_r-U_l))*R(:,1);
        elseif x(j)/t(n)<lambda(3)
            U2(:,j)=U_l+(R(:,1)'*(U_r-U_l))*R(:,1)+(R(:,2)'*(U_r-U_l))*R(:,2);
        else
            U2(:,j)=U_r;
        end
    end
    plot(x,U(1,:),'-',x,U2(1,:),':')
   
    shg
    pause(dt)
end
            
       
                