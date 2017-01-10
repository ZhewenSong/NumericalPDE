clear;
h=0.5e-2;
k=5e-4;
x=(-1:h:10)';
J=size(x,1);
u_LLF=zeros(J,1);
u_LW=zeros(J,1);
u_God=zeros(J,1);
u_exact=zeros(J,1);

%%%% Initial Conditions --- Prob.2 %%%%
for j=1:J
    if x(j)>0 && x(j)<1
        u_LLF(j)=2; u_LW(j)=2; u_God(j)=2; u_exact(j)=2;
    else
        u_LLF(j)=0; u_LW(j)=0; u_God(j)=0; u_exact(j)=0;
    end
end

%%%% Initial Conditions --- Prob.3 %%%%
%{
for j=1:J
    if x(j)>0 && x(j)<1
        u_LLF(j)=2; u_LW(j)=2; u_God(j)=2; u_exact(j)=2;
    else
        u_LLF(j)=4; u_LW(j)=4; u_God(j)=4; u_exact(j)=4;
    end
end
%}

figure
plt=1;
for n=1:1.5/k
    for j=2:J-1
        %%%% Local Lax-Friedrichs %%%%
        Fp=(u_LLF(j+1)^2/2+u_LLF(j)^2/2)/2-max(abs(u_LLF(j:j+1)))/2*(u_LLF(j+1)-u_LLF(j));
        Fm=(u_LLF(j-1)^2/2+u_LLF(j)^2/2)/2-max(abs(u_LLF(j-1:j)))/2*(u_LLF(j)-u_LLF(j-1));
        u_LLF(j)=u_LLF(j)-k/h*(Fp-Fm);
        
        %%%% Lax-Wendroff %%%%
        
        up=(u_LW(j)+u_LW(j+1))/2-k/2/h*(u_LW(j+1)^2/2-u_LW(j)^2/2);
        um=(u_LW(j-1)+u_LW(j))/2-k/2/h*(u_LW(j)^2/2-u_LW(j-1)^2/2);
        u_LW(j)=u_LW(j)-k/h*(up^2/2-um^2/2);
        
        %%%% Godunov %%%%
        if u_God(j-1)<=u_God(j)
            Fm=min(u_God(j-1:j).^2/2);
        else
            Fm=max(u_God(j-1:j).^2/2);
        end
        
        if u_God(j)<=u_God(j+1)
            Fp=min(u_God(j:j+1).^2/2);
        else
            Fp=max(u_God(j:j+1).^2/2);
        end        
        u_God(j)=u_God(j)-k*(Fp-Fm)/h;
        
        %%%% Exact Solution --- Prob.2 %%%%
        t=n*k;
        if t<1
            if x(j)<0
                u_exact(j)=0;
            elseif x(j)<2*t
                u_exact(j)=x(j)/t;
            elseif x(j)<1+t
                u_exact(j)=2;
            else
                u_exact(j)=0;
            end
        else
            if x(j)<0
                u_exact(j)=0;
            elseif x(j)<2*sqrt(t)
                u_exact(j)=x(j)/t;
            else
                u_exact(j)=0;
            end
        end
        
        %%%% Exact Solution --- Prob.3 %%%%
        %{
        t=n*k;
        if t<1
            if x(j)<3*t
                u_exact(j)=4;
            elseif x(j)<2*t+1
                u_exact(j)=2;
            elseif x(j)<4*t+1
                u_exact(j)=(x(j)-1)/t;
            else
                u_exact(j)=4;
            end
        else
            if x(j)<4*t+1-2*sqrt(t)
                u_exact(j)=4;
            elseif x(j)<4*t+1
                u_exact(j)=(x(j)-1)/t;
            else
                u_exact(j)=4;
            end
        end
        %}
    end
end


