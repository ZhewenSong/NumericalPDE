function out = u_hat(u_l,u_r,vbar,x,t)
    r1=[1,vbar-1]';
    r2=[1,vbar+1]';
    if x/t<=vbar-1
        out=u_l;
    elseif x/t<vbar+1
        out=u_l+(r1'*(u_r-u_l))*r1;
    else
        out=u_r;
    end
end

