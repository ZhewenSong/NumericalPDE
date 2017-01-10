function U = TridiagonalMatrixSolver(a,b,c,d)
   J=size(a,1);
   cc=zeros(J-1,1);
   cc(1)=c(1)/b(1);
   for j=2:J-1
       cc(j)=c(j)/(b(j)-a(j)*cc(j-1));
   end
   dd=zeros(J,1);
   dd(1)=d(1)/b(1);
   for j=2:J
       dd(j)=(d(j)-a(j)*dd(j-1))/(b(j)-a(j)*cc(j-1));
   end
   U=zeros(J,1);
   U(J)=dd(J);
   for j=J-1:-1:1
       U(j)=dd(j)-cc(j)*U(j+1);
   end