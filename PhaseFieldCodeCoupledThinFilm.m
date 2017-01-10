%% Latest version ----> function PhaseField()
% Spinodal Decomposition code
% An unstable homogeneous solution of species A and B decomposes into two stable phases of different compositions and different structure
% In this code, only the conserved phase field variable (concentration) is considered, while the order parameter isn't.
% Elastic energy caused by lattice misfit bewteen matrix and precipitate is considered.
clear all
clc

TimeStart=fix(clock);
save TimeStart.mat TimeStart

global Nmesh average_composition e0
Nmesh=50;   % mesh size.
average_composition=0.5000;%0.2484;%0.2484; %0.1261;%;
e0=0.0; %0.05; --?    % composition expansion coef.
%C11=232.0;
%C12=153.0;   % with possion effect
%C44=116.0;  % 117.0

%% -----------------------MAIN FUNCTION BEGINS--------------------------*/
% 
%   int i,j;
%   int NT=0;
%   double C[Nmesh+4][Nmesh+4][2];
%   double noise=0.0;
%   double pro=0.0;
%   double F=0.0;
%   double freeE=0.0;
%   double gradE=0.0;
%   double miu_el[Nmesh+4][Nmesh+4]; //elastic energy term
%   double e_el[Nmesh+4][Nmesh+4]; //elastic energy term
C(:,:,1)=0.00*(2*rand(Nmesh+4)-1);
C(1:0.4*Nmesh+2,:,1)=1;
% for i=(0.4*Nmesh+3):(0.5*Nmesh+2)
%     C(i,:,1)=(-10*(i-2)/Nmesh+5)*ones(1,Nmesh+4);
% end
C(0.4*Nmesh+3:0.5*Nmesh+2,:,1)=0.2*rand(0.1*Nmesh,Nmesh+4)+0.4;
% C(:,:,1)=average_composition*zeros(Nmesh+4);
% C(:,:,1)=zeros(Nmesh+4,Nmesh+4);
% for i=1:Nmesh
%     for j=1:Nmesh
%         if (i-25)*(i-25)+(j-25)*(j-25)<200 %=Nmesh*Nmesh/2/pi
%             C(i+2,j+2,1)=1;
% %         else if (i-50)*(i-50)+(j-50)*(j-50)<400
% %                 C(i+2,j+2,1)=0.9;
% %             
% %             end
%         end
%     end
% end 

%C=CC;
miu_el=zeros(Nmesh+4,Nmesh+4);
e_el=zeros(Nmesh+4,Nmesh+4);
dx=1;
dy=1;
DX=1/(dx*dx);
dt=0.01;% 0.005; % NT=5,000:10,000
Energy_Counter=0;
counter=1;
x0=zeros(2*Nmesh*Nmesh,1);
x1=zeros(2*Nmesh*Nmesh,1);
x2=zeros(2*Nmesh*Nmesh,1);
x3=zeros(2*Nmesh*Nmesh,1);
% load x0.mat
% x0=x0;

%% Main Loop
for NT=0:1000000 %?5?0???5?0
  if mod(NT,100)==0
     process=100*NT/4000000.00;
  end
  for i=1:Nmesh+4
      ii=i;
      if(i<3)
          ii=6-i; % i+Nmesh-1; % Be careful that the period is (Nmesh-1)...
      end
      if(i>Nmesh+2)
          ii=2*Nmesh+4-i; % i-Nmesh+1;
      end
      for j=1:Nmesh+4
          if((i>=3)&&(i<=Nmesh+2)&&(j>=3)&&(j<=Nmesh+2)) 
               continue;
          end
          jj=j;
          if(j<3)
            jj=j+Nmesh-1;
          end
          if(j>Nmesh+2)
            jj=j-Nmesh+1;
          end
          C(i,j,1)=C(ii,jj,1);
      end
  end
   % h[i][j] denotes df/dc, f is the bulk free energy density. So the sole term h expresses the type of f.
   %h=4.*(C(:,:,1).*C(:,:,1)-0.9.*C(:,:,1)+0.18).*(2.*C(:,:,1)-0.9);
   %h=-2.*(C(:,:,1)-0.5)+10.*(C(:,:,1)-0.5).*(C(:,:,1)-0.5).*(C(:,:,1)-0.5);
   h=4.*C(:,:,1).*(C(:,:,1)-1).*(2.*C(:,:,1)-1);
   
   %h=50.*h; % increase initial free energy density by 50 times...
   
   u=zeros(Nmesh+4);
   v=zeros(Nmesh+4);
   uu=zeros(Nmesh+4);
   
    
%%  % /*-------------------Elastic_Equilibrium Equations BEGIN----------------------*/
  if (NT==10)||(NT==30)||(NT==50)||(NT==70)||(mod(NT,25)==0)
      w=Nmesh*Nmesh;
      Dimension=2*Nmesh*Nmesh; %dimension of the coef matrix in A*x=b.
      A=zeros(Dimension);
      b=zeros(Dimension,1);
    epsilon_f=0.1;
      for i=1:Nmesh
          
          m1=i+1;
%           if(m1>Nmesh) 
%              m1=m1-(Nmesh-1);
%           end
          m2=i-1;
%           if(m2<1) 
%              m2=m2+(Nmesh-1);
%           end
          k=(i-1)*Nmesh;
          
          for j=1:Nmesh
%               k=k+1;
              m3=j+1;
	          if(m3>Nmesh) 
                 m3=m3-(Nmesh-1);
              end
	          m4=j-1;
	          if(m4<1) 
                  m4=m4+(Nmesh-1);
              end
                            
              if(i==1)  % &&(j==1)
                  A(k+j,k+j)=1;
                  A(k+j+w,k+j+w)=1;
                  continue;
              end  

                            
              
              % Calculate elastic constant
              interp=C(i+2,j+2,1)*C(i+2,j+2,1)*C(i+2,j+2,1)*(6*C(i+2,j+2,1)*C(i+2,j+2,1)-15*C(i+2,j+2,1)+10);
              
              C11=2+102.0*interp;
              C12=1+43.5*interp;
              C44=1+29.0*interp;
              
%               if C(i+2,j+2,1)<0.01
%                   A(k+j,k+j)=1;
%                   A(k+j+w,k+j+w)=1;
%                   continue;
% %                   C11=2.0;
% %                   C12=2.0;
% %                   C44=2.0;
%               end
              
              % P denotes the derivative of interp w.r.t composition
              P=30*C(i+2,j+2,1)*C(i+2,j+2,1)*(C(i+2,j+2,1)-1)*(C(i+2,j+2,1)-1);
              
%               if(i==Nmesh)
%                   A(k+j,k+j)=1;
%                   A(k+j+w,k+j+w)=1;
%                   continue;
%               end
              if(i==Nmesh)
                   % /*-------first equation begins-------*/ u displacement
                 A(k+j,m2*Nmesh-Nmesh+j)=C11;
                 A(k+j,m2*Nmesh-Nmesh+j)=C11+A(k+j,m2*Nmesh-Nmesh+j);
                 A(k+j,i*Nmesh-Nmesh+j)=-2*C11;
                 A(k+j,i*Nmesh-Nmesh+m3)=0;
                 A(k+j,i*Nmesh-Nmesh+m4)=0;
                       % v displacement
                  A(k+j,w+m2*Nmesh-Nmesh+m3)=0;
                  A(k+j,w+m2*Nmesh-Nmesh+m4)=0;
                  A(k+j,w+m2*Nmesh-Nmesh+m3)=0+A(k+j,w+m2*Nmesh-Nmesh+m3);
                  A(k+j,w+m2*Nmesh-Nmesh+m4)=0+A(k+j,w+m2*Nmesh-Nmesh+m4);
                  A(k+j,w+i*Nmesh-Nmesh+m3)=0;
                  A(k+j,w+i*Nmesh-Nmesh+m4)=0;
                  A(k+j,w+m2*Nmesh-Nmesh+j)=50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
                  A(k+j,w+m2*Nmesh-Nmesh+j)=-50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0+A(k+j,w+m2*Nmesh-Nmesh+j);
	
                  b(k+j)=e0*(C11+C12)*(C(i+3,j+2,1)-C(i+1,j+2,1))/2.0+100*P*(C(i+3,j+2,1)-C(i+1,j+2,1))*e0*(C(i+2,j+2,1)-average_composition); %// Be cautious here...
                 % /*--------first equation ends---------*/
	
                     % /*-------second equation begins-------*/ 
                     % u displacement
                  A(w+k+j,m2*Nmesh-Nmesh+m3)=0;
                  A(w+k+j,m2*Nmesh-Nmesh+m4)=0;
                  A(w+k+j,m2*Nmesh-Nmesh+m3)=0+A(w+k+j,m2*Nmesh-Nmesh+m3);
                  A(w+k+j,m2*Nmesh-Nmesh+m4)=0+A(w+k+j,m2*Nmesh-Nmesh+m4);
                  A(w+k+j,m2*Nmesh-Nmesh+j)=0;
                  A(w+k+j,m2*Nmesh-Nmesh+j)=0+A(w+k+j,m2*Nmesh-Nmesh+j);
                  A(w+k+j,i*Nmesh-Nmesh+m3)=0;
                  A(w+k+j,i*Nmesh-Nmesh+m4)=0;
                        % v displacement
                  A(w+k+j,w+i*Nmesh-Nmesh+m3)=0;
                  A(w+k+j,w+i*Nmesh-Nmesh+m4)=0;
                  A(w+k+j,w+i*Nmesh-Nmesh+j)=-2*C44;
                  A(w+k+j,w+m2*Nmesh-Nmesh+j)=C44;
                  A(w+k+j,w+m2*Nmesh-Nmesh+j)=C44+A(w+k+j,w+m2*Nmesh-Nmesh+j);
	
                  b(w+k+j)=e0*(C12+C11)*(C(i+2,j+3,1)-C(i+2,j+1,1))/2.0+100*P*(C(i+2,j+3,1)-C(i+2,j+1,1))*e0*(C(i+2,j+2,1)-average_composition); %// Be cautious here...
                  % /*--------second equation ends--------*/
                  
                  continue;
              end
              if(i==1)
                   % /*-------first equation begins-------*/ v displacement
                 A(k+j,m1*Nmesh-Nmesh+j)=C11;
                 A(k+j,m1*Nmesh-Nmesh+j)=C11+A(k+j,m1*Nmesh-Nmesh+j);
                 A(k+j,i*Nmesh-Nmesh+j)=-2*C11;%0;%%
                 A(k+j,i*Nmesh-Nmesh+m3)=0;
                 A(k+j,i*Nmesh-Nmesh+m4)=0;
                       % u displacement
                  A(k+j,w+m1*Nmesh-Nmesh+m3)=0;
                  A(k+j,w+m1*Nmesh-Nmesh+m4)=0;
                  A(k+j,w+m1*Nmesh-Nmesh+m3)=0+A(k+j,w+m1*Nmesh-Nmesh+m3);
                  A(k+j,w+m1*Nmesh-Nmesh+m4)=0+A(k+j,w+m1*Nmesh-Nmesh+m4);
                  A(k+j,w+i*Nmesh-Nmesh+m3)=0;
                  A(k+j,w+i*Nmesh-Nmesh+m4)=0;
                  A(k+j,w+m1*Nmesh-Nmesh+j)=50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
                  A(k+j,w+m1*Nmesh-Nmesh+j)=-50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0+A(k+j,w+m1*Nmesh-Nmesh+j);
	
                  b(k+j)=0.1*(C11+150*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0)+e0*(C11+C12)*(C(i+3,j+2,1)-C(i+1,j+2,1))/2.0+100*P*(C(i+3,j+2,1)-C(i+1,j+2,1))*e0*(C(i+2,j+2,1)-average_composition); %// Be cautious here...
                 % /*--------first equation ends---------*/
	
                     % /*-------second equation begins-------*/ 
                     % v displacement
                  A(w+k+j,m1*Nmesh-Nmesh+m3)=0;
                  A(w+k+j,m1*Nmesh-Nmesh+m4)=0;
                  A(w+k+j,m1*Nmesh-Nmesh+m3)=0+A(w+k+j,m1*Nmesh-Nmesh+m3);
                  A(w+k+j,m1*Nmesh-Nmesh+m4)=0+A(w+k+j,m1*Nmesh-Nmesh+m4);
                  A(w+k+j,m1*Nmesh-Nmesh+j)=0;
                  A(w+k+j,m1*Nmesh-Nmesh+j)=0+A(w+k+j,m1*Nmesh-Nmesh+j);
                  A(w+k+j,i*Nmesh-Nmesh+m3)=0;
                  A(w+k+j,i*Nmesh-Nmesh+m4)=0;
                        % u displacement
                  A(w+k+j,w+i*Nmesh-Nmesh+m3)=0;
                  A(w+k+j,w+i*Nmesh-Nmesh+m4)=0;
                  A(w+k+j,w+i*Nmesh-Nmesh+j)=-2*C44;%0;%%
                  A(w+k+j,w+m1*Nmesh-Nmesh+j)=C44;
                  A(w+k+j,w+m1*Nmesh-Nmesh+j)=C44+A(w+k+j,w+m1*Nmesh-Nmesh+j);
	
                  b(w+k+j)=2*epsilon_f*(50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0)-2*epsilon_f*((C12+C44)/4.0)+2*epsilon_f*((C12+C44)/4.0)+e0*(C12+C11)*(C(i+2,j+3,1)-C(i+2,j+1,1))/2.0+100*P*(C(i+2,j+3,1)-C(i+2,j+1,1))*e0*(C(i+2,j+2,1)-average_composition); %// Be cautious here...
                  % /*--------second equation ends--------*/
                  
                  continue;
              end
              
%               if(i==Nmesh-1)
%                   % /*-------first equation begins-------*/
%                   A(k+j,m1*Nmesh-Nmesh+j)=C11+150*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
%                   A(k+j,m2*Nmesh-Nmesh+j)=C11-150*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
%                   A(k+j,i*Nmesh-Nmesh+j)=-2*C11-2*C44;
%                   A(k+j,i*Nmesh-Nmesh+m3)=C44+50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
%                   A(k+j,i*Nmesh-Nmesh+m4)=C44-50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
%                   
%                   A(k+j,w+m1*Nmesh-Nmesh+m3)=0;
%                   A(k+j,w+m2*Nmesh-Nmesh+m4)=(C12+C44)/4.0;
%                   A(k+j,w+m2*Nmesh-Nmesh+m3)=-(C12+C44)/4.0;
%                   A(k+j,w+m1*Nmesh-Nmesh+m4)=0;
%                   A(k+j,w+i*Nmesh-Nmesh+m3)=50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
%                   A(k+j,w+i*Nmesh-Nmesh+m4)=-50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
%                   A(k+j,w+m1*Nmesh-Nmesh+j)=50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
%                   A(k+j,w+m2*Nmesh-Nmesh+j)=-50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
% 	
%                   b(k+j)=e0*(C11+C12)*(C(i+3,j+2,1)-C(i+1,j+2,1))/2.0+100*P*(C(i+3,j+2,1)-C(i+1,j+2,1))*e0*(C(i+2,j+2,1)-average_composition); %// Be cautious here...
% 	              % /*--------first equation ends---------*/
% 	
% 	              % /*-------second equation begins-------*/
%                   A(w+k+j,m1*Nmesh-Nmesh+m3)=0;
%                   A(w+k+j,m2*Nmesh-Nmesh+m4)=(C12+C44)/4.0;
%                   A(w+k+j,m2*Nmesh-Nmesh+m3)=-(C12+C44)/4.0;
%                   A(w+k+j,m1*Nmesh-Nmesh+m4)=0;
%                   A(w+k+j,m1*Nmesh-Nmesh+j)=50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
%                   A(w+k+j,m2*Nmesh-Nmesh+j)=-50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
%                   A(w+k+j,i*Nmesh-Nmesh+m3)=50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
%                   A(w+k+j,i*Nmesh-Nmesh+m4)=-50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
%                   
%                   A(w+k+j,w+i*Nmesh-Nmesh+m3)=C11+150*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
%                   A(w+k+j,w+i*Nmesh-Nmesh+m4)=C11-150*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
%                   A(w+k+j,w+i*Nmesh-Nmesh+j)=-2*C11-2*C44;
%                   A(w+k+j,w+m1*Nmesh-Nmesh+j)=C44+50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
%                   A(w+k+j,w+m2*Nmesh-Nmesh+j)=C44-50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
% 	
%                   b(w+k+j)=e0*(C12+C11)*(C(i+2,j+3,1)-C(i+2,j+1,1))/2.0+100*P*(C(i+2,j+3,1)-C(i+2,j+1,1))*e0*(C(i+2,j+2,1)-average_composition); %// Be cautious here...
% 	              % /*--------second equation ends--------*/
%               end
%               if(i==Nmesh-1)
%                   continue;
%               end
              
	   % /*-------first equation begins-------*/
              A(k+j,m1*Nmesh-Nmesh+j)=C11+150*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
	          A(k+j,m2*Nmesh-Nmesh+j)=C11-150*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
              A(k+j,i*Nmesh-Nmesh+j)=-2*C11-2*C44;
              A(k+j,i*Nmesh-Nmesh+m3)=C44+50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
              A(k+j,i*Nmesh-Nmesh+m4)=C44-50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
	
              A(k+j,w+m1*Nmesh-Nmesh+m3)=(C12+C44)/4.0;
              A(k+j,w+m2*Nmesh-Nmesh+m4)=(C12+C44)/4.0;
              A(k+j,w+m2*Nmesh-Nmesh+m3)=-(C12+C44)/4.0;
              A(k+j,w+m1*Nmesh-Nmesh+m4)=-(C12+C44)/4.0;
              A(k+j,w+i*Nmesh-Nmesh+m3)=50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
              A(k+j,w+i*Nmesh-Nmesh+m4)=-50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
              A(k+j,w+m1*Nmesh-Nmesh+j)=50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
              A(k+j,w+m2*Nmesh-Nmesh+j)=-50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
	
	          b(k+j)=e0*(C11+C12)*(C(i+3,j+2,1)-C(i+1,j+2,1))/2.0+100*P*(C(i+3,j+2,1)-C(i+1,j+2,1))*e0*(C(i+2,j+2,1)-average_composition); %// Be cautious here...
	   % /*--------first equation ends---------*/
	
	   % /*-------second equation begins-------*/
	          A(w+k+j,m1*Nmesh-Nmesh+m3)=(C12+C44)/4.0;
              A(w+k+j,m2*Nmesh-Nmesh+m4)=(C12+C44)/4.0;
              A(w+k+j,m2*Nmesh-Nmesh+m3)=-(C12+C44)/4.0;
              A(w+k+j,m1*Nmesh-Nmesh+m4)=-(C12+C44)/4.0;
              A(w+k+j,m1*Nmesh-Nmesh+j)=50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
              A(w+k+j,m2*Nmesh-Nmesh+j)=-50*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
              A(w+k+j,i*Nmesh-Nmesh+m3)=50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
              A(w+k+j,i*Nmesh-Nmesh+m4)=-50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
	
              A(w+k+j,w+i*Nmesh-Nmesh+m3)=C11+150*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
              A(w+k+j,w+i*Nmesh-Nmesh+m4)=C11-150*P*(C(i+2,j+3,1)-C(i+2,j+1,1))/4.0;
              A(w+k+j,w+i*Nmesh-Nmesh+j)=-2*C11-2*C44;
              A(w+k+j,w+m1*Nmesh-Nmesh+j)=C44+50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
              A(w+k+j,w+m2*Nmesh-Nmesh+j)=C44-50*P*(C(i+3,j+2,1)-C(i+1,j+2,1))/4.0;
	
              b(w+k+j)=e0*(C12+C11)*(C(i+2,j+3,1)-C(i+2,j+1,1))/2.0+100*P*(C(i+2,j+3,1)-C(i+2,j+1,1))*e0*(C(i+2,j+2,1)-average_composition); %// Be cautious here...
	   % /*--------second equation ends--------*/
       
              if(i==2)  % &&(j==1)
                  A(k+j,:)=zeros(1,Dimension);
                  A(k+j,k+j)=1;
                  b(k+j)=epsilon_f;
              end
              if(j==1)||(j==Nmesh)
                  A(k+j+w,:)=zeros(1,Dimension);
                  A(k+j+w,k+j+w)=1;
              end
          end
      end

      %Anor=sum(sum((A)))
      b=reshape(b,Dimension,1);
      
%      A(Nmesh+1:Dimension,1:Nmesh)=zeros(Dimension-Nmesh,Nmesh);
  
%         A(:,1)=zeros(Dimension,1);
%           A(1,1)=1;
%         A(:,1+w)=zeros(Dimension,1);
%           A(1+w,1+w)=1;
%         A(:,Nmesh)=zeros(Dimension,1);
%           A(Nmesh,Nmesh)=1;
%         A(:,Nmesh+w)=zeros(Dimension,1);
%           A(Nmesh+w,Nmesh+w)=1;
%         A(:,(Nmesh-1)*Nmesh+1)=zeros(Dimension,1);
%           A((Nmesh-1)*Nmesh+1,(Nmesh-1)*Nmesh+1)=1;
%         A(:,(Nmesh-1)*Nmesh+1+w)=zeros(Dimension,1);
%           A((Nmesh-1)*Nmesh+1+w,(Nmesh-1)*Nmesh+1+w)=1;
%         A(:,Nmesh*Nmesh)=zeros(Dimension,1);
%           A(Nmesh*Nmesh,Nmesh*Nmesh)=1;
%         A(:,Dimension)=zeros(Dimension,1);
%           A(Dimension,Dimension)=1;

      setup.type='ilutp';
      setup.milu='row';
      setup.droptol=1.0e-2;
      
      [L,U] = ilu(sparse(A),setup);
      %x=gmres(sparse(A),b,[],1.0e-8,50,L,U,x0);%,w,[],800);
%       if NT==0
        
%         U=cholinc(A,0.01);
%         L=U';
         % conditionNo=cond(A)
          %x=bicgstab(sparse(A),b,1.0e-8,50,[],[],x1-3*x2+3*x3);
          %x=cgs(sparse(A),b,1.0e-8,50,L,U,x1-3*x2+3*x3);
          x=gmres(sparse(A),b,10,1.0e-8,100,L,U,x1-3*x2+3*x3);%,w,[],800);
          %x=A\b;
          %x1=x2;
          %x2=x;
          x1=x2;
          x2=x3;
          x3=x;

      u1=reshape(x(1:w),Nmesh,Nmesh);
      u2=reshape(x(w+1:Dimension),Nmesh,Nmesh);
      u1=u1';
      u2=u2';
   
      %/*------calculate elastic energy miu_el[Nmesh+4][Nmesh+4] and e_el[Nmesh+4][Nmesh+4]----*/
      [u1_y,u1_x]=gradient(u1,1); % Du1/Dj, Du1/Di
      [u2_y,u2_x]=gradient(u2,1); % Du2/Dj, Du2/Di
      
      for i=1:Nmesh
          for j=1:Nmesh
              
              e1=u1_x(i,j)-e0*(C(i+2,j+2,1)-average_composition);
              e2=u2_y(i,j)-e0*(C(i+2,j+2,1)-average_composition);
              e6=u1_y(i,j)+u2_x(i,j);
              
              strain=[e1;e2;e6];
              
              % Calculate elastic constant
              interp=C(i+2,j+2,1)*C(i+2,j+2,1)*C(i+2,j+2,1)*(6*C(i+2,j+2,1)*C(i+2,j+2,1)-15*C(i+2,j+2,1)+10);
              
%               C11=225+150*interp;
%               C12=75+50*interp;
%               C44=75+50*interp;

%               C11=324-92*interp;
%               C12=60.8+92.2*interp;
%               C44=390-273*interp;

              C11=1+102.0*interp;
              C12=0.4+43.5*interp;
              C44=0.3+29.0*interp;
              
              stiffness=[C11 C12 0; C12 C11 0; 0 0 C44];
              
              % P denotes the derivative of interp w.r.t composition
              P=30*C(i+2,j+2,1)*C(i+2,j+2,1)*(C(i+2,j+2,1)-1)*(C(i+2,j+2,1)-1);
              
%               stiffnessderivative=50*P*[3 1 0; 1 3 0; 0 0 1];
              stiffnessderivative=P*[102 43.5 0;43.5 102 0;0 0 29];
              
              e_el(i+2,j+2)=0.5*strain'*stiffness*strain;
              miu_el(i+2,j+2)=-e0*(C11+C12)*(e1+e2)+0.5*strain'*stiffnessderivative*strain;
          end
      end

      for i=1:Nmesh+4
          ii=i;
%           if(i<3)
%             ii=i+Nmesh-1;
%           end
%           if(i>Nmesh+2)
%             ii=i-Nmesh+1;
%           end
          if(i<3)
             ii=6-i; % i+Nmesh-1; % Be careful that the period is (Nmesh-1)...
          end
          if(i>Nmesh+2)
             ii=2*Nmesh+4-i; % i-Nmesh+1;
          end
          for j=1:Nmesh+4
             if((i>=3)&&(i<=Nmesh+2)&&(j>=3)&&(j<=Nmesh+2)) 
               continue;
             end
             jj=j;
             if(j<3)
               jj=j+Nmesh-1;
             end
             if(j>Nmesh+2)
               jj=j-Nmesh+1;
             end
             miu_el(i,j)=miu_el(ii,jj);
          end
      end
              
  end  % /*-------------------Elastic Equilibrium Equations ENDS----------------------*/ 
%%    
  h=h+miu_el;
  
     % *calculate the whole energy of the system, of which "freeE" term is actually the bulk energy density, and e/2==1 due to e==2 here.Also,dx=dy=1*/
  if mod(NT,100)==0
     F=0.0;
     Out_freeE=0.0;
     Out_gradE=0.0;
     Out_e_el=0.0;
     for i=3:Nmesh+2
         for j=3:Nmesh+2
             %freeE=2*50*(C(i,j,1)-0.3)*(C(i,j,1)-0.3)*(C(i,j,1)-0.6)*(C(i,j,1)-0.6);
             %freeE=-(C(i,j,1)-0.5).*(C(i,j,1)-0.5)+2.5.*(C(i,j,1)-0.5).*(C(i,j,1)-0.5).*(C(i,j,1)-0.5).*(C(i,j,1)-0.5);
             freeE=2*C(i,j,1)*C(i,j,1)*(1-C(i,j,1))*(1-C(i,j,1));
             % //(WRONG!!)-- gradE=C[i-1][j][0]+C[i+1][j][0]+C[i][j-1][0]+C[i][j+1][0]-4*C[i][j][0];
             gradE=0.5*((C(i+1,j,1)-C(i,j,1))*(C(i+1,j,1)-C(i,j,1))+(C(i,j+1,1)-C(i,j,1))*(C(i,j+1,1)-C(i,j,1)));
             F=F+freeE+gradE+e_el(i,j);
             Out_freeE=Out_freeE+freeE;
             Out_gradE=Out_gradE+gradE;
             Out_e_el=Out_e_el+e_el(i,j);
         end
     end
    % OUT_Tot_freeE_gradE_elaE=[Out_freeE+Out_gradE Out_freeE Out_gradE Out_e_el]
     Energy_Counter=Energy_Counter+1;
     Time_Energy(Energy_Counter,1)=NT;  % store the evotion time
     Time_Energy(Energy_Counter,2)=F;   % store the system energy
     OUT=[NT F Out_freeE Out_gradE Out_e_el sum(sum(C(3:Nmesh+2,3:Nmesh+2,1)))/(Nmesh*Nmesh)]
     OUTdata(counter,:)=OUT;
     counter=counter+1;
  end
     %% /*-------------------EVOLUTION BEGINS----------------------*/
  u(2:Nmesh+3,2:Nmesh+3)=h(3:Nmesh+4,2:Nmesh+3)+h(1:Nmesh+2,2:Nmesh+3)+h(2:Nmesh+3,3:Nmesh+4)+h(2:Nmesh+3,1:Nmesh+2)-4.0*h(2:Nmesh+3,2:Nmesh+3);
  v(2:Nmesh+3,2:Nmesh+3)=C(3:Nmesh+4,2:Nmesh+3,1)+C(1:Nmesh+2,2:Nmesh+3,1)+C(2:Nmesh+3,3:Nmesh+4,1)+C(2:Nmesh+3,1:Nmesh+2,1)-4.0*C(2:Nmesh+3,2:Nmesh+3,1);
     
  uu(3:Nmesh+2,3:Nmesh+2)=v(4:Nmesh+3,3:Nmesh+2)+v(2:Nmesh+1,3:Nmesh+2)+v(3:Nmesh+2,4:Nmesh+3)+v(3:Nmesh+2,2:Nmesh+1)-4.0*v(3:Nmesh+2,3:Nmesh+2); % *DX==1;
  
     % Main loop
  C(3:Nmesh+2,3:Nmesh+2,2)=C(3:Nmesh+2,3:Nmesh+2,1)+dt*(u(3:Nmesh+2,3:Nmesh+2)-1*uu(3:Nmesh+2,3:Nmesh+2)); % The last red number is gradient energy coef, here it's 1.
  C(3:Nmesh+2,3:Nmesh+2,1)=C(3:Nmesh+2,3:Nmesh+2,2);                                                       % Note that we use k/2 and e/2 in the energy expression.
  
     %%  /*-------------------EVOLUTION ENDS----------------------*/
%%
  if NT==0
     C_0=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_0.mat C_0;
  end 
  if NT==100
     C_100=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_100.mat C_100;
  end 
  if NT==1000
     C_1000=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_1000.mat C_1000;
  end 
  if NT==5000
     C_5000=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_5000.mat C_5000;
  end 
  if NT==10000
     C_10000=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_10000.mat C_10000;
  end 
  if NT==50000
     C_50000=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_50000.mat C_50000;
  end 
  if NT==100000
     C_100000=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_100000.mat C_100000;
  end 
  if NT==200000
     C_200000=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_200000.mat C_200000;
  end 
  if NT==500000
     C_500000=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_500000.mat C_500000;
  end 
  if NT==1000000
     C_1000000=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_1000000.mat C_1000000;
  end 
  if NT==2000000
     C_2000000=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_2000000.mat C_2000000;
  end 
  if NT==3999999
     C_3999999=C(3:Nmesh+2,3:Nmesh+2,2);
     save C_3999999.mat C_3999999;
  end
 
end
save Time_Energy.mat Time_Energy
TimeEnd=fix(clock);
save TimeEnd.mat TimeEnd
% end
%% /*-----------------------MAIN FUNCTION ENDS--------------------------*/

