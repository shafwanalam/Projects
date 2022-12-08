%Steady one dimensional heat conduction program
clear all

%Input data
L=1.5; %[m]
N=20; %number of grid points
N_list=4:5:200;
alpha=20.0; %thermal diffusivity [m^2/s]
T_left=350.0; %Temperature left wall [K]
T_right=300.0; %Temperature right wall [K]
Qt=50; %Heat source proportional to temperature, i.e.Qt*T, [1/s]
maxIter=1e6; %Maximum number of iterations for the solver
maxRes=1e-6; %Solver residual |Ax-b|<maxRes

Error=zeros(1,length(N_list));
cas=1;
for N=N_list
%Mesh
%posX=0:L/(N-1):L;
x=linspace(0,L,N);

%Preallocate coefficients & mesh distances
ap=ones(1,N);
ae=zeros(1,N);
aw=zeros(1,N);
b=zeros(1,N);
Dxe=zeros(1,N);
Dxw=zeros(1,N);
Dx=zeros(1,N);

%Preallocate vector for local residual of solver
loc_res=zeros(1,N);

%Inner coefficients
for iX=2:N-1
  Dxe(iX)=(x(iX+1)-x(iX));
  Dxw(iX)=(x(iX)-x(iX-1));
  Dx(iX)=(Dxe(iX)+Dxw(iX))/2;
  ae(iX)=alpha/Dxe(iX);
  aw(iX)=alpha/Dxw(iX);
  ap(iX)=-(ae(iX)+aw(iX)+Qt*Dx(iX));
  b(iX)=0;
end

%Boundary coefficients
b(1)=T_left;
b(end)=T_right;

%Initialize temperature vector
%T=ones(size(x))*((Tw+Te)*0.5);
T=zeros(size(x));
T(1)=T_left;
T(end)=T_right;

%Solver
res=maxRes+1;
ite=0;
tic;
Jacobi = [];
while res>maxRes && ite<maxIter
  %Gauss-Seidel iteration
  T(1)=(-ae(1)*T(2)+b(1))/ap(1);
  for iX=2:numel(T)-1
    T(iX)=(-aw(iX)*T(iX-1)-ae(iX)*T(iX+1)+b(iX))/ap(iX);
  end
  T(end)=(-aw(end)*T(end-1)+b(end))/ap(end);
  
  %Calculation of the solver residual res=|Ax-b|
  for iX=2:numel(T)-1
    loc_res(iX)=aw(iX)*T(iX-1)+ap(iX)*T(iX)+ae(iX)*T(iX+1)-b(iX);
  end
  res=max(abs(loc_res));

  ite=ite+1;
  if mod(ite,10000)==0
    fprintf('ite: %d solver residual: %e\n',ite,res);
  end
end
if ite==maxIter
  warning(['Maximum number of iterations reached (lastRes=',num2str(res),').']);
end

%Matlab solver for systems of linear equations
%A is the full matrix. Try full(A) to see the whole matrix
A=spdiags([[aw(2:end)'; 0] ap(:) [0; ae(1:end-1)']],[-1 0 1],N,N);
T2=A\b(:);
toc;

%Tana=@(x) (-0.5*Q/lambda)*x.^2+((T_right-T_left+0.5*Q/lambda*L^2)/L)*x+T_left; %analytic solution
k=sqrt(Qt/alpha);
C2 = (T_right-T_left*exp(k*L))/(exp(-k*L)-exp(k*L));
C1 = T_left-C2;
Tana=@(x) C1*exp(k*x)+C2*exp(-k*x); %analytic solution

%plot(x,T);
%hold on
%plot(x,T2);
[AX,H1,H2] = plotyy(x,T,x,feval(Tana,x)-T);
legend('T numerical','Error')
ylabel(AX(1),'Temperature') % left y-axis 
ylabel(AX(2),'Error') % right y-axis
title('Temperature distribution');
drawnow

Error(cas)=max(abs(feval(Tana,x)-T));
cas=cas+1;
end

figure(2);
loglog(L./N_list,Error,'-o','LineWidth',3);
ylabel('Error');
xlabel('Dx')
title('Truncation error analysis');

