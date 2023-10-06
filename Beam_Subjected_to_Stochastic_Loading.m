clc;
clear all;
%Parameters
m=30000;
l=10;
EJ=13748*10^4;
c=0.91;
n=1;
omega_n=(n*pi/l)^2*sqrt(EJ/m);
alpha=c/(2*m);
v0=30;
sigma_v=0.025;

%Time Steps
T=100;
dt=0.01;
t=0:dt:T-dt;
Nt=length(t);
Nsim=1;

deltamat = [sqrt(dt)            0;
            dt^1.5/2    dt^1.5/(2*sqrt(3))];

for j=1:Nsim
    X=0.1*randn;
    forcing_func=10*sin(n*pi*X/l);
    y=zeros(4,Nt);
    for i=1:Nt-1
        dZ1 = deltamat(2,:)*randn(2,1);
        dW1 = (deltamat(1,:)*randn(2,1))';
        %Drift Coeficients drift is doubt
         z3=sin(n*pi*X/l);
         z4=cos(n*pi*X/l);
         a1 = y(2,i);
         a2=-omega_n^2*y(1,i)-2*alpha*y(2,i)+2/(m*l)*z3;
         a3=n*pi*v0/l*z4-0.5*(n*pi/l)^2*sigma_v^2*z3;
         a4=-n*pi*v0/l*z3-0.5*(n*pi/l)^2*sigma_v^2*z4;
         a_coeff=[a1;a2;a3;a4];
         b_coeff=[0;0;n*pi*sigma_v*z4/l;-n*pi*sigma_v*z3/l];

         %Kolmogorov Moments (time differntials consired 0)
         L0a1=a2;
         L0a2=-a1*omega_n^2-2*a2*alpha+2*a3/(m*l);
         L0a3=-a3*(n*pi*v0/l*tan(n*pi*X/l)+0.5*(n*pi/l)^2*sigma_v^2)+a4*(n*pi*v0/l+0.5*(n*pi/l)^2*sigma_v^2*z4)-1/2*(n*pi*sigma_v/l)^2*n*pi*a3*v0/(z4^2*l);
         L0a4=a3*(0.5*(n*pi/l)^2*sigma_v^2*tan(n*pi*X/l)-n*pi*v0/l)+a4*(n*pi*v0/l*cot(n*pi*X/l)-0.5*(n*pi/l)^2*sigma_v^2)-n*pi*sigma_v*0.5*(n*pi/l)^2*sigma_v^2/(z4*l);

         L1a1=0;
         L1a2=2*n*pi*sigma_v*z4/(m*l^2);
         L1a3=n*pi*sigma_v*z4/l*(-n*pi*sigma_v*v0*tan(n*pi*X/l)/l-0.5*(n*pi/l)^2*sigma_v^2)-n*pi*sigma_v*z3/l*(n*pi*sigma_v*v0+0.5*(n*pi/l)^2*sigma_v^2*cot(n*pi*X/l));
         L1a4=n*pi*sigma_v*z4/l*(-n*pi*sigma_v*v0+0.5*(n*pi/l)^2*sigma_v^2*tan(n*pi*X/l))+n*pi*sigma_v*z3/l*(n*pi*sigma_v*v0*cot(n*pi*X/l)-0.5*(n*pi/l)^2*sigma_v^2);

         DW1 = dW1;   % 1 x 1
         DZ1 = dZ1;

         L1a = [L1a1 L1a2 L1a3 L1a4]';
         L0a = [L0a1 L0a2 L0a3 L0a4]';

         y(:,i+1) = y(:,i) + a_coeff*dt + b_coeff(:,1).*DW1 +...
               L1a.*DZ1 + L0a*0.5*dt^2 + forcing_func*dt;
         end

% Storage of variables,
Y1(j,:) = y(1,:);
Y2(j,:) = y(2,:);
end
u1= mean(Y1,1);
u1d= mean(Y2,1);
figure(2),plot(u1,u1d);hold on
%plot(t,Y1);hold on;
%plot(t,Y2)

       

         









