% Bidirectional duffing SDOF system simulated using correlated white noise base excitation
%
clc;
clear;
% close all;
% % 
%% System parameters
m = 1;
omegax = 2*pi*0.2;
zetax = 0.05;
alphax = 1;
S0x = 0.005;

T = 100; % Total time of integration
dt = 0.01; % Time step
t = 0:dt:T-dt;    % Time vector
Nt = length(t);
Nsim = 1; % MC simulations

deltamat = [sqrt(dt)            0;
            dt^1.5/2    dt^1.5/(2*sqrt(3))];
%
sigx = 0.0001;%sqrt(2*pi*S0x);

% load DATA.mat
for j = 1:Nsim
    j
    %
    y = zeros(2,Nt);
    ff = 10*sin(2*pi*0.6*t);

    for i = 1:Nt-1
    FORCE = [0;ff(i)];
        
%     dZ1 = deltamat(2,:)*[DATA(1,i,j);DATA(1,i,j)]; % 1 x 1
%     dW1 = (deltamat(1,:)*[DATA(2,i,j);DATA(2,i,j)])'; % 1 x 1

    dZ1 = deltamat(2,:)*randn(2,1); % 1 x 1
    dW1 = (deltamat(1,:)*randn(2,1))'; % 1 x 1

    % b matrix
    b11 = 0;
    b21 = -sigx;

    b_coef = [b11;
              b21];    % 4 x 2
    
    % Drift Coefficients
    
    a1 = y(2,i);
    a2 = -2-zetax*omegax*y(2,i)-omegax^2*y(1,i)-alphax/m*(y(1,i))^3;
    
    a_coef = [a1; a2];
    
    % Higher order Moments
    L1a1=b21; 
    L1a2 = -2*zetax*omegax*b21;
    L2a1=0; 
    L2a2 = -2*zetax*omegax*0;
    
    L0a1=a2;
    L0a2 = -omegax^2*a1-2*zetax*omegax*a2-3*alphax/m*(y(1,i))^2*a1;
    
    DW1 = dW1;   % 1 x 1
    DZ1 = dZ1;   % 1 x 1

       L1a = [L1a1 L1a2]';
       L2a = [L2a1 L2a2]';
       L0a = [L0a1 L0a2]';
   
    y(:,i+1) = y(:,i) + a_coef.*dt + b_coef(:,1).*DW1 +...
               L1a.*DZ1 + L0a*0.5*dt^2 + FORCE*dt; 
    end

% Storage of variables,
Y1(j,:) = y(1,:);
Y2(j,:) = y(2,:);
end

%
% Sample solution of the Monte Carlo prediction :
u1= mean(Y1,1);

% figure,plot(t,u1,'b')
% figure(1),plot(t,u1,'b');hold on
%
% [f,s] = fftplot(dt,u1);
% figure(2),plot(f,s); xlim([0 5])
% hold on

u1d= mean(Y2,1);
figure(2),plot(u1,u1d);hold on

% End