clc; clear; close all;
%% 1. Load data
data = load('C:\Users\Admin\Downloads\Kobe_earthquake.txt');
time = data(:,1);
xgdd = data(:,2);
dt = time(2) -  time(1);

%% System Properties
m = 0.5;     % kN-sec^2/cm
k = 18;         % kN/cm
c = 0.3;     % kN-sec/cm

% Initial conditions

u0 = 0;         % cm
ud0 = 0;        % cm/s
udd0 = xgdd(1) - ud0*c/m - u0*k/m;

% case: constant average acceleration method
% gamma = 1/2 and beta = 1/4
gamma = 1/2;
beta = 1/4;

% calculation of a1, a2 and a3

a1 = m/(beta*dt^2) + gamma*c/(beta*dt);
a2 = m/(beta*dt) + (gamma/beta  - 1)*c;
a3 = (1/(2*beta) - 1)*m + dt*(gamma/(2*beta) - 1)*c;

khat = k + a1;
%% Calculation for each time step
u = zeros(size(time));
ud = zeros(size(time));
udd = zeros(size(time));
phat = zeros(size(time));

u(1) = u0;
ud(1) = ud0;
udd(1) = udd0;

for i = 1:length(time)-1
    phat(i+1) = m*xgdd(i+1) + a1*u(i) + a2*ud(i) + a3*udd(i);
    u(i+1) = phat(i+1)/khat;
    ud(i+1) = gamma/(beta*dt)*(u(i+1) - u(i)) + (1 - gamma/beta)*ud(i) + dt*(1 - gamma/(2*beta))*udd(i);
    udd(i+1) = 1/(beta*dt^2)*(u(i+1) - u(i)) - 1/(beta*dt)*ud(i) - (1/(2*beta) - 1) *udd(i);
end

%% Plots
figure,
set(gcf,'color','w')
plot(time,xgdd,'r','LineWidth',2);
xlabel('t [sec]');
ylabel('Ground Accln');
set(gca,'fontname','Latin Modern Roman','fontweight','bold','fontsize',18,'linewidth',1.2)

figure,
set(gcf,'position',[200 120 1000 800],'color','w')
subplot(3,1,1)
plot(time,u,'k','LineWidth',2);
xlabel('t [sec]');
ylabel('Displacement, u');
set(gca,'fontname','Latin Modern Roman','fontweight','bold','fontsize',18,'linewidth',1.2)

subplot(3,1,2)
plot(time,ud,'k','LineWidth',2);
xlabel('t [sec]');
ylabel('Velocity, ud');
set(gca,'fontname','Latin Modern Roman','fontweight','bold','fontsize',18,'linewidth',1.2)

subplot(3,1,3)
plot(time,udd,'k','LineWidth',2);
xlabel('t [sec]');
ylabel('Acceleration, udd');
set(gca,'fontname','Latin Modern Roman','fontweight','bold','fontsize',18,'linewidth',1.2)






