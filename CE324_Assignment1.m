g=9.81;
m=70;
Cd=0.2;
v(1)=0;
t=0:0.5:20;
dt=0.5;

for i=1:length(t)-1
    v(i+1)=v(i)+(g-((Cd/m)*(v(i)^2)))*dt;
end
plot(t,v,'r')
xlabel('Time(s)')
ylabel('Velocity(m/s)')
title('Velocity of bungee jumper')

Va=sqrt(g*m/Cd)*tanh(sqrt(g*Cd/m)*t);
hold on
plot(t,Va,'g')
legend('Numerical Solution','Analytical Solution')


