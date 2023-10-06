x=[0,1,2,3,4,5];
y=[2.1,7.7,13.6,27.2,40.9,61.1];
p=polyfit(x,y,2);
x_fit=0:.2:5;
y_fit=polyval(p,x_fit);
plot(x,y,'ro',x_fit,y_fit)