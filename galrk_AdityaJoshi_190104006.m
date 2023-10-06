syms u(x)
ode=(1+x)*diff(u,x,2)+diff(u,x)+x^2==0;
cond1=u(0)==0;
cond2=u(1)==0;
conds=[cond1 cond2];

uSol(x)=dsolve(ode,conds);

syms x
p1=sin(2*pi*x);
p2=sin(2*pi*x)^2*cos(x);
p3=sin(2*pi*x)^3*cos(x)^2;
p4=sin(2*pi*x)^2*cos(1-x);
k=1+x;
f=-x^2;

p=[p1 p2 p3 p4 ];

K=zeros(4,4);

for i=1:4
    for j=1:4
        K(i,j)=val(p(i),p(j),k);
    end
end

F=zeros(4,1);

for i=1:4
    F(i)=l2(p(i),f);
end

c=K\F;

syms P(x)
P=0;
for j=1:4
    P=P+c(j)*p(j);
end

fplot(uSol,'g',[0 1])
hold on
fplot(P,'r',[0 1])


syms d(x)
syms g(x)
syms a(x)

function I=val(d,g,a)
    d=diff(d);
    g=diff(g);
    I=int(d*g*a,[0 1]);
end

function Z=l2(d,g)
    Z=int(d*g,[0 1]);
end
       










