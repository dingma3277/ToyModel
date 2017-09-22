%grid 
yl=13.3;
ny=1000;
dy=2*yl/1000;
x=[-yl+0.5*dy:dy:yl-0.5*dy]*1500;
t=[0:0.8:160];
k=fopen('h.dat','r');
n=fopen('evap.dat','r');
o=fopen('prec.dat','r');
p=fopen('u.dat','r');
q=fopen('pt.dat','r');
r=fopen('atendh.dat','r');
s=fopen('v.dat','r');
h=zeros(200,1000);
prec=h;evap=h;u=h;
pt=h;atendh=h;v=h;
dt=zeros(10,1000);
dt2=zeros(size(dt));
for i=1:200
[a,count]=fread(k,1000,'float');
%[b,count]=fread(m,1000,'float');
[c,count]=fread(n,1000,'float');
[d,count]=fread(o,1000,'float');
[dd,count]=fread(p,1000,'float');
[bb,count]=fread(q,1000,'float');
[cc,count]=fread(r,1000,'float');
[kk,count]=fread(s,1000,'float');
h(i,:)=a';
%xm(i,:)=b';
evap(i,:)=c';
prec(i,:)=d';
u(i,:)=dd';
pt(i,:)=bb';
atendh(i,:)=cc';
v(i,:)=kk';
dt2(i,2:999)=h(i,1:998)-2.*h(i,2:999)+h(i,3:1000);
dt1(i,1)=h(i,2)-2.*h(i,1)+h(i,999);
dt1(i,1000)=h(i,1)-2.*h(i,1000)+h(i,999);
end

gf=zeros(1000,1000);
k=fopen('gfunc.dat','r')

for i=1:1000
[a,count]=fread(k,1000,'float');
gf(i,:)=a';
end

figure(1)
[cc,hh]=contour(x,t(2:200),u(2:200,:)*17);clabel(cc,hh);
xlabel('x (km)')
ylabel('t (days)')
title('u (m/s)')
%print -depsc u_contour.eps

figure(2);
hold off;
plot(x,prec(160,:));
hold
plot(x,u(160,:)*17,'r-.');
legend('precip','u');
xlabel('x (km)')
%print -depsc u_and_p.eps

figure(3)
%[cc,hh]=contour(x,t(2:200),pt(2:200,:)*17);clabel(cc,hh);
plot(x,h(160,:)/70.);
xlabel('x (km)')
ylabel('t (days)')
title('h')
%print -depsc u_contour.eps

%figure(3);
%[cc,hh]=contour(x,t(2:200),pt(2:200,:));clabel(cc,hh);
