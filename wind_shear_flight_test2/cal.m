function df = cal(t,x,u,alpha,h)
% h_ref=30.50;
df=zeros(6,1);
df(1)=0.04937-16.22*u*cos(alpha)*x(3)^2 - 0.207*x(1)*sqrt(x(1)^2+x(2)^2);
df(2)=16.22*u*x(3)^2*sin(alpha)-0.207*x(2)*sqrt(x(1)^2+x(2)^2);
df(3)=-349.22*(1.156*x(4) + (x(2)*sin(alpha)-x(1)*cos(alpha))/x(3))*u*x(3)^2 - 0.01822*x(3)^2*(1+4.65*((x(1)*sin(alpha)+x(2)*cos(alpha))/x(3))^2);
df(4)=33.55*x(3)*u-77.50*x(3)*x(4)*sqrt((1.156*x(4)+(x(2)*sin(alpha)-x(1)*cos(alpha))/x(3))^2+((x(2)*cos(alpha)+x(1)*sin(alpha))/x(3))^2 + 0.5*u*(1-2*((x(2)*cos(alpha)+x(1)*sin(alpha))/(x(3)*sqrt(0.5*u)))^2)*(1/(2 + (1.156*x(4)+(x(2)*sin(alpha)-x(1)*cos(alpha))/x(3))/sqrt(0.5*u))^2 - ((1.156*x(4)+(x(2)*sin(alpha)-x(1)*cos(alpha))/x(3))/sqrt(0.5*u))^2 + (1+(1.156*x(4)+(x(2)*sin(alpha)-x(1)*cos(alpha))/x(3))/sqrt(0.5*u))*(0.109+0.217*(((1.156*x(4)+(x(2)*sin(alpha)-x(1)*cos(alpha))/x(3))/sqrt(0.5*u))-0.15)^2)));
df(5)=-(37*5.37/h)*x(1);
df(6)=(37*5.37/h)*x(2);
end