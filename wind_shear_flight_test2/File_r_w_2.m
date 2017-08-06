clear all;close all;clc;
addpath('apm');
s='http://byu.apmonitor.com';

                          %%  Initial condition with time delay %%
vel=5.5;          % forward velocity
h=100;            % altitude
h1=h;
m=1360;                 % mass of helicopter
rho=1.2256;            % density of air

fe=2.8;                % drag 
g=9.81;                 % gravitational 
A=90.61;                % area of disk
R=5.37;                 % Radius of rotor
omega_0=37;             % rpm of rotor
for ii=1:length(h)
    for jj=1:length(vel)
E_init(ii,jj)=m*g*h(ii) + 0.5*m*(vel(jj).^2)+ 0.5*973.58*(omega_0)^2;
    end
end
format long
for j=1:length(h)
        H(j)=h(j);
for i=1:length(vel)
   alpha(i)=atan(rho*fe*vel(i).^2/(2*m*g)) ;
   Ct(i)=(m*g)/(rho*A*R^2*omega_0^2*cos(alpha(i)));
   nu(i)=vel(i)/(omega_0*R);
   ff=@induced;
    L(i)=fsolve(@(lambda)ff(lambda,nu(i),alpha(i),Ct(i)),0.0036);
    a_i=(L- nu.*sin(alpha))/1.154;
        tspan=linspace(0,0.75,100);
        x0(j,i,:)=[0 nu(i) 1 a_i(i) 1 0];
       [t,x(j,i,:,:)]=ode45(@(t,x)cal(t,x,Ct(i),alpha(i),h(j)),tspan,x0(j,i,:)); 
       aa=x(j,i,:,1);bb=x(j,i,:,2);cc=x(j,i,:,3);dd=x(j,i,:,4);ee=x(j,i,:,5);ff=x(j,i,:,6);
       subplot(3,2,1);plot(tspan,aa(:,:));subplot(3,2,2);plot(tspan,bb(:,:));subplot(3,2,3);plot(tspan,cc(:,:));
       subplot(3,2,4);plot(tspan,dd(:,:));subplot(3,2,5);plot(tspan,ee(:,:));subplot(3,2,6);plot(tspan,ff(:,:))
% %        figure
% %        plot(tspan,Ct(i)')
% %        figure
AA=x(:,:,end,1);BB=x(:,:,end,2);CC=x(:,:,end,3);DD=x(:,:,end,4);EE=x(:,:,end,5);FF=x(:,:,end,6);
x_int=[AA(:,:,end) BB(:,:,end) CC(:,:,end) DD(:,:,end) EE(:,:,end) FF(:,:,end) ]
aa=AA(j,i,end);
end
end

                           %% reading text of existing file  and writting a new apm file%%
fid =fopen('wind_shear_flight_test2.apm');
C=textscan(fid,'%s','delimiter','\n');
[n,m]=size(AA);
for k=1:n
for i=1:m


C{1,1}{3,1}=['u =' num2str(Ct(i)) ' , >0.0001 , <=1                          % set Upper and Lower limit '];
C{1,1}{4,1}=['alpha =' num2str(alpha(i)) '  , >=-0.2094 , <= 0.08726              % -12 degree lower limit and 5 degree upper limit'];    % num2str(alpha(i))
C{1,1}{10,1}=['href=' num2str(h(j))];
C{1,1}{14,1}=['x1 =' num2str(AA(j,i)) '  ,>=0 ,<=0.04826 '];
C{1,1}{15,1}=['x2 =' num2str(BB(j,i)) ',>0'];
C{1,1}{16,1}=['x3 =' num2str(CC(j,i)) ', >=0.65 ,<=1.15 '];
C{1,1}{17,1}=['x4 =' num2str(DD(j,i))];
C{1,1}{18,1}=['x5 =' num2str(EE(j,i)) ' ,>=0 ,< 1 '];
C{1,1}{19,1}=['x6 =' num2str(FF(j,i)) '  ,>=0 '];

                            %% print new file %%
fName = 'new_file.apm';
fid = fopen(fName,'w');            % Open the file
for k=1:numel(C{1,1}) 
  fprintf(fid,'%s\r\n',C{1,1}{k,1}); 
end
a='File_r_w_2';
                       %% updating the .txt file to apm %%
apm(s,a,'clear all');

apm_load(s,a,'new_file.apm');
csv_load(s,a,'data111.csv');


apm_option(s,a,'nlc.nodes',6)           % No of collocation points
apm_option(s,a,'nlc.solver',3)          % 1= APOPT & 3= IPOPT
apm_option(s,a,'nlc.imode',6)            % 6= Dynamic Optimization

apm_option(s,a,'nlc.max_iter',1000)
apm_option(s,a,'nlc.mv_type',1)
apm_option(s,a,'nlc.cv_type',1)
apm_option(s,a,'nlc.atol',1e-8)
apm_option(s,a,'nlc.rtol',1e-8)
% apm_option(s,a,'nlc.ctrl_hor',0.25)
% apm_option(s,a,'nlc.ctrl_time',0.01)
% apm_option(s,a,'nlc.diaglevel',0)
% apm_option(s,a,'nlc.otol',1e-4)
% apm_option(s,a,'nlc.pred_hor',0.25)
% apm_option(s,a,'nlc.reqctrlmode',3)
% apm_option(s,a,'nlc.scaling',0)


apm_info(s,a,'FV','tf')
apm_option(s,a,'tf.status',1)

apm_info(s,a,'SV','x1')
% apm_option(s,a,'x1.status',1)
% apm_option(s,a,'x1.dcost',1e-1)
% apm_option(s,a,'x1.cost',1e-1)
% apm_option(s,a,'x1.dmax',1e-1)
apm_option(s,a,'x1.lower',0)
apm_option(s,a,'x1.upper',0.04602)
apm_info(s,a,'SV','x2')
apm_option(s,a,'x2.lower',0)
% apm_option(s,a,'x2.status',1)
% apm_option(s,a,'x2_type',1)
% apm_option(s,a,'x2.dcost',1e-1)
% apm_option(s,a,'x2.cost',1e-1)
% apm_option(s,a,'x2.dmax',1e-1)
% apm_option(s,a,'x2.tr_init',1)
apm_info(s,a,'SV','x6')
% apm_option(s,a,'x6.status',1)
% apm_option(s,a,'x6_type',1)
% apm_option(s,a,'x6.dcost',1e-1)
% apm_option(s,a,'x6.cost',1e-1)
% apm_option(s,a,'x6.dmax',1e-1)
% apm_option(s,a,'x6.lower',0)
apm_info(s,a,'SV','x5')

% apm_option(s,a,'x5.sp',0)
apm_info(s,a,'SV','x4')
apm_info(s,a,'SV','x3')
% apm_info(s,a,'x3.lower',0.65)
% apm_info(s,a,'x3.upper',1.15)
% apm_info(s,a,'CV','x7')


apm_info(s,a,'MV','u')
apm_option(s,a,'u.status',1)
apm_option(s,a,'u_type',1)
apm_option(s,a,'u.dmax',1e-5)
 apm_option(s,a,'u.cost',1e-1)
 apm_option(s,a,'u.dcost',1e-6)
 apm_info(s,a,'u.lower',0.00002)
apm_info(s,a,'u.upper',0.0072)

apm_info(s,a,'MV','alpha')
apm_option(s,a,'alpha.status',1)
apm_option(s,a,'alpha_type',1)
apm_option(s,a,'alpha.dcost',1e-6)
apm_option(s,a,'alpha.cost',1e-3)
apm_option(s,a,'alpha.dmax',1e-3)
% apm_info(s,a,'alpha.lower',-0.2094)
% apm_info(s,a,'alpha.upper',0.1047)

%apm_info(s,a,'CV','obstacle');
%apm_option(s,a,'obstacle.sphi',1);
%apm_option(s,a,'obstacle.splo',0);
%apm_option(s,a,'obstacle.wsphi',-100);
%apm_option(s,a,'obstacle.wsplo',0);
%apm_option(s,a,'obstacle.tr_init',0);
%apm_option(s,a,'obstacle.status',1);

output=apm(s,a,'solve')
disp(output)
y=apm_sol(s,a)
apm_web_var(s,a)
apm_get(s,a,'infeasibilities.txt')
z=y.x;
tf=z.tf(end)
t=z.time;
t_flight=t*tf;
t_flight1=t_flight;
figure
subplot(3,2,1)
plot(t,z.x1,'--','Linewidth',2)
xlabel('time');ylabel('sink rate');grid on;
subplot(3,2,2)
plot(t,z.x2,'--','Linewidth',2)
xlabel('time');ylabel('forward veocity');grid on ;
subplot(3,2,3)
plot(t,z.x3,'--','Linewidth',2)
xlabel('time');ylabel('\omega/\omega_0');grid on;
subplot(3,2,4)
plot(t,z.x5,'-.','Linewidth',2);%set(gca,'Xdir','reverse')
xlabel('time')
ylabel('height')
grid on 
subplot(3,2,5)
plot(t,z.alpha*57.29);xlabel('Time(sec)');ylabel('Alpha')
subplot(3,2,6)
plot(t,z.u,'Linewidth',2);xlabel('time(sec)');ylabel('Ct'); %set(gca,'Xdir','reverse')
figure
subplot(2,2,1)
plot(z.x5,z.x1);xlabel('height');ylabel('SInk rate');set(gca,'Xdir','reverse')
subplot(2,2,2)
plot(z.x5,z.x2);xlabel('height');ylabel('forward velocity');set(gca,'Xdir','reverse')
subplot(2,2,3)
plot(z.x5,z.x3);xlabel('height');ylabel('\omega/\omega_0');set(gca,'Xdir','reverse')
subplot(2,2,4)
plot(z.x5,z.u);xlabel('height');ylabel('C_T');set(gca,'Xdir','reverse')
figure
CT=z.u;
Ct_z=CT.*cos(z.alpha);
Ct_x=CT.*sin(z.alpha);
plot(t,Ct_x,'-.',t,Ct_z,'-.k')
omega_0=37;
R=5.37;
nd2d=(omega_0*R);r2d=180/pi;m=1360;
w=z.x1*nd2d;u=z.x2*nd2d;omega=z.x3*omega_0;V=sqrt(u.^2+w.^2);rho=1.2256;fe=2.32;sigma=0.048;CD_0=0.0087;R=5.37;A=90.61;a=z.x4;I=1822;
mu=(u.*cos(z.alpha)+w.*sin(z.alpha))./(omega*R);
P_parasite=0.5*rho*fe*(V).^3;
C_profile=sigma*CD_0*(1+4.65*mu.^2)/8;
P_profile=rho*A*R^3*omega.^3.*C_profile;
nu_i=0.866*a;
P_induced=rho*A*(omega.*R).^3.*Ct.*nu_i;
P_total=P_parasite+P_induced+P_profile;
figure
subplot(1,2,1)
plot(t,P_total,t,P_parasite,'--',t,P_profile,'-.',t,P_induced,':','Linewidth',2);xlabel('time(sec)')
subplot(1,2,2)
plot(mu,P_total,mu,P_parasite,'--',mu,P_profile,'-.',mu,P_induced,':','Linewidth',2);xlabel('\mu')
legend('Total Power','Parasite Power','Profile Power','Induced Power')
figure
KE=0.5*m*V.^2;                                                          % Kinetic Energy
PE=m*9.81*(z.x5*100);                                                 % Potential Energy
RE=0.5*I*(omega.^2);                                                   % Ritational Energy 
TE=KE+PE+RE;
plot(t,TE,t,RE,'--r',t,PE,'-.k',t,KE,':c')
legend('Total Energy','Rotational Energy','Potential Energy','Kinetic Energy')
xlabel('Time(sec)')
figure
subplot(1,2,1)
n_z_forward=(z.u.*z.x3.^2.*cos(z.alpha)+0.01280*z.x1.*sqrt(z.x1.^2+z.x2.^2))/0.003043;
plot(t,n_z_forward);xlabel('Time (Sec)');ylabel('n_z(Load Factor)')
subplot(1,2,2)
Th=rho*A*omega.^2*R^2.*z.u;
plot(t,Th);xlabel('Time (Sec)');ylabel('Thrust')
figure
subplot(1,2,1)
plot(t,z.x4);xlabel('time(sec)');ylabel('\lambda_i')
subplot(1,2,2)
lambda=z.x4+(z.x2.*sin(z.alpha)-z.x1.*cos(z.alpha))./z.x3;
plot(t,lambda);xlabel('time(sec)');ylabel('\lambda')
figure
subplot(1,2,1)
plot(t,z.u);xlabel('time(sec)');ylabel('C_T')
subplot(1,2,2)
plot(z.x5,z.u);xlabel('height');ylabel('C_T');set(gca,'Xdir','reverse')
figure
height=z.x5*h1; 
range=z.x6*h1;
plot(range,height);xlabel('Range');ylabel('Height')
grid on
hold on
a=2.5; % horizontal radius
b=35; % vertical radius
x0=40; % x0,y0 ellipse centre coordinates
y0=35;
tt=-pi:0.01:pi;
xxx=x0+a*cos(tt);
yyy=y0+b*sin(tt);
plot(xxx,yyy,'Linewidth',2)
hold off
% axis equal
figure
subplot(3,2,1)
plot(t_flight1,z.x1*(37*5.37),'--','Linewidth',2)
xlabel('time');ylabel('sink rate');grid on;
subplot(3,2,2)
plot(t_flight1,z.x2*(37*5.37),'--','Linewidth',2)
xlabel('time');ylabel('forward veocity');grid on ;
subplot(3,2,3)
plot(t_flight1,z.x3*37,'--','Linewidth',2)
xlabel('time');ylabel('\omega/\omega_0');grid on;
subplot(3,2,4)
plot(t_flight1,z.x5*h1,'-.','Linewidth',2);%set(gca,'Xdir','reverse')
xlabel('time')
ylabel('height')
grid on 
subplot(3,2,5)
plot(t_flight1,z.alpha*57.29);xlabel('Time(sec)');ylabel('Alpha')
subplot(3,2,6)
plot(t_flight1,z.u,'Linewidth',2);xlabel('time(sec)');ylabel('Ct'); %set(gca,'Xdir','reverse')

end
end
%  (((z.x6-(40/h)).^2/(1/h)^2) +( (z.x5-(35/h)).^2/(35/h)^2))
