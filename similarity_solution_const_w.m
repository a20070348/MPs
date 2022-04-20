clear;
clc;

Sc=0.5;us=0.5;k=0.41;%constant paramters
L=-30;%Obukhov length L for stability of boundary layer
Cr1=2;Cr2=2;Cr0=2;%reference concentration
z=0:1:100;
ws=0.0641;ws0=0.45;%settling velocity for case with extremly low dissipation rate and spherical particle
alpha=ws*Sc/(k*us);%Rosse number for case with extremly low dissipation rate
alpha0=ws0*Sc/(k*us);%Rosse number for spherical particle
dissip=((1+0.5*(-z/L).^(2/3)).^1.5*us^3/k)./z;%dissipation rate formulated from observations
FL=0.005;el=1.5*10^(-5);asp=250;%parameters for fibers
SF=5*0.1*ws^2/(8*el*log(asp))*((FL^2)./dissip).^(1/3);
if SF<0.1
   ws1=ws*4/3; 
elseif SF>5.0
    ws1=ws+ws*2/15/(SF.^2);
else
    ws1=ws+ws*(0.07531*SF.^(-0.6692)-0.0188);
end
alpha1=ws1*Sc/(k*us);
tau=16*z/L*(0.5-alpha)-(1+alpha);
y=2*alpha./(sqrt(tau.^2-64*alpha*z/L*(0.5+alpha))-tau);
r=y.^2/alpha+(((1-y).^2-(y.^2).*(1-y).^2*0.5).*(16*z/L).^2)./(alpha*(1-(16*z/L).*y).^2);
C=-((((1+alpha)^(0.5+alpha)*r.^(-0.5)).*(y/alpha).^alpha).*(1-y)).*(1-(16*z/L).*y).^(-0.5);
C00=(C).*(1./(z.^alpha))-C;
C=(C+Cr1).*(1./(z.^alpha))-C;
tau0=16*z/L*(0.5-alpha0)-(1+alpha0);
y0=2*alpha0./(sqrt(tau0.^2-64*alpha0*z/L.*(0.5+alpha0))-tau0);
r0=y.^2/alpha0+(((1-y0).^2-(y0.^2).*(1-y0).^2*0.5).*(16*z/L).^2)./(alpha0*(1-(16*z/L).*y0).^2);
C0=-((((1+alpha0)^(0.5+alpha0)*r0.^(-0.5)).*(y/alpha0).^alpha0).*(1-y0)).*(1-(16*z/L).*y0).^(-0.5);
C01=(C0).*(1./(z.^alpha0))-C0;
C0=(C0+Cr0).*(1./(z.^alpha0))-C0;
tau1=(16*z/L).*(0.5-alpha1)-(1+alpha1);
y1=2*alpha1./(sqrt(tau1.^2-(64*alpha1.*z/L).*(0.5+alpha1))-tau1);
r1=(y.^2)./alpha1+(1-y1).^2-((((y1.^2).*(1-y1)).^2*0.5).*((16*z/L)).^2)./((alpha1.*(1-(16*z/L).*y1)).^2);
C1=-(((1+alpha1).^(0.5+alpha1)).*(r1.^(-0.5))).*((((y1./alpha1).^alpha1).*(1-y1)).*(1-(16*z/L).*y1).^(-0.5));
C02=(C1).*(1./(z.^alpha1))-C1;
C1=(C1+Cr2).*(1./(z.^alpha1))-C1;



% plot(z,C01,'red')
% hold on;
% text(50,0.1,'sphere')
% hold on;
% plot(z,C00,'blue')
% hold on
% text(50,1.15,'fiber')
% hold on;
% plot(z,C02,'black')
% hold on
% text(50,0.85,'fiber local')
plot(z,log(dissip),'green')