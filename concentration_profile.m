Sc=0.5;us=0.5;k=0.41;L=-30;Cr=15;Cr0=2;Cr1=15;zr=1;
z=0:1:100;
ws=0.0641;ws0=0.45;
alpha=ws*Sc/(k*us);
alpha0=ws0*Sc/(k*us);
dissip=((1+0.5*(-z/L).^(2/3)).^1.5*us^3/k)./z;
FL=0.005;el=1.5*10^(-5);asp=250;
SF=5*0.1*ws^2/(8*el*log(asp))*((FL^2)./dissip).^(1/3);
if SF<0.1
   ws1=ws*4/3; 
elseif SF>5.0
    ws1=ws+ws*2/15/(SF.^2);
else
    ws1=ws+ws*(0.07531*SF.^(-0.6692)-0.0188);
end
alpha1=ws1*Sc/(k*us);

tau=(16*z/L).*(0.5-alpha)-(1+alpha);
y=2*alpha./(sqrt(tau.^2-64*alpha*z/L*(0.5+alpha))-tau);
r=y.^2/alpha+(1-y).^2-(((y.^2).*(1-y).^2*0.5).*(16*z/L).^2)./(alpha*(1-(16*z/L).*y).^2);
wm=-((((1+alpha)^(0.5+alpha)*r.^(-0.5)).*(y/alpha).^alpha).*(1-y)).*(1-(16*z/L).*y).^(-0.5);
C=(-wm(2)/ws+(k/ws+(1-16*zr/L)^(-0.5))/(k+ws)).*((zr./z).^alpha)+wm/ws;

tau0=(16*z/L).*(0.5-alpha0)-(1+alpha0);
y0=2*alpha0./(sqrt(tau0.^2-64*alpha0*z/L*(0.5+alpha0))-tau0);
r0=y0.^2/alpha0+(1-y0).^2-(((y0.^2).*(1-y0).^2*0.5).*(16*z/L).^2)./(alpha0*(1-(16*z/L).*y0).^2);
wm0=-((((1+alpha0)^(0.5+alpha0)*r0.^(-0.5)).*(y0/alpha0).^alpha0).*(1-y0)).*(1-(16*z/L).*y0).^(-0.5);
C0=(-wm0(2)/ws0+(k/ws0+(1-16*zr/L)^(-0.5))/(k+ws0)).*((zr./z).^alpha0)+wm0/ws0;

tau1=(16*z/L).*(0.5-alpha1)-(1+alpha1);
y1=2*alpha1./(sqrt(tau1.^2-(64*alpha1.*z/L).*(0.5+alpha1))-tau1);
r1=(y1.^2)./alpha1+(1-y1).^2-(((y1.^2).*(1-y1).^2*0.5).*(16*z/L).^2)./(alpha1.*(1-(16*z/L).*y1).^2);
wm1=-(((1+alpha1).^(0.5+alpha1)).*r1.^(-0.5)).*((((y1./alpha1).^alpha1).*(1-y1)).*(1-(16*z/L).*y1).^(-0.5));
C1=(-wm1(2)/ws1(2)+(k/ws1(2)+(1-16*zr/L)^(-0.5))/(k+ws1(2))).*((zr./z).^alpha1)+wm1./ws1;

plot(z,C0,'red')
hold on;
text(50,0,'sphere')
hold on;
plot(z,C,'blue')
hold on
text(50,5.5,'fiber')
hold on;
plot(z,C1,'black')
hold on
text(50,2.5,'fiber local')