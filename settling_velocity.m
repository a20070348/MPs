clear;
clc;

nu = 1.5* 10^(-5);mu = 1.8* 10^(-5);eulerc = 0.5772156649;dissip = 0.0001;
rho = 1000.0;g = 9.8;n1=100;n2=150;n21=50;
beta1 =1:1:n1; beta2 =n21:1:n2;
% D = 1.0* 10^(-5):1.0* 10^(-6):5.0* 10^(-5);
D1 = 1.0* 10^(-5);
D2 = 2.0* 10^(-5);
Dsp = D2*(1.5*beta2).^(1.0/3.0);
Dspd = D2*(1.5*beta2).^(1.0/3.0)*10^(6);
% b=1.0*10^(-5);l=5.0*10^(-5):1.0*10^(-6):2.5*10^(-3);

syms x
f = exp(-x)/x;
Rel = 0.0*beta2;Relf = 0.0*beta2;
wslRe = 0.0*beta2;wslRemin = 0.0*beta2;
Fv = 0.0*beta2;Fh = 0.0*beta2;Mv = 0.0*beta2;Mh = 0.0*beta2;
Fic = 0.0*beta2;
SF = 0.0*beta2;costheta = 0.0*beta2+0.5;theta = 0.0*beta2;
i=1;
while (i<=n2-n21+1)
    Rel(i) = beta2(i)*(log(2.0*beta2(i))+log(4)-0.5)*D2^3*rho*g/(32*mu*nu);
    Relf(i) = 0.0;
    while (abs(Rel(i)-Relf(i))/Rel(i)>0.05)
        Relf(i) = Rel(i); 
        Fv(i) = int(f,x,Rel(i),inf) + log(Rel(i))-(exp(-Rel(i))-1)/Rel(i) + eulerc -0.5 -log(4);
        Fh(i) = 0.5*((int(f,x,2*Rel(i),inf)+log(2*Rel(i))-exp(-2*Rel(i))+eulerc+1)/(2.0*Rel(i))+int(f,x,2*Rel(i),inf)+log(Rel(i))+eulerc-3*log(2)+1);
        Mv(i) = log(2.0*beta2(i))-Fv(i);
        Mh(i) = 2.0*log(2.0*beta2(i))-2.0*Fh(i);
        wslRemin(i) = D2^2*rho*g/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i))); 
        
        theta(i) = acos(costheta(i)^0.5);
        
        Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
        
        if beta2(i)*D2<=(nu^3/dissip)^0.25
           SF(i) = 5.0*wslRemin(i)^2/(8.0*nu^0.5*dissip^0.5*log(2.0*beta2(i)));
        else
           SF(i) = 5.0*wslRemin(i)^2*(beta2(i)*D2)^(2.0/3.0)/(8.0*nu*dissip^(1.0/3.0)*log(2.0*beta2(i)));
        end
        
        if SF(i)<=0.1
           costheta(i) = 1.0/3.0;
        elseif SF(i)>=5.0
           costheta(i) = 2.0/(15.0*SF(i)^2);
        else
           costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
        end
        
        wslRe(i) = rho*g/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
        Rel(i) = (wslRe(i)*beta2(i))*D2^3/(2.0*nu);
    end
    i=i+1;
end
wslRe = wslRe * D2^2;

figure(1)
plot(beta2,wslRe,'--black','linewidth',3)
figure(2)
plot(beta2,SF,'-.green','linewidth',3)
figure(3)
plot(beta2,costheta,'magenta','linewidth',3)
figure(4)
plot(beta2,Fic,'red','linewidth',3)
figure(5)
plot(beta2,costheta.*(Mh-Mv)./Mv,'blue','linewidth',3)
figure(6)
plot(beta2,(Mh-Mv)./Mv,'blue','linewidth',3)
% set(gca,'FontSize',24);
% title(' Normalized Settling Velocity','fontsize',18)
% xlabel('$\beta$','fontsize',36,'Interpreter','latex') 
% ylabel('$\frac{w_s}{\sqrt {g*D}}$','fontsize',36,'Interpreter','latex')
% legend({'horizontally aligned fiber','randomly aligned fiber','spherical particle'},'Location','northwest','fontsize',18)
