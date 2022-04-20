clear;
clc;

nu = 1.5* 10^(-5);mu = 1.8* 10^(-5);eulerc = 0.5772156649;dissip = 0.001;
D = 1.5* 10^(-5);
D1 = 1.5;
rho = 1140.0;g = 9.8;
% beta =5:1:250;
beta = [1.75/D1,3.75/D1,7.5/D1,17.5/D1,37.5/D1,75/D1,125/D1,175/D1,225/D1,275/D1,350/D1,450/D1];

Dsp = D*(1.5*beta).^(1.0/3.0);

syms x
f = exp(-x)/x;
Rel = 0.0*beta;Relf = 0.0*beta;
wslRe = 0.0*beta;wslRemin = 0.0*beta;
Fv = 0.0*beta;Fh = 0.0*beta;Mv = 0.0*beta;Mh = 0.0*beta;
Fic = 0.0*beta;
SF = 0.0*beta;costheta = 0.0*beta;theta = 0.0*beta;

i=1;
while (i<=12)
    Rel(i) = beta(i)*(log(2.0*beta(i))+log(4)-0.5)*D^3*rho*g/(32*mu*nu);
    Relf(i) = 0.0;
    while (abs(Rel(i)-Relf(i))/Rel(i)>0.01)
        Relf(i) = Rel(i); 
        Fv(i) = int(f,x,Rel(i),inf) + log(Rel(i))-(exp(-Rel(i))-1)/Rel(i) + eulerc -0.5 -log(4);
        Fh(i) = 0.5*((int(f,x,2*Rel(i),inf)+log(2*Rel(i))-exp(-2*Rel(i))+eulerc+1)/(2.0*Rel(i))+int(f,x,2*Rel(i),inf)+log(Rel(i))+eulerc-3*log(2)+1);
        Mv(i) = log(2.0*beta(i))-Fv(i);
        Mh(i) = 2.0*log(2.0*beta(i))-2.0*Fh(i);
        wslRemin(i) = rho*g*D^2/(16*mu)*Mv(i); %incomning flow velocity as U in COX
        Rel(i) = (wslRemin(i)*beta(i))*D/(2.0*nu);
    end
    i=i+1;
end

i=1;
while (i<=12)
    
    theta(i) = pi/4.0;
    Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
    
    if beta(i)*D<=(nu^3/dissip)^0.25
        SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip^0.5*log(2.0*beta(i)));
    else
        SF(i) = 5.0*wslRemin(i)^2*(beta(i)*D)^(2.0/3.0)*Fic(i)/(8.0*nu*dissip^(1.0/3.0)*log(2.0*beta(i)));
    end
        
    if SF(i)<=0.1
        costheta(i) = 1.0/3.0;
    elseif SF(i)>=5.0
        costheta(i) = 2.0/(15.0*SF(i)^2);
    else
        costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
    end
    wslRe(i) = rho*g*D^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
 
    i=i+1;
end    

Resp = Dsp.^3*rho*g/(18*mu*nu);Respf = 0.0*beta;wssp = 0.0*beta;
i=1;
while (i<=12)
    while (abs(Resp(i)-Respf(i))/Resp(i)>0.01)
        Respf(i) = Resp(i);  
        wssp(i) = Dsp(i)^2*rho*g/(18*mu)/(1+0.15*Resp(i)^0.687);
        Resp(i) = wssp(i)*Dsp(i)/nu;
    end
    i=i+1;
end

Dz1=sqrt(wslRe*18*nu/rho/g);
Dz2=sqrt(wssp*18*nu/rho/g);

scatter(beta*D1*10,Dz1*10^6,'o','LineWidth',14)
hold on
scatter(beta*D1*10,Dz2*10^6,'+','LineWidth',14)
hold on
set(gca,'FontSize',24);
title('Equavelent Spherical Size','fontsize',18)
xlabel('fiber length $(\mu m)$','fontsize',36,'Interpreter','latex') 
ylabel('equivalent spherical size${(\mu m)}$','fontsize',36,'Interpreter','latex')
legend({'MPF','spherical particle with equivalent volume'},'Location','northwest','fontsize',24,'Interpreter','latex')
% plot(Dsp*10^6,wslRe,'-r','linewidth',3)
% hold on
% plot(Dsp*10^6,wssp,'-.b','linewidth',3)
% xlim([30 100])
% % ylim([0 0.1])
% set(gca,'FontSize',24);
% title(' Settling Velocity','fontsize',18)
% xlabel('equaivalent diameter $D_{ev} (\mu m)$ based on MPF size','fontsize',36,'Interpreter','latex') 
% ylabel('${w_s (m/s)}$','fontsize',36,'Interpreter','latex')
% legend({'MPF','spherical particle with equivalent volume'},'Location','northwest','fontsize',24,'Interpreter','latex')