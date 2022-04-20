clear;
clc;

nu = 1.5* 10^(-5);mu = 1.8* 10^(-5);eulerc = 0.5772156649;dissip = 0.001;
rho = 1000.0;g = 9.8; costheta1 = 0; costheta2 = 1.0/3.0;n1=100;n2=250;n21=1;
beta1 =1:1:n1; beta2 =n21:1:n2;
% D = 1.0* 10^(-5):1.0* 10^(-6):5.0* 10^(-5);
D1 = 1.0* 10^(-5);
D2 = 2.5* 10^(-5);
Dsp = D2*(1.5*beta2).^(1.0/3.0);
Dspd = D2*(1.5*beta2).^(1.0/3.0)*10^(6);
% b=1.0*10^(-5);l=5.0*10^(-5):1.0*10^(-6):2.5*10^(-3);
wshsRe = rho*g/(16*mu)*((log(beta1)+log(4)-0.5)+costheta1*(log(beta1)+log(4)-2.5));
wsransRe = rho*g/(16*mu)*((log(beta1)+log(4)-0.5)+costheta2*(log(beta1)+log(4)-2.5));
% wshsRe = wshsRe * D^2;
% wsransRe = wsransRe * D^2;

Res = (wshsRe.*beta1)*D1^3/nu;Resf = 0.0*beta1;% this may not be right,missiong factor of 2 since l instead of L
i=1;
while (i<=n1)
    while (abs(Res(i)-Resf(i))>0.1)
        Resf(i) = Res(i);  
        wshsRe(i) = rho*g/(16*mu)*((log(beta1(i))+log(4)-0.5-0.5*Res(i))+costheta1*(log(beta1(i))+log(4)-2.5));
        wsransRe(i) = rho*g/(16*mu)*((log(beta1(i))+log(4)-0.5-0.5*Res(i))+costheta2*(log(beta1(i))+log(4)-2.5));
        Res(i) = (wshsRe(i)*beta1(i))*D1^3/nu;
    end
    i=i+1;
end
wshsRe = wshsRe * D1^2;
wsransRe = wsransRe * D1^2;

syms x
f = exp(-x)/x;
Rehl = 0.0*beta2;Rehlf = 0.0*beta2;Reranl = 0.0*beta2;Reranlf = 0.0*beta2;
wshlRe = 0.0*beta2;wsranlRe = 0.0*beta2;
Fvh = 0.0*beta2;Fhh = 0.0*beta2;Mvh = 0.0*beta2;Mhh = 0.0*beta2;
Fvran = 0.0*beta2;Fhran = 0.0*beta2;Mvran = 0.0*beta2;Mhran = 0.0*beta2;
i=1;
while (i<=n2-n21+1)
%     wshlRe(i) = rho*g/(16*mu)*((log(beta2(i))+log(4)-0.5)+costheta1*(log(beta2(i))+log(4)-2.5));
%     wsranlRe(i) = rho*g/(16*mu)*((log(beta2(i))+log(4)-0.5)+costheta2*(log(beta2(i))+log(4)-2.5));
%     Rel(i) = (wshlRe(i)*beta2(i))*D2^3/mu;
    Rehl(i) = (beta2(i)*log(beta2(i)))*D2^3*rho*g/(32*mu*nu);
    Rehlf(i) = 0.0;
    Reranl(i) = (beta2(i)*log(beta2(i)))*D2^3*rho*g/(32*mu*nu)*4.0/3.0;
    Reranlf(i) = 0.0;
    while ((abs(Rehl(i)-Rehlf(i))/Rehl(i)>0.05)||(abs(Reranl(i)-Reranlf(i))/Reranl(i)>0.05))
        Rehlf(i) = Rehl(i); 
        Fvh(i) = int(f,x,Rehl(i),inf) + log(Rehl(i))-(exp(-Rehl(i))-1)/Rehl(i) + eulerc -0.5 -log(4);
        Fhh(i) = 0.5*((int(f,x,2*Rehl(i),inf)+log(2*Rehl(i))-exp(-2*Rehl(i))+eulerc+1)/(2.0*Rehl(i))+int(f,x,2*Rehl(i),inf)+log(Rehl(i))+eulerc-3*log(2)+1);
        Mvh(i) = log(beta2(i))-Fvh(i);
        Mhh(i) = 2.0*log(beta2(i))-2.0*Fhh(i);
        wshlRe(i) = rho*g/(16*mu)*(Mvh(i)+costheta1*(Mhh(i)-Mvh(i)));        
        Rehl(i) = (wshlRe(i)*beta2(i))*D2^3/(2.0*nu);
        
        Reranlf(i) = Reranl(i);
        Fvran(i) = int(f,x,Reranl(i),inf) + log(Reranl(i))-(exp(-Reranl(i))-1)/Reranl(i) + eulerc -0.5 -log(4);
        Fhran(i) = 0.5*((int(f,x,2*Reranl(i),inf)+log(2*Reranl(i))-exp(-2*Reranl(i))+eulerc+1)/(2.0*Reranl(i))+int(f,x,2*Reranl(i),inf)+log(Reranl(i))+eulerc-3*log(2)+1);
        Mvran(i) = log(beta2(i))-Fvran(i);
        Mhran(i) = 2.0*log(beta2(i))-2.0*Fhran(i);
        wsranlRe(i) = rho*g/(16*mu)*(Mvran(i)+costheta2*(Mhran(i)-Mvran(i)));
        Reranl(i) = (wsranlRe(i)*beta2(i))*D2^3/(2.0*nu);
    end
    i=i+1;
end
wshlRe = wshlRe * D2^2;
wsranlRe = wsranlRe * D2^2;

wshlRed = wshlRe;
wsranlRed = wsranlRe;

wshlRe = wshlRe/(sqrt(D2*g));
wsranlRe = wsranlRe/(sqrt(D2*g));

Mhratio = (Mhh-Mvh)./Mvh;
Mranratio = (Mhran-Mvran)./Mvran;

plot(beta2,Mhratio,'--black','linewidth',3)
hold on
plot(beta2,Mranratio,'-.green','linewidth',3)
% while (abs(Re1-Ref)>0.00000001)
%    Ref = Re1; 
%   
%    Fv = int(f,x,Re1,inf) + log(Re1)-(exp(-Re1)-1)./Re1 + eulerc -0.5 -log(4);
%    Fh = 0.5*((int(f,x,2*Re1,inf)+log(2*Re1)-exp(-2*Re1)-eulerc+1)./(2.0*Re1)+int(f,x,2*Re1,inf)+log(Re1)+eulerc-3*log(2)+1);
%    Mv = log(beta)-Fv;
%    Mh = 2.0*log(beta)-2.0*Fh;
%    wsh2 = rho*g/(16*mu)*(Mv+costheta1*(Mh-Mv));
%    wsran2 = rho*g/(16*mu)*(Mv+costheta2*(Mh-Mv));
%    Re1 = (wsh2.*beta).*D.^3/mu;  
% end
% plot(beta1,wshsRe,'red')
% hold on
% plot(beta1,wsransRe,'blue')
% hold on
% figure(1)
% plot(beta2,wshlRe,'--black','linewidth',3)
% hold on
% plot(beta2,wsranlRe,'-.green','linewidth',3)
% hold on
% plot(beta2,wssp,'magenta','linewidth',3)
% set(gca,'FontSize',24);
% title(' Normalized Settling Velocity','fontsize',18)
% xlabel('$\beta$','fontsize',36,'Interpreter','latex') 
% ylabel('$\frac{w_s}{\sqrt {g*D}}$','fontsize',36,'Interpreter','latex')
% legend({'horizontally aligned fiber','randomly aligned fiber','spherical particle'},'Location','northwest','fontsize',18)
% 
% figure(2)
% 
% plot(Dspd,wshlRed,'--black','linewidth',3)
% hold on
% plot(Dspd,wsranlRed,'-.green','linewidth',3)
% hold on
% plot(Dspd,wsspd,'magenta','linewidth',3)
% set(gca,'FontSize',24);
% title('Settling Velocity','fontsize',18)
% xlabel('$D (\mu m)$','fontsize',36,'Interpreter','latex') 
% ylabel('$w_s (m/s)$','fontsize',36,'Interpreter','latex')
% legend({'horizontally aligned fiber','randomly aligned fiber','spherical particle'},'Location','northwest','fontsize',18)
