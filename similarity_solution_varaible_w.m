clear;
clc;

Sc=0.5;us=0.5;k=0.41;%constant paramters
L=-186.5;%Obukhov length L for stability of boundary layer
Chr=171;Cspr=170;Clr=169;%reference concentration
ns=1;ne=100;
z=ns:1:ne;
z1=ns:0.001:ne;
wsh=0.0641;wssp=0.45;%settling velocity for case with extremly low dissipation rate and spherical particle
deltaw=0.0634;

dd=105.0;
dissip=((1+0.5*(-z/L).^(dd)).^(1.0/dd)*us^3/k)./z;%dissipation rate formulated from observations
dissip1=-3.319e-10*z1.^3 + 7.473e-08*z1.^2 + -5.654e-06*z1 + 0.0002844;
dissip0=((1+0.5*(-z1/L).^(dd)).^(1.0/dd)*us^3/k)./z1;
FL=0.005;el=1.5*10^(-5);asp=250;%parameters for fibers
SF=5*0.1*wsh^2/(8*el*log(asp))*((FL^2)./dissip).^(1/3);
SF1=5*0.1*wsh^2/(8*el*log(asp))*((FL^2)./dissip1).^(1/3);
if SF<0.1
   wsl=wsh*4/3; 
elseif SF>5.0
    wsl=wsh+deltaw*2/15/(SF.^2);
else
    wsl=wsh+deltaw*(0.07531*SF.^(-0.6692)-0.0188);
end

if SF1<0.1
   wsl1=wsh*4/3; 
elseif SF1>5.0
    wsl1=wsh+deltaw*2/15/(SF1.^2);
else
    wsl1=wsh+deltaw*(0.07531*SF1.^(-0.6692)-0.0188);
end

alphah=wsh*Sc/(k*us);%Rosse number for case with extremly low dissipation rate
tauh=16*z/L*(0.5-alphah)-(1+alphah);
yh=2*alphah./(sqrt(tauh.^2-64*alphah*z/L*(0.5+alphah))-tauh);
rh=yh.^2/alphah+(1-yh).^2-(((yh.^2).*(1-yh).^2*0.5).*(16*z/L).^2)./(alphah*(1-(16*z/L).*yh).^2);
Wh=((((1+alphah)^(0.5+alphah)*rh.^(-0.5)).*(yh/alphah).^alphah).*(1-yh)).*(1-(16*z/L).*yh).^(-0.5);
Wh=Wh/wsh;
Ch=(Wh(1)+Chr).*(1./(z.^alphah))-Wh;

alphasp=wssp*Sc/(k*us);%Rosse number for spherical particle
tausp=16*z/L*(0.5-alphasp)-(1+alphasp);
ysp=2*alphasp./(sqrt(tausp.^2-64*alphasp*z/L.*(0.5+alphasp))-tausp);
rsp=ysp.^2/alphasp+(1-ysp).^2-(((ysp.^2).*(1-ysp).^2*0.5).*(16*z/L).^2)./(alphasp*(1-(16*z/L).*ysp).^2);
Wsp=((((1+alphasp)^(0.5+alphasp)*rsp.^(-0.5)).*(ysp/alphasp).^alphasp).*(1-ysp)).*(1-(16*z/L).*ysp).^(-0.5);
Wsp=Wsp/wssp;
Csp=(Wsp(1)+Cspr).*(1./(z.^alphasp))-Wsp;

alphal=wsl*Sc/(k*us);
taul=(16*z/L).*(0.5-alphal)-(1+alphal);
yl=2*alphal./(sqrt(taul.^2-(64*alphal.*z/L).*(0.5+alphal))-taul);
rl=(yl.^2)./alphal+(1-yl).^2-((((yl.^2).*(1-yl)).^2*0.5).*((16*z/L)).^2)./((alphal.*(1-(16*z/L).*yl)).^2);
Wl=(((1+alphal).^(0.5+alphal)).*(rl.^(-0.5))).*((((yl./alphal).^alphal).*(1-yl)).*(1-(16*z/L).*yl).^(-0.5));
Wl=Wl./wsl;
Cl=(Wl(1)+Clr).*(1./(z.^alphal))-Wl;

phic=(1-16.0*z1/L).^(-0.5);
alphal1=wsl1*Sc/(k*us);
CFDl=0.0*z1;
CFDl(1)=Clr;
i=1;
while (i<=(ne-ns)*1000)   
    CFDl(i+1)=CFDl(i)-(z1(2)-z1(1))*(phic(i)*Sc/(k*us)/z1(i)+CFDl(i)/z1(i)*alphal1(i));
    i=i+1;
end


phic=(1-16.0*z1/L).^(-0.5);
alphah1=wsh*Sc/(k*us);
CFDh=0.0*z1;
CFDh(1)=Chr;
i=1;
while (i<=(ne-ns)*1000)   
    CFDh(i+1)=CFDh(i)-(z1(2)-z1(1))*(phic(i)*Sc/(k*us)/z1(i)+CFDh(i)/z1(i)*alphah1);
    i=i+1;
end

phic=(1-16.0*z1/L).^(-0.5);
alphasp1=wssp*Sc/(k*us);
CFDsp=0.0*z1;
CFDsp(1)=Cspr;
i=1;
while (i<=(ne-ns)*1000)   
    CFDsp(i+1)=CFDsp(i)-(z1(2)-z1(1))*(phic(i)*Sc/(k*us)/z1(i)+CFDsp(i)/z1(i)*alphasp1);
    i=i+1;
end

CLESl=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/CZCClV3.txt');
CLESh=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/CZCChV3.txt');
DissipLES=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/dissipationzcc.txt');
DissipLES1 = DissipLES(:,1);
DissipLES2 = DissipLES(:,2);

% plot(Csp,z,'-red','linewidth',3)
% hold on;
% plot(Ch,z,'-.blue','linewidth',3)
% hold on;
% plot(Cl,z,'black')
% hold on

plot(CFDsp,z1,':r','linewidth',3)
hold on;
% plot(CFDl,z1,'-green','linewidth',3)
% hold on;
plot(CFDh,z1,'--m','linewidth',3)
hold on;
% 

scatter(CLESh(:,1),CLESh(:,2),'o','LineWidth',14)
hold on
scatter(CLESl(:,1),CLESl(:,2),'+','LineWidth',14)
hold on
xlim([0 170])
ylim([0 100])

set(gca,'FontSize',24);
title('Fiber Concentration','fontsize',18)
xlabel('$C (kg/m^3)$','fontsize',36,'Interpreter','latex')
ylabel('$z (m)$','fontsize',36,'Interpreter','latex') 
legend({'spherical particle with equivalent volume as MPF','MPF with settling velocity for regional model',sprintf('LES result for MPF with settling velocity for \nregional model'),sprintf('LES result for MPF with settling velocity denpending \non local dissipation rate')},'Location','northeast','fontsize',24)

% plot(log10(dissip0),z1,'green','linewidth',3)
% hold on
% scatter(log10(DissipLES(:,1)),DissipLES(:,2),'o','LineWidth',3)
% set(gca,'FontSize',24);
% title('Vertical Dissipation Rate','fontsize',18)
% xlabel('$log(\epsilon) (m^2/s^3)$','fontsize',36,'Interpreter','latex')
% ylabel('$z (m)$','fontsize',36,'Interpreter','latex') 
% legend({sprintf('dissipation rate based on the formula  \nfitted from Kansas experiment'),sprintf('dissipation rate from LES results')},'Location','northeast','fontsize',24)
