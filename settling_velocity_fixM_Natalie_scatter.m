clear;
clc;

nu = 1.5* 10^(-5);mu = 1.8* 10^(-5);eulerc = 0.5772156649;dissip = 0.0001;
rho = 1000.0;g = 9.8;

geodata1=load('/Users/shuolinxiao/Downloads/PostDoc/Cornell_Qi/projects/MPF/code/Beaver_Mountain_Fiber.dat');
geodata2=load('/Users/shuolinxiao/Downloads/PostDoc/Cornell_Qi/projects/MPF/code/CO98_Fiber.dat');
geodata3=load('/Users/shuolinxiao/Downloads/PostDoc/Cornell_Qi/projects/MPF/code/UT95_Fiber.dat');
% geodata4=[geodata1;geodata2;geodata3];
geodata4=load('/Users/shuolinxiao/Downloads/PostDoc/Cornell_Qi/projects/MPF/code/sample.dat');

geodata1=10^(-6)*geodata1;
geodata2=10^(-6)*geodata2;
geodata3=10^(-6)*geodata3;
geodata4=10^(-6)*geodata4;

width=geodata4(:,2);height = 0.2* 10^(-5);
radius=(0.25*width.^2+height^2)/(2.0*height);
De = sqrt(asin(0.5*width./radius)*4.0.*radius.^2/pi-2.0*width.*(radius-height)/pi);

betafl = geodata4(:,1)./De;
betacy = geodata4(:,1)./width;
Dfl = De;
Dcy = width;

Dim = size(width);

syms x
f = exp(-x)/x;
Rel = 0.0*width;Relf = 0.0*width;
wslRe = 0.0*width;wslRemin = 0.0*width;
wslRecfl = 0.0*width;wslReccy = 0.0*width;
Fv = 0.0*width;Fh = 0.0*width;Mv = 0.0*width;Mh = 0.0*width;
Fic = 0.0*width;
SF = 0.0*width;costheta = 0.0*width+0.5;theta = 0.0*width;thetaf = 0.0*width;

i=1;
while (i<=Dim(:,1))
    Rel(i) = betafl(i)*(log(2.0*betafl(i))+log(4)-0.5)*Dfl(i)^3*rho*g/(32*mu*nu);
    Relf(i) = 0.0;
    while (abs(Rel(i)-Relf(i))/Rel(i)>0.01)
        Relf(i) = Rel(i); 
        Fv(i) = int(f,x,Rel(i),inf) + log(Rel(i))-(exp(-Rel(i))-1)/Rel(i) + eulerc -0.5 -log(4);
        Fh(i) = 0.5*((int(f,x,2*Rel(i),inf)+log(2*Rel(i))-exp(-2*Rel(i))+eulerc+1)/(2.0*Rel(i))+int(f,x,2*Rel(i),inf)+log(Rel(i))+eulerc-3*log(2)+1);
        Mv(i) = log(2.0*betafl(i))-Fv(i);
        Mh(i) = 2.0*log(2.0*betafl(i))-2.0*Fh(i);
        wslRemin(i) = rho*g*Dfl(i)^2/(16*mu)*Mv(i); %incomning flow velocity as U in COX
        Rel(i) = wslRemin(i)*betafl(i)*Dfl(i)/(2.0*nu);
    end
    i=i+1;
end    

i=1;
while (i<=Dim(:,1))
    
    theta(i) = pi/4.0;
    
    Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
    
    if betafl(i)*Dfl(i)<=(nu^3/dissip)^0.25
        SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip^0.5*log(2.0*betafl(i)));
    else
        SF(i) = 5.0*wslRemin(i)^2*(betafl(i)*Dfl(i))^(2.0/3.0)*Fic(i)/(8.0*nu*dissip^(1.0/3.0)*log(2.0*betafl(i)));
    end
        
    if SF(i)<=0.1
        costheta(i) = 1.0/3.0;
    elseif SF(i)>=5.0
        costheta(i) = 2.0/(15.0*SF(i)^2);
    else
        costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
    end
       
    wslRecfl(i) = rho*g*Dfl(i)^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
 
    i=i+1;
end     

i=1;
while (i<=Dim(:,1))
    Rel(i) = betacy(i)*(log(2.0*betacy(i))+log(4)-0.5)*Dcy(i)^3*rho*g/(32*mu*nu);
    Relf(i) = 0.0;
    while (abs(Rel(i)-Relf(i))/Rel(i)>0.01)
        Relf(i) = Rel(i); 
        Fv(i) = int(f,x,Rel(i),inf) + log(Rel(i))-(exp(-Rel(i))-1)/Rel(i) + eulerc -0.5 -log(4);
        Fh(i) = 0.5*((int(f,x,2*Rel(i),inf)+log(2*Rel(i))-exp(-2*Rel(i))+eulerc+1)/(2.0*Rel(i))+int(f,x,2*Rel(i),inf)+log(Rel(i))+eulerc-3*log(2)+1);
        Mv(i) = log(2.0*betacy(i))-Fv(i);
        Mh(i) = 2.0*log(2.0*betacy(i))-2.0*Fh(i);
        wslRemin(i) = rho*g*Dcy(i)^2/(16*mu)*Mv(i); %incomning flow velocity as U in COX
        Rel(i) = wslRemin(i)*betacy(i)*Dcy(i)/(2.0*nu);
    end
    i=i+1;
end    

i=1;
while (i<=Dim(:,1))
    
    theta(i) = pi/4.0;
    
    Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
    
    if betacy(i)*Dcy(i)<=(nu^3/dissip)^0.25
        SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip^0.5*log(2.0*betacy(i)));
    else
        SF(i) = 5.0*wslRemin(i)^2*(betacy(i)*Dcy(i))^(2.0/3.0)*Fic(i)/(8.0*nu*dissip^(1.0/3.0)*log(2.0*betacy(i)));
    end
        
    if SF(i)<=0.1
        costheta(i) = 1.0/3.0;
    elseif SF(i)>=5.0
        costheta(i) = 2.0/(15.0*SF(i)^2);
    else
        costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
    end
       
    wslReccy(i) = rho*g*Dcy(i)^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
 
    i=i+1;
end 

Resultfl0 = [width,wslRecfl];
Resultfl = sortrows(Resultfl0);

Resultcy0 = [width,wslReccy];
Resultcy = sortrows(Resultcy0);

Result=[wslRecfl*100,wslReccy*100];

wratio=wslReccy./wslRecfl;

h1 = histfit(wslRecfl*100,10,'kernel');
delete(h1(1))
set(h1(2),'color','m')
set(h1(2),'LineStyle','--')
hold on
h2 = histfit(wslReccy*100,100,'kernel');
delete(h2(1))
set(h2(2),'color','r')

xlim([0 15])
set(gca,'FontSize',24);
title(' PDF of Fiber Settling Velocity','fontsize',18)
xlabel('${w_s (cm/s)}$','fontsize',36,'Interpreter','latex') 
ylabel('$PDF$','fontsize',36,'Interpreter','latex')
legend({'flat fiber','cylindrical fiber'},'Location','northeast','fontsize',18)

figure(2)
gscatter(width*10^6,geodata4(:,1)*10^6,wratio,36);
% colormap(subfig1,'jet');
set(gca,'FontSize',24);
title('scatter plot for horizontally aligned fiber','fontsize',18,'position',[100,4.5])
xlabel('width $(\mu m)$','fontsize',36,'Interpreter','latex') 
ylabel('length $(\mu m)$','fontsize',36,'Interpreter','latex')
% xlim([-100 300])

% figure(1)
% h1 = histfit(wslRecfl*100,10,'kernel');
% delete(h1(1))
% set(h1(2),'color','m')
% set(h1(2),'LineStyle','--')
% xlim([0 1.5])
% set(gca,'FontSize',24);
% title(' Curved Histograms of Flat Fiber Settling Velocity','fontsize',18)
% xlabel('${w_s (cm/s)}$','fontsize',36,'Interpreter','latex') 
% ylabel('$count$','fontsize',36,'Interpreter','latex')
% 
% figure(2)
% h2 = histfit(wslReccy*100,10,'kernel');
% delete(h2(1))
% set(h2(2),'color','r')
% hold on
% 
% xlim([0 15])
% set(gca,'FontSize',24);
% title(' Curved Histograms of Cylindrical Fiber Settling Velocity','fontsize',18)
% xlabel('${w_s (cm/s)}$','fontsize',36,'Interpreter','latex') 
% ylabel('$count$','fontsize',36,'Interpreter','latex')

% plot(D,wslRe,'--green','linewidth',3)
% hold on
% scatter(Resultfl(:,1)*10^6,Resultfl(:,2),'o','LineWidth',14)
% hold on
% scatter(Resultcy(:,1)*10^6,Resultcy(:,2),'+','LineWidth',14)
% set(gca,'FontSize',24);
% title(' Fiber Settling Velocity','fontsize',18)
% xlabel('${w (\mu m)}$','fontsize',36,'Interpreter','latex') 
% ylabel('${w_s (m/s)}$','fontsize',36,'Interpreter','latex')
% legend({'high dissipation rate','low dissipation rate'},'Location','northwest','fontsize',18)