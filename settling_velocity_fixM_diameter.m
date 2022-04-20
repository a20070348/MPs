clear;
clc;

nu = 1.5* 10^(-5);mu = 1.8* 10^(-5);eulerc = 0.5772156649;dissip = 0.001;
rho = 1000.0;g = 9.8;
beta =10;
D = 1.0* 10^(-6):1.0* 10^(-6):5.0* 10^(-5);
D1 = D*10^6;

syms x
f = exp(-x)/x;
Rel = 0.0*D;Relf = 0.0*D;
wslRe = 0.0*D;wslRec = 0.0*D;wslRemin = 0.0*D;
Fv = 0.0*D;Fh = 0.0*D;Mv = 0.0*D;Mh = 0.0*D;
Fic = 0.0*D;
SF = 0.0*D;costheta = 0.0*D+0.5;theta = 0.0*D;thetaf = 0.0*D;
i=1;
while (i<=50)
    Rel(i) = beta*(log(2.0*beta)+log(4)-0.5)*D(i)^3*rho*g/(32*mu*nu);
    Relf(i) = 0.0;
    while (abs(Rel(i)-Relf(i))/Rel(i)>0.01)
        Relf(i) = Rel(i); 
        Fv(i) = int(f,x,Rel(i),inf) + log(Rel(i))-(exp(-Rel(i))-1)/Rel(i) + eulerc -0.5 -log(4);
        Fh(i) = 0.5*((int(f,x,2*Rel(i),inf)+log(2*Rel(i))-exp(-2*Rel(i))+eulerc+1)/(2.0*Rel(i))+int(f,x,2*Rel(i),inf)+log(Rel(i))+eulerc-3*log(2)+1);
        Mv(i) = log(2.0*beta)-Fv(i);
        Mh(i) = 2.0*log(2.0*beta)-2.0*Fh(i);
        wslRemin(i) = rho*g*D(i)^2/(16*mu)*Mv(i); %incomning flow velocity as U in COX
        Rel(i) = wslRemin(i)*beta*D(i)/(2.0*nu);
    end
    i=i+1;
end

% i=1;
% while (i<=50)
%     
%     theta(i) = pi/4.0;
%     while (abs(theta(i)-thetaf(i))/theta(i)>0.001)
%         Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
%     
%         if beta*D(i)<=(nu^3/dissip)^0.25
%             SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip^0.5*log(2.0*beta));
%         else
%             SF(i) = 5.0*wslRemin(i)^2*(beta*D(i))^(2.0/3.0)*Fic(i)/(8.0*nu*dissip^(1.0/3.0)*log(2.0*beta));
%         end
%         
%         if SF(i)<=0.1
%             costheta(i) = 1.0/3.0;
%         elseif SF(i)>=5.0
%             costheta(i) = 2.0/(15.0*SF(i)^2);
%         else
%             costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
%         end
%         thetaf(i) = theta(i);
%         theta(i) = acos(costheta(i)^0.5);
%     end   
%     wslRe(i) = rho*g*D(i)^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
%  
%     i=i+1;
% end     

i=1;
while (i<=50)
    
    theta(i) = pi/4.0;
    
    Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
    
    if beta*D(i)<=(nu^3/dissip)^0.25
        SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip^0.5*log(2.0*beta));
    else
        SF(i) = 5.0*wslRemin(i)^2*(beta*D(i))^(2.0/3.0)*Fic(i)/(8.0*nu*dissip^(1.0/3.0)*log(2.0*beta));
    end
        
    if SF(i)<=0.1
        costheta(i) = 1.0/3.0;
    elseif SF(i)>=5.0
        costheta(i) = 2.0/(15.0*SF(i)^2);
    else
        costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
    end
       
    wslRec(i) = rho*g*D(i)^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
 
    i=i+1;
end     

% plot(D,wslRe,'--green','linewidth',3)
% hold on
plot(D*10^6,wslRec,'-black','linewidth',3)
title(' Settling Velocity','fontsize',18)
xlabel('$D (\mu m)$','fontsize',36,'Interpreter','latex') 
ylabel('${w_s (m/s)}$','fontsize',36,'Interpreter','latex')
% legend({'high dissipation rate','low dissipation rate'},'Location','northwest','fontsize',18)