clear;
clc;

nu = 1.5* 10^(-5);mu = 1.8* 10^(-5);eulerc = 0.5772156649;dissip1 = 0.01;dissip2=0.001;dissip3=0.0001;
rho = 1140.0;g = 9.8;
beta =5:1:250;
D1 = 0.2* 10^(-5);
D2 = 1.5* 10^(-5);

syms x
f = exp(-x)/x;
Rel = 0.0*beta;Relf = 0.0*beta;
wslRemin = 0.0*beta;
wslRec11 = 0.0*beta;wslRec12 = 0.0*beta;wslRec13 = 0.0*beta;wslRec21 = 0.0*beta;wslRec22 = 0.0*beta;wslRec23 = 0.0*beta;
Fv = 0.0*beta;Fh = 0.0*beta;Mv = 0.0*beta;Mh = 0.0*beta;
Fic = 0.0*beta;
SF = 0.0*beta;costheta = 0.0*beta;theta = 0.0*beta;thetaf = 0.0*beta;

i=1;
while (i<=250-5+1)
    Rel(i) = beta(i)*(log(2.0*beta(i))+log(4)-0.5)*D1^3*rho*g/(32*mu*nu);
    Relf(i) = 0.0;
    while (abs(Rel(i)-Relf(i))/Rel(i)>0.01)
        Relf(i) = Rel(i); 
        Fv(i) = int(f,x,Rel(i),inf) + log(Rel(i))-(exp(-Rel(i))-1)/Rel(i) + eulerc -0.5 -log(4);
        Fh(i) = 0.5*((int(f,x,2*Rel(i),inf)+log(2*Rel(i))-exp(-2*Rel(i))+eulerc+1)/(2.0*Rel(i))+int(f,x,2*Rel(i),inf)+log(Rel(i))+eulerc-3*log(2)+1);
        Mv(i) = log(2.0*beta(i))-Fv(i);
        Mh(i) = 2.0*log(2.0*beta(i))-2.0*Fh(i);
        wslRemin(i) = rho*g*D1^2/(16*mu)*Mv(i); %incomning flow velocity as U in COX
        Rel(i) = wslRemin(i)*beta(i)*D1/(2.0*nu);
    end
    i=i+1;
end    

i=1;
while (i<=250-5+1)
    
    theta(i) = pi/4.0;
    
    Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
    
    if beta(i)*D1<=(nu^3/dissip1)^0.25
        SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip1^0.5*log(2.0*beta(i)));
    else
        SF(i) = 5.0*wslRemin(i)^2*(beta(i)*D1)^(2.0/3.0)*Fic(i)/(8.0*nu*dissip1^(1.0/3.0)*log(2.0*beta(i)));
    end
        
    if SF(i)<=0.1
        costheta(i) = 1.0/3.0;
    elseif SF(i)>=5.0
        costheta(i) = 2.0/(15.0*SF(i)^2);
    else
        costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
    end
       
    wslRec11(i) = rho*g*D1^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
 
    i=i+1;
end     

i=1;
while (i<=250-5+1)
    
    theta(i) = pi/4.0;
    
    Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
    
    if beta(i)*D1<=(nu^3/dissip2)^0.25
        SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip2^0.5*log(2.0*beta(i)));
    else
        SF(i) = 5.0*wslRemin(i)^2*(beta(i)*D1)^(2.0/3.0)*Fic(i)/(8.0*nu*dissip2^(1.0/3.0)*log(2.0*beta(i)));
    end
        
    if SF(i)<=0.1
        costheta(i) = 1.0/3.0;
    elseif SF(i)>=5.0
        costheta(i) = 2.0/(15.0*SF(i)^2);
    else
        costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
    end
       
    wslRec12(i) = rho*g*D1^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
 
    i=i+1;
end 

i=1;
while (i<=250-5+1)
    
    theta(i) = pi/4.0;
    
    Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
    
    if beta(i)*D1<=(nu^3/dissip3)^0.25
        SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip3^0.5*log(2.0*beta(i)));
    else
        SF(i) = 5.0*wslRemin(i)^2*(beta(i)*D1)^(2.0/3.0)*Fic(i)/(8.0*nu*dissip3^(1.0/3.0)*log(2.0*beta(i)));
    end
        
    if SF(i)<=0.1
        costheta(i) = 1.0/3.0;
    elseif SF(i)>=5.0
        costheta(i) = 2.0/(15.0*SF(i)^2);
    else
        costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
    end
       
    wslRec13(i) = rho*g*D1^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
 
    i=i+1;
end 

i=1;
while (i<=250-5+1)
    Rel(i) = beta(i)*(log(2.0*beta(i))+log(4)-0.5)*D2^3*rho*g/(32*mu*nu);
    Relf(i) = 0.0;
    while (abs(Rel(i)-Relf(i))/Rel(i)>0.01)
        Relf(i) = Rel(i); 
        Fv(i) = int(f,x,Rel(i),inf) + log(Rel(i))-(exp(-Rel(i))-1)/Rel(i) + eulerc -0.5 -log(4);
        Fh(i) = 0.5*((int(f,x,2*Rel(i),inf)+log(2*Rel(i))-exp(-2*Rel(i))+eulerc+1)/(2.0*Rel(i))+int(f,x,2*Rel(i),inf)+log(Rel(i))+eulerc-3*log(2)+1);
        Mv(i) = log(2.0*beta(i))-Fv(i);
        Mh(i) = 2.0*log(2.0*beta(i))-2.0*Fh(i);
        wslRemin(i) = rho*g*D2^2/(16*mu)*Mv(i); %incomning flow velocity as U in COX
        Rel(i) = wslRemin(i)*beta(i)*D2/(2.0*nu);
    end
    i=i+1;
end    

i=1;
while (i<=250-5+1)
    
    theta(i) = pi/4.0;
    
    Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
    
    if beta(i)*D2<=(nu^3/dissip1)^0.25
        SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip1^0.5*log(2.0*beta(i)));
    else
        SF(i) = 5.0*wslRemin(i)^2*(beta(i)*D2)^(2.0/3.0)*Fic(i)/(8.0*nu*dissip1^(1.0/3.0)*log(2.0*beta(i)));
    end
        
    if SF(i)<=0.1
        costheta(i) = 1.0/3.0;
    elseif SF(i)>=5.0
        costheta(i) = 2.0/(15.0*SF(i)^2);
    else
        costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
    end
       
    wslRec21(i) = rho*g*D2^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
 
    i=i+1;
end     

i=1;
while (i<=250-5+1)
    
    theta(i) = pi/4.0;
    
    Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
    
    if beta(i)*D2<=(nu^3/dissip2)^0.25
        SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip2^0.5*log(2.0*beta(i)));
    else
        SF(i) = 5.0*wslRemin(i)^2*(beta(i)*D2)^(2.0/3.0)*Fic(i)/(8.0*nu*dissip2^(1.0/3.0)*log(2.0*beta(i)));
    end
        
    if SF(i)<=0.1
        costheta(i) = 1.0/3.0;
    elseif SF(i)>=5.0
        costheta(i) = 2.0/(15.0*SF(i)^2);
    else
        costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
    end
       
    wslRec22(i) = rho*g*D2^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
 
    i=i+1;
end 

i=1;
while (i<=250-5+1)
    
    theta(i) = pi/4.0;
    
    Fic(i)=-12.0*sin(2.0*theta(i))/(5.0*Rel(i)^2)*(0.5/(1-cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1-cos(theta(i))))-2)/((1-cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1-cos(theta(i))),inf)-log((Rel(i)*(1-cos(theta(i)))))-eulerc)+0.5/(1+cos(theta(i)))*(2+(2.0*exp(-Rel(i)*(1+cos(theta(i))))-2)/((1+cos(theta(i)))*Rel(i))-int(f,x,Rel(i)*(1+cos(theta(i))),inf)-log((Rel(i)*(1+cos(theta(i)))))-eulerc)-1.0/cos(theta(i))/(1-cos(theta(i)))*(1-(1-exp(-Rel(i)*(1-cos(theta(i)))))/(Rel(i)*(1-cos(theta(i)))))+1.0/cos(theta(i))/(1+cos(theta(i)))*(1-(1-exp(-Rel(i)*(1+cos(theta(i)))))/(Rel(i)*(1+cos(theta(i))))));
    
    if beta(i)*D2<=(nu^3/dissip3)^0.25
        SF(i) = 5.0*wslRemin(i)^2*Fic(i)/(8.0*nu^0.5*dissip3^0.5*log(2.0*beta(i)));
    else
        SF(i) = 5.0*wslRemin(i)^2*(beta(i)*D2)^(2.0/3.0)*Fic(i)/(8.0*nu*dissip3^(1.0/3.0)*log(2.0*beta(i)));
    end
        
    if SF(i)<=0.1
        costheta(i) = 1.0/3.0;
    elseif SF(i)>=5.0
        costheta(i) = 2.0/(15.0*SF(i)^2);
    else
        costheta(i) = 0.0753*SF(i)^(-0.6692)-0.0188;
    end
       
    wslRec23(i) = rho*g*D2^2/(16*mu)*(Mv(i)+costheta(i)*(Mh(i)-Mv(i)));        
 
    i=i+1;
end 

figure(1)
plot(beta,wslRec11,'-r','linewidth',3)
hold on
plot(beta,wslRec12,'--b','linewidth',3)
hold on
plot(beta,wslRec13,':g','linewidth',3)
xlim([5 250])
ylim([0 0.002])
set(gca,'FontSize',24);
title(' Settling Velocity','fontsize',18)
xlabel('aspect ratio $\beta $','fontsize',36,'Interpreter','latex') 
ylabel('${w_s (m/s)}$','fontsize',36,'Interpreter','latex')
legend({'$D=2 \mu m$ surrounded by turbulence with high dissipation rate','$D=2 \mu m$ surrounded by turbulence with intermediate dissipation rate','$D=2 \mu m$ surrounded by turbulence with low dissipation rate'},'Location','northwest','fontsize',24,'Interpreter','latex')

figure(2)
plot(beta,wslRec21,':m','linewidth',3)
hold on
plot(beta,wslRec22,'--c','linewidth',3)
hold on
plot(beta,wslRec23,'-.y','linewidth',3)
xlim([5 250])
ylim([0 0.1])
set(gca,'FontSize',24);
title(' Settling Velocity','fontsize',18)
xlabel('aspect ratio $\beta $','fontsize',36,'Interpreter','latex') 
ylabel('${w_s (m/s)}$','fontsize',36,'Interpreter','latex')
legend({'$D=15 \mu m$ surrounded by turbulence with high dissipation rate','$D=15 \mu m$ surrounded by turbulence with intermediate dissipation rate','$D=15 \mu m$ surrounded by turbulence with low dissipation rate'},'Location','northwest','fontsize',24,'Interpreter','latex')