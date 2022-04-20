clear;
clc;

L = 1.0* 10^(-4):0.25* 10^(-4):25.0 * 10^(-4);
abundance = 10^(1.1*3.4)*(L.*10^6).^(-1.1);

CSF =0.00:0.005:0.65;
Freq = 0.06/(sqrt(2*pi*0.03^2))*exp(-(CSF-0.08).^2/(2*0.03^2))+0.94/(sqrt(2*pi*0.19^2))*exp(-(CSF-0.44).^2/(2*0.19^2));

nu = 1.5* 10^(-5);mu = 1.8* 10^(-5);eulerc = 0.5772156649;dissip = 0.001;
rho = 1000.0;g = 9.8;
beta =1./(CSF.^2);
D = L'.*CSF.^2;
Dim = size(D);

Cfreq0 = abundance'.*Freq;
Cfreq = Cfreq0(:);

syms x
f = exp(-x)/x;
Rel = 0.0*D;Relf = 0.0*D;
wslRe = 0.0*D;wslRec = 0.0*D;wslRemin = 0.0*D;
Fv = 0.0*D;Fh = 0.0*D;Mv = 0.0*D;Mh = 0.0*D;
Fic = 0.0*D;
SF = 0.0*D;costheta = 0.0*D+0.5;theta = 0.0*D;thetaf = 0.0*D;
i=1;j=1;
while (i<=Dim(1,:))
    while (j<=Dim(:,2))
        Rel(i,j) = beta(j)*(log(2.0*beta(j))+log(4)-0.5)*D(i,j)^3*rho*g/(32*mu*nu);
        Relf(i,j) = 0.0;
        while (abs(Rel(i,j)-Relf(i,j))/Rel(i,j)>0.01)
            Relf(i,j) = Rel(i,j); 
            Fv(i,j) = int(f,x,Rel(i,j),inf) + log(Rel(i,j))-(exp(-Rel(i,j))-1)/Rel(i,j) + eulerc -0.5 -log(4);
            Fh(i,j) = 0.5*((int(f,x,2*Rel(i,j),inf)+log(2*Rel(i,j))-exp(-2*Rel(i,j))+eulerc+1)/(2.0*Rel(i,j))+int(f,x,2*Rel(i,j),inf)+log(Rel(i,j))+eulerc-3*log(2)+1);
            Mv(i,j) = log(2.0*beta(j))-Fv(i,j);
            Mh(i,j) = 2.0*log(2.0*beta(j))-2.0*Fh(i,j);
            wslRemin(i,j) = rho*g*D(i,j)^2/(16*mu)*Mv(i,j); %incomning flow velocity as U in COX
            Rel(i,j) = wslRemin(i,j)*beta(j)*D(i,j)/(2.0*nu);
        end
        j=j+1;
    end
    i=i+1;
end

i=1;j=1;
while (i<=Dim(1,:))
    while (j<=Dim(:,2))
    
        theta(i,j) = pi/4.0;
    
        Fic(i,j)=-12.0*sin(2.0*theta(i,j))/(5.0*Rel(i,j)^2)*(0.5/(1-cos(theta(i,j)))*(2+(2.0*exp(-Rel(i,j)*(1-cos(theta(i,j))))-2)/((1-cos(theta(i,j)))*Rel(i,j))-int(f,x,Rel(i,j)*(1-cos(theta(i,j))),inf)-log((Rel(i,j)*(1-cos(theta(i,j)))))-eulerc)+0.5/(1+cos(theta(i,j)))*(2+(2.0*exp(-Rel(i,j)*(1+cos(theta(i,j))))-2)/((1+cos(theta(i,j)))*Rel(i,j))-int(f,x,Rel(i,j)*(1+cos(theta(i,j))),inf)-log((Rel(i,j)*(1+cos(theta(i,j)))))-eulerc)-1.0/cos(theta(i,j))/(1-cos(theta(i,j)))*(1-(1-exp(-Rel(i,j)*(1-cos(theta(i,j)))))/(Rel(i,j)*(1-cos(theta(i,j)))))+1.0/cos(theta(i,j))/(1+cos(theta(i,j)))*(1-(1-exp(-Rel(i,j)*(1+cos(theta(i,j)))))/(Rel(i,j)*(1+cos(theta(i,j))))));
    
        if beta(j)*D(i,j)<=(nu^3/dissip)^0.25
            SF(i,j) = 5.0*wslRemin(i,j)^2*Fic(i,j)/(8.0*nu^0.5*dissip^0.5*log(2.0*beta(j)));
        else
            SF(i,j) = 5.0*wslRemin(i,j)^2*(beta(j)*D(i,j))^(2.0/3.0)*Fic(i,j)/(8.0*nu*dissip^(1.0/3.0)*log(2.0*beta(j)));
        end
        
        if SF(i,j)<=0.1
            costheta(i,j) = 1.0/3.0;
        elseif SF(i)>=5.0
            costheta(i,j) = 2.0/(15.0*SF(i,j)^2);
        else
            costheta(i,j) = 0.0753*SF(i,j)^(-0.6692)-0.0188;
        end
       
        wslRec(i,j) = rho*g*D(i,j)^2/(16*mu)*(Mv(i,j)+costheta(i,j)*(Mh(i,j)-Mv(i,j))); 
        
        j=j+1;
    end
    i=i+1;
end  

wslRec0 = wslRec(:);
Result0 = [wslRec0,Cfreq];
Result = sortrows(Result0);

% plot(D,wslRe,'--green','linewidth',3)
% hold on
plot(Result(:,1),Result(:,2),'-black','linewidth',3)
set(gca,'FontSize',24);
title(' Fiber Settling Velocity Frequency','fontsize',18)
xlabel('${w_s (m/s)}$','fontsize',36,'Interpreter','latex') 
ylabel('Frequency','fontsize',36,'Interpreter','latex')
% legend({'high dissipation rate','low dissipation rate'},'Location','northwest','fontsize',18)