clear;
clc;

nu = 1.5* 10^(-5);mu = 1.8* 10^(-5);eulerc = 0.5772156649;dissip = 0.0001;
rho = 1000.0;g = 9.8;

geodata1=load('/Users/shuolinxiao/Downloads/PostDoc/Cornell_Qi/projects/MPF/code/Beaver_Mountain_Fiber.dat');
geodata2=load('/Users/shuolinxiao/Downloads/PostDoc/Cornell_Qi/projects/MPF/code/CO98_Fiber.dat');
geodata3=load('/Users/shuolinxiao/Downloads/PostDoc/Cornell_Qi/projects/MPF/code/UT95_Fiber.dat');
geodata4=load('/Users/shuolinxiao/Downloads/PostDoc/Cornell_Qi/projects/MPF/code/Beaver_Mountain_Fiber_Wet.dat');
% geodata4=load('/Users/shuolinxiao/Downloads/PostDoc/Cornell_Qi/projects/MPF/code/sample.dat');

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
    if (wslRecfl(i)<=0.0)
        wslRecfl(i) = 0.0;
    end
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
    if (wslReccy(i)<=0.0)
        wslReccy(i) = 0.0;
    end
    i=i+1;
end 


Dzenderfl=sqrt(wslRecfl*18*nu/rho/g)*10^6;
Dzendercy=sqrt(wslReccy*18*nu/rho/g)*10^6;
Dzender=[Dzenderfl';Dzendercy']';

Weightfl = zeros([1 6]);Weightcy = zeros([1 6]);

i=1;
while (i<=Dim(:,1))
    if Dzenderfl(i)<=0.5
        Weightfl(1) = Weightfl(1)+1.0/6.0*pi*rho*Dzenderfl(i)^3*10^(-18);
    elseif Dzenderfl(i)<=5
        Weightfl(2) = Weightfl(2)+1.0/6.0*pi*rho*Dzenderfl(i)^3*10^(-18);
    elseif Dzenderfl(i)<=10
        Weightfl(3) = Weightfl(3)+1.0/6.0*pi*rho*Dzenderfl(i)^3*10^(-18);
    elseif Dzenderfl(i)<=20
        Weightfl(4) = Weightfl(4)+1.0/6.0*pi*rho*Dzenderfl(i)^3*10^(-18);
    elseif Dzenderfl(i)<=50
        Weightfl(5) = Weightfl(5)+1.0/6.0*pi*rho*Dzenderfl(i)^3*10^(-18);
    elseif Dzenderfl(i)<=100
        Weightfl(6) = Weightfl(6)+1.0/6.0*pi*rho*Dzenderfl(i)^3*10^(-18);
    end
    i=i+1;
end 
Weightfl=Weightfl/(sum(Weightfl));

i=1; 
while (i<=Dim(:,1))
    if Dzendercy(i)<=0.5
        Weightcy(1) = Weightcy(1)+1.0/6.0*pi*rho*Dzendercy(i)^3*10^(-18);
    elseif Dzendercy(i)<=5
        Weightcy(2) = Weightcy(2)+1.0/6.0*pi*rho*Dzendercy(i)^3*10^(-18);
    elseif Dzendercy(i)<=10
        Weightcy(3) = Weightcy(3)+1.0/6.0*pi*rho*Dzendercy(i)^3*10^(-18);
    elseif Dzendercy(i)<=20
        Weightcy(4) = Weightcy(4)+1.0/6.0*pi*rho*Dzendercy(i)^3*10^(-18);
    elseif Dzendercy(i)<=50
        Weightcy(5) = Weightcy(5)+1.0/6.0*pi*rho*Dzendercy(i)^3*10^(-18);
    elseif Dzendercy(i)<=100
        Weightcy(6) = Weightcy(6)+1.0/6.0*pi*rho*Dzendercy(i)^3*10^(-18);  
    end
    i=i+1;
end 

Weightcy=Weightcy/(sum(Weightcy));



