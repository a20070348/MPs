clear;
clc;

nu = 1.5* 10^(-5);
mu = 1.8* 10^(-5);
eulerc = 0.5772156649;
dissip = 0.0001;
g = 9.8;

% the shape parameters for flat fiber
width=0.2* 10^(-6); 
height = 0.11* 10^(-6);
length = 9* 10^(-6);
rho = 3400.0;

B=0.5*width;
C=0.5*height; % B and C are the parameters for flat fiber

Defl=2.0*sqrt(B*C);
betafl = length/Defl;

% no touch for the following if no background on the theory
f = @(x) exp(-x)./x;

fun = @(y) sqrt(1+C^2*y.^2./(B^2*(B^2-y.^2)));
Rel = betafl*(log(2.0*betafl)+log(4)-0.5-log(0.5*pi*(B+C)/integral(fun,0,B))-0.5*(B-C)/(B+C))*Defl^3*rho*g/(32*mu*nu);
Relf = 0.0;

while (abs(Rel-Relf)/Rel > 0.01)
    Relf = Rel; 
    Fv = integral(f,Rel,inf) + log(Rel)-(exp(-Rel)-1)/Rel + eulerc -0.5 -log(4)+log(0.5*pi*(B+C)/integral(fun,0,B))+0.5*(B-C)/(B+C);
    Fh = 0.5*((integral(f,2*Rel,inf)+log(2*Rel)-exp(-2*Rel)+eulerc+1)/(2.0*Rel)+integral(f,2*Rel,inf)+log(Rel)+eulerc-3*log(2)+1+6*log(0.25*pi*(B+C)/integral(fun,0,B))+2*(B-C)/(B+C));
    Mv = log(2.0*betafl)-Fv;
    Mh = 2.0*log(2.0*betafl)-2.0*Fh;
    wslRemin = rho*g*Defl^2/(16*mu)*Mv; %incoming flow velocity as U in COX
    Rel = wslRemin*betafl*Defl/(2.0*nu);
end 

theta = pi/4.0;

Fic_fun = @(x) -12.0*sin(2.0*theta)/(5.0*Rel^2)*(0.5/(1-cos(theta))*(2+(2.0*exp(-Rel*(1-cos(theta)))-2)/((1-cos(theta))*Rel)-integral(f,Rel*(1-cos(theta)),inf)-log(Rel*(1-cos(theta)))-eulerc)+0.5/(1+cos(theta))*(2+(2.0*exp(-Rel*(1+cos(theta)))-2)/((1+cos(theta))*Rel)-integral(f,Rel*(1+cos(theta)),inf)-log(Rel*(1+cos(theta)))-eulerc)-1.0/cos(theta)/(1-cos(theta))*(1-(1-exp(-Rel*(1-cos(theta))))/(Rel*(1-cos(theta))))+1.0/cos(theta)/(1+cos(theta))*(1-(1-exp(-Rel*(1+cos(theta))))/(Rel*(1+cos(theta)))));

Fic = Fic_fun(0);

if betafl*Defl<=(nu^3/dissip)^0.25
    SF = 5.0*wslRemin^2*Fic/(8.0*nu^0.5*dissip^0.5*log(2.0*betafl));
else
    SF = 5.0*wslRemin^2*(betafl*Defl)^(2.0/3.0)*Fic/(8.0*nu*dissip^(1.0/3.0)*log(2.0*betafl));
end

if SF<=0.1
    costheta = 1.0/3.0;
elseif SF>=5.0
    costheta = 2.0/(15.0*SF^2);
else
    costheta = 0.0753*SF^(-0.6692)-0.0188;
end

% the settling velocity and aerodynamic diameter for flat fiber
wslRecfl = rho*g*Defl^2/(16*mu)*(Mv+costheta*(Mh-Mv));
Dzenderfl=sqrt(wslRecfl*18*nu/rho/g)*10^6;