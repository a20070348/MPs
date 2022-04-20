nu = 1.5* 10^(-5);mu = 1.8* 10^(-5);eulerc = 0.5772156649;dissip = 0.001;
b=1.0*10^(-5);l=5.0*10^(-5):1.0*10^(-6):2.5*10^(-3);beta=l/b;

Re1 = 10.0;Re2 = 0.0;
while (abs(Re1-Re2)>0.0001)
   Re1 = Re2;
   Fb = integral(exp(-x)/x,Re1,Inf) + log(Re1)-(exp(-Re1)-1)/Re1 + eulerc -0.5 -log(4);
   Fp = 0.5*((integral(exp(-x)/x,2*Re1,Inf)+log(2*Re1)-exp(-2*Re1)-eulerc+1)/2/Re1+integral(exp(-x)/x,2*Re1,Inf)+log(Re1)+eulerc-3*log(2)+1);
   Mb = log(beta)/(8*pi*mu*l)*(1-Fb/log(beta));
   Mp = log(beta)/(4*pi*mu*l)*(1-Fp/log(beta));
    
   
end