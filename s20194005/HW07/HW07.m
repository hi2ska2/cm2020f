q = 1.602192e-19; 
k = 1.380662e-23;
ni = 1e16;
T=300;
VT=(k*T)/q;
N1=zeros(9,1);
N2=zeros(9,1);
for i=1:9
     N1(i,1)=1*10^(15+i);
     N2(i,1)=-1*10^(15+i);
end

for i=1:9
    phi = 0.5;
    for newton = 1:10000
        Jaco = ni*(-1/VT)*exp(-phi/VT)-ni*(1/VT)*exp(phi/VT);
        res = N1(i,1)+ni*exp(-phi/VT)-ni*exp(phi/VT); 
        update= Jaco\(-res);
        phi=phi+update;
    end
    phi_numerical(i,1)=phi;
end

for i=1:9
    phi = -0.5;
    for newton = 1:10000
        Jaco = ni*(-1/VT)*exp(-phi/VT)-ni*(1/VT)*exp(phi/VT);
        res = N2(i,1)+ni*exp(-phi/VT)-ni*exp(phi/VT);
        update= Jaco\(-res);
        phi=phi+update;
    end
    phi_numerical2(i,1)=phi;
end

for i = 1:9
    phi_analytic(i,1) = VT*asinh(N1(i,1)/(2*ni));
    phi_analytic2(i,1) = VT*asinh(N2(i,1)/(2*ni));
    error(i,1) = (phi_analytic(i,1)-phi_numerical(i,1))/phi_analytic(i,1)*100;
    error2(i,1) = (phi_analytic2(i,1)-phi_numerical2(i,1))/phi_analytic2(i,1)*100;
end
