run('config.m')
cte = constants_and_units.constants;

R = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/R').*100;
Z = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/Z').*100;
omega_pe = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/omega_pe');
omega_ce = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/omega_ce'); 

N = (2*pi*13.56e6*omega_pe.r./cte.qe).^2*cte.me*cte.eps0;
B = 2*pi*13.56e6*omega_ce.r*cte.me/cte.qe;

z = Z(1,:);
l = 0*z;

for i = 1:length(z)
    
    [l(i),~,~] = lambda(N(1,i), B(1,i));
    
end

plot(z,l);xlim([20,67])

i = find(flip(l)>=0.24); i = i(1); i = length(z)-i;

nref  = N(1,i);
Bref  = B(1,i);
[lref,wceref,wperef] = lambda(nref,Bref);

[l,wce,wpe] = lambda(5e17,100e-4);

function [l,wce,wpe] = lambda(n,Ba)

    cte = constants_and_units.constants;

    w  = 2*pi * 13.56e6;
%     n  = 5e17;
%     Ba = 0.01;

    wpe = sqrt(n * cte.qe^2 / cte.eps0 / cte.me);
    wce = cte.qe * Ba / cte.me;

    de = cte.c0 / wpe;
    l  = 2*pi * de * sqrt(wce/w);

end