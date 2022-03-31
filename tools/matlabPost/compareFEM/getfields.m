R = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/R');
Z = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/Z');

omega_ce = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/omega_ce'); omega_ce = omega_ce.r;
omega_pe = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/omega_pe'); omega_pe = omega_pe.r;
nu = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/nu'); nu = nu.r;
thetaB = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/thetaB'); thetaB = thetaB.r;

P = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/P');
F = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/F');

jz = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/jz'); jz = jz.r + 1i*jz.i;
jy = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/jy'); jy = jy.r + 1i*jy.i;
jx = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/jx'); jx = jx.r + 1i*jx.i;

Ez = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/Ez'); Ez = (Ez.r + 1i*Ez.i);
Ey = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/Ey'); Ey = (Ey.r + 1i*Ey.i);
Ex = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/Ex'); Ex = (Ex.r + 1i*Ex.i);

Bz = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/Bz'); Bz = 1e4*(Bz.r + 1i*Bz.i);
By = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/By'); By = 1e4*(By.r + 1i*By.i);
Bx = h5read([WORKFOLDER,'cases\',study,'\results.h5'],'/Bx'); Bx = 1e4*(Bx.r + 1i*Bx.i);