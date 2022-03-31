close all; clear
PYEE_PATH = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\pyee\';
FIGS_PATH = 'C:\Users\pedro\OneDrive\Escritorio\figs\';
Pinput = 350;

fig = figure;

study = 'SPC2020'; FILE_NAME = '\results.h5'; 
subplot(2,1,1)
getfields; 
P(P<1e-10) = 1e-10; 
LOGZ  = Zr(1,:)>=20 & Zr(1,:)<=44; LOGR = Rr(:,1)<=3; 
PointFlag = PointFlag(LOGR,LOGZ);
P1 = P(LOGR,LOGZ); Z1 = Zr(LOGR,LOGZ); R1 = Rr(LOGR,LOGZ); j1 = jz(LOGR,LOGZ);
plotField(Z1, R1, P1, PointFlag*0, 'title', '$Q_a [W/m^{3}]$','scale','log')
caxis([1e2,1e8]);c = colorbar;c.Ticks = [0.1,10,1000,1e5,1e7,1e8];
patch([20, 32, 37, 20], [1.25, 1.25, 3, 3],'white')
cont = contour(Z1, R1, double(abs(j1)>max(abs(j1(:)))*0.99),'LineWidth',10,'LineColor','magenta');
xlabel('');set(gca, 'xtick',[])

study = 'SPC2020'; FILE_NAME = '\results_small.h5';
subplot(2,1,2)
getfields; P(P<1e-10) = 1e-10; Zr = Zr + 20;
plotField(Zr, Rr, P, PointFlag*0, 'title', '','scale','log')
caxis([1e2,1e8]);c = colorbar;c.Ticks = [0.1,10,1000,1e5,1e7,1e8];
patch([20, 32, 37, 20], [1.25, 1.25, 3, 3],'white')
cont = contour(Z1, R1, double(abs(j1)>max(abs(j1(:)))*0.99),'LineWidth',10,'LineColor','magenta');

save_fig(fig,'Pcomp', FIGS_PATH, true)


figure;
study = 'SPC2020'; FILE_NAME = '\results.h5'; 
getfields; E = sqrt(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2);
E(E<1e-10) = 1e-10; 
LOGZ  = Zr(1,:)>=20 & Zr(1,:)<=44; LOGR = Rr(:,1)<=3; 
PointFlag = PointFlag(LOGR,LOGZ);
E1 = E(LOGR,LOGZ); Z1 = Zr(LOGR,LOGZ); R1 = Rr(LOGR,LOGZ); j1 = jz(LOGR,LOGZ);
plotField(Z1, R1, E1, PointFlag*0, 'title', '$\mid\mathbf{E}\mid [V/m]$','scale','log')
caxis([1e4,1e8]);c = colorbar;c.Ticks = [1e4,1e6,1e8];
patch([20, 32, 37, 20], [1.25, 1.25, 3, 3],'white')
cont = contour(Z1, R1, double(abs(j1)>max(abs(j1(:)))*0.99),'LineWidth',10,'LineColor','magenta');