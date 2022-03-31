clear;clc;close all;

%% Load FD fields and mesh
<<<<<<< HEAD
study = 'vacuum';
% WORKFOLDER = 'C:\Users\Pedro Jimenez\Documents\GitHub\pyee\';
WORKFOLDER = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\pyee\';

getfields
Ez = abs(Ez)./max(abs(Ez(:)));
Ey = abs(Ey)./max(abs(Ey(:)));
Ex = abs(Ex)./max(abs(Ex(:)));
nP = P./max(P(:));

%% Interpolate and load FEM fields
EzFEM = getFEM([study, '/Ez'], Z, R); EzFEM = EzFEM./max(EzFEM(:));
EyFEM = getFEM([study, '/Ey'], Z, R); EyFEM = EyFEM./max(EyFEM(:));
ExFEM = getFEM([study, '/Ex'], Z, R); ExFEM = ExFEM./max(ExFEM(:));
PFEM = getFEM([study, '/P'], Z, R); PFEM = PFEM./max(PFEM(:));

%% Plots
Zr = Z*100;
Rr = R*100;

colorE = [0 1];
colormapT = 'parula';
figure('units','normalized','outerposition',[0 0 1 1])

subplot('Position',[0.1 0.63 0.4 0.3])
pcolor(Zr, Rr, Ez);  shading flat;hold on;
ylabel('\begin{tabular}{c} $\mid\mathbf{E_z}\mid$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
set(gca,'FontSize',18);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
set(gca, 'xtick',[],'ytick',[20, 50, 80]);
c = colorbar('northoutside');
c.Label.String = 'Normalized Field FD'; c.Label.Interpreter = 'LaTex'; c.FontSize=20;
caxis(colorE)
colormap(colormapT)

subplot('Position',[0.53 0.63 0.4 0.3])
pcolor(Zr, Rr, EzFEM); view(2); shading flat;hold on;
set(gca,'FontSize',18);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
set(gca,'xtick',[], 'ytick',[])
set(gca, 'YAxisLocation','right')
c = colorbar('northoutside');
c.Label.String = 'Normalized Field FEM'; c.Label.Interpreter = 'LaTex'; c.FontSize=20;
caxis(colorE)

subplot('Position',[0.1 0.37 0.4 0.22])
pcolor(Zr, Rr, Ex); view(2); shading flat;hold on;
ylabel('\begin{tabular}{c} $\mid\mathbf{E_r}\mid$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
set(gca,'FontSize',18);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
set(gca,'xtick',[],'ytick',[20, 50, 80])
caxis(colorE)
colormap(colormapT)

subplot('Position',[0.53 0.37 0.4 0.22])
pcolor(Zr, Rr, ExFEM); shading flat;hold on;
set(gca, 'YAxisLocation','right')
set(gca,'FontSize',18);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
set(gca,'xtick',[], 'ytick',[])
caxis(colorE)
colormap(colormapT)

subplot('Position',[0.1 0.1 0.4 0.22])
pcolor(Zr, Rr, Ey); shading flat;hold on;
xlabel('z (cm)','Interpreter','latex')
ylabel('\begin{tabular}{c} $\mid\mathbf{E_\theta}\mid$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
set(gca,'FontSize',18);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
set(gca,'ytick',[20, 50, 80])
caxis(colorE)
colormap(colormapT)

subplot('Position',[0.53 0.1 0.4 0.22])
pcolor(Zr, Rr, EyFEM); shading flat;hold on;
xlabel('z (cm)','Interpreter','latex')
set(gca, 'YAxisLocation','right')
set(gca,'FontSize',18);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
set(gca,'ytick',[])
caxis(colorE)
colormap(colormapT)

figure('units','normalized','outerposition',[0 0 1 1])
normFD = sqrt(Ex.^2 + Ez.^2 + Ey.^2); normFD = normFD./max(normFD(:));
normFEM = sqrt(ExFEM.^2 + EzFEM.^2 + EyFEM.^2); normFEM = normFEM./max(normFEM(:));
pcolor(Zr,Rr,abs(normFD - normFEM)*100); shading interp; 
set(gca,'FontSize',24)
ylabel('r(cm)','Interpreter','latex');
xlabel('z(cm)','Interpreter','latex');
c=colorbar;c.Label.String = 'Relative Error %'; c.Label.Interpreter='Latex';c.FontSize=22;
c.Location = 'northoutside';
colormap('jet')
daspect([1 2 1])
caxis([0,3])
=======
study = 'Helicon';
% WORKFOLDER = 'C:\Users\Pedro Jimenez\Documents\GitHub\pyee\';
WORKFOLDER = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\pyee\';
getfields

%% Interpolate and load FEM fields
EzFEM = getFEM('Ez', Z, R);
EyFEM = getFEM('Ey', Z, R);
ExFEM = getFEM('Ex', Z, R);
P = getFEM('P', Z, R);

surf(Z,R,EzFEM);shading interp;view(2);colorbar;
>>>>>>> parent of 85236e0... 	new file:   cases/rwave/config_sim.py
