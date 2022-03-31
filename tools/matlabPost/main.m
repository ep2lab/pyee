run('config.m')
%% Run Pyee
% cd(PYEE_PATH)
% 
% if computeFlag
%     errorFlag = dos(['activate base & python run_pyee.py ',study], '-echo');
%     if errorFlag
%         error('Python Computation Failed')
%     end
% end
% 
% cd(MAIN_PATH)

%% Postprocess
f = getfields(PYEE_PATH, STUDY, FILENAME);
c = constants_and_units.constants;

%% Plots
CASE_PATH = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\pyee\cases\05MextP';
SIM_PATH = 'C:\Users\pedro\OneDrive\Documentos\LocalSimData\HPT05MC0';
ss_path = [SIM_PATH,'\CORE\out\SimState.hdf5'];
mf_array = Coil_HPT05M(ss_path);
[~,Bz0,Br0] = mf_array.field_2d(f.Zr/100,f.Rr/100);

%%
fact = 0.90*(1 - abs(32-f.Zr)./100).^0.08;
Br = Br0 .* fact;
Bz = Bz0 ./ fact;

figure
% str = streamslice(f.Zr,f.Rr,f.omega_ce.*cos(f.thetaB),f.omega_ce.*sin(f.thetaB),3); hold on;
% plot3(0, 0, 1e30, 'g--','LineWidth',2);hold on
contourf(f.Zr, f.Rr, f.omega_ce/c.qe*c.me*2*pi*13.56*1e6*10000,power(10,-2:0.5:2.5)/c.qe*c.me*2*pi*13.56*1e6*10000);hold on
% contourf(f.Zr, f.Rr,B*10000,power(10,-2:0.5:2)/c.qe*c.me*2*pi*13.56*1e6*10000); hold on
str = streamslice(f.Zr,f.Rr,Bz,Br,3); 
[~,i] = find(f.Zr(1,:)>40);
contour(f.Zr(:,i),f.Rr(:,i),f.PointFlag(:,i),'LineWidth',3,'LineColor','red','LineStyle','--');
lineobj = streamline(f.Zr,f.Rr,Bz,Br,40,3.7,[.5,10000]);
lineobj.LineWidth = 2;
lineobj.LineStyle = "--";
lineobj.Color     = "blue";
set(str,'Color','black');
view(2); co=colorbar; shading flat;set(gca,'colorscale','log') 
co.Ticks = 10.^(-2:1:3);
pbaspect([2, 1, 1]);caxis([1e-2,2e3])
title('$B [G]$','FontSize',28,'Interpreter','latex')
xlabel('z (cm)','FontSize',28,'Interpreter','latex')
ylabel('r (cm)','FontSize',28,'Interpreter','latex')
set(gca,'FontSize',28);box on;ylim([min(f.Rr(:)) max(f.Rr(:))]); xlim([min(f.Zr(:)) max(f.Zr(:))]);

%%
figure
spcFigs;
nu = (f.nu * 2*pi*13.56*1e6);
plotField(f.Zr, f.Rr, nu,  'title', '$\nu (Hz)$','scale','log')
contour(f.Zr,f.Rr,f.PointFlag,'LineWidth',3,'LineColor','red');
l = line([40,40],[0,3.8]); set(l, 'LineWidth',4,'Color','red','LineStyle','--')
rectangle('Position',[10 0 10 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
rectangle('Position',[28.5 3 7 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
c = colorbar; caxis([1e6,1e10]);
set(gca, 'FontSize',30)

figure
spcFigs; c = constants_and_units.constants;
ne = (f.omega_pe * 2*pi*13.56*1e6).^2 / c.qe.^2 * c.me * c.eps0;
ne(~f.PointFlag) = 0;
plotField(f.Zr, f.Rr, ne,  'title', '$n_e [1/m^3]$','scale','log')
contour(f.Zr,f.Rr,f.PointFlag,'LineWidth',3,'LineColor','red');
l = line([40,40],[0,3.8]); set(l, 'LineWidth',3,'Color','red','LineStyle','--')
rectangle('Position',[10 0 10 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
rectangle('Position',[28.5 3 7 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
c = colorbar; c.Ticks = 10.^linspace(14,19,6);caxis([1e14,1.5e19])
set(gca, 'FontSize',35)

figure
flag = f.PointFlag*0 + 1;
flag(or(f.Z.Qa<20 | f.Z.Qa>32.5, f.R.Qa>1.25)) = 0;
flag(and(f.Z.Qa<40, f.R.Qa <= (1.25 + (3.8-1.25)*(f.Z.Qa - 32.5)/7.5)))=1;
flag(~flag)=NaN;
rectangle('Position',[10 0 10 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
rectangle('Position',[28.5 3 7 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
c = constants_and_units.constants;
ne = (f.omega_pe * 2*pi*13.56*1e6).^2 / c.qe.^2 * c.me * c.eps0;
plotField(f.Zr, f.Rr, ne,'scale','log')
line([20 20 32.5],[0 1.25 1.25],'LineWidth',3,'Color','red');
line([32.5 40 40],[1.25 3.8 0],'LineWidth',3,'Color','red','LineStyle','--');
str = streamslice(f.Zr,f.Rr,f.omega_ce.*cos(f.thetaB),f.omega_ce.*sin(f.thetaB),15,'noarrows'); 
set(str,'Color','black');
contour(f.Zr, f.Rr, double(abs(f.jz)>max(abs(f.jz(:)))*0.95),'LineWidth',25,'LineColor','magenta');
c = colorbar;
caxis([1e18,1.5e19]);c.Ticks = [0.1,0.2,0.5,1]*10.^19;
c.TickLabels = {'10^{18}','2\cdot10^{18}','5\cdot10^{18}','10^{19}'};
xlim([20 40]);ylim([0 3.8])
pbaspect([2 1 1]);set(gca,'FontSize',35)
box on

figure
flag = f.PointFlag*0 + 1;
flag(or(f.Z.Qa<20 | f.Z.Qa>32.5, f.R.Qa>1.25)) = 0;
flag(and(f.Z.Qa<40, f.R.Qa <= (1.25 + (3.8-1.25)*(f.Z.Qa - 32.5)/7.5)))=1;
flag(~flag)=NaN;
hold on;
c = constants_and_units.constants;
B = (f.omega_ce * 2*pi*13.56*1e6)*c.me / c.qe * 10000;
contourf(f.Zr, f.Rr, B)
c = colorbar; caxis([0 1500]);set(gca,'colorscale','log') 
c.Ticks = [0.2,0.3,0.5,0.7,1,1.5]*1000;
c.TickLabels = {'200','300','500','700','1000','1500'};
set(gca,'FontSize',35)
% title('$B [G]$','FontSize',28,'Interpreter','latex')
xlabel('z (cm)','FontSize',35,'Interpreter','latex')
ylabel('r (cm)','FontSize',35,'Interpreter','latex')
line([20 20 32.5],[0 1.25 1.25],'LineWidth',3,'Color','red');
line([32.5 40 40],[1.25 3.8 0],'LineWidth',3,'Color','red','LineStyle','--');
str = streamslice(f.Zr,f.Rr,f.omega_ce.*cos(f.thetaB),f.omega_ce.*sin(f.thetaB),15,'noarrows'); 
set(str,'Color','black');
contour(f.Zr, f.Rr, double(abs(f.jz)>max(abs(f.jz(:)))*0.95),'LineWidth',25,'LineColor','magenta');
xlim([20 40]);ylim([0 3.8])
pbaspect([2 1 1])
box on
% rectangle('Position',[10 0 10 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
rectangle('Position',[28.5 3 7 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)


% save_fig(gcf,'oce', FIGS_PATH, saveFlag)
% spcFigs
% 
% figure
% spcFigs
% plotField(f.Zr, f.Rr, f.nu*2*pi*13.56*1e6,'title', '$\nu_{e} [Hz]$','pointflag',f.PointFlag)
% caxis([1e5,1e9])
% c = colorbar; c.Ticks = 10.^linspace(5,9,5); 
% spcFigs
% figure
% plotField(f.Zr, f.Rr, f.F,  'title', 'Propagation Regions')
% save_fig(gcf,'oce', FIGS_PATH, saveFlag)

% figure
% spcFigs
% normE = sqrt(abs(f.Ex).^2 + abs(f.Ey).^2 + abs(f.Ez).^2);
% plotField(f.Zr, f.Rr, normE,  'title', '$\mid\mathbf{E}\mid [V/m]$', 'scale','log')
% caxis([1e2,1e8])
% save_fig(gcf,'oce', FIGS_PATH, saveFlag)

% figure
% mesh(f.Z.Ez, f.R.Ez,angle(f.Ez));hold on;
% colormap('hsv')
% caxis([-pi pi]);view(2);colorbar
% spcFigs;shading flat

% clims = [1e-2 1e6];
% clims = [0 300];
% alt = 1;
% % 
% figure
% mesh(f.Z.Ez, f.R.Ez, abs(f.Ez)*alt, abs(f.Ez));hold on; 
% view(2);set(gca,'colorscale','log');colorbar
% contour3(f.Zr,f.Rr,f.PointFlag*1e6,[1 1]*1e6,'LineWidth',1,'LineColor','magenta');
% colormap('jet'); caxis(clims); shading flat
% 
% figure
% mesh(f.Z.Ex, f.R.Ex, abs(f.Ex)*alt, abs(f.Ex));hold on; 
% view(2);set(gca,'colorscale','log');colorbar
% contour3(f.Zr,f.Rr,f.PointFlag*1e20,[1 1]*1e20,'LineWidth',1,'LineColor','magenta');
% colormap('jet'); caxis(clims); shading flat
% 
% figure
% mesh(f.Z.Ey, f.R.Ey, abs(f.Ey)*alt, abs(f.Ey));hold on; 
% view(2);set(gca,'colorscale','log');colorbar
% contour3(f.Zr,f.Rr,f.PointFlag*1e20,[1 1]*1e20,'LineWidth',1,'LineColor','magenta');
% colormap('jet'); caxis(clims); shading flat

% figure
% plotField(f.Zr, f.Rr, f.Qa, 'scale', 'log','title','[$W/m^3 ]$');
% caxis([1e-4 1e4])
% spcFigs;
% 
% figure
% plotField(f.Z.Ey, f.R.Ey, f.Ey,'type','phase','title','$[rad]$');
% spcFigs;
% 
% res = 2*pi*trapz(f.Zr(1,:)*0.01, trapz(f.Rr(:,1)*0.01, f.Qa.*f.Rr(:,1)*0.01));
% disp(sum(f.nu(:)))
% figure
% mesh(f.Z.y, f.R.y,angle(f.Ey));hold on;
% colormap('hsv')
% caxis([-pi pi]);view(2);colorbar
