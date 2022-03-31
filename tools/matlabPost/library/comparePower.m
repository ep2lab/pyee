getfields
cd(['P_', study])

P1 = load('P1.mat'); %m = 1
P2 = load('P2.mat'); %m = -1
P3 = load('P3.mat'); %m = 3
P4 = load('P4.mat'); %m = -3

color = [1e-2 1e6];
% color = [1e-2 1e4];
colormapT = 'jet';

figure()

subplot('Position',[0.1 0.52 0.4 0.4])
pcolor(Zr, Rr, P1.P);  shading flat;hold on;
contour(Zr, Rr, bound, 1, 'r', 'LineWidth',1)
ylabel('\begin{tabular}{c} $m=1$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
set(gca, 'xtick',[],'ytick',[0.2, 1 ,1.8] ,'colorscale','log', 'FontSize',24);
c = colorbar('northoutside');
c.Label.String = '$[W/m^3]$'; c.Label.Interpreter = 'LaTex';
caxis(color)
colormap(colormapT)

subplot('Position',[0.53 0.52 0.4 0.4])
pcolor(Zr, Rr, P2.P); shading flat;hold on;
contour(Zr, Rr, bound, 1, 'r', 'LineWidth',1)
ylabel('$m=-1$',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle')
box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
set(gca, 'YAxisLocation','right')
set(gca, 'xtick',[],'ytick',[] ,'colorscale','log', 'FontSize',24);
c = colorbar('northoutside');
c.Label.String = '$[W/m^3]$'; c.Label.Interpreter = 'LaTex';
caxis(color)
colormap(colormapT)

subplot('Position',[0.1 0.1 0.4 0.4])
pcolor(Zr, Rr, P3.P); shading flat;hold on;
contour(Zr, Rr, bound, 1, 'r', 'LineWidth',1)
xlabel('z (cm)','FontSize',24,'Interpreter','latex')
ylabel('\begin{tabular}{c} $m=3$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
set(gca, 'ytick',[0.2, 1 ,1.8], 'colorscale','log','FontSize',24);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
caxis(color)
colormap(colormapT)

subplot('Position',[0.53 0.1 0.4 0.4])
pcolor(Zr, Rr, P4.P); shading flat;hold on;
contour(Zr, Rr, bound,'r')
set(gca, 'YAxisLocation','right')
xlabel('z (cm)','FontSize',24,'Interpreter','latex')
ylabel('$m=-3$',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle')
set(gca, 'ytick',[],'colorscale','log','FontSize',24);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
caxis(color)
colormap(colormapT)

cd(MAIN_PATH)
