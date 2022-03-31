norm150 = load('norm150.mat');
P150 = load('P150.mat');
norm300 = load('norm300.mat'); 
P300 = load('P300.mat'); 

color = [1e-2 1e4];
colorE = [0 500];
colormapE = 'parula';
colormapT = 'jet';
bound = zeros(size(Zr));
bound(~nu==0) = 1;

figure()

subplot('Position',[0.13 0.52 0.4 0.4])
pcolor(Zr, Rr, norm150.field);  shading flat;hold on;
contour(Zr, Rr, bound, 1, 'r', 'LineWidth',1)
ylabel('\begin{tabular}{c} $\mathbf{B}_{0,max} = 150\; G$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
set(gca, 'xtick',[],'ytick',[0.2, 1 ,1.8] , 'FontSize',18);
c = colorbar('northoutside');
c.Label.String = 'Electric field norm $\mid\mathbf{E}\mid$ $[V/m]$'; c.Label.Interpreter = 'LaTex'; c.FontSize=20;
caxis(colorE)
colormap(colormapE)

subplot('Position',[0.58 0.52 0.4 0.4])
pcolor(Zr, Rr, P150.P); shading flat;hold on;
contour(Zr, Rr, bound, 1, 'r', 'LineWidth',1)
box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
set(gca, 'YAxisLocation','right')
set(gca, 'xtick',[],'ytick',[] ,'colorscale','log', 'FontSize',18);
c = colorbar('northoutside');
c.Label.String = 'Power Deposition $[W/m^3]$'; c.Label.Interpreter = 'LaTex'; c.FontSize=20;
caxis(color)
colormap(colormapT)

subplot('Position',[0.13 0.1 0.4 0.4])
pcolor(Zr, Rr, norm300.field); shading flat;hold on;
contour(Zr, Rr, bound, 1, 'r', 'LineWidth',1)
xlabel('z (cm)','FontSize',24,'Interpreter','latex')
ylabel('\begin{tabular}{c} $\mathbf{B}_{0,max} = 300\; G$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
set(gca, 'ytick',[0.2, 1 ,1.8], 'FontSize',18);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
caxis(colorE)
colormap(colormapE)

subplot('Position',[0.58 0.1 0.4 0.4])
pcolor(Zr, Rr, P300.P); shading flat;hold on;
contour(Zr, Rr, bound,'r')
set(gca, 'YAxisLocation','right')
xlabel('z (cm)','FontSize',24,'Interpreter','latex')
set(gca, 'ytick',[],'colorscale','log','FontSize',18);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
caxis(color)
colormap(colormapT)

