run('config.m')
f = getfields(PYEE_PATH, STUDY, FILENAME);

colorE = [1e-2 1e4]; ticksE = 10.^(-2:2:4);
colorB = [1e-6 1];   ticksB = 10.^(-6:2:1);
% colormapT = 'jet';

fig = figure
% set(fig,'units','normalized','outerposition',[0 0 1 1],'Menu','none','ToolBar','none');

subplot('Position',[0.1 0.63 0.4 0.28])
[s,ax,c] = plotField(f.Z.Ez, f.R.Ez, f.Ez,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c);
ylabel('\begin{tabular}{c} $\mid\mathbf{E_z}\mid$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
ax.XTick = [];
ax.YTick = [5 10 15];
c = colorbar('northoutside');
c.Label.String = '$[V/m]$'; c.Label.Interpreter = 'LaTex'; c.Ticks = ticksE;
caxis(colorE)
spcFigs

subplot('Position',[0.53 0.63 0.4 0.28])
[s,ax,c] = plotField(f.Z.Ey, f.R.Ey, f.Bz,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c);
ylabel('$\mid\mathbf{B_z}\mid$',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle')
ax.YTick = [];
ax.XTick = [];
set(ax, 'YAxisLocation','right')
c = colorbar('northoutside');
c.Label.String = '$[G]$'; c.Label.Interpreter = 'LaTex'; c.Ticks = ticksB;
caxis(colorB)
spcFigs

subplot('Position',[0.1 0.37 0.4 0.22])
[s,ax,c] = plotField(f.Z.Ex, f.R.Ex, f.Ex,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c)
ylabel('\begin{tabular}{c} $\mid\mathbf{E_r}\mid$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
ax.XTick = [];
caxis(colorE)
spcFigs

subplot('Position',[0.53 0.37 0.4 0.22])
[s,ax,c] = plotField(f.Z.Ey, f.R.Ey, f.Bx,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c);
set(gca, 'YAxisLocation','right')
ylabel('$\mid\mathbf{B_r}\mid$',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle')
ax.YTick = [];
ax.XTick = [];
caxis(colorB)
spcFigs

subplot('Position',[0.1 0.11 0.4 0.22])
[s,ax,c] = plotField(f.Z.Ey, f.R.Ey, f.Ey,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c);
xlabel('z (cm)','Interpreter','latex')
ylabel('\begin{tabular}{c} $\mid\mathbf{E_\theta}\mid$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
caxis(colorE)
spcFigs

subplot('Position',[0.53 0.11 0.4 0.22])
[s,ax,c] = plotField(f.Z.Ey, f.R.Ey, f.By,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c);
ylabel('$\mid\mathbf{B_\theta}\mid$',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle')
xlabel('z (cm)','Interpreter','latex')
set(ax, 'YAxisLocation','right');
ax.YTick = [];
caxis(colorB)
spcFigs

figure;plotField(f.Z.Ey,f.R.Ey,f.Ey,'type','phase','title','[rad]');spcFigs;