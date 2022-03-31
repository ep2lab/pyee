run('config.m')
f = getfields(PYEE_PATH, STUDY, FILENAME);


Ez = interp2(f.Z.Ez,f.R.Ez,f.Ez,f.Z.Qa,f.R.Qa,'spline');
Ex = interp2(f.Z.Ex,f.R.Ex,f.Ex,f.Z.Qa,f.R.Qa,'spline');


Epar = Ez .* cos(f.thetaB)  + Ex .* sin(f.thetaB);
Eper = -Ez .* sin(f.thetaB) + Ex .* cos(f.thetaB);

figure

subplot('Position',[0.1 0.63 0.4 0.28])
[s,ax,c] = plotField(f.Z.Qa, f.R.Qa, Epar,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c); caxis([1e-2 1e4])
ylabel('\begin{tabular}{c} $\mid\mathbf{E_\parallel}\mid$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
ax.XTick = [];
ax.YTick = [5 10 15];
c = colorbar('northoutside');
c.Label.String = '$[V/m]$'; c.Label.Interpreter = 'LaTex'; 
spcFigs

subplot('Position',[0.53 0.63 0.4 0.28])
[s,ax,c] = plotField(f.Z.Qa, f.R.Qa, Epar,'type','phase');delete(ax.YLabel);delete(ax.XLabel);delete(c);
ax.YTick = [];
ax.XTick = [];
set(ax, 'YAxisLocation','right')
c = colorbar('northoutside');
c.Label.String = '$[rad]$'; c.Label.Interpreter = 'LaTex'; 
spcFigs

subplot('Position',[0.1 0.37 0.4 0.22])
[s,ax,c] = plotField(f.Z.Qa, f.R.Qa, Eper,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c); caxis([1e-2 1e4])
ylabel('\begin{tabular}{c} $\mid\mathbf{E_\perp}\mid$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
ax.XTick = [];
spcFigs

subplot('Position',[0.53 0.37 0.4 0.22])
[s,ax,c] = plotField(f.Z.Qa, f.R.Qa, Eper,'type','phase');delete(ax.YLabel);delete(ax.XLabel);delete(c);
set(gca, 'YAxisLocation','right')
ax.YTick = [];
ax.XTick = [];
spcFigs

subplot('Position',[0.1 0.11 0.4 0.22])
[s,ax,c] = plotField(f.Z.Ey, f.R.Ey, f.Ey,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c); caxis([1e-2 1e4])
xlabel('z (cm)','Interpreter','latex')
ylabel('\begin{tabular}{c} $\mid\mathbf{E_\theta}\mid$ \\ r(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
spcFigs

subplot('Position',[0.53 0.11 0.4 0.22])
[s,ax,c] = plotField(f.Z.Ey, f.R.Ey, f.Ey,'type','phase');delete(ax.YLabel);delete(ax.XLabel);delete(c);
xlabel('z (cm)','Interpreter','latex')
set(ax, 'YAxisLocation','right');
ax.YTick = [];
spcFigs