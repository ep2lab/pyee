run('config.m')
f1 = getfields(PYEE_PATH, 'paper/param_overA3_clean', FILENAME);
f = getfields(PYEE_PATH, 'paper/param_overA3', FILENAME);
f2 = getfields(PYEE_PATH, 'paper/param_vacA3', FILENAME);

%%
colorRange = [1e-1 1e3];
figure
subplot('Position',[0.1 0.63 0.4 0.28])
[s,ax,c] = plotField(f.Z.Ey, f.R.Ey, f.Ey,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c); caxis(colorRange)
ylabel('\begin{tabular}{c} $\mathrm{R}$. \\ $r$(cm)  \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
ax.XTick = [];
ax.YTick = [5 10 15];
c = colorbar('northoutside');
c.Label.String = '(V/m)'; c.Label.Interpreter = 'LaTex'; 
c.Ticks = [1e-1,1e0,1e1,1e2,1e3];
rectangle('Position',[10 0 10 4], 'FaceColor', 'black')
rectangle('Position',[28.5 3 7 4], 'FaceColor', 'black')
pbaspect([2.5, 1, 1])
spcFigs

subplot('Position',[0.53 0.63 0.4 0.28])
[s,ax,c] = plotField(f.Z.Ey, f.R.Ey, f.Ey,'type','phase');delete(ax.YLabel);delete(ax.XLabel);delete(c);
ax.YTick = [];
ax.XTick = [];
set(ax, 'YAxisLocation','right')
c = colorbar('northoutside');
c.Label.String = '(deg)'; c.Label.Interpreter = 'LaTex'; c.Ticks = [-180,-90,0,90,180];
rectangle('Position',[10 0 10 4], 'FaceColor', 'black')
rectangle('Position',[28.5 3 7 4], 'FaceColor', 'black')
pbaspect([2.5, 1, 1])

subplot('Position',[0.1 0.37 0.4 0.22])
[s,ax,c] = plotField(f1.Z.Ey, f1.R.Ey, f1.Ey,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c); caxis(colorRange)
ylabel('\begin{tabular}{c} $\mathrm{T}.$ \\ $r$(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
ax.XTick = [];
pbaspect([2.5, 1, 1])
spcFigs

subplot('Position',[0.53 0.37 0.4 0.22])
[s,ax,c] = plotField(f1.Z.Ey, f1.R.Ey, f1.Ey,'type','phase');delete(ax.YLabel);delete(ax.XLabel);delete(c);
set(gca, 'YAxisLocation','right')
ax.YTick = [];
ax.XTick = [];
pbaspect([2.5, 1, 1])

subplot('Position',[0.1 0.11 0.4 0.22])
[s,ax,c] = plotField(f2.Z.Ey, f2.R.Ey, f2.Ey,'scale','log');delete(ax.YLabel);delete(ax.XLabel);delete(c); caxis(colorRange)
xlabel('$z$ (cm)','Interpreter','latex')
ylabel('\begin{tabular}{c} $\mathrm{V}.$ \\ $r$(cm) \end{tabular}',...
    'Interpreter','latex','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
rectangle('Position',[10 0 10 4], 'FaceColor', 'black')
rectangle('Position',[28.5 3 7 4], 'FaceColor', 'black')
pbaspect([2.5, 1, 1])
spcFigs

subplot('Position',[0.53 0.11 0.4 0.22])
[s,ax,c] = plotField(f2.Z.Ey, f2.R.Ey, f2.Ey,'type','phase');delete(ax.YLabel);delete(ax.XLabel);delete(c);
xlabel('$z$ (cm)','Interpreter','latex')
set(ax, 'YAxisLocation','right');
ax.YTick = [];
rectangle('Position',[10 0 10 4], 'FaceColor', 'black')
rectangle('Position',[28.5 3 7 4], 'FaceColor', 'black')
pbaspect([2.5, 1, 1])