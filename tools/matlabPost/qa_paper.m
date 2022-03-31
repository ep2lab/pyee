run("config.m")
f1 = getfields(PYEE_PATH,'paper\param_overA3_clean',getname('\results',7));
f = getfields(PYEE_PATH,'paper\param_overA3',getname('\results',7));
f2 = getfields(PYEE_PATH,'paper\param_vacA3',getname('\results',7));


%%
colorRange = [1e-3 1e4];
fig = figure;
subplot(4,1,1)
[s,ax,c] = plotField(f.Zr, f.Rr, f.Qa, 'title', '$Q_a^1 (W/m^{3})$','scale','log');
ylabel("R. r(cm)", 'FontSize', 18)
set(ax.XLabel,'visible','off')
set(ax,'FontSize',18)
delete(c) 
delete(ax.XLabel); ax.XTick = [];
caxis(colorRange)
rectangle('Position',[10 0 10 4], 'FaceColor', 'black')
rectangle('Position',[28.5 3 7 4], 'FaceColor', 'black')
pbaspect([3, 1, 1])

subplot(4,1,2)
spcFigs
[s,ax,c] = plotField(f1.Zr, f1.Rr, f1.Qa,'scale','log');
set(ax.Title,'visible','off')
set(ax,'FontSize',18)
delete(c)
delete(ax.XLabel); ax.XTick = [];
caxis(colorRange)
ylabel("T. r(cm)", 'FontSize', 18)
pbaspect([3, 1, 1])

subplot(4,1,3)
spcFigs
[s,ax,c] = plotField(f2.Zr, f2.Rr, f2.Qa,'scale','log');
set(ax.Title,'visible','off')
set(ax,'FontSize',18)
delete(c)
caxis(colorRange)
ylabel("V. r(cm)", 'FontSize', 18)
pbaspect([3, 1, 1])
rectangle('Position',[10 0 10 4], 'FaceColor', 'black')
rectangle('Position',[28.5 3 7 4], 'FaceColor', 'black')

subplot(4,1,4)
spcFigs
[s,ax,c] = plotField(f2.Zr, f2.Rr, f2.Qa,'scale','log');
set(ax.Title,'visible','off')
set(ax,'FontSize',18)
delete(c)
caxis(colorRange)
ylabel("V. r(cm)", 'FontSize', 18)
pbaspect([3, 1, 1])
rectangle('Position',[10 0 10 4], 'FaceColor', 'black')
rectangle('Position',[28.5 3 7 4], 'FaceColor', 'black')

h = axes(fig,'visible','off');
set(h,'colorscale','log'); 
c = colorbar(h); 
c.Ticks = 10.^(-2:2:4);
caxis(colorRange)
set(h,'FontSize',20)