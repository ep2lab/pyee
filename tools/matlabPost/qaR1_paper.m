run("config.m")
f  = getfields(PYEE_PATH,'paper\param_overA3',getname('\results',7));
f1 = getfields(PYEE_PATH,'paper\param_overA3_clean',getname('\results',7));
f2 = getfields(PYEE_PATH,'paper\param_vacA3',getname('\results',7));
f3 = getfields(PYEE_PATH,'paper\overA3_small',getname('\results',2));

%%
fig = figure;
subplot(4,1,1)
[s,ax,c] = plotField(f.Zr, f.Rr, f.Qa, 'title', '$Q_a^1 (W/m^{3})$','scale','log','pointflag',f.PointFlag & f.Zr<40 & f.Rr<1.25);
ylabel("R. r(cm)", 'FontSize', 16)
set(ax.XLabel,'visible','off')
set(ax,'FontSize',16)
delete(c) 
delete(ax.XLabel); ax.XTick = [];
caxis([1e-2 1e4])
contour(f.Zr,f.Rr,f.PointFlag,'LineWidth',3,'LineColor','red');
pbaspect([3, 1, 1])

subplot(4,1,2)
spcFigs
[s,ax,c] = plotField(f1.Zr, f1.Rr, f1.Qa,'scale','log','pointflag',f1.PointFlag & f1.Zr<40 & f1.Rr<1.25);
set(ax.Title,'visible','off')
set(ax,'FontSize',16)
delete(c)
delete(ax.XLabel); ax.XTick = [];
caxis([1e-2 1e4])
ylabel("T. r(cm)", 'FontSize', 16)
contour(f1.Zr,f1.Rr,f1.PointFlag,'LineWidth',3,'LineColor','red');
pbaspect([3, 1, 1])

subplot(4,1,3)
spcFigs
[s,ax,c] = plotField(f2.Zr, f2.Rr, f2.Qa,'scale','log','pointflag',f2.PointFlag & f2.Zr<40 & f2.Rr<1.25);
set(ax.Title,'visible','off')
delete(ax.XLabel); ax.XTick = [];
set(ax,'FontSize',16)
delete(c)
caxis([1e-2 1e4])
ylabel("V. r(cm)", 'FontSize', 16)
contour(f2.Zr,f2.Rr,f2.PointFlag,'LineWidth',3,'LineColor','red');
pbaspect([3, 1, 1])

subplot(4,1,4)
spcFigs
[s,ax,c] = plotField(f3.Zr+20, f3.Rr, f3.Qa,'scale','log','pointflag',f3.Zr<40 & f3.Rr<1.25);
set(ax.Title,'visible','off')
set(ax,'FontSize',16)
delete(c)
caxis([1e-2 1e4])
ylabel("S. r(cm)", 'FontSize', 16)
contour(f2.Zr+0.01,f2.Rr,f2.PointFlag,'LineWidth',3,'LineColor','red');
pbaspect([3, 1, 1])

h = axes(fig,'visible','off');
set(h,'colorscale','log'); 
c = colorbar(h); 
c.Ticks = 10.^(-2:2:4);
caxis([1e-2 1e4])
set(h,'FontSize',20)