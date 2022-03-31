c = constants_and_units.constants;
fig = figure;
ne = (f.omega_pe * 2*pi*13.56*1e6).^2 / c.qe.^2 * c.me * c.eps0;
ne(~f.PointFlag) = 0;

subplot(2,1,1)
[s,ax,c1] = plotField(f.Zr, f.Rr, ne,  'title', '$n_e (1/m^3)$','scale','log');
caxis([1e14,3e19])
delete(c1)
rectangle('Position',[10 0 10 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
rectangle('Position',[28.5 3 7 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
set(gca, 'FontSize',18)
set(ax.XLabel,'visible','off')
spcFigs
contour(f.Zr,f.Rr,f.PointFlag,'LineWidth',3,'LineColor','red','LineStyle','--');

subplot(2,1,2)
[s,ax,c2] = plotField(f.Zr, f.Rr, ne,  'title', '$n_e (1/m^3)$','scale','log');
set(ax.Title,'visible','off')
caxis([1e14,3e19])
delete(c2);
rectangle('Position',[28.5 3 7 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
set(gca, 'FontSize',18)
xlim([20 40]);ylim([0 3.8])
spcFigs
line([20 20 32.5],[0 1.18 1.18],'LineWidth',3,'Color','red');
line([32.5 40 40],[1.18 3.8 0],'LineWidth',3,'Color','red','LineStyle','--');

h = axes(fig,'visible','off');
set(h,'colorscale','log'); 
set(gca,'FontSize',18);
co = colorbar; 
co.Ticks = 10.^linspace(14,19,6);
caxis([1e14,3e19])

