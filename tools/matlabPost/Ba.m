c = constants_and_units.constants;
fig = figure;
B = f.omega_ce/c.qe*c.me*2*pi*13.56*1e6*10000;
levels = power(10,-2:0.5:2.5)/c.qe*c.me*2*pi*13.56*1e6*10000;

fact = 0.81*(1 - abs(32-f.Zr)./100).^0.16;
thetaB = atan2(sin(f.thetaB).*fact, cos(f.thetaB));

subplot(2,1,1)
str = streamslice(f.Zr,f.Rr,f.omega_ce.*cos(thetaB),f.omega_ce.*sin(thetaB),3); hold on;
set(str,'Color','black');
plot3(0, 0, 1e30, 'g--','LineWidth',2);hold on
contour(f.Zr, f.Rr, f.omega_ce/c.qe*c.me*2*pi*13.56*1e6*10000,power(10,-1.0:0.5:2.0)/c.qe*c.me*2*pi*13.56*1e6*10000,'k');
hold on
s = surface(f.Zr, f.Rr, f.Rr*0-10,  f.omega_ce/c.qe*c.me*2*pi*13.56*1e6*10000);
s.LineStyle = "none";
set(gca,'colorscale','log');
caxis([1e-2,2e3]);
pbaspect([2, 1, 1]);
title('$B_a (G)$','FontSize',18,'Interpreter','latex')
% xlabel('z (cm)','FontSize',18,'Interpreter','latex')
ylabel('r (cm)','FontSize',18,'Interpreter','latex')
set(gca,'FontSize',18);
box on;ylim([min(f.Rr(:)) max(f.Rr(:))]); xlim([min(f.Zr(:)) max(f.Zr(:))]);
contour(f.Zr, f.Rr, double(abs(f.jz)>max(abs(f.jz(:)))*0.95),'LineWidth',5,'LineColor','magenta');
l = line([32.5 40 40],[1.18 3.75 0]); set(l, 'LineWidth',3,'Color','red','LineStyle','--')
line([20 20 32.5],[0 1.18 1.18],'LineWidth',3,'Color','red');
contour(f.Zr, f.Rr, f.omega_ce, [1,1],'LineWidth',3,'LineColor','green');
% % contour(f.Zr,f.Rr,f.PointFlag,'LineWidth',3,'LineColor','red','LineStyle','--');
lineobj = streamline(f.Zr,f.Rr,f.omega_ce.*cos(thetaB),f.omega_ce.*sin(thetaB),40,3.75,[.5,10000]);
lineobj.LineWidth = 3;
lineobj.LineStyle = "--";
lineobj.Color     = "red";
rectangle('Position',[10 0 10 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
rectangle('Position',[28.5 3 7 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
caxis([1e-2,2e3])

subplot(2,1,2)
contour(f.Zr, f.Rr, f.omega_ce/c.qe*c.me*2*pi*13.56*1e6*10000,power(10,-1.0:0.5:2.0)/c.qe*c.me*2*pi*13.56*1e6*10000,'k');
hold on
s = surface(f.Zr, f.Rr, f.Rr*0-10,  f.omega_ce/c.qe*c.me*2*pi*13.56*1e6*10000);
s.LineStyle = "none";
set(gca,'colorscale','log');
caxis([1e-2,2e3]);
xlabel('z (cm)','FontSize',18,'Interpreter','latex')
ylabel('r (cm)','FontSize',18,'Interpreter','latex')
set(gca,'FontSize',18);
line([20 20 32.5],[0 1.18 1.18],'LineWidth',3,'Color','red');
line([32.5 40 40],[1.18 3.8 0],'LineWidth',3,'Color','red','LineStyle','--');
str = streamslice(f.Zr,f.Rr,f.omega_ce.*cos(thetaB),f.omega_ce.*sin(thetaB),15,'noarrows'); 
set(str,'Color','black');
contour(f.Zr, f.Rr, double(abs(f.jz)>max(abs(f.jz(:)))*0.95),'LineWidth',18,'LineColor','magenta');
xlim([20 40]);ylim([0 3.8])
rectangle('Position',[28.5 3 7 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
pbaspect([2 1 1])
box on

h = axes(fig,'visible','off');
set(h,'colorscale','log'); 
set(gca,'FontSize',18);
co = colorbar; 
co.Ticks = 10.^(-2:1:3);
caxis([1e-2,2e3])

