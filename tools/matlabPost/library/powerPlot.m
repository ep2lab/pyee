%Scaled to give 350 W in the domain%

P0 = P; 
P0(~PointFlag) = nan;

P1 = P0;
P1(isnan(P0))  = 0;

dep = 2*pi*trapz(Zr(1,:)*0.01, trapz(Rr(:,1)*0.01, P.*Rr(:,1)*0.01));

fac = 350/dep;

P0 = P0*fac;

Zr0 = Zr - 5;
figure
pcolor(Zr0, Rr, P0); view(2);  colorbar; set(gca,'colorscale','log'); shading flat;hold on;
contour(Zr0, Rr, bound, 1, 'k', 'LineWidth',1)
pbaspect([2, 1, 1])
xlabel('z [cm]','FontSize',30,'Interpreter','latex')
ylabel('r [cm]','FontSize',30,'Interpreter','latex')
set(gca,'FontSize',30);box on; ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr0(:)) max(Zr0(:))]);
title('$Q_a [W/m^3]$','FontSize',30,'Interpreter','latex')
pbaspect([1.5 1 1]); set(gca, 'YTick', [0 0.5 1 1.5 2 2.5], 'XTick', [0 10 20])
grid off
caxis([1e4 1e8])

figure
pcolor(Zr0, Rr, abs(Ez)*sqrt(fac)); view(2);  colorbar; set(gca,'colorscale','linear'); shading flat;hold on;
contour(Zr0, Rr, bound, 1, 'k', 'LineWidth', 1)
pbaspect([2, 1, 1])
xlabel('z [cm]','FontSize',30,'Interpreter','latex')
ylabel('r [cm]','FontSize',30,'Interpreter','latex')
set(gca,'FontSize',30);box on; ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr0(:)) max(Zr0(:))]);
title('$|\hat{E}_z| [V/m]$','FontSize',30,'Interpreter','latex')
pbaspect([1.5 1 1]); set(gca, 'YTick', [0 0.5 1 1.5 2 2.5], 'XTick', [0 10 20])
% caxis([1e-2 1e5])
grid off

figure
pcolor(Zr0, Rr, abs(Ex)*sqrt(fac)); view(2);  colorbar; set(gca,'colorscale','log'); shading flat;hold on;
contour(Zr0, Rr, bound, 1, 'k', 'LineWidth', 1)
pbaspect([2, 1, 1])
xlabel('z [cm]','FontSize',30,'Interpreter','latex')
ylabel('r [cm]','FontSize',30,'Interpreter','latex')
set(gca,'FontSize',30);box on; ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr0(:)) max(Zr0(:))]);
title('$|\hat{E}_r| [V/m]$','FontSize',30,'Interpreter','latex')
pbaspect([1.5 1 1]); set(gca, 'YTick', [0 0.5 1 1.5 2 2.5], 'XTick', [0 10 20])
% caxis([1e-2 1e5])
grid off

figure
pcolor(Zr0, Rr, abs(Ey)*sqrt(fac)); view(2);  colorbar; set(gca,'colorscale','log'); shading flat;hold on;
contour(Zr0, Rr, bound, 1, 'k', 'LineWidth', 1)
pbaspect([2, 1, 1])
xlabel('z [cm]','FontSize',30,'Interpreter','latex')
ylabel('r [cm]','FontSize',30,'Interpreter','latex')
set(gca,'FontSize',30);box on; ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr0(:)) max(Zr0(:))]);
title('$|\hat{E}_\theta| [V/m]$','FontSize',30,'Interpreter','latex')
pbaspect([1.5 1 1]); set(gca, 'YTick', [0 0.5 1 1.5 2 2.5], 'XTick', [0 10 20])
% caxis([1e-2 1e5])
grid off
