figure
str = streamslice(Zr,Rr,omega_ce.*cos(thetaB),omega_ce.*sin(thetaB),3); hold on;
set(str,'Color','black');
plot3(0, 0, 1e30, 'g--','LineWidth',2);hold on
surf(Zr, Rr, omega_ce); view(2); colorbar; shading flat;
pbaspect([2, 1, 1])
% title('Electron Cyclotron Frequency $\omega_{ce}/\omega$','FontSize',20,'Interpreter','latex')
xlabel('z (cm)','FontSize',20,'Interpreter','latex')
ylabel('r (cm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);

figure
surf(Zr, Rr, nu); view(2);  colorbar; shading flat;hold on;
pbaspect([2, 1, 1])
title('Collision Frequency $\nu/\omega$','FontSize',20,'Interpreter','latex')
xlabel('z (cm)','FontSize',20,'Interpreter','latex')
ylabel('r (cm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);

figure
surf(Zr, Rr, omega_pe); view(2);  colorbar; shading flat;hold on;%set(gca,'colorscale','log');
pbaspect([2, 1, 1])
% title('Plasma Frequency $\omega_{pe}/\omega$','FontSize',20,'Interpreter','latex')
xlabel('z (cm)','FontSize',20,'Interpreter','latex')
ylabel('r (cm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);

figure
field = sqrt(abs(Ez).^2 + abs(Ey).^2 + abs(Ex).^2); 
% streamslice(Zm*0.001 + 0.05,Rm*0.001,Bz,Br,2); hold on;
% plot3([-3,-3,57]*0.001+0.05, [0,25, 25]*0.001, [1e30,1e30, 1e30], 'g--','LineWidth',2);hold on
% plot3([57, 80, z(end)]*0.001+0.05, [25, r(end), r(end)]*0.001,[1e40,1e40, 1e40],'m--','LineWidth',2); hold on;
surf(Zr, Rr, field); view(2);  colorbar; shading flat;hold on;  set(gca,'colorscale','log');
% [C,h] = contour(Zr,Rr,F, [1 2 3 4 5 6 7 8],'k'); clabel(C,h); clabel(C,h, 'FontSize',15,'Color','red');
pbaspect([2, 1, 1])
% title('Electric field Norm $[V/m]$','FontSize',20,'Interpreter','latex')
xlabel('z (cm)','FontSize',20,'Interpreter','latex')
ylabel('r (cm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20);box on; ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);

figure
P(Rr>1.5)=0;
surf(Zr, Rr, P); view(2);  colorbar; set(gca,'colorscale','log'); shading flat;hold on;
pbaspect([2, 1, 1])
xlabel('z (cm)','FontSize',20,'Interpreter','latex')
ylabel('r (cm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20);box on; ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);

figure
colormap('jet');set(gca,'CLim',[1 8])
contourf(Zr,Rr,F);c = colorbar; c.Ticks = [1 2 3 4 5 6 7 8];hold on;
% plot3([-3,-3,57]*0.001+0.05, [0,25, 25]*0.001, [1e30,1e30, 1e30], 'g--','LineWidth',2);hold on
% plot3([57, 80, z(end)]*0.001+0.05, [25, r(end), r(end)]*0.001,[1e40,1e40, 1e40],'m--','LineWidth',2);
pbaspect([2, 1, 1])
% title('CMA Propagation regions','FontSize',20,'Interpreter','latex')
xlabel('z (cm)','FontSize',20,'Interpreter','latex')
ylabel('r (cm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20);box on; ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);

figure
surf(Zr, Rr, abs(Ey)); view(2);  colorbar; shading flat;hold on;
pbaspect([2, 1, 1])
xlabel('z (cm)','FontSize',20,'Interpreter','latex')
ylabel('r (cm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);

figure
surf(Zr, Rr, angle(Ez)); view(2);  colorbar; shading flat;hold on; colormap hsv;caxis([-pi pi]);
pbaspect([2, 1, 1])
xlabel('z (cm)','FontSize',20,'Interpreter','latex')
ylabel('r (cm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20);box on;ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);