tag = 'm1pos_h05neg_Bneg';

%% Plot fields
allFields
figHandle = gcf;
set(figHandle, 'units','normalized','outerposition',[0 0 1 1])
saveas(figHandle, ['local_data/polarity_HPT05M/', tag, '_fields.png'])

%% Plot power
% P(Rr>0.999) = nan; 
figure
colormap('jet')
pcolor(Zr, Rr, P); view(2);  colorbar; set(gca,'colorscale','log'); shading flat;hold on;
xlabel('z [cm]','FontSize',30,'Interpreter','latex')
ylabel('r [cm]','FontSize',30,'Interpreter','latex')
set(gca,'FontSize',30);box on; ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
title('$Q_a [W/m^3]$','FontSize',30,'Interpreter','latex')
pbaspect([2 1 1]); set(gca, 'YTick', [0 0.5 1 1.5 2 2.5], 'XTick', [0 10 20])
grid off
caxis([1e-1 1e4])
figHandle = gcf;
set(figHandle, 'units','normalized','outerposition',[0 0 1 1])
saveas(figHandle, ['local_data/polarity_HPT05M/', tag, '_power.png'])

close all
