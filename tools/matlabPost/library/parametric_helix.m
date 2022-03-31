%% Paths
% close all; clc; clear;
MAIN_PATH = pwd;

study = 'HPT05M';

PYEE_PATH = 'C:\Users\Pedro Jimenez\Documents\GitHub\pyee\';
% PYEE_PATH = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\pyee\';

addpath('C:\Users\Pedro Jimenez\OneDrive - Universidad Carlos III de Madrid\Codes\magnetic\HPT03 FEMM\hpt03');
addpath('C:\Users\Pedro Jimenez\OneDrive - Universidad Carlos III de Madrid\Codes\utilities-master')
addpath('C:\Users\Pedro Jimenez\OneDrive - Universidad Carlos III de Madrid\Codes\constants_and_units-master')
addpath('C:\Users\Pedro Jimenez\OneDrive - Universidad Carlos III de Madrid\Codes\magnetic_field-master')
% 
% addpath('C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\utilities-master')
% addpath('C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\constants_and_units-master')
% addpath('C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\magnetic_field-master')

%% Run Pyee
i = 0;
hvec = -0.8:0.05:-0.6;
for h = hvec
    
    disp([num2str(100*i/length(hvec)),'%'])
    i = i + 1;
    
    cd(PYEE_PATH)
    flag = dos(['activate base & python run_parametric.py ', num2str(h)]);
    
    if flag == 1
        break
    end
    %% Postprocess
    cd(MAIN_PATH)
    getfields
    bound = PointFlag;
    LogScale = false;

    if h==0
        tag = ['m+1_+B_h', num2str(h*10)];
    elseif h>0
        tag = ['m+1_+B_h+', num2str(abs(h)*10)];
    else
        tag = ['m+1_+B_h-', num2str(abs(h)*10)];
    end

    %% Plot power
    P1 = P;
    P1(isnan(P1))  = 0;
    dep(i) = 2*pi*trapz(Zr(1,:)*0.01, trapz(Rr(:,1)*0.01, P1.*Rr(:,1)*0.01));
        
    figure
    set(gcf,'visible','off')
    colormap('jet')
    pcolor(Zr, Rr, P); view(2);  colorbar; set(gca,'colorscale','log'); shading flat;hold on;
    xlabel('z [cm]','FontSize',30,'Interpreter','latex')
    ylabel('r [cm]','FontSize',30,'Interpreter','latex')
    set(gca,'FontSize',30);box on; ylim([min(Rr(:)) max(Rr(:))]); xlim([min(Zr(:)) max(Zr(:))]);
    title(['$Q_a [W/m^3]$, $helix = $', num2str(h)],'FontSize',30,'Interpreter','latex')
    pbaspect([2 1 1]); set(gca, 'YTick', [0 0.5 1 1.5 2 2.5], 'XTick', [0 10 20])
    grid off
    caxis([1e-1 1e4])
    figHandle = gcf;
    set(figHandle, 'units','normalized','outerposition',[0 0 1 1])
    saveas(figHandle, ['local_data/parametric_helix/', tag, '_power.png'])
    close all
end

figure
semilogy(hvec,dep,'*-','LineWidth',2)
set(gca,'FontSize',30);box on;
xlabel('Helix number','FontSize',30,'Interpreter','latex')
ylabel('Power Deposition','FontSize',30,'Interpreter','latex')

%% Make gif
filename = 'powerHPT05.gif';
for n = 1:length(hvec)
        
    if hvec(n)==0
        image_name = ['local_data/parametric_helix/m+1_+B_h', num2str(hvec(n)*10),'_power','.png'];
    elseif hvec(n)>0
        image_name = ['local_data/parametric_helix/m+1_+B_h+', num2str(abs(hvec(n))*10),'_power','.png'];
    else
        image_name = ['local_data/parametric_helix/m+1_+B_h-', num2str(abs(hvec(n))*10),'_power','.png'];
    end
    
    im = imread(image_name);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end