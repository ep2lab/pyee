%% Paths
% close all; clc; clear;
MAIN_PATH = pwd;

study = 'helicon';

% PYEE_PATH = 'C:\Users\Pedro Jimenez\Documents\GitHub\pyee\';
PYEE_PATH = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\pyee\';

%% Run Pyee
i = 0;
nodesvec = [25, 50, 100, 150, 200, 250, 300];
depIn  = 0*nodesvec;
depTot = 0*nodesvec;
for nodes = nodesvec
    
    disp([num2str(100*i/length(nodesvec)),'%'])
    i = i + 1;
    
    cd(PYEE_PATH)
    flag = dos(['activate base & python run_parametric.py ', num2str(nodes)]);
    
    if flag == 1
        break
    end
    %% Postprocess
    cd(MAIN_PATH)
    getfields

    %% Compute power
    P1 = P;
    P1(isnan(P1))  = 0;
    depTot(i) = 2*pi*trapz(Zr(1,:)*0.01, trapz(Rr(:,1)*0.01, P1.*Rr(:,1)*0.01));
    P1(~PointFlag) = 0;
    depIn(i) = 2*pi*trapz(Zr(1,:)*0.01, trapz(Rr(:,1)*0.01, P1.*Rr(:,1)*0.01));   
end

figure
loglog(nodesvec.^2,depTot,'*-','LineWidth',2,'DisplayName','Total $Q_a$'); hold on
loglog(nodesvec.^2,depIn ,'*-','LineWidth',2,'DisplayName','Tube $Q_a$')
set(gca,'FontSize',30);box on;
xlabel('Number of nodes','FontSize',30,'Interpreter','latex')
ylabel('Power Deposition','FontSize',30,'Interpreter','latex')
l = legend; l.Interpreter = 'Latex';