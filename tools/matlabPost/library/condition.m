close all; clc; clear;
MAIN_PATH = pwd;

PYEEPATH = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\pyee\';
FDWAVESPATH = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\fdwaves\';

%% Run Python
cd(PYEEPATH)
dos('activate base & python condition.py ')

condPyee = h5read([PYEEPATH,'condition.h5'],'/cond');
freqPyee = h5read([PYEEPATH,'condition.h5'],'/freq');

%% Run Fdwaves with double nodes
cd(FDWAVESPATH)
condFdwaves = condPyee*0;
for i=1:length(freqPyee)

    userdata = [];
    userdata.wave.omega = double(freqPyee(i))*1e6*2*pi;
    
    simrcfile='config_condition';     %If no file is set the default configuration is used
    data=fdwaves.preprocessor.preprocessor(simrcfile,userdata);

    [A,~]=fdwaves.assemble.(data.simType).coeff(data);
    
    condFdwaves(i) = condest(A);
end

%% Run Fdwaves with same nodes
condFdwaves2 = condPyee*0;
for i=1:length(freqPyee)

    userdata = [];
    userdata.wave.omega = double(freqPyee(i))*1e6*2*pi;
    
    simrcfile='config_condition2';     %If no file is set the default configuration is used
    data=fdwaves.preprocessor.preprocessor(simrcfile,userdata);

    [A,~]=fdwaves.assemble.(data.simType).coeff(data);
    
    condFdwaves2(i) = condest(A);
end

cd(MAIN_PATH)
save('condition.mat', 'freqPyee','condPyee','condFdwaves','condFdwaves2')

%% Plots
load condition.mat

figure
loglog(double(freqPyee)*1e6, condPyee, '+-', 'LineWidth', 2,'DisplayName', 'Pyee $10\times 10$');
hold on;
loglog(double(freqPyee)*1e6, condFdwaves, '*-', 'LineWidth', 2,'DisplayName', 'Sofd5 $10\times 10$');
hold on;
loglog(double(freqPyee)*1e6, condFdwaves2,'o-', 'LineWidth', 2,'DisplayName', 'Sofd5 $20\times 20$');
lgd = legend;
set(lgd, 'Interpreter', 'LaTex')
set(gca,'FontSize',24)
xlabel('Frequency [Hz]','Interpreter','Latex')
ylabel('Condition Number', 'Interpreter','Latex')

% figure
% loglog(double(freqPyee)*1e6,condFdwaves,'-*', double(freqPyee)*1e6, condFdwaves2,'-*')
% lgd = legend('Sofd5 $70\times 70$', 'Sofd5 $35\times 35$'); set(lgd,'Interpreter','Latex')
% set(gca,'FontSize',24)
% set(gca,'FontSize',24)
% xlabel('Frequency [Hz]','Interpreter','Latex')
% ylabel('Condition Number', 'Interpreter','Latex')
