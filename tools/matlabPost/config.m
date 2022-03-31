clear;
MAIN_PATH = pwd;

%% Config
% study = '05MextP2';
STUDY = 'fem_test';
FILENAME = getname('\results',1);
computeFlag = 0;

saveName = FILENAME(9:end-3);
saveFlag = false;

PYEE_PATH = 'C:\Users\Pedro Jimenez\Documents\GitHub\pyee\';
% PYEE_PATH = 'C:\Users\pedro\GitHub\pyee\';
FIGS_PATH = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\SPC 2020\figs\';

Pinput = 350;

%% Setup Profiles
% pre
% drawnow;
