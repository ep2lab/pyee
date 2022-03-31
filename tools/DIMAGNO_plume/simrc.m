function data = simrc(data)

%% General
data.dimagno.simdir = fullfile(pwd,'sims'); % directory where simulation files will be saved

data.dimagno.plotfront = 1; % whether to plot current front position during simulation

data.dimagno.max_iter = 1000; % maximum number of iterations
data.dimagno.tol = 1e-6; % Tolerance to exit iterations
data.dimagno.corrector = 1; % corrector in predictor-corrector on?

data.dimagno.B_on_ions = 0; % 1/0 whether to activate magnetic force on ions (if 0, B only affects electrons)
data.dimagno.beta0 = 0; % plasma beta at the origin, beta0 =  mu0 n T / Ba0^2, for iterations of the self-induced field
data.dimagno.beta_0_tolerance = 1e-4;

%% Logger 
data.logger.filedebuglevel = 3; % file debug level. A higher number prints less messages
data.logger.screendebuglevel = 3; % screen debug level. A higher number prints less messages
data.logger.linelength = 80; % maximum line length in the logs

%% Solver
data.solver.solver = 'dMoC_advance_front'; % Type of front advancer to use               
data.exit.exit = 'number_of_steps'; % Type of exit condition to use
data.exit.parameters = struct(... % Parameters for exit condition. 
            'n_steps', 1000); 

%% Postprocessor
data.postprocessor.postfunctions = {'fronts_to_arrays','induced_field'}; % Cell array with the names of postprocessor functions to run after iteration process
data.postprocessor.postparameters.induced_field = struct(...
            'Xlim',[0,0.35],...
            'Ylim',[0,1.1],...
            'xpoints',20,...
            'ypoints',20);
        
