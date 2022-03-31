%% Parameters 
omega_pe = 2.5; % Electron plasma frequency
omega_pi = 0; % Ion plasma frequency
omega_ce = 2.; % Electron cyclotron frequency
omega_ci = 0.; % Ion cyclotron frequency
nu_e = 0.5; % Effective electron collisional frequency
nu_i = 0; % Effective ion collisional frequency

thetaB = 0*pi/2; % azimuth angle [rad] of B0
phiB = 0*pi/180; % elevation angle [rad] of B0

%% Kappa tensor
t.TestData.kappa = wave_explorer.permittivity.coldplasma(omega_pe,omega_pi,omega_ce,omega_ci,nu_e,nu_i,thetaB,phiB);  

[KX,KZ] = meshgrid(linspace(-5,5,100),linspace(-5,5,100));
    for i=1:length(KX(:,1))
        for j=1:length(KX(1,:))
            D(i,j) = wave_explorer.dispersion.continuous.dispersion(t.TestData.kappa,KX(i,j),0,KZ(i,j)); 
        end
    end
    figure
    contour(KX,KZ,real(D),[0,0],'color','k');
    hold on
    contour(KX,KZ,imag(D),[0,0],'color','k','linestyle','--'); 
    hold on