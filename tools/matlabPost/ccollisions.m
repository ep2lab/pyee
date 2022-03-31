SIM_PATH = 'C:\Users\pedro\OneDrive\Documentos\LocalSimData\HPT05MC0';
ss_path = [SIM_PATH,'\CORE\out\SimState.hdf5'];
pd_path = [SIM_PATH,'\CORE\out\PostData.hdf5'];
cte = constants_and_units.constants;

%% Load data
Rpic = double(h5read(ss_path,'/picM/rs')); 
Zpic = double(h5read(ss_path,'/picM/zs'));
Npic = double(h5read(ss_path,'/ssD_picM_acc/n'));
NNpic = double(h5read(ss_path,'/ssD_picM_acc/nn'));
NNpic = reshape(NNpic,size(Npic));

elem_geom = h5read(ss_path,'/eFldM/element_geom');
Re = double(elem_geom(:,2));
Ze = double(elem_geom(:,1));

Te   = double(h5read(ss_path,'/ssD_eFld_e_acc/Te'));
Te = scatteredInterpolant(Ze,Re,Te);
Tpic = Te(Zpic,Rpic);

% Nu   = double(h5read(ss_path,'/ssD_eFld_e_acc/freq_e_tot'));
% Nu = scatteredInterpolant(Ze,Re,Nu);
% Nupic = Nu(Zpic,Rpic);
% 
% Nu   = double(h5read(ss_path,'/ssD_eFld_e_acc/freq_e_tot_eff_ine'));
% Nu = scatteredInterpolant(Ze,Re,Nu);
% nu_ine = Nu(Zpic,Rpic);

%% Comparing collisions
nu_en  = NNpic.*sqrt(8*cte.qe*Tpic/pi/cte.me)*27e-20;
% nu_cex = NNpic.*sqrt(cte.kB.*Tpic/22e-26)*81e-20.*(1-0.2*log10(sqrt(cte.kB.*Tpic/22e-26)/1000)).^2;
nu_ei  = Npic.*(1./Tpic).^1.5 .* (9 + 0.5*log((1e18./Npic).*Tpic.^3)) * 2.9*1e-12;

figure
surf(Zpic,Rpic,nu_ei);shading flat; view(2); colorbar;
title('$\nu_{ei}$','Interpreter','LaTex','FontSize',20)
set(gca,'ColorScale','log');

figure
surf(Zpic,Rpic,nu_en);shading flat; view(2); colorbar;
title('$\nu_{en}$','Interpreter','LaTex','FontSize',20)
set(gca,'ColorScale','log');caxis([1e4,1e8])

% figure
% surf(Zpic,Rpic,nu_cex);shading flat; view(2); colorbar;
% title('$\nu_{cex}$','Interpreter','LaTex','FontSize',20)
% set(gca,'ColorScale','log');


%% Actual collision frequencies

nu_en_pd = double(h5read(pd_path,'/picM_data/freq_e_acc1')); 
nu_ei_pd = double(h5read(pd_path,'/picM_data/freq_e_acc4'));
nu_et_pd = double(h5read(pd_path,'/picM_data/freq_e_tot_acc'));

nu_en_pd = reshape(nu_en_pd(end,:,:),size(Rpic)); 
nu_ei_pd = reshape(nu_ei_pd(end,:,:),size(Rpic)); 
nu_et_pd = reshape(nu_et_pd(end,:,:),size(Rpic));

figure
surf(Zpic,Rpic,nu_ei_pd);shading flat; view(2); colorbar;
title('$\nu_{ei}$ PD','Interpreter','LaTex','FontSize',20)
set(gca,'ColorScale','log');

figure
surf(Zpic,Rpic,nu_en_pd);shading flat; view(2); colorbar;
title('$\nu_{en}$ PD','Interpreter','LaTex','FontSize',20)
set(gca,'ColorScale','log'); caxis([1e4,1e8])

figure
surf(Zpic,Rpic,(nu_et_pd - nu_ei_pd)./nu_et_pd);shading flat; view(2); colorbar;
title('Rel err $\nu_{en}-\nu_{ei}$ PD','Interpreter','LaTex','FontSize',20)
set(gca,'ColorScale','log'); %caxis([1e6,1e10])

figure
surf(Zpic,Rpic,nu_et_pd*20);shading flat; view(2); colorbar;
title('$\nu_{et}$ wave','Interpreter','LaTex','FontSize',20)
set(gca,'ColorScale','log'); %caxis([1e6,1e10])

lnC = (9 + 0.5*log((1e18./Npic).*Tpic.^3));
figure
surf(Zpic,Rpic,Tpic);shading flat; view(2); colorbar;
title('Coulomb log','Interpreter','LaTex','FontSize',20)
% set(gca,'ColorScale','log'); %caxis([1e6,1e10])

nu_en_ref = mean(nu_en_pd(end,:));
nu_ei_ref = mean(nu_ei_pd(end,:));
% figure
% surf(Zpic,Rpic,Nupic);shading flat; view(2); colorbar;
% title('$\nu_{ine}$','Interpreter','LaTex','FontSize',20)
% set(gca,'ColorScale','log');