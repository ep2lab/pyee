clear; close all; 
% rmdir('sims','s')
%% Config
CASE_PATH = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\pyee\cases\05MextP';
SIM_PATH = 'C:\Users\pedro\OneDrive\Documentos\LocalSimData\HPT05MC0';
ss_path = [SIM_PATH,'\CORE\out\SimState.hdf5'];
zend = 0.2;
rend = 0.038;
gamma_DIM = 1.2;
gamma_post = 1.2;

%% Applied and total magnetic field
[Zpic, Rpic, Npic, Nupic, Tpic, r, M, n, plasma, r0, I0, n0, T0] = load_plume_data(ss_path, zend, rend, gamma_DIM);
mf_array = Coil_HPT05M(ss_path);

% Move the coil -> Dimagno only accepts z=0 (bug) & dimensionless
for i=1:mf_array.n_generators

mf_array.generators{i}.ZL = (mf_array.generators{i}.ZL - zend)/r0;
mf_array.generators{i}.RL = (mf_array.generators{i}.RL)/r0;
mf_array.generators{i}.I  =  mf_array.generators{i}.I / I0;

end

userdata.applied = mf_array; % Only the applied field
userdata.field = userdata.applied; % Total (applied plus plasma-induced) field

%% Plasma (Non-dimensional, plasma var is used for cs)
userdata.plasma = fluid_plasma.plasma; 
userdata.plasma.electrons{1}.gamma = gamma_DIM;

%% Initial conditions
r_ = linspace(0,1,1000);

warning off
% f  = fit(r',M','smoothingspline'); M_ = f(r_);
ft = fittype('a*exp(-(x/s)^2)', 'independent', 'x', 'dependent', 'y' );
f  = fit(r',n',ft); n_ = f(r_);
warning on

% M_ = 0*M_ + min(M_);
% M_ = 0*M_ + M_(1);
f  = fit(r',M'./1.5,'poly2'); M_ = f(r_);

figure
plot(r,M); hold on; plot(r_,M_);
figure
plot(r,n); hold on; plot(r_,n_);
hold off
userdata.ic = dimagno.ic('plasma',userdata.plasma,'field',userdata.field,'r_',r_,'ne_',{n_},'M_',M_);

%% Initial front
userdata.initialfront.i = 1; % Index for the initial front
userdata.initialfront.front = dimagno.front('r_',linspace(0,userdata.ic.r_max,200).','z_',zeros(200,1));
    userdata.initialfront.front = userdata.ic.initial_front(userdata.initialfront.front);

%% Run DIMAGNO
% dimagno.dimagno('simrc.m', userdata)

%% Post
% load('sims_back3/post.mat')

zvec  = [];
rvec  = [];
nvec  = [];
tvec  = [];

for i=1:size(Zpic,1)
    for j=1:size(Zpic,2)
        zvec(end+1) = Zpic(i,j);
        rvec(end+1) = Rpic(i,j);
        nvec(end+1) = Npic(i,j);
        tvec(end+1) = Tpic(i,j);
    end
end

for i=2:size(z,2)  % Do not duplicate the first front
    for j=1:size(z,1)
        zvec(end+1) = z(j,i)*r0 + zend;
        rvec(end+1) = r(j,i)*r0;
        nvec(end+1) = ne1(j,i)*n0;
        tvec(end+1) = T0*ne1(j,i)^((gamma_post-1)/gamma_post);
    end
end

N  = scatteredInterpolant(zvec',rvec',nvec','linear','none');
Te = scatteredInterpolant(zvec',rvec',tvec','linear','none');
Nu = scatteredInterpolant(Zpic(:),Rpic(:),Nupic(:),'linear','none');

[Z,R] = meshgrid(linspace(-0.2,0.48,500),linspace(0,max(rvec),500));
N  = N(Z,R);
Te = Te(Z,R);
Nu = Nu(Z,R);

edge.z = [0; Zpic(:,end); z(end, 2:end)'*r0 + zend; flip(z(:,end))*r0 + zend]; 
edge.r = [0; Rpic(:,end); r(end, 2:end)'*r0; flip(r(:,end))*r0];
pointF = inpolygon(Z,R,edge.z,edge.r);

N(~pointF) = min(ne1(:)) * n0;
N(N<1e14)  = 1e14;
Te(~pointF)= 0;
Tmin       = T0 * (1e14/n0).^((gamma_post-1)/gamma_post);
Te(Te<Tmin)= Tmin;

[~, TF] = rmoutliers(N,'movmedian',10);
TF(R(:,1)<0.04) = false;
ZT = Z(~TF,:);
RT = R(~TF,:);
NT = N(~TF,:);
NT = scatteredInterpolant(ZT(:),RT(:),NT(:),'linear','nearest');
N = NT(Z,R);
% N = imgaussfilt(N,1);

% Te = imgaussfilt(Te,20);
% [~, TF] = rmoutliers(Te,'movmedian',10);
% TF(R(:,1)<0.04) = false;
% ZT = Z(~TF,:);
% RT = R(~TF,:);
% TT = Te(~TF,:);
% TT = scatteredInterpolant(ZT(:),RT(:),TT(:),'linear','nearest');
% Te = TT(Z,R);

% Te = imgaussfilt(Te,5);
TUBE = Z>=0 & Z<=0.12 & R<=0.0125;
Rei = (1./Te(~TUBE)).^1.5 .* (9 + 0.5*log((1e18./N(~TUBE)).*Te(~TUBE).^3)) * 2.9*1e-12;
Nu(~TUBE)  =  N(~TUBE) .* Rei;

% [~, TF] = rmoutliers(Nu,'movmedian',10,'ThresholdFactor',2);
% TF(R(:,1)<0.1) = false;
% ZT = Z(~TF,:);
% RT = R(~TF,:);
% NuT = Nu(~TF,:);
% NuT = scatteredInterpolant(ZT(:),RT(:),NuT(:),'linear','nearest');
% Nu = NuT(Z,R);

N  = imgaussfilt(N, 0.5);
% Nu = imgaussfilt(Nu,3);

mf_array = Coil_HPT05M(ss_path);
for i=1:mf_array.n_generators

mf_array.generators{i}.ZL = mf_array.generators{i}.ZL;
mf_array.generators{i}.RL = mf_array.generators{i}.RL;
mf_array.generators{i}.I  = mf_array.generators{i}.I;

end

[~,Bz,Br] = mf_array.field_2d(Z,R);
B         = sqrt(Bz.^2+Br.^2);

figure
N2 = N; N2(~pointF) = nan;
surf(Z,R,N*0,N2);shading flat; view(2); colorbar;title('$n_e$','Interpreter','LaTex','FontSize',20)
set(gca,'ColorScale','log');
streamslice(Z,R,Bz,Br,3)

Z = Z + 0.2;
figure
surf(Z,R,N);shading flat; view(2); colorbar;title('$n_e$','Interpreter','LaTex','FontSize',20)
set(gca,'ColorScale','log');

figure
surf(Z,R,imgaussfilt(Nu, 2));shading flat; view(2); colorbar;title('$\nu_e$','Interpreter','LaTex','FontSize',20)
set(gca,'ColorScale','log');
spcFigs

figure
mesh(Zpic*100 + 20, Rpic*100, -1 + Rpic*0' +1); hold on;
mesh((z*r0*100 + zend*100) + 20, r*r0*100, -1 + 0*r); view(2); 
streamslice(Z*100,R*100,Bz,Br,3)
colormap('winter')
% title('HYPHEN (green) and DIMAGNO (blue) meshes','FontSize',20,'Interpreter','latex')
xlabel('z (cm)','FontSize',20,'Interpreter','latex')
ylabel('r (cm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20);box on;


% %% Write to file
% PATH = pwd;
% cd(CASE_PATH)
% 
% delete('profile_data.h5')
% if ~isfile('profile_data.h5')
%     h5create('profile_data.h5','/Z',size(Z))
%     h5create('profile_data.h5','/R',size(R))
%     h5create('profile_data.h5','/Bz',size(Bz))
%     h5create('profile_data.h5','/Br',size(Br))
%     h5create('profile_data.h5','/N',size(N))
%     h5create('profile_data.h5','/Nu',size(Nu))
%     h5create('profile_data.h5','/PointFlag',size(pointF))
% end
% 
% h5write('profile_data.h5','/Z', Z)
% h5write('profile_data.h5','/R', R)
% h5write('profile_data.h5','/Bz', Bz)
% h5write('profile_data.h5','/Br', Br)
% h5write('profile_data.h5','/N', N)
% h5write('profile_data.h5','/Nu', Nu)
% h5write('profile_data.h5','/PointFlag', double(pointF))
% 
% cd(PATH)
% 
% if any(isnan(Nu(:))) || any(isnan(B(:))) || any(isnan(N(:)))
%     warning('nans on matrices, check interpolators or change to nearest extrapolation')
% end