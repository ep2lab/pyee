close all;
run('config.m')
STUDY = 'paper\nom_2000filt_vac';
f = getfields(PYEE_PATH, STUDY, FILENAME);

f.Qa(f.Qa<=0) = 1e-16;

%% Filter
tic

threshold = 20;
window = [10 10];
gfilt = 0.05;

Qafilt = field_filter(f.Qa, threshold, window, 0.05);

toc

%% PLOTS
% Total resistance
res1 = 2*pi*trapz(f.Zr(1,:)*0.01, trapz(f.Rr(:,1)*0.01, f.Qa.*f.Rr*0.01));
res2 = 2*pi*trapz(f.Zr(1,:)*0.01, trapz(f.Rr(:,1)*0.01, Qafilt.*f.Rr*0.01));

figure
plotField(f.Zr, f.Rr, f.Qa, f.PointFlag*0,'title', '$Q_a [W/m^{3}]$','scale','log');
% caxis([1e-4 1])
caxis([1e-2 5e1])
clim = get(gca,'CLim');
zlim = get(gca,'zlim');
spcFigs

figure;
plotField(f.Zr, f.Rr, Qafilt, f.PointFlag*0,'title', '$Q_a [W/m^{3}]$','scale','log');
caxis(clim)
set(gca,'zlim',zlim)
spcFigs


%% Plot Noisy region 
zind = find(f.Zr(1,:)>=20 & f.Zr(1,:)<=40); rind = find(f.Rr(:,1)>=1 & f.Rr(:,1)<=3);

figure
plotField(f.Zr(rind,zind), f.Rr(rind,zind), Qafilt(rind,zind), f.PointFlag(rind,zind)*0,'title', '$Q_a [W/m^{3}]$','scale','log')
view(3)
caxis([1e-2 5e1])
clim = get(gca,'CLim');
zlim = get(gca,'zlim');
spcFigs

figure;
plotField(f.Zr(rind,zind), f.Rr(rind,zind), f.Qa(rind,zind), f.PointFlag(rind,zind)*0,'title', '$Q_a [W/m^{3}]$','scale','log');
view(3)
caxis(clim)
set(gca,'zlim',zlim)
spcFigs