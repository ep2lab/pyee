% run('config.m')
STUDY = 'paper/param_overA3';
FILENAME = getname('\results',7);
f = getfields(PYEE_PATH, STUDY, FILENAME);


% Field and Region def
close all;
field  = 'Ez';
m = 1;
Nth = 157;
L = 10;
w = 5;

myZ     = [f.Z.(field)(2:end,:); f.Z.(field)];
myR     = [-flip(f.R.(field)(2:end,:)); f.R.(field)];
myfield = [-flip(exp(1i*pi*m)*f.(field)(2:end,:)); f.(field)];
% myfield = exp(-20*1i*(myZ + 0*myR)); 

close all;

fig1 = figure;
plotField(myZ, myR, myfield, 'type','magnitude','scale','log');
set(gcf,'units','normalized','outerposition',[0 0.3 0.6 0.6]);

fig2 = figure;
set(gcf,'units','normalized','outerposition',[0.6 0.3 0.4 0.6]);

theta = linspace(0,pi,Nth);

d = max(myZ(1,2)-myZ(1,1),myR(2,1)-myR(1,1));

[Ztemp, Rtemp] = ndgrid(-L/2:d:L/2,-w/2:0.3*d:w/2);

cir1  = []; 
cir2  = [];
pnt   = [];
val   = [];
kz    = [];
kp    = [];

F = griddedInterpolant(myZ',myR',myfield');

Epar = interp2(f.Z.('Ez'), f.R.('Ez'), f.('Ez'), f.Zr, f.Rr).*cos(f.thetaB) + interp2(f.Z.('Ex'), f.R.('Ex'), f.('Ex'), f.Zr, f.Rr).*sin(f.thetaB);
Eper = -interp2(f.Z.('Ez'), f.R.('Ez'), f.('Ez'), f.Zr, f.Rr).*sin(f.thetaB) + interp2(f.Z.('Ex'), f.R.('Ex'), f.('Ex'), f.Zr, f.Rr).*cos(f.thetaB);

figure
plotField(f.Zr, f.Rr, Eper, 'type','magnitude','scale','log');
figure
plotField(f.Zr, f.Rr, Eper, 'type','phase');
figure
plotField(f.Z.('Ex'), f.R.('Ex'), f.('Ex'),'scale','log');caxis([1e-3,1e3])


% Test it (uses ndgrid)
% [Z,R] = ndgrid(linspace(0,67,500),linspace(0,20,500));
% plotField(Z, R, F(Z,R), 'type','phase');

while 1

    figure(fig1);
    hold on
    [zo,ro] = getpts;
    delete(cir1);delete(cir2); delete(pnt);
    cir1 = rectangle('Position',[zo-L/2 ro-L/2 L L],'Curvature',[1 1],'LineWidth',2);
    cir2 = rectangle('Position',[zo-L/2 ro-w/2 L w],'LineWidth',2);
    pnt = plot(zo,ro,'+','MarkerSize',5);
    
    for i = 1:Nth
        
    % Rotate the rectangle
    Zrot = zo + cos(theta(i)) .* Ztemp - sin(theta(i)) .* Rtemp;
    Rrot = ro + cos(theta(i)) .* Rtemp + sin(theta(i)) .* Ztemp;
    
    samp = sum(F(Zrot,Rrot),2);
    % Test it (uses ndgrid)
%     figure
%     plotField(Zrot, Rrot, F(Zrot,Rrot), 'type','phase');

    
    [wavelet,fs] = cwt(samp);
    
%     disp(size(fs))
    
    val(i,:) = wavelet(:,round(length(samp)/2),1);
    val(2*Nth - i + 1,:) = wavelet(:,round(length(samp)/2),2);
        
    kz(i,:) = cos(theta(i)) .* fs / d * 100 * 2*pi;
    kp(i,:) = sin(theta(i)) .* fs / d * 100 * 2*pi;
    
    kz(2*Nth - i + 1,:) = cos(pi + theta(i)) .* fs / d * 100 * 2*pi;
    kp(2*Nth - i + 1,:) = sin(pi + theta(i)) .* fs / d * 100 * 2*pi;
    
    end

    figure(fig2)
    p1 = surf(kz,kp,0.*kp,abs(val));hold on; axis equal; colorbar;
    shading flat;  view(2); set(gca,'ColorScale','Log');
    
    mag_angle = interp2(f.Zr,f.Rr,f.thetaB,zo,ro); zmax = max(kz(:));
    quiver(0,0,zmax*cos(mag_angle), zmax*sin(mag_angle),'LineWidth',2,'Color','magenta')
    
    hold off;
     
end
