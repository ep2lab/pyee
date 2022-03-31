%% Geometry
rc = 12.5;
Lc = 120;

zp  = 130;
zpu = 160;
rpu = 30;

%% Profiles
Bflag = Bflag + 0*Zm;
Nu_flag = Nu_flag + 0*Zm;
N = 0*Zm + Nout./refN;
sr = 3;
sz = 3;
for iz = 1:size(N, 2)
    for jr = 1:size(N, 1)
        zv = Zm(1,iz);
        rv = Rm(jr,1);
        if zv>=0
            if zv<=zp
                if rv<=rc
                    N(jr, iz) = exp(-abs(zv)*sz/Lc)*exp(-rv*sr/rc);
                    Nu_flag(jr, iz) = 1;
                    Bflag(jr, iz) = 1;
                end
            elseif zv<=zpu
                Rmax = rc + (rpu-rc)/(zpu-zp)*(zv-zp);
                if rv<Rmax
                    N(jr, iz) = exp(-(abs(zv))*sz/Lc)*exp(-rv*sr/rc); 
                    Nu_flag(jr, iz) = 1;
                    Bflag(jr, iz) = 1;
                end  
            else
                if rv<rpu
                    N(jr, iz) = exp(-abs(zv)*sz/Lc)*exp(-rv*sr/rc);
                    Nu_flag(jr, iz) = 1;
                    Bflag(jr, iz) = 1;
                end
            end
        end
    end
end


N = imgaussfilt(N, 2); N = N./max(N(:));

[~, ind] = min(abs(Zm(1,:) - 120));
Baxis    = BFEM(1, ind);

N = N*refN;
% Nu = Nu*refNu;
Nu = N./max(N(:))*refNu.*Nu_flag;
Bz = BzFEM./Baxis*refB.*Bflag;
Br = BrFEM./Baxis*refB.*Bflag;
B  = BFEM./Baxis*refB.*Bflag;

figure
pcolor(Zm,Rm,B*10000.*Bflag);hold on;
shading flat
colormap('jet')
c = colorbar; set(gca,'CLim',[0,3000])
caxis([0 1500])
h = streamslice(Zm,Rm,Bz,Br,3); set(h, 'LineWidth',1.5);set(h, 'Color', 'white'); hold on; 
plot3([0, 0, Lc], [0, rc, rc], [1e40,1e40, 1e40], 'k-','LineWidth',3);hold on
plot3([Lc, zp, zpu, Zm(1,end)], [rc, rc, rpu, rpu],[1e40,1e40, 1e40, 1e40],'r-','LineWidth',3);
title('Magnetic field [G]','FontSize',20,'Interpreter','latex')
xlabel('z (mm)','FontSize',20,'Interpreter','latex')
ylabel('r (mm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20)
% daspect([1, 1, 1])

figure
set(gca,'FontSize',20)
surface(Zm,Rm,N);view(2);shading flat; colorbar; set(gca,'colorscale','log');hold on;
plot3([0, 0, Lc], [0, rc, rc], [1e40,1e40, 1e40], 'k--','LineWidth',2);hold on
plot3([Lc, zp, zpu, Zm(1,end)], [rc, rc, rpu, rpu],[1e40,1e40, 1e40, 1e40],'r--','LineWidth',2);
title('Plasma Density $[m^{-3}]$','FontSize',20,'Interpreter','latex')
xlabel('z (mm)','FontSize',20,'Interpreter','latex')
ylabel('r (mm)','FontSize',20,'Interpreter','latex')
daspect([1, 1, 1])

figure
surface(Zm,Rm,Nu);view(2);shading flat; colorbar; set(gca,'colorscale','log'); hold on;
plot3([0, 0, Lc], [0, rc, rc], [1e40,1e40, 1e40], 'k--','LineWidth',4);hold on
plot3([Lc, zp, zpu, Zm(1,end)], [rc, rc, rpu, rpu],[1e40,1e40, 1e40, 1e40],'r--','LineWidth',4);
title('Collison Frequency $[\nu/\omega]$','FontSize',20,'Interpreter','latex')
xlabel('z (mm)','FontSize',20,'Interpreter','latex')
ylabel('r (mm)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',20)
daspect([1, 1, 1])

%% Write to file
cd([PYEE_PATH,'/cases/',study])

delete('magnetic_data.h5')
if ~isfile('magnetic_data.h5')
    h5create('magnetic_data.h5','/Z',size(Zm))
    h5create('magnetic_data.h5','/R',size(Rm))
    h5create('magnetic_data.h5','/BZ',size(Bz))
    h5create('magnetic_data.h5','/BR',size(Br))
    h5create('magnetic_data.h5','/N',size(N))
    h5create('magnetic_data.h5','/Nu',size(Nu))
end

h5write('magnetic_data.h5','/Z', Zm)
h5write('magnetic_data.h5','/R', Rm)
h5write('magnetic_data.h5','/BZ', Bz)
h5write('magnetic_data.h5','/BR', Br)
h5write('magnetic_data.h5','/N', N)
h5write('magnetic_data.h5','/Nu', Nu)
% 
% movefile('magnetic_data.h5',PYEE_PATH);

cd(MAIN_PATH)