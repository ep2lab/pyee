function [Zpic, Rpic, Npic, Nupic, Tpic, r, M, n, plasma, r0, I0, n0, T0] = load_plume_data(ss_path, zend, rend, gamma)

    c = constants_and_units.constants;
    
    Rpic = double(h5read(ss_path,'/picM/rs')); 
    Zpic = double(h5read(ss_path,'/picM/zs')); 
    Fi_z = double(h5read(ss_path,'/ssD_picM_acc/fi_z'));
    Npic = double(h5read(ss_path,'/ssD_picM_acc/n'));
        
    iend = find(and(Zpic(2:end,1) > zend, Zpic(1:end-1,1) < zend));  
    iend = round((zend - Zpic(iend,1))/(Zpic(iend+1,1) - Zpic(iend,1))) + iend;
    
%     jend = find(and(Rpic(iend,2:end) > rend, Rpic(iend,1:end-1) < rend));  
%     jend = round((rend - Rpic(iend,jend))/(Rpic(iend,jend+1) - Rpic(iend,jend))) + jend;
%     
    Zpic = Zpic(1:iend,:);
    Rpic = Rpic(1:iend,:);
   
    Fi_z = reshape(Fi_z(:,1:iend,:),[numel(Zpic),1]);
    N    = reshape(Npic(1:iend,:),[numel(Zpic),1]);
    
    Nf    = scatteredInterpolant(Zpic(:),Rpic(:),N,'linear','nearest');
    Fi_zf = scatteredInterpolant(Zpic(:),Rpic(:),Fi_z,'linear','nearest');
    
    for i = 1:size(Zpic,1)
        for j = 1:size(Zpic,2)        
            if Zpic(i,j) >= 0.12
                Rpic(i,j) = (0.0125+(Zpic(i,j) - 0.12)/(Zpic(iend,j) - 0.12)*(rend-0.0125))*Rpic(1,j)/Rpic(1,end);
            end    
        end
    end
    
    Npic = Nf(Zpic,Rpic);
    Fi_z = Fi_zf(Zpic,Rpic);
    
    r  = Rpic(end,:);
    r0 = max(r);
    r = r/r0;
    
    elem_geom = h5read(ss_path,'/eFldM/element_geom');
    Re = double(elem_geom(:,2));
    Ze = double(elem_geom(:,1));
 
    Te   = double(h5read(ss_path,'/ssD_eFld_e_acc/Te'));
    Te = scatteredInterpolant(Ze,Re,Te);
    Tpic = Te(Zpic,Rpic);
    
    Nu   = double(h5read(ss_path,'/ssD_eFld_e_acc/freq_e_tot'));
    Nu = scatteredInterpolant(Ze,Re,Nu);
    Nupic = Nu(Zpic,Rpic);
    
    T0 = Tpic(end,1);
    n0 = Npic(end,1);
    n  = Npic(end,:)/n0;
    
    % Needed to compute sound speed
    plasma = fluid_plasma.plasma;
    plasma.ions.m = double(h5read(ss_path,'/ssIons/ssIons1/mass'));
    plasma.ions.q = c.qe;
    plasma.electrons{1}.q = -c.qe; 
    plasma.electrons{1}.T0 = c.eV2J(T0);
    plasma.electrons{1}.m  = 0;
    plasma.electrons{1}.gamma = gamma;
    plasma.electrons{1}.n0 = n0;
    
    imass = double(h5read(ss_path,'/ssIons/ssIons1/mass'));
    cs = plasma.cs(n0);    
    I0 =  imass * cs^2 / r0 /(c.qe * cs / r0); % N/A absorbs mu0

    vi_z = 0 * n;
    M    = 0 * n;
    for i = 1:size(Npic,2)
        vi_z(i) = Fi_z(end,i)/(n(i)*n0);
        M(i)    = vi_z(i)/sqrt(c.eV2J(Tpic(end,i))/imass);
    end
        
end

