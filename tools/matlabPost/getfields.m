
function f = getfields(PYEE_PATH, STUDY, FILENAME)

f.Rr = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/R').*100;
f.Zr = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/Z').*100;

PointFlag = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/PointFlag');
f.PointFlag = reshape(cellfun(@(entry)isequal(entry,'TRUE'), PointFlag), size(f.Zr));

omega_ce = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/omega_ce'); 
omega_pe = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/omega_pe'); 
nu = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/nu'); f.nu = nu.r;
thetaB = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/thetaB'); 
f.omega_ce = omega_ce.r;
f.omega_pe = omega_pe.r;
f.thetaB = thetaB.r;

% f.nthreads    = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/nthreads');
f.simTime   = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/time');
f.solveTime = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/solve_time');

f.F = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/CMA');

jz = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/Jz'); 
jy = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/Jy'); 
jx = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/Jx'); 
f.jz = jz.r + 1i*jz.i;
f.jy = jy.r + 1i*jy.i;
f.jx = jx.r + 1i*jx.i;

Ez = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/Ez'); 
Ey = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/Ey'); 
Ex = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/Ex'); 

f.Ez = (Ez.r + 1i*Ez.i);
f.Ex = (Ex.r + 1i*Ex.i);
f.Ey = (Ey.r + 1i*Ey.i);

if numel(f.Zr)>numel(f.Ez)                          % Legacy FD stagged
    f.Z.Ez = f.Zr(1:2:end,2:2:end); 
    f.R.Ez = f.Rr(1:2:end,2:2:end);
    f.Z.Ex = f.Zr(2:2:end,1:2:end);
    f.R.Ex = f.Rr(2:2:end,1:2:end);
    f.Z.Ey = f.Zr(1:2:end,1:2:end);
    f.R.Ey = f.Rr(1:2:end,1:2:end);
    f.Qa = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/P');
else                                                 % New FEM meshes
    f.Z.Ez = f.Zr; 
    f.R.Ez = f.Rr;
    f.Z.Ex = f.Zr;
    f.R.Ex = f.Rr;
    f.Z.Ey = f.Zr;
    f.R.Ey = f.Rr;
    Qa = h5read([PYEE_PATH,'cases\',STUDY,FILENAME],'/Qa');
    f.Qa = Qa.r;
end
    f.Z.Qa = f.Zr;
    f.R.Qa = f.Rr;
end