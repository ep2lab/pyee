%% Load data
clear; close all

run('config.m')
filename = '\results'; 
f_c   = getfields(PYEE_PATH, 'paper\param_vac', getname(filename,6));
fA3_c = getfields(PYEE_PATH, 'paper\param_vacA3', getname(filename,6));
f_f   = getfields(PYEE_PATH, 'paper\param_vac', getname(filename,7));
fA3_f = getfields(PYEE_PATH, 'paper\param_vacA3', getname(filename,7));

field = "Qa";

%% Get an smooth curve for the contour

M = contour(f_f.Zr,f_f.Rr,f_f.PointFlag,[1,1]);

zb = M(1,2:end);
rb = M(2,2:end);
rb = rb(zb>=32.5);
zb = zb(zb>=32.5);

f1 = fit(zb(zb<40)', rb(zb<40)', 'exp2' );
z1 = unique(zb(zb<40)');
r1 = f1(z1);

f2 = fit(zb(zb>=40)', rb(zb>=40)', 'exp2' );
z2 = unique(zb(zb>=40)');
r2 = f2(z2);

z = [z1;z2];
r = [r1;r2]-0.1;

s = 0*r;
for i=2:length(s)
    s(i) = s(i-1) + norm([z(i),r(i)] - [z(i-1),r(i-1)]);  
end
% plot(zb,rb);hold on; plot(z,r);

%% Create interpolators

o_ce  = interp2(f_f.Zr,f_f.Rr,f_f.omega_ce,z,r);


vz = -0.04:0.02:0.04;
vr = -0.02:0.01:0.02;

[Zm,Rm] = meshgrid(vz,vr);


for i=1:length(z)
    q_f(i)   = sum(sum(interp2(f_f.Z.(field),f_f.R.(field),f_f.(field),z(i) + Zm, r(i) + Rm)))/numel(Zm);
    qA3_f(i) = sum(sum(interp2(fA3_f.Z.(field),fA3_f.R.(field),fA3_f.(field),z(i) + Zm,r(i) + Rm)))/numel(Zm);
    q_c(i)   = sum(sum(interp2(f_c.Z.(field),f_c.R.(field),f_c.(field),z(i) + Zm,r(i) + Rm)))/numel(Zm);
    qA3_c(i) = sum(sum(interp2(fA3_c.Z.(field),fA3_c.R.(field),fA3_c.(field),z(i) + Zm,r(i) + Rm)))/numel(Zm);
end


%% Error metrics

INT = @(Z,R,F) 2*pi*trapz(Z(1,:)*0.01, trapz(R(:,1)*0.01, abs(F).*R(:,1)*0.01));

err   = abs(qA3_f-q_c)/INT(fA3_f.Z.(field),fA3_f.R.(field),fA3_f.(field)) * INT(f_f.Z.(field),f_f.R.(field),1 + 0*f_f.(field));
errA3 = abs(qA3_f-qA3_c)/INT(fA3_f.Z.(field),fA3_f.R.(field),fA3_f.(field)) * INT(fA3_f.Z.(field),fA3_f.R.(field),1 + 0*fA3_f.(field));

% err   = abs(1-q_c./qA3_f);
% errA3 = abs(1-qA3_c./qA3_f);

er = err./errA3;


% s = s(o_ce>1);
% err = err(o_ce>1); 
% errA3 = errA3(o_ce>1);

ecr = s(o_ce>1);
ecr = ecr(end);
regs = s(z>=40);
regs = regs(1);

figure
semilogy(s,err,'LineWidth',2); hold on; semilogy(s, errA3,'LineWidth',2)
xline(ecr,'LineStyle','--','Color','g','LineWidth',2)
xline(regs,'LineStyle','--','Color','k','LineWidth',2)
set(gca,'FontSize',35)
ylabel('Relative Error','FontSize',35, 'Interpreter','latex')
xlabel('Arc Length (cm)','FontSize',35, 'Interpreter','latex')

figure
semilogy(s,er,'LineWidth',2);
xline(ecr,'LineStyle','--','Color','g','LineWidth',2)
xline(regs,'LineStyle','--','Color','k','LineWidth',2)