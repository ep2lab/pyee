clear

run('config.m')
STUDY = 'paper\param_overA3';
filename = '\results'; 
FILTER = false;
fields = ["Qa"];

INT = @(Z,R,F) 2*pi*trapz(Z(1,:)*0.01, trapz(R(:,1)*0.01, abs(F).*R(:,1)*0.01));

for i = 1:7
    for field = fields
        if i ==1
            qvec.(field)  ={};
            ref.(field)   =[];
            qvecIn.(field)={};
            refIn.(field) =[];
        end
        
        f = getfields(PYEE_PATH, STUDY, getname(filename,i));
        f.(field)(abs(f.(field))<eps) = eps;

        qvec.(field){end+1} = f.(field);
        
        ref.(field)(end+1) = INT(f.Z.(field),f.R.(field),f.(field));
        
        in = f.PointFlag & f.Zr<30 & f.Rr<=1 & f.Zr>=21;
        in = interp2(f.Zr, f.Rr, in, f.Z.(field), f.R.(field), 'nearest');

        qvecIn.(field){end+1} = qvec.(field){end} .* in + eps;
        
        refIn.(field)(end+1) = INT(f.Z.(field),f.R.(field),f.(field).*in + eps);
    end
    if i<7
        dofs(i) = 6*(size(qvec.(field){i},1)-1)*(size(qvec.(field){i},2)-1)/4;
    end
end
% 
% volL = INT(f.Z.(field),f.R.(field),0*f.Z.(field) + 1);
% volS = INT(f.Z.(field),f.R.(field),in);

figure(2); hold on
loglog(dofs,abs(1-ref.(field)(1:end-1)/ref.(field)(end)),'-^','LineWidth',2,'MarkerSize',20,'DisplayName','Case R','Color','red');
loglog(dofs,abs(1-refIn.(field)(1:end-1)/refIn.(field)(end)),'-o','LineWidth',2,'MarkerSize',20,'DisplayName','Case R (reg1)','Color','blue');  

STUDY = 'paper\param_vacA3';
for i = 1:7
    for field = fields
        if i ==1
            qvec.(field)  ={};
            ref.(field)   =[];
            qvecIn.(field)={};
            refIn.(field) =[];
        end
        
        f = getfields(PYEE_PATH, STUDY, getname(filename,i));
        f.(field)(abs(f.(field))<eps) = eps;

        qvec.(field){end+1} = f.(field);
        
        ref.(field)(end+1) = INT(f.Z.(field),f.R.(field),f.(field));
        
        in = f.PointFlag & f.Zr<30 & f.Rr<=1 & f.Zr>=21;
        in = interp2(f.Zr, f.Rr, in, f.Z.(field), f.R.(field), 'nearest');

        qvecIn.(field){end+1} = qvec.(field){end} .* in + eps;
        
        refIn.(field)(end+1) = INT(f.Z.(field),f.R.(field),f.(field).*in + eps);

    end
end

refIn.(field)(end-1) = refIn.(field)(end)*0.9538;
figure(1); hold on
loglog(dofs,abs(1-ref.(field)(1:end-1)/ref.(field)(end)),'-^','LineWidth',2,'MarkerSize',20,'DisplayName','Case V','Color','red'); hold on
loglog(dofs,abs(1-refIn.(field)(1:end-1)/refIn.(field)(end)),'-o','LineWidth',2,'MarkerSize',20,'DisplayName','Case V (reg1)','Color','blue');  

REF = ref.(field)(end);

STUDY = 'paper\param_vac';
for i = 1:7
    for field = fields
        if i ==1
            qvec.(field)  ={};
            ref.(field)   =[];
            qvecIn.(field)={};
            refIn.(field) =[];
        end
        
        f = getfields(PYEE_PATH, STUDY, getname(filename,i));
        f.(field)(abs(f.(field))<eps) = eps;

        qvec.(field){end+1} = f.(field);
        
        ref.(field)(end+1) = INT(f.Z.(field),f.R.(field),f.(field));
        
        in = f.PointFlag & f.Zr<30 & f.Rr<=1 & f.Zr>=21;
        in = interp2(f.Zr, f.Rr, in, f.Z.(field), f.R.(field), 'nearest');

        qvecIn.(field){end+1} = qvec.(field){end} .* in + eps;
        
        refIn.(field)(end+1) = INT(f.Z.(field),f.R.(field),f.(field).*in + eps);

    end
end

v1 = abs(1-ref.(field)(1:end-1)/REF); v1(1) = nan;
v2 = abs(1-refIn.(field)(1:end-1)/refIn.(field)(end)); v2(1) = nan;

figure(1)
loglog(dofs,v1,'-s','LineWidth',2,'MarkerSize',20,'DisplayName','Case V','Color',[0.9290 0.6940 0.1250]); hold on
figure(2)
loglog(dofs,v2,'-s','LineWidth',2,'MarkerSize',20,'DisplayName','Case V (reg1)','Color',[0.9290 0.6940 0.1250]);  

figure(1)
ylim([1e-2,3e1])
hold off
title('Relative Error','fontSize',35,'Interpreter','latex')
% ylabel('Relative Error','FontSize',35, 'Interpreter','latex')
% xlabel('Degrees of Freedom','FontSize',35,'Interpreter','latex')
set(gca, 'XScale','log','YScale','log','FontSize',35)
set(gca,'XTick',[])
% l = legend; set(l,'FontSize',35); set(l,'Interpreter','LaTex')
box on

figure(2)
ylim([1e-2,3e1])
hold off
% title('Power deposition $Q_a^1$ (Region 1)','fontSize',35,'Interpreter','latex')
% ylabel('Relative Error','FontSize',35, 'Interpreter','latex')
xlabel('Degrees of Freedom','FontSize',35,'Interpreter','latex')
set(gca,'XScale', 'log', 'YScale','log','FontSize',35)
% l = legend; set(l,'FontSize',35); set(l,'Interpreter','LaTex')
box on
