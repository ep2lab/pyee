clear

run('config.m')
STUDY = 'paper\param_vacA3';
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
        
%         if FILTER 
%             window      = [1 1]*ceil(size(f.(field),2)/100);
%             threshold   = (window(1)*2 + 1)/2;
%             qvec.(field){end+1} = field_filter(f.(field), threshold, window, 0.05);
%         else
%             qvec.(field){end+1} = f.(field);
%         end
        
        
        ref.(field)(end+1) = INT(f.Z.(field),f.R.(field),f.(field));
        
        in = f.PointFlag & f.Zr<30 & f.Rr<=1 & f.Zr>=21;
        in = interp2(f.Zr, f.Rr, in, f.Z.(field), f.R.(field), 'nearest');

        qvecIn.(field){end+1} = qvec.(field){end} .* in + eps;
        
        refIn.(field)(end+1) = INT(f.Z.(field),f.R.(field),f.(field).*in + eps);

    %     figure
    %     spcFigs
    %     plotField(1:size(qvec{end},2), 1:size(qvec{end},1), qvec{end}, 'title', '$Q_a [W/m^{3}]$','scale','log')
    %     caxis([1e-2 1e4])

    end
end

volL = INT(f.Z.(field),f.R.(field),0*f.Z.(field) + 1);
volS = INT(f.Z.(field),f.R.(field),in);

n = 7;
    
for i=1:n-1 
    for field = fields
        
        if i==1
            err_rms.(field) =[];
            err_rmsIn.(field) =[];
        end
        
        q1  = qvec.(field){i};
        q2  = qvec.(field){end}(1:2^(n-i):end,1:2^(n-i):end);
          
        err = (q2(:) - q1(:))/ref.(field)(end);
        err_rms.(field)(i) = rms(err);

        q1In = qvecIn.(field){i};
        q2In = qvecIn.(field){end}(1:2^(n-i):end,1:2^(n-i):end);

        err  = (q2In(:) - q1In(:))/refIn.(field)(end);
        err_rmsIn.(field)(i) = rms(err) * sqrt(numel(in)/ sum(in(:)));     
    end
    dofs(i) = 6*(size(q1,1)-1)*(size(q1,2)-1)/4;
end

color = ["red","blue","green","magenta"]; i = 1;
figure 
for field = fields
    loglog(dofs,err_rms.(field),'-+','LineWidth',2,'DisplayName',field,'Color',color(i)); hold on
    loglog(dofs,err_rmsIn.(field),'-*','LineWidth',2,'DisplayName',[char(field),' tube'],'Color',color(i));  
i = i + 1;
end
ylim([1e-3,Inf])
hold off
ylabel('Error to most accurate simulation','FontSize',16)
xlabel('DOFS','FontSize',16)
set(gca, 'YScale', 'log','FontSize',16)
l = legend; set(l,'FontSize',16)
set(gca,'FontSize',18)

figure; i = 1;
for field = fields
    loglog(dofs,abs(1-ref.(field)(1:end-1)/ref.(field)(end)),'-+','LineWidth',2,'DisplayName',field,'Color',color(i)); hold on
    loglog(dofs,abs(1-refIn.(field)(1:end-1)/refIn.(field)(end)),'-*','LineWidth',2,'DisplayName',[char(field),' tube'],'Color',color(i));  
i = i + 1;
end
ylim([1e-3,Inf])
hold off
ylabel('Error','FontSize',25)
xlabel('DOFS','FontSize',25)
set(gca, 'YScale', 'log','FontSize',25)
l = legend; set(l,'FontSize',25)
set(gca,'FontSize',25)


% figure;
% plotField(1:size(q1,2), 1:size(q1,1), abs(q1-q2)/Qaref, 'scale', 'log','topo','on','title','Normalized Error');
% caxis([1e-4,100])
% 
% disp(sum(f.nu(:)))



