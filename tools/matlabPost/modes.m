clear; close all;
PYEE_PATH = 'C:\Users\pedro\OneDrive - Universidad Carlos III de Madrid\Codes\pyee\';
FIGS_PATH = 'C:\Users\pedro\OneDrive\Escritorio\figs\';
study = 'SPC2020\';

tags = {'results_mm3.h5', 'results_mm1.h5', 'results.h5', 'results_m3.h5'};
    
Ptv = [];
Rpv = [];
dofs= [];

Pinput = 350;
    
for t = tags
    
    FILE_NAME = t{1};
    
    getfields;
    pfigs;
    
    dofs(end+1) = 36 * (size(P,1)-1)/2 * (size(P,2)-1)/2;
    Ptv(end+1) = Ptube;
    Rpv(end+1) = Rp;

end

% figure
% a=[Rpv' zeros(4,1)];
% b=[zeros(4,1) Ptv'/350];
% [AX,H1,H2] =plotyy([-3,-1,1,3],a,[-3,-1,1,3],b, 'bar', 'bar');
% set(H1,'FaceColor','r') % a
% set(H2,'FaceColor','b') % b
% set(AX(1),'Ycolor','r','FontSize',30);% yyaxis left; ylabel('$[\Omega]$', 'Interpreter', 'latex','Rotation',0)
% set(AX(2),'Ycolor','b','FontSize',30)
% xlabel('Azimuthal mode number $m$','Interpreter','latex','FontSize',30)
figure
a=[Rpv'/sum(Rpv) zeros(4,1)];
b=[zeros(4,1) (Ptv.*Rpv/sum(Rpv))'/350];
[AX,H1,H2] =plotyy([-3,-1,1,3],a,[-3,-1,1,3],b, 'bar', 'bar');
set(H1,'FaceColor','r') % a
set(H2,'FaceColor','b') % b
set(AX(1),'Ycolor','r','FontSize',30);% yyaxis left; ylabel('$[\Omega]$', 'Interpreter', 'latex','Rotation',0)
set(AX(2),'Ycolor','b','FontSize',30)
xlabel('Azimuthal mode number $m$','Interpreter','latex','FontSize',30)

save_fig(gcf,'modes', FIGS_PATH, true)