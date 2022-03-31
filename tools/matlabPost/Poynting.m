%% Poynting Flux
% arrDens = size(f.Ez,1);
% 
% S = NaN*ones(size(f.Ez,1),size(f.Ez,1),3);
% Snorm = NaN*ones(size(f.Ez,1),size(f.Ez,1));
% for i = 1:round(size(f.Ez,1)/arrDens):size(f.Ez,1)
%     for j = 1:round(size(f.Ez,2)/arrDens):size(f.Ez,2)
%     val   = 0.5 * real(cross([f.Ez(i,j),f.Ex(i,j),f.Ey(i,j)],conj([Bz(i,j),Bx(i,j),By(i,j)])));
%     S(i,j,:) = val;%/norm(val(1:2));
%     Snorm(i,j) = norm(val);
%     end
% end
% 
% 
% PF = (f.Zr >= 0.15*100) & (f.Zr <= 0.45*100) & (f.Rr <= 0.1*100);
% 
% f.Zr(~PF) = nan;
% f.Rr(~PF) = nan;
% 
% P(P<1e-10) = 1e-10; 
% Snorm(Snorm<1e-10) = 1e-10;
% plotField(f.Zr, f.Rr, Snorm, f.PointFlag*0, 'title', 'Power','scale','log')
% figure
% str = streamslice(f.Zr,f.Rr,S(:,:,1),S(:,:,2),7.5); hold on;
% set(str,'Color','magenta');set(str,'LineWidth',2);