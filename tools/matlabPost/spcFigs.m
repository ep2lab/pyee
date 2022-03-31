%% Specific options for figures

% rectangle('Position',[10 0 10 4], 'FaceColor', 'black')
% rectangle('Position',[28.5 3 7 4], 'FaceColor', 'black')
% contour(f.Zr, f.Rr, double(abs(f.jz)>max(abs(f.jz(:)))*0.95),'LineWidth',7,'LineColor','magenta');
% contour(f.Zr, f.Rr, f.omega_ce, [1,1],'LineWidth',1,'LineColor','green','LineStyle','--');
% l = line([40,40],[0,3.8]); set(l, 'LineWidth',3,'Color','red','LineStyle','--')
% contour(f.Zr,f.Rr,f.PointFlag,'LineWidth',3,'LineColor','red','LineStyle','--');
% line([20 20 32.5],[0 1.18 1.18],'LineWidth',3,'Color','red');
% line([40,40,32.5],[0,3.8,1.25], 'LineWidth',3,'Color','red','LineStyle','--');
% 
% rectangle('Position',[10 0 10 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
% rectangle('Position',[28.5 3 7 4], 'EdgeColor', 'black','LineStyle','-.','LineWidth',4)
% 
% contour(f.Zr, f.Rr, double(abs(f.jz)>max(abs(f.jz(:)))*0.95),'LineWidth',5,'LineColor','magenta');
% contour(f.Zr, f.Rr, f.omega_ce, [1,1],'LineWidth',3,'LineColor','green');
% 
% lineobj = streamline(f.Zr,f.Rr,Bz,Br,40,3.6,[.5,10000]);
% lineobj.LineWidth = 2;
% lineobj.LineStyle = "--";
% lineobj.Color     = "red";