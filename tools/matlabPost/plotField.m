function [s,ax,c] = plotField(Z, R, field,varargin)

FSize = 18;

%% Parser
% Defaults
options = struct('type','magnitude', 'scale','none',...
    'title',string(),'colorlims',[0, 100], 'pointflag',0.*Z,'topo','off');
optionNames = fieldnames(options);

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = lower(pair{1}); %# make case insensitive

   if any(strcmp(inpName,optionNames))
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

%% Options
if strcmp(options.topo,'on')
    alt = 1;
else
    alt = 0;
end

if ~sum(options.pointflag(:)) == 0
    Z(~options.pointflag) = nan;
    R(~options.pointflag) = nan;
    field(~options.pointflag) = nan;
end

if strcmp(options.type, 'magnitude')
    field = abs(field);
    cmap = parula;
    clims = min(field(:)) + (max(field(:)) - min(field(:))).*options.colorlims./100;    
elseif strcmp(options.type, 'phase')
%     field = angle(field);
    cmap = hsv;
%     clims = [-pi pi];
    
    field = angle(field) * 360/(2*pi);
    clims = [-180 180];
    
    alt = 0;
end

%% Plot
% figure
if strcmp(options.scale, 'log')
%     field = abs(field);
    field(field<1e-10) = 1e-10; 
    clims = min(field(:)) + (max(field(:)) - min(field(:))).*options.colorlims./100; % Update
end
s = surf(Z, R, -1 + field*alt, field); c = colorbar;hold on; view(2);
set(s, 'FaceColor', 'flat','EdgeColor','none');
pbaspect([2, 1, 1])
xlabel('$z$(cm)','FontSize',FSize,'Interpreter','latex')
ylabel('$r$(cm)','FontSize',FSize,'Interpreter','latex')
title(options.title,'FontSize',FSize,'Interpreter','latex')
set(gca,'FontSize',FSize);
box on;ylim([min(R(:)) max(R(:))]); xlim([min(Z(:)) max(Z(:))]);
set(gca,'Colormap',cmap)
if strcmp(options.scale, 'log')
    set(gca,'colorscale','log')    
end
if clims(1) == clims(2)
    clims(2) = clims(1)*1.0001;
end
caxis(clims)
ax = gca;
end