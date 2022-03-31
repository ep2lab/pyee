function save_fig(fig, title, PATH, saveFlag)
if saveFlag
    set(fig,'units','normalized','outerposition',[0 0 1 1],'Menu','none','ToolBar','none');
    exportgraphics(gca,[PATH, title,'.png'])
end
end