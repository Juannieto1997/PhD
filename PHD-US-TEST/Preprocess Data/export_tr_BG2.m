fig = gcf;
ax = gca;
ax.FontSize = 15;
set(fig, 'color', 'none');    
set(ax, 'color', 'none');
exportgraphics(fig,'temp.eps',...   % since R2020a
    'ContentType','vector',...
    'BackgroundColor','none')