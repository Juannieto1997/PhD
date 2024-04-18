function f_boxplot(v_Elec,v_US,str_y,str_name)
    v_Data = [v_Elec,v_US];
    v_ID = [zeros(1,length(v_Elec)),ones(1,length(v_US))];
    catName = categorical(v_ID,[0 1],{'Electrical' 'Ultrasound'});
    v_Green = [48 153 34]/255;
    v_Brown = [109 79 58]/255;
    figure();
    % Left axes
    tiledlayout(1,2)
    ax1 = nexttile;
    h1 = boxchart(ax1,v_Data(1,:),'GroupByColor',catName);
    h1(1).BoxFaceColor =v_Brown;
    h1(2).BoxFaceColor =v_Green;
    ylabel(str_y)
    ax1.XAxis.Color = 'none';
   ax2 = nexttile;
    h1 = boxchart(ax2,v_Data(2,:)*1000,'GroupByColor',catName);
    h1(1).BoxFaceColor =v_Brown;
    h1(2).BoxFaceColor =v_Green;
    ylabel('Delay (ms)')
    ax2.XAxis.Color = 'none';
    saveas(gcf,str_name)