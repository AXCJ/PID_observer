function plotboundary(x, xy, style, color, legendname)
    legendonoff = 'on';
    if(nargin < 5)
        legendonoff = 'off';
        legendname = '';
    end
    if(xy == 'y')
        xlimits = get(gca,'XLim');
        meanLine = line([xlimits(1) xlimits(2)],...
		 [x x],'Color', color, 'LineStyle',style, 'Displayname', legendname);
    else
        ylimits = get(gca,'YLim');
        meanLine = line([x x],...
		 [ylimits(1) ylimits(2)],'Color', color, 'LineStyle',style, 'Displayname', legendname);
    end
    
     set(get(get(meanLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle',legendonoff);
end