%% Standard figure
function figstdAD(box)

movegui('onscreen')

set(gcf,'Color','w')

set(gca,'xcolor','k','ycolor','k'...
    ,'FontName','Helvetica'...
    ,'FontSize',12 ...
    ,'FontWeight','normal' ...
    ,'Linewidth',1 ...
    ,'Box','off' ...
    ,'Layer','top')

if nargin==1
    if box
        set(gca,'box','on')
    end
end

end