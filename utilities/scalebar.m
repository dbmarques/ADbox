% SCALEBAR() - Insert scale bars in plot
%
% 
% 
% Author: Danilo Benette Marques, 2024

function [h,htxt] = scalebar(x,xlabel,y,ylabel)

if nargin<1 | isempty(x)
    x = 1;
end

if nargin<2 | isempty(xlabel)
    xlabel = '';
end

if nargin<3 | isempty(y)
    y = 1;
end

if nargin<4 | isempty(ylabel)
    ylabel = '';
end
   
%Handles
ax = gca;
h = [];
htxt = [];

% %Set X scale bar coordinates
xpos = [ax.XLim(2)-x ax.XLim(2) ; ...
        ax.YLim(1) ax.YLim(1)];

%Plot X scale bar
hold on,h(end+1)=line(xpos(1,:),xpos(2,:),'color','k','linewidth',1);
%Plot X scale bar label
hold on,htxt(end+1)=text(mean(xpos(1,:)),xpos(2,2)-.025*(ax.YLim(2)-ax.YLim(1)),xlabel...
    ,'color','k','horizontalalignment','center');

% %Set Y scale bar coordinates
ypos = [ax.XLim(2) ax.XLim(2) ; ...
        ax.YLim(1) ax.YLim(1)+y];

%Plot Y scale bar
hold on,h(end+1)=line(ypos(1,:),ypos(2,:),'color','k','linewidth',1);
%Plot Y scale bar label
hold on,htxt(end+1)=text(ypos(1,1)+.025*(ax.XLim(2)-ax.XLim(1)),mean(ypos(2,:)),ylabel...
    ,'color','k','horizontalalignment','left');

end