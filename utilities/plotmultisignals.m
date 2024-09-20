% PLOTMULTISIGNALS() - Plot multiple signals spaced apart
% 
%   Usage:
%       [h,ax] = plotmultisignals(t,X,space,map)
% 
%   Inputs: 
%       t       = time vector. if  empty, uses indices
%       X       = signals (channels,samples)
%       space   = amp. space between signals. if empty, uses maximum amp.
%                   across all signals
%       map     = plot complete signals comparison (0/1). [default: 0]
% 
%   Outputs: 
%       h       = plot handles(channels,subplots)
%       ax      = axes handles(subplots)
% 
% Author: Danilo Benette Marques, 2022

function [h,ax] = plotmultisignals(t,X,space,map)

%time vector
if isempty(t)
    t = 1:size(X,2); %sampling rate=1
end

%space between signals
if nargin < 3 | isempty(space)
    space = -2*max(max(abs(X),[],1));
end

%Number of channels
NX = size(X,1);
cmap = colormap('lines');

%plot complete
if NX==1
    map = 0;
end

if nargin<4
    map = 0;
end

% % figure('color','w')
fig = gcf;
% set(fig,'position',...
%     [fig.Position(1)/2 ...
%     fig.Position(2)/2 ...
%     fig.Position(3)*2 ...
%     fig.Position(4)*2]);
movegui(fig,'onscreen')

%Plot spaced signals
if map
ax(1) = subplot(2,2,1);
end

iter = 0;
for idx_X = 1:NX
    hold on
    h(idx_X,1) = plot(t,X(idx_X,:)+iter,'color',cmap(idx_X,:));
%     text(t(end),iter,num2str(idx_X))
    
    iter = iter+space;
end
% title('Signals')
xlabel('Time')
ylabel('Amp.')
axis tight
ax(1) = gca;
pause(0)

if map
%Plot PCA
[coeff,score] = pca(X,'numcomponents',3);

ax(2) = subplot(2,2,2);
for idx_X = 1:NX
    hold on
    if size(score,2)>=3
        h(idx_X,2)=scatter3(score(idx_X,1),score(idx_X,2),score(idx_X,3),'filled','sizedata',200,'markerfacecolor',cmap(idx_X,:)); 
        text(score(idx_X,1),score(idx_X,2),score(idx_X,3),num2str(idx_X))
    else
        h(idx_X,2)=scatter(score(idx_X,1),score(idx_X,2),'filled','sizedata',200,'markerfacecolor',cmap(idx_X,:));
        text(score(idx_X,1),score(idx_X,2),num2str(idx_X))
    end
end
title('PCA')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
pause(0)

%Plot stacked signals
ax(3) = subplot(2,2,3);
for idx_X = 1:NX
    hold on
    h(idx_X,3)=plot(t,X(idx_X,:),'color',cmap(idx_X,:));
end
xlabel('Time')
ylabel('Amp.')
axis tight
pause(0)

%Plot interleaved signals
ax(4) = subplot(2,2,4);
L = size(X,2);

iter=0;
for idx_X = 1:NX
    hold on
    h(idx_X,4)=plot(t+iter,X(idx_X,:),'color',cmap(idx_X,:)); 
    plot([t(1)+iter ; t(1)+iter],[min(X,[],'all') ; max(X,[],'all')],'--k'); %grid between signals
    text(t(round(L/2))+iter,mean(X(idx_X,:)),num2str(idx_X));

    iter = iter+t(end);
end
xlabel('Time')
ylabel('Amp.')
axis tight
pause(0)

% %Plot PCA coeffs
% ax(3) = subplot(2,2,3);
% iter=0;
% for idx_pc = 1:3
%     hold on
%     h(idx_X,5)=plot(t,coeff(:,idx_pc)+iter,'color',cmap(idx_X,:));
%     
%     iter = iter-max(abs(coeff),[],'all');
% end
% title('PC coeff')
% xlabel('Time')
% ylabel('Coeff.')
% axis tight
% pause(0)

linkaxes(ax([1 3]),'x')
linkaxes(ax([3 4]),'y')
    
end