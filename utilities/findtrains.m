% FINDTRAINS() - Find indices of a specified train of pulses by the
%   convolution of an artifact-containing signal with a simulated train.
% 
%   Usage:
%       [idx_train,C] = findtrains(X,Fs,ftrain,durtrain,mininterval)    
% 
%   Inputs:
%       X = signal
%       Fs = sampling rate
%       ftrain = frequency of train (Hz)
%       durtrain = duration of train (s)
%       mininterval = minimum interval between trains (s). [default:
%           2*durtrain]
%       thr = minimum amplitude of convoluted signal (ratio of max)
%           [default: 3*std of modulus of convoluted signal]
% 
%   Outputs
%       idx_train = indices of identified trains
%       C = convolution with simulated train
% 
% Author: Danilo Benette Marques, 2023

function [idx_train,C] = findtrains(X,Fs,ftrain,durtrain,mininterval,thr)

if exist('eegfilt')
    X = eegfilt(shiftdim(X)',Fs,200,0);
end

durtrain = durtrain*Fs;

%Simulated train
simtrain = zeros(durtrain,1);
simtrain(1:(Fs/ftrain):end) = 1;

%Convolution with simulated train
disp('Applying convolution with simulated train')
C = conv(X,simtrain,'same');
C = abs(C);

if nargin<5 | isempty(mininterval)
    mininterval = 2*durtrain;
else
    mininterval = mininterval*Fs;
end

if nargin<6 | isempty(thr)
    thr = 3*std(C);
else
    thr = thr*max(C);
end

disp('Finding train')
[~,idx_train] = findpeaks(C,'minpeakdist',mininterval,'minpeakheight',thr);
%Fix peak location from the center of train
idx_train = idx_train-(durtrain/2); %start
% idx_train = idx_train+(durtrain/2); %stop

% figure
% plot(C)
% hold on,plot(idx_train,C(idx_train),'sr')
end