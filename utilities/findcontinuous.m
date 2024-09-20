% FINDCONTINUOUS() - Find index of starts of squared pulses (e.g. TTL) in
%   signal.
% 
%   Usage
%       [idx_startstop,idx_continuous] = findcontinuous(X,thr,minduration,mininterval)
% 
%   Inputs
%       X = signal vector with squared pulses
%       thr = amplitude threshold to find continuous data above. [default: 1.5*iqr]
%       minduration = minimum duration (samples) of pulses. [default: 1]
%       maxinterval = maximum interval (samples) between end and next start to join. [default: 1]
%   Outputs
%       idx_pulses = index of starts of pulses
% 
% Autor: Danilo Benette Marques, 2024

function [idx_startstop,idx_continuous] = findcontinuous(X,thr,minduration,maxinterval)

if ~isvector(X)
    error('X must be a vector')
end
if ~isscalar(thr) | ~isscalar(minduration) | ~isscalar(maxinterval)
    error('thr, minduration, and maxinterval must be scalars')
end    
if nargin<2 | isempty(thr)
    thr = 1.5*iqr(X);
end
if nargin<3 | isempty(minduration)
    minduration = 1;
end
if nargin<4 | isempty(maxinterval)
    maxinterval = 1;
end

X = shiftdim(X);

idx_continuous = zeros(length(X),1); %allocate

% turn signal to 0/1
idx_thr = find(X>thr);
idx_thr01 = zeros(size(X));
idx_thr01(idx_thr) = 1;

idx_thr01(1)=0; idx_thr01(end)=0; %set initial and end as 0

% find start and end of continuous data above thr
diff01 = diff([0 ; idx_thr01]);
idx_thr1 = find(diff01==1); %start (0-->1)
idx_thr0 = find(diff01==-1); %end (1-->0)

% get continuous as starts to ends of every data above thr
idx_startstop = [idx_thr1 idx_thr0];

% get continuous with minimum duration
if nargin>=2 & ~isempty(minduration)
    diff10 = idx_startstop(:,2) - idx_startstop(:,1); %pulses durations (0-->1 1-->0)
    idx_thrdur = find(diff10>minduration); %find pulses longer than minduration
    
    idx_startstop = idx_startstop(idx_thrdur,:); % get pulses as starts of pulses with minimum duration
end

if isempty(idx_startstop)
    return
end

% join continuouns with short (maximum interval) intervals between
if nargin>=3 & ~isempty(maxinterval)
interval = [idx_startstop(2:end,1) - idx_startstop(1:end-1,2)]; % intervals between next start and previous end
idx_thrinter = find(interval<maxinterval);
for iinter = 1:numel(idx_thrinter)
    idx_startstop(idx_thrinter(iinter),2) = idx_startstop(idx_thrinter(iinter)+1,2); %replace end with next end
    idx_startstop(idx_thrinter(iinter)+1,:) = NaN; %erase continuous
end
idx_startstop(any(isnan(idx_startstop),2),:)=[];
end

%get indices of continuous data
for icont = 1:size(idx_startstop,1)
   idx_continuous(idx_startstop(icont,1):idx_startstop(icont,2)) = 1;
end
idx_continuous = logical(idx_continuous);


% %get pulses with minimum interval between starts
% if nargin>=3 & ~isempty(mininterval)
%     diff1 = diff(idx_pulses); %intervals between starts
%     idx_thrinter = [1 ; find(diff1>mininterval)+1]; %intervals longer than minimum interval
%     
%     idx_pulses = idx_pulses(idx_thrinter); % get pulses as starts of pulses with minimum interval
% end

end




