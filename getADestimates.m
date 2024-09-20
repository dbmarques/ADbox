% GETADESTIMATES() - Get conventional estimates of afterdischarges (AD),
%   per train, such as quantities, latency to first, durations, and 
%   amplitude. Also retrieves the amplitude of LFP and the LFP only during
%   ADs.
% 
%   Usage:
%       [ADqnt, ADlat, ADdur, ADamp, Lfp_amp, Lfp_ad] = getADestimates(Lfp,Idx_ad,Fs,pretrain,traindur,posttrain)
% 
%   Inputs:
%       Lfp = signals to locate AD (trains x samples x regions x days x
%           subjects). Works with up to 5 dimensions
%       Idx_ad = cell array of indices of ADs retrieved by findAD
%       Fs = sampling rate
%       pretrain = last sample of pre-train epoch (sample)
%       traindur = duration of train to cut (samples)
%       posttrain = duration of post-train epoch (samples)
% 
%   Outputs
%       ADqnt = cell array of AD quantities per train
%       ADlat = cell array of latency to first AD per train. If none, 
%           returns NaN
%       ADdur = cell array of AD durations per train
%       ADamp = cell array of AD amplitudes per train. Amplitudes are
%           calculated by 1-sec moving maximum of LFP modulus.
%       Lfp_amp = array of cut LFP amplitudes. Amplitudes are calculated by 
%           1-sec moving maximum of LFP modulus.
%       Lfp_ad = array of cut LFP during ADs with else hidden (NaN). Tip: 
%           use ~isnan(Lfp_ad) to get logical indices
% 
% Author: Danilo Benette Marques, 2024

function [ADqnt, ADlat, ADdur, ADamp, Lfp_amp, Lfp_ad] = getADestimates(Lfp,Idx_ad,Fs,pretrain,traindur,posttrain)

Ntrain = size(Lfp,1); %number of trains

%Allocate all with NaN (if empty or missing keeps NaN)
ADqnt(1:Ntrain,1,size(Lfp,3),size(Lfp,4),size(Lfp,5)) = {[NaN]};                
ADlat(1:Ntrain,1,size(Lfp,3),size(Lfp,4),size(Lfp,5)) = {[NaN]};               
ADdur(1:Ntrain,1,size(Lfp,3),size(Lfp,4),size(Lfp,5)) = {[NaN]};       
ADamp(1:Ntrain,1,size(Lfp,3),size(Lfp,4),size(Lfp,5)) = {[NaN]};  

Lfp_amp = nan(size(Lfp,1),pretrain+posttrain,size(Lfp,3),size(Lfp,4),size(Lfp,5));
Lfp_ad = nan(size(Lfp,1),pretrain+posttrain,size(Lfp,3),size(Lfp,4),size(Lfp,5));

%Run for every dimension above trains
for idx_subj = 1:size(Lfp,5)
    for idx_day = 1:size(Lfp,4)
        for idx_region = 1:size(Lfp,3)
        
        lfp = Lfp(:,:,idx_region,idx_day,idx_subj);
        
        if isempty(lfp) | all(isnan(lfp),'all') %obs: actually no need for isempty anymore. not cell
           continue
        end
        
        disp(['Getting AD estimates of subject: ' num2str(idx_subj) ', day: ' num2str(idx_day) ', region: ' num2str(idx_region)])

        %Sepasubje pre- and post-
        lfp_pre = lfp(:,1:pretrain); %first 1 min+10,1 seg
        lfp_pos = lfp(:,pretrain+traindur+1 : (pretrain+traindur+posttrain)); %>1 min +10.1(train) 
        
        %Turn NaNs to 0 (for filtering)
        isnan_lfp_pre = isnan(lfp_pre);
        isnan_lfp_pos = isnan(lfp_pos);
        lfp_pre(isnan(lfp_pre))=0;
        lfp_pos(isnan(lfp_pos))=0;
        
        %Filter
        locutoff = 2;
        hicutoff = 20;
        lfp_pre = eegfilt(lfp_pre,Fs,locutoff,0);
            lfp_pre = eegfilt(lfp_pre,Fs,0,hicutoff);
        lfp_pos = eegfilt(lfp_pos,Fs,locutoff,0);
            lfp_pos = eegfilt(lfp_pos,Fs,0,hicutoff);
        
        %Reconcatenate
        lfp = [lfp_pre lfp_pos];
        t = linspace(0,size(lfp,2)/Fs,size(lfp,2)); %time vector (s) as in app
        
        %Return NaN
%         idx_nan = Idx_noise{idx_subj,idx_day};
%         lfp(idx_nan) = NaN;

        %AD detected
        idx_ad = Idx_ad(:,1,idx_region,idx_day,idx_subj);
        idx_ad = cellfun(@(x) x+pretrain,idx_ad,'un',0);

        %LFP amplitude
%         movingwin = 10*Fs;
%         lfp_amp = abs(lfp);
%         lfp_amp = movmean(lfp.*lfp,movingwin,2); %energy mean
        lfp_amp = movmax(abs(lfp),1*Fs,2);

         %Run for each AD of every train
        lfp_ad = nan(size(lfp));
        for idx_train = 1:size(lfp,1) %each train
                     
            %AD quantities
            ADqnt{idx_train,1,idx_region,idx_day,idx_subj} = size(idx_ad{idx_train},1); %0 if none
            
            %if no AD or no info in train
            if isempty(idx_ad{idx_train}) %if no AD in train, fill with 0 or NaN
                ADlat{idx_train,1,idx_region,idx_day,idx_subj} = [NaN]; %"infinite" latency            
                ADdur{idx_train,1,idx_region,idx_day,idx_subj} = [0]; %"zero" duration                
                ADamp{idx_train,1,idx_region,idx_day,idx_subj} = [0]; %"zero" amplitude
            elseif all(isnan(idx_ad{idx_train})) %if no AD info, fill with NaN (don't know)
                idx_ad{idx_train} = []; %make empty to not run ADs loop
                ADqnt{idx_train,1,idx_region,idx_day,idx_subj} = [NaN]; %NaN (replaces 1)
                ADlat{idx_train,1,idx_region,idx_day,idx_subj} = [NaN]; %NaN          
                ADdur{idx_train,1,idx_region,idx_day,idx_subj} = [NaN]; %NaN         
                ADamp{idx_train,1,idx_region,idx_day,idx_subj} = [NaN]; %NaN
            end
            
            for iad = 1:size(idx_ad{idx_train},1) %each train's AD
                idx_ad_ad = idx_ad{idx_train}(iad,:); %indices of current AD
                idx01_ad_ad = idx2logical(idx_ad_ad(1):idx_ad_ad(2),size(lfp,2)); %indices 0/1 of current AD

                lfp_ad(idx_train,idx_ad_ad(1):idx_ad_ad(2)) = movmax(abs(lfp(idx_train,idx_ad_ad(1):idx_ad_ad(2))),1*Fs);
%                 lfp_ad(idx_train,idx01_ad_ad==0) = NaN;     

                lfp_ad_ad = lfp(idx_train,idx_ad_ad(1):idx_ad_ad(2)); %lfp of specific AD

                %Amplitude
%                 lfp_ad_ad_amp = abs(lfp_ad_ad); %modulus
%                 lfp_ad_ad_amp = lfp_ad_ad.*lfp_ad_ad; %energy
                lfp_ad_ad_amp = movmax(abs(lfp_ad_ad),1*Fs); %smoothed max

                %Calculate AD estimates
                %Latency of first
                ADlat{idx_train,1,idx_region,idx_day,idx_subj} = (idx_ad_ad(1)-pretrain)/Fs;

                %Dusubjion
                ADdur{idx_train,1,idx_region,idx_day,idx_subj}(iad,:) = (idx_ad_ad(2)-idx_ad_ad(1))/Fs;

                %Amplitude
                ADamp{idx_train,1,idx_region,idx_day,idx_subj}(iad,:) = nanmean(lfp_ad_ad_amp);

            end %ad
        end %train
        
        Lfp_amp(:,:,idx_region,idx_day,idx_subj) = lfp_amp;
        Lfp_ad(:,:,idx_region,idx_day,idx_subj) = lfp_ad;
        
        end %region
    end %day
end %subj


end