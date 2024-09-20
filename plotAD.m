% PLOTAD() - Find indices of afterdischarges (AD) based on a semiautomatic 
%   detection based on baseline-referenced energy amplitude threshold
%   followed by manual fixes.
% 
%   Usage:
%       [Idx_ad] = findAD(Lfp,Fs,pretrain,traindur,posttrain,showfigs)    
% 
%   Inputs:
%       Lfp = signals to locate AD (trains x samples x regions x days x
%           subjects)
%       Fs = sampling rate
%       pretrain = last sample of pre-train epoch (sample)
%       traindur = duration of train to cut (samples)
%       posttrain = duration of post-train epoch (samples)
% 
%   Outputs
%       Idx_ad = cell array containing [start stop] indices for every AD of 
%           each train for every dimension
% 
% Author: Danilo Benette Marques, 2024

function [h] = plotAD(Lfp,Idx_ad,Fs,pretrain,traindur,posttrain)

for idx_subj = 1:size(Lfp,5)
    for idx_day = 1:size(Lfp,4)
        for idx_region = 1:size(Lfp,3)

            lfp = Lfp(:,:,idx_region,idx_day,idx_subj);
            if isempty(lfp) | all(isnan(lfp),'all')
                continue %next loop itesubjion if no data
            end
            
            %Separate pre- and post-
            lfp_pre = lfp(:,1:pretrain); %first 1 min+10,1 seg
            lfp_pos = lfp(:,pretrain+traindur+1 : pretrain+traindur+posttrain); %>1 min +10.1(train) 
            
            %Turn NaNs to 0
            isnan_lfp_pre = isnan(lfp_pre);
            isnan_lfp_pos = isnan(lfp_pos);
            lfp_pre(isnan(lfp_pre))=0;
            lfp_pos(isnan(lfp_pos))=0;
            
            %Filter
            locutoff = 2;
            hicutoff = 50;
            lfp_pre = eegfilt(lfp_pre,Fs,locutoff,0);
                lfp_pre = eegfilt(lfp_pre,Fs,0,hicutoff);
            lfp_pos = eegfilt(lfp_pos,Fs,locutoff,0);
                lfp_pos = eegfilt(lfp_pos,Fs,0,hicutoff);

            %Return NaN
            lfp_pre(isnan_lfp_pre) = NaN;
            lfp_pos(isnan_lfp_pos) = NaN;
            
            %Reconcantenate
            lfp = [lfp_pre lfp_pos];
            t = linspace(0,size(lfp,2)/Fs,size(lfp,2)); %time vector (s) as in app
    
            %AD detection measure
            movingwin = 10*Fs;
            lfp_amp = movmean(lfp.*lfp,movingwin,2); %energy mean
            lfp_amp_pre = movmean(lfp_pre.*lfp_pre,movingwin,2);
            
            %Threshold for AD detection
            thr_amp_ad = 1.5*median(lfp_amp_pre(:));
    %         thr_amp_ad = median(lfp_amp_pre(:))% + 3*mad(lfp_amp_pre,[1],'all');
    %         thr_amp_ad = mean(lfp_amp_pre(:)) + 1*std(lfp_amp_pre(:));
    %         thr_amp_ad = 3*iqr(lfp_amp_pre(:));
    %         thr_amp_ad = quantile(lfp_amp_pre(:),.75) + 1.5*iqr(lfp_amp_pre(:));
    %         thr_amp_ad = 2*rms(lfp_amp_pre(:));
            
            %General session figure
            lfp_ad = abs(lfp);
            for idx_train = 1:size(lfp,1)
                idx_ad = Idx_ad{idx_train,1,idx_region,idx_day,idx_subj} + pretrain;
                idx_ad01 = zeros(size(lfp,2),1);
                for iad = 1:size(idx_ad,1)
                    idx_ad01(idx_ad(iad,1):idx_ad(iad,2))=1;
                end
                
            lfp_ad(idx_train,idx_ad01==0)=NaN; %hide non-AD
            end %train                

            %All trains figure
            figure('color','w')
            hold on,h(:,1)=plotmultisignals(t,abs(lfp),-2); set(h(:,1),'color',[.8 .8 .8]) %lfp modulus
            hold on,h(:,2)=plotmultisignals(t,lfp,-2); set(h(:,2),'color','k') %lfp
%             hold on,h(:,3)=plotmultisignals(t,lfp_amp,[]); set(h(:,3),'color','b') %lfp energy    
            hold on,h(:,4)=plotmultisignals(t,movmax(abs(lfp_ad),1*Fs,2),-2); set(h(:,4),'color','r','linewidth',2); %lfp amp
            g0=gridxy([pretrain/Fs (pretrain+traindur)/Fs],[],'linestyle','--','color','k'); uistack(g0,'top')    
            scalebar(60,'1 min',1,'1 mV'); figstdAD; set(gca,'xcolor','none','ycolor','none');           
            
            title(['subject: ' num2str(idx_subj) ', day: ' num2str(idx_day) ', region: ' num2str(idx_region)])
            pause(0)

            
        end %region
    end %day
end %subj


