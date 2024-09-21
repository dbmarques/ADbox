% FINDAD() - Find indices of afterdischarges (AD) based on a semiautomatic 
%   detection based on baseline-referenced energy amplitude threshold
%   followed by manual fixes.
% 
%   Usage:
%       [Idx_ad] = findAD(Lfp,Fs,pretrain,traindur,posttrain,method,showfigs)    
% 
%   Inputs:
%       Lfp = signals to locate AD (trains x samples x regions x days x
%           subjects)
%       Fs = sampling rate
%       pretrain = last sample of pre-train epoch (sample)
%       traindur = duration of train to cut (samples)
%       posttrain = duration of post-train epoch (samples)
%       method = method to estimate lfp amplitude {'energy'
%           [default],'amplitude','neo'}. Energy estimates the 10-sec
%           moving mean of the signal energy (mV^2). Amplitude estimates
%           the 1-sec moving maximum of the modulus of the signal. NEO
%           estimates the energy and the Nonlinear Energy Operator and its
%           amplitude and asymmetry, then gets a general score.
% 
%   Outputs
%       Idx_ad = cell array containing [start stop] indices for every AD of 
%           each train for every dimension
% 
%   Requirements:
%       eegfilt
%       findcontinuous
%       idx2logical
%       gridxy
%       plotmultisignals
%       signaltsmarker
% 
% Author: Danilo Benette Marques, 2024
% Last update: 2024-09-19

function [Idx_ad] = findAD(Lfp,Idx_ad,Fs,pretrain,traindur,posttrain,method,showfigs)

%Use presaved input Idx_ad
if ~isempty(Idx_ad) 
    presaved = 1;
else 
    presaved = 0;
end

if nargin<7 | isempty(method)
    method = 'energy';
end

for idx_subj = 1:size(Lfp,5)
    for idx_day = 1:size(Lfp,4)
        for idx_region = 1:size(Lfp,3)

            lfp = Lfp(:,:,idx_region,idx_day,idx_subj);
            if isempty(lfp) | all(isnan(lfp),'all')
                continue %next loop itesubjion if no data
            end
    %         t = linspace(0,size(lfp,2)/Fs,size(lfp,2)); %time vector (s)
            
            %Separate pre- and post-
            lfp_pre = lfp(:,1:pretrain); %first 1 min+10,1 seg
            lfp_pos = lfp(:,pretrain+traindur+1 : pretrain+traindur+posttrain); %>1 min +10.1(train) 
            
            %Turn NaNs to 0
            isnan_lfp_pre = isnan(lfp_pre);
            isnan_lfp_pos = isnan(lfp_pos);
            lfp_pre(isnan(lfp_pre))=0;
            lfp_pos(isnan(lfp_pos))=0;
            
            %Automatic detection of ADs
            switch method
                case 'energy'

                    %Filter
                    locutoff = 2;
                    hicutoff = 20; %overall amplitude
                    lfp_pre = eegfilt(lfp_pre,Fs,locutoff,0);
                        lfp_pre = eegfilt(lfp_pre,Fs,0,hicutoff);
                    lfp_pos = eegfilt(lfp_pos,Fs,locutoff,0);
                        lfp_pos = eegfilt(lfp_pos,Fs,0,hicutoff);
    
                    %Reconcantenate
                    lfp = [lfp_pre lfp_pos];
                    t = linspace(0,size(lfp,2)/Fs,size(lfp,2)); %time vector (s) as in app
    
                    %Energy
                    lfp_amp = lfp.*lfp;
                    lfp_amp_pre = lfp_pre.*lfp_pre;

                    %Smooth
                    movingwin = 10*Fs;
                    lfp_amp = movmean(lfp_amp,movingwin,2); 
                    lfp_amp_pre = movmean(lfp_amp_pre,movingwin,2);

                case 'amplitude'
                    
                    %Filter
                    locutoff = 2;
                    hicutoff = 20; %overall amplitude
                    lfp_pre = eegfilt(lfp_pre,Fs,locutoff,0);
                        lfp_pre = eegfilt(lfp_pre,Fs,0,hicutoff);
                    lfp_pos = eegfilt(lfp_pos,Fs,locutoff,0);
                        lfp_pos = eegfilt(lfp_pos,Fs,0,hicutoff);
    
                    %Reconcantenate
                    lfp = [lfp_pre lfp_pos];
                    t = linspace(0,size(lfp,2)/Fs,size(lfp,2)); %time vector (s) as in app
    
                    %Amplitude
                    movingwin = 1*Fs;
                    lfp_amp = movmax(abs(lfp),movingwin,2);
                    lfp_amp_pre = movmax(abs(lfp_pre),movingwin,2);

                case 'neo'

                    %Filter
                    locutoff = 2;
                    hicutoff = 50; %preserves spike waveforms
                    lfp_pre = eegfilt(lfp_pre,Fs,locutoff,0);
                        lfp_pre = eegfilt(lfp_pre,Fs,0,hicutoff);
                    lfp_pos = eegfilt(lfp_pos,Fs,locutoff,0);
                        lfp_pos = eegfilt(lfp_pos,Fs,0,hicutoff);
    
                    %Reconcantenate
                    lfp = [lfp_pre lfp_pos];
                    t = linspace(0,size(lfp,2)/Fs,size(lfp,2)); %time vector (s) as in app
    
                    %Energy
                    movingwin = 1*Fs;
                    lfp_eng = movmean(lfp.*lfp,movingwin,2);
                    lfp_eng_pre = movmean(lfp_pre.*lfp_pre,movingwin,2);
        
                    %NEO
                    d = 50;
                    lfp_neo = neo(lfp',d)';
                    lfp_neo_pre = neo(lfp_pre',d)';
        
                    %Amplitude
                    movingwin = 1*Fs;
                    lfp_amp = movmax(abs(lfp_neo),movingwin,2);
                    lfp_amp_pre = movmax(abs(lfp_neo_pre),movingwin,2);
    
                    %Asymmetry
                    movingwin = 1*Fs;
                    lfp_asym = movmax(lfp_neo,movingwin,2) - abs(movmin(lfp_neo,movingwin,2));
                    lfp_asym_pre = movmax(lfp_neo_pre,movingwin,2) - abs(movmin(lfp_neo_pre,movingwin,2));
        
                    %General AD score
                    lfp_eng = normalize(lfp_eng,2);
                    lfp_eng_pre = normalize(lfp_eng_pre,2);
                    lfp_amp = normalize(lfp_amp,2);
                    lfp_amp_pre = normalize(lfp_amp_pre,2);
                    lfp_asym = normalize(lfp_asym,2);
                    lfp_asym_pre = normalize(lfp_asym_pre,2);
    
                    lfp_amp = + lfp_amp + lfp_asym;
                    lfp_amp_pre = + lfp_amp_pre + lfp_asym_pre;

                    %Smooth score
                    movingwin = 10*Fs;
                    lfp_amp = movmean(lfp_amp,movingwin,2); 
                    lfp_amp_pre = movmean(lfp_amp_pre,movingwin,2);
            end
            
            %Threshold for AD detection
            thr_amp_ad = 1.5*median(lfp_amp_pre(:));
%             thr_amp_ad = median(lfp_amp_pre(:)) + 3*mad(lfp_amp_pre,1,'all');
%             thr_amp_ad = mean(lfp_amp_pre(:)) + 2*std(lfp_amp_pre(:)):
%             thr_amp_ad = 3*iqr(lfp_amp_pre(:));
%             thr_amp_ad = quantile(lfp_amp_pre(:),.75) + 1.5*iqr(lfp_amp_pre(:));
%             thr_amp_ad = 2*rms(lfp_amp_pre(:));
            
            minduration = 10*Fs; %min duration above threshold
            maxinterval = 5*Fs; %max interval between events to join
            
            %Gets continuous data above threshold or presaved indices
            lfp_ad = abs(lfp);
            for idx_train = 1:size(lfp,1)
                switch presaved
                    case 0
                        [idx_ad,idx_ad01] = findcontinuous(lfp_amp(idx_train,:),thr_amp_ad,minduration,maxinterval); %based on 'findpulses' to define minimum duration
                    case 1
                        idx_ad = Idx_ad{idx_train,1,idx_region,idx_day,idx_subj} + pretrain;
                        idx_ad01 = zeros(size(lfp,2),1);
                        for iad = 1:size(idx_ad,1)
                            idx_ad01(idx_ad(iad,1):idx_ad(iad,2))=1;
                        end
                end                   
                        
            %do not detect in pre-train
            if any(idx_ad<pretrain,'all') %if detected before train
                idx_ad_pre = idx_ad<pretrain;
                idx_ad(find(sum(idx_ad_pre,2)==2),:)=[]; %both pre
                idx_ad_pre = idx_ad<pretrain;
                idx_ad(idx_ad_pre(:,1),1) = pretrain+1; %turn before to first pos
                idx_ad01(1:pretrain)=false;
                clear idx_ad_pre
            end
                
            lfp_ad(idx_train,idx_ad01==0)=NaN; %hide non-AD
            
            %referentece to train end
            idx_ad = idx_ad - pretrain; 
            
            %Concatenate AD indices and LFP     
            Idx_ad{idx_train,1,idx_region,idx_day,idx_subj} = idx_ad;
            end
            
            %All trains figure
            if showfigs
            figure
            space = -2;
            hold on,h=plotmultisignals(t,abs(lfp),space); set(h,'color',[.8 .8 .8]) 
            hold on,h=plotmultisignals(t,lfp,space); set(h,'color','k')
%             hold on,h=plotmultisignals(t,lfp_amp,space); set(h,'color','b')    
            hold on,h=plotmultisignals(t,movmax(abs(lfp_ad),1*Fs,2),space); set(h,'color','r','linewidth',2);
            g0=gridxy([pretrain/Fs (pretrain+traindur)/Fs],[],'linestyle','--','color','k'); uistack(g0,'top')    
            title(['subject: ' num2str(idx_subj) ', day: ' num2str(idx_day) ', region: ' num2str(idx_region)])
            pause(0)
            end

            %---Perform manual fix of AD classification?---
            disp('Do you want to perform manual fixes of AD classification?');
            trials_interest = input('No: Enter; Yes: Inform the trials you want to fix (e.g: [1 2 4 7]): ');

            %Semi-autimatic detection by train
            for idx_train = trials_interest %1:size(lfp,1)
                
                %Automatic detection of ADs
                [idx_ad,idx_ad01] = findcontinuous(lfp_amp(idx_train,:),thr_amp_ad,minduration,maxinterval); %based on 'findpulses' to define minimum duration
                switch presaved
                    case 0
                        [idx_ad,idx_ad01] = findcontinuous(lfp_amp(idx_train,:),thr_amp_ad,minduration,maxinterval); %based on 'findpulses' to define minimum duration
                    case 1
                        idx_ad = Idx_ad{idx_subj,idx_day}{idx_train,1} + pretrain;
                        idx_ad01 = zeros(size(lfp,2),1);
                        for iad = 1:size(idx_ad,1)
                            idx_ad01(idx_ad(iad,1):idx_ad(iad,2))=1;
                        end
                end
                
                %do not detect in pre-train
                if any(idx_ad<pretrain,'all') %if detected before train
                    idx_ad_pre = idx_ad<pretrain;
                    idx_ad(find(sum(idx_ad_pre,2)==2),:)=[]; %both pre
                    idx_ad_pre = idx_ad<pretrain;
                    idx_ad(idx_ad_pre(:,1),1) = pretrain+1; %turn before to first pos
                    idx_ad01(1:pretrain)=false;
                    clear idx_ad_pre
                end
                
                %Single train figure
                figure,
                hold on,plot(t,abs(lfp(idx_train,:)),'color',[.8 .8 .8]);
                hold on,plot(t,lfp(idx_train,:),'k');
                hold on,plot(t,lfp_amp(idx_train,:),'b');
                for iad = 1:size(idx_ad,1)
    %             hold on,scatter(find(idx_ad01),lfp_amp(idx_train,idx_ad01),'filled','sb');
                hold on,plot(t(idx_ad(iad,1):idx_ad(iad,2)),movmax(abs(lfp(idx_train,idx_ad(iad,1):idx_ad(iad,2))),1*Fs,2),'r','linewidth',2);
                end
                axad(1)=gca;
                
                ylim([-1 1])
                
                g0=gridxy([pretrain/Fs (pretrain+traindur)/Fs],[],'linestyle','--','color','k'); uistack(g0,'top')
                g0=gridxy([],thr_amp_ad,'linestyle','--','color','b','linewidth',2); uistack(g0,'top')
    
                ylabel('mV')
                xlabel('Time (s)')
                
                title(['subject: ' num2str(idx_subj) ', day: ' num2str(idx_day) ', region: ' num2str(idx_region) ', train: ' num2str(idx_train)])
                
                %Signal Timestamps Marker App 
                signalmarkers=[]; open signalmarkers
                sigapp = signaltsmarker(lfp(idx_train,:),Fs); 
                axad(2)=gca;
                
                %Add to 'signaltsmarker' app
                idx_ad = transpose(idx_ad); idx_ad = idx_ad(:);
                sigapp.Markers{2,1} = idx_ad;
                sigapp.Markers{3,1} = idx_ad/Fs;
    
                %Add located trains to 'signaltsmarker' app for manual fix
                sigapp.plotSignal; 
    
                ylim([-1 1])
                linkaxes(axad,'xy')
                
                g0=gridxy([pretrain/Fs (pretrain+traindur)/Fs],[],'linestyle','--','color','k'); uistack(g0,'top')
                g0=gridxy([],thr_amp_ad,'linestyle','--','color','b','linewidth',2); uistack(g0,'top')
                
                uiwait(sigapp.SignalTimestampMarkerappUIFigure) %waits until close app
                signalmarkers = evalin('base','signalmarkers');
                
                close
                close
    
                %Get classified ADs
                idx_ad = signalmarkers{2,1};
                
                %concatenate in start-stop
                clear idx_ad_ss
                idx_ad_ss(:,1) = idx_ad(1:2:end); 
                idx_ad_ss(:,2) = idx_ad(2:2:end); 
                idx_ad = idx_ad_ss; 
    
                %fix if detected in pre-train
                if any(idx_ad<pretrain,'all') %if detected before train
                    idx_ad_pre = idx_ad<pretrain;
                    idx_ad(find(sum(idx_ad_pre,2)==2),:)=[]; %both pre
                    idx_ad_pre = idx_ad<pretrain;
                    idx_ad(idx_ad_pre(:,1),1) = pretrain+1; %turn before to first pos
                    idx_ad01(1:pretrain)=false;
                    clear idx_ad_pre
                end
                
                %referentece to train end
                idx_ad = idx_ad - pretrain; 

                %Concatenate fixed AD indices and LFP     
%                 Lfp_amp(idx_train,:,idx_region,idx_day,idx_subj) = lfp_amp; %does not need AD indices
%                 Lfp_ad(idx_train,:,idx_region,idx_day,idx_subj) = lfp_ad; %needs fixing!

%                 Idx_ad{idx_subj,idx_day}{idx_train,1} = idx_ad;
                Idx_ad{idx_train,1,idx_region,idx_day,idx_subj} = idx_ad;
                          
            end %train
        
        end %region
    end %day
end %subj



