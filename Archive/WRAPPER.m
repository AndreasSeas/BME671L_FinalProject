%% clean slate
close all;clear;clc;


%% load in the dataset of interest
% ptsznum=1;
% edfname="PN05-2.edf";
% txtname="Seizures-list-PN05.txt";

ptsznum=1;
edfname="PN00-1.edf";
txtname="Seizures-list-PN00.txt";
cdrawdata='/Users/as822/Box/!PhD/Courses/BME671L/BME671L_FinalProject/BME671L_FinalProject/RAW DATA';

disp("loading and processing " + edfname)
[t,rawmat,colnames]=edf2eeg(cdrawdata,edfname);
dt=t(2)-t(1);
Fs=1/dt;
 
eeg_raw=rawmat(:,1:29);% change eventually to be "intelligent"

%% load in the 2d organization data (chanmap)
chanmap=readtable("net_10_20_xy.xlsx");

%% get indices of sz start and end
[Start_Seiz, Dur_Seiz]=seizureTime(cdrawdata,txtname);


startidx=find(t==Start_Seiz(ptsznum));
endidx=find(t==Start_Seiz(ptsznum)+Dur_Seiz(ptsznum));

%% bandpower - all leads
freqlo=[1,4,8,12,30];
freqhi=[4,8,12,30,100];
freqnames={"delta","theta","alpha","beta","gamma"};
L=length(eeg_raw);
ilo=[1,startidx,endidx];
ihi=[startidx,endidx,L];
inames={"before","during","after"};
% 
% aF=[];% freq
% aE=[];% epoch
% % aL=[];% lead
% aLobe=[];% lobe
% aSide=[];% side
% Pval=[];
% for i=1:numel(chanmap.idx)
%     lead=chanmap.idx(i);
%     
% %     Phold=zeros(numel(freqlo),numel(ilo));
%     
%     for m=1:numel(freqlo)
%         
%         for n=1:numel(ilo)
%             
%             Ptemp=bandpower(eeg_raw(ilo(n):ihi(n),lead),Fs,[freqlo(m),freqhi(m)]);
% %             disp("done lead "+ num2str(lead) + ", " + freqnames{m} + ", " + inames{n});
%             %         disp(["done "  freqnames{m} ", "  inames{n}],'Interpreter','latex');
%             
%             aF=[aF;freqnames{m}];
%             aE=[aE;inames{n}];
%             aLobe=[aLobe;chanmap.Lobe(i)];
%             aSide=[aSide;chanmap.Side(i)];
%             Pval=[Pval;Ptemp];
%             
%         end
%         
%     end
% end

% 
% %% anova to prove that we aren't *totally* insane and stats protected
% [p,tbl,stats]=anovan(Pval,{aF,aE,aLobe,aSide});
% % [p,tbl,stats]=anovan(Pval,{aF,aE,aLobe,aSide},'model','interaction');
% %%
% results = multcompare(stats,'Dimension',[2,4]);

%% for each channel, get bandpower for each second

window_dt=1;% second... make sure inverse is always divisible by 2 if !=1
window_di=window_dt.*Fs;
i_band=1:window_di:numel(t);
t_band=t(i_band);
% t_start=t_band(1:end-1);
% t_end=t_band(2:end);
t_center=t_band(1:end-1)+window_dt/2;
BPwT=zeros(numel(t_band)-1,numel(chanmap.idx),numel(freqlo));

for i_freq=1:numel(freqlo)
   
    for i_chan=1:numel(chanmap.idx)
       disp("frequency = " + freqnames{i_freq} + " and chan = " + chanmap.Sensor{i_chan})
       tic
        for i_t=1:numel(t_band)-1
            
            BPwT(i_t,i_chan,i_freq)=bandpower(...
                eeg_raw(i_band(i_t):i_band(i_t+1),chanmap.idx(i_chan))...
                ,Fs,[freqlo(i_freq),freqhi(i_freq)]);
            
        end
        tempt=toc;
        disp("      took " + num2str(tempt) + " seconds")
    end
    
end

%% plot average bandpower over brain with time
freqstart=2;
avgpr=mean(BPwT(:,:,freqstart:end),2);
[r,~,c]=size(avgpr);
f=figure; hold on; grid on;
maxval=max(max(avgpr));
fill([Start_Seiz(ptsznum),Start_Seiz(ptsznum)+Dur_Seiz(ptsznum),...
    Start_Seiz(ptsznum)+Dur_Seiz(ptsznum),Start_Seiz(ptsznum)],...
    [0,0,maxval,maxval],'r','facealpha',0.3)
plot(t_center,reshape(avgpr,r,c))
legend(['sztime', freqnames(freqstart:end)])
xlim([Start_Seiz(ptsznum)-Dur_Seiz(ptsznum),Start_Seiz(ptsznum)+2*Dur_Seiz(ptsznum)])

%% new highpass
delidx=endidx-startidx;
eeg_raw2=eeg_raw(startidx-10*delidx:startidx+11*delidx,chanmap.idx);
tval2=t(startidx-10*delidx:startidx+11*delidx);
V=highpass(eeg_raw2',4,Fs)';

nustart=(startidx-(startidx-10*delidx));
nuend=nustart+delidx;

[n,nchan]=size(V);
V2=V;
mV=mean(V,2);
for i=1:nchan
   
    V2(:,i)=V2(:,i)-mV;
    
end

%%
% soundsc

return
%% topogram plot

ch_list = chanmap.Sensor;
ch_list(8)=[];
ch_list(end)=[];
values=V2;
values(:,8)=[];
values(:,end)=[];
valminmax=values(nustart-Fs*3:nuend+Fs*3,:);
cbds=[min(min(valminmax)),max(max(valminmax))];
plot(values)

f=figure; hold on;
daspect([1,1,1]);
timestr=1;
for i=nustart-Fs*3:512/16:nuend+Fs*3%:length(values)
    h=plot_topography(ch_list,values(i,:),false,'10-20',true,false,100);
    colormap("hot");
    caxis(cbds);
    
    
    if i<nustart 
        epoch="before";
    elseif i>nuend
        epoch="after";
    else
        epoch="seizure";
    end
        
        
    title("t = " + timestr + " seconds... epoch =  " + epoch);
    timestr=timestr+1/16;
    drawnow
    
%     pause(0.001);
    
end

%% topographic representation of power for each band

bandid=2;
ch_list = chanmap.Sensor;
ch_list(8)=[];
ch_list(end)=[];
values=BPwT(:,:,bandid);
maxvaltot=max(max(log(values)));
minvaltot=min(min(log(values)));


f=figure; hold on;
for i=1100:1500%length(values)
    values=BPwT(:,:,bandid);
    values(:,8)=[];
    values(:,end)=[];
    h=plot_topography(ch_list,log(values(i,:)),false,'10-20',true,false,100);
    colormap("hot");
    caxis([minvaltot,maxvaltot]);
    
    
    if i>1143 && i<1213
        namemeyo='seizure';
    else
        namemeyo='not';
    end
    
    title(namemeyo);
    drawnow
        
%     pause(0.01);
    
end




% Víctor Martínez-Cagigal (2021). Topographic EEG/MEG plot
% (https://www.mathworks.com/matlabcentral/fileexchange/72729-topographic-eeg-meg-plot),
% MATLAB Central File Exchange. Retrieved December 4, 2021.


% [X,Y] = meshgrid(-1:0.025:1, -1:0.025:1);
% xq = X(:);
% yq = Y(:);
% rTest = sqrt(xq.^2 + yq.^2);
% ckUnit = rTest <= 1;
% xq = xq(ckUnit);
% yq = yq(ckUnit);
% nInt = length(xq);
% 
% [ntimes,nchan]=size(meanvolts);
% 
% EEG_int = zeros(5025, ntimes); 
% % interpolate at each time time point per case
% for jj = 1:ntimes
%     EEG_int(:,jj) = griddata(xy(:,1), xy(:,2), ...
%         meanvolts(1:end-1,jj), xq, yq);
% end


return





%% could also potentially do this to find start and end unsupervised

% results = multcompare(stats,'Dimension',[1 2,3])

% %% bandpower - lead 1 only
% 
% lead=1;
% freqlo=[0.5,4,8,12,30];
% freqhi=[4,8,12,30,100];
% freqnames={"\delta","\theta","\alpha","\beta","\gamma"};
% L=length(eeg_raw);
% ilo=[1,startidx,endidx];
% ihi=[startidx,endidx,L];
% inames={"before","during","after"};
% 
% Phold=zeros(numel(freqlo),numel(ilo));
% 
% for m=1:numel(freqlo)
%    
%     for n=1:numel(ilo)
%        
%         Phold(m,n)=bandpower(eeg_raw(ilo(n):ihi(n),lead),Fs,[freqlo(m),freqhi(m)]);
%         disp("done " + freqnames{m} + ", " + inames{n});
% %         disp(["done "  freqnames{m} ", "  inames{n}],'Interpreter','latex');
%         
%     end
%     
% end

% BP_D=bandpower(eeg_raw(1:startidx,lead),Fs,[1,4]);


% y_D=bandpass(eeg_raw(:,lead),[1,4],Fs);

%% 
figure; hold on;
plot(eeg_raw(:,lead))
plot(y_D)
return

%% FFT
lead=1;
L=length(eeg_raw);
Y=fft(eeg_raw(:,lead));
P2=abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,(P1));
xlim([1,100])

return
%% only filter into the bands we want
% eeg_filt

%% filter in outside function
eeg_filt=filterEEG(eeg_raw,Fs);

%% design custom filter to remove baseline noise

fvtool(d)

% rather than a notch, assume the data is free of 60Hz noise to avoid
% ringing

%% get only specific bands of signal
Y=fft(eeg_raw(1,:));
    
return
%% highpass filter raw data

% in order to minimize edge artifacts, we can do some cool manipulation of
% the data to create an artificial edge
padlength=10000;% add these points to edge
eeg_pad=[flipud(eeg_raw(1:padlength,:));...
    eeg_raw;
    flipud(eeg_raw(length(eeg_raw)-padlength:end,:))]';
% eeg_hp=eeg_pad;
% for i=1:29
%     eeg_hp(i,:)=filter(d,eeg_pad(i,:))';
% end
% 
eeg_hp=filter(d,eeg_pad)';
eeg_hp=eeg_hp(padlength+1:end-padlength,:);

f=figure; 
subplot(2,1,1);
plot(eeg_raw);
subplot(2,1,2);
plot(eeg_hp);

return
%% high pass filter to remove moving baseline
disp("hpf of dataset");
eeg_hp = highpass(eeg',0.4,Fs)';

%% notch filter to remove 60 Hz noise
disp("notch dataset");
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);
% fvtool(d,'fs',fs)
eeg_notch = filtfilt(d,eeg_hp')';

%save
% temp={eeg_notch,t,fs,rawmat,colnames};
% save("temp.mat","temp");

%% rereference
iCz=15;
fCz = -sum(eeg_notch(:,~15))/29;%

% eeg_rr=

%% sumpower for all
sumpower=sum((eeg_notch.^2),2);
figure;hold on;grid on;
% Sstart=1143;Send=1213;
Sstart=1220;Send=1274;
fill([Sstart,Send,Send,Sstart],[0,0,1,1],'r','facealpha',0.5);
plot(t,sumpower./max(sumpower));




% Seizure n 1
% File name: PN00-1.edf
% Registration start time: 19.39.33
% Registration end time:  20.22.58
% Seizure start time: 19.58.36
% Seizure end time: 19.59.46

%% sum of power
% cheatrange=[1100,1250];
% cheatidx=[find(t==cheatrange(1)),find(t==cheatrange(2))];
% eeg_ready=eeg_notch(cheatidx(1):cheatidx(2),:);
% t_ready=t(cheatidx(1):cheatidx(2));
% 
% sumpower=sum((eeg_ready.^2),2);
% plot(t_ready,sumpower);

% %% trying shit with windows and correl
% cheatrange=[1100,1250];
% cheatidx=[find(t==cheatrange(1)),find(t==cheatrange(2))];
% eeg_ready=eeg_notch(cheatidx(1):cheatidx(2),:);
% t_ready=t(cheatidx(1):cheatidx(2));
% 
% % plot(t_ready,eeg_ready);
% windowlength=fs;
% numberwindows=floor(numel(t_ready)/windowlength);%
% i_start=1:windowlength:windowlength*(numberwindows-1);
% i_end=i_start+(windowlength-1);
% 
% corrmat=[];
% tcentr=[];
% for i=1:numel(i_start)
%     y=eeg_ready(i_start(i):i_end(i),:);
%     x=t_ready(i_start(i):i_end(i));
%     ycor=corr(y);
%     ycor
% %     corrmat=[corrmat;mean(mean(corr(y)))];
%     tcentr=[tcentr;(max(x)+min(x))/2];
% end
% 
% figure;hold on;grid on;
% plot(tcentr,corrmat);


% %% attempt filtration/ID of initiation of SZ
% % do combo of:
% %	[X] power analysis, change in freq content to hi
% %	[X] correlation b/w channels
% %	[X] derivative/movvar with time
% %   [ ] measures of synchrony
% cheatrange=[1100,1250];
% cheatidx=[find(t==cheatrange(1)),find(t==cheatrange(2))];
% eeg_ready=eeg_notch(cheatidx(1):cheatidx(2),:);
% t_ready=t(cheatidx(1):cheatidx(2));
% 
% % https://www.sciencedirect.com/science/article/pii/S0010482518304037
% eeg1=eeg_ready(:,1);
% eeg2=eeg_ready(:,2);
% 
% % c=[];
% % for i=1:numel(eeg1)-10
% % 	c=[c;corr(eeg1(i:i+10),eeg2(i:i+10))];
% % end
% 
% % corr(eeg1,eeg2)
% 
% 
% 
% % c=dt*conv(eeg1,fliplr(eeg2));
% 
% %spectrogram
% pspectrum(eeg_ready(:,1),fs,'spectrogram');
% 
% 
% %plot(t_ready,eeg_ready)
% 
% % plot(t_ready,eeg_ready)
% % hold on; 
% % plot([1143,1143],[min(min(eeg_ready)),max(max(eeg_ready))])
% % plot([1213,1213],[min(min(eeg_ready)),max(max(eeg_ready))])
% 
% % %% plot reference electrodes
% 
% % if makeplots==1
% % 	f=figure;hold on; grid on;
% % 	plot(t,eeg_hp(:,1),'-b');
% % 	plot(t,eeg(:,1),'-r');
% % end
% 
% 
% %% topomaker - for later
% 
