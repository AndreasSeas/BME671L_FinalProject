%% clean slate
close all;clear;clc;


%% save data? save plots?
savedata=0;
makeplots=1;
saveplots=0;

%% load in the dataset of interest
edfname="PN00-2.edf";
directorystring='/Users/as822/Box/!PhD/Courses/BME671L/BME671L_FinalProject/RAW DATA';

disp("loading and processing " + edfname)
[t,rawmat,colnames]=edf2eeg(directorystring,edfname);
dt=t(2)-t(1);
fs=1/dt;
 
eeg=rawmat(:,1:29);% change eventually to be "intelligent"

%% high pass filter to remove moving baseline
disp("hpf of dataset");
eeg_hp = highpass(eeg',0.4,fs)';

%% notch filter to remove 60 Hz noise
disp("notch dataset");
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);
% fvtool(d,'fs',fs)
eeg_notch = filtfilt(d,eeg_hp')';

%save
% temp={eeg_notch,t,fs,rawmat,colnames};
% save("temp.mat","temp");

%% rereference
iCz=15;
fCz = -sum(eeg_notch(:,~15))/29;%

eeg_rr=

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
