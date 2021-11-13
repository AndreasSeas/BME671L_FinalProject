%% clean slate
close all;clear;clc;


%% save data? save plots?
savedata=0;
makeplots=1;
saveplots=0;

%% load in the dataset of interest
edfname="PN00-1.edf";
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

%% attempt filtration/ID of initiation of SZ
% do combo of:
%	[X] power analysis, change in freq content to hi
%	[X] correlation b/w channels
%	[X] derivative/movvar with time
%   [ ] measures of synchrony
cheatrange=[1100,1250];
cheatidx=[find(t==cheatrange(1)),find(t==cheatrange(2))];
eeg_ready=eeg_notch(cheatidx(1):cheatidx(2),:);
t_ready=t(cheatidx(1):cheatidx(2));

% https://www.sciencedirect.com/science/article/pii/S0010482518304037
eeg1=eeg_ready(:,1);
eeg2=eeg_ready(:,2);

% c=[];
% for i=1:numel(eeg1)-10
% 	c=[c;corr(eeg1(i:i+10),eeg2(i:i+10))];
% end

corr(eeg1,eeg2)



% c=dt*conv(eeg1,fliplr(eeg2));

%spectrogram
% pspectrum(eeg_ready(:,1),fs,'spectrogram');


%plot(t_ready,eeg_ready)

% plot(t_ready,eeg_ready)
% hold on; 
% plot([1143,1143],[min(min(eeg_ready)),max(max(eeg_ready))])
% plot([1213,1213],[min(min(eeg_ready)),max(max(eeg_ready))])

% %% plot reference electrodes

% if makeplots==1
% 	f=figure;hold on; grid on;
% 	plot(t,eeg_hp(:,1),'-b');
% 	plot(t,eeg(:,1),'-r');
% end
