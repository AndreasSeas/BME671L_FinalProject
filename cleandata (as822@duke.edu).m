close all
clear
clc

%% make directory set
prgdir=cd;
cd ..
homedir=cd;
cd("RAW DATA");
datadir=cd;

%% load dataset
cd(datadir)
raw=edfread("PN00-1.edf");
cd(prgdir)

%% process

[numtimes,numelem]=size(raw);
fs=numel(raw{1,1}{1});
t=[0:1/fs:numtimes(end)-1/fs]';
vars=raw.Properties.VariableNames';
rawmat=cell2mat(raw{:,:});

%% plot
f=figure; hold on;
tiledlayout(29,1,'TileSpacing','none');
pf=[];
for i=1:29
    nexttile;
    plot(t,rawmat(:,i)); grid on;
    ylabel(vars{i});
    set(gca,'visible','off');
    
    
    y0 = fftshift(rawmat(:,i));         % shift y values
    n=length(t);
    f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
    power0 = abs(y0).^2/n;    % 0-centered power
    pf=[pf;power0'];
    
    %     
%     plot(f0,power0)
%     xlabel('Frequency')
%     ylabel('Power')
    
end

%% ***


%% etc

% figure;
% temp=downsample(rawmat,10000);
% surf(temp');view(2);
% T=timetable2table(raw);
% cell2mat(raw.EEGC3)
return
%% temp

eeg=[];
for C=1:numelem
%     tempstream=[];
%     for R=1:numtimes
%         tempstream=[tempstream;raw{R,C}{1}];
%     end
    
    disp(num2str(R)+","+num2str(C))
    
    eeg=[eeg, tempstream];
end
