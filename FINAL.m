% final code for BME 671L Project 
% Yining, Emily, Andreas 
% Saturday, December 4th

%% clean slate
close all;clear;clc;

%% set up directories
%establish data directory
dir_home=cd;
cd("RAW DATA");
dir_raw=cd;
cd(dir_home);
cd("OUTPUTS");
dir_out=cd;
cd(dir_home);

%% load in the dataset of interest
ptsznum=1; % this patient had several seizures recorded, we want the first
edfname="PN00-1.edf";% name of the .edf file
txtname="Seizures-list-PN00.txt";% name of timing/info file
disp("loading and processing " + edfname)% display what we are working on
[t,rawmat,colnames]=edf2eeg(dir_raw,edfname); % use our edf2eeg to load eeg
dt=t(2)-t(1);% get the dt independently
Fs=1/dt;% get the sampling frequency%

eeg_raw=rawmat(:,1:29);% get the raw EEG data from the raw matrix

%% load in the 2d organization data (chanmap)
disp("loading channel map for 10-20 array");
chanmap=readtable("net_10_20_xy.xlsx");% load the 10-20 EEG array locations

f=figure; hold on; grid on;
fill([-1,0,1],[0,1.6,0],[0.9,0.9,0.9])
theta=linspace(0,2*pi,1000);
r=1.3;
x=r*cos(theta);y=r*sin(theta);
fill(x,y,[0.9,0.9,0.9])
f.Position=[440 160 768 638];
scatter(chanmap.X,chanmap.Y,1900,[1,0.4,0.4],'filled')
daspect([1,1,1])
a=gca;a.Visible='off';
text(chanmap.X,chanmap.Y,chanmap.Sensor,'FontSize',20,"HorizontalAlignment","center")
cd(dir_out);
exportgraphics(f,'10-20array.pdf','Resolution',300);
cd(dir_home);

%% get indices of sz start and end
[Start_Seiz, Dur_Seiz]=seizureTime(dir_raw,txtname);% extract the seizure
% times for this specific edf
startidx=find(t==Start_Seiz(ptsznum));% put this in terms of index
endidx=find(t==Start_Seiz(ptsznum)+Dur_Seiz(ptsznum));

%% highpass filtration


% we need to do high-pass filtration here to better understand interpret
% our signal

% however, there is a TON of data... so we want to focus on the area around
% the seizure itself (about 10 min before and 10 min after

delidx=endidx-startidx;% get length of sz in index
% now define new eeg_raw2 that has only region of interest
eeg_raw2=eeg_raw(startidx-60*Fs:endidx+60*Fs,:);% get new eeg
tval2=t(startidx-60*Fs:endidx+60*Fs);% get corresponding tvals

% use in-build highpass filtration
V=highpass(eeg_raw2,4,Fs);
% note that we are high-pass filtering above 4 Hz. Note that this is by no
% means perfect.... there is important data between 1-4 Hz called the delta
% band. However, we are leaving that out here, mainly because this band
% also comes with a TON of noise from movement, eye blinking, and being
% human.... When we go and look at our bandpower, you'll see how delta has
% a TON of power (no not COVID delta)

nustart=(startidx-(startidx-60*Fs));% get a new starting index for the SZ
nuend=nustart+delidx;

%% rereference the data
% the scalp is not very well grounded usually... and no we dont mean in
% reality..... but our scalps have a ton of problems. There is sweat, hair,
% and a ton of other *crap* that makes it so charge can accumulate or sink
% in weird ways. Therefore we need to "re-reference" all of our electrodes
% in order to help get the most accurate and realistic data out of each one

% note that, while 10-20 only has 19 electrodes, we are keeping all 29
% original electrodes to help wiuth re-referencing... the more data on the
% conductivity of the scalp, the better

[~,nchan]=size(V);% get the size of V
V2=V;% save a new copy of Voltages to rereference
meanV=mean(V,2);%
for i=1:nchan
    
    V2(:,i)=V2(:,i)-meanV;
    
end

V2_mem=V2;% just saving a version here for fun

% lets now cut down to only the voltages we want for the 10-20 array
V2=V2(:,chanmap.idx);%

%% plot the different steps of processing
f=figure;
subplot(3,1,1);grid on; hold on;
maxval=max(max(eeg_raw2));minval=min(min(eeg_raw2));
fill([tval2(nustart),tval2(nuend),tval2(nuend),tval2(nustart)],...
    [minval,minval,maxval,maxval],'r','facealpha',0.3)
plot(tval2,eeg_raw2);
xlabel('time (seconds)'); ylabel('voltage (\muV)'); title('raw data')
set(gca,'FontSize',20)
subplot(3,1,2);grid on; hold on;
maxval=max(max(V));minval=min(min(V));
fill([tval2(nustart),tval2(nuend),tval2(nuend),tval2(nustart)],...
    [minval,minval,maxval,maxval],'r','facealpha',0.3)
plot(tval2,V);
xlabel('time (seconds)'); ylabel('voltage (\muV)'); title('high-pass filtered')
set(gca,'FontSize',20)
subplot(3,1,3);grid on; hold on;
maxval=max(max(V2_mem));minval=min(min(V2_mem));
fill([tval2(nustart),tval2(nuend),tval2(nuend),tval2(nustart)],...
    [minval,minval,maxval,maxval],'r','facealpha',0.3)
plot(tval2,V2_mem);
xlabel('time (seconds)'); ylabel('voltage (\muV)'); title('re-referenced and filtered')
set(gca,'FontSize',20)

f.Position=[26 1 1415 804];

cd(dir_out);
exportgraphics(f,'raw_EEG_processing.pdf','Resolution',300);
cd(dir_home);

%% visualize voltage topogram

ch_list = chanmap.Sensor;
ch_list(8)=[];ch_list(end)=[];% remove electrodes not in 10-20
values=V2;values(:,8)=[];values(:,end)=[];% remove electrodes not in 10-20

rng=mean(std(values))*5;
cbds=[-rng,rng];
% cbds=[min(min(values)),max(max(values))];

f=figure('Visible','on','Position',[26 1 943 804]); hold on;
daspect([1,1,1]);
subplot(4,1,1);hold on; grid on;
fill([tval2(nustart),tval2(nuend),tval2(nuend),tval2(nustart)],...
    3*[cbds(1),cbds(1),cbds(2),cbds(2)],'r','facealpha',0.3)
plot(tval2,values);
i=1;
tempplot=plot([tval2(i),tval2(i)],cbds,'--k','linewidth',1);

for i=1:Fs:length(values)
    
    subplot(4,1,1);hold on; grid on;
    delete(tempplot)
    tempplot=plot([tval2(i),tval2(i)],3*cbds,'-k','linewidth',1);
    hold off
    subplot(4,1,2:4);
    h=plot_topography(ch_list,values(i,:),false,'10-20',true,false,100);
    % this is an already-existing function we got from mathworks:
    %     Víctor Martínez-Cagigal (2021). Topographic EEG/MEG plot
    %     (https://www.mathworks.com/matlabcentral/fileexchange/72729-topographic-eeg-meg-plot),
    %     MATLAB Central File Exchange. Retrieved December 5, 2021.
    colormap("hot");
    caxis(cbds);
    hold off
    drawnow
end

%% Now, lets look at the band power for fun - 1 second window
% that way we can look at bandpower and how it changes with time

freqlo=4;%[1,4,8,12,30];% identify the start frequency of major bands
freqhi=8;%[4,8,12,30,100];% end freq of bands
freqnames={"theta"};%{"delta","theta","alpha","beta","gamma"};% names of bands
L=length(eeg_raw);% length of EEG data
ilo=[1,startidx,endidx];% start index of each epoch (before/during/after)
ihi=[startidx,endidx,L];% end index
inames={"before","during","after"};% name of each epoch

window_dt=1;% second... make sure inverse is always divisible by 2 if !=1
window_di=window_dt.*Fs;
i_band=1:window_di:numel(t);
t_band=t(i_band);
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

%% visualize bandpower in theta band (already associated with sz)
values=BPwT(1080:1280,:);values(:,8)=[];values(:,end)=[];% remove electrodes not in 10-20
cbds=[0,max(max(values))];

f=figure('Visible','on','Position',[26 1 943 804]); hold on;
daspect([1,1,1]);
subplot(4,1,1);hold on; grid on;
fill([tval2(nustart),tval2(nuend),...
    tval2(nuend),tval2(nustart)],...
    [cbds(1),cbds(1),cbds(2),cbds(2)],'r','facealpha',0.3)
plot(t_center(1080:1280),values(:,:));
xlim([1080,1280]);
tfornow=t_center(1080:1280);
tempplot=plot([tval2(i),tval2(i)],cbds,'--k','linewidth',1);

for i=1:length(values)% no Fs since already in seconds
    
    subplot(4,1,1);hold on; grid on;
    delete(tempplot)
    tempplot=plot(tfornow(i)*ones(1,2),cbds,'--k','linewidth',1);
    hold off
    subplot(4,1,2:4);
    h=plot_topography(ch_list,values(i,:),false,'10-20',true,false,100);
    % this is an already-existing function we got from mathworks:
    %     Víctor Martínez-Cagigal (2021). Topographic EEG/MEG plot
    %     (https://www.mathworks.com/matlabcentral/fileexchange/72729-topographic-eeg-meg-plot),
    %     MATLAB Central File Exchange. Retrieved December 5, 2021.
    colormap("hot");
    caxis(cbds);
    hold off
    drawnow
end

%% functions
function [t,rawmat,colnames]=edf2eeg(directorystring,edfname)
prgdir=pwd;
cd(directorystring)
raw=edfread(edfname);
cd(prgdir)

%% process
[numtimes,~]=size(raw);
fs=numel(raw{1,1}{1});
t=[0:1/fs:numtimes(end)-1/fs]';
colnames=raw.Properties.VariableNames';
rawmat=cell2mat(raw{:,:});

end

function [Start_Seiz, Dur_Seiz]=seizureTime(dirname,filename)
homedir=cd;
cd(dirname);
list = fopen(filename);

tline = fgetl(list);

lfile = [];
rV = [];
sSV = [];
sEV = [];

while ischar(tline)
    if strfind(tline, 'File name:') ~= 0
        lfile = [lfile ; tline];
    else if strfind(tline, 'Registration') ~= 0
            tline = tline(find(~isspace(tline)));
            deleteC = isletter(tline);
            tline(deleteC) = [];
            tline = erase(tline, ':');
            rV = [rV ; tline];
        else if strfind(tline, 'Seizure start') ~= 0
                tline = tline(find(~isspace(tline)));
                deleteC = isletter(tline);
                tline(deleteC) = [];
                tline = erase(tline, ':');
                sSV = [sSV ; tline];
            else if strfind(tline, 'Seizure end') ~= 0
                    tline = tline(find(~isspace(tline)));
                    deleteC = isletter(tline);
                    tline(deleteC) = [];
                    tline = erase(tline, ':');
                    sEV = [sEV ; tline];
                end
            end
        end
    end
    tline = fgetl(list);
end

[krow, kcolumn] = size(rV);
regS_secV =[];
sSV_secV = [];
sES_secV = [];

for step = 1:krow
    if mod(step, 2) ~= 0
        regS = rV(step, :);
        regS_sec = str2num(regS(1:2))*60*60 + str2num(regS(4:5))*60 + str2num(regS(7:8));
        regS_secV = [regS_secV regS_sec];
    end
    if step < krow/2+1
        sSS = sSV(step, :);
        sSS_sec = str2num(sSS(1:2))*60*60 + str2num(sSS(4:5))*60 + str2num(sSS(7:8));
        sSV_secV = [sSV_secV sSS_sec];
        
        sES = sEV(step, :);
        sES_sec = str2num(sES(1:2))*60*60 + str2num(sES(4:5))*60 + str2num(sES(7:8));
        sES_secV = [sES_secV sES_sec];
    end
end

% display(lfile)
Start_Seiz =sSV_secV - regS_secV;  %in seconds
Dur_Seiz = sES_secV - sSV_secV;    %in seconds

cd(homedir)
end


% we have a lot more stuff to show, because EEG is complicated.... but
% thats all for now folks


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % BELOW ARE JUST SCRAPS N FUN IDEAS N STATS % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% % % 
% % % % below are some stats for SnG
% % % %% we also want to look at the band power of different types of rhythms
% % % freqlo=[1,4,8,12,30];% identify the start frequency of major bands
% % % freqhi=[4,8,12,30,100];% end freq of bands
% % % freqnames={"delta","theta","alpha","beta","gamma"};% names of bands
% % % L=length(eeg_raw);% length of EEG data
% % % ilo=[1,startidx,endidx];% start index of each epoch (before/during/after)
% % % ihi=[startidx,endidx,L];% end index
% % % inames={"before","during","after"};% name of each epoch
% % % 
% % % aF=[];% freq aE=[];% epoch % aL=[];% lead aLobe=[];% lobe aSide=[];% side
% % % Pval=[]; 
% % % for i=1:numel(chanmap.idx)
% % %     lead=chanmap.idx(i);
% % %     
% % %     %     Phold=zeros(numel(freqlo),numel(ilo));
% % %     
% % %     for m=1:numel(freqlo)
% % %         
% % %         for n=1:numel(ilo)
% % %             
% % %             Ptemp=bandpower(eeg_raw(ilo(n):ihi(n),lead),Fs,[freqlo(m),freqhi(m)]);
% % %             %             disp("done lead "+ num2str(lead) + ", " + freqnames{m} + "," + inames{n});
% % %             %         disp(["done "  freqnames{m} ", "inames{n}],'Interpreter','latex');
% % %             
% % %             aF=[aF;freqnames{m}]; aE=[aE;inames{n}];
% % %             aLobe=[aLobe;chanmap.Lobe(i)]; aSide=[aSide;chanmap.Side(i)];
% % %             Pval=[Pval;Ptemp];
% % %             
% % %         end
% % %         
% % %     end
% % % end
% % % 
% % % 
% % % %% anova to prove that we aren't *totally* insane and stats protected
% % % [p,tbl,stats]=anovan(Pval,{aF,aE,aLobe,aSide}); %
% % % [p,tbl,stats]=anovan(Pval,{aF,aE,aLobe,aSide},'model','interaction'); %%
% % % results = multcompare(stats,'Dimension',[2,4]);
% % % 
% % % return
% % % %% bandpower - all leads
% % % 
% % % 
% % % %% for each channel, get bandpower for each second
% % % 
% % % window_dt=1;% second... make sure inverse is always divisible by 2 if !=1
% % % window_di=window_dt.*Fs;
% % % i_band=1:window_di:numel(t);
% % % t_band=t(i_band);
% % % % t_start=t_band(1:end-1); t_end=t_band(2:end);
% % % t_center=t_band(1:end-1)+window_dt/2;
% % % BPwT=zeros(numel(t_band)-1,numel(chanmap.idx),numel(freqlo));
% % % 
% % % for i_freq=1:numel(freqlo)
% % %     
% % %     for i_chan=1:numel(chanmap.idx)
% % %         disp("frequency = " + freqnames{i_freq} + " and chan = " + chanmap.Sensor{i_chan})
% % %         tic
% % %         for i_t=1:numel(t_band)-1
% % %             
% % %             BPwT(i_t,i_chan,i_freq)=bandpower(...
% % %                 eeg_raw(i_band(i_t):i_band(i_t+1),chanmap.idx(i_chan))...
% % %                 ,Fs,[freqlo(i_freq),freqhi(i_freq)]);
% % %             
% % %         end
% % %         tempt=toc;
% % %         disp("      took " + num2str(tempt) + " seconds")
% % %     end
% % %     
% % % end
% % % 
% % % %% plot average bandpower over brain with time
% % % freqstart=2;
% % % avgpr=mean(BPwT(:,:,freqstart:end),2);
% % % [r,~,c]=size(avgpr);
% % % f=figure; hold on; grid on;
% % % maxval=max(max(avgpr));
% % % fill([Start_Seiz(ptsznum),Start_Seiz(ptsznum)+Dur_Seiz(ptsznum),...
% % %     Start_Seiz(ptsznum)+Dur_Seiz(ptsznum),Start_Seiz(ptsznum)],...
% % %     [0,0,maxval,maxval],'r','facealpha',0.3)
% % % plot(t_center,reshape(avgpr,r,c))
% % % legend(['sztime', freqnames(freqstart:end)])
% % % xlim([Start_Seiz(ptsznum)-Dur_Seiz(ptsznum),Start_Seiz(ptsznum)+2*Dur_Seiz(ptsznum)])
% % % 
% % % 
% % % %%
% % % % soundsc
% % % 
% % % return
% % % %% topogram plot
% % % 
% % % ch_list = chanmap.Sensor;
% % % ch_list(8)=[];
% % % ch_list(end)=[];
% % % values=V2;
% % % values(:,8)=[];
% % % values(:,end)=[];
% % % valminmax=values(nustart-Fs*3:nuend+Fs*3,:);
% % % cbds=[min(min(valminmax)),max(max(valminmax))];
% % % plot(values)
% % % 
% % % f=figure; hold on;
% % % daspect([1,1,1]);
% % % timestr=1;
% % % for i=nustart-Fs*3:512/16:nuend+Fs*3%:length(values)
% % %     h=plot_topography(ch_list,values(i,:),false,'10-20',true,false,100);
% % %     colormap("hot");
% % %     caxis(cbds);
% % %     
% % %     
% % %     if i<nustart
% % %         epoch="before";
% % %     elseif i>nuend
% % %         epoch="after";
% % %     else
% % %         epoch="seizure";
% % %     end
% % %     
% % %     
% % %     title("t = " + timestr + " seconds... epoch =  " + epoch);
% % %     timestr=timestr+1/16;
% % %     drawnow
% % %     
% % % end
% % % 
% % % %% topographic representation of power for each band
% % % 
% % % bandid=2;
% % % ch_list = chanmap.Sensor;
% % % ch_list(8)=[];
% % % ch_list(end)=[];
% % % values=BPwT(:,:,bandid);
% % % maxvaltot=max(max(log(values)));
% % % minvaltot=min(min(log(values)));
% % % 
% % % 
% % % f=figure; hold on;
% % % for i=1100:1500%length(values)
% % %     values=BPwT(:,:,bandid);
% % %     values(:,8)=[];
% % %     values(:,end)=[];
% % %     h=plot_topography(ch_list,log(values(i,:)),false,'10-20',true,false,100);
% % %     colormap("hot");
% % %     caxis([minvaltot,maxvaltot]);
% % %     
% % %     
% % %     if i>1143 && i<1213
% % %         namemeyo='seizure';
% % %     else
% % %         namemeyo='not';
% % %     end
% % %     
% % %     title(namemeyo);
% % %     drawnow
% % %     
% % %     %     pause(0.01);
% % %     
% % % end
% % % 
% % % 
% % % 
% % % 
% % % % Víctor Martínez-Cagigal (2021). Topographic EEG/MEG plot
% % % % (https://www.mathworks.com/matlabcentral/fileexchange/72729-topographic-eeg-meg-plot),
% % % % MATLAB Central File Exchange. Retrieved December 4, 2021.
% % % 
% % % 
% % % % [X,Y] = meshgrid(-1:0.025:1, -1:0.025:1); xq = X(:); yq = Y(:); rTest =
% % % % sqrt(xq.^2 + yq.^2); ckUnit = rTest <= 1; xq = xq(ckUnit); yq =
% % % % yq(ckUnit); nInt = length(xq);
% % % %
% % % % [ntimes,nchan]=size(meanvolts);
% % % %
% % % % EEG_int = zeros(5025, ntimes); % interpolate at each time time point per
% % % % case for jj = 1:ntimes
% % % %     EEG_int(:,jj) = griddata(xy(:,1), xy(:,2), ...
% % % %         meanvolts(1:end-1,jj), xq, yq);
% % % % end
% % % 
% % % 
% % % return
% % % 
% % % 
% % % 
% % % 
% % % 
% % % %% could also potentially do this to find start and end unsupervised
% % % 
% % % % results = multcompare(stats,'Dimension',[1 2,3])
% % % 
% % % % %% bandpower - lead 1 only
% % % %
% % % % lead=1; freqlo=[0.5,4,8,12,30]; freqhi=[4,8,12,30,100];
% % % % freqnames={"\delta","\theta","\alpha","\beta","\gamma"};
% % % % L=length(eeg_raw); ilo=[1,startidx,endidx]; ihi=[startidx,endidx,L];
% % % % inames={"before","during","after"};
% % % %
% % % % Phold=zeros(numel(freqlo),numel(ilo));
% % % %
% % % % for m=1:numel(freqlo)
% % % %
% % % %     for n=1:numel(ilo)
% % % %
% % % %         Phold(m,n)=bandpower(eeg_raw(ilo(n):ihi(n),lead),Fs,[freqlo(m),freqhi(m)]);
% % % %         disp("done " + freqnames{m} + ", " + inames{n});
% % % % %         disp(["done "  freqnames{m} ", "
% % % % inames{n}],'Interpreter','latex');
% % % %
% % % %     end
% % % %
% % % % end
% % % 
% % % % BP_D=bandpower(eeg_raw(1:startidx,lead),Fs,[1,4]);
% % % 
% % % 
% % % % y_D=bandpass(eeg_raw(:,lead),[1,4],Fs);
% % % 
% % % %%
% % % figure; hold on;
% % % plot(eeg_raw(:,lead))
% % % plot(y_D)
% % % return
% % % 
% % % %% FFT
% % % lead=1;
% % % L=length(eeg_raw);
% % % Y=fft(eeg_raw(:,lead));
% % % P2=abs(Y/L);
% % % P1 = P2(1:L/2+1);
% % % P1(2:end-1) = 2*P1(2:end-1);
% % % f = Fs*(0:(L/2))/L;
% % % figure;
% % % plot(f,(P1));
% % % xlim([1,100])
% % % 
% % % return
% % % %% only filter into the bands we want
% % % % eeg_filt
% % % 
% % % %% filter in outside function
% % % eeg_filt=filterEEG(eeg_raw,Fs);
% % % 
% % % %% design custom filter to remove baseline noise
% % % 
% % % fvtool(d)
% % % 
% % % % rather than a notch, assume the data is free of 60Hz noise to avoid
% % % % ringing
% % % 
% % % %% get only specific bands of signal
% % % Y=fft(eeg_raw(1,:));
% % % 
% % % return
% % % %% highpass filter raw data
% % % 
% % % % in order to minimize edge artifacts, we can do some cool manipulation of
% % % % the data to create an artificial edge
% % % padlength=10000;% add these points to edge
% % % eeg_pad=[flipud(eeg_raw(1:padlength,:));...
% % %     eeg_raw;
% % %     flipud(eeg_raw(length(eeg_raw)-padlength:end,:))]';
% % % % eeg_hp=eeg_pad; for i=1:29
% % % %     eeg_hp(i,:)=filter(d,eeg_pad(i,:))';
% % % % end
% % % %
% % % eeg_hp=filter(d,eeg_pad)';
% % % eeg_hp=eeg_hp(padlength+1:end-padlength,:);
% % % 
% % % f=figure;
% % % subplot(2,1,1);
% % % plot(eeg_raw);
% % % subplot(2,1,2);
% % % plot(eeg_hp);
% % % 
% % % return
% % % %% high pass filter to remove moving baseline
% % % disp("hpf of dataset");
% % % eeg_hp = highpass(eeg',0.4,Fs)';
% % % 
% % % %% notch filter to remove 60 Hz noise
% % % disp("notch dataset");
% % % d = designfilt('bandstopiir','FilterOrder',2, ...
% % %     'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
% % %     'DesignMethod','butter','SampleRate',Fs);
% % % % fvtool(d,'fs',fs)
% % % eeg_notch = filtfilt(d,eeg_hp')';
% % % 
% % % %save
% % % % temp={eeg_notch,t,fs,rawmat,colnames}; save("temp.mat","temp");
% % % 
% % % %% rereference
% % % iCz=15;
% % % fCz = -sum(eeg_notch(:,~15))/29;%
% % % 
% % % % eeg_rr=
% % % 
% % % %% sumpower for all
% % % sumpower=sum((eeg_notch.^2),2);
% % % figure;hold on;grid on;
% % % % Sstart=1143;Send=1213;
% % % Sstart=1220;Send=1274;
% % % fill([Sstart,Send,Send,Sstart],[0,0,1,1],'r','facealpha',0.5);
% % % plot(t,sumpower./max(sumpower));
% % % 
% % % 

