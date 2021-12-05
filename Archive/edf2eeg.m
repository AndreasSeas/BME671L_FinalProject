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