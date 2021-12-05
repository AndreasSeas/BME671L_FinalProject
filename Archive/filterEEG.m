function [a]=filterEEG(EEG,Fs,padlength)

Fstop = 0.4;
Fpass = 1;
Astop = 65;
Apass = 1;

d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
  'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
  'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','equiripple');

