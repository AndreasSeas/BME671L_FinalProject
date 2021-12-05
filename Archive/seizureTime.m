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
