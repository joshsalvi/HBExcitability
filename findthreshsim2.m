function [threshrec, noisefloor, totpower, frdetnoise, frstonoise, frdettot, frstotot, const] = findthreshsim2(filename,dwnspl,biftype,stiffind,offset)
% This function finds the noise floor for simulation data in order to
% calculate the threshold required for a peak-finding algorithm. The
% function also calculates the total power in the signal as the upper bound
% for the threshold for both deterministic and stochastic cases.
%
% [noisefloor totpower frdetnoise frstonoise frdettot frstotot] = findthreshsim(filename,dwnspl,biftype,offset)
%
% Thresholds are recommended as the noisefloor plus divisions (0.25%) the
% difference between the total power and noise floor. That is, threshrec is
% equal to max(noisefloor) + n*[max(totpower)/max(noisefloor)]/4,
% where n = 1, 2, 3.
%
% filename: name of MAT file with directory
% dwnspl: downsampling factor
% biftype: 1=supercritical Hopf, 2=SNIC, 3=subcritical Hopf, 4=HB model;
% 5=Hopf with offset
% stiffind: index for HB model; 1 if biftype~=4
% offset: offset of the signal in the time domain
%
% threshrec: recommended thresholds
% noisefloor: noise floor for deterministic(2) and stochastic(1) cases
% totpower: total power, det(2), sto(1)
% frstonoise,frdetnoise: all noise floor values
% frdettot,frstotot: all total power values
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%

% Import data
load(filename);

% Type of bifurcation
        % 1=supercritical Hopf; 2=SNIC; 3=subcritical Hopf, % 4=HB model
if biftype == 1
    Xdet1 = Xdet; Xsto1 = Xsto;
    clear i;
    for j = 1:length(Xdet1)
        for k = 1:length(Xdet1{1})
            Xdet{j}{k} = Xdet1{j}{k}(1,1:dwnspl:end) + offset;  % use only real part
            Xsto{j}{k} = Xsto1{j}{k}(1,1:dwnspl:end) + offset;
        end
    end
    clear Xdet1 Xsto1 Hopfdet Hopfsto
elseif biftype == 2
    for j = 1:length(Xdet)
        for k = 1:length(Xdet{1})
            SNICdet=Xdet;SNICsto=Xsto;
            Xdet{j}{k} = SNICdet{j}{k}(1:dwnspl:end) + offset;
            Xsto{j}{k} = SNICsto{j}{k}(1:dwnspl:end) + offset;
        end
    end
elseif biftype == 3
    Xdet1 = Xdet; Xsto1 = Xsto;
    for j = 1:length(Xdet1)
        for k = 1:length(Xdet1{1})
            Xdet{j}{k} = Xdet1{j}{k}(1,1:dwnspl:end).*sin(Xdet1{j}{k}(2,1:dwnspl:end)) + offset;
            Xsto{j}{k} = Xsto1{j}{k}(1,1:dwnspl:end).*sin(Xsto1{j}{k}(2,1:dwnspl:end)) + offset;
        end
    end
    clear Xdet1 Xsto1
elseif biftype == 4
    Xsto1=Xsto; Xdet1=Xdet;clear Xdet Xsto;
    sizeX = size(Xdet1);
for j = 1:sizeX(1)       % Isolate the appropriate index
    for m = 1:sizeX(3)
        Xsto{m}{j} = Xsto1{j,stiffind,m}(1,1:dwnspl:end) + offset;
        Xdet{m}{j} = Xdet1{j,stiffind,m}(1,1:dwnspl:end) + offset;
    end
end
elseif biftype == 5
    Xsto1=Xsto;Xdet1=Xdet;
    clear Xdet Xsto
    sizeX = size(Xdet1);
    for j = 1:sizeX(1)
        for k = 1:sizeX(2)
            Xsto{k}{j} = Xsto1{j,k}(1,1:dwnspl:end) + offset;
            Xdet{k}{j} = Xdet1{j,k}(1,1:dwnspl:end) + offset;
        end
    end
    
else
    disp('No bifurcation type chosen');
end




% Define time vector
tvec=t(1:dwnspl:end);
tvec=tvec(1:length(Xdet{1}{1}));
% Define the sample rate if not already done
Fs=1/(tvec(2)-tvec(1));
Fnyq = Fs/2;
% High-pass filter
%fsmooth = 0.5; %Hz

% Time range
tmin=(tvec(round(0.1*length(tvec)))); 
tmax=tvec(end);
minindex = find(abs(tvec-tmin)==min(abs(tvec-tmin)));
maxindex = find(abs(tvec-tmax)==min(abs(tvec-tmax)))-1;
Ttotal = tmax-tmin;


% Find amplitudes and frequencies using FFT
T = 1/Fs;
L = length(Xsto{1}{1});
NFFT = (2^2)*2^nextpow2(L);
nw=10;
XsegL = floor(length(Xsto{1}{1})/nw);
welchwin = round(XsegL);
NPSD = floor(NFFT/nw);
noverlap = 0;
winfunc = hamming(welchwin);
freq = 0.005;
Xsine = sin(2*pi*freq.*t);
[Xsinepsd,~] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
for j = 1:length(Xsto)
    for k = 1:length(Xsto{1})
        clear Xstoscaled Xdetscaled frsto frdet
        [Xstofft{j,k}, fstofft{j,k}]= pwelch(Xsto{j}{k},winfunc,noverlap,NPSD,Fs);
        fnyq = findnearest(fstofft{j,k},Fnyq);fnyq=fnyq(1);
        L2=length(fstofft{j,k});
        df=fstofft{j,k}(2)-fstofft{j,k}(1);
        freqend=fstofft{j,k}(end);
        fstoind = find(fstofft{j,k} > 0.001);
        [Xdetfft{j,k}, fdetfft{j,k}] = pwelch(Xdet{j}{k},winfunc,noverlap,NPSD,Fs);
        fdetind = find(fdetfft{j,k} > 0.001);
        fscale = 1e3;
        Xstofft{j,k} = Xstofft{j,k}./fscale;
        Xdetfft{j,k} = Xdetfft{j,k}./fscale;
        if max(2*abs(Xstofft{j,k})) > 1e-10
        Xstofftmaxind = find(Xstofft{j,k}(fstoind)==max(Xstofft{j,k}(fstoind)));
        Xdetfftmaxind = find(Xdetfft{j,k}(fdetind)==max(Xdetfft{j,k}(fdetind)));
        Xstomax=max(Xstofft{j,k}(fstoind));Xstomin=max(Xstofft{j,k}(fstoind));
        Xdetmax=max(Xdetfft{j,k}(fdetind));Xdetmin=max(Xdetfft{j,k}(fdetind));
        xnsd = 1:Xdetfftmaxind;
        xnss = 1:Xstofftmaxind;
        xnsdL = length(xnsd);
        xnssL = length(xnss);
        frsto = Xstofft{j,k};
        frdet = Xdetfft{j,k};
        fsmooth = fstofft{j,k}(Xstofftmaxind);
        % Find mean and standard deviation of noise far away from peak
        %{
        fftstomean = mean(Xstofft{j,k}(round(length(Xstofft{j,k})/2):end));
        fftdetmean = mean(Xdetfft{j,k}(round(length(Xdetfft{j,k})/2):end));
        fftstostd = std(Xstofft{j,k}(round(length(Xstofft{j,k})/2):end));
        fftdetstd = std(Xdetfft{j,k}(round(length(Xdetfft{j,k})/2):end));
        %}
        %{
        fftstomean = mean(Xstofft{j,k}(1:round(xnssL/2)));
        fftdetmean = mean(Xdetfft{j,k}(1:round(xnsdL/2)));
        fftstostd = std(Xstofft{j,k}(1:round(xnssL/2)));
        fftdetstd = std(Xdetfft{j,k}(1:round(xnsdL/2)));
        %}
        
        % Include a window from 0-50 Hz that excludes a 3-Hz band around
        % the peak -> note that 1Hz in points is (1/df)
        %{
        if fstoind(Xstofftmaxind(1)) > 1
        fftstomean = mean([mean(Xstofft{j,k}(1:fstoind(Xstofftmaxind(1))-round(1*freqend/df))) mean(Xstofft{j,k}(fstoind(Xstofftmaxind(1))+round(1*freqend/df):end))]);
        fftstostd = std([mean(Xstofft{j,k}(1:fstoind(Xstofftmaxind(1))-round(1*freqend/df))) mean(Xstofft{j,k}(fstoind(Xstofftmaxind(1))+round(1*freqend/df):end))]);
        else
            fftstomean = mean(Xstofft{j,k}(fstoind(Xstofftmaxind(1))+round(1*freqend/df):end));
            fftstostd = std(Xstofft{j,k}(fstoind(Xstofftmaxind(1))+round(1*freqend/df):end));
        end
        if fdetind(Xdetfftmaxind(1)) > 1
            fftdetmean = mean([mean(Xdetfft{j,k}(1:fdetind(Xdetfftmaxind(1))-round(1*freqend/df))) mean(Xdetfft{j,k}(fdetind(Xdetfftmaxind(1))+round(1*freqend/df):end))]);
            fftdetstd = std([mean(Xdetfft{j,k}(1:fdetind(Xdetfftmaxind(1))-round(1*freqend/df))) mean(Xdetfft{j,k}(fdetind(Xdetfftmaxind(1))+round(1*freqend/df):end))]);
        else
            fftdetmean = mean(Xdetfft{j,k}(fdetind(Xdetfftmaxind(1))+round(1*freqend/df):end));
            fftdetstd = std(Xdetfft{j,k}(fdetind(Xdetfftmaxind(1))+round(1*freqend/df):end));
        end
        %}
        
        
        %}
        
        % Take median and mad up to half of the Nyquist frequency
        fftstomed = median(Xstofft{j,k}(1:fnyq));
        fftstomad = mad(Xstofft{j,k}(1:fnyq));
        fftdetmed = median(Xdetfft{j,k}(1:fnyq));
        fftdetmad = mad(Xdetfft{j,k}(1:fnyq));
        stomax(j,k) = Xstomax;
        detmax(j,k) = Xdetmax;
        % Remove peak
        if biftype == 4 || biftype == 1
        frsto(frsto>(fftstomed+2*fftstomad))=0;
        frdet(frdet>(fftdetmed+2*fftdetmad))=0;
        %frsto=sqrt(frsto);
        %frdet=sqrt(frdet);
        else
        frsto(frsto>(fftstomed+2*fftstomad))=0;
        frdet(frdet>(fftdetmed+2*fftdetmad))=0; 
        end
        
        %}
        %{
        frsto(frsto>(Xstomax*0.5))=0;
        frdet(frdet>(Xdetmax*0.5))=0;
        %}
        frstonoise(j,k) = sum(frsto(1:fnyq));               % find noise floor (Parseval's Theorem)
        frdetnoise(j,k) = sum(frdet(1:fnyq));
        frstotot(j,k) = sum(Xstofft{j,k}(1:fnyq));          % find total variance (Parseval's Theorem)
        frdettot(j,k) = sum(Xdetfft{j,k}(1:fnyq));
        else
        frstonoise(j,k) = 0;
        frdetnoise(j,k) = 0;
        end
        clear xstosmooth xdetsmooth
        xstosmooth = smooth(Xsto{j}{k}(1,:),ceil(Fs/fsmooth));
        xdetsmooth = smooth(Xdet{j}{k}(1,:),ceil(Fs/fsmooth));
        xstovar(j,k) = std(xstosmooth(round(L/4):round(3*L/4))-mean(xstosmooth(round(L/4):round(3*L/4))));
        xdetvar(j,k) = std(xdetsmooth(round(L/4):round(3*L/4))-mean(xdetsmooth(round(L/4):round(3*L/4))));
    end
end

for j = 1:length(Xsto)
    noisefloor(j,1) = sqrt(max(max(frstonoise(j,:))));
    noisefloor(j,2) = sqrt(max(max(frdetnoise(j,:))));
    noisefloorSTD(j,1) = sqrt(max(max(std(frstonoise(j,:)))));
    totpower(j,1) = sqrt(max(max(frstotot(j,:))));
    totpower(j,2) = sqrt(max(max(frdettot(j,:))));


powdiff(j) = totpower(j,1)-noisefloor(j,1);
M = 3;
for m = 1:M
    threshrec(j,m) = m*powdiff(j)/4+max(noisefloor(j,1));
end
end
c1=0;c2=0;
const = [c1 c2];
end
