function [threshrec noisefloor totpower frdetnoise frstonoise frdettot frstotot] = findthreshsim(filename,dwnspl,biftype,stiffind,offset)
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
% biftype: 1=supercritical Hopf, 2=SNIC, 3=subcritical Hopf, 4=HB model
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
    sizeX = size(Xdet);
    Xdet1 = Xdet; Xsto1 = Xsto;
    clear Xdet Xsto
    for j = 1:sizeX(1)
        for k = 1:sizeX(2)
            for l = 1:sizeX(3)
                Xdet{j}{k}{l} = Xdet1{j,k,l}(1,1:dwnspl:end) + offset;
                Xsto{j}{k}{l} = Xsto1{j,k,l}(1,1:dwnspl:end) + offset;
            end
        end
    end
    clear Xdet1 Xsto1
else
    disp('No bifurcation type chosen');
end

if biftype == 1 || biftype == 2 || biftype == 3
    stiffind=1;
elseif biftype == 4
        Xsto1=Xsto; Xdet1=Xdet; pksto1=pksto;pkdet1=pkdet;trsto1=trsto;trdet1=trdet;clear Xdet Xsto pkdet pksto trdet trsto;
for j = 1:length(Xdet1)       % Isolate the appropriate index
    for m = 1:length(Xdet1{1}{1})
        Xsto{m}{j} = Xsto1{j}{stiffind}{m};
        Xdet{m}{j} = Xdet1{j}{stiffind}{m};
        pksto{m,j} = pksto1{j,stiffind,m};pkdet{m,j} = pkdet1{j,stiffind,m};
        trsto{m,j} = trsto1{j,stiffind,m};trdet{m,j} = pkdet1{j,stiffind,m};
    end
end
end 



% Define time vector
tvec=t(1:dwnspl:end);
tvec=tvec(1:length(Xdet{1}{1}));
% Define the sample rate if not already done
Fs=1/(tvec(2)-tvec(1));

% Time range
tmin=(tvec(round(0.1*length(tvec)))); 
tmax=tvec(end);
minindex = find(abs(tvec-tmin)==min(abs(tvec-tmin)));
maxindex = find(abs(tvec-tmax)==min(abs(tvec-tmax)))-1;
Ttotal = tmax-tmin;


% Find amplitudes and frequencies using FFT
T = 1/Fs;
L = length(Xsto{1}{1});
NFFT = (2^4)*2^nextpow2(L);
nw=10;
XsegL = floor(length(Xsto{1}{1})/nw);
welchwin = round(XsegL);
NPSD = floor(NFFT/nw);
noverlap = 0;
winfunc = hamming(welchwin);
freq = 0.005;
Xsine = sin(2*pi*freq.*t);
[Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
for j = 1:length(Xsto)
    for k = 1:length(Xsto{1})
        clear Xstoscaled Xdetscaled frsto frdet
        [Xstofft{j,k}, fstofft{j,k}]= pwelch(Xsto{j}{k},winfunc,noverlap,NPSD,Fs);
        fstoind = find(fstofft{j,k} > 0.0001);
        [Xdetfft{j,k}, fdetfft{j,k}] = pwelch(Xdet{j}{k},winfunc,noverlap,NPSD,Fs);
        fdetind = find(fdetfft{j,k} > 0.0001);
        fscale = 1e3;
        Xstofft{j,k} = Xstofft{j,k}./fscale;
        Xdetfft{j,k} = Xdetfft{j,k}./fscale;
        if max(2*abs(Xstofft{j,k})) > 1e-2
        Xstofftmaxind = find(Xstofft{j,k}(fstoind)==max(Xstofft{j,k}(fstoind)));
        Xdetfftmaxind = find(Xdetfft{j,k}(fdetind)==max(Xdetfft{j,k}(fdetind)));
        frsto = Xstofft{j,k};
        frdet = Xdetfft{j,k};
        % Find mean and standard deviation of noise far away from peak
        fftstomean = mean(Xstofft{j,k}(fstoind(Xstofftmaxind(1))+1000:fstoind(Xstofftmaxind(1))+2000));
        fftstostd = std(Xstofft{j,k}(fstoind(Xstofftmaxind(1))+1000:fstoind(Xstofftmaxind(1))+2000));
        fftdetmean = mean(Xdetfft{j,k}(fstoind(Xdetfftmaxind(1))+1000:fstoind(Xdetfftmaxind(1))+2000));
        fftdetstd = std(Xdetfft{j,k}(fstoind(Xdetfftmaxind(1))+1000:fstoind(Xdetfftmaxind(1))+2000));
        % Remove peak
        frsto(frsto>(fftstomean+fftstostd))=0;
        frdet(frdet>(fftdetmean+fftdetstd))=0;
        frstonoise(j,k) = sum(frsto);               % find noise floor (Parseval's Theorem)
        frdetnoise(j,k) = sum(frdet);
        frstotot(j,k) = sum(Xstofft{j,k});          % find total variance (Parseval's Theorem)
        frdettot(j,k) = sum(Xdetfft{j,k});
        else
        frstonoise(j,k) = 0;
        frdetnoise(j,k) = 0;
        end
    end
end

noisefloor(1) = max(max(frstonoise));
noisefloor(2) = max(max(frdetnoise));
totpower(1) = max(max(frstotot));
totpower(2) = max(max(frdettot));

powdiff = max(totpower)-max(noisefloor);
powdiff = powdiff/4;
for j = 1:3
    threshrec(j) = max(noisefloor)+j*powdiff;
end


end
