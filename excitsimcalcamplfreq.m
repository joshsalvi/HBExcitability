function excitsimcalcamplfreq(filename1,biftype,stiffind)
% Calculates the time for individual pulses, coefficients of variation for
% them, amplitudes from the peaks, and amplitudes and frequencies from the
% FFT (normalized by the length of the signal).
%
% excitsimcalcamplfreq(filename,biftype,stiffind)
%
% biftype = 1,2,3 for SNIC/Hopf/sHopf; biftype = 4 for HB model
% stiffind : stiffness/force index for biftype==4
%
% jsalvi@rockefeller.edu
load(filename1)

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
% Calculate the time for each peak and trough
sizeP = size(pksto);
nperc = 0.5; % percentage for threshold
for j = 1:sizeP(1)
    for k = 1:sizeP(2)
      if isempty(pksto{j,k}) == 0 && isempty(trsto{j,k}) == 0
       loopmaxtr = length(trsto{j,k});loopmaxpk = length(pksto{j,k});
       if pksto{j,k}(1) < trsto{j,k}(1) && pksto{j,k}(end) > trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
                   trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==loopmaxpk
               pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs; 
               end
           end
       elseif pksto{j,k}(1) < trsto{j,k}(1) && pksto{j,k}(end) < trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
                   trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       elseif pksto{j,k}(1) > trsto{j,k}(1) && pksto{j,k}(end) < trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
                   pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs; 
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs;
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       elseif pksto{j,k}(1) > trsto{j,k}(1) && pksto{j,k}(end) > trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
                   pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==min([loopmaxpk loopmaxtr])-1
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       end
      else
          trtime{j,k}=0;pktime{j,k}=0;
      end
    end
end

sizeP = size(pkdet);
nperc = 0.5; % percentage for threshold
for j = 1:sizeP(1)
    for k = 1:sizeP(2)
      if isempty(pkdet{j,k}) == 0 && isempty(trdet{j,k}) == 0
       loopmaxtr = length(trdet{j,k});loopmaxpk = length(pkdet{j,k});
       if pkdet{j,k}(1) < trdet{j,k}(1) && pkdet{j,k}(end) > trdet{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   trthresh = nperc*(mean([Xdet{j}{k}(pkdet{j,k}(l)) Xdet{j}{k}(pkdet{j,k}(l+1))])-Xdet{j}{k}(trdet{j,k}(l))) + Xdet{j}{k}(trdet{j,k}(l));
                   trtimedet{j,k}(l) = sum(Xdet{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xdet{j}{k}(pkdet{j,k}(l)) Xdet{j}{k}(pkdet{j,k}(l+1))])-Xdet{j}{k}(trdet{j,k}(l))) + Xdet{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xdet{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xdet{j}{k}(trdet{j,k}(l)) Xdet{j}{k}(trdet{j,k}(l+1))])-Xdet{j}{k}(pkdet{j,k}(l))) + Xdet{j}{k}(pkdet{j,k}(l));
               pktimedet{j,k}(l) = sum((Xdet{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==loopmaxpk
               pkthresh = nperc*(mean([Xdet{j}{k}(trdet{j,k}(l)) Xdet{j}{k}(trdet{j,k}(l+1))])-Xdet{j}{k}(pkdet{j,k}(l))) + Xdet{j}{k}(pkdet{j,k}(l));
               pktimedet{j,k}(l) = sum((Xdet{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs; 
               end
           end
       elseif pkdet{j,k}(1) < trdet{j,k}(1) && pkdet{j,k}(end) < trdet{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   trthresh = nperc*(mean([Xdet{j}{k}(pkdet{j,k}(l)) Xdet{j}{k}(pkdet{j,k}(l+1))])-Xdet{j}{k}(trdet{j,k}(l))) + Xdet{j}{k}(trdet{j,k}(l));
                   trtimedet{j,k}(l) = sum(Xdet{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xdet{j}{k}(pkdet{j,k}(l)) Xdet{j}{k}(pkdet{j,k}(l+1))])-Xdet{j}{k}(trdet{j,k}(l))) + Xdet{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xdet{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xdet{j}{k}(trdet{j,k}(l)) Xdet{j}{k}(trdet{j,k}(l+1))])-Xdet{j}{k}(pkdet{j,k}(l))) + Xdet{j}{k}(pkdet{j,k}(l));
               pktimedet{j,k}(l) = sum((Xdet{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xdet{j}{k}(pkdet{j,k}(l)) Xdet{j}{k}(pkdet{j,k}(l+1))])-Xdet{j}{k}(trdet{j,k}(l))) + Xdet{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xdet{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       elseif pkdet{j,k}(1) > trdet{j,k}(1) && pkdet{j,k}(end) < trdet{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xdet{j}{k}(trdet{j,k}(l)) Xdet{j}{k}(trdet{j,k}(l+1))])-Xdet{j}{k}(pkdet{j,k}(l))) + Xdet{j}{k}(pkdet{j,k}(l));
                   pktimedet{j,k}(l) = sum((Xdet{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs; 
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xdet{j}{k}(pkdet{j,k}(l)) Xdet{j}{k}(pkdet{j,k}(l+1))])-Xdet{j}{k}(trdet{j,k}(l))) + Xdet{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xdet{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xdet{j}{k}(trdet{j,k}(l)) Xdet{j}{k}(trdet{j,k}(l+1))])-Xdet{j}{k}(pkdet{j,k}(l))) + Xdet{j}{k}(pkdet{j,k}(l));
               pktimedet{j,k}(l) = sum((Xdet{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs;
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xdet{j}{k}(pkdet{j,k}(l)) Xdet{j}{k}(pkdet{j,k}(l+1))])-Xdet{j}{k}(trdet{j,k}(l))) + Xdet{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xdet{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       elseif pkdet{j,k}(1) >= trdet{j,k}(1) && pkdet{j,k}(end) >= trdet{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xdet{j}{k}(trdet{j,k}(l)) Xdet{j}{k}(trdet{j,k}(l+1))])-Xdet{j}{k}(pkdet{j,k}(l))) + Xdet{j}{k}(pkdet{j,k}(l));
                   pktimedet{j,k}(l) = sum((Xdet{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xdet{j}{k}(pkdet{j,k}(l)) Xdet{j}{k}(pkdet{j,k}(l+1))])-Xdet{j}{k}(trdet{j,k}(l))) + Xdet{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xdet{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xdet{j}{k}(trdet{j,k}(l)) Xdet{j}{k}(trdet{j,k}(l+1))])-Xdet{j}{k}(pkdet{j,k}(l))) + Xdet{j}{k}(pkdet{j,k}(l));
               pktimedet{j,k}(l) = sum((Xdet{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==min([loopmaxpk loopmaxtr])-1
               trthresh = nperc*(mean([Xdet{j}{k}(pkdet{j,k}(l)) Xdet{j}{k}(pkdet{j,k}(l+1))])-Xdet{j}{k}(trdet{j,k}(l))) + Xdet{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xdet{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       end
      else
          trtimedet{j,k}=0;pktimedet{j,k}=0;
      end
    end
end

% Find coefficients of variation for peaks and troughs
for j = 1:sizeP(1)
    for k = 1:sizeP(2)
        meantrtime(j,k) = mean(trtime{j,k});
        meanpktime(j,k) = mean(pktime{j,k});
        meantrtimedet(j,k) = mean(trtimedet{j,k});
        meanpktimedet(j,k) = mean(pktimedet{j,k});
        CVtrtime(j,k) = std(trtime{j,k})/mean(trtime{j,k});
        CVpktime(j,k) = std(pktime{j,k})/mean(pktime{j,k});
        CVtrtimedet(j,k) = std(trtimedet{j,k})/mean(trtimedet{j,k});
        CVpktimedet(j,k) = std(pktimedet{j,k})/mean(pktimedet{j,k});
        CDtrtime(j,k) = var(trtime{j,k})/mean(trtime{j,k});
        CDpktime(j,k) = var(pktime{j,k})/mean(pktime{j,k});
        CDtrtimedet(j,k) = var(trtimedet{j,k})/mean(trtimedet{j,k});
        CDpktimedet(j,k) = var(pktimedet{j,k})/mean(pktimedet{j,k});
    end
end

% Find amplitudes using peaks
for j = 1:sizeP(1)
    for k = 1:sizeP(2)
        if isempty(pksto{j,k}) == 0 && isempty(trsto{j,k}) == 0
            loopmaxtr = length(trsto{j,k});loopmaxpk = length(pksto{j,k});
            if min([loopmaxpk loopmaxtr])>2
            for l = 1:min([loopmaxpk loopmaxtr])-1
                amplsto{j,k}(l) = Xsto{j}{k}(pksto{j,k}(l)) - Xsto{j}{k}(trsto{j,k}(l));
            end
            else
                amplsto{j,k} = 0;
            end
        else
            amplsto{j,k} = 0;
        end
        if isempty(pkdet{j,k}) == 0 && isempty(trdet{j,k}) == 0
            loopmaxtr = length(trdet{j,k});loopmaxpk = length(pkdet{j,k});
            if min([loopmaxpk loopmaxtr])>2
            for l = 1:min([loopmaxpk loopmaxtr])-1
                ampldet{j,k}(l) = Xdet{j}{k}(pkdet{j,k}(l)) - Xdet{j}{k}(trdet{j,k}(l));
            end
            else
                ampldet{j,k}=0;
            end
        else
            ampldet{j,k}=0;
        end
        PKamplmeansto(j,k) = mean(amplsto{j,k});PKamplsemsto(j,k)=std(amplsto{j,k})/sqrt(length(amplsto{j,k}));
        PKamplmeandet(j,k) = mean(ampldet{j,k});PKamplsemdet(j,k)=std(ampldet{j,k})/sqrt(length(ampldet{j,k}));
    end
end

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
        [Xstofft{j,k}, fstofft{j,k}]= pwelch(Xsto{j}{k},winfunc,noverlap,NPSD,Fs);
        fstoind = find(fstofft{j,k} > 0.001);
        [Xdetfft{j,k}, fdetfft{j,k}] = pwelch(Xdet{j}{k},winfunc,noverlap,NPSD,Fs);
        fdetind = find(fdetfft{j,k} > 0.001);
        fscale = 1e3;
        Xstofft{j,k} = Xstofft{j,k}./fscale;
        Xdetfft{j,k} = Xdetfft{j,k}./fscale;
        % Is the peak a minimum height? If not, set ampl/freq to zero.
        if max(2*abs(Xstofft{j,k})) > 1e-3
        Xstofftmaxind = find(Xstofft{j,k}(fstoind)==max(Xstofft{j,k}(fstoind)));
        Xdetfftmaxind = find(Xdetfft{j,k}(fdetind)==max(Xdetfft{j,k}(fdetind)));
        fftamplsto(j,k) = Xstofft{j,k}(fstoind(Xstofftmaxind(1)));
        fftampldet(j,k) = Xdetfft{j,k}(fdetind(Xdetfftmaxind(1)));
        fftamplsto(j,k) = (sqrt(fscale.*fftamplsto(j,k).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
        fftampldet(j,k) = (sqrt(fscale.*fftampldet(j,k).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
        
        fftfreqsto(j,k) = fstofft{j,k}(fstoind(Xstofftmaxind(1)));
        fftfreqdet(j,k) = fdetfft{j,k}(fdetind(Xdetfftmaxind(1)));
        % Find mean and standard deviation of noise far away from peak
        fftstomean = mean(Xstofft{j,k}(fstoind(Xstofftmaxind(1))+1000:fstoind(Xstofftmaxind(1))+2000));
        fftstostd = std(Xstofft{j,k}(fstoind(Xstofftmaxind(1))+1000:fstoind(Xstofftmaxind(1))+2000));
        fftdetmean = mean(Xdetfft{j,k}(fstoind(Xdetfftmaxind(1))+1000:fstoind(Xdetfftmaxind(1))+2000));
        fftdetstd = std(Xdetfft{j,k}(fstoind(Xdetfftmaxind(1))+1000:fstoind(Xdetfftmaxind(1))+2000));
        % Is peak greater than mean plus std? If not, set ampl/freq to
        % zero.
        if fftamplsto(j,k) < (fftstomean+fftstostd)
            fftamplsto(j,k) = 0;
            fftfreqsto(j,k) = 0;
        end
        if fftampldet(j,k) < (fftdetmean+fftdetstd)
            fftampldet(j,k) = 0;
            fftfreqdet(j,k) = 0;
        end
        
        else
        fftfreqsto(j,k) = 0;
        fftfreqdet(j,k) = 0;
        fftamplsto(j,k) = 0;
        fftampldet(j,k) = 0;  
        end
    end
end

    disp('Saving...');
    fnind2=find(filename1=='/');fdir2 = filename1(1:fnind2(end));fprefix2=filename1(fnind2(end)+1:end-4);
    save(sprintf('%s%s%s%s%s%s',fdir2,'AmplFreq-',fprefix2,'-StiffForceInd-',num2str(stiffind),'-analyzed.mat'),'fftfreqsto','fftfreqdet','fftamplsto','fftampldet','PKamplmeansto','PKamplmeandet','PKamplsemsto','PKamplsemdet','CVtrtime','CVpktime','trtime','pktime','meantrtime','meanpktime','CVtrtimedet','CVpktimedet','trtimedet','pktimedet','meantrtimedet','meanpktimedet','CDpktime','CDtrtime','CDpktimedet','CDtrtimedet','Xstofft','Xdetfft','fstofft','fdetfft','-v7.3')

    disp('Finished.');

end
        
