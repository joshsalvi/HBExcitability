function excitsimcalcamplfreq(filename,biftype)
% Calculates the time for individual pulses, coefficients of variation for
% them, amplitudes from the peaks, and amplitudes and frequencies from the
% FFT (normalized by the length of the signal).
%
% excitsimcalcamplfreq(filename,biftype)
%
% biftype = 1,2,3 for SNIC/Hopf/sHopf; biftype = 4 for HB model
%
% jsalvi@rockefeller.edu

load(filename)

if biftype == 1 || biftype == 2 || biftype == 3
    
% Calculate the time for each peak and trough
sizeP = size(pksto);
nperc = 0.7; % percentage for threshold
for j = 1:sizeP(1)
    for k = 1:sizeP(2)
      if isempty(pksto{j,k}) == 0 && isempty(trsto{j,k}) == 0
       loopmaxtr = length(trsto{j,k});loopmaxpk = length(pksto{j,k});
       if pksto{j,k}(1) < trsto{j,k}(1) && pksto{j,k}(end) > trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
                   trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))>trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))>trthresh)/Fs;
               pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>pkthresh))/Fs; 
               end
               if l==loopmaxpk
               pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>pkthresh))/Fs; 
               end
           end
       elseif pksto{j,k}(1) < trsto{j,k}(1) && pksto{j,k}(end) < trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
                   trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))>trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))>trthresh)/Fs;
               pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>pkthresh))/Fs; 
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))>trthresh)/Fs;
               end
           end
       elseif pksto{j,k}(1) > trsto{j,k}(1) && pksto{j,k}(end) < trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
                   pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>pkthresh))/Fs; 
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))>trthresh)/Fs;
               pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>pkthresh))/Fs;
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))>trthresh)/Fs;
               end
           end
       elseif pksto{j,k}(1) > trsto{j,k}(1) && pksto{j,k}(end) > trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
                   pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>pkthresh))/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))>trthresh)/Fs;
               pkthresh = nperc*(mean([Xsto{j}{k}(trsto{j,k}(l)) Xsto{j}{k}(trsto{j,k}(l+1))])-Xsto{j}{k}(pksto{j,k}(l))) + Xsto{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xsto{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>pkthresh))/Fs; 
               end
               if l==min([loopmaxpk loopmaxtr])-1
               trthresh = nperc*(mean([Xsto{j}{k}(pksto{j,k}(l)) Xsto{j}{k}(pksto{j,k}(l+1))])-Xsto{j}{k}(trsto{j,k}(l))) + Xsto{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xsto{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))>trthresh)/Fs;
               end
           end
       end
       end
    end
end

% Find coefficients of variation for peaks and troughs
for j = 1:sizeP(1)
    for k = 1:sizeP(2)
        CVtrtime(j,k) = std(trtime{j,k})/mean(trtime{j,k});
        CVpktime(j,k) = std(pktime{j,k})/mean(pktime{j,k});
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
f = Fs/2*linspace(0,1,NFFT/2+1);
for j = 1:length(Xsto)
    for k = 1:length(Xsto{1})
        clear Xstofft Xdetfft
        Xstofft = fft(Xsto{j}{k},NFFT)/L;
        Xdetfft = fft(Xdet{j}{k},NFFT)/L;
        if max(2*abs(Xstofft)) > 1e-2
        Xstofftmaxind = find(2*abs(Xstofft(1:NFFT/2+1))==max(2*abs(Xstofft(1:NFFT/2+1))));
        Xdetfftmaxind = find(2*abs(Xdetfft(1:NFFT/2+1))==max(2*abs(Xdetfft(1:NFFT/2+1))));
        fftfreqsto(j,k) = f(Xstofftmaxind(1));
        fftfreqdet(j,k) = f(Xdetfftmaxind(1));
        fftamplsto(j,k) = 2*abs(Xstofft(Xstofftmaxind(1)));
        fftampldet(j,k) = 2*abs(Xdetfft(Xdetfftmaxind(1)));
        else
        fftfreqsto(j,k) = 0;
        fftfreqdet(j,k) = 0;
        fftamplsto(j,k) = 0;
        fftampldet(j,k) = 0;  
        end
    end
end
    disp('Saving...');
    fnind=find(filename=='/');fdir = filename(1:fnind(end));fprefix=filename(fnind(end)+1:end-4);
    save(sprintf('%s%s%s%s',fdir,'AmplFreq-',fprefix,'-analyzed.mat'),'fftfreqsto','fftfreqdet','fftamplsto','fftampldet','PKamplmeansto','PKamplmeandet','PKamplsemsto','PKamplsemdet','CVtrtime','CVpktime','trtime','pktime')
    disp('Finished.');
    
elseif biftype==4
    disp('Saving...');
    fnind=find(filename=='/');fdir = filename(1:fnind(end));fprefix=filename(fnind(end)+1:end-4);
    save(sprintf('%s%s%s%s',fdir,'AmplFreq-',fprefix,'-analyzed.mat'),'fftfreqsto','fftfreqdet','fftamplsto','fftampldet','PKamplmeansto','PKamplmeandet','PKamplsemsto','PKamplsemdet','CVtrtime','CVpktime','trtime','pktime')
    disp('Finished.');
    
end

end
        
