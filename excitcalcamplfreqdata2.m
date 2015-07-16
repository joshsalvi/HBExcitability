function excitcalcamplfreqdata2(filename1,biftype,stiffind)
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
% Import data
importfile=dir(sprintf('%s%s',filepath,'Extracted Data.mat'));
load(sprintf('%s%s',filepath,importfile.name));

tstart=1;
if exist('Xd_pulse')==1
    Xd=Xd_pulse;
end
% Downsample data
for j = 1:a
    for i = 1:(logdata.data(1,8))
        Xd_dwnspl{i}{j} = Xd(tstart:dwnspl:length(Xd),i,j);  % downsample
        Xd(:,i,j) = Xd(:,i,j) - smooth(Xd(:,i,j),length(Xd(:,i,j))) + offset;   % remove drift
        Xd_dwnspl{i}{j} = Xd_dwnspl{i}{j} - smooth(Xd_dwnspl{i}{j},length(Xd_dwnspl{i}{j})/10) + offset;    % remove drift
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
                   trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pksto{j,k}(l)) Xd_dwnspl{j}{k}(pksto{j,k}(l+1))])-Xd_dwnspl{j}{k}(trsto{j,k}(l))) + Xd_dwnspl{j}{k}(trsto{j,k}(l));
                   trtime{j,k}(l) = sum(Xd_dwnspl{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pksto{j,k}(l)) Xd_dwnspl{j}{k}(pksto{j,k}(l+1))])-Xd_dwnspl{j}{k}(trsto{j,k}(l))) + Xd_dwnspl{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xd_dwnspl{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trsto{j,k}(l)) Xd_dwnspl{j}{k}(trsto{j,k}(l+1))])-Xd_dwnspl{j}{k}(pksto{j,k}(l))) + Xd_dwnspl{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xd_dwnspl{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==loopmaxpk
               pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trsto{j,k}(l)) Xd_dwnspl{j}{k}(trsto{j,k}(l+1))])-Xd_dwnspl{j}{k}(pksto{j,k}(l))) + Xd_dwnspl{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xd_dwnspl{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs; 
               end
           end
       elseif pksto{j,k}(1) < trsto{j,k}(1) && pksto{j,k}(end) < trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pksto{j,k}(l)) Xd_dwnspl{j}{k}(pksto{j,k}(l+1))])-Xd_dwnspl{j}{k}(trsto{j,k}(l))) + Xd_dwnspl{j}{k}(trsto{j,k}(l));
                   trtime{j,k}(l) = sum(Xd_dwnspl{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pksto{j,k}(l)) Xd_dwnspl{j}{k}(pksto{j,k}(l+1))])-Xd_dwnspl{j}{k}(trsto{j,k}(l))) + Xd_dwnspl{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xd_dwnspl{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trsto{j,k}(l)) Xd_dwnspl{j}{k}(trsto{j,k}(l+1))])-Xd_dwnspl{j}{k}(pksto{j,k}(l))) + Xd_dwnspl{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xd_dwnspl{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pksto{j,k}(l)) Xd_dwnspl{j}{k}(pksto{j,k}(l+1))])-Xd_dwnspl{j}{k}(trsto{j,k}(l))) + Xd_dwnspl{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xd_dwnspl{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       elseif pksto{j,k}(1) > trsto{j,k}(1) && pksto{j,k}(end) < trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trsto{j,k}(l)) Xd_dwnspl{j}{k}(trsto{j,k}(l+1))])-Xd_dwnspl{j}{k}(pksto{j,k}(l))) + Xd_dwnspl{j}{k}(pksto{j,k}(l));
                   pktime{j,k}(l) = sum((Xd_dwnspl{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs; 
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pksto{j,k}(l)) Xd_dwnspl{j}{k}(pksto{j,k}(l+1))])-Xd_dwnspl{j}{k}(trsto{j,k}(l))) + Xd_dwnspl{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xd_dwnspl{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trsto{j,k}(l)) Xd_dwnspl{j}{k}(trsto{j,k}(l+1))])-Xd_dwnspl{j}{k}(pksto{j,k}(l))) + Xd_dwnspl{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xd_dwnspl{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs;
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pksto{j,k}(l)) Xd_dwnspl{j}{k}(pksto{j,k}(l+1))])-Xd_dwnspl{j}{k}(trsto{j,k}(l))) + Xd_dwnspl{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xd_dwnspl{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       elseif pksto{j,k}(1) > trsto{j,k}(1) && pksto{j,k}(end) > trsto{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trsto{j,k}(l)) Xd_dwnspl{j}{k}(trsto{j,k}(l+1))])-Xd_dwnspl{j}{k}(pksto{j,k}(l))) + Xd_dwnspl{j}{k}(pksto{j,k}(l));
                   pktime{j,k}(l) = sum((Xd_dwnspl{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pksto{j,k}(l)) Xd_dwnspl{j}{k}(pksto{j,k}(l+1))])-Xd_dwnspl{j}{k}(trsto{j,k}(l))) + Xd_dwnspl{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xd_dwnspl{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trsto{j,k}(l)) Xd_dwnspl{j}{k}(trsto{j,k}(l+1))])-Xd_dwnspl{j}{k}(pksto{j,k}(l))) + Xd_dwnspl{j}{k}(pksto{j,k}(l));
               pktime{j,k}(l) = sum((Xd_dwnspl{j}{k}(trsto{j,k}(l):trsto{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==min([loopmaxpk loopmaxtr])-1
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pksto{j,k}(l)) Xd_dwnspl{j}{k}(pksto{j,k}(l+1))])-Xd_dwnspl{j}{k}(trsto{j,k}(l))) + Xd_dwnspl{j}{k}(trsto{j,k}(l));
               trtime{j,k}(l) = sum(Xd_dwnspl{j}{k}(pksto{j,k}(l):pksto{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       end
      else
          trtime{j,k}=0;pktime{j,k}=0;
      end
    end
end

% Find coefficients of variation for peaks and troughs
for j = 1:sizeP(1)
    for k = 1:sizeP(2)
        meantrtime(j,k) = mean(trtime{j,k});
        meanpktime(j,k) = mean(pktime{j,k});
        CVtrtime(j,k) = std(trtime{j,k})/mean(trtime{j,k});
        CVpktime(j,k) = std(pktime{j,k})/mean(pktime{j,k});
        CDtrtime(j,k) = var(trtime{j,k})/mean(trtime{j,k});
        CDpktime(j,k) = var(pktime{j,k})/mean(pktime{j,k});
    end
end

% Find amplitudes using peaks
for j = 1:sizeP(1)
    for k = 1:sizeP(2)
        if isempty(pksto{j,k}) == 0 && isempty(trsto{j,k}) == 0
            loopmaxtr = length(trsto{j,k});loopmaxpk = length(pksto{j,k});
            if min([loopmaxpk loopmaxtr])>2
            for l = 1:min([loopmaxpk loopmaxtr])-1
                amplsto{j,k}(l) = Xd_dwnspl{j}{k}(pksto{j,k}(l)) - Xd_dwnspl{j}{k}(trsto{j,k}(l));
            end
            else
                amplsto{j,k} = 0;
            end
        else
            amplsto{j,k} = 0;
        end
        PKamplmeansto(j,k) = mean(amplsto{j,k});PKamplsemsto(j,k)=std(amplsto{j,k})/sqrt(length(amplsto{j,k}));
    end
end

% Find amplitudes and frequencies using FFT
T = 1/Fs;
L = length(Xd_dwnspl{1}{1});
NFFT = (2^4)*2^nextpow2(L);
nw=10;
XsegL = floor(length(Xd_dwnspl{1}{1})/nw);
welchwin = round(XsegL);
NPSD = floor(NFFT/nw);
noverlap = 0;
winfunc = hamming(welchwin);
freq = 0.005;
Xsine = sin(2*pi*freq.*t);
[Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
for j = 1:length(Xd_dwnspl)
    for k = 1:length(Xd_dwnspl{1})
        [Xd_dwnsplfft{j,k}, fstofft{j,k}]= pwelch(Xd_dwnspl{j}{k},winfunc,noverlap,NPSD,Fs);
        fstoind = find(fstofft{j,k} > 0.001);
        fscale = 1e3;
        Xd_dwnsplfft{j,k} = Xd_dwnsplfft{j,k}./fscale;
        % Is the peak a minimum height? If not, set ampl/freq to zero.
        if max(2*abs(Xd_dwnsplfft{j,k})) > 1e-3
        Xd_dwnsplfftmaxind = find(Xd_dwnsplfft{j,k}(fstoind)==max(Xd_dwnsplfft{j,k}(fstoind)));
        fftamplsto(j,k) = Xd_dwnsplfft{j,k}(fstoind(Xd_dwnsplfftmaxind(1)));
        fftamplsto(j,k) = (sqrt(fscale.*fftamplsto(j,k).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
        
        fftfreqsto(j,k) = fstofft{j,k}(fstoind(Xd_dwnsplfftmaxind(1)));
        % Find mean and standard deviation of noise far away from peak
        fftstomean = mean(Xd_dwnsplfft{j,k}(fstoind(Xd_dwnsplfftmaxind(1))+1000:fstoind(Xd_dwnsplfftmaxind(1))+2000));
        fftstostd = std(Xd_dwnsplfft{j,k}(fstoind(Xd_dwnsplfftmaxind(1))+1000:fstoind(Xd_dwnsplfftmaxind(1))+2000));
        % Is peak greater than mean plus std? If not, set ampl/freq to
        % zero.
        if fftamplsto(j,k) < (fftstomean+fftstostd)
            fftamplsto(j,k) = 0;
            fftfreqsto(j,k) = 0;
        end

        else
        fftfreqsto(j,k) = 0;
        fftamplsto(j,k) = 0;
        end
    end
end



if exist('threshcheck')==0
    threshcheck=1;
end
if exist('fftpkcheck')==0
    fftpkcheck=1;
end



  if threshcheck==1  
    for n = 1:length(threshrange)
        clear IEIpksto pkd pks
        IEIpksto = cell(length(Xd_dwnspl),length(Xd_dwnspl{1}));
        for j = 1:length(Xd_dwnspl)
            for k = 1:length(Xd_dwnspl{1})
                % Peak finding and CD
                pks{j,k} = PTDetect(Xd_dwnspl{j}{k},threshrange(n));
                for l = 2:length(pks{j,k})
                    IEIpksto{j,k}(l-1) = (pks{j,k}(l) - pks{j,k}(l-1));
                end
                if isempty(IEIpksto{j,k})==0
                    CDstopk{m,n}(j,k) = var(IEIpksto{j,k})/mean(IEIpksto{j,k});
                else
                    CDstopk{m,n}(j,k) = 0;
                end
                % Phase distribution
                clear Xd_dwnsplhilb Xd_dwnsplhilb
                Xd_dwnsplhilb = hilbert(Xd_dwnspl{j}{k}(round(length(Xd_dwnspl{j}{k})/2:end)));
                Xd_dwnsplphase = atan2(imag(Xd_dwnsplhilb),real(Xd_dwnsplhilb));
                clear i;
                VSsto(m,j,k) = abs(1/length(Xd_dwnsplphase)*sum(exp(i*Xd_dwnsplphase))); % vector strength (sto)
            end
        end
        if mod(n,5)==0
            disp(['n = ' num2str(n) '/' num2str(length(threshrange))]);
        end
    end
  end
    if exist('t') == 0
        Fs=50;
    else
        Fs = 1/t(2);
    end
    
    L = length(Xd_dwnspl{1}{1});
    NFFT = (2^2)*2^nextpow2(L);
    nw=10;
    XsegL = floor(length(Xd_dwnspl{1}{1})/nw);
    welchwin = round(XsegL);
    NPSD = floor(NFFT/nw);
    noverlap = 0;
    winfunc = hamming(welchwin);
    freq = 0.005;
    Xsine = sin(2*pi*freq.*t);
    [Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
    winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
    for j = 1:length(Xd_dwnspl)
        for k = 1:length(Xd_dwnspl{1})
            if fftpkcheck==1
            clear fftpkspkheight fftpkdpkheight
            [Xd_dwnsplfft, fstofft, Xd_dwnsplfftCI]= pwelch(Xd_dwnspl{j}{k},winfunc,noverlap,NPSD,Fs,'ConfidenceLevel',0.999);
            [a1s b1s] = twoclass(Xd_dwnsplfft,1e-4);
            pkfs = PTDetect(Xd_dwnsplfft,max([a1s b1s]));
            mxs = 1; 
            for q = 1:length(pkfs)  % test each peak (stochastic) -> is it statistically distinct from peaks with less then half of its power on both sides?
                fhps = find(Xd_dwnsplfft<=(Xd_dwnsplfft(pkfs(q))/2));
                testL = sign(Xd_dwnsplfftCI(pkfs(q),1)-Xd_dwnsplfftCI(fhps(fhps<pkfs(q)),2));
                testH = sign(Xd_dwnsplfftCI(pkfs(q),1)-Xd_dwnsplfftCI(fhps(fhps>pkfs(q)),2));
                if sum(testL==1) >= 1 && sum(testH==1) >= 1 
                    fftpkspkheight(q) = Xd_dwnsplfft(pkfs(q));
                    fftpkspkfreq(q) = fstofft(pkfs(q));
                    fhpsL = fhps(fhps<pkfs(q));fhpsH = fhps(fhps>pkfs(q));
                    Lind = fhpsL(testL==1);Lind=Lind(end);
                    Hind = fhpsH(testH==1);Hind=Hind(1);
                    fftpkswidth(m,q,1:2) = [fstofft(Lind) fstofft(Hind)];
                    mxs = mxs+1;
                end
            end
            if mxs > 1
                clear ind1
                ind1 = find(fftpkspkheight(:)==max(fftpkspkheight(:))); ind1=ind1(1);
                if numel(ind1)==0;ind1=1;
                end
                fftmaxpkstoind(m,j,k) = pkfs(ind1);
                fftmaxpkheightsto(m,j,k) = Xd_dwnsplfft(fftmaxpkstoind(m,j,k));
                fftmaxpkfreqtsto(m,j,k) =  fstofft(fftmaxpkstoind(m,j,k));
                fftmaxpkwidthtstoL(m,j,k) = fftpkswidth(m,ind1,1);
                fftmaxpkwidthtstoR(m,j,k) = fftpkswidth(m,ind1,2);
                Qsto(m,j,k) = fftmaxpkfreqtsto(m,j,k)/(fftmaxpkwidthtstoR(m,j,k)-fftmaxpkwidthtstoL(m,j,k));
            else
                fftmaxpkheightsto(m,j,k) = 0;
                fftmaxpkfreqtsto(m,j,k) = 0;
                fftmaxpkwidthtstoL(m,j,k) = 0;
                fftmaxpkwidthtstoR(m,j,k) = 0;
                Qsto(m,j,k) = 0;
            end
            
          
            end
            if k == 1
                disp('Dip test...');
                disp([num2str(j) '/' num2str(length(Xd_dwnspl))])
            end
            % Find the dip statistic and associated p-value
            Nboot = 10^1;   % number of bootstraps
            %[as bs] = hist(Xd_dwnspl{j}{k},freedmandiaconis(Xd_dwnspl{j}{k}));as=as./sum(as);
            %[ad bd] = hist(Xd_dwnspl{j}{k},freedmandiaconis(Xd_dwnspl{j}{k}));ad=ad./sum(ad);
            [dipsto(m,j,k), punisto(m,j,k), Xlowsto(m,j,k), Xupsto(m,j,k)]=HartigansDipSignifTest(Xd_dwnspl{j}{k},Nboot);
            if threshcheck ~= 1
                if mod(k,25)==0
                    disp([num2str(k) '/' num2str(length(Xd_dwnspl{1}))])
                end
            end
           
        end
    end
    if threshcheck==1 && fftpkcheck==1
disp(['m = ' num2str(m) '/5']);
disp(['Saving' num2str(m) '...']);

    end

%Pre-allocate memory
mu1s = zeros(length(threshrange),length(Xd_dwnspl));
CD1s = zeros(length(threshrange),length(Xd_dwnspl));
sizeC1 = size(CDstopk);
sizeC2 = size(CDstopk{1,1});

if threshcheck ==1 && fftpkcheck ==1
for n = 1:sizeC1(2)
    for j = 1:sizeC2(1)
        for k = 1:sizeC2(2)
            CDstomean{n}(j,k) = mean([CDstopk{1,n}(j,k) CDstopk{2,n}(j,k) CDstopk{3,n}(j,k) CDstopk{4,n}(j,k) CDstopk{5,n}(j,k)]);
            CDstosem{n}(j,k) = std([CDstopk{1,n}(j,k) CDstopk{2,n}(j,k) CDstopk{3,n}(j,k) CDstopk{4,n}(j,k) CDstopk{5,n}(j,k)])/sqrt(5);
            Qstomean(j,k) = mean([Qsto(1,j,k) Qsto(2,j,k) Qsto(3,j,k) Qsto(4,j,k) Qsto(5,j,k)]);
            Qstosem(j,k) = std([Qsto(1,j,k) Qsto(2,j,k) Qsto(3,j,k) Qsto(4,j,k) Qsto(5,j,k)])/sqrt(5);
            VSstomean(j,k) = mean([VSsto(1,j,k) VSsto(2,j,k) VSsto(3,j,k) VSsto(4,j,k) VSsto(5,j,k)]);
            VSstosem(j,k) = std([VSsto(1,j,k) VSsto(2,j,k) VSsto(3,j,k) VSsto(4,j,k) VSsto(5,j,k)])/sqrt(5);
        end
        CD1s1 = findnearest(CDstomean{n}(j,:),1);
        if isempty(CD1s1)==0;
            CD1s(n,j)=CDstomean{n}(j,CD1s1(1));
        else
            CD1s(n,j)=NaN;
        end
        mu1s(n,j) = mu(CD1s1(1));       % control parameter at which CD ~ 1 (stochastic)
    end
end
end

