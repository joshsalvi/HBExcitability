function [pkspikeratedet CDdetpk fftampldet fftfreqdet PKamplmeandet c1det c2det] = excitabilityanalysisdata(filepath,dwnspl,offset,nnthresh,c12thresh,scaleyn,saveyn)
% This function imports raw data and performs a peak-finding algorithm
% analysis. 
%
% [pkspikeratedet CDdetpk fftampldet fftfreqdet PKamplmeandet]= excitabilityanalysisdata(filepath,dwnspl,offset,nnthresh,c12thresh,saveyn)
%
% Ex. [pks, CD, PSDfreqampl] = excitabilityanalysisdata('/Users/joshsalvi/Documents/Lab/Lab/Clamp
% Data/2014-09-30.01/Ear 1/Cell 2/', 10,0,0);
%
% filepath : directory for your file
% dwnspl : downsampling factor (1= no downsampling)
% offset : offset applied to time trace
% nnthresh : nearest-neighbor clustering on a defined point or absolute
% threshold? (1= nearest neighbor; 2= absolute threshold)
% c12thresh : threshold based upon nnthresh (index number or threshold
% value) - if index number use the form [a b]
% scaleyn : use power spectra to rescale the data? (1=yes)
% saveyn : save? (1=yes)
%
% If you want to perform a proper histogram analysis, modify the code to
% change the thresholds. The raw data file must be a matlab MAT file called
% Extracted Data.mat. This is generated from importVLCdata().
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

% Define time vector
tvec=time(tstart:dwnspl:length(time));
clear time;
% Define the sample rate if not already done
Fs=1/(tvec(2)-tvec(1))*1e3;



% Time range
tmin=tvec(round(0.1*length(tvec)));
tmax=tvec(end);
minindex = find(abs(tvec-tmin)==min(abs(tvec-tmin)));
maxindex = length(Xd_dwnspl{1}{1});
%maxindex = find(abs(tvec-tmax)==min(abs(tvec-tmax)))-1;
Ttotal = tmax-tmin;

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
Xsine = sin(2*pi*freq.*tvec);
[Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
for j = 1:length(Xd_dwnspl)
    for k = 1:length(Xd_dwnspl{1})
        clear Xstoscaled Xd_dwnsplscaled
        [Xd_dwnsplfft{j,k}, fdetfft{j,k}] = pwelch(Xd_dwnspl{j}{k},winfunc,noverlap,NPSD,Fs);
        g = findnearest(fdetfft{1,1},100);g=g(1);
        g2 = findnearest(fdetfft{1,1},1000);g2=g2(1);
        fdetind = find(fdetfft{j,k} > 0.0001);
        fscale = 10^3;
        Xd_dwnsplscaled = Xd_dwnsplfft{j,k}./fscale;
        if max(2*abs(Xd_dwnsplfft{j,k})) > 1e-2
        Xd_dwnsplfftmaxind = find(Xd_dwnsplfft{j,k}(fdetind)==max(Xd_dwnsplfft{j,k}(fdetind)));
        fftampldet(j,k) = Xd_dwnsplfft{j,k}(fdetind(Xd_dwnsplfftmaxind(1)));
        fftampldet(j,k) = (sqrt(fscale.*fftampldet(j,k).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
        fftfreqdet(j,k) = fdetfft{j,k}(fdetind(Xd_dwnsplfftmaxind(1)));
        else
        fftfreqdet(j,k) = 0;
        fftampldet(j,k) = 0;  
        end
        if scaleyn == 1
        Xscale2(j,k) = rms(Xd_dwnsplfft{1,1}(g:g2))/rms(Xd_dwnsplfft{j,k}(g:g2)); %rescale the signals from the noise
        else
            Xscale2(j,k) = 1;
        end
    end
end

% Prepare to remove the mean from time traces
fdrift = 0.001; %kHz
fsmooth = 3*fdrift;deltat=1/Fs;ws = (1/fsmooth)/deltat;
if rem(floor(ws),2) == 1 %odd
    ws = floor(ws);
else %even
    ws = floor(ws)+1; 
end
% Find the running mean and std for a time greater than the longest period (1/fdrift)
tss = 1.5*(1/fsmooth);
ssws = round(tss/deltat);
nstdss = 2;
stdfrac = 0.1;%Get to 90% of std ss
%Discount end discontinuities
ndev = 2;

% Pre-allocate memory
detmean = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
dethista = cell(length(Xd_dwnspl),length(Xd_dwnspl{1})); dethistb = cell(length(Xd_dwnspl),length(Xd_dwnspl{1}));
detvar = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1})); 


% Generate histograms and calculate statistics
for k = 1:size(Xd_dwnspl{1})
    for j = 1:size(Xd_dwnspl)
        Xd_dwnspl{j}{k} = bsxfun(@minus, Xd_dwnspl{j}{k}, mean(Xd_dwnspl{j}{k}(minindex:maxindex)));
        [dethista{j,k}, dethistb{j,k}] = hist(Xd_dwnspl{j}{k},freedmandiaconis(Xd_dwnspl{j}{k}));
        bincount = sum(dethista{j,k}); dethista{j,k} = dethista{j,k}./bincount;
        detmean(j,k) = mean(Xd_dwnspl{j}{k}); detvar(j,k) = var(Xd_dwnspl{j}{k});
    end
end

% Pre-allocate memory
hsymmdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));psymmdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
Kstatdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
KSstatdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
hsymmsto = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
pKdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
dipdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
punidet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
Xlowdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
Xupdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));


%{
% Perform three statistical tests on histograms
for k = 1:a
    for j = 1:(logdata.data(1,8))
        
        % Symmetry test [Asai '06]
        upperdet = Xd_dwnspl{j}{k}(Xd_dwnspl{j}{k} >= detmean(j,k)) - detmean(j,k);   % divide the distributions
        lowerdet = detmean(j,k) - Xd_dwnspl{j}{k}(Xd_dwnspl{j}{k} <= detmean(j,k));
        if iscolumn(upperdet) == 1
        upperdet = [upperdet' -upperdet']';lowerdet = [lowerdet' -lowerdet']'; % symmetrize
        else
        upperdet = [upperdet -upperdet]';lowerdet = [lowerdet -lowerdet]'; % symmetrize
        end
        [hsymmdet(j,k),psymmdet(j,k),KSstatdet(j,k)] = kstest2(upperdet(1,:),lowerdet(1,:),10^-2);    % kstest between the distributions
        
        % Kurtosis test, Reference: Basic Statistics for Social Research (1997)
        NK = length(Xd_dwnspl{j}{k});
        SES = sqrt(6*NK*(NK-1)/((NK-2)*(NK+1)*(NK+3))); % standard error in skewness
        SEK = 2*SES*sqrt((NK^2-1)/((NK-3)*(NK+5))); % standard error in kurtosis = std of kurtosis
        Kstatdet(j,k) = (kurtosis(Xd_dwnspl{j}{k},0)-3)/SEK; % find kurtosis in each case
        pKdet(j,k) = cdf('Normal',Kstatdet(j,k),0,1);   % find p-value for kurtosis, where Kstat is approximately normal for NK>10^3
        
        % Unimodality test
        % null is unimodality (no dip) - sort X and look for concave/convex
        % change ... dip~0.1 => not unimodal
        Nboot = 2*10^2; % number of bootstraps to find puni
        [dipdet(j,k), punidet(j,k), Xlowdet(j,k), Xupdet(j,k)]=hartigansdipsigniftest(Xd_dwnspl{j}{k},Nboot);
        
    end
end

    
% Exclude cases that go unstable or remain at zero (including those that fail
% the above statistical tests)
asymthresh = 0.0001;      % threshold for asymmetry test (KS>thresh -> asymmetric)
kurtthresh = -18;       % threshold for kurtosis (K<thresh -> fat)
hartthresh = 0.008;       % threshold for unimodality (dip>thresh -> multimodal)
%}
% Preallocate memory
detincl = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));detinclstat = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));

%{
% detinclstat/stoinclstat correspond to the number of tests that particular point
% passes (multimodal/fat/asymmetric). If zero, the point is quiescent.
% detincl/stoincl are the operating points that will be included in future
% analyses.
for k = 1:a
    for j = 1:(logdata.data(1,8))
        % Check each of the three statistical tests and include those that
        % pass.
        if KSstatdet(j,k) <= asymthresh
            detincl(j,k) = detincl(j,k) + 1;
        end
        if Kstatdet(j,k) <= kurtthresh  
            detincl(j,k) = detincl(j,k) + 1;
        end
        if dipdet(j,k) >= hartthresh
            detincl(j,k) = detincl(j,k) + 1;
        end
        detinclstat(j,k) = detincl(j,k); 
        
        % Exclude traces that go to infinity or have NaNs
        if sum(isinf(Xd_dwnspl{j}{k})) > 0 || sum(isnan(Xd_dwnspl{j}{k})) > 0
            detincl(j,k) = 0;
        end
        
        % Exclude traces that reside at zero
        if detvar(j,k) <= 10^-6
            detincl(j,k) = 0;
        end
    end
end
    
% All indices equal to 1 are included and others are not. Ensure a simply
% connected region.
detincl(detincl>=1) = 1;
for k = 1:a
    detind = find(diff(detincl(j,:))==1);
    if numel(detind)>0
        detincl(j,detind(1):end) = 1; 
    end
end
%}

% Nearest-neighbor clustering and peak finding

% Preallocate memory
c1det = cell(length(Xd_dwnspl),length(Xd_dwnspl{1}));c2det = cell(length(Xd_dwnspl),length(Xd_dwnspl{1}));
pkdet = cell(length(Xd_dwnspl),length(Xd_dwnspl{1}));trdet = cell(length(Xd_dwnspl),length(Xd_dwnspl{1}));
IEIpkdet = cell(length(Xd_dwnspl),length(Xd_dwnspl{1}));IEItrdet = cell(length(Xd_dwnspl),length(Xd_dwnspl{1}));

NNerr = 10^-6;   % error threshold for nearest-neighbor clustering
for k = 1:a
    for j = 1:(logdata.data(1,8))
        if j == 1 && k == 1
            if nnthresh == 1
                [c1det,c2det]=twoclass(Xd_dwnspl{c12thresh(1)}{c12thresh(2)},NNerr);  % nearest-neighbor clustering
                [pkdet{j,k},trdet{j,k}] = PTDetect(Xd_dwnspl{j}{k}, max([c1det c2det]));
            elseif nnthresh==2
                [pkdet{j,k},trdet{j,k}] = PTDetect(Xd_dwnspl{j}{k}, c12thresh);
            end
            for l = 2:length(pkdet{j,k})
                IEIpkdet{j,k}(l-1) = (pkdet{j,k}(l) - pkdet{j,k}(l-1))/Fs;
            end
            for l = 2:length(trdet{j,k})
                IEItrdet{j,k}(l-1) = (trdet{j,k}(l) - trdet{j,k}(l-1))/Fs;
            end
        else
            if nnthresh == 1
                [c1det,c2det]=twoclass(Xd_dwnspl{c12thresh(1)}{c12thresh(2)},NNerr);  % nearest-neighbor clustering
                [pkdet{j,k},trdet{j,k}] = PTDetect(Xd_dwnspl{j}{k}, max([c1det c2det]));
            elseif nnthresh==2
                [pkdet{j,k},trdet{j,k}] = PTDetect(Xd_dwnspl{j}{k}, c12thresh);
            end
            for l = 2:length(pkdet{j,k})
                IEIpkdet{j,k}(l-1) = (pkdet{j,k}(l) - pkdet{j,k}(l-1))/Fs;
            end
            for l = 2:length(trdet{j,k})
                IEItrdet{j,k}(l-1) = (trdet{j,k}(l) - trdet{j,k}(l-1))/Fs;
            end
            %pkdet{j,k} = 0; trdet{j,k} = 0;
            %IEIpkdet{j,k} = 0; IEItrdet{j,k} = 0;
        end
    end
end

% Pre-allocate memory
pkspikeratedet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
trspikeratedet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
CDdettime = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
CDdetpk = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
CDdettr = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
pkdiffusiondet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
trdiffusiondet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
IEIvarpkdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));IEIvartrdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
IEImeanpkdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));IEImeantrdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
pktcorrdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
trtcorrdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
meanIEIpkdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
meanIEItrdet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
pkIEIspikeratiodet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));
trIEIspikeratiodet = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}));

% Calculate statistics: coefficient of dispersion, spike rate, diffusion
% coefficient, correlation time
for k = 1:a
    for j = 1:(logdata.data(1,8))
        pkspikeratedet(j,k) = length(pkdet{j,k})/Ttotal*1e3;      % spike rate, peaks
        trspikeratedet(j,k) = length(trdet{j,k})/Ttotal*1e3;      % spike rate, troughs
        CDdettime(j,k) = detvar(j,k)/detmean(j,k);        % coefficient of dispersion
        CDdetpk(j,k) = var(IEIpkdet{j,k})/mean(IEIpkdet{j,k});
        CDdettr(j,k) = var(IEItrdet{j,k})/mean(IEItrdet{j,k});
        IEIvarpkdet(j,k) = var(IEIpkdet{j,k});IEIvartrdet(j,k) = var(IEItrdet{j,k});
        IEImeanpkdet(j,k) = mean(IEIpkdet{j,k});IEImeantrdet(j,k) = mean(IEItrdet{j,k});
        pkdiffusiondet(j,k) = 0.5*CDdetpk(j,k)^2*pkspikeratedet(j,k); % diffusion coefficient, peaks
        trdiffusiondet(j,k) = 0.5*CDdettr(j,k)^2*trspikeratedet(j,k); % diffusion coefficient, troughs
        pktcorrdet(j,k) = 0.5*pkdiffusiondet(j,k) - 1/pkspikeratedet(j,k);       % correlation time, peaks
        trtcorrdet(j,k) = 0.5*trdiffusiondet(j,k) - 1/trspikeratedet(j,k);       % correlation time, troughs
        meanIEIpkdet(j,k) = mean(IEIpkdet{j,k});meanIEItrdet(j,k) = mean(IEItrdet{j,k});    % mean IEI, should be equal to 1/(spikerate)
        pkIEIspikeratiodet(j,k) = meanIEIpkdet(j,k)*pkspikeratedet(j,k);    % ratio of IEI to 1/spikerate
        trIEIspikeratiodet(j,k) = meanIEItrdet(j,k)*trspikeratedet(j,k);
    end
end

    
% Generate data for a Poisson process by stepping through window sizes,
% sliding these windows, and counting the number of peaks within each
% window (if the peaks are given as indices, simply step through different
% number ranges and count the number of values within those index ranges.
poisswin = round(0.0001*length(Xd_dwnspl{1}{1}):(0.1*length(Xd_dwnspl{1}{1})-0.001*length(Xd_dwnspl{1}{1}))/9:0.1*length(Xd_dwnspl{1}{1}));
poisswin(poisswin<3) = 3;

clear poisscountsdetpk poisscountsstopk poisscountsdettr poisscountsstotr
% Pre-allocate memory
%cumsumdetpk = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}),length(poisswin));cumsumstopk = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}),length(poisswin));
%cumsumdettr = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}),length(poisswin));cumsumstotr = zeros(length(Xd_dwnspl),length(Xd_dwnspl{1}),length(poisswin));

for k = 1:a
    for j = 1:(logdata.data(1,8))
        Xd_dwnsplpk{j}{k} = zeros(1,length(Xd_dwnspl{1}{1}));
        Xd_dwnspltr{j}{k} = zeros(1,length(Xd_dwnspl{1}{1}));
        if length(pkdet{j,k})>1
            Xd_dwnsplpk{j}{k}(pkdet{j,k}) = 1;
        end
        if length(trdet{j,k})>1   
            Xd_dwnspltr{j}{k}(trdet{j,k}) = 1;
        end 
        cumsumdetpk(j,k,:) = cumsum(Xd_dwnsplpk{j}{k});
        cumsumdettr(j,k,:) = cumsum(Xd_dwnsplpk{j}{k});
        for l = 1:length(poisswin)
            numwindet = floor(length(Xd_dwnspl{j}{k})/poisswin(l));
            for m = 1:numwindet
                poisscountsdetpk{j,k,l}(m) = -cumsumdetpk(j,k,1+(m-1)*poisswin(l))+cumsumdetpk(j,k,m*poisswin(l));
                poisscountsdettr{j,k,l}(m) = -cumsumdettr(j,k,1+(m-1)*poisswin(l))+cumsumdettr(j,k,m*poisswin(l));
            end
            CDpoisscountsdetpk(j,k,l) = var(poisscountsdetpk{j,k,l})/mean(poisscountsdetpk{j,k,l});
            CDpoisscountsdettr(j,k,l) = var(poisscountsdettr{j,k,l})/mean(poisscountsdettr{j,k,l,:});
        end
    end
end


% Perform statistics on the distributions to see if they follow a Poisson
% process.
for k = 1:a
    for j = 1:(logdata.data(1,8))
        for l = 1:length(poisswin)
           warning off
           [poissHdetpk{j,k,l}, poissPdetpk(j,k,l), poissSTATdetpk{j,k,l}] = chi2gof(squeeze(poisscountsdetpk{j,k,l}),'cdf',@(z)poisscdf(z,mean(squeeze(poisscountsdetpk{j,k,l}))),'nparams',1);
           [poissHdettr{j,k,l}, poissPdettr(j,k,l), poissSTATdettr{j,k,l}] = chi2gof(squeeze(poisscountsdettr{j,k,l}),'cdf',@(z)poisscdf(z,mean(squeeze(poisscountsdettr{j,k,l}))),'nparams',1);
        end
    end
end

for j = 1:length(Xd_dwnspl)
    for k = 1:length(Xd_dwnspl{1})
        for l = 1:length(poisswin)
            poissdetpkchi2stat(j,k,l) = poissSTATdetpk{j,k,l}.chi2stat;
            poissdettrchi2stat(j,k,l) = poissSTATdettr{j,k,l}.chi2stat;
        end
    end
end


    
% Generate a power spectrum for each case, find the peak in the power
% spectrum, and calculate the quality of the peak (wmax/dW) and the degree
% of coherence by taking the ratio of the peak height to the width at which
% the spectrum decays to sqrt(2) of its peak (S(wmax)/[dW/wmax]).


% Calculate the time for each peak and trough
sizeP = size(pkdet);
nperc = 0.5; % percentage for threshold

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
                   trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pkdet{j,k}(l)) Xd_dwnspl{j}{k}(pkdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(trdet{j,k}(l))) + Xd_dwnspl{j}{k}(trdet{j,k}(l));
                   trtimedet{j,k}(l) = sum(Xd_dwnspl{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pkdet{j,k}(l)) Xd_dwnspl{j}{k}(pkdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(trdet{j,k}(l))) + Xd_dwnspl{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xd_dwnspl{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trdet{j,k}(l)) Xd_dwnspl{j}{k}(trdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(pkdet{j,k}(l))) + Xd_dwnspl{j}{k}(pkdet{j,k}(l));
               pktimedet{j,k}(l) = sum((Xd_dwnspl{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==loopmaxpk
               pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trdet{j,k}(l)) Xd_dwnspl{j}{k}(trdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(pkdet{j,k}(l))) + Xd_dwnspl{j}{k}(pkdet{j,k}(l));
               pktimedet{j,k}(l) = sum((Xd_dwnspl{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs; 
               end
           end
       elseif pkdet{j,k}(1) < trdet{j,k}(1) && pkdet{j,k}(end) < trdet{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pkdet{j,k}(l)) Xd_dwnspl{j}{k}(pkdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(trdet{j,k}(l))) + Xd_dwnspl{j}{k}(trdet{j,k}(l));
                   trtimedet{j,k}(l) = sum(Xd_dwnspl{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pkdet{j,k}(l)) Xd_dwnspl{j}{k}(pkdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(trdet{j,k}(l))) + Xd_dwnspl{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xd_dwnspl{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trdet{j,k}(l)) Xd_dwnspl{j}{k}(trdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(pkdet{j,k}(l))) + Xd_dwnspl{j}{k}(pkdet{j,k}(l));
               pktimedet{j,k}(l) = sum((Xd_dwnspl{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pkdet{j,k}(l)) Xd_dwnspl{j}{k}(pkdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(trdet{j,k}(l))) + Xd_dwnspl{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xd_dwnspl{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       elseif pkdet{j,k}(1) > trdet{j,k}(1) && pkdet{j,k}(end) < trdet{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trdet{j,k}(l)) Xd_dwnspl{j}{k}(trdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(pkdet{j,k}(l))) + Xd_dwnspl{j}{k}(pkdet{j,k}(l));
                   pktimedet{j,k}(l) = sum((Xd_dwnspl{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs; 
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pkdet{j,k}(l)) Xd_dwnspl{j}{k}(pkdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(trdet{j,k}(l))) + Xd_dwnspl{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xd_dwnspl{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trdet{j,k}(l)) Xd_dwnspl{j}{k}(trdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(pkdet{j,k}(l))) + Xd_dwnspl{j}{k}(pkdet{j,k}(l));
               pktimedet{j,k}(l) = sum((Xd_dwnspl{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs;
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pkdet{j,k}(l)) Xd_dwnspl{j}{k}(pkdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(trdet{j,k}(l))) + Xd_dwnspl{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xd_dwnspl{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               end
           end
       elseif pkdet{j,k}(1) >= trdet{j,k}(1) && pkdet{j,k}(end) >= trdet{j,k}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trdet{j,k}(l)) Xd_dwnspl{j}{k}(trdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(pkdet{j,k}(l))) + Xd_dwnspl{j}{k}(pkdet{j,k}(l));
                   pktimedet{j,k}(l) = sum((Xd_dwnspl{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pkdet{j,k}(l)) Xd_dwnspl{j}{k}(pkdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(trdet{j,k}(l))) + Xd_dwnspl{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xd_dwnspl{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xd_dwnspl{j}{k}(trdet{j,k}(l)) Xd_dwnspl{j}{k}(trdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(pkdet{j,k}(l))) + Xd_dwnspl{j}{k}(pkdet{j,k}(l));
               pktimedet{j,k}(l) = sum((Xd_dwnspl{j}{k}(trdet{j,k}(l):trdet{j,k}(l+1))>=pkthresh))/Fs; 
               end
               if l==min([loopmaxpk loopmaxtr])-1
               trthresh = nperc*(mean([Xd_dwnspl{j}{k}(pkdet{j,k}(l)) Xd_dwnspl{j}{k}(pkdet{j,k}(l+1))])-Xd_dwnspl{j}{k}(trdet{j,k}(l))) + Xd_dwnspl{j}{k}(trdet{j,k}(l));
               trtimedet{j,k}(l) = sum(Xd_dwnspl{j}{k}(pkdet{j,k}(l):pkdet{j,k}(l+1))<=trthresh)/Fs;
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
        meantrtimedet(j,k) = mean(trtimedet{j,k});
        meanpktimedet(j,k) = mean(pktimedet{j,k});
        CVtrtimedet(j,k) = std(trtimedet{j,k})/mean(trtimedet{j,k});
        CVpktimedet(j,k) = std(pktimedet{j,k})/mean(pktimedet{j,k});
        CDtrtimedet(j,k) = var(trtimedet{j,k})/mean(trtimedet{j,k});
        CDpktimedet(j,k) = var(pktimedet{j,k})/mean(pktimedet{j,k});
    end
end

% Downsample data
for j = 1:a
    for i = 1:(logdata.data(1,8))
        Xd_dwnspl2{i}{j} = Xd_dwnspl{i}{j}.*Xscale2(i,j);
    end
end
% Find amplitudes using peaks
for j = 1:sizeP(1)
    for k = 1:sizeP(2)
        if isempty(pkdet{j,k}) == 0 && isempty(trdet{j,k}) == 0
            loopmaxtr = length(trdet{j,k});loopmaxpk = length(pkdet{j,k});
            if min([loopmaxpk loopmaxtr])>2
            for l = 1:min([loopmaxpk loopmaxtr])-1
                ampldet{j,k}(l) = Xd_dwnspl2{j}{k}(pkdet{j,k}(l)) - Xd_dwnspl2{j}{k}(trdet{j,k}(l));
            end
            else
                ampldet{j,k}=0;
            end
        else
            ampldet{j,k}=0;
        end
        PKamplmeandet(j,k) = mean(ampldet{j,k});PKamplsemdet(j,k)=std(ampldet{j,k})/sqrt(length(ampldet{j,k}));
    end
end


% Find amplitudes and frequencies using FFT
T = 1/Fs;
L = length(Xd_dwnspl2{1}{1});
NFFT = (2^4)*2^nextpow2(L);
nw=10;
XsegL = floor(length(Xd_dwnspl2{1}{1})/nw);
welchwin = round(XsegL);
NPSD = floor(NFFT/nw);
noverlap = 0;
winfunc = hamming(welchwin);
freq = 0.005;
Xsine = sin(2*pi*freq.*tvec);
[Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
for j = 1:length(Xd_dwnspl2)
    for k = 1:length(Xd_dwnspl2{1})
        clear Xstoscaled Xd_dwnsplscaled
        [Xd_dwnsplfft{j,k}, fdetfft{j,k}] = pwelch(Xd_dwnspl2{j}{k},winfunc,noverlap,NPSD,Fs);
        fdetind = find(fdetfft{j,k} > 0.0001);
        fscale = 10^3;
        Xd_dwnsplscaled = Xd_dwnsplfft{j,k}./fscale;
        if max(2*abs(Xd_dwnsplfft{j,k})) > 1e-2
        Xd_dwnsplfftmaxind = find(Xd_dwnsplfft{j,k}(fdetind)==max(Xd_dwnsplfft{j,k}(fdetind)));
        fftampldet(j,k) = Xd_dwnsplfft{j,k}(fdetind(Xd_dwnsplfftmaxind(1)));
        fftampldet(j,k) = (sqrt(fscale.*fftampldet(j,k).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
        fftfreqdet(j,k) = fdetfft{j,k}(fdetind(Xd_dwnsplfftmaxind(1)));
        else
        fftfreqdet(j,k) = 0;
        fftampldet(j,k) = 0;  
        end
    end
end

if saveyn == 1
%savefile=file;  % Overwrite old file
save(sprintf('%s%s%s%s%s%s',filepath,'Dwnspl',num2str(dwnspl),'-Thresh',num2str(c12thresh),'-output.mat'),'-v7.3');
end
end
