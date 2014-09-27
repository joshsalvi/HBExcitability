%% (1) DATA IMPORT

% Import data
%file='/Users/joshsalvi/Desktop/SNICstochoutput.mat';
%load(file);

% Type of bifurcation
biftype = 2;        % 1=supercritical Hopf; 2=SNIC; 3=subcritical Hopf
if biftype == 1
    Xdet1 = Hopfdet; Xsto1 = Hopfsto;
    clear i;
    for j = 1:length(Xdet)
        for k = 1:length(Xdet{1})
            Xdet{j,k} = Xdet1{j,k}(1,:) + i*Xdet1{j,k}(2,:);
            Xsto{j,k} = Xsto1{j,k}(1,:) + i*Xsto1{j,k}(2,:);
        end
    end
    clear Xdet1 Xsto1
elseif biftype == 2
    Xdet = SNICdet; Xsto = SNICsto;
elseif biftype == 3
    Xdet1 = Hopfsubdet; Xsto1 = Hopfsubsto;
    for j = 1:length(Xdet)
        for k = 1:length(Xdet{1})
            Xdet{j,k} = Xdet1{j,k}(1,:).*sin(Xdet1{j,k}(2,:));
            Xsto{j,k} = Xsto1{j,k}(1,:).*sin(Xsto1{j,k}(2,:));
        end
    end
    clear Xdet1 Xsto1
else
    disp('No bifurcation type chosen');
end

% Define time vector
tvec=t;
% Define the sample rate if not already done
Fs=1/(tvec(2)-tvec(1));

analysisstep=1;
%% (2) HISTOGRAM ANALYSIS

if analysisstep == 1
% Time range
tmin=tvec(0.1*length(tvec));
tmax=tvec(end);
minindex = find(abs(tvec-tmin)==min(abs(tvec-tmin)));
maxindex = find(abs(tvec-tmax)==min(abs(tvec-tmax)))-1;
Ttotal = tmax-tmin;

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
detmean = zeros(length(Xdet),length(Xdet{1}));
stomean = zeros(length(Xdet),length(Xdet{1}));
dethista = cell(length(Xdet),length(Xdet{1})); dethistb = cell(length(Xdet),length(Xdet{1}));
stohista = cell(length(Xdet),length(Xdet{1})); stohistb = cell(length(Xdet),length(Xdet{1}));
detvar = zeros(length(Xdet),length(Xdet{1})); stovar = zeros(length(Xdet),length(Xdet{1}));


% Generate histograms and calculate statistics
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        Xdet{j}{k} = bsxfun(@minus, Xdet{j}{k}, mean(Xdet{j}{k}(minindex:maxindex)));
        Xsto{j}{k} = bsxfun(@minus, Xsto{j}{k}, mean(Xsto{j}{k}(minindex:maxindex)));
        [dethista{j,k}, dethistb{j,k}] = hist(Xdet{j}{k},freedmandiaconis(Xdet{j}{k}));
        bincount = sum(dethista{j,k}); dethista{j,k} = dethista{j,k}./bincount;
        bincount = sum(stohista{j,k}); stohista{j,k} = stohista{j,k}./bincount;
        [stohista{j,k}, stohistb{j,k}] = hist(Xsto{j}{k},freedmandiaconis(Xsto{j}{k}));
        detmean(j,k) = mean(Xdet{j}{k}); detvar(j,k) = var(Xdet{j}{k});
        stomean(j,k) = mean(Xsto{j}{k}); stovar(j,k) = var(Xsto{j}{k});
    end
end

% Pre-allocate memory
hsymmdet = zeros(length(Xdet),length(Xdet{1}));psymmdet = zeros(length(Xdet),length(Xdet{1}));
Kstatdet = zeros(length(Xdet),length(Xdet{1}));Kstatsto = zeros(length(Xdet),length(Xdet{1}));
KSstatdet = zeros(length(Xdet),length(Xdet{1}));KSstatsto = zeros(length(Xdet),length(Xdet{1}));
hsymmsto = zeros(length(Xdet),length(Xdet{1}));psymmsto = zeros(length(Xdet),length(Xdet{1}));
pKdet = zeros(length(Xdet),length(Xdet{1}));pKsto = zeros(length(Xdet),length(Xdet{1}));
dipdet = zeros(length(Xdet),length(Xdet{1}));dipsto = zeros(length(Xdet),length(Xdet{1}));
punidet = zeros(length(Xdet),length(Xdet{1}));punisto = zeros(length(Xdet),length(Xdet{1}));
Xlowdet = zeros(length(Xdet),length(Xdet{1}));Xlowsto = zeros(length(Xdet),length(Xdet{1}));
Xupdet = zeros(length(Xdet),length(Xdet{1}));Xupsto = zeros(length(Xdet),length(Xdet{1}));

% Perform three statistical tests on histograms
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        
        % Symmetry test [Asai '06]
        upperdet = Xdet{j}{k}(Xdet{j}{k} >= detmean(j,k)) - detmean(j,k);   % divide the distributions
        lowerdet = detmean(j,k) - Xdet{j}{k}(Xdet{j}{k} <= detmean(j,k));
        uppersto = Xsto{j}{k}(Xsto{j}{k} >= stomean(j,k)) - stomean(j,k);   
        lowersto = stomean(j,k) - Xsto{j}{k}(Xsto{j}{k} <= stomean(j,k));
        if iscolumn(upperdet) == 1
        upperdet = [upperdet' -upperdet']';lowerdet = [lowerdet' -lowerdet']'; % symmetrize
        uppersto = [uppersto' -uppersto']';lowersto = [lowersto' -lowersto']';
        else
        upperdet = [upperdet -upperdet]';lowerdet = [lowerdet -lowerdet]'; % symmetrize
        uppersto = [uppersto -uppersto]';lowersto = [lowersto -lowersto]';
        end
        [hsymmdet(j,k),psymmdet(j,k),KSstatdet(j,k)] = kstest2(upperdet,lowerdet,10^-2);    % kstest between the distributions
        [hsymmsto(j,k),psymmsto(j,k),KSstatsto(j,k)] = kstest2(uppersto,lowersto,10^-2);
        
        % Kurtosis test, Reference: Basic Statistics for Social Research (1997)
        NK = length(Xdet{j}{k});
        SES = sqrt(6*NK*(NK-1)/((NK-2)*(NK+1)*(NK+3))); % standard error in skewness
        SEK = 2*SES*sqrt((NK^2-1)/((NK-3)*(NK+5))); % standard error in kurtosis = std of kurtosis
        Kstatdet(j,k) = (kurtosis(Xdet{j}{k},0)-3)/SEK; % find kurtosis in each case
        Kstatsto(j,k) = (kurtosis(Xsto{j}{k},0)-3)/SEK;
        pKdet(j,k) = cdf('Normal',Kstatdet(j,k),0,1);   % find p-value for kurtosis, where Kstat is approximately normal for NK>10^3
        pKsto(j,k) = cdf('Normal',Kstatsto(j,k),0,1);
        
        % Unimodality test
        % null is unimodality (no dip) - sort X and look for concave/convex
        % change ... dip~0.1 => not unimodal
        Nboot = 2*10^2; % number of bootstraps to find puni
        [dipdet(j,k), punidet(j,k), Xlowdet(j,k), Xupdet(j,k)]=hartigansdipsigniftest(Xdet{j}{k},Nboot);
        [dipsto(j,k), punisto(j,k), Xlowsto(j,k), Xupsto(j,k)]=hartigansdipsigniftest(Xsto{j}{k},Nboot);
        
    end
end
analysisstep = 2;
else
    disp('Run previous cell.');
end


%% (3) EXCLUDE SELECTED EXAMPLES        
if analysisstep == 2
    
% Exclude cases that go unstable or remain at zero (including those that fail
% the above statistical tests)
asymthresh = 0.06;      % threshold for asymmetry test (KS>thresh -> asymmetric)
kurtthresh = -18;       % threshold for kurtosis (K<thresh -> fat)
hartthresh = 0.02;       % threshold for unimodality (dip>thresh -> multimodal)

% Preallocate memory
detincl = zeros(length(Xdet),length(Xdet{1}));detinclstat = zeros(length(Xdet),length(Xdet{1}));
stoincl = zeros(length(Xdet),length(Xdet{1}));stoinclstat = zeros(length(Xdet),length(Xdet{1}));

% detinclstat/stoinclstat correspond to the number of tests that particular point
% passes (multimodal/fat/asymmetric). If zero, the point is quiescent.
% detincl/stoincl are the operating points that will be included in future
% analyses.
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        % Check each of the three statistical tests and include those that
        % pass.
        if KSstatdet(j,k) >= asymthresh(j,k) && psymmdet(j,k) <= 10^-2
            detincl(j,k) = detincl(j,k) + 1;
        end
        if KSstatsto(j,k) >= asymthresh(j,k) && psymmsto(j,k) <= 10^-2
            stoincl(j,k) = stoincl(j,k) + 1;
        end
        if Kstatdet(j,k) <= kurtthresh && pKdet(j,k) <= 10^-3   
            detincl(j,k) = detincl(j,k) + 1;
        end
        if Kstatsto(j,k) <= kurtthresh && pKsto(j,k) <= 10^-3   
            stoincl(j,k) = stoincl(j,k) + 1;
        end
        if dipdet(j,k) >= hartthresh && punidet(j,k) <= 10^-4
            detincl(j,k) = detincl(j,k) + 1;
        end
        if dipsto(j,k) >= hartthresh && punisto(j,k) <= 10^-4
            stoincl(j,k) = stoincl(j,k) + 1;
        end
        detinclstat(j,k) = detincl(j,k); stoinclstat(j,k) = stoincl(j,k);
        
        % Exclude traces that go to infinity or have NaNs
        if sum(isinf(Xdet{j}{k})) > 0 || sum(isnan(Xdet{j}{k})) > 0
            detincl(j,k) = 0;
        end
        if sum(isinf(Xsto{j}{k})) > 0 || sum(isnan(Xsto{j}{k})) > 0
            stoincl(j,k) = 0;
        end 
        
        % Exclude traces that reside at zero
        if detvar(j,k) <= 10^-6
            detincl(j,k) = 0;
        end
        if stovar(j,k) <= 10^-6
            stoincl(j,k) = 0;
        end
        
    end
end
    
detincl(detincl>1) = 1;
stoincl(stoincl>1) = 1;

analysisstep = 3;
else
    disp('Run previous cell.');
end
%% (4) PEAK-FINDING ALGORITHM
if analysisstep == 3
    
% Nearest-neighbor clustering and peak finding

% Preallocate memory
c1det = cell(length(Xdet),length(Xdet{1}));c2det = cell(length(Xdet),length(Xdet{1}));
pkdet = cell(length(Xdet),length(Xdet{1}));trdet = cell(length(Xdet),length(Xdet{1}));
c1sto = cell(length(Xdet),length(Xdet{1}));c2sto = cell(length(Xdet),length(Xdet{1}));
pksto = cell(length(Xdet),length(Xdet{1}));trsto = cell(length(Xdet),length(Xdet{1}));
IEIpkdet = cell(length(Xdet),length(Xdet{1}));IEItrdet = cell(length(Xdet),length(Xdet{1}));
IEIpksto = cell(length(Xdet),length(Xdet{1}));IEItrsto = cell(length(Xdet),length(Xdet{1}));

NNerr = 10^-6   % error threshold for nearest-neighbor clustering
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        if detincl(j,k) == 1
            [c1det(j,k),c2det(j,k)]=twoclass(Xdet{j}{k},NNerr);  % nearest-neighbor clustering
            [pkdet{j,k},trdet{j,k}] = PTDetect(Xdet{j}{k}, max([c1det{j,k} c2det{j,k}]));
            for l = 2:length(pkdet{j,k})
                IEIpkdet{j,k}(l-1) = (pkdet{j,k}(l) - pkdet{j,k}(l-1))/Fs;
            end
            for l = 2:length(trdet{j,k})
                IEItrdet{j,k}(l-1) = (trdet{j,k}(l) - trdet{j,k}(l-1))/Fs;
            end
        else
            pkdet{j,k} = 0; trdet{j,k} = 0;
            IEIpkdet{j,k} = 0; IEItrdet{j,k} = 0;
        end
        if stoincl(j,k) == 1
            [c1sto(j,k),c2sto(j,k)]=twoclass(Xsto{j}{k},NNerr);  % nearest-neighbor clustering
            [pksto{j,k},trsto{j,k}] = PTDetect(Xsto{j}{k}, max([c1sto{j,k} c2sto{j,k}]));
            for l = 2:length(pksto{j,k})
                IEIpksto{j,k}(l-1) = (pksto{j,k}(l) - pksto{j,k}(l-1))/Fs;
            end
            for l = 2:length(trsto{j,k})
                IEItrsto{j,k}(l-1) = (trsto{j,k}(l) - trsto{j,k}(l-1))/Fs;
            end
        else
            pksto{j,k} = 0; trsto{j,k} = 0;
            IEIpksto{j,k} = 0; IEItrsto{j,k} = 0;
        end
    end
end

% Pre-allocate memory
pkspikeratedet = zeros(length(Xdet),length(Xdet{1}));pkspikeratesto = zeros(length(Xdet),length(Xdet{1}));
trspikeratedet = zeros(length(Xdet),length(Xdet{1}));trspikeratesto = zeros(length(Xdet),length(Xdet{1}));
CDdet = zeros(length(Xdet),length(Xdet{1}));CDsto = zeros(length(Xdet),length(Xdet{1}));
pkdiffusiondet = zeros(length(Xdet),length(Xdet{1}));pkdiffusionsto = zeros(length(Xdet),length(Xdet{1}));
trdiffusiondet = zeros(length(Xdet),length(Xdet{1}));trdiffusionsto = zeros(length(Xdet),length(Xdet{1}));
pktcorrdet = zeros(length(Xdet),length(Xdet{1}));pktcorrsto = zeros(length(Xdet),length(Xdet{1}));
trtcorrdet = zeros(length(Xdet),length(Xdet{1}));trtcorrsto = zeros(length(Xdet),length(Xdet{1}));
meanIEIpkdet = zeros(length(Xdet),length(Xdet{1}));meanIEItrdet = zeros(length(Xdet),length(Xdet{1}));
meanIEIpksto = zeros(length(Xdet),length(Xdet{1}));meanIEItrsto = zeros(length(Xdet),length(Xdet{1}));
pkIEIspikeratiodet = zeros(length(Xdet),length(Xdet{1}));pkIEIspikeratiosto = zeros(length(Xdet),length(Xdet{1}));
trIEIspikeratiodet = zeros(length(Xdet),length(Xdet{1}));trIEIspikeratiosto = zeros(length(Xdet),length(Xdet{1}));

% Calculate statistics: coefficient of dispersion, spike rate, diffusion
% coefficient, correlation time
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        pkspikeratedet(j,k) = length(pkdet{j,k})/Ttotal;      % spike rate, peaks
        pkspikeratesto(j,k) = length(pksto{j,k})/Ttotal;
        trspikeratedet(j,k) = length(trdet{j,k})/Ttotal;      % spike rate, troughs
        trspikeratesto(j,k) = length(trsto{j,k})/Ttotal;
        CDdet(j,k) = detvar(j,k)/detmean(j,k);        % coefficient of dispersion
        CDsto(j,k) = stovar(j,k)/stomean(j,k);
        pkdiffusiondet(j,k) = 0.5*CDdet(j,k)^2*pkspikeratedet(j,k); % diffusion coefficient, peaks
        pkdiffusionsto(j,k) = 0.5*CDsto(j,k)^2*pkspikeratesto(j,k);
        trdiffusiondet(j,k) = 0.5*CDdet(j,k)^2*trspikeratedet(j,k); % diffusion coefficient, troughs
        trdiffusionsto(j,k) = 0.5*CDsto(j,k)^2*trspikeratesto(j,k);
        pktcorrdet(j,k) = 0.5*pkdiffusiondet(j,k) - 1/pkspikeratedet(j,k);       % correlation time, peaks
        pktcorrsto(j,k) = 0.5*pkdiffusionsto(j,k) - 1/pkspikeratesto(j,k);
        trtcorrdet(j,k) = 0.5*trdiffusiondet(j,k) - 1/trspikeratedet(j,k);       % correlation time, troughs
        trtcorrsto(j,k) = 0.5*trdiffusionsto(j,k) - 1/trspikeratesto(j,k);
        meanIEIpkdet(j,k) = mean(IEIpkdet{j,k});meanIEItrdet(j,k) = mean(IEItrdet{j,k});    % mean IEI, should be equal to 1/(spikerate)
        meanIEIpksto(j,k) = mean(IEIpksto{j,k});meanIEItrsto(j,k) = mean(IEItrsto{j,k});    
        pkIEIspikeratiodet(j,k) = meanIEIpkdet(j,k)*pkspikeratedet(j,k);    % ratio of IEI to 1/spikerate
        trIEIspikeratiodet(j,k) = meanIEItrdet(j,k)*trspikeratedet(j,k);
        pkIEIspikeratiosto(j,k) = meanIEIpksto(j,k)*pkspikeratesto(j,k);    % ratio of IEI to 1/spikerate
        trIEIspikeratiosto(j,k) = meanIEItrsto(j,k)*trspikeratesto(j,k);
    end
end

analysisstep = 4;
else
    disp('Run previous cell.');
end
%% (5) POISSON PROCESS STATISTICS
if analysisstep == 4
    
% Generate data for a Poisson process by stepping through window sizes,
% sliding these windows, and counting the number of peaks within each
% window (if the peaks are given as indices, simply step through different
% number ranges and count the number of values within those index ranges.

% Generate distributions from the above method.

% Perform statistics on the distributions to see if they follow a Poisson
% process.

analysisstep = 5;
else
    disp('Run previous cell.');
end
%% (6) POWER SPECTRA STATISTICS
if analysisstep == 5
    
% Generate a power spectrum for each case, find the peak in the power
% spectrum, and calculate the quality of the peak (wmax/dW) and the degree
% of coherence by taking the ratio of the peak height to the width at which
% the spectrum decays to sqrt(2) of its peak (S(wmax)/[dW/wmax]).

analysisstep = 6;
else
    disp('Run previous cell.');
end

%% (X) SAVE DATA

savefile='';    % New file
savefile=file;  % Overwrite old file
save(savefile);
