%% (1) DATA IMPORT

% Downsample rate
dwnspl = 10;


% Import data
%file='/Users/joshsalvi/Desktop/SNICstochoutput.mat';
%load(file);

% Type of bifurcation
biftype = 2;        % 1=supercritical Hopf; 2=SNIC; 3=subcritical Hopf
if biftype == 1
    Xdet1 = Hopfdet; Xsto1 = Hopfsto;
    clear i;
    for j = 1:length(Xdet1)
        for k = 1:length(Xdet1{1})
            Xdet{j}{k} = Xdet1{j}{k}(1,1:dwnspl:end) + 1*Xdet1{j}{k}(2,1:dwnspl:end);
            Xsto{j}{k} = Xsto1{j}{k}(1,1:dwnspl:end) + 1*Xsto1{j}{k}(2,1:dwnspl:end);
        end
    end
    clear Xdet1 Xsto1 Hopfdet Hopfsto
elseif biftype == 2
    for j = 1:length(SNICdet)
        for k = 1:length(SNICdet{1})
            Xdet{j}{k} = SNICdet{j}{k}(1:dwnspl:end);
            Xsto{j}{k} = SNICsto{j}{k}(1:dwnspl:end);
        end
    end
elseif biftype == 3
    Xdet1 = Hopfsubdet; Xsto1 = Hopfsubsto;
    for j = 1:length(Xdet1)
        for k = 1:length(Xdet1{1})
            Xdet{j}{k} = Xdet1{j}{k}(1,1:dwnspl:end).*sin(Xdet1{j}{k}(2,1:dwnspl:end));
            Xsto{j}{k} = Xsto1{j}{k}(1,1:dwnspl:end).*sin(Xsto1{j}{k}(2,1:dwnspl:end));
        end
    end
    clear Xdet1 Xsto1
else
    disp('No bifurcation type chosen');
end

% Define time vector
tvec=t(1:dwnspl:end);
% Define the sample rate if not already done
Fs=1/(tvec(2)-tvec(1));

analysisstep=1;

%% (1B) OPTIONAL: PARTITION DATA


% Separate the data into N traces.
% Repeat this process manually for all traces 1:N
Xdet3=Xdet;Xsto3=Xsto;

N=5;    % divide into N traces
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        for l = 1:N
            Xdet2{l}{j}{k} = Xdet3{j}{k}(round(1+(l-1)*length(Xdet{1}{1})/N):round((l)*length(Xdet{1}{1})/N));
            Xsto2{l}{j}{k} = Xsto3{j}{k}(round(1+(l-1)*length(Xsto{1}{1})/N):round((l)*length(Xsto{1}{1})/N));
        end
    end
end
%clear Xdet Xsto

% Which sample would you like to analyze first?
sample = 5;
Xdet = Xdet2{sample};
Xsto = Xsto2{sample};

tvec = tvec(1:length(Xdet{1}));

%clear Xdet2 Xsto2 Xdet3 Xsto3


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
        [hsymmdet(j,k),psymmdet(j,k),KSstatdet(j,k)] = kstest2(upperdet(1,:),lowerdet(1,:),10^-2);    % kstest between the distributions
        [hsymmsto(j,k),psymmsto(j,k),KSstatsto(j,k)] = kstest2(uppersto(1,:),lowersto(1,:),10^-2);
        
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
analysisstep=2;
if analysisstep == 2
    
% Exclude cases that go unstable or remain at zero (including those that fail
% the above statistical tests)
asymthresh = 0.001;      % threshold for asymmetry test (KS>thresh -> asymmetric)
kurtthresh = -18;       % threshold for kurtosis (K<thresh -> fat)
hartthresh = 0.005;       % threshold for unimodality (dip>thresh -> multimodal)

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
        if KSstatdet(j,k) <= asymthresh
            detincl(j,k) = detincl(j,k) + 1;
        end
        if KSstatsto(j,k) <= asymthresh
            stoincl(j,k) = stoincl(j,k) + 1;
        end
        if Kstatdet(j,k) <= kurtthresh  
            detincl(j,k) = detincl(j,k) + 1;
        end
        if Kstatsto(j,k) <= kurtthresh  
            stoincl(j,k) = stoincl(j,k) + 1;
        end
        if dipdet(j,k) >= hartthresh
            detincl(j,k) = detincl(j,k) + 1;
        end
        if dipsto(j,k) >= hartthresh
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
    
% All indices equal to 1 are included and others are not. Ensure a simply
% connected region.
detincl(detincl>=1) = 1;
stoincl(stoincl>=1) = 1;
for j = 1:length(Xdet)
    detind = find(diff(detincl(j,:))==1); stoind = find(diff(stoincl(j,:))==1);
    if numel(detind)>0
        detincl(j,detind(1):end) = 1; 
    end
    if numel(stoind)>0
        stoincl(j,stoind(1):end) = 1;
    end
end

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

NNerr = 10^-6;   % error threshold for nearest-neighbor clustering
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        if detincl(j,k) == 1
            [c1det{j,k},c2det{j,k}]=twoclass(Xdet{j}{k},NNerr);  % nearest-neighbor clustering
            [pkdet{j,k},trdet{j,k}] = PTDetect(Xdet{j}{k}, max([c1det{j,k} c2det{j,k}]));
            for l = 2:length(pkdet{j,k})
                IEIpkdet{j,k}(l-1) = (pkdet{j,k}(l) - pkdet{j,k}(l-1))/Fs;
            end
            for l = 2:length(trdet{j,k})
                IEItrdet{j,k}(l-1) = (trdet{j,k}(l) - trdet{j,k}(l-1))/Fs;
            end
        else
            [c1det{j,k},c2det{j,k}]=twoclass(Xdet{j}{k},NNerr);  % nearest-neighbor clustering
            [pkdet{j,k},trdet{j,k}] = PTDetect(Xdet{j}{k}, max([c1det{j,k} c2det{j,k}]));
            for l = 2:length(pkdet{j,k})
                IEIpkdet{j,k}(l-1) = (pkdet{j,k}(l) - pkdet{j,k}(l-1))/Fs;
            end
            for l = 2:length(trdet{j,k})
                IEItrdet{j,k}(l-1) = (trdet{j,k}(l) - trdet{j,k}(l-1))/Fs;
            end
            %pkdet{j,k} = 0; trdet{j,k} = 0;
            %IEIpkdet{j,k} = 0; IEItrdet{j,k} = 0;
        end
        if stoincl(j,k) == 1
            [c1sto{j,k},c2sto{j,k}]=twoclass(Xsto{j}{k},NNerr);  % nearest-neighbor clustering
            [pksto{j,k},trsto{j,k}] = PTDetect(Xsto{j}{k}, max([c1sto{j,k} c2sto{j,k}]));
            for l = 2:length(pksto{j,k})
                IEIpksto{j,k}(l-1) = (pksto{j,k}(l) - pksto{j,k}(l-1))/Fs;
            end
            for l = 2:length(trsto{j,k})
                IEItrsto{j,k}(l-1) = (trsto{j,k}(l) - trsto{j,k}(l-1))/Fs;
            end
        else
            %pksto{j,k} = 0; trsto{j,k} = 0;
            %IEIpksto{j,k} = 0; IEItrsto{j,k} = 0;
            [c1sto{j,k},c2sto{j,k}]=twoclass(Xsto{j}{k},NNerr);  % nearest-neighbor clustering
            [pksto{j,k},trsto{j,k}] = PTDetect(Xsto{j}{k}, max([c1sto{j,k} c2sto{j,k}]));
            for l = 2:length(pksto{j,k})
                IEIpksto{j,k}(l-1) = (pksto{j,k}(l) - pksto{j,k}(l-1))/Fs;
            end
            for l = 2:length(trsto{j,k})
                IEItrsto{j,k}(l-1) = (trsto{j,k}(l) - trsto{j,k}(l-1))/Fs;
            end
        end
    end
end

% Pre-allocate memory
pkspikeratedet = zeros(length(Xdet),length(Xdet{1}));pkspikeratesto = zeros(length(Xdet),length(Xdet{1}));
trspikeratedet = zeros(length(Xdet),length(Xdet{1}));trspikeratesto = zeros(length(Xdet),length(Xdet{1}));
CDdettime = zeros(length(Xdet),length(Xdet{1}));CDstotime = zeros(length(Xdet),length(Xdet{1}));
CDdetpk = zeros(length(Xdet),length(Xdet{1}));CDstopk = zeros(length(Xdet),length(Xdet{1}));
CDdettr = zeros(length(Xdet),length(Xdet{1}));CDstotr = zeros(length(Xdet),length(Xdet{1}));
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
        CDdettime(j,k) = detvar(j,k)/detmean(j,k);        % coefficient of dispersion
        CDstotime(j,k) = stovar(j,k)/stomean(j,k);
        CDdetpk(j,k) = var(IEIpkdet{j,k})/mean(IEIpkdet{j,k});
        CDdettr(j,k) = var(IEItrdet{j,k})/mean(IEItrdet{j,k});
        CDstopk(j,k) = var(IEIpksto{j,k})/mean(IEIpksto{j,k});
        CDstotr(j,k) = var(IEItrsto{j,k})/mean(IEItrsto{j,k});
        pkdiffusiondet(j,k) = 0.5*CDdetpk(j,k)^2*pkspikeratedet(j,k); % diffusion coefficient, peaks
        pkdiffusionsto(j,k) = 0.5*CDstopk(j,k)^2*pkspikeratesto(j,k);
        trdiffusiondet(j,k) = 0.5*CDdettr(j,k)^2*trspikeratedet(j,k); % diffusion coefficient, troughs
        trdiffusionsto(j,k) = 0.5*CDstotr(j,k)^2*trspikeratesto(j,k);
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
poisswin = round(0.0001*length(Xdet{1}{1}):(0.1*length(Xdet{1}{1})-0.001*length(Xdet{1}{1}))/9:0.1*length(Xdet{1}{1}));
poisswin(poisswin<3) = 3;


% Pre-allocate memory
%cumsumdetpk = zeros(length(Xdet),length(Xdet{1}),length(poisswin));cumsumstopk = zeros(length(Xdet),length(Xdet{1}),length(poisswin));
%cumsumdettr = zeros(length(Xdet),length(Xdet{1}),length(poisswin));cumsumstotr = zeros(length(Xdet),length(Xdet{1}),length(poisswin));

for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        Xdetpk{j}{k} = zeros(1,length(Xdet{1}{1}));Xstopk{j}{k} = zeros(1,length(Xdet{1}{1}));
        Xdettr{j}{k} = zeros(1,length(Xdet{1}{1}));Xstotr{j}{k} = zeros(1,length(Xdet{1}{1}));
        if length(pkdet{j,k})>1
            Xdetpk{j}{k}(pkdet{j,k}) = 1;
        end
        if length(trdet{j,k})>1   
            Xdettr{j}{k}(trdet{j,k}) = 1;
        end 
        if length(pksto{j,k})>1
            Xstopk{j}{k}(pksto{j,k}) = 1;
        end
        if length(trsto{j,k})>1   
            Xstotr{j}{k}(trsto{j,k}) = 1;
        end
        cumsumdetpk(j,k,:) = cumsum(Xdetpk{j}{k});
        cumsumstopk(j,k,:) = cumsum(Xstopk{j}{k});
        cumsumdettr(j,k,:) = cumsum(Xdetpk{j}{k});
        cumsumstotr(j,k,:) = cumsum(Xstotr{j}{k});
        for l = 1:length(poisswin)
            numwindet = floor(length(Xdet{j}{k})/poisswin(l));
            numwinsto = floor(length(Xsto{j}{k})/poisswin(l));
            for m = 1:numwindet
                poisscountsdetpk(j,k,l,m) = -cumsumdetpk(j,k,1+(m-1)*poisswin(l))+cumsumdetpk(j,k,m*poisswin(l));
                poisscountsdettr(j,k,l,m) = -cumsumdettr(j,k,1+(m-1)*poisswin(l))+cumsumdettr(j,k,m*poisswin(l));
            end
            for m = 2:numwinsto
                poisscountsstopk(j,k,l,m) = -cumsumstopk(j,k,1+(m-1)*poisswin(l))+cumsumstopk(j,k,m*poisswin(l));
                poisscountsstotr(j,k,l,m) = -cumsumstotr(j,k,1+(m-1)*poisswin(l))+cumsumstotr(j,k,m*poisswin(l));                       
            end
        end
    end
end


% Perform statistics on the distributions to see if they follow a Poisson
% process.
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        for l = 1:length(poisswin)
           warning off
           [poissHdetpk{j,k,l}, poissPdetpk(j,k,l), poissSTATdetpk{j,k,l}] = chi2gof(squeeze(poisscountsdetpk(j,k,l,:)),'cdf',@(z)poisscdf(z,mean(squeeze(poisscountsdetpk(j,k,l,:)))),'nparams',1);
           [poissHdettr{j,k,l}, poissPdettr(j,k,l), poissSTATdettr{j,k,l}] = chi2gof(squeeze(poisscountsdettr(j,k,l,:)),'cdf',@(z)poisscdf(z,mean(squeeze(poisscountsdettr(j,k,l,:)))),'nparams',1);
           [poissHstopk{j,k,l}, poissPstopk(j,k,l), poissSTATstopk{j,k,l}] = chi2gof(squeeze(poisscountsstopk(j,k,l,:)),'cdf',@(z)poisscdf(z,mean(squeeze(poisscountsstopk(j,k,l,:)))),'nparams',1);
           [poissHstotr{j,k,l}, poissPstotr(j,k,l), poissSTATstotr{j,k,l}] = chi2gof(squeeze(poisscountsstotr(j,k,l,:)),'cdf',@(z)poisscdf(z,mean(squeeze(poisscountsstotr(j,k,l,:)))),'nparams',1);        
        end
    end
end

for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        for l = 1:length(poisswin)
            poissdetpkchi2stat(j,k,l) = poissSTATdetpk{j,k,l}.chi2stat;
            poissdettrchi2stat(j,k,l) = poissSTATdettr{j,k,l}.chi2stat;
            poissstopkchi2stat(j,k,l) = poissSTATstopk{j,k,l}.chi2stat;
            poissstotrchi2stat(j,k,l) = poissSTATstotr{j,k,l}.chi2stat;
        end
    end
end
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

% Find peaks in power spectra.
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        [pxxdet,fdet]=pwelch(Xdet{j}{k},[],[],[],Fs);
        [pxxsto,fsto]=pwelch(Xsto{j}{k},[],[],[],Fs);
        PSDdetpk_ampl(j,k)=pxxdet(pxxdet==max(pxxdet));
        PSDdetpk_freq(j,k)=fdet(pxxdet==max(pxxdet));
        PSDstopk_ampl(j,k)=pxxsto(pxxsto==max(pxxsto));
        PSDstopk_freq(j,k)=fsto(pxxsto==max(pxxsto));
    end
end
analysisstep = 6;
else
    disp('Run previous cell.');
end
%% (7) PLOT DATA
biftype2='SNIC';

for i = 1:6
    detinclind = find(detincl(i,:)==0); detinclind=detinclind(end);
    stoinclind = find(stoincl(i,:)==0); 
    if length(stoinclind) > 1
        stoinclind=stoinclind(end);
    else
        stoinclind=1;
    end
figure(1);subplot(3,2,i);patch([I(1) I(1) I(detinclind) I(detinclind)],[0 max(pkspikeratedet(i,:)) max(pkspikeratedet(i,:)) 0],[0.85 0.85 0.85],'EdgeColor','none');hold on;patch([I(1) I(1) I(stoinclind) I(stoinclind)],[0 max(pkspikeratedet(i,:)) max(pkspikeratedet(i,:)) 0],[0.6 0.3 0.3],'EdgeColor','none');title(sprintf('%s%s %s%s','Peaks, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Quiescent Region');
figure(2);subplot(3,2,i);patch([I(1) I(1) I(detinclind) I(detinclind)],[0 max(pkspikeratedet(i,:)) max(pkspikeratedet(i,:)) 0],[0.85 0.85 0.85],'EdgeColor','none');hold on;patch([I(1) I(1) I(stoinclind) I(stoinclind)],[0 max(pkspikeratedet(i,:)) max(pkspikeratedet(i,:)) 0],[0.6 0.3 0.3],'EdgeColor','none');title(sprintf('%s%s %s%s','Troughs, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Quiescent Region');
end;clear i

for i = 1:6
    detinclind = find(detincl(i,:)==0); detinclind=detinclind(end);
    stoinclind = find(stoincl(i,:)==0); 
    if length(stoinclind) > 1
        stoinclind=stoinclind(end);
    else
        stoinclind=1;
    end
figure(3);subplot(3,2,i);patch([I(1) I(1) I(detinclind) I(detinclind)],[0 1 1 0],[0.85 0.85 0.85],'EdgeColor','none');hold on;patch([I(1) I(1) I(stoinclind) I(stoinclind)],[0 1 1 0],[0.6 0.3 0.3],'EdgeColor','none');plot(I,pkspikeratedet(i,:),'k--');hold on;plot(I,pkspikeratesto(i,:),'r');title(sprintf('%s%s %s%s','Peaks, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Spike Rate');
figure(4);subplot(3,2,i);patch([I(1) I(1) I(detinclind) I(detinclind)],[0 1 1 0],[0.85 0.85 0.85],'EdgeColor','none');hold on;patch([I(1) I(1) I(stoinclind) I(stoinclind)],[0 1 1 0],[0.6 0.3 0.3],'EdgeColor','none');plot(I,trspikeratedet(i,:),'k--');hold on;plot(I,trspikeratesto(i,:),'r');title(sprintf('%s%s %s%s','Troughs, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Spike Rate');
end;clear i

for i = 1:6
    detinclind = find(detincl(i,:)==0); detinclind=detinclind(end);
    stoinclind = find(stoincl(i,:)==0); 
    if length(stoinclind) > 1
        stoinclind=stoinclind(end);
    else
        stoinclind=1;
    end
figure(5);subplot(3,2,i);patch([I(1) I(1) I(detinclind) I(detinclind)],[0 max(CDstopk(i,:)) max(CDstopk(i,:)) 0],[0.85 0.85 0.85],'EdgeColor','none');hold on;patch([I(1) I(1) I(stoinclind) I(stoinclind)],[0 max(CDstopk(i,:)) max(CDstopk(i,:)) 0],[0.6 0.3 0.3],'EdgeColor','none');plot(I,CDdetpk(i,:),'k--');plot(I,CDstopk(i,:),'r');plot(I,ones(1,length(I)),'g--');title(sprintf('%s%s %s%s','Peaks, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Coeff. Dispersion');
figure(6);subplot(3,2,i);patch([I(1) I(1) I(detinclind) I(detinclind)],[0 max(CDstotr(i,:)) max(CDstotr(i,:)) 0],[0.85 0.85 0.85],'EdgeColor','none');hold on;patch([I(1) I(1) I(stoinclind) I(stoinclind)],[0 max(CDstotr(i,:)) max(CDstotr(i,:)) 0],[0.6 0.3 0.3],'EdgeColor','none');plot(I,CDdettr(i,:),'k--');plot(I,CDstotr(i,:),'r');plot(I,ones(1,length(I)),'g--');title(sprintf('%s%s %s%s','Troughs, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Coeff. Dispersion');
end;clear i

for i = 1:6
    detinclind = find(detincl(i,:)==0); detinclind=detinclind(end);
    stoinclind = find(stoincl(i,:)==0); 
    if length(stoinclind) > 1
        stoinclind=stoinclind(end);
    else
        stoinclind=1;
    end
figure(7);subplot(3,2,i);patch([I(1) I(1) I(detinclind) I(detinclind)],[0 max(pkdiffusiondet(i,:)) max(pkdiffusiondet(i,:)) 0],[0.85 0.85 0.85],'EdgeColor','none');hold on;patch([I(1) I(1) I(stoinclind) I(stoinclind)],[0 max(pkdiffusiondet(i,:)) max(pkdiffusiondet(i,:)) 0],[0.6 0.3 0.3],'EdgeColor','none');plot(I,pkdiffusiondet(i,:),'k--');plot(I,pkdiffusionsto(i,:),'r');title(sprintf('%s%s %s%s','Peaks, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Diffusion Coeff.');
figure(8);subplot(3,2,i);patch([I(1) I(1) I(detinclind) I(detinclind)],[0 max(trdiffusiondet(i,:)) max(trdiffusiondet(i,:)) 0],[0.85 0.85 0.85],'EdgeColor','none');hold on;patch([I(1) I(1) I(stoinclind) I(stoinclind)],[0 max(trdiffusiondet(i,:)) max(trdiffusiondet(i,:)) 0],[0.6 0.3 0.3],'EdgeColor','none');plot(I,trdiffusiondet(i,:),'k--');plot(I,trdiffusionsto(i,:),'r');title(sprintf('%s%s %s%s','Troughs, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Diffusion Coeff.');
end;clear i

for i = 1:6
    detinclind = find(detincl(i,:)==0); detinclind=detinclind(end);
    stoinclind = find(stoincl(i,:)==0); 
    if length(stoinclind) > 1
        stoinclind=stoinclind(end);
    else
        stoinclind=1;
    end
figure(9);subplot(3,2,i);patch([I(1) I(1) I(detinclind) I(detinclind)],[0 max(pkIEIspikeratiodet(i,:)) max(pkIEIspikeratiodet(i,:)) 0],[0.85 0.85 0.85],'EdgeColor','none');hold on;patch([I(1) I(1) I(stoinclind) I(stoinclind)],[0 max(pkIEIspikeratiodet(i,:)) max(pkIEIspikeratiodet(i,:)) 0],[0.6 0.3 0.3],'EdgeColor','none');plot(I,pkIEIspikeratiodet(i,:),'k--');plot(I,pkIEIspikeratiosto(i,:),'r');title(sprintf('%s%s %s%s','Peaks, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('MeanIEI/(1/spikerate)');
figure(10);subplot(3,2,i);patch([I(1) I(1) I(detinclind) I(detinclind)],[0 max(trIEIspikeratiodet(i,:)) max(trIEIspikeratiodet(i,:)) 0],[0.85 0.85 0.85],'EdgeColor','none');hold on;patch([I(1) I(1) I(stoinclind) I(stoinclind)],[0 max(trIEIspikeratiodet(i,:)) max(trIEIspikeratiodet(i,:)) 0],[0.6 0.3 0.3],'EdgeColor','none');plot(I,trIEIspikeratiodet(i,:),'k--');plot(I,trIEIspikeratiosto(i,:),'r');title(sprintf('%s%s %s%s','Troughs, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('MeanIEI/(1/spikerate)');
end;clear i



%% (X) SAVE DATA

savefile='/Users/joshsalvi/Desktop/SNICstochoutput2-analyzed-4.mat';    % New file
%savefile=file;  % Overwrite old file
save(savefile,'-v7.3','I','detincl','stoincl','pkspikeratedet','pkspikeratesto','trspikeratedet','trspikeratedet','CDstopk','CDstotr','CDdetpk','CDdettr','pkdiffusiondet','pkdiffusionsto','trdiffusiondet','trdiffusionsto','pkIEIspikeratiodet','pkIEIspikeratiosto','trIEIspikeratiodet','trIEIspikeratiosto','noiselevel');
