function [pkspikeratedet, pkspikeratesto, CDdetpk, CDstopk,c12det] = excitabilityanalysissim(filename,dwnspl,biftype,c12constrain,offset,saveyn)
% This function imports simulation data and performs a peak-finding algorithm
% analysis. 
%
% [pkspikeratedet, pkspikeratesto, CDdetpk, CDstopk,c12det] = excitabilityanalysissim(filename,dwnspl,biftype,c12constrain,offset,saveyn)
% 
% Ex. excitabilityanalysissim('/Users/joshsalvi/Desktop/Hopf/Hopfstochoutput4-finemu-analyzed-dwnspl10.mat',10,1,1)
%
% filename: name of MAT file with directory
% dwnspl: downsampling rate
% biftype: 1=supercritical Hopf, 2=SNIC, 3=subcritical Hopf, 4=hb model
% c12constrain: index of threshold from twoclass
%
%
% jsalvi@rockefeller.edu

% Import data
load(filename);

% Import data
%file='/Users/joshsalvi/Desktop/SNICstochoutput.mat';
%load(file);

% Type of bifurcation
        % 1=supercritical Hopf; 2=SNIC; 3=subcritical Hopf; 4=hair-bundle
        % model
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

% Define time vector
tvec=t(1:dwnspl:end);
tvec=tvec(1:length(Xdet{1}{1}));
% Define the sample rate if not already done
Fs=1/(tvec(2)-tvec(1));

% Time range
tmin=tvec(0.1*length(tvec));
tmax=tvec(end);
minindex = find(abs(tvec-tmin)==min(abs(tvec-tmin)));
maxindex = find(abs(tvec-tmax)==min(abs(tvec-tmax)))-1;
Ttotal = tmax-tmin;

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
        [stohista{j,k}, stohistb{j,k}] = hist(Xsto{j}{k},freedmandiaconis(Xsto{j}{k}));
        bincount = sum(stohista{j,k}); stohista{j,k} = stohista{j,k}./bincount;
        detmean(j,k) = mean(Xdet{j}{k}); detvar(j,k) = var(Xdet{j}{k});
        stomean(j,k) = mean(Xsto{j}{k}); stovar(j,k) = var(Xsto{j}{k});
    end
end


% Nearest-neighbor clustering and peak finding

% Preallocate memory
pkdet = cell(length(Xdet),length(Xdet{1}));trdet = cell(length(Xdet),length(Xdet{1}));
pksto = cell(length(Xdet),length(Xdet{1}));trsto = cell(length(Xdet),length(Xdet{1}));
IEIpkdet = cell(length(Xdet),length(Xdet{1}));IEItrdet = cell(length(Xdet),length(Xdet{1}));
IEIpksto = cell(length(Xdet),length(Xdet{1}));IEItrsto = cell(length(Xdet),length(Xdet{1}));

NNerr = 10^-6;   % error threshold for nearest-neighbor clustering
for j = 1:length(Xdet)
    [c1det,c2det]=twoclass(Xdet{1}{c12constrain},NNerr);  % nearest-neighbor clustering
    for k = 1:length(Xdet{1})
            [pkdet{j,k},trdet{j,k}] = PTDetect(Xdet{j}{k}, max([c1det c2det]));
            for l = 2:length(pkdet{j,k})
                IEIpkdet{j,k}(l-1) = (pkdet{j,k}(l) - pkdet{j,k}(l-1))/Fs;
            end
            for l = 2:length(trdet{j,k})
                IEItrdet{j,k}(l-1) = (trdet{j,k}(l) - trdet{j,k}(l-1))/Fs;
            end
            
            [pksto{j,k},trsto{j,k}] = PTDetect(Xsto{j}{k}, max([c1det c2det]));
            for l = 2:length(pksto{j,k})
                IEIpksto{j,k}(l-1) = (pksto{j,k}(l) - pksto{j,k}(l-1))/Fs;
            end
            for l = 2:length(trsto{j,k})
                IEItrsto{j,k}(l-1) = (trsto{j,k}(l) - trsto{j,k}(l-1))/Fs;
            end
    end
end
c12det=max([c1det c2det]);

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
IEIvarpkdet = zeros(length(Xdet),length(Xdet{1}));IEIvartrdet = zeros(length(Xdet),length(Xdet{1}));
IEIvarpksto = zeros(length(Xdet),length(Xdet{1}));IEIvartrsto = zeros(length(Xdet),length(Xdet{1}));
IEImeanpkdet = zeros(length(Xdet),length(Xdet{1}));IEImeantrdet = zeros(length(Xdet),length(Xdet{1}));
IEImeanpksto = zeros(length(Xdet),length(Xdet{1}));IEImeantrsto = zeros(length(Xdet),length(Xdet{1}));

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
        IEIvarpkdet(j,k) = var(IEIpkdet{j,k});IEIvartrdet(j,k) = var(IEItrdet{j,k});
        IEIvarpksto(j,k) = var(IEIpksto{j,k});IEIvartrsto(j,k) = var(IEItrsto{j,k});
        IEImeanpkdet(j,k) = mean(IEIpkdet{j,k});IEImeantrdet(j,k) = mean(IEItrdet{j,k});
        IEImeanpksto(j,k) = mean(IEIpksto{j,k});IEImeantrsto(j,k) = mean(IEItrsto{j,k});
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

else
    

% Define time vector
tvec=t(1:dwnspl:end);
tvec=tvec(1:length(Xdet{1}{1}{1}));
% Define the sample rate if not already done
Fs=1/(tvec(2)-tvec(1));

% Time range
tmin=tvec(0.1*length(tvec));
tmax=tvec(end);
minindex = find(abs(tvec-tmin)==min(abs(tvec-tmin)));
maxindex = find(abs(tvec-tmax)==min(abs(tvec-tmax)))-1;
Ttotal = tmax-tmin;

% Pre-allocate memory
detmean = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
stomean = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
dethista = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1})); dethistb = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
stohista = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1})); stohistb = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
detvar = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1})); stovar = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));


% Generate histograms and calculate statistics
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        for l = 1:length(Xdet{1}{1})
            Xdet{j}{k}{l} = bsxfun(@minus, Xdet{j}{k}{l}, mean(Xdet{j}{k}{l}(minindex:maxindex)));
            Xsto{j}{k}{l} = bsxfun(@minus, Xsto{j}{k}{l}, mean(Xsto{j}{k}{l}(minindex:maxindex)));
            [dethista{j,k,l}, dethistb{j,k,l}] = hist(Xdet{j}{k}{l},freedmandiaconis(Xdet{j}{k}{l}));
            bincount = sum(dethista{j,k,l}); dethista{j,k,l} = dethista{j,k,l}./bincount;
            [stohista{j,k,l}, stohistb{j,k,l}] = hist(Xsto{j}{k}{l},freedmandiaconis(Xsto{j}{k}{l}));
            bincount = sum(stohista{j,k,l}); stohista{j,k,l} = stohista{j,k,l}./bincount;
            detmean(j,k,l) = mean(Xdet{j}{k}{l}); detvar(j,k,l) = var(Xdet{j}{k}{l});
            stomean(j,k,l) = mean(Xsto{j}{k}{l}); stovar(j,k,l) = var(Xsto{j}{k}{l});
    end
    end
end

% Nearest-neighbor clustering and peak finding

% Preallocate memory
pkdet = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));trdet = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
pksto = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));trsto = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
IEIpkdet = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));IEItrdet = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
IEIpksto = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));IEItrsto = cell(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));

NNerr = 10^-6;   % error threshold for nearest-neighbor clustering
for j = 1:length(Xdet)
    [c1det,c2det]=twoclass(Xdet{c12constrain}{1}{1},NNerr);  % nearest-neighbor clustering
    for k = 1:length(Xdet{1})
        for l = 1:length(Xdet{1}{1})
            [pkdet{j,k,l},trdet{j,k,l}] = PTDetect(Xdet{j}{k}{l}, max([c1det c2det]));
            for m = 2:length(pkdet{j,k,l})
                IEIpkdet{j,k,l}(m-1) = (pkdet{j,k,l}(m) - pkdet{j,k,l}(m-1))/Fs;
            end
            for m = 2:length(trdet{j,k,l})
                IEItrdet{j,k,l}(m-1) = (trdet{j,k,l}(m) - trdet{j,k,l}(m-1))/Fs;
            end
            
            [pksto{j,k,l},trsto{j,k,l}] = PTDetect(Xsto{j}{k}{l}, max([c1det c2det]));
            for m = 2:length(pksto{j,k,l})
                IEIpksto{j,k,l}(m-1) = (pksto{j,k,l}(m) - pksto{j,k,l}(m-1))/Fs;
            end
            for m = 2:length(trsto{j,k,l})
                IEItrsto{j,k,l}(m-1) = (trsto{j,k,l}(m) - trsto{j,k,l}(m-1))/Fs;
            end
    end
    end
end
c12det=max([c1det c2det]);
  
% Pre-allocate memory
pkspikeratedet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));pkspikeratesto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
trspikeratedet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));trspikeratesto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
CDdettime = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));CDstotime = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
CDdetpk = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));CDstopk = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
CDdettr = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));CDstotr = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
pkdiffusiondet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));pkdiffusionsto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
trdiffusiondet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));trdiffusionsto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
pktcorrdet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));pktcorrsto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
trtcorrdet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));trtcorrsto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
meanIEIpkdet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));meanIEItrdet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
meanIEIpksto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));meanIEItrsto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
pkIEIspikeratiodet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));pkIEIspikeratiosto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
trIEIspikeratiodet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));trIEIspikeratiosto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
IEIvarpkdet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));IEIvartrdet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
IEIvarpksto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));IEIvartrsto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
IEImeanpkdet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));IEImeantrdet = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));
IEImeanpksto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));IEImeantrsto = zeros(length(Xdet),length(Xdet{1}),length(Xdet{1}{1}));

% Calculate statistics: coefficient of dispersion, spike rate, diffusion
% coefficient, correlation time
for j = 1:length(Xdet)
    for k = 1:length(Xdet{1})
        for l = 1:length(Xdet{1}{1})
        pkspikeratedet(j,k,l) = length(pkdet{j,k,l})/Ttotal;      % spike rate, peaks
        pkspikeratesto(j,k,l) = length(pksto{j,k,l})/Ttotal;
        trspikeratedet(j,k,l) = length(trdet{j,k,l})/Ttotal;      % spike rate, troughs
        trspikeratesto(j,k,l) = length(trsto{j,k,l})/Ttotal;
        CDdettime(j,k,l) = (detvar(j,k,l))/detmean(j,k,l);        % coefficient of dispersion
        CDstotime(j,k,l) = (stovar(j,k,l))/stomean(j,k,l);
        IEIvarpkdet(j,k,l) = var(IEIpkdet{j,k,l});IEIvartrdet(j,k,l) = var(IEItrdet{j,k,l});
        IEIvarpksto(j,k,l) = var(IEIpksto{j,k,l});IEIvartrsto(j,k,l) = var(IEItrsto{j,k,l});
        IEImeanpkdet(j,k,l) = mean(IEIpkdet{j,k,l});IEImeantrdet(j,k,l) = mean(IEItrdet{j,k,l});
        IEImeanpksto(j,k,l) = mean(IEIpksto{j,k,l});IEImeantrsto(j,k,l) = mean(IEItrsto{j,k,l});
        CDdetpk(j,k,l) = var(IEIpkdet{j,k,l})/mean(IEIpkdet{j,k,l});
        CDdettr(j,k,l) = var(IEItrdet{j,k,l})/mean(IEItrdet{j,k,l});
        CDstopk(j,k,l) = var(IEIpksto{j,k,l})/mean(IEIpksto{j,k,l});
        CDstotr(j,k,l) = var(IEItrsto{j,k,l})/mean(IEItrsto{j,k,l});
        pkdiffusiondet(j,k,l) = 0.5*CDdetpk(j,k,l)^2*pkspikeratedet(j,k,l); % diffusion coefficient, peaks
        pkdiffusionsto(j,k,l) = 0.5*CDstopk(j,k,l)^2*pkspikeratesto(j,k,l);
        trdiffusiondet(j,k,l) = 0.5*CDdettr(j,k,l)^2*trspikeratedet(j,k,l); % diffusion coefficient, troughs
        trdiffusionsto(j,k,l) = 0.5*CDstotr(j,k,l)^2*trspikeratesto(j,k,l);
        pktcorrdet(j,k,l) = 0.5*pkdiffusiondet(j,k,l) - 1/pkspikeratedet(j,k,l);       % correlation time, peaks
        pktcorrsto(j,k,l) = 0.5*pkdiffusionsto(j,k,l) - 1/pkspikeratesto(j,k,l);
        trtcorrdet(j,k,l) = 0.5*trdiffusiondet(j,k,l) - 1/trspikeratedet(j,k,l);       % correlation time, troughs
        trtcorrsto(j,k,l) = 0.5*trdiffusionsto(j,k,l) - 1/trspikeratesto(j,k,l);
        meanIEIpkdet(j,k,l) = mean(IEIpkdet{j,k,l});meanIEItrdet(j,k,l) = mean(IEItrdet{j,k,l});    % mean IEI, should be equal to 1/(spikerate)
        meanIEIpksto(j,k,l) = mean(IEIpksto{j,k,l});meanIEItrsto(j,k,l) = mean(IEItrsto{j,k,l});    
        pkIEIspikeratiodet(j,k,l) = meanIEIpkdet(j,k,l)*pkspikeratedet(j,k,l);    % ratio of IEI to 1/spikerate
        trIEIspikeratiodet(j,k,l) = meanIEItrdet(j,k,l)*trspikeratedet(j,k,l);
        pkIEIspikeratiosto(j,k,l) = meanIEIpksto(j,k,l)*pkspikeratesto(j,k,l);    % ratio of IEI to 1/spikerate
        trIEIspikeratiosto(j,k,l) = meanIEItrsto(j,k,l)*trspikeratesto(j,k,l);
    end
end
end

end

if saveyn==1
    disp('Saving...');
filename2 = filename(1:end-4);
save(sprintf('%s%s%s%s%s%s%s%s',filename2,'-dwnspl',num2str(dwnspl),'-thresh',num2str(c12constrain),'-biftype',num2str(biftype),'.mat'))
disp('Finished.');
else
    disp('Not saved.');
end
end
