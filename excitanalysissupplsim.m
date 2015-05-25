function excitanalysissupplsim(directory,biftype,threshrange,threshcheck,fftpkcheck,stiffind)
%
% Finds additional metrics on simulated data from bifurcation normal forms.
%
% excitanalysissupplsim(directory,biftype,threshrange,threshcheck,fftpkcheck)
%
% Directory: The directory of the files to be searched e.g. '/Users/.../'
% biftype: bifurcation type to be analyzed e.g. 1
% threshrange: range of thresholds e.g. linspace(0.01,1,100)
% threshcheck: check range of thresholds? (1=yes;0=no) default=1
% fftpkcheck: check FFT? (1=yes;0=no) default=1
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%

if exist('threshcheck')==0
    threshcheck=1;
end
if exist('fftpkcheck')==0
    fftpkcheck=1;
end
    
%Fs = 1/t(2);
dwnspl=1;offset=0;
if biftype ~= 1 && isempty(dir([directory '*-rect-pt' num2str(1) '.mat']))==1
    for j = 1:5
        files{j} = dir([directory '*-all-pt' num2str(j) '.mat']);
    end
elseif biftype == 4
    for j = 1:5
        files{j} = dir([directory '*-pt' num2str(j) '.mat']);
    end
elseif isempty(dir([directory '*-rect-pt' num2str(1) '.mat']))==0
    for j = 1:5
        files{j} = dir([directory '*-rect-pt' num2str(j) '.mat']);
    end
elseif isempty(findstr(directory,'offset')) == 1
    for j = 1:5
        files{j} = dir([directory '*-pt' num2str(j) '.mat']);
    end
else
    for j = 1:5
        files{j} = dir([directory '*-all-pt' num2str(j) '.mat']);
    end
end
%{
if biftype == 4
    stiffind = input('Which index in stiffness? (1,2,3,4,5)...   ');
end
%}

stiffind2=stiffind;

for m = 1:5
    load([directory files{m}.name]);
    stiffind=stiffind2;
    if m == 1
        % Pre-allocate memory
        dipsto = zeros(5,length(Xdet),length(Xdet{1}));dipdet = zeros(5,length(Xdet),length(Xdet{1}));
        punisto = zeros(5,length(Xdet),length(Xdet{1}));punidet = zeros(5,length(Xdet),length(Xdet{1}));
        fftmaxpkstoind = zeros(5,length(Xdet),length(Xdet{1}));fftmaxpkdetind = zeros(5,length(Xdet),length(Xdet{1}));
        fftmaxpkheightsto = zeros(5,length(Xdet),length(Xdet{1}));fftmaxpkheightdet = zeros(5,length(Xdet),length(Xdet{1}));
        fftmaxpkfreqtsto = zeros(5,length(Xdet),length(Xdet{1}));fftmaxpkfreqtdet = zeros(5,length(Xdet),length(Xdet{1}));
        fftmaxpkwidthtdetL = zeros(5,length(Xdet),length(Xdet{1}));fftmaxpkwidthtstoL = zeros(5,length(Xdet),length(Xdet{1}));
        fftmaxpkwidthtstoR = zeros(5,length(Xdet),length(Xdet{1}));fftmaxpkwidthtdetR = zeros(5,length(Xdet),length(Xdet{1}));
        Xlowsto = zeros(5,length(Xdet),length(Xdet{1}));Xlowdet = zeros(5,length(Xdet),length(Xdet{1}));
        Qdet = zeros(5,length(Xdet),length(Xdet{1}));Qsto = zeros(5,length(Xdet),length(Xdet{1}));
        Xupsto = zeros(5,length(Xdet),length(Xdet{1}));Xupdet = zeros(5,length(Xdet),length(Xdet{1}));
        CDdetpk = cell(5,length(threshrange));CDstopk = cell(5,length(threshrange));
    end
    
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
    if exist('Xsto')==1
        Xsto1=Xsto; Xdet1=Xdet;clear Xdet Xsto;
    end
    if isempty (findstr([directory files{m}.name],'force')) == 0
        mu = F; rt=mu; I=rt;
    else
        mu = ke; rt=mu; I=rt;
    end
    sizeX = size(Xdet1);
for j = 1:sizeX(1)       % Isolate the appropriate index
    for qp = 1:sizeX(3)
        Xsto{qp}{j} = Xsto1{j,stiffind,qp}(1,1:dwnspl:end) + offset;
        Xdet{qp}{j} = Xdet1{j,stiffind,qp}(1,1:dwnspl:end) + offset;
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
    if exist('mu') == 0
    if exist('I') == 1
        mu = I;
    elseif exist('rt') == 1
        mu = rt;
    else 
        disp('No control parameter exists');
        break;
    end
    end
    if m == 1
        for j = 1:5
            for k = 1:length(threshrange)
                CDdetpk{j,k} = zeros(length(Xdet),length(Xdet{1}));
                CDstopk{j,k} = zeros(length(Xdet),length(Xdet{1}));
            end
        end
    end
  if threshcheck==1  
    for n = 1:length(threshrange)
        clear IEIpkdet IEIpksto pkd pks
        IEIpkdet = cell(length(Xdet),length(Xdet{1}));IEIpksto = cell(length(Xdet),length(Xdet{1}));
        for j = 1:length(Xdet)
            for k = 1:length(Xdet{1})
                % Peak finding and CD
                pkd{j,k} = PTDetect(Xdet{j}{k},threshrange(n));
                pks{j,k} = PTDetect(Xsto{j}{k},threshrange(n));
                for l = 2:length(pkd{j,k})
                    IEIpkdet{j,k}(l-1) = (pkd{j,k}(l) - pkd{j,k}(l-1));
                end
                for l = 2:length(pks{j,k})
                    IEIpksto{j,k}(l-1) = (pks{j,k}(l) - pks{j,k}(l-1));
                end
                if isempty(IEIpkdet{j,k})==0
                    CDdetpk{m,n}(j,k) = var(IEIpkdet{j,k})/mean(IEIpkdet{j,k});
                else
                    CDdetpk{m,n}(j,k) = 0;
                end
                if isempty(IEIpksto{j,k})==0
                    CDstopk{m,n}(j,k) = var(IEIpksto{j,k})/mean(IEIpksto{j,k});
                else
                    CDstopk{m,n}(j,k) = 0;
                end
                % Phase distribution
                clear xdethilb xstohilb
                xdethilb = hilbert(Xdet{j}{k}(round(length(Xdet{j}{k})/2:end)));
                xstohilb = hilbert(Xsto{j}{k}(round(length(Xsto{j}{k})/2:end)));
                xdetphase = atan2(imag(xdethilb),real(xdethilb));
                xstophase = atan2(imag(xstohilb),real(xstohilb));
                clear i;
                VSdet(m,j,k) = abs(1/length(xdetphase)*sum(exp(i*xdetphase))); % vector strength (det)
                VSsto(m,j,k) = abs(1/length(xstophase)*sum(exp(i*xstophase))); % vector strength (sto)
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
    [Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
    winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
    for j = 1:length(Xsto)
        for k = 1:length(Xsto{1})
            if fftpkcheck==1
            clear fftpkspkheight fftpkdpkheight
            [Xstofft, fstofft, XstofftCI]= pwelch(Xsto{j}{k},winfunc,noverlap,NPSD,Fs,'ConfidenceLevel',0.999);
            [Xdetfft, fdetfft, XdetfftCI]= pwelch(Xdet{j}{k},winfunc,noverlap,NPSD,Fs,'ConfidenceLevel',0.999);
            [a1s b1s] = twoclass(Xstofft,1e-4);
            [a1d b1d] = twoclass(Xdetfft,1e-4);
            pkfs = PTDetect(Xstofft,max([a1s b1s]));
            pkfd = PTDetect(Xstofft,max([a1d b1d]));
            mxs = 1; mxd = 1;
            for q = 1:length(pkfs)  % test each peak (stochastic) -> is it statistically distinct from peaks with less then half of its power on both sides?
                fhps = find(Xstofft<=(Xstofft(pkfs(q))/2));
                testL = sign(XstofftCI(pkfs(q),1)-XstofftCI(fhps(fhps<pkfs(q)),2));
                testH = sign(XstofftCI(pkfs(q),1)-XstofftCI(fhps(fhps>pkfs(q)),2));
                if sum(testL==1) >= 1 && sum(testH==1) >= 1 
                    fftpkspkheight(q) = Xstofft(pkfs(q));
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
                fftmaxpkheightsto(m,j,k) = Xstofft(fftmaxpkstoind(m,j,k));
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
            
            for q = 1:length(pkfd)  % test each peak (deterministic) -> is it statistically distinct from peaks with less then half of its power on both sides?
                fhpd = find(Xdetfft<=(Xdetfft(pkfd(q))/2));
                testL = sign(Xdetfft(pkfd(q),1)-XdetfftCI(fhpd(fhpd<pkfd(q)),2));
                testH = sign(XdetfftCI(pkfd(q),1)-XdetfftCI(fhpd(fhpd>pkfd(q)),2));
                if sum(testL==1) >= 1 && sum(testH==1) >= 1
                    fftpkdpkheight(q) = Xdetfft(pkfd(q));
                    fftpkdpkfreq(q) = fdetfft(pkfd(q));
                    fhpdL = fhpd(fhpd<pkfd(q));fhpdH = fhpd(fhpd>pkfd(q));
                    Lind = fhpdL(testL==1);Lind=Lind(end);
                    Hind = fhpdH(testH==1);Hind=Hind(1);
                    fftpkdwidth(m,q,1:2) = [fdetfft(Lind) fdetfft(Hind)];
                    mxd = mxd+1;
                end
            end
            if mxd > 1
                clear ind1
                ind1 = find(fftpkdpkheight(:)==max(fftpkdpkheight(:))); ind1=ind1(1);
                if numel(ind1)==0;ind1=1;
                end
                fftmaxpkdetind(m,j,k) = pkfd(ind1);
                fftmaxpkheightdet(m,j,k) = Xdetfft(fftmaxpkdetind(m,j,k));
                fftmaxpkfreqtdet(m,j,k) =  fdetfft(fftmaxpkdetind(m,j,k));
                fftmaxpkwidthtdetL(m,j,k) = fftpkdwidth(m,ind1,1);
                fftmaxpkwidthtdetR(m,j,k) = fftpkdwidth(m,ind1,2);
                Qdet(m,j,k) = fftmaxpkfreqtdet(m,j,k)/(fftmaxpkwidthtdetR(m,j,k)-fftmaxpkwidthtdetL(m,j,k));
            else
                fftmaxpkheightdet(m,j,k) = 0;
                fftmaxpkfreqtdet(m,j,k) = 0;
                fftmaxpkwidthtdetL(m,j,k) = 0;
                fftmaxpkwidthtdetR(m,j,k) = 0;
                Qdet(m,j,k) = 0;
            end  
            end
            if k == 1
                disp('Dip test...');
                disp([num2str(j) '/' num2str(length(Xsto))])
            end
            % Find the dip statistic and associated p-value
            Nboot = 10^1;   % number of bootstraps
            %[as bs] = hist(Xsto{j}{k},freedmandiaconis(Xsto{j}{k}));as=as./sum(as);
            %[ad bd] = hist(Xsto{j}{k},freedmandiaconis(Xdet{j}{k}));ad=ad./sum(ad);
            [dipsto(m,j,k), punisto(m,j,k), Xlowsto(m,j,k), Xupsto(m,j,k)]=HartigansDipSignifTest(Xsto{j}{k},Nboot);
            [dipdet(m,j,k), punidet(m,j,k), Xlowdet(m,j,k), Xupdet(m,j,k)]=HartigansDipSignifTest(Xdet{j}{k},Nboot);
            if threshcheck ~= 1
                if mod(k,25)==0
                    disp([num2str(k) '/' num2str(length(Xsto{1}))])
                end
            end
           
        end
    end
    if threshcheck==1 && fftpkcheck==1
disp(['m = ' num2str(m) '/5']);
disp(['Saving' num2str(m) '...']);
if biftype == 4
    save([directory 'stiffind-' num2str(stiffind) '-output' num2str(m) '.mat'],'Qsto','Qdet','VSdet','VSsto','CDdetpk','CDstopk','punisto','punidet','fftmaxpkheightdet','fftmaxpkheightsto','fftmaxpkfreqtdet','fftmaxpkfreqtsto','fstofft','mu','threshrange','dipsto','dipdet');
    disp('Finished.');
else
    save([directory 'output' num2str(m) '.mat'],'Qsto','Qdet','VSdet','VSsto','CDdetpk','CDstopk','punisto','punidet','fftmaxpkheightdet','fftmaxpkheightsto','fftmaxpkfreqtdet','fftmaxpkfreqtsto','fstofft','mu','threshrange','dipsto','dipdet');
    disp('Finished.');
end
    end
end

%Pre-allocate memory
mu1d = zeros(length(threshrange),length(Xdet));mu1s = zeros(length(threshrange),length(Xdet));
CD1d = zeros(length(threshrange),length(Xdet));CD1s = zeros(length(threshrange),length(Xdet));
sizeC1 = size(CDdetpk);
sizeC2 = size(CDdetpk{1,1});

if threshcheck ==1 && fftpkcheck ==1
for n = 1:sizeC1(2)
    for j = 1:sizeC2(1)
        for k = 1:sizeC2(2)
            CDdetmean{n}(j,k) = mean([CDdetpk{1,n}(j,k) CDdetpk{2,n}(j,k) CDdetpk{3,n}(j,k) CDdetpk{4,n}(j,k) CDdetpk{5,n}(j,k)]);
            CDdetsem{n}(j,k) = std([CDdetpk{1,n}(j,k) CDdetpk{2,n}(j,k) CDdetpk{3,n}(j,k) CDdetpk{4,n}(j,k) CDdetpk{5,n}(j,k)])/sqrt(5);
            CDstomean{n}(j,k) = mean([CDstopk{1,n}(j,k) CDstopk{2,n}(j,k) CDstopk{3,n}(j,k) CDstopk{4,n}(j,k) CDstopk{5,n}(j,k)]);
            CDstosem{n}(j,k) = std([CDstopk{1,n}(j,k) CDstopk{2,n}(j,k) CDstopk{3,n}(j,k) CDstopk{4,n}(j,k) CDstopk{5,n}(j,k)])/sqrt(5);
            Qdetmean(j,k) = mean([Qdet(1,j,k) Qdet(2,j,k) Qdet(3,j,k) Qdet(4,j,k) Qdet(5,j,k)]);
            Qdetsem(j,k) = std([Qdet(1,j,k) Qdet(2,j,k) Qdet(3,j,k) Qdet(4,j,k) Qdet(5,j,k)])/sqrt(5);
            Qstomean(j,k) = mean([Qsto(1,j,k) Qsto(2,j,k) Qsto(3,j,k) Qsto(4,j,k) Qsto(5,j,k)]);
            Qstosem(j,k) = std([Qsto(1,j,k) Qsto(2,j,k) Qsto(3,j,k) Qsto(4,j,k) Qsto(5,j,k)])/sqrt(5);
            VSdetmean(j,k) = mean([VSdet(1,j,k) VSdet(2,j,k) VSdet(3,j,k) VSdet(4,j,k) VSdet(5,j,k)]);
            VSstomean(j,k) = mean([VSsto(1,j,k) VSsto(2,j,k) VSsto(3,j,k) VSsto(4,j,k) VSsto(5,j,k)]);
            VSdetsem(j,k) = std([VSdet(1,j,k) VSdet(2,j,k) VSdet(3,j,k) VSdet(4,j,k) VSdet(5,j,k)])/sqrt(5);
            VSstosem(j,k) = std([VSsto(1,j,k) VSsto(2,j,k) VSsto(3,j,k) VSsto(4,j,k) VSsto(5,j,k)])/sqrt(5);
        end
        CD1d1 = findnearest(CDdetmean{n}(j,:),1);if isempty(CD1d1)==0;CD1d(n,j)=CDdetmean{n}(j,CD1d1(1));else CD1d(n,j)=NaN;end
        mu1d(n,j) = mu(CD1d1(1));       % control parameter at which CD ~ 1 (deterministic)
        CD1s1 = findnearest(CDstomean{n}(j,:),1);if isempty(CD1s1)==0;CD1s(n,j)=CDstomean{n}(j,CD1s1(1));else CD1s(n,j)=NaN;end
        mu1s(n,j) = mu(CD1s1(1));       % control parameter at which CD ~ 1 (stochastic)
    end
end
end
if threshcheck==1 && fftpkcheck ==1
if biftype == 4
    disp('Saving...');
    save([directory 'outputALL-stiffind-' num2str(stiffind) '.mat'],'Qsto','Qdet','VSdet','VSsto','mu1s','mu1d','CDdetmean','CDdetsem','CDstomean','CDstosem','Qstomean','Qstosem','Qdetmean','Qdetsem','VSdetmean','VSdetsem','VSstomean','VSstosem','dipsto','dipdet','punisto','punidet','fftmaxpkheightdet','fftmaxpkheightsto','fftmaxpkfreqtdet','fftmaxpkfreqtsto','fstofft','mu','threshrange','CDdetpk','CDstopk','CD1d','CD1s');
    disp('Finished.');
else
    disp('Saving...');
    save([directory 'outputALL.mat'],'Qsto','Qdet','VSdet','VSsto','mu1s','mu1d','CDdetmean','CDdetsem','CDstomean','CDstosem','Qstomean','Qstosem','Qdetmean','Qdetsem','VSdetmean','VSdetsem','VSstomean','VSstosem','dipsto','dipdet','punisto','punidet','fftmaxpkheightdet','fftmaxpkheightsto','fftmaxpkfreqtdet','fftmaxpkfreqtsto','fstofft','mu','threshrange','CDdetpk','CDstopk','CD1d','CD1s');
    disp('Finished.');
end
elseif threshcheck==1 && fftpkcheck ==0
    if biftype == 4
    disp('Saving...');
    save([directory 'output-NOTHRESHCHECK-stiffind-' num2str(stiffind) '.mat'],'VSdet','VSsto','mu1s','mu1d','VSdetmean','VSdetsem','VSstomean','VSstosem','dipsto','dipdet','punisto','punidet','mu','threshrange');
    disp('Finished.');
else
    disp('Saving...');
    save([directory 'output-NOTHRESHCHECK.mat'],'VSdet','VSsto','mu1s','mu1d','VSdetmean','VSdetsem','VSstomean','VSstosem','dipsto','dipdet','punisto','punidet','mu','threshrange');
    disp('Finished.');
    end
else
    if biftype == 4
    disp('Saving...');
    save([directory 'output-NOTHRESHCHECK-stiffind-' num2str(stiffind) '.mat'],'dipsto','dipdet','punisto','punidet','mu','threshrange');
    disp('Finished.');
else
    disp('Saving...');
    save([directory 'output-NOTHRESHCHECK.mat'],'dipsto','dipdet','punisto','punidet','mu','threshrange');
    disp('Finished.');
    end
end
