j=2;
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs{j}, 'Type');  %type of low-level graphics object
xdata = get(dataObjs{j}, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs{j}, 'YData');
udata = get(dataObjs{j}, 'uData');

%%

j=1;
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
zdata = get(dataObjs, 'ZData');

%% Extract error bars

j=2;
h = gca; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
udata = get(dataObjs, 'UData');

figure;
setfiguredefaults(length(ydata))
m = 1:length(ydata);m=fliplr(m);
for j = 1:length(ydata)
    a=find(ydata{m(j)}==max(ydata{m(j)}));
    ha = errorbar(xdata{m(j)}(a:end),ydata{m(j)}(a:end),udata{m(j)}(a:end));hold all;
    hb = get(gca,'children');        
    maxy(j) = max(ydata{m(j)});
end
for j = 1:length(ydata)
    errorbar_tick(hb(j),1000);
end
if max(maxy) >= 1.5
    axis([-1 1 0.2 1.1*max(maxy)]);
else
    axis([-1 1 0.2 1.8]);
end
hold on;plot(-5:1:5,ones(1,length(-5:1:5)),'g--')
legend([{'Thresh 1'};{'Thresh 2'};{'Thresh 3'};{'CV = 1'}])
xlabel('Control parameter');ylabel('Coefficient of variation');

%%
m=1; % plot number (index/threshold)

sizeX = size(IEIpkdet);

% Calculate SEM for spike rate
for j = 1:sizeX(1)
    for k = 1:a
        meanIEIpkdet(j,k)=mean(IEIpkdet{j,k});
        pkspikerateIEIpkdet{j,k}=1./(IEIpkdet{j,k});
        pkspikeratedetsem{m}(j,k)=std(pkspikerateIEIpkdet{j,k})/sqrt(length(pkspikerateIEIpkdet{j,k}));
    end
end

% Calculate errors for CD
for j = 1:sizeX(1)
    for k = 1:a
        N = length(IEIpkdet{j,k});
        SES = sqrt(6*N*(N-1)/((N-2)*(N+1)*(N+3)));
        SEK = 2*SES*sqrt((N^2-1)/((N-3)*(N+5)));
        varIEIpkdet{m}(j,k) = sqrt(var(IEIpkdet{j,k})^2*(2/(N-1)+(kurtosis(IEIpkdet{j,k})-3)/(N*SEK)));
        semIEIpkdet{m}(j,k) = std(IEIpkdet{j,k})/sqrt(length(IEIpkdet{j,k}));
        CDdetpkvar{m}(j,k) = sqrt((sqrt(var(IEIpkdet{j,k})^2*(2/(N-1)+(kurtosis(IEIpkdet{j,k})-3)/N))/CDdetpk(j,k))^2+(std(IEIpkdet{j,k})/sqrt(length(IEIpkdet{j,k})/mean(IEIpkdet{j,k})))^2);
    end
end


% Calculate dip statistic
dipyn=0;
%dipyn = input('Calculate dip p-value? (1=yes):  ');
if dipyn == 1
    Nboot = input('Number of bootstraps: ');
end

for k = 1:a
    for j = 1:(logdata.data(1,8))
   
        % Symmetry test [Asai '06]
        detmean(j,k) = mean(Xd_dwnspl{j}{k});
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
        if dipyn ~=1
            Nboot = 10^1; % number of bootstraps to find puni
        end
        [dipdet(j,k), punidet(j,k), Xlowdet(j,k), Xupdet(j,k)]=HartigansDipSignifTest(Xd_dwnspl{j}{k},Nboot);
        %disp(['j = ' num2str(j) ' k = ' num2str(k)]);
    end
end

% Calculate RMS magnitude
for k = 1:a
    for j = 1:(logdata.data(1,8))
        Xdcenter=Xd_dwnspl{j}{k}-mean(Xd_dwnspl{j}{k});
        RMSmag{m}(j,k) = std(Xdcenter);
        N=length(Xdcenter);
        RMSmagstd{m}(j,k) = sqrt(var(Xdcenter)^2*(2/(N-1)+(kurtosis(Xdcenter)-3)/N));
    end
end

PKampldetSEM{m}=PKamplsemdet;

% Calculate mutual information (real vs imaginary)
for k = 1:length(Xd_dwnspl)
        for l = 1:length(Xd_dwnspl{1})
            MIrealimag(k,l) = mutualinfo(real(Xd_dwnspl{k}{l}),imag(hilbert(Xd_dwnspl{k}{l})));
        end
end

%% Plot CD data

figure;
for j = 1:length(xdata)
        errorbar(xdata{j},ydata{j},reshape(CDdetpkvar{j},1,length(xdata{1}))-reshape(CDdetpkvar{j},1,length(xdata{1})),reshape(CDdetpkvar{j},1,length(xdata{1}))); hold all;
end

%% Plot spike rate data

figure;
for j = 1:length(xdata)
    if j ==length(xdata)
        errorbar(xdata{j},ydata{j},reshape(pkspikeratedetsem{j},1,length(xdata{1}))-reshape(pkspikeratedetsem{j},1,length(xdata{1})),reshape(pkspikeratedetsem{j},1,length(xdata{1}))); hold all;
    elseif j==1
        errorbar(xdata{j},ydata{j},reshape(pkspikeratedetsem{j},1,length(xdata{1})),reshape(pkspikeratedetsem{j},1,length(xdata{1}))-reshape(pkspikeratedetsem{j},1,length(xdata{1}))); hold all;
    else
        plot(xdata{j},ydata{j});hold all;
    end
end

%% Plot amplitude data

figure;
for j = 1:length(xdata)
    if j == 5
        errorbar(xdata{j},ydata{j},reshape(PKampldetSEM{j},1,length(xdata{1}))-reshape(PKampldetSEM{j},1,length(xdata{1})),reshape(PKampldetSEM{j},1,length(xdata{1}))); hold all;
    elseif j==10
        errorbar(xdata{j},ydata{j},reshape(PKampldetSEM{j},1,length(xdata{1})),reshape(PKampldetSEM{j},1,length(xdata{1}))-reshape(PKampldetSEM{j},1,length(xdata{1}))); hold all;
    else
         errorbar(xdata{j},ydata{j},reshape(PKampldetSEM{j},1,length(xdata{1})));hold all;
       % plot(xdata{j},ydata{j});hold all;
    end
end

%% Plot dip statistic and RMS magnitude

figure;
subplot(2,2,1);plot(xdata{1},reshape(dipdet,1,length(xdata{1})));xlabel('Control parameter');ylabel('Dip statistic');set(gca,'xdir','reverse')
subplot(2,2,2);errorbar(xdata{1},reshape(RMSmag{1},1,length(xdata{1})),reshape(RMSmagstd{1},1,length(xdata{1})));xlabel('Control parameter');ylabel('RMS magnitude');set(gca,'xdir','reverse')
subplot(2,2,3);plot(xdata{1},reshape(Kstatdet,1,length(xdata{1})));xlabel('Control parameter');ylabel('Excess kurtosis');set(gca,'xdir','reverse')
subplot(2,2,4);plot(xdata{1},reshape(KSstatdet,1,length(xdata{1})));xlabel('Control parameter');ylabel('KS-test statistic');set(gca,'xdir','reverse')

%% Plot mutual information
sizeM=size(MIrealimag);
MIrealimag=reshape(MIrealimag,1,sizeM(1)*sizeM(2));
figure;plot(xdata{1},MIrealimag(1:length(xdata{1})));
xlabel('Control parameter');ylabel('Mutual information (real/imag)');
set(gca,'xdir','reverse')


%% Bootstrap data
sizeI = size(IEIpkdet);
Nboot = 1e4;
m=2;
clear CDpkdetSEM pkspikeratedetSEM
for j = 1:sizeI(1)
    for k = 1:sizeI(2)
        x = IEIpkdet{j,k};
        if isempty(x) == 0 && length(x) > 1
            [a]=bootstrp(Nboot,@(x) var(x)/mean(x),x);
            [a2]=bootstrp(Nboot,'mean',1./x);
            CDpkdetSEM(j,k) = std(a);
            pkspikeratedetSEM(j,k) = std(a2);
        end
    end
    disp(num2str(j));
end
%CDpkdetSEM = reshape(CDpkdetSEM,1,size(CDpkdetSEM,1)*size(CDpkdetSEM,2));   
%pkspikeratedetSEM = reshape(pkspikeratedetSEM,1,size(pkspikeratedetSEM,1)*size(pkspikeratedetSEM,2)); 
%L2 = length(CDpkdetSEM);L3 = length(pkspikeratedetSEM);scaling=ydata{m}(1)./(var(IEIpkdet{1,1})/mean(IEIpkdet{1,1}));
%CDpkdetSEM(L2+1:length(ydata{1}))=0;
%pkspikeratedetSEM(L3+1:length(ydata{1}))=0;

%% Bootstrap the dip statistic
clc;
sizeI = size(IEIpkdet);
Nboot = 1e4;
m=2;
clear CDpkdetSEM pkspikeratedetSEM
disp('Running...');
tic;
for j = 1:sizeI(1)
    if j ==1
    disp(['Completed:  '  num2str(0) '/' num2str(sizeI(1)) sprintf('\nTime elapsed:  ')  num2str(toc)]);  
    end
    for k = 1:sizeI(2)
        x = Xd_dwnspl{j}{k};
        if isempty(x) == 0 && length(x) > 1
            [dipdet(j,k) dipdetpval(j,k)] = HartigansDipSignifTest(x,1);
            [a]=bootstrp(Nboot,@(x) HartigansDipTest(x),x);
            %[p1]=bootstrp(Nboot,@(x) HartigansDipSignifTestpvalue(x,10),x);
            dipdetSEM(j,k) = std(a);
            %dipdetpvalSEM(j,k) = std(p1);
        end
    end
    disp(['Completed:  '  num2str(j) '/' num2str(sizeI(1)) sprintf('\nTime elapsed:  ')  num2str(toc)]);
end

%% Plot the result - CD with errors
m=2;
figure;errorbar(xdata{m},ydata{m},CDpkdetSEM-CDpkdetSEM,CDpkdetSEM.*scaling/2)

%% Plot the result - spike rate with errors
m=2;mm=2;
if mm==2
figure;errorbar(xdata{m},ydata{m},pkspikeratedetSEM,pkspikeratedetSEM-pkspikeratedetSEM)
elseif mm==1
    figure;errorbar(xdata{m},ydata{m},pkspikeratedetSEM-pkspikeratedetSEM,pkspikeratedetSEM)
end

%% bootstrap the coefficient of variation
files{1} = '/Volumes/Promise Pegasus/Manual Backup/Lab/Clamp Data/2013-04-19.01/Ear 1/Cell 8/Extracted Data2.matDwnspl20-Thresh6-output.mat';
files{2} = '/Volumes/Promise Pegasus/Manual Backup/Lab/Clamp Data/2013-04-19.01/Ear 1/Cell 8/Extracted Data2.matDwnspl20-Thresh8-output.mat';
files{3} = '/Volumes/Promise Pegasus/Manual Backup/Lab/Clamp Data/2013-04-19.01/Ear 1/Cell 8/Extracted Data2.matDwnspl20-Thresh10-output.mat';

Nboot = 1e3;

for qr = 1:length(files)
    disp('loading')
    load(files{qr});
for j = 1:size(pktimedet,1)
    for k = 1:size(pktimedet,2)
        disp(num2str(k));
        xpk = pktimedet{j,k};
        xtr = trtimedet{j,k};
        if length(xpk) >1;
        cvpktimedet{qr}(j,k) = std(xpk)/mean(xpk);
        else
             cvpktimedet{qr}(j,k) = NaN;
        end
        if length(xtr) >1
        cvtrtimedet{qr}(j,k) = std(xtr)/mean(xtr);
        else
            cvtrtimedet{qr}(j,k) = NaN;
        end
        mpk{qr}(j,k) = isnan(cvpktimedet{qr}(j,k));
        mtr{qr}(j,k) = isnan(cvtrtimedet{qr}(j,k));
        if isnan(cvpktimedet{qr}(j,k)) ==0
        [apk] = bootstrp(Nboot,@(x) std(x)/mean(x),xpk);
        
        cvpktimedetSEM{qr}(j,k) = std(apk(isnan(apk)==0));
        
        else
        cvpktimedetSEM{qr}(j,k) = NaN;
        
        end
        if isnan(cvtrtimedet{qr}(j,k)) ==0 
            [atr] = bootstrp(Nboot,@(x) std(x)/mean(x),xtr);
            cvtrtimedetSEM{qr}(j,k) = std(atr(isnan(atr)==0));
        else
            cvtrtimedetSEM{qr}(j,k) = NaN;
        end
    end
end
end

