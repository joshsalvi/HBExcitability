biftype =2;
forcestiff=1;
hbcalc=1;

% Plot excitability data (and perform statistics)
% plotexcitabilitydata2(biftypemforcestiff)
% biftype: 1=supercritical Hopf, 2=SNIC, 3=subcritical Hopf, 4=Cusp/Fold ,5=HB model
% forcestiff: index for HB model
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%
if exist('I') == 0 && exist('mu') == 1
    I=mu;
end

if biftype ==2
% SNIC
close all;
maxiter=4;
noiselevel = [0.05 0.1 0.2 0.4];

for j = 1:maxiter
    for k = 1:500
        pkspikeratestoAVG(j,k)= mean([pkspikeratesto1(j,k) pkspikeratesto2(j,k) pkspikeratesto3(j,k) pkspikeratesto4(j,k) pkspikeratesto5(j,k)]);
        pkspikeratestoSEM(j,k)= std([pkspikeratesto1(j,k) pkspikeratesto2(j,k) pkspikeratesto3(j,k) pkspikeratesto4(j,k) pkspikeratesto5(j,k)])/sqrt(5);
    end
end


% coefficient of variation 1 = variance/mean -- between peak/peak and
% trough/trough
for j = 1:maxiter
    for k = 1:500
        CDstopkAVG(j,k)= mean([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)]);
        CDstopkSEM(j,k)= std([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)])/sqrt(5);
    end
end

% averages of coefficients for peak/trough, amplitudes and frequencies from
% FFT and amplitudes from peak/trough
for j = 1:maxiter
    for k = 1:500
        [~,CDpsto(j,k),~,CDstatsto(j,k)] =  ttest([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)],1);CDtstat(j,k)=CDstatsto(j,k).tstat;
        CDpktimeAVG(j,k)= mean([CDpktime1(j,k) CDpktime2(j,k) CDpktime3(j,k) CDpktime4(j,k) CDpktime5(j,k)]);
        CDpktimeSEM(j,k)= std([CDpktime1(j,k) CDpktime2(j,k) CDpktime3(j,k) CDpktime4(j,k) CDpktime5(j,k)])/sqrt(5);
        CDtrtimeAVG(j,k)= mean([CDtrtime1(j,k) CDtrtime2(j,k) CDtrtime3(j,k) CDtrtime4(j,k) CDtrtime5(j,k)]);
        CDtrtimeSEM(j,k)= std([CDtrtime1(j,k) CDtrtime2(j,k) CDtrtime3(j,k) CDtrtime4(j,k) CDtrtime5(j,k)])/sqrt(5);
        CVpktimeAVG(j,k)= mean([CVpktime1(j,k) CVpktime2(j,k) CVpktime3(j,k) CVpktime4(j,k) CVpktime5(j,k)]);
        CVpktimeSEM(j,k)= std([CVpktime1(j,k) CVpktime2(j,k) CVpktime3(j,k) CVpktime4(j,k) CVpktime5(j,k)])/sqrt(5);
        CVtrtimeAVG(j,k)= mean([CVtrtime1(j,k) CVtrtime2(j,k) CVtrtime3(j,k) CVtrtime4(j,k) CVtrtime5(j,k)]);
        CVtrtimeSEM(j,k)= std([CVtrtime1(j,k) CVtrtime2(j,k) CVtrtime3(j,k) CVtrtime4(j,k) CVtrtime5(j,k)])/sqrt(5);
        fftfreqstoAVG(j,k)= mean([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)]);
        fftfreqstoSEM(j,k)= std([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)])/sqrt(5);
        fftamplstoAVG(j,k)= mean([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)]);
        fftamplstoSEM(j,k)= std([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)])/sqrt(5);
        meanpktimeAVG(j,k)= mean([meanpktime1(j,k) meanpktime2(j,k) meanpktime3(j,k) meanpktime4(j,k) meanpktime5(j,k)]);
        meanpktimeSEM(j,k)= std([meanpktime1(j,k) meanpktime2(j,k) meanpktime3(j,k) meanpktime4(j,k) meanpktime5(j,k)])/sqrt(5);
        meantrtimeAVG(j,k)= mean([meantrtime1(j,k) meantrtime2(j,k) meantrtime3(j,k) meantrtime4(j,k) meantrtime5(j,k)]);
        meantrtimeSEM(j,k)= std([meantrtime1(j,k) meantrtime2(j,k) meantrtime3(j,k) meantrtime4(j,k) meantrtime5(j,k)])/sqrt(5);
        PKamplmeanstoAVG(j,k)= mean([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)]);
        PKamplmeanstoSEM(j,k)= std([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)])/sqrt(5);
    end
    %CDtstat(j,isinf(CDtstat(j,:)))=0;
    klp = findnearest(((CDtstat(j,:))),0);     transitionend(j)=klp(end);clear klp;
    transitionend_p(j) = CDpsto(j,transitionend(j));
end
%}

% coefficient of dispersion = std/mean -- between peak/peak and
% trough/trough

for j = 1:maxiter
    for k = 1:500
        CVstopk1(j,k)=sqrt(IEIvarpksto1(j,k))/IEImeanpksto1(j,k);
        CVstopk2(j,k)=sqrt(IEIvarpksto2(j,k))/IEImeanpksto2(j,k);
        CVstopk3(j,k)=sqrt(IEIvarpksto3(j,k))/IEImeanpksto3(j,k);
        CVstopk4(j,k)=sqrt(IEIvarpksto4(j,k))/IEImeanpksto4(j,k);
        CVstopk5(j,k)=sqrt(IEIvarpksto5(j,k))/IEImeanpksto5(j,k);
        CVstopkAVG(j,k)= mean([CVstopk1(j,k) CVstopk2(j,k) CVstopk3(j,k) CVstopk4(j,k) CVstopk5(j,k)]);
        CVstopkSEM(j,k)= std([CVstopk1(j,k) CVstopk2(j,k) CVstopk3(j,k) CVstopk4(j,k) CVstopk5(j,k)])/sqrt(5);
        IEIvarpkstoAVG(j,k) = mean([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]);
        IEIstdpkstoAVG(j,k) = sqrt(mean([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]));
        IEIvarpkstoSEM(j,k) = std([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)])/sqrt(5);
        IEIstdpkstoSEM(j,k) = sqrt(std([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]))/sqrt(5);
        IEImeanpkstoAVG(j,k) =  mean([IEImeanpksto1(j,k) IEImeanpksto2(j,k) IEImeanpksto3(j,k) IEImeanpksto4(j,k) IEImeanpksto5(j,k)]);
        IEImeanpkstoSEM(j,k) =  std([IEImeanpksto1(j,k) IEImeanpksto2(j,k) IEImeanpksto3(j,k) IEImeanpksto4(j,k) IEImeanpksto5(j,k)])/sqrt(5);
    end
end


% Find the averages for all of the deterministic cases
pkspikeratedet=mean(pkspikeratedet,1);trspikeratedet=mean(trspikeratedet,1);
CDdetpk=mean(CDdetpk,1);CDdettr=mean(CDdettr,1);CDpktimedet=mean(CDpktimedet,1);
CDtrtimedet=mean(CDtrtimedet,1);CVpktimedet=mean(CVpktimedet,1);CVtrtimedet=mean(CVtrtimedet,1);
fftampldet=mean(fftampldet,1);fftfreqdet=mean(fftfreqdet,1);IEImeanpkdet=mean(IEImeanpkdet,1);IEImeantrdet=mean(IEImeantrdet,1);
IEIvarpkdet=mean(IEIvarpkdet,1);IEIvartrdet=mean(IEIvartrdet,1);meanIEIpkdet=mean(meanIEIpkdet,1);
meanIEItrdet=mean(meanIEItrdet,1);meanpktimedet=mean(meanpktimedet,1);meantrtimedet=mean(meantrtimedet,1);
PKamplmeandet=mean(PKamplmeandet,1);pkdiffusiondet=mean(pkdiffusiondet,1);pkIEIspikeratiodet=mean(pkIEIspikeratiodet,1);
pktcorrdet=mean(pktcorrdet,1);trdiffusiondet=mean(trdiffusiondet,1);trtcorrdet=mean(trtcorrdet,1);

tr2=find(pkspikeratedet(1,:)==0);
tr2=max(tr2);
for j = 1:maxiter
    pkspikshift(j) = I(transitionend(j));
    I_shifted(j,:) = I - pkspikshift(j);
end

% PLOT THE SPIKE RATES
figure(1);
subplot(2,1,1);
plot(I,pkspikeratedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,pkspikeratestoAVG(j,:),pkspikeratestoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 max(max(pkspikeratestoAVG))]);
xlabel('Control Parameter');ylabel('Spike Rate (spikes/sec)');
title('SNIC Bifurcation; Averages = 5');

subplot(2,1,2);
plot(I,pkspikeratedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
     ha=errorbar(I_shifted(j,:),pkspikeratestoAVG(j,:),pkspikeratestoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 max(max(pkspikeratestoAVG))]);
xlabel('Control Parameter (shifted)');ylabel('Spike Rate (spikes/sec)');
title('SNIC Bifurcation; Averages = 5');


%{
ft=fittype('a*x^b+c');
[fitdet gofdet] = fit(I(find(I>0))',pkspikeratedet(1,find(I>0))',ft);
plot(fitdet,'k');
for j = 1:maxiter
    [fitsto{j} gofsto{j}] = fit(I(I>0)',pkspikeratestoAVG(j,I>0)',ft);
    plot(fitsto{j});hold all;;
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4',sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Deterministic','y=a*x^b+c','a=',num2str(fitdet.a),'b=',num2str(fitdet.b),'c=',num2str(fitdet.c),'R^2=',num2str(gofdet.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.05','y=a*x^b+c','a=',num2str(fitsto{1}.a),'b=',num2str(fitsto{1}.b),'c=',num2str(fitsto{1}.c),'R^2=',num2str(gofsto{1}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.1','y=a*x^b+c','a=',num2str(fitsto{2}.a),'b=',num2str(fitsto{2}.b),'c=',num2str(fitsto{2}.c),'R^2=',num2str(gofsto{2}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.2','y=a*x^b+c','a=',num2str(fitsto{3}.a),'b=',num2str(fitsto{3}.b),'c=',num2str(fitsto{3}.c),'R^2=',num2str(gofsto{3}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.4','y=a*x^b+c','a=',num2str(fitsto{4}.a),'b=',num2str(fitsto{4}.b),'c=',num2str(fitsto{4}.c),'R^2=',num2str(gofsto{4}.rsquare)));
axis([mu(1) mu(end) 0 max(max(pkspikeratestoAVG))]);
xlabel('Control Parameter');ylabel('Spike Rate (spikes/sec)');
title('SNIC Bifurcation; Averages = 5');
%}

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(2);
subplot(2,1,1);
plot(I,CDdetpk(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    hold all;
    plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVstopkAVG(j,:),CVstopkSEM(j,:));
    hb = get(ha,'children'); 
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(IEImeanpkstoAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - peak/peak');
title('SNIC Bifurcation - Between Peaks; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR PEAK/PEAK and TROUGH/TROUGH
figure(2)
subplot(2,1,2);
plot(I,CDdetpk(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDstopkAVG(j,:),CDstopkSEM(j,:));
    hb = get(ha,'children'); 
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - peak/peak');
title('SNIC Bifurcation - Between Peaks; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR EACH PEAK
figure(4)
subplot(2,2,1);
plot(I,CDpktimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(2,2,1);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDpktimeAVG(j,:),CDpktimeSEM(j,:));
    hb = get(ha,'children');
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH PEAK');
title('SNIC Bifurcation - PEAKS; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR EACH TROUGH
figure(4)
subplot(2,2,2);
plot(I,CDtrtimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    subplot(2,2,2);
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDtrtimeAVG(j,:),CDtrtimeSEM(j,:));
    hb = get(ha,'children');
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH TROUGH');
title('SNIC Bifurcation - TROUGHS; Averages = 5');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(4);
subplot(2,2,3);
plot(I,CVpktimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    subplot(2,2,3);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVpktimeAVG(j,:),CVpktimeSEM(j,:));
    hb = get(ha,'children');
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(meanpktimeAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH PEAK');
title('SNIC Bifurcation - Each Peak; Averages = 5');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(4);
subplot(2,2,4);
plot(I,CVtrtimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    subplot(2,2,4);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVtrtimeAVG(j,:),CVtrtimeSEM(j,:));
    hb = get(ha,'children');  
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(meantrtimeAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH TROUGH');
title('SNIC Bifurcation - Each Trough; Averages = 5');

% PLOT THE AMPLITUDES FROM FFT
figure(5);
subplot(3,1,2);
plot(I,fftampldet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,2);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,fftamplstoAVG(j,:),fftamplstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.2*fftampldet(1,end)]);
xlabel('Control Parameter');ylabel('Amplitude from FFT');
title('SNIC Bifurcation; Averages = 5');

% PLOT THE FREQUENCIES FROM FFT
figure(5);
subplot(3,1,3);
plot(I,fftfreqdet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,3);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,fftfreqstoAVG(j,:),fftfreqstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.2*fftfreqdet(1,end)]);
xlabel('Control Parameter');ylabel('Frequency from FFT');
title('SNIC Bifurcation; Averages = 5');

% PLOT THE PEAK-TO-TROUGH AMPLITUDES
figure(5);
subplot(3,1,1);
plot(I,PKamplmeandet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,1);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,PKamplmeanstoAVG(j,:),PKamplmeanstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.2*PKamplmeandet(1,end)]);
xlabel('Control Parameter');ylabel('Peak-to-Trough Amplitude');
title('SNIC Bifurcation; Averages = 5');


Iselected = [-0.15 -0.1 -0.05 -0.01 -0.005 0.005 0.01 0.05 0.1 0.15]; % CHOOSE
clear Iselect
for k = 1:length(Iselected)
    Iselect(k) = max(findnearest(I,Iselected(k)));
end
tmin=10000;tmax=15000;    % CHOOSE
ymin=-2;ymax=2;
yimin=ymin*2;yimax=ymax*2;
dwnsplquiver=20;quiverstart=2000;quiverend=4000;quiverscale=1;
dwnsplrealimag=10;realimagstart=2000;realimagend=4000;
fmin = 0.05; fmax = 0.3;

% PLOT EXAMPLE TIME TRACES AND PHASE PORTRAITS
% Deterministic
xdr=real(hilbert(Xdet{1}{Iselect(end)}));xdi=imag(hilbert(Xdet{1}{Iselect(end)}));
[bw dens mx1 my1]=kde2d([xdr',xdi']);
for k = 1:length(Iselect);
    figure(6);
    sph=subplot(6,length(Iselect),k);plot(Xdet{1}{Iselect(k)},'k');axis([tmin tmax ymin ymax]);
    set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s %s%s','Deterministic','I = ',num2str(I(Iselect(k)))));
    spp = get(sph, 'pos');
    set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xdet{1}{Iselect(k)}));xdi=imag(hilbert(Xdet{1}{Iselect(k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1
        sph=subplot(6,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        if length(my(:,1)) < size(dens,1)
            my(end:size(dens,1),:)=0;
        elseif length(my(:,1)) > size(dens,1)
            dens(end:length(my(:,1)),:)=0;
        end
        if length(mx(1,:)) < size(dens,2)
            mx(end:size(dens,2))=0;
        elseif length(mx(1,:)) > size(dens,2)
            dens(:,end:length(mx(1,:)))=0;
        end
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(6,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 1.1*max(dens1)]);axis([mx(1,1) mx(1,end) 0 0.001]);
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+5*length(Iselect));plot(fdetfft{j,Iselect(k)},Xdetfft{1,Iselect(k)},'k');axis([fmin fmax 0 1.05*max(Xdetfft{1,Iselect(end)})]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %tt=findnearest(fdetfft{j,Iselect(k)},fftfreqdet(1,Iselect(k)));tt=tt(1);hold on;scatter(fdetfft{j,Iselect(k)}(tt),Xdetfft{1,Iselect(k)}(tt),'b.');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([yimin yimax yimin yimax]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([yimin yimax yimin yimax]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end

% Stochastic
xdr=real(hilbert(Xsto{j}{Iselect(end)}));xdi=imag(hilbert(Xsto{j}{Iselect(end)}));
[bw dens mx1 my1]=kde2d([xdr',xdi']);
for j = 1:maxiter
for k = 1:length(Iselect);
    figure(6+j);
    sph=subplot(6,length(Iselect),k);plot(Xsto{j}{Iselect(k)},'r');axis([tmin+500*j tmax+500*j ymin ymax]);
    set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s%s %s%s','Noise = ',num2str(noiselevel(j)),', I = ',num2str(I(Iselect(k)))));
    spp = get(sph, 'pos');
    set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xsto{j}{Iselect(k)}));xdi=imag(hilbert(Xsto{j}{Iselect(k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1
        sph=subplot(5,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        if length(my(:,1)) < size(dens,1)
            my(end:size(dens,1),:)=0;
        elseif length(my(:,1)) > size(dens,1)
            dens(end:length(my(:,1)),:)=0;
        end
        if length(mx(1,:)) < size(dens,2)
            mx(end:size(dens,2))=0;
        elseif length(mx(1,:)) > size(dens,2)
            dens(:,end:length(mx(1,:)))=0;
        end
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(6,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 0.001]);axis([mx(1,1) mx(1,end) 0 0.001]);
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+5*length(Iselect));plot(fstofft{j,Iselect(k)},Xstofft{j,Iselect(k)},'r');axis([fmin fmax 0 1.05*max(Xstofft{j,Iselect(end)})]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %tt=findnearest(fstofft{j,Iselect(k)},fftfreqsto(1,Iselect(k)));tt=tt(1);hold on;scatter(fstofft{j,Iselect(k)}(tt),Xstofft{1,Iselect(k)}(tt),'b.');        
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([yimin yimax yimin yimax]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([yimin yimax yimin yimax]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end
end

elseif biftype ==1 
% SUPERCRITICAL HOPF
close all;
maxiter = 4;
noiselevel = [0.05 0.1 0.2 0.4];

for j = 1:maxiter
    for k = 1:500
        pkspikeratestoAVG(j,k)= mean([pkspikeratesto1(j,k) pkspikeratesto2(j,k) pkspikeratesto3(j,k) pkspikeratesto4(j,k) pkspikeratesto5(j,k)]);
        pkspikeratestoSEM(j,k)= std([pkspikeratesto1(j,k) pkspikeratesto2(j,k) pkspikeratesto3(j,k) pkspikeratesto4(j,k) pkspikeratesto5(j,k)])/sqrt(5);
    end
end


% coefficient of variation 1 = variance/mean
for j = 1:maxiter
    for k = 1:500
        CDstopkAVG(j,k)= mean([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)]);
        CDstopkSEM(j,k)= std([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)])/sqrt(5);
    end
end
%}

% FFT and amplitudes from peak/trough
for j = 1:maxiter
    for k = 1:500
        [~,CDpsto(j,k),~,CDstatsto(j,k)] =  ttest([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)],1);CDtstat(j,k)=CDstatsto(j,k).tstat;
        CDpktimeAVG(j,k)= mean([CDpktime1(j,k) CDpktime2(j,k) CDpktime3(j,k) CDpktime4(j,k) CDpktime5(j,k)]);
        CDpktimeSEM(j,k)= std([CDpktime1(j,k) CDpktime2(j,k) CDpktime3(j,k) CDpktime4(j,k) CDpktime5(j,k)])/sqrt(5);
        CDtrtimeAVG(j,k)= mean([CDtrtime1(j,k) CDtrtime2(j,k) CDtrtime3(j,k) CDtrtime4(j,k) CDtrtime5(j,k)]);
        CDtrtimeSEM(j,k)= std([CDtrtime1(j,k) CDtrtime2(j,k) CDtrtime3(j,k) CDtrtime4(j,k) CDtrtime5(j,k)])/sqrt(5);
        CVpktimeAVG(j,k)= mean([CVpktime1(j,k) CVpktime2(j,k) CVpktime3(j,k) CVpktime4(j,k) CVpktime5(j,k)]);
        CVpktimeSEM(j,k)= std([CVpktime1(j,k) CVpktime2(j,k) CVpktime3(j,k) CVpktime4(j,k) CVpktime5(j,k)])/sqrt(5);
        CVtrtimeAVG(j,k)= mean([CVtrtime1(j,k) CVtrtime2(j,k) CVtrtime3(j,k) CVtrtime4(j,k) CVtrtime5(j,k)]);
        CVtrtimeSEM(j,k)= std([CVtrtime1(j,k) CVtrtime2(j,k) CVtrtime3(j,k) CVtrtime4(j,k) CVtrtime5(j,k)])/sqrt(5);
        fftfreqstoAVG(j,k)= mean([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)]);
        fftfreqstoSEM(j,k)= std([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)])/sqrt(5);
        fftamplstoAVG(j,k)= mean([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)]);
        fftamplstoSEM(j,k)= std([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)])/sqrt(5);
        meanpktimeAVG(j,k)= mean([meanpktime1(j,k) meanpktime2(j,k) meanpktime3(j,k) meanpktime4(j,k) meanpktime5(j,k)]);
        meanpktimeSEM(j,k)= std([meanpktime1(j,k) meanpktime2(j,k) meanpktime3(j,k) meanpktime4(j,k) meanpktime5(j,k)])/sqrt(5);
        meantrtimeAVG(j,k)= mean([meantrtime1(j,k) meantrtime2(j,k) meantrtime3(j,k) meantrtime4(j,k) meantrtime5(j,k)]);
        meantrtimeSEM(j,k)= std([meantrtime1(j,k) meantrtime2(j,k) meantrtime3(j,k) meantrtime4(j,k) meantrtime5(j,k)])/sqrt(5);
        PKamplmeanstoAVG(j,k)= mean([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)]);
        PKamplmeanstoSEM(j,k)= std([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)])/sqrt(5);
    end
    %CDtstat(j,isinf(CDtstat(j,:)))=0;
    klp = findnearest(((CDtstat(j,:))),0);
    transitionend(j)=klp(end);clear klp
    transitionend_p(j) = CDpsto(j,transitionend(j));
end

% coefficient of dispersion = std/mean

for j = 1:maxiter
    for k = 1:500
        CVstopk1(j,k)=sqrt(IEIvarpksto1(j,k))/IEImeanpksto1(j,k);
        CVstopk2(j,k)=sqrt(IEIvarpksto2(j,k))/IEImeanpksto2(j,k);
        CVstopk3(j,k)=sqrt(IEIvarpksto3(j,k))/IEImeanpksto3(j,k);
        CVstopk4(j,k)=sqrt(IEIvarpksto4(j,k))/IEImeanpksto4(j,k);
        CVstopk5(j,k)=sqrt(IEIvarpksto5(j,k))/IEImeanpksto5(j,k);
        CVstopkAVG(j,k)= mean([CVstopk1(j,k) CVstopk2(j,k) CVstopk3(j,k) CVstopk4(j,k) CVstopk5(j,k)]);
        CVstopkSEM(j,k)= std([CVstopk1(j,k) CVstopk2(j,k) CVstopk3(j,k) CVstopk4(j,k) CVstopk5(j,k)])/sqrt(5);
        IEIvarpkstoAVG(j,k) = mean([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]);
        IEIstdpkstoAVG(j,k) = sqrt(mean([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]));
        IEIvarpkstoSEM(j,k) = std([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)])/sqrt(5);
        IEIstdpkstoSEM(j,k) = sqrt(std([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]))/sqrt(5);
        IEImeanpkstoAVG(j,k) =  mean([IEImeanpksto1(j,k) IEImeanpksto2(j,k) IEImeanpksto3(j,k) IEImeanpksto4(j,k) IEImeanpksto5(j,k)]);
        IEImeanpkstoSEM(j,k) =  std([IEImeanpksto1(j,k) IEImeanpksto2(j,k) IEImeanpksto3(j,k) IEImeanpksto4(j,k) IEImeanpksto5(j,k)])/sqrt(5);
    end
end

% Find the averages for all of the deterministic cases
pkspikeratedet=mean(pkspikeratedet,1);trspikeratedet=mean(trspikeratedet,1);
CDdetpk=mean(CDdetpk,1);CDdettr=mean(CDdettr,1);CDpktimedet=mean(CDpktimedet,1);
CDtrtimedet=mean(CDtrtimedet,1);CVpktimedet=mean(CVpktimedet,1);CVtrtimedet=mean(CVtrtimedet,1);
fftampldet=mean(fftampldet,1);fftfreqdet=mean(fftfreqdet,1);IEImeanpkdet=mean(IEImeanpkdet,1);IEImeantrdet=mean(IEImeantrdet,1);
IEIvarpkdet=mean(IEIvarpkdet,1);IEIvartrdet=mean(IEIvartrdet,1);meanIEIpkdet=mean(meanIEIpkdet,1);
meanIEItrdet=mean(meanIEItrdet,1);meanpktimedet=mean(meanpktimedet,1);meantrtimedet=mean(meantrtimedet,1);
PKamplmeandet=mean(PKamplmeandet,1);pkdiffusiondet=mean(pkdiffusiondet,1);pkIEIspikeratiodet=mean(pkIEIspikeratiodet,1);
pktcorrdet=mean(pktcorrdet,1);trdiffusiondet=mean(trdiffusiondet,1);trtcorrdet=mean(trtcorrdet,1);

tr2=find(pkspikeratedet(1,:)==0);
tr2=max(tr2);
for j = 1:maxiter
    pkspikshift(j) = I(transitionend(j));
    I_shifted(j,:) = I - pkspikshift(j);
end

figure(1);
subplot(2,1,1);
plot(mu,pkspikeratedet(1,:),'k');I=mu;
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    %plot(mu,pkspikeratestoAVG(j,:));
    ha=errorbar(mu,pkspikeratestoAVG(j,:),pkspikeratestoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.2;
    Xdata(xright) = Xdata(xright)-0.2;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 2*max(max(pkspikeratesto))]);
xlabel('Control Parameter');ylabel('Spike Rate (spikes/sec)');
title('Hopf Bifurcation; Averages = 5');

subplot(2,1,2);
plot(I,pkspikeratedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
     ha=errorbar(I_shifted(j,:),pkspikeratestoAVG(j,:),pkspikeratestoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 max(max(pkspikeratestoAVG))]);
xlabel('Control Parameter (shifted)');ylabel('Spike Rate (spikes/sec)');
title('Hopf Bifurcation; Averages = 5');

%{
ft=fittype('a*x^b+c');
[fitdet gofdet] = fit(I(I>0)',pkspikeratedet(1,I>0)',ft);
plot(fitdet,'k');
for j = 1:maxiter
    [fitsto{j} gofsto{j}] = fit(I(I>0)',pkspikeratestoAVG(j,I>0)',ft);
    plot(fitsto{j});hold all;;
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4',sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Deterministic','y=a*x^b+c','a=',num2str(fitdet.a),'b=',num2str(fitdet.b),'c=',num2str(fitdet.c),'R^2=',num2str(gofdet.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.05','y=a*x^b+c','a=',num2str(fitsto{1}.a),'b=',num2str(fitsto{1}.b),'c=',num2str(fitsto{1}.c),'R^2=',num2str(gofsto{1}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.1','y=a*x^b+c','a=',num2str(fitsto{2}.a),'b=',num2str(fitsto{2}.b),'c=',num2str(fitsto{2}.c),'R^2=',num2str(gofsto{2}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.2','y=a*x^b+c','a=',num2str(fitsto{3}.a),'b=',num2str(fitsto{3}.b),'c=',num2str(fitsto{3}.c),'R^2=',num2str(gofsto{3}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.4','y=a*x^b+c','a=',num2str(fitsto{4}.a),'b=',num2str(fitsto{4}.b),'c=',num2str(fitsto{4}.c),'R^2=',num2str(gofsto{4}.rsquare)));
axis([mu(1) mu(end) 0 2*max(max(pkspikeratesto))]);
xlabel('Control Parameter');ylabel('Spike Rate (spikes/sec)');
title('Hopf Bifurcation; Averages = 5');
%}

figure(2);
subplot(2,1,1);
plot(mu,CDdetpk(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(mu,CVstopkAVG(j,:),CVstopkSEM(j,:));
    set(get(ha,'Parent'),'YScale','log');
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(mu,1./sqrt(IEImeanpkstoAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(mu,ones(1,length(mu)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 1e-2 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation');
title('Hopf Bifurcation; Averages = 5');

figure(2);
subplot(2,1,2);
plot(mu,CDdetpk(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(mu,CDstopkAVG(j,:),CDstopkSEM(j,:));
    set(get(ha,'Parent'),'YScale','log');
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
plot(mu,ones(1,length(mu)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 1e-2 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion');
title('Hopf Bifurcation; Averages = 5');


% PLOT THE COEFFICIENT OF DISPERSION FOR EACH PEAK
figure(4)
subplot(2,2,1);
plot(I,CDpktimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(2,2,1);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDpktimeAVG(j,:),CDpktimeSEM(j,:));
    set(get(ha,'Parent'),'YScale','log');
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 1e-2 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH PEAK');
title('Hopf Bifurcation - PEAKS; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR EACH TROUGH
figure(4)
subplot(2,2,2);
plot(I,CDtrtimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    subplot(2,2,2);
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDtrtimeAVG(j,:),CDtrtimeSEM(j,:));
    set(get(ha,'Parent'),'YScale','log');
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 1e-2 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH TROUGH');
title('Hopf Bifurcation - TROUGHS; Averages = 5');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(4);
subplot(2,2,3);
plot(I,CVpktimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    subplot(2,2,3);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVpktimeAVG(j,:),CVpktimeSEM(j,:));
    set(get(ha,'Parent'),'YScale','log');
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(meanpktimeAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 1e-2 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH PEAK');
title('Hopf Bifurcation - Each Peak; Averages = 5');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(4);
subplot(2,2,4);
plot(I,CVtrtimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    subplot(2,2,4);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVtrtimeAVG(j,:),CVtrtimeSEM(j,:));
    set(get(ha,'Parent'),'YScale','log');
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(meantrtimeAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 1e-2 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH TROUGH');
title('Hopf Bifurcation - Each Trough; Averages = 5');

% PLOT THE AMPLITUDES FROM FFT
figure(5);
subplot(3,1,2);
plot(I,fftampldet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,2);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,fftamplstoAVG(j,:),fftamplstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.1*max(max(fftamplsto))]);
xlabel('Control Parameter');ylabel('Amplitude from FFT');
title('Hopf Bifurcation; Averages = 5');

% PLOT THE FREQUENCIES FROM FFT
figure(5);
subplot(3,1,3);
plot(I,fftfreqdet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,3);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,fftfreqstoAVG(j,:),fftfreqstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.1*max(max(fftfreqsto))]);
xlabel('Control Parameter');ylabel('Frequency from FFT');
title('Hopf Bifurcation; Averages = 5');

% PLOT THE PEAK-TO-TROUGH AMPLITUDES
figure(5);
subplot(3,1,1);
plot(I,PKamplmeandet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,1);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,PKamplmeanstoAVG(j,:),PKamplmeanstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.1*max(max(PKamplmeansto))]);
xlabel('Control Parameter');ylabel('Peak-to-Trough Amplitude');
title('Hopf Bifurcation; Averages = 5');



Iselected = [17 18.5 19.5 20 20.5 20.7 20.8 20.9 21.5 22]; % CHOOSE
clear Iselect
for k = 1:length(Iselected)
    Iselect(k) = max(findnearest(I,Iselected(k)));
    maxx(k) = max([max(abs(Xsto{1}{Iselect(k)})) max(abs(Xsto{2}{Iselect(k)})) max(abs(Xsto{3}{Iselect(k)})) max(abs(Xsto{4}{Iselect(k)}))]);
end
maxy=max(maxx);
tmin=10000;tmax=15000;    % CHOOSE
ymin=-1.2*maxy;ymax=1.2*maxy;
yimin=ymin*2;yimax=ymax*2;
dwnsplquiver=1;quiverstart=2000;quiverend=2050;quiverscale=1;
dwnsplrealimag=1;realimagstart=2000;realimagend=4000;
fmin = 0.05; fmax = 3;

% PLOT EXAMPLE TIME TRACES AND PHASE PORTRAITS
% Deterministic
xdr=real(hilbert(Xdet{1}{Iselect(end)}));xdi=imag(hilbert(Xdet{1}{Iselect(end)}));
[bw dens mx1 my1]=kde2d([xdr',xdi']);
for k = 1:length(Iselect);
    figure(6);
    sph=subplot(6,length(Iselect),k);plot(Xdet{1}{Iselect(k)},'k');axis([tmin tmax ymin ymax]);set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s %s%s','Deterministic','I = ',num2str(I(Iselect(k)))));
    spp = get(sph, 'pos');
    set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xdet{1}{Iselect(k)}));xdi=imag(hilbert(Xdet{1}{Iselect(k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1
        sph=subplot(6,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        if length(my(:,1)) < size(dens,1)
            my(end:size(dens,1),:)=0;
        elseif length(my(:,1)) > size(dens,1)
            dens(end:length(my(:,1)),:)=0;
        end
        if length(mx(1,:)) < size(dens,2)
            mx(end:size(dens,2))=0;
        elseif length(mx(1,:)) > size(dens,2)
            dens(:,end:length(mx(1,:)))=0;
        end
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(6,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 1.1*max(dens1)]);axis([mx(1,1) mx(1,end) 0 0.001]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+5*length(Iselect));plot(fdetfft{j,Iselect(k)},Xdetfft{1,Iselect(k)},'k');axis([fmin fmax 0 1.05*max(Xdetfft{1,Iselect(end)})]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %tt=findnearest(fdetfft{j,Iselect(k)},fftfreqdet(1,Iselect(k)));tt=tt(1);hold on;scatter(fdetfft{j,Iselect(k)}(tt),Xdetfft{1,Iselect(k)}(tt),'b.');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([yimin yimax yimin yimax]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([yimin yimax yimin yimax]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end

% Stochastic
xdr=real(hilbert(Xsto{j}{Iselect(end)}));xdi=imag(hilbert(Xsto{j}{Iselect(end)}));
[bw dens mx1 my1]=kde2d([xdr',xdi']);
for j = 1:maxiter
for k = 1:length(Iselect);
    figure(6+j);
    sph=subplot(6,length(Iselect),k);plot(Xsto{j}{Iselect(k)},'r');axis([tmin+500*j tmax+500*j ymin ymax]);
    set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s%s %s%s','Noise = ',num2str(noiselevel(j)),', I = ',num2str(I(Iselect(k)))));
    spp = get(sph, 'pos');
    set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xsto{j}{Iselect(k)}));xdi=imag(hilbert(Xsto{j}{Iselect(k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1
        sph=subplot(5,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        if length(my(:,1)) < size(dens,1)
            my(end:size(dens,1),:)=0;
        elseif length(my(:,1)) > size(dens,1)
            dens(end:length(my(:,1)),:)=0;
        end
        if length(mx(1,:)) < size(dens,2)
            mx(end:size(dens,2))=0;
        elseif length(mx(1,:)) > size(dens,2)
            dens(:,end:length(mx(1,:)))=0;
        end
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(6,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 0.001]);axis([mx(1,1) mx(1,end) 0 0.0005]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+5*length(Iselect));plot(fstofft{j,Iselect(k)},Xstofft{j,Iselect(k)},'r');axis([fmin fmax 0 1.05*max(Xstofft{j,Iselect(end)})]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %tt=findnearest(fstofft{j,Iselect(k)},fftfreqsto(1,Iselect(k)));tt=tt(1);hold on;scatter(fstofft{j,Iselect(k)}(tt),Xstofft{1,Iselect(k)}(tt),'b.');        
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([yimin yimax yimin yimax]);
        set(gca, 'XTickLabel', []);set(gca, 'YTickLabel', []);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end
end

elseif biftype ==3
% SUBCRITICAL HOPF
close all;
maxiter = 4;
noiselevel = [0.05 0.1 0.2 0.4];

for j = 1:maxiter
    for k = 1:500
        pkspikeratestoAVG(j,k)= mean([pkspikeratesto1(j,k) pkspikeratesto2(j,k) pkspikeratesto3(j,k) pkspikeratesto4(j,k) pkspikeratesto5(j,k)]);
        pkspikeratestoSEM(j,k)= std([pkspikeratesto1(j,k) pkspikeratesto2(j,k) pkspikeratesto3(j,k) pkspikeratesto4(j,k) pkspikeratesto5(j,k)])/sqrt(5);
    end
end


% coefficient of variation 1 = variance/mean
for j = 1:maxiter
    for k = 1:500
        CDstopkAVG(j,k)= mean([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)]);
        CDstopkSEM(j,k)= std([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)])/sqrt(5);
    end
end
%}

% FFT and amplitudes from peak/trough
for j = 1:maxiter
    for k = 1:500
        [~,CDpsto(j,k),~,CDstatsto(j,k)] =  ttest([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)],1);CDtstat(j,k)=CDstatsto(j,k).tstat;
        CDpktimeAVG(j,k)= mean([CDpktime1(j,k) CDpktime2(j,k) CDpktime3(j,k) CDpktime4(j,k) CDpktime5(j,k)]);
        CDpktimeSEM(j,k)= std([CDpktime1(j,k) CDpktime2(j,k) CDpktime3(j,k) CDpktime4(j,k) CDpktime5(j,k)])/sqrt(5);
        CDtrtimeAVG(j,k)= mean([CDtrtime1(j,k) CDtrtime2(j,k) CDtrtime3(j,k) CDtrtime4(j,k) CDtrtime5(j,k)]);
        CDtrtimeSEM(j,k)= std([CDtrtime1(j,k) CDtrtime2(j,k) CDtrtime3(j,k) CDtrtime4(j,k) CDtrtime5(j,k)])/sqrt(5);
        CVpktimeAVG(j,k)= mean([CVpktime1(j,k) CVpktime2(j,k) CVpktime3(j,k) CVpktime4(j,k) CVpktime5(j,k)]);
        CVpktimeSEM(j,k)= std([CVpktime1(j,k) CVpktime2(j,k) CVpktime3(j,k) CVpktime4(j,k) CVpktime5(j,k)])/sqrt(5);
        CVtrtimeAVG(j,k)= mean([CVtrtime1(j,k) CVtrtime2(j,k) CVtrtime3(j,k) CVtrtime4(j,k) CVtrtime5(j,k)]);
        CVtrtimeSEM(j,k)= std([CVtrtime1(j,k) CVtrtime2(j,k) CVtrtime3(j,k) CVtrtime4(j,k) CVtrtime5(j,k)])/sqrt(5);
        fftfreqstoAVG(j,k)= mean([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)]);
        fftfreqstoSEM(j,k)= std([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)])/sqrt(5);
        fftamplstoAVG(j,k)= mean([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)]);
        fftamplstoSEM(j,k)= std([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)])/sqrt(5);
        meanpktimeAVG(j,k)= mean([meanpktime1(j,k) meanpktime2(j,k) meanpktime3(j,k) meanpktime4(j,k) meanpktime5(j,k)]);
        meanpktimeSEM(j,k)= std([meanpktime1(j,k) meanpktime2(j,k) meanpktime3(j,k) meanpktime4(j,k) meanpktime5(j,k)])/sqrt(5);
        meantrtimeAVG(j,k)= mean([meantrtime1(j,k) meantrtime2(j,k) meantrtime3(j,k) meantrtime4(j,k) meantrtime5(j,k)]);
        meantrtimeSEM(j,k)= std([meantrtime1(j,k) meantrtime2(j,k) meantrtime3(j,k) meantrtime4(j,k) meantrtime5(j,k)])/sqrt(5);
        PKamplmeanstoAVG(j,k)= mean([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)]);
        PKamplmeanstoSEM(j,k)= std([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)])/sqrt(5);
    end
    %CDtstat(j,isinf(CDtstat(j,:)))=0;
    klp = findnearest(((CDtstat(j,:))),0);
    transitionend(j)=klp(end);clear klp
    transitionend_p(j) = CDpsto(j,transitionend(j));
end

% coefficient of dispersion = std/mean

for j = 1:maxiter
    for k = 1:500
        CVstopk1(j,k)=sqrt(IEIvarpksto1(j,k))/IEImeanpksto1(j,k);
        CVstopk2(j,k)=sqrt(IEIvarpksto2(j,k))/IEImeanpksto2(j,k);
        CVstopk3(j,k)=sqrt(IEIvarpksto3(j,k))/IEImeanpksto3(j,k);
        CVstopk4(j,k)=sqrt(IEIvarpksto4(j,k))/IEImeanpksto4(j,k);
        CVstopk5(j,k)=sqrt(IEIvarpksto5(j,k))/IEImeanpksto5(j,k);
        CVstopkAVG(j,k)= mean([CVstopk1(j,k) CVstopk2(j,k) CVstopk3(j,k) CVstopk4(j,k) CVstopk5(j,k)]);
        CVstopkSEM(j,k)= std([CVstopk1(j,k) CVstopk2(j,k) CVstopk3(j,k) CVstopk4(j,k) CVstopk5(j,k)])/sqrt(5);
        IEIvarpkstoAVG(j,k) = mean([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]);
        IEIstdpkstoAVG(j,k) = sqrt(mean([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]));
        IEIvarpkstoSEM(j,k) = std([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)])/sqrt(5);
        IEIstdpkstoSEM(j,k) = sqrt(std([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]))/sqrt(5);
        IEImeanpkstoAVG(j,k) =  mean([IEImeanpksto1(j,k) IEImeanpksto2(j,k) IEImeanpksto3(j,k) IEImeanpksto4(j,k) IEImeanpksto5(j,k)]);
        IEImeanpkstoSEM(j,k) =  std([IEImeanpksto1(j,k) IEImeanpksto2(j,k) IEImeanpksto3(j,k) IEImeanpksto4(j,k) IEImeanpksto5(j,k)])/sqrt(5);
    end
end

% Find the averages for all of the deterministic cases
pkspikeratedet=mean(pkspikeratedet,1);trspikeratedet=mean(trspikeratedet,1);
CDdetpk=mean(CDdetpk,1);CDdettr=mean(CDdettr,1);CDpktimedet=mean(CDpktimedet,1);
CDtrtimedet=mean(CDtrtimedet,1);CVpktimedet=mean(CVpktimedet,1);CVtrtimedet=mean(CVtrtimedet,1);
fftampldet=mean(fftampldet,1);fftfreqdet=mean(fftfreqdet,1);IEImeanpkdet=mean(IEImeanpkdet,1);IEImeantrdet=mean(IEImeantrdet,1);
IEIvarpkdet=mean(IEIvarpkdet,1);IEIvartrdet=mean(IEIvartrdet,1);meanIEIpkdet=mean(meanIEIpkdet,1);
meanIEItrdet=mean(meanIEItrdet,1);meanpktimedet=mean(meanpktimedet,1);meantrtimedet=mean(meantrtimedet,1);
PKamplmeandet=mean(PKamplmeandet,1);pkdiffusiondet=mean(pkdiffusiondet,1);pkIEIspikeratiodet=mean(pkIEIspikeratiodet,1);
pktcorrdet=mean(pktcorrdet,1);trdiffusiondet=mean(trdiffusiondet,1);trtcorrdet=mean(trtcorrdet,1);

tr2=find(pkspikeratedet(1,:)==0);
tr2=max(tr2);
for j = 1:maxiter
    pkspikshift(j) = I(transitionend(j));
    I_shifted(j,:) = I - pkspikshift(j);
end

figure(1);
subplot(2,1,1);
plot(mu,pkspikeratedet(1,:),'k');I=mu;
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    %plot(mu,pkspikeratestoAVG(j,:));
    ha=errorbar(mu,pkspikeratestoAVG(j,:),pkspikeratestoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 2*max(max(pkspikeratesto))]);
xlabel('Control Parameter');ylabel('Spike Rate (spikes/sec)');
title('Subcritical Hopf Bifurcation; Averages = 5');

subplot(2,1,2);
plot(I,pkspikeratedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
     ha=errorbar(I_shifted(j,:),pkspikeratestoAVG(j,:),pkspikeratestoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 max(max(pkspikeratestoAVG))]);
xlabel('Control Parameter (shifted)');ylabel('Spike Rate (spikes/sec)');
title('Subcritical Hopf Bifurcation; Averages = 5');

%{
ft=fittype('a*x^b+c');
[fitdet gofdet] = fit(I((I>0))',pkspikeratedet(1,I>0)',ft);
plot(fitdet,'k');
for j = 1:maxiter
    [fitsto{j} gofsto{j}] = fit(I(I>0)',pkspikeratestoAVG(j,(I>0))',ft);
    plot(fitsto{j});hold all;
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4',sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Deterministic','y=a*x^b+c','a=',num2str(fitdet.a),'b=',num2str(fitdet.b),'c=',num2str(fitdet.c),'R^2=',num2str(gofdet.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.05','y=a*x^b+c','a=',num2str(fitsto{1}.a),'b=',num2str(fitsto{1}.b),'c=',num2str(fitsto{1}.c),'R^2=',num2str(gofsto{1}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.1','y=a*x^b+c','a=',num2str(fitsto{2}.a),'b=',num2str(fitsto{2}.b),'c=',num2str(fitsto{2}.c),'R^2=',num2str(gofsto{2}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.2','y=a*x^b+c','a=',num2str(fitsto{3}.a),'b=',num2str(fitsto{3}.b),'c=',num2str(fitsto{3}.c),'R^2=',num2str(gofsto{3}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.4','y=a*x^b+c','a=',num2str(fitsto{4}.a),'b=',num2str(fitsto{4}.b),'c=',num2str(fitsto{4}.c),'R^2=',num2str(gofsto{4}.rsquare)));
axis([mu(1) mu(end) 0 1.5*pkspikeratedet(1,findnearest(I,1))]);
xlabel('Control Parameter');ylabel('Spike Rate (spikes/sec)');
title('Subcritical Hopf Bifurcation; Averages = 5');
%}

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(2);
subplot(2,1,1);
plot(I,CDdetpk(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    hold all;
    plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVstopkAVG(j,:),CVstopkSEM(j,:));
    hb = get(ha,'children'); 
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(IEImeanpkstoAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - peak/peak');
title('Subcritical Hopf Bifurcation - Between Peaks; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR PEAK/PEAK and TROUGH/TROUGH
figure(2)
subplot(2,1,2);
plot(I,CDdetpk(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDstopkAVG(j,:),CDstopkSEM(j,:));
    hb = get(ha,'children'); 
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - peak/peak');
title('Subcritical Hopf Bifurcation - Between Peaks; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR EACH PEAK
figure(4)
subplot(2,2,1);
plot(I,CDpktimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(2,2,1);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDpktimeAVG(j,:),CDpktimeSEM(j,:));
    hb = get(ha,'children');
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH PEAK');
title('Subcritical Hopf Bifurcation - PEAKS; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR EACH TROUGH
figure(4)
subplot(2,2,2);
plot(I,CDtrtimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    subplot(2,2,2);
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDtrtimeAVG(j,:),CDtrtimeSEM(j,:));
    hb = get(ha,'children');
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH TROUGH');
title('Subcritical Hopf Bifurcation - TROUGHS; Averages = 5');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(4);
subplot(2,2,3);
plot(I,CVpktimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    subplot(2,2,3);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVpktimeAVG(j,:),CVpktimeSEM(j,:));
    hb = get(ha,'children');
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(meanpktimeAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH PEAK');
title('Subcritical Hopf Bifurcation - Each Peak; Averages = 5');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(4);
subplot(2,2,4);
plot(I,CVtrtimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    subplot(2,2,4);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVtrtimeAVG(j,:),CVtrtimeSEM(j,:));
    hb = get(ha,'children');  
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(meantrtimeAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH TROUGH');
title('Subcritical Hopf Bifurcation - Each Trough; Averages = 5');

% PLOT THE AMPLITUDES FROM FFT
figure(5);
subplot(3,1,2);
plot(I,fftampldet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,2);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,fftamplstoAVG(j,:),fftamplstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.2*fftampldet(1,end)]);
xlabel('Control Parameter');ylabel('Amplitude from FFT');
title('Subcritical Hopf Bifurcation; Averages = 5');

% PLOT THE FREQUENCIES FROM FFT
figure(5);
subplot(3,1,3);
plot(I,fftfreqdet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,3);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,fftfreqstoAVG(j,:),fftfreqstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.2*fftfreqdet(1,end)]);
xlabel('Control Parameter');ylabel('Frequency from FFT');
title('Subcritical Hopf Bifurcation; Averages = 5');

% PLOT THE PEAK-TO-TROUGH AMPLITUDES
figure(5);
subplot(3,1,1);
plot(I,PKamplmeandet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,1);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,PKamplmeanstoAVG(j,:),PKamplmeanstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.2*PKamplmeandet(1,end)]);
xlabel('Control Parameter');ylabel('Peak-to-Trough Amplitude');
title('Subcritical Hopf Bifurcation; Averages = 5');


Iselected = [-0.4 -0.3 -0.2 -0.1 -0.05 -0.01 0.05 0.1 0.2 0.3 0.4]; % CHOOSE
clear Iselect
for k = 1:length(Iselected)
    Iselect(k) = max(findnearest(I,Iselected(k)));
end
tmin=10000;tmax=15000;    % CHOOSE
ymin=-2;ymax=2;
yimin=ymin*2;yimax=ymax*2;
dwnsplquiver=20;quiverstart=2000;quiverend=4000;quiverscale=1;
dwnsplrealimag=10;realimagstart=2000;realimagend=4000;
fmin = 0.05; fmax = 0.5;

% PLOT EXAMPLE TIME TRACES AND PHASE PORTRAITS
% Deterministic
xdr=real(hilbert(Xdet{1}{Iselect(end)}));xdi=imag(hilbert(Xdet{1}{Iselect(end)}));
[bw dens mx1 my1]=kde2d([xdr',xdi']);
for k = 1:length(Iselect);
    figure(6);
    sph=subplot(6,length(Iselect),k);plot(Xdet{1}{Iselect(k)},'k');axis([tmin tmax ymin ymax]);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s %s%s','Deterministic','I = ',num2str(I(Iselect(k)))));
    spp = get(sph, 'pos');
    set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xdet{1}{Iselect(k)}));xdi=imag(hilbert(Xdet{1}{Iselect(k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1 && (max(xdr)-min(xdr)) > 1e-14 && (max(xdi)-min(xdi)) > 1e-14
        sph=subplot(6,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);
        if length(my(:,1)) < size(dens,1)
            my(end:size(dens,1),:)=0;
        elseif length(my(:,1)) > size(dens,1)
            dens(end:length(my(:,1)),:)=0;
        end
        if length(mx(1,:)) < size(dens,2)
            mx(end:size(dens,2))=0;
        elseif length(mx(1,:)) > size(dens,2)
            dens(:,end:length(mx(1,:)))=0;
        end
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(6,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 1.1*max(dens1)]);axis([mx(1,1) mx(1,end) 0 0.001]);
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+5*length(Iselect));plot(fdetfft{j,Iselect(k)},Xdetfft{1,Iselect(k)},'k');axis([fmin fmax 0 1.05*max(Xdetfft{1,Iselect(end)})]);
        %tt=findnearest(fdetfft{j,Iselect(k)},fftfreqdet(1,Iselect(k)));tt=tt(1);hold on;scatter(fdetfft{j,Iselect(k)}(tt),Xdetfft{1,Iselect(k)}(tt),'b.');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end

% Stochastic
xdr=real(hilbert(Xsto{j}{Iselect(end)}));xdi=imag(hilbert(Xsto{j}{Iselect(end)}));
[bw dens mx1 my1]=kde2d([xdr',xdi']);
for j = 1:maxiter
for k = 1:length(Iselect);
    figure(6+j);
    sph=subplot(6,length(Iselect),k);plot(Xsto{j}{Iselect(k)},'r');axis([tmin+500*j tmax+500*j ymin ymax]);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s%s %s%s','Noise = ',num2str(noiselevel(j)),', I = ',num2str(I(Iselect(k)))));
    spp = get(sph, 'pos');
    set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xsto{j}{Iselect(k)}));xdi=imag(hilbert(Xsto{j}{Iselect(k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1
        sph=subplot(5,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);
        if length(my(:,1)) < size(dens,1)
            my(end:size(dens,1),:)=0;
        elseif length(my(:,1)) > size(dens,1)
            dens(end:length(my(:,1)),:)=0;
        end
        if length(mx(1,:)) < size(dens,2)
            mx(end:size(dens,2))=0;
        elseif length(mx(1,:)) > size(dens,2)
            dens(:,end:length(mx(1,:)))=0;
        end
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(6,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 1.1*max(dens1)]);axis([mx(1,1) mx(1,end) 0 0.001]);
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+5*length(Iselect));plot(fstofft{j,Iselect(k)},Xstofft{j,Iselect(k)},'r');axis([fmin fmax 0  1.05*max(Xstofft{j,Iselect(end)})]);
        %tt=findnearest(fstofft{j,Iselect(k)},fftfreqsto(1,Iselect(k)));tt=tt(1);hold on;scatter(fstofft{j,Iselect(k)}(tt),Xstofft{1,Iselect(k)}(tt),'b.');        
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end
end

elseif biftype ==4 
% Fold
close all;
maxiter = 4;
noiselevel = [0.05 0.1 0.2 0.4];

for j = 1:maxiter
    for k = 1:500
        pkspikeratestoAVG(j,k)= mean([pkspikeratesto1(j,k) pkspikeratesto2(j,k) pkspikeratesto3(j,k) pkspikeratesto4(j,k) pkspikeratesto5(j,k)]);
        pkspikeratestoSEM(j,k)= std([pkspikeratesto1(j,k) pkspikeratesto2(j,k) pkspikeratesto3(j,k) pkspikeratesto4(j,k) pkspikeratesto5(j,k)])/sqrt(5);
    end
end


% coefficient of variation 1 = variance/mean
for j = 1:maxiter
    for k = 1:500
        CDstopkAVG(j,k)= mean([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)]);
        CDstopkSEM(j,k)= std([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)])/sqrt(5);
    end
end
%}

% FFT and amplitudes from peak/trough
for j = 1:maxiter
    for k = 1:500
        [~,CDpsto(j,k),~,CDstatsto(j,k)] =  ttest([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)],1);CDtstat(j,k)=CDstatsto(j,k).tstat;
        CDpktimeAVG(j,k)= mean([CDpktime1(j,k) CDpktime2(j,k) CDpktime3(j,k) CDpktime4(j,k) CDpktime5(j,k)]);
        CDpktimeSEM(j,k)= std([CDpktime1(j,k) CDpktime2(j,k) CDpktime3(j,k) CDpktime4(j,k) CDpktime5(j,k)])/sqrt(5);
        CDtrtimeAVG(j,k)= mean([CDtrtime1(j,k) CDtrtime2(j,k) CDtrtime3(j,k) CDtrtime4(j,k) CDtrtime5(j,k)]);
        CDtrtimeSEM(j,k)= std([CDtrtime1(j,k) CDtrtime2(j,k) CDtrtime3(j,k) CDtrtime4(j,k) CDtrtime5(j,k)])/sqrt(5);
        CVpktimeAVG(j,k)= mean([CVpktime1(j,k) CVpktime2(j,k) CVpktime3(j,k) CVpktime4(j,k) CVpktime5(j,k)]);
        CVpktimeSEM(j,k)= std([CVpktime1(j,k) CVpktime2(j,k) CVpktime3(j,k) CVpktime4(j,k) CVpktime5(j,k)])/sqrt(5);
        CVtrtimeAVG(j,k)= mean([CVtrtime1(j,k) CVtrtime2(j,k) CVtrtime3(j,k) CVtrtime4(j,k) CVtrtime5(j,k)]);
        CVtrtimeSEM(j,k)= std([CVtrtime1(j,k) CVtrtime2(j,k) CVtrtime3(j,k) CVtrtime4(j,k) CVtrtime5(j,k)])/sqrt(5);
        fftfreqstoAVG(j,k)= mean([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)]);
        fftfreqstoSEM(j,k)= std([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)])/sqrt(5);
        fftamplstoAVG(j,k)= mean([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)]);
        fftamplstoSEM(j,k)= std([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)])/sqrt(5);
        meanpktimeAVG(j,k)= mean([meanpktime1(j,k) meanpktime2(j,k) meanpktime3(j,k) meanpktime4(j,k) meanpktime5(j,k)]);
        meanpktimeSEM(j,k)= std([meanpktime1(j,k) meanpktime2(j,k) meanpktime3(j,k) meanpktime4(j,k) meanpktime5(j,k)])/sqrt(5);
        meantrtimeAVG(j,k)= mean([meantrtime1(j,k) meantrtime2(j,k) meantrtime3(j,k) meantrtime4(j,k) meantrtime5(j,k)]);
        meantrtimeSEM(j,k)= std([meantrtime1(j,k) meantrtime2(j,k) meantrtime3(j,k) meantrtime4(j,k) meantrtime5(j,k)])/sqrt(5);
        PKamplmeanstoAVG(j,k)= mean([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)]);
        PKamplmeanstoSEM(j,k)= std([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)])/sqrt(5);
    end
    %CDtstat(j,isinf(CDtstat(j,:)))=0;
    klp = findnearest(((CDtstat(j,:))),0);
    if isempty(klp)==0 && sum(klp)~=0
    transitionend(j)=klp(end);clear klp;
    transitionend_p(j) = CDpsto(j,transitionend(j));
    else
        transitionend(j)=500;
        transitionend_p(j)=NaN;
    end
end

% coefficient of dispersion = std/mean

for j = 1:maxiter
    for k = 1:500
        CVstopk1(j,k)=sqrt(IEIvarpksto1(j,k))/IEImeanpksto1(j,k);
        CVstopk2(j,k)=sqrt(IEIvarpksto2(j,k))/IEImeanpksto2(j,k);
        CVstopk3(j,k)=sqrt(IEIvarpksto3(j,k))/IEImeanpksto3(j,k);
        CVstopk4(j,k)=sqrt(IEIvarpksto4(j,k))/IEImeanpksto4(j,k);
        CVstopk5(j,k)=sqrt(IEIvarpksto5(j,k))/IEImeanpksto5(j,k);
        CVstopkAVG(j,k)= mean([CVstopk1(j,k) CVstopk2(j,k) CVstopk3(j,k) CVstopk4(j,k) CVstopk5(j,k)]);
        CVstopkSEM(j,k)= std([CVstopk1(j,k) CVstopk2(j,k) CVstopk3(j,k) CVstopk4(j,k) CVstopk5(j,k)])/sqrt(5);
        IEIvarpkstoAVG(j,k) = mean([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]);
        IEIstdpkstoAVG(j,k) = sqrt(mean([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]));
        IEIvarpkstoSEM(j,k) = std([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)])/sqrt(5);
        IEIstdpkstoSEM(j,k) = sqrt(std([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]))/sqrt(5);
        IEImeanpkstoAVG(j,k) =  mean([IEImeanpksto1(j,k) IEImeanpksto2(j,k) IEImeanpksto3(j,k) IEImeanpksto4(j,k) IEImeanpksto5(j,k)]);
        IEImeanpkstoSEM(j,k) =  std([IEImeanpksto1(j,k) IEImeanpksto2(j,k) IEImeanpksto3(j,k) IEImeanpksto4(j,k) IEImeanpksto5(j,k)])/sqrt(5);
    end
end

% Find the averages for all of the deterministic cases
pkspikeratedet=mean(pkspikeratedet,1);trspikeratedet=mean(trspikeratedet,1);
CDdetpk=mean(CDdetpk,1);CDdettr=mean(CDdettr,1);CDpktimedet=mean(CDpktimedet,1);
CDtrtimedet=mean(CDtrtimedet,1);CVpktimedet=mean(CVpktimedet,1);CVtrtimedet=mean(CVtrtimedet,1);
fftampldet=mean(fftampldet,1);fftfreqdet=mean(fftfreqdet,1);IEImeanpkdet=mean(IEImeanpkdet,1);IEImeantrdet=mean(IEImeantrdet,1);
IEIvarpkdet=mean(IEIvarpkdet,1);IEIvartrdet=mean(IEIvartrdet,1);meanIEIpkdet=mean(meanIEIpkdet,1);
meanIEItrdet=mean(meanIEItrdet,1);meanpktimedet=mean(meanpktimedet,1);meantrtimedet=mean(meantrtimedet,1);
PKamplmeandet=mean(PKamplmeandet,1);pkdiffusiondet=mean(pkdiffusiondet,1);pkIEIspikeratiodet=mean(pkIEIspikeratiodet,1);
pktcorrdet=mean(pktcorrdet,1);trdiffusiondet=mean(trdiffusiondet,1);trtcorrdet=mean(trtcorrdet,1);


tr2=find(pkspikeratedet(1,:)==0);
tr2=max(tr2);
for j = 1:maxiter
    pkspikshift(j) = I(transitionend(j));
    I_shifted(j,:) = I - pkspikshift(j);
end

figure(1);
subplot(2,1,1);
set(0,'DefaultAxesColorOrder',ametrine(4));
plot(mu,pkspikeratedet(1,:),'k');I=mu;
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    %plot(mu,pkspikeratestoAVG(j,:));
    ha=errorbar(mu,pkspikeratestoAVG(j,:),pkspikeratestoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) =  .01;
    Xdata(xright) =  .01;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.5*max(max(pkspikeratestoAVG))]);
xlabel('Control Parameter');ylabel('Spike Rate (spikes/sec)');
title('Fold Bifurcation; Averages = 5');

subplot(2,1,2);
plot(I,pkspikeratedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
     ha=errorbar(I_shifted(j,:),pkspikeratestoAVG(j,:),pkspikeratestoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = .01;
    Xdata(xright) =  .01;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 max(max(pkspikeratestoAVG))]);
xlabel('Control Parameter (shifted)');ylabel('Spike Rate (spikes/sec)');
title('Fold Bifurcation; Averages = 5');

%{
ft=fittype('a*x^b+c');
[fitdet gofdet] = fit(I(I>0)',pkspikeratedet(1,I>0)',ft);
plot(fitdet,'k');
for j = 1:maxiter
    [fitsto{j} gofsto{j}] = fit(I(I>0)',pkspikeratestoAVG(j,I>0)',ft);
    plot(fitsto{j});hold all;;
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4',sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Deterministic','y=a*x^b+c','a=',num2str(fitdet.a),'b=',num2str(fitdet.b),'c=',num2str(fitdet.c),'R^2=',num2str(gofdet.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.05','y=a*x^b+c','a=',num2str(fitsto{1}.a),'b=',num2str(fitsto{1}.b),'c=',num2str(fitsto{1}.c),'R^2=',num2str(gofsto{1}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.1','y=a*x^b+c','a=',num2str(fitsto{2}.a),'b=',num2str(fitsto{2}.b),'c=',num2str(fitsto{2}.c),'R^2=',num2str(gofsto{2}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.2','y=a*x^b+c','a=',num2str(fitsto{3}.a),'b=',num2str(fitsto{3}.b),'c=',num2str(fitsto{3}.c),'R^2=',num2str(gofsto{3}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.4','y=a*x^b+c','a=',num2str(fitsto{4}.a),'b=',num2str(fitsto{4}.b),'c=',num2str(fitsto{4}.c),'R^2=',num2str(gofsto{4}.rsquare)));
axis([mu(1) mu(end) 0 1.5*max(max(pkspikeratestoAVG))]);
xlabel('Control Parameter');ylabel('Spike Rate (spikes/sec)');
title('Fold Bifurcation; Averages = 5');
%}

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(2);
subplot(2,1,1);
plot(I,CDdetpk(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    hold all;
    plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVstopkAVG(j,:),CVstopkSEM(j,:));
    hb = get(ha,'children'); 
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = 0.01;
    Xdata(xright) = 0.01;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(IEImeanpkstoAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 0.1 100000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - peak/peak');
title('Cusp Bifurcation - Between Peaks; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR PEAK/PEAK and TROUGH/TROUGH
figure(2)
subplot(2,1,2);
plot(I,CDdetpk(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDstopkAVG(j,:),CDstopkSEM(j,:));
    hb = get(ha,'children'); 
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = 0.01;
    Xdata(xright) = 0.01;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0.1 100000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - peak/peak');
title('Cusp Bifurcation - Between Peaks; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR EACH PEAK
figure(4)
subplot(2,2,1);
plot(I,CDpktimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(2,2,1);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDpktimeAVG(j,:),CDpktimeSEM(j,:));
    hb = get(ha,'children');
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = 0.01;
    Xdata(xright) = 0.01;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0.1 100000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH PEAK');
title('Cusp Bifurcation - PEAKS; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR EACH TROUGH
figure(4)
subplot(2,2,2);
plot(I,CDtrtimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    subplot(2,2,2);
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDtrtimeAVG(j,:),CDtrtimeSEM(j,:));
    hb = get(ha,'children');
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = 0.01;
    Xdata(xright) = 0.01;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0.1 100000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH TROUGH');
title('Cusp Bifurcation - TROUGHS; Averages = 5');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(4);
subplot(2,2,3);
plot(I,CVpktimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    subplot(2,2,3);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVpktimeAVG(j,:),CVpktimeSEM(j,:));
    hb = get(ha,'children');
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = 0.01;
    Xdata(xright) = 0.01;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(meanpktimeAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 0.1 100000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH PEAK');
title('Cusp Bifurcation - Each Peak; Averages = 5');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(4);
subplot(2,2,4);
plot(I,CVtrtimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    subplot(2,2,4);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVtrtimeAVG(j,:),CVtrtimeSEM(j,:));
    hb = get(ha,'children');  
    set(get(ha,'Parent'),'YScale','log');
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = 0.01;
    Xdata(xright) = 0.01;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(meantrtimeAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)','Noise = 0.4; 1/sqrt(mean)');
axis([mu(1) mu(end) 0.1 100000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH TROUGH');
title('Cusp Bifurcation - Each Trough; Averages = 5');

% PLOT THE AMPLITUDES FROM FFT
figure(5);
subplot(3,1,2);
plot(I,fftampldet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,2);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,fftamplstoAVG(j,:),fftamplstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = 0.01;
    Xdata(xright) = 0.01;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.2*fftamplstoAVG(4,end)]);
xlabel('Control Parameter');ylabel('Amplitude from FFT');
title('Cusp Bifurcation; Averages = 5');

% PLOT THE FREQUENCIES FROM FFT
figure(5);
subplot(3,1,3);
plot(I,fftfreqdet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,3);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,fftfreqstoAVG(j,:),fftfreqstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = 0.01;
    Xdata(xright) = 0.01;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.05*max(max(fftfreqstoAVG(:,:)))]);
xlabel('Control Parameter');ylabel('Frequency from FFT');
title('Cusp Bifurcation; Averages = 5');

% PLOT THE PEAK-TO-TROUGH AMPLITUDES
figure(5);
subplot(3,1,1);
plot(I,PKamplmeandet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,1);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,PKamplmeanstoAVG(j,:),PKamplmeanstoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = 0.01;
    Xdata(xright) = 0.01;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 1.2*PKamplmeanstoAVG(4,findnearest(mu,mu(end)))]);
xlabel('Control Parameter');ylabel('Peak-to-Trough Amplitude');
title('Cusp Bifurcation; Averages = 5');


Iselected = [-0.25 -0.2 -0.15 -0.1 -0.05 -0.04 -0.03 -0.1 -0.01 0]; % CHOOSE
clear Iselect
for k = 1:length(Iselected)
    Iselect(k) = max(findnearest(I,Iselected(k)));
end
tmin=10000;tmax=15000;    % CHOOSE
ymin=-2;ymax=2;
yimin=ymin*2;yimax=ymax*2;
dwnsplquiver=20;quiverstart=2000;quiverend=4000;quiverscale=1;
dwnsplrealimag=10;realimagstart=2000;realimagend=4000;
fmin = 0; fmax = 0.5;

% PLOT EXAMPLE TIME TRACES AND PHASE PORTRAITS
% Deterministic
xdr=real(hilbert(Xdet{1}{Iselect(end)}));xdi=imag(hilbert(Xdet{1}{Iselect(end)}));
%[bw dens mx1 my1]=kde2d([xdr',xdi']);
for k = 1:length(Iselect);
    figure(6);
    sph=subplot(6,length(Iselect),k);plot(Xdet{1}{Iselect(k)},'k');axis([tmin tmax ymin ymax]);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s %s%s','Deterministic','I = ',num2str(I(Iselect(k)))));
    spp = get(sph, 'pos');
    set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xdet{1}{Iselect(k)}));xdi=imag(hilbert(Xdet{1}{Iselect(k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1 && (max(xdr)-min(xdr)) > 1e-14 && (max(xdi)-min(xdi)) > 1e-14
        sph=subplot(6,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);
        if length(my(:,1)) < size(dens,1)
            my(end:size(dens,1),:)=0;
        elseif length(my(:,1)) > size(dens,1)
            dens(end:length(my(:,1)),:)=0;
        end
        if length(mx(1,:)) < size(dens,2)
            mx(end:size(dens,2))=0;
        elseif length(mx(1,:)) > size(dens,2)
            dens(:,end:length(mx(1,:)))=0;
        end
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(6,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 1.1*max(dens1)]);axis([mx(1,1) mx(1,end) 0 0.0005]);
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+5*length(Iselect));plot(fdetfft{j,Iselect(k)},Xdetfft{1,Iselect(k)},'k');axis([fmin fmax 0 1.1*max(Xdetfft{1,Iselect(end)})]);
        %tt=findnearest(fdetfft{j,Iselect(k)},fftfreqdet(1,Iselect(k)));tt=tt(1);hold on;scatter(fdetfft{j,Iselect(k)}(tt),Xdetfft{1,Iselect(k)}(tt),'b.');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end

% Stochastic
xdr=real(hilbert(Xsto{j}{Iselect(end)}));xdi=imag(hilbert(Xsto{j}{Iselect(end)}));
[bw dens mx1 my1]=kde2d([xdr',xdi']);
for j = 1:maxiter
for k = 1:length(Iselect);
    figure(6+j);
    sph=subplot(6,length(Iselect),k);plot(Xsto{j}{Iselect(k)},'r');axis([tmin+500*j tmax+500*j ymin ymax]);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s%s %s%s','Noise = ',num2str(noiselevel(j)),', I = ',num2str(I(Iselect(k)))));
    spp = get(sph, 'pos');
    set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xsto{j}{Iselect(k)}));xdi=imag(hilbert(Xsto{j}{Iselect(k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1
        sph=subplot(5,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);
        if length(my(:,1)) < size(dens,1)
            my(end:size(dens,1),:)=0;
        elseif length(my(:,1)) > size(dens,1)
            dens(end:length(my(:,1)),:)=0;
        end
        if length(mx(1,:)) < size(dens,2)
            mx(end:size(dens,2))=0;
        elseif length(mx(1,:)) > size(dens,2)
            dens(:,end:length(mx(1,:)))=0;
        end
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(6,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 1.1*max(dens1)]);axis([mx(1,1) mx(1,end) 0 0.0004]);
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+5*length(Iselect));plot(fstofft{j,Iselect(k)},Xstofft{j,Iselect(k)},'r');axis([fmin fmax 0 1.1*max(Xstofft{j,Iselect(end)})]);
        %tt=findnearest(fstofft{j,Iselect(k)},fftfreqsto(1,Iselect(k)));tt=tt(1);hold on;scatter(fstofft{j,Iselect(k)}(tt),Xstofft{1,Iselect(k)}(tt),'b.');        
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end
end

elseif biftype ==5
    if hbcalc==1
% HAIR BUNDLE MODEL - GENERATE AVERAGES

index = 3;          % Which force or stiffness index would you like to use?
forcestiff = 1;     % stiffness or force scan? (1=force scan; 2=stiffness scan);
Xsto1=Xsto; Xdet1=Xdet; pkspikeratedet1=pkspikeratedet; IEIvarpkdet1=IEIvarpkdet;IEImeanpkdet1=IEImeanpkdet;CDdetpk1=CDdetpk;
clear Xsto Xdet pkspikeratedet IEImeanpkdet IEIvarpkdet CDdetpk
noiselevel = [0.05 0.1 0.2];
maxiter = length(noiselevel);

if forcestiff==1
    mu = F;
else
    mu = ke;
end

for j = 1:500       % Isolate the appropriate index
    for m = 1:length(noiselevel)
        Xsto{m}{j} = Xsto1{j}{index}{m};
        Xdet{m}{j} = Xdet1{j}{index}{m};
        pkspikeratedet(m,j) = pkspikeratedet1(j,index,m);
        pkspikeratesto1(m,j) = pkspikeratesto_1(j,index,m);
        pkspikeratesto2(m,j) = pkspikeratesto_2(j,index,m);
        pkspikeratesto3(m,j) = pkspikeratesto_3(j,index,m);
        pkspikeratesto4(m,j) = pkspikeratesto_4(j,index,m);
        pkspikeratesto5(m,j) = pkspikeratesto_5(j,index,m);
        CDdetpk(m,j) = CDdetpk1(j,index,m);
        CDstopk1(m,j) = CDstopk_1(j,index,m);
        CDstopk2(m,j) = CDstopk_2(j,index,m);
        CDstopk3(m,j) = CDstopk_3(j,index,m);
        CDstopk4(m,j) = CDstopk_4(j,index,m);
        CDstopk5(m,j) = CDstopk_5(j,index,m);
        IEIvarpkdet(m,j) = IEIvarpkdet1(j,index,m);
        IEIvarpksto1(m,j) = IEIvarpksto_1(j,index,m);
        IEIvarpksto2(m,j) = IEIvarpksto_2(j,index,m);
        IEIvarpksto3(m,j) = IEIvarpksto_3(j,index,m);
        IEIvarpksto4(m,j) = IEIvarpksto_4(j,index,m);
        IEIvarpksto5(m,j) = IEIvarpksto_5(j,index,m);
        IEImeanpkdet(m,j) = IEImeanpkdet1(j,index,m);
        IEImeanpksto1(m,j) = IEImeanpksto_1(j,index,m);
        IEImeanpksto2(m,j) = IEImeanpksto_2(j,index,m);
        IEImeanpksto3(m,j) = IEImeanpksto_3(j,index,m);
        IEImeanpksto4(m,j) = IEImeanpksto_4(j,index,m);
        IEImeanpksto5(m,j) = IEImeanpksto_5(j,index,m);
    end
end

if forcestiff == 1
    mu = F;
elseif forcestiff == 2
    mu = ke;
end
I=mu;

for j = 1:maxiter
    for k = 1:500
        pkspikeratestoAVG(j,k)= mean([pkspikeratesto1(j,k) pkspikeratesto2(j,k) pkspikeratesto3(j,k) pkspikeratesto4(j,k) pkspikeratesto5(j,k)]);
        pkspikeratestoSEM(j,k)= std([pkspikeratesto1(j,k) pkspikeratesto2(j,k) pkspikeratesto3(j,k) pkspikeratesto4(j,k) pkspikeratesto5(j,k)])/sqrt(5);
    end
end


% FFT and amplitudes from peak/trough
for j = 1:maxiter
    for k = 1:500
        [~,CDpsto(j,k),~,CDstatsto(j,k)] =  ttest([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)],1);CDtstat(j,k)=CDstatsto(j,k).tstat;
        CDpktimeAVG(j,k)= mean([CDpktime1(j,k) CDpktime2(j,k) CDpktime3(j,k) CDpktime4(j,k) CDpktime5(j,k)]);
        CDpktimeSEM(j,k)= std([CDpktime1(j,k) CDpktime2(j,k) CDpktime3(j,k) CDpktime4(j,k) CDpktime5(j,k)])/sqrt(5);
        CDtrtimeAVG(j,k)= mean([CDtrtime1(j,k) CDtrtime2(j,k) CDtrtime3(j,k) CDtrtime4(j,k) CDtrtime5(j,k)]);
        CDtrtimeSEM(j,k)= std([CDtrtime1(j,k) CDtrtime2(j,k) CDtrtime3(j,k) CDtrtime4(j,k) CDtrtime5(j,k)])/sqrt(5);
        CVpktimeAVG(j,k)= mean([CVpktime1(j,k) CVpktime2(j,k) CVpktime3(j,k) CVpktime4(j,k) CVpktime5(j,k)]);
        CVpktimeSEM(j,k)= std([CVpktime1(j,k) CVpktime2(j,k) CVpktime3(j,k) CVpktime4(j,k) CVpktime5(j,k)])/sqrt(5);
        CVtrtimeAVG(j,k)= mean([CVtrtime1(j,k) CVtrtime2(j,k) CVtrtime3(j,k) CVtrtime4(j,k) CVtrtime5(j,k)]);
        CVtrtimeSEM(j,k)= std([CVtrtime1(j,k) CVtrtime2(j,k) CVtrtime3(j,k) CVtrtime4(j,k) CVtrtime5(j,k)])/sqrt(5);
        fftfreqstoAVG(j,k)= mean([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)]);
        fftfreqstoSEM(j,k)= std([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)])/sqrt(5);
        fftamplstoAVG(j,k)= mean([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)]);
        fftamplstoSEM(j,k)= std([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)])/sqrt(5);
        meanpktimeAVG(j,k)= mean([meanpktime1(j,k) meanpktime2(j,k) meanpktime3(j,k) meanpktime4(j,k) meanpktime5(j,k)]);
        meanpktimeSEM(j,k)= std([meanpktime1(j,k) meanpktime2(j,k) meanpktime3(j,k) meanpktime4(j,k) meanpktime5(j,k)])/sqrt(5);
        meantrtimeAVG(j,k)= mean([meantrtime1(j,k) meantrtime2(j,k) meantrtime3(j,k) meantrtime4(j,k) meantrtime5(j,k)]);
        meantrtimeSEM(j,k)= std([meantrtime1(j,k) meantrtime2(j,k) meantrtime3(j,k) meantrtime4(j,k) meantrtime5(j,k)])/sqrt(5);
        PKamplmeanstoAVG(j,k)= mean([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)]);
        PKamplmeanstoSEM(j,k)= std([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)])/sqrt(5);
    end
    %CDtstat(j,isinf(CDtstat(j,:)))=0;
    klp = findnearest(((CDtstat(j,:))),0);     transitionend(j)=klp(end);clear klp;
    transitionend_p(j) = CDpsto(j,transitionend(j));
end

% coefficient of variation 1 = variance/mean
for j = 1:maxiter
    for k = 1:500
        CDstopkAVG(j,k)= mean([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)]);
        CDstopkSEM(j,k)= std([CDstopk1(j,k) CDstopk2(j,k) CDstopk3(j,k) CDstopk4(j,k) CDstopk5(j,k)])/sqrt(5);
    end
end
%}

% coefficient of dispersion = std/mean

for j = 1:maxiter
    for k = 1:500
        CVstopk1(j,k)=sqrt(IEIvarpksto1(j,k))/IEImeanpksto1(j,k);
        CVstopk2(j,k)=sqrt(IEIvarpksto2(j,k))/IEImeanpksto2(j,k);
        CVstopk3(j,k)=sqrt(IEIvarpksto3(j,k))/IEImeanpksto3(j,k);
        CVstopk4(j,k)=sqrt(IEIvarpksto4(j,k))/IEImeanpksto4(j,k);
        CVstopk5(j,k)=sqrt(IEIvarpksto5(j,k))/IEImeanpksto5(j,k);
        CVstopkAVG(j,k)= mean([CVstopk1(j,k) CVstopk2(j,k) CVstopk3(j,k) CVstopk4(j,k) CVstopk5(j,k)]);
        CVstopkSEM(j,k)= std([CVstopk1(j,k) CVstopk2(j,k) CVstopk3(j,k) CVstopk4(j,k) CVstopk5(j,k)])/sqrt(5);
        IEIvarpkstoAVG(j,k) = mean([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]);
        IEIstdpkstoAVG(j,k) = sqrt(mean([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]));
        IEIvarpkstoSEM(j,k) = std([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)])/sqrt(5);
        IEIstdpkstoSEM(j,k) = sqrt(std([IEIvarpksto1(j,k) IEIvarpksto2(j,k) IEIvarpksto3(j,k) IEIvarpksto4(j,k) IEIvarpksto5(j,k)]))/sqrt(5);
        IEImeanpkstoAVG(j,k) =  mean([IEImeanpksto1(j,k) IEImeanpksto2(j,k) IEImeanpksto3(j,k) IEImeanpksto4(j,k) IEImeanpksto5(j,k)]);
        IEImeanpkstoSEM(j,k) =  std([IEImeanpksto1(j,k) IEImeanpksto2(j,k) IEImeanpksto3(j,k) IEImeanpksto4(j,k) IEImeanpksto5(j,k)])/sqrt(5);
    end
end

% Find the averages for all of the deterministic cases
pkspikeratedet=mean(pkspikeratedet,1);trspikeratedet=mean(trspikeratedet,1);
CDdetpk=mean(CDdetpk,1);CDdettr=mean(CDdettr,1);CDpktimedet=mean(CDpktimedet,1);
CDtrtimedet=mean(CDtrtimedet,1);CVpktimedet=mean(CVpktimedet,1);CVtrtimedet=mean(CVtrtimedet,1);
fftampldet=mean(fftampldet,1);fftfreqdet=mean(fftfreqdet,1);IEImeanpkdet=mean(IEImeanpkdet,1);IEImeantrdet=mean(IEImeantrdet,1);
IEIvarpkdet=mean(IEIvarpkdet,1);IEIvartrdet=mean(IEIvartrdet,1);meanIEIpkdet=mean(meanIEIpkdet,1);
meanIEItrdet=mean(meanIEItrdet,1);meanpktimedet=mean(meanpktimedet,1);meantrtimedet=mean(meantrtimedet,1);
PKamplmeandet=mean(PKamplmeandet,1);pkdiffusiondet=mean(pkdiffusiondet,1);pkIEIspikeratiodet=mean(pkIEIspikeratiodet,1);
pktcorrdet=mean(pktcorrdet,1);trdiffusiondet=mean(trdiffusiondet,1);trtcorrdet=mean(trtcorrdet,1);

tr2=find(pkspikeratedet(1,:)==0);
tr2=min(tr2);
for j = 1:maxiter
    pkspikshift(j) = I(transitionend(j));
    I_shifted(j,:) = I - pkspikshift(j);
end

    else
% HB MODEL PLOTS

close all;
mustart=350;muend=500;

figure(1);
subplot(2,1,1);
plot(mu,pkspikeratedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    %plot(mu,pkspikeratestoAVG(j,:));
    ha=errorbar(mu,pkspikeratestoAVG(j,:),pkspikeratestoSEM(j,:));
    set(gca,'xdir','reverse')
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft) + .02;
    Xdata(xright) = Xdata(xright) - .02;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2');
axis([mu(mustart) mu(muend) 0 1.1*max(max(pkspikeratestoAVG))]);
if forcestiff==1
    xlabel('Force');ylabel('Spike Rate (spikes/sec)');
elseif forcestiff==2
    xlabel('Stiffness');ylabel('Spike Rate (spikes/sec)');
end
title('Hair Bundle Model; Averages = 5');

subplot(2,1,2);
plot(I,pkspikeratedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
     ha=errorbar(I_shifted(j,:),pkspikeratestoAVG(j,:),pkspikeratestoSEM(j,:));
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft)+0.1;
    Xdata(xright) = Xdata(xright)-0.1;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4');
axis([mu(1) mu(end) 0 max(max(pkspikeratestoAVG))]);
xlabel('Control Parameter (shifted)');ylabel('Spike Rate (spikes/sec)');
title('Hair Bundle Model; Averages = 5');


%{
ft=fittype('a*x^b+c');
[fitdet gofdet] = fit(I(2:g)',fliplr(pkspikeratedet(1,2:g))',ft);
plot(fitdet,'k');
for j = 1:maxiter
    [fitsto{j} gofsto{j}] = fit(I(2:g)',fliplr(pkspikeratestoAVG(j,2:g))',ft);
    plot(fitsto{j});hold all;
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2',sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Deterministic','y=a*x^b+c','a=',num2str(fitdet.a),'b=',num2str(fitdet.b),'c=',num2str(fitdet.c),'R^2=',num2str(gofdet.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.05','y=a*x^b+c','a=',num2str(fitsto{1}.a),'b=',num2str(fitsto{1}.b),'c=',num2str(fitsto{1}.c),'R^2=',num2str(gofsto{1}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.1','y=a*x^b+c','a=',num2str(fitsto{2}.a),'b=',num2str(fitsto{2}.b),'c=',num2str(fitsto{2}.c),'R^2=',num2str(gofsto{2}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.2','y=a*x^b+c','a=',num2str(fitsto{3}.a),'b=',num2str(fitsto{3}.b),'c=',num2str(fitsto{3}.c),'R^2=',num2str(gofsto{3}.rsquare)));
axis([mu(2) mu(g) 0 1.1*max(max(pkspikeratestoAVG))]);
if forcestiff==1
    xlabel('Force');ylabel('Spike Rate (spikes/sec)');
elseif forcestiff==2
    xlabel('Stiffness');ylabel('Spike Rate (spikes/sec)');
end
title('Hair Bundle Model; Averages = 5 --- REVERSED FOR FITTING');
%}

figure(2);
plot(mu,CDdetpk(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(mu,CVstopkAVG(j,:),CVstopkSEM(j,:));set(gca,'xdir','reverse')
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft) + .02;
    Xdata(xright) = Xdata(xright) - .02;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(mu,1./sqrt(IEImeanpkstoAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(mu,ones(1,length(mu)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)');
axis([mu(mustart) mu(muend) 0 1.5]);
xlabel('Control Parameter');ylabel('Coefficient of Variation');
title('Hair Bundle Model; Averages = 5');

figure(3)
plot(mu,CDdetpk(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(mu,CDstopkAVG(j,:),CDstopkSEM(j,:));set(gca,'xdir','reverse')
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft) + .02;
    Xdata(xright) = Xdata(xright) - .02;
    set(hb(2),'Xdata',Xdata)
end
plot(mu,ones(1,length(mu)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2');
axis([mu(mustart) mu(muend) 0 1.5]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion');
title('Hair Bundle Model; Averages = 5');


% PLOT THE COEFFICIENT OF DISPERSION FOR EACH PEAK
figure(4)
subplot(2,2,1);
plot(I,CDpktimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(2,2,1);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDpktimeAVG(j,:),CDpktimeSEM(j,:));set(gca,'xdir','reverse')
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft) + .02;
    Xdata(xright) = Xdata(xright) - .02;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2');
axis([mu(mustart) mu(muend) 0 1.5]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH PEAK');
title('Hair Bundle Model - PEAKS; Averages = 5');

% PLOT THE COEFFICIENT OF DISPERSION FOR EACH TROUGH
figure(4)
subplot(2,2,2);
plot(I,CDtrtimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    hold all;
    subplot(2,2,2);
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CDtrtimeAVG(j,:),CDtrtimeSEM(j,:));set(gca,'xdir','reverse')
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft) + .02;
    Xdata(xright) = Xdata(xright) - .02;
    set(hb(2),'Xdata',Xdata)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2');
axis([mu(mustart) mu(muend) 0 1.5]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH TROUGH');
title('Hair Bundle Model - TROUGHS; Averages = 5');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(4);
subplot(2,2,3);
plot(I,CVpktimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    subplot(2,2,3);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVpktimeAVG(j,:),CVpktimeSEM(j,:));set(gca,'xdir','reverse')
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft) + .02;
    Xdata(xright) = Xdata(xright) - .02;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(meanpktimeAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)');
axis([mu(mustart) mu(muend) 0 1.5]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH PEAK');
title('Hair Bundle Model - Each Peak; Averages = 5');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
figure(4);
subplot(2,2,4);
plot(I,CVtrtimedet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
colors2=ametrine(4);
for j = 1:maxiter
    subplot(2,2,4);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,CVtrtimeAVG(j,:),CVtrtimeSEM(j,:));set(gca,'xdir','reverse')
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft) + .02;
    Xdata(xright) = Xdata(xright) - .02;
    set(hb(2),'Xdata',Xdata)
end
for j = 1:maxiter
        plot(I,1./sqrt(meantrtimeAVG(j,:)),'Color',colors2(j,:),'LineStyle','--','LineWidth',0.3)
end
plot(I,ones(1,length(I)),'g--');
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.05; 1/sqrt(mean)','Noise = 0.1; 1/sqrt(mean)','Noise = 0.2; 1/sqrt(mean)');
axis([mu(mustart) mu(muend) 0 1.5]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH TROUGH');
title('Hair Bundle Model - Each Trough; Averages = 5');

% PLOT THE AMPLITUDES FROM FFT
figure(5);
subplot(3,1,2);
plot(I,smooth(fftampldet(1,:),5),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,2);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,fftamplstoAVG(j,:),fftamplstoSEM(j,:));set(gca,'xdir','reverse')
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft) + .02;
    Xdata(xright) = Xdata(xright) - .02;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2');
axis([mu(mustart) mu(muend) 0 1.1*fftampldet(1,findnearest(mu,0))]);
xlabel('Control Parameter');ylabel('Amplitude from FFT');
title('Hair Bundle Model; Averages = 5');

% PLOT THE FREQUENCIES FROM FFT
figure(5);
subplot(3,1,3);
plot(I,fftfreqdet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,3);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,fftfreqstoAVG(j,:),fftfreqstoSEM(j,:));set(gca,'xdir','reverse')
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft) + .02;
    Xdata(xright) = Xdata(xright) - .02;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2');
axis([mu(mustart) mu(muend) 0 1.1*fftfreqsto(1,findnearest(mu,0))]);
xlabel('Control Parameter');ylabel('Frequency from FFT');
title('Hair Bundle Model; Averages = 5');

% PLOT THE PEAK-TO-TROUGH AMPLITUDES
figure(5);
subplot(3,1,1);
plot(I,PKamplmeandet(1,:),'k');
set(0,'DefaultAxesColorOrder',ametrine(4));
for j = 1:maxiter
    subplot(3,1,1);
    hold all;
    %plot(I,pkspikeratestoAVG(j,:));
    ha=errorbar(I,PKamplmeanstoAVG(j,:),PKamplmeanstoSEM(j,:));set(gca,'xdir','reverse')
    hb = get(ha,'children');  
    Xdata = get(hb(2),'Xdata');
    temp = 4:3:length(Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1; 
    % Change line length
    Xdata(xleft) = Xdata(xleft) + .02;
    Xdata(xright) = Xdata(xright) - .02;
    set(hb(2),'Xdata',Xdata)
end
%legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2');
axis([mu(mustart) mu(muend) 0 1.1*PKamplmeansto(3,findnearest(mu,0))]);
xlabel('Control Parameter');ylabel('Peak-to-Trough Amplitude');
title('Hair Bundle Model; Averages = 5');

% PLOT EXAMPLE TIME TRACES AND PHASE PORTRAITS % Deterministics for k = 1:length(Iselect);     figure(6);     sph=subplot(5,length(Iselect),k);plot(Xdet{1}{Iselect(k)},'k');axis([tmin tmax ymin ymax]);     %xlabel('Time');ylabel('Position');     title(sprintf('%s %s%s','Deterministic','I = ',num2str(I(Iselect(k)))));     spp = get(sph, 'pos');     set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);     xdr=real(hilbert(Xdet{1}{Iselect(k)}));xdi=imag(hilbert(Xdet{1}{Iselect(k)}));     if length(unique(xdr)) > 1 && length(unique(xdi)) > 1         sph=subplot(5,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi',2^8,[ymin,ymin],[ymax,ymax]]]);         if length(my(:,1)) < size(dens,1)             my(end:size(dens,1),:)=0;         elseif length(my(:,1)) > size(dens,1)             dens(end:length(my(:,1)),:)=0;         end         if length(mx(1,:)) < size(dens,2)             mx(end:size(dens,2))=0;         elseif length(mx(1,:)) > size(dens,2)             dens(:,end:length(mx(1,:)))=0;         end         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         sph=subplot(5,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([mx(1,1) mx(1,end) my(1,1) my(end,1)]);         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);         sph=subplot(5,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([mx(1,1) mx(1,end) my(1,1) my(end,1)]);         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);         sph=subplot(5,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);         [bw1,dens1,xmesh]=kde1d(xdr);         sph=subplot(5,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'r');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);     else         sph=subplot(5,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([yimin yimax yimin yimax]);         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);         sph=subplot(5,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([yimin yimax yimin yimax]);         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);     end end % Stochastic for j = 1:maxiter for k = 1:length(Iselect);     figure(6+j);     sph=subplot(5,length(Iselect),k);plot(Xsto{j}{Iselect(k)},'r');axis([tmin+500*j tmax+500*j ymin ymax]);     %xlabel('Time');ylabel('Position');     title(sprintf('%s%s %s%s','Noise = ',num2str(noiselevel(j)),', I = ',num2str(I(Iselect(k)))));     spp = get(sph, 'pos');     set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);     xdr=real(hilbert(Xsto{j}{Iselect(k)}));xdi=imag(hilbert(Xsto{j}{Iselect(k)}));     if length(unique(xdr)) > 1 && length(unique(xdi)) > 1         sph=subplot(5,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi',2^8,[ymin,ymin],[ymax,ymax]]]);         if length(my(:,1)) < size(dens,1)             my(end:size(dens,1),:)=0;         elseif length(my(:,1)) > size(dens,1)             dens(end:length(my(:,1)),:)=0;         end         if length(mx(1,:)) < size(dens,2)             mx(end:size(dens,2))=0;         elseif length(mx(1,:)) > size(dens,2)             dens(:,end:length(mx(1,:)))=0;         end         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         sph=subplot(5,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([mx(1,1) mx(1,end) my(1,1) my(end,1)]);         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);         sph=subplot(5,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([mx(1,1) mx(1,end) my(1,1) my(end,1)]);         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);         sph=subplot(5,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);         [bw1,dens1,xmesh]=kde1d(xdr);         sph=subplot(5,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'r');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);     else         sph=subplot(5,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([yimin yimax yimin yimax]);         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);         sph=subplot(5,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([yimin yimax yimin yimax]);         %xlabel('Real');ylabel('Imaginary');         spp = get(sph, 'pos');         set(sph, 'Position', [spp(1) 0.95*spp(2) 1.3*spp(3) 1.4*spp(4)]);     end end end
Iselected = [3.5 3.4 3.3 3.2 3.1 3.0 2.9 2.8 2.7 2.6]; % CHOOSE
Iselected = fliplr(Iselected);

clear Iselect
for k = 1:length(Iselected)
    Iselect(k) = max(findnearest(mu,Iselected(k)));
end
tmin=3000;tmax=5000;    % CHOOSE
ymin=-4;ymax=4;
yimin=ymin*2;yimax=ymax*2;
dwnsplquiver=1;quiverstart=3000;quiverend=4000;quiverscale=2;
dwnsplrealimag=1;realimagstart=2000;realimagend=4000;
fmin = 0; fmax = 0.1;

% PLOT EXAMPLE TIME TRACES AND PHASE PORTRAITS
% Deterministic
xdr=real(hilbert(Xdet{1}{Iselect(end)}));xdi=imag(hilbert(Xdet{1}{Iselect(end)}));
%[bw dens mx1 my1]=kde2d([xdr',xdi']);
for k = 1:length(Iselect);
    figure(6);
    sph=subplot(6,length(Iselect),k);plot(Xdet{1}{Iselect(k)},'k');axis([tmin tmax ymin ymax]);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s %s%s','Deterministic','I = ',num2str(I(Iselect(k)))));
    spp = get(sph, 'pos');
    set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xdet{1}{Iselect(k)}));xdi=imag(hilbert(Xdet{1}{Iselect(k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1 && (max(xdr)-min(xdr)) > 1e-14 && (max(xdi)-min(xdi)) > 1e-14
        sph=subplot(6,length(Iselect),k+3*length(Iselect));[bw dens mx1 my1]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);[bw dens mx my]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);
        if length(my(:,1)) < size(dens,1)
            my(end:size(dens,1),:)=0;
        elseif length(my(:,1)) > size(dens,1)
            dens(end:length(my(:,1)),:)=0;
        end
        if length(mx(1,:)) < size(dens,2)
            mx(end:size(dens,2))=0;
        elseif length(mx(1,:)) > size(dens,2)
            dens(:,end:length(mx(1,:)))=0;
        end
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(6,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 1.1*max(dens1)]);axis([mx(1,1) mx(1,end) 0 0.001]);
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+5*length(Iselect));plot(fdetfft{j,Iselect(k)},Xdetfft{1,Iselect(k)},'k');axis([fmin fmax 0 1.1*max(Xdetfft{1,Iselect(end)})]);
        %tt=findnearest(fdetfft{j,Iselect(k)},fftfreqdet(1,Iselect(k)));tt=tt(1);hold on;scatter(fdetfft{j,Iselect(k)}(tt),Xdetfft{1,Iselect(k)}(tt),'b.');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end

% Stochastic
xdr=real(hilbert(Xsto{j}{Iselect(end)}));xdi=imag(hilbert(Xsto{j}{Iselect(end)}));
[bw dens mx1 my1]=kde2d([xdr',xdi']);
for j = 1:maxiter
for k = 1:length(Iselect);
    figure(6+j);
    sph=subplot(6,length(Iselect),k);plot(Xsto{j}{Iselect(k)},'r');axis([tmin+500*j tmax+500*j ymin ymax]);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s%s %s%s','Noise = ',num2str(noiselevel(j)),', I = ',num2str(I(Iselect(k)))));
    spp = get(sph, 'pos');
    set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xsto{j}{Iselect(k)}));xdi=imag(hilbert(Xsto{j}{Iselect(k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1
        sph=subplot(5,length(Iselect),k+3*length(Iselect));[bw dens mx my]=kde2d([xdr',xdi'],2^8,[ymin,ymin],[ymax,ymax]);
        if length(my(:,1)) < size(dens,1)
            my(end:size(dens,1),:)=0;
        elseif length(my(:,1)) > size(dens,1)
            dens(end:length(my(:,1)),:)=0;
        end
        if length(mx(1,:)) < size(dens,2)
            mx(end:size(dens,2))=0;
        elseif length(mx(1,:)) > size(dens,2)
            dens(:,end:length(mx(1,:)))=0;
        end
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([mx1(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+3*length(Iselect));pcolor(mx(1,:)',my(:,1)',dens./sum(sum(dens)));shading interp;load jetnew.mat;colormap(cjetnew);caxis([-1e-5 0.001]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(6,length(Iselect),k+4*length(Iselect));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 1.1*max(dens1)]);axis([mx(1,1) mx(1,end) 0 0.001]);
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+5*length(Iselect));plot(fstofft{j,Iselect(k)},Xstofft{j,Iselect(k)},'r');axis([fmin fmax 0 1.1*max(Xstofft{j,Iselect(k)})]);
        %tt=findnearest(fstofft{j,Iselect(k)},fftfreqsto(1,Iselect(k)));tt=tt(1);hold on;scatter(fstofft{j,Iselect(k)}(tt),Xstofft{1,Iselect(k)}(tt),'b.');        
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(6,length(Iselect),k+length(Iselect));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'r');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(6,length(Iselect),k+2*length(Iselect));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'r');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end
end
    end
end


