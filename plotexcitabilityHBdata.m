function plotexcitabilityHBdata(filename,mu,saveyn)
% This function plots data for hair-bundle excitability experiments. The
% input file should be one generated from excitabilityanalysisdata().
%
% plotexcitabilityHBdata(filename,mu,saveyn)
%
% filename: file name from excitabilityanalysisdata()
% mu: array of control parameter values
% saveyn: save the figures? (1=yes)
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%


load(filename)
close all;


% Reshape all vectors
pkspikeratedet = reshape(pkspikeratedet,1,a*(logdata.data(1,8)));
CDdetpk = reshape(CDdetpk,1,a*(logdata.data(1,8)));
fftampldet = reshape(fftampldet,1,a*(logdata.data(1,8)));
detmean = reshape(detmean,1,a*(logdata.data(1,8)));
detvar = reshape(detvar,1,a*(logdata.data(1,8)));
trspikeratedet = reshape(trspikeratedet,1,a*(logdata.data(1,8)));
CDdettime = reshape(CDdettime,1,a*(logdata.data(1,8)));
CDdettr = reshape(CDdettr,1,a*(logdata.data(1,8)));
IEIvarpkdet = reshape(IEIvarpkdet,1,a*(logdata.data(1,8)));
IEIvartrdet = reshape(IEIvartrdet,1,a*(logdata.data(1,8)));
IEImeanpkdet = reshape(IEImeanpkdet,1,a*(logdata.data(1,8)));
IEImeantrdet = reshape(IEImeantrdet,1,a*(logdata.data(1,8)));
trtcorrdet = reshape(trtcorrdet,1,a*(logdata.data(1,8)));
pktcorrdet = reshape(pktcorrdet,1,a*(logdata.data(1,8)));
meanIEIpkdet = reshape(meanIEIpkdet,1,a*(logdata.data(1,8)));
meanIEItrdet = reshape(meanIEItrdet,1,a*(logdata.data(1,8)));
pkIEIspikeratiodet = reshape(pkIEIspikeratiodet,1,a*(logdata.data(1,8)));
trIEIspikeratiodet = reshape(trIEIspikeratiodet,1,a*(logdata.data(1,8)));
meantrtimedet = reshape(meantrtimedet,1,a*(logdata.data(1,8)));
meanpktimedet = reshape(meanpktimedet,1,a*(logdata.data(1,8)));
CDtrtimedet = reshape(CDtrtimedet,1,a*(logdata.data(1,8)));
CDpktimedet = reshape(CDpktimedet,1,a*(logdata.data(1,8)));
CVtrtimedet = reshape(CVtrtimedet,1,a*(logdata.data(1,8)));
CVpktimedet = reshape(CVpktimedet,1,a*(logdata.data(1,8)));
PKamplmeandet = reshape(PKamplmeandet,1,a*(logdata.data(1,8)));
PKamplsemdet = reshape(PKamplsemdet,1,a*(logdata.data(1,8)));
fftampldet = reshape(fftampldet,1,a*(logdata.data(1,8)));
fftfreqdet = reshape(fftfreqdet,1,a*(logdata.data(1,8)));


% PLOT THE SPIKE RATES
I=mu;
hrt(1)=figure(1);
plot(mu,pkspikeratedet(1,:),'k');
set(0,'DefaultAxesColorOrder',jet(4));
    %plot(I,pkspikeratestoAVG(j,:));
    ha=plot(I,pkspikeratedet);
    set(gca,'xdir','reverse')
axis([mu(1) mu(end) 0 max(max(pkspikeratedet))]);
xlabel('Control Parameter');ylabel('Spike Rate (spikes/sec)');
title('Hair Bundle Data');


%{
ft=fittype('a*x^b+c');
[fitdet gofdet] = fit(I(find(I>0))',(1,find(I>0))',ft);
plot(fitdet,'k');
for j = 1:maxiter
    [fitsto{j} gofsto{j}] = fit(I(I>0)',pkspikeratestoAVG(j,I>0)',ft);
    plot(fitsto{j});hold all;;
end
legend('Deterministic','Noise = 0.05','Noise = 0.1','Noise = 0.2','Noise = 0.4',sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Deterministic','y=a*x^b+c','a=',num2str(fitdet.a),'b=',num2str(fitdet.b),'c=',num2str(fitdet.c),'R^2=',num2str(gofdet.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.05','y=a*x^b+c','a=',num2str(fitsto{1}.a),'b=',num2str(fitsto{1}.b),'c=',num2str(fitsto{1}.c),'R^2=',num2str(gofsto{1}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.1','y=a*x^b+c','a=',num2str(fitsto{2}.a),'b=',num2str(fitsto{2}.b),'c=',num2str(fitsto{2}.c),'R^2=',num2str(gofsto{2}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.2','y=a*x^b+c','a=',num2str(fitsto{3}.a),'b=',num2str(fitsto{3}.b),'c=',num2str(fitsto{3}.c),'R^2=',num2str(gofsto{3}.rsquare)),sprintf('%s\n%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s','Noise=0.4','y=a*x^b+c','a=',num2str(fitsto{4}.a),'b=',num2str(fitsto{4}.b),'c=',num2str(fitsto{4}.c),'R^2=',num2str(gofsto{4}.rsquare)));
axis([mu(1) mu(end) 0 max(max(pkspikeratestoAVG))]);
xlabel('Control Parameter');ylabel('Spike Rate (spikes/sec)');
%}

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
hrt(2)=figure(2);
subplot(2,1,1);
ha=plot(I,CDdetpk(1,:),'k');
hb = get(ha,'children'); 
set(get(ha,'Parent'),'YScale','log');
set(gca,'xdir','reverse');hold on;
plot(I,1./sqrt(IEImeanpkdet),'k--')
plot(I,ones(1,length(I)),'g--');
legend('Data','1/sqrt(mean','CV=1');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - peak/peak');
title('Bundle Data - Between Peaks');

% PLOT THE COEFFICIENT OF DISPERSION FOR PEAK/PEAK and TROUGH/TROUGH
hrt(2)=figure(2);
subplot(2,1,2);
ha=plot(I,CDdettr(1,:),'k');
hb = get(ha,'children'); 
set(get(ha,'Parent'),'YScale','log');
set(gca,'xdir','reverse'); hold on;
plot(I,ones(1,length(I)),'g--');
legend('Data','CD=1');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - trough/trough');
title('Bundle Data - Between Peaks');

% PLOT THE COEFFICIENT OF DISPERSION FOR EACH PEAK
hrt(3)=figure(3); 
subplot(2,2,1);
ha=plot(I,CDpktimedet(1,:),'k');
hb = get(ha,'children');
set(get(ha,'Parent'),'YScale','log');
set(gca,'xdir','reverse'); hold on;
plot(I,ones(1,length(I)),'g--');
legend('Data','CD=1');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH PEAK');
title('Bundle Data - PEAKS');

% PLOT THE COEFFICIENT OF DISPERSION FOR EACH TROUGH
hrt(3)=figure(3); 
subplot(2,2,2);
ha=plot(I,CDtrtimedet(1,:),'k');
hb = get(ha,'children');
set(get(ha,'Parent'),'YScale','log');
set(gca,'xdir','reverse'); hold on;
plot(I,ones(1,length(I)),'g--');
legend('Data','CD=1');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Dispersion - EACH TROUGH');
title('Bundle Data - TROUGHS');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
hrt(3)=figure(3); 
subplot(2,2,3);
ha=plot(I,CVpktimedet(1,:),'k');
hb = get(ha,'children');
set(get(ha,'Parent'),'YScale','log');
set(gca,'xdir','reverse'); hold on;
plot(I,1./sqrt(meanpktimedet),'k--')
plot(I,ones(1,length(I)),'g--');
legend('Data','1/sqrt(mean','CV=1');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH PEAK');
title('Bundle Data - Each Peak');

% PLOT THE COEFFICIENT OF VARIATION BETWEEN PEAK/PEAK and TROUGH/TROUGH
hrt(3)=figure(3); 
subplot(2,2,4);
ha=plot(I,CVtrtimedet(1,:),'k');
hb = get(ha,'children');
set(get(ha,'Parent'),'YScale','log');
set(gca,'xdir','reverse'); hold on;
plot(I,1./sqrt(meantrtimedet),'k--')
plot(I,ones(1,length(I)),'g--');
legend('Data','1/sqrt(mean','CV=1');
axis([mu(1) mu(end) 0.1 1000]);
xlabel('Control Parameter');ylabel('Coefficient of Variation - EACH TROUGH');
title('Bundle Data - Each Trough');

% PLOT THE AMPLITUDES FROM FFT
hrt(4)=figure(4);
subplot(3,1,2);
ha=plot(I,fftampldet(1,:),'k');
hb = get(ha,'children'); 
set(gca,'xdir','reverse'); hold on;
axis([mu(1) mu(end) 0 1.2*max(fftampldet)]);
xlabel('Control Parameter');ylabel('Amplitude from FFT');
title('Bundle Data');

% PLOT THE FREQUENCIES FROM FFT
hrt(4)=figure(4);
subplot(3,1,3);
ha=plot(I,fftfreqdet(1,:),'k');
hb = get(ha,'children'); 
set(gca,'xdir','reverse'); hold on;
axis([mu(1) mu(end) 0 1.2*max(fftfreqdet)]);
xlabel('Control Parameter');ylabel('Frequency from FFT');
title('Bundle Data');

% PLOT THE PEAK-TO-TROUGH AMPLITUDES
hrt(4)=figure(4);
subplot(3,1,1);
ha=plot(I,PKamplmeandet(1,:),'k');
hb = get(ha,'children'); 
set(gca,'xdir','reverse'); hold on;
axis([mu(1) mu(end) 0 1.2*max(PKamplmeandet)]);
xlabel('Control Parameter');ylabel('Peak-to-trough amplitude');
title('Bundle Data');



Ilength = ceil(length(I)/5);
for j = 1:Ilength
    if j < Ilength
        Iselect2(j,:) = I(1+(j-1)*5:5*j);
        Iselect(j,:) = 1+(j-1)*5:5*j;
    else
        Iselect2(j,:) = zeros(1,5);
        Iselect(j,:) = zeros(1,5);
        for k = 1:length(I(1+(j-1)*5:end))
            Iselect2(j,k) = I(1+(j-1)*5+(k-1));
            Iselect(j,k) = 1+(j-1)*5+(k-1);
        end
    end
end
tmin=round((length(Xd_dwnspl{1}{1}))*0.8);tmax=(length(Xd_dwnspl{1}{1}));    % CHOOSE
ymin=3*min(Xd_dwnspl{1}{1});ymax=3*max(Xd_dwnspl{1}{1});
yimin=2*ymin;yimax=2*ymax;
dwnsplquiver=20;quiverstart=tmin;quiverend=tmax;quiverscale=1;
dwnsplrealimag=10;realimagstart=tmin;realimagend=tmax;dwnspl2=dwnsplrealimag;
fmin = 0; fmax = 5*max(max(fftfreqdet));


% PLOT EXAMPLE TIME TRACES AND PHASE PORTRAITS
% Deterministic
n=1;

for k = 1:length(Xd_dwnspl{1})
    for j = 1:length(Xd_dwnspl)
        Xd_dwnspl2{n} = Xd_dwnspl{j}{k};
        fdetfft2{n} = fdetfft{j,k};
        Xd_dwnsplfft2{n} = Xd_dwnsplfft{j,k};
        n=n+1;
    end
end
Xd_dwnspl=Xd_dwnspl2;fdetfft=fdetfft2;Xd_dwnsplfft=Xd_dwnsplfft2;clear Xd_dwnspl2 Xd_dwnsplfft2 fdetfft2;
for j = 1:Ilength
for k = 1:sum((Iselect(j,:)~=0))
    hrt(4+j)=figure(4+j);
    sph=subplot(8,length(Iselect(1,:)),k);plot(Xscale2(k+(j-1)*5).*Xd_dwnspl{Iselect(j,k)},'k');axis([tmin tmax ymin ymax]);
    %xlabel('Time');ylabel('Position');
    title(sprintf('%s %s%s','Bundle Data','mu = ',num2str((Iselect2(j,k)))));
    spp = get(sph, 'pos');
    axis off;set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    xdr=real(hilbert(Xscale2(k+(j-1)*5).*Xd_dwnspl{Iselect(j,k)}));xdi=imag(hilbert(Xscale2(k+(j-1)*5).*Xd_dwnspl{Iselect(j,k)}));
    if length(unique(xdr)) > 1 && length(unique(xdi)) > 1
        sph=subplot(8,length(Iselect(1,:)),k+3*length(Iselect(1,:)));[bw dens mx my]=kde2d([xdr,xdi],2^8,[ymin,ymin],[ymax,ymax]);
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
        sph=subplot(8,length(Iselect(1,:)),k+length(Iselect(1,:)));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([mx(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        axis off;set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(8,length(Iselect(1,:)),k+2*length(Iselect(1,:)));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([mx(1,1) mx(1,end) my(1,1) my(end,1)]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        axis off;set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(8,length(Iselect(1,:)),k+3*length(Iselect(1,:)));pcolor(mx(1,:)',my(:,1)',dens);shading interp;load jetnew.mat;colormap(cjetnew);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        axis off;set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        [bw1,dens1,xmesh]=kde1d(xdr);dens1=dens1./sum(dens1);
        sph=subplot(8,length(Iselect(1,:)),k+4*length(Iselect(1,:)));plot(xmesh,dens1,'k');axis([mx(1,1) mx(1,end) 0 1.1*max(dens1)]);
        spp = get(sph, 'pos');
        axis off;set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);Xdpos=Xd_dwnspl{Iselect(j,k)}(tmin:dwnspl2:tmax);Xdvel=gradient(Xd_dwnspl{Iselect(j,k)}(tmin:dwnspl2:tmax));[bw dens mxn myn]=kde2d([Xdpos,Xdvel],2^8,[ymin,ymin],[ymax,ymax]);
        sph=subplot(8,length(Iselect(1,:)),k+5*length(Iselect(1,:)));plot(Xdpos,Xdvel);axis([mxn(1,1) mxn(1,end) myn(1,1) myn(end,1)]);
        %tt=findnearest(fdetfft{j,Iselect(j,k)},fftfreqdet(1,Iselect(j,k)));tt=tt(1);hold on;scatter(fdetfft{j,Iselect(j,k)}(tt),Xscale2(k+(j-1)*5).*Xd_dwnsplfft{1,Iselect(j,k)}(tt),'b.');
        spp = get(sph, 'pos');
        axis off;set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(8,length(Iselect(1,:)),k+6*length(Iselect(1,:)));pcolor(mx(1,:)',my(:,1)',dens);shading interp;load jetnew.mat;colormap(cjetnew);
        %tt=findnearest(fdetfft{j,Iselect(j,k)},fftfreqdet(1,Iselect(j,k)));tt=tt(1);hold on;scatter(fdetfft{j,Iselect(j,k)}(tt),Xscale2(k+(j-1)*5).*Xd_dwnsplfft{1,Iselect(j,k)}(tt),'b.');
        spp = get(sph, 'pos');
        axis off;set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(8,length(Iselect(1,:)),k+7*length(Iselect(1,:)));plot(fdetfft{Iselect(j,k)},Xscale2(k+(j-1)*5).*Xd_dwnsplfft{Iselect(j,k)},'k');axis([fmin fmax 0 2*max(Xd_dwnsplfft{Iselect(j,k)})]);
        %tt=findnearest(fdetfft{j,Iselect(j,k)},fftfreqdet(1,Iselect(j,k)));tt=tt(1);hold on;scatter(fdetfft{j,Iselect(j,k)}(tt),Xscale2(k+(j-1)*5).*Xd_dwnsplfft{1,Iselect(j,k)}(tt),'b.');
        spp = get(sph, 'pos');
        axis off;set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    else
        sph=subplot(8,length(Iselect(1,:)),k+length(Iselect(1,:)));plot(xdr(realimagstart:dwnsplrealimag:realimagend),xdi(realimagstart:dwnsplrealimag:realimagend),'k');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        axis off;set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
        sph=subplot(8,length(Iselect(1,:)),k+2*length(Iselect(1,:)));quiver(xdr(quiverstart:dwnsplquiver:quiverend),xdi(quiverstart:dwnsplquiver:quiverend),gradient(xdr(quiverstart:dwnsplquiver:quiverend)),gradient(xdi(quiverstart:dwnsplquiver:quiverend)),quiverscale,'k');axis([yimin yimax yimin yimax]);
        %xlabel('Real');ylabel('Imaginary');
        spp = get(sph, 'pos');
        axis off;set(sph, 'Position', [spp(1) 1*spp(2) 1.3*spp(3) 1.4*spp(4)]);
    end
end
end

if saveyn == 1
    disp('Saving...');
    filename2 = filename(1:end-4);
    savefig(hrt,sprintf('%s%s',filename2,'-Figures.fig'));
    disp('Finished.');
else
    disp('Not saved.');
end

end
