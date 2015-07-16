% Initialize
clear pkspikeratestoAVG pkspikeratestoSEM pkspikeratedetAVG pkspikeratedetSEM CDstopkAVG CDdetpkAVG CDdetpkSEM CDstopkSEM RMSmagdetAVG RMSmagstoAVG RMSmagdetSEM RMSmagstoSEM
clear IEIvarpkstoAVG IEIvarpkstoSEM IEImeanpkstoAVG IEImeanpkstoSEM
pkspikeratedet1=pkspikeratedet2;CDdetpk1=CDdetpk2;PKamplmeandet1=PKamplmeandet2;fftampldet1=fftampldet2;fftfreqdet1=fftfreqdet2;

% Find averages
for k = 1:5
    for l = 1:3
        for j = 1:500
            pkspikeratestoAVG(j,k,l)= mean([pkspikeratesto_1(j,k,l) pkspikeratesto_2(j,k,l) pkspikeratesto_3(j,k,l) pkspikeratesto_4(j,k,l) pkspikeratesto_5(j,k,l)]);
            pkspikeratestoSEM(j,k,l)= std([pkspikeratesto_1(j,k,l) pkspikeratesto_2(j,k,l) pkspikeratesto_3(j,k,l) pkspikeratesto_4(j,k,l) pkspikeratesto_5(j,k,l)])/sqrt(5);
            pkspikeratedetAVG(j,k,l)= mean([pkspikeratedet1(j,k,l) pkspikeratedet2(j,k,l) pkspikeratedet3(j,k,l) pkspikeratedet4(j,k,l) pkspikeratedet5(j,k,l)]);
            pkspikeratedetSEM(j,k,l)= std([pkspikeratedet1(j,k,l) pkspikeratedet2(j,k,l) pkspikeratedet3(j,k,l) pkspikeratedet4(j,k,l) pkspikeratedet5(j,k,l)])/sqrt(5);
            CDstopkAVG(j,k,l)= mean([CDstopk_1(j,k,l) CDstopk_2(j,k,l) CDstopk_3(j,k,l) CDstopk_4(j,k,l) CDstopk_5(j,k,l)]);
            CDstopkSEM(j,k,l)= std([CDstopk_1(j,k,l) CDstopk_2(j,k,l) CDstopk_3(j,k,l) CDstopk_4(j,k,l) CDstopk_5(j,k,l)])/sqrt(5);
            CDdetpkAVG(j,k,l)= mean([CDdetpk1(j,k,l) CDdetpk2(j,k,l) CDdetpk3(j,k,l) CDdetpk4(j,k,l) CDdetpk5(j,k,l)]);
            CDdetpkSEM(j,k,l)= std([CDdetpk1(j,k,l) CDdetpk2(j,k,l) CDdetpk3(j,k,l) CDdetpk4(j,k,l) CDdetpk5(j,k,l)])/sqrt(5);
            %PKamplmeandetAVG(j,k,l)= mean([PKamplmeandet1(j,k,l) PKamplmeandet2(j,k,l) PKamplmeandet3(j,k,l) PKamplmeandet4(j,k,l) PKamplmeandet5(j,k,l)]);
            %PKamplmeandetSEM(j,k,l)= std([PKamplmeandet1(j,k,l) PKamplmeandet2(j,k,l) PKamplmeandet3(j,k,l) PKamplmeandet4(j,k,l) PKamplmeandet5(j,k,l)])/sqrt(5);
            %PKamplmeanstoAVG(j,k,l)= mean([PKamplmeansto1(j,k,l) PKamplmeansto2(j,k,l) PKamplmeansto3(j,k,l) PKamplmeansto4(j,k,l) PKamplmeansto5(j,k,l)]);
            %PKamplmeanstoSEM(j,k,l)= std([PKamplmeansto1(j,k,l) PKamplmeansto2(j,k,l) PKamplmeansto3(j,k,l) PKamplmeansto4(j,k,l) PKamplmeansto5(j,k,l)])/sqrt(5);
            if exist('Xdet1')==1
            L = length(Xdet1{1,1,1});
            for m = 1:5
                RMSmagdet(m) = std(Xdet1{j,k,l}(1,1+(m-1)*floor(L/5):m*floor(L/5))-mean(Xdet1{j,k,l}(1,1+(m-1)*floor(L/5):m*floor(L/5))));
                RMSmagsto(m) = std(Xsto1{j,k,l}(1,1+(m-1)*floor(L/5):m*floor(L/5))-mean(Xsto1{j,k,l}(1,1+(m-1)*floor(L/5):m*floor(L/5))));
            end
            RMSmagdetAVG(j,k,l) = mean(RMSmagdet);RMSmagdetSEM(j,k,l) = std(RMSmagdet)/sqrt(5);
            RMSmagstoAVG(j,k,l) = mean(RMSmagsto);RMSmagstoSEM(j,k,l) = std(RMSmagsto)/sqrt(5);
            else
                RMSmagdetAVG(j,k,l) = NaN;
                RMSmagstoAVG(j,k,l) = NaN;
            end
            clear RMSmagdet RMSmagsto
            IEIvarpkstoAVG(j,k,l)= mean([IEIvarpksto_1(j,k,l) IEIvarpksto_2(j,k,l) IEIvarpksto_3(j,k,l) IEIvarpksto_4(j,k,l) IEIvarpksto_5(j,k,l)]);
            IEIvarpkstoSEM(j,k,l)= std([IEIvarpksto_1(j,k,l) IEIvarpksto_2(j,k,l) IEIvarpksto_3(j,k,l) IEIvarpksto_4(j,k,l) IEIvarpksto_5(j,k,l)])/sqrt(5);
            IEImeanpkstoAVG(j,k,l)= mean([IEImeanpksto_1(j,k,l) IEImeanpksto_2(j,k,l) IEImeanpksto_3(j,k,l) IEImeanpksto_4(j,k,l) IEImeanpksto_5(j,k,l)]);
            IEImeanpkstoSEM(j,k,l)= std([IEImeanpksto_1(j,k,l) IEImeanpksto_2(j,k,l) IEImeanpksto_3(j,k,l) IEImeanpksto_4(j,k,l) IEImeanpksto_5(j,k,l)])/sqrt(5);
    
        end
    end
end

for j = 1:3
    for k = 1:500
        PKamplmeandetAVG(j,k)= mean([PKamplmeandet1(j,k) PKamplmeandet2(j,k) PKamplmeandet3(j,k) PKamplmeandet4(j,k) PKamplmeandet5(j,k)]);
        PKamplmeandetSEM(j,k)= std([PKamplmeandet1(j,k) PKamplmeandet2(j,k) PKamplmeandet3(j,k) PKamplmeandet4(j,k) PKamplmeandet5(j,k)])/sqrt(5);
        PKamplmeanstoAVG(j,k)= mean([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)]);
        PKamplmeanstoSEM(j,k)= std([PKamplmeansto1(j,k) PKamplmeansto2(j,k) PKamplmeansto3(j,k) PKamplmeansto4(j,k) PKamplmeansto5(j,k)])/sqrt(5);
        fftampldetAVG(j,k)= mean([fftampldet1(j,k) fftampldet2(j,k) fftampldet3(j,k) fftampldet4(j,k) fftampldet5(j,k)]);
        fftampldetSEM(j,k)= std([fftampldet1(j,k) fftampldet2(j,k) fftampldet3(j,k) fftampldet4(j,k) fftampldet5(j,k)])/sqrt(5);
        fftamplstoAVG(j,k)= mean([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)]);
        fftamplstoSEM(j,k)= std([fftamplsto1(j,k) fftamplsto2(j,k) fftamplsto3(j,k) fftamplsto4(j,k) fftamplsto5(j,k)])/sqrt(5);
        fftfreqstoAVG(j,k)= mean([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)]);
        fftfreqstoSEM(j,k)= std([fftfreqsto1(j,k) fftfreqsto2(j,k) fftfreqsto3(j,k) fftfreqsto4(j,k) fftfreqsto5(j,k)])/sqrt(5);
        fftfreqdetAVG(j,k)= mean([fftfreqdet1(j,k) fftfreqdet2(j,k) fftfreqdet3(j,k) fftfreqdet4(j,k) fftfreqdet5(j,k)]);
        fftfreqdetSEM(j,k)= std([fftfreqdet1(j,k) fftfreqdet2(j,k) fftfreqdet3(j,k) fftfreqdet4(j,k) fftfreqdet5(j,k)])/sqrt(5);
    end
end


%% Plot the data
thr=3;
setfiguredefaults(4);
stiffind=1;
mu=F;
addplot=1;
ctrlparam='Control parameter (Force)';
noiselevel=[0.05 0.1 0.2 0.4];
%noiselevel=[0.5 1];
colors1={'r','m','b'};


% Spike rate
for j = 1:3
    for k = 1:5
        hh(j)=figure(j);
        if addplot==1
            an=findobj(gcf,'type','axes');anm=1:length(an);anm=fliplr(anm);
            hold on;subplot(3,2,k,an(anm(k)));hold on;
        else
            subplot(3,2,k);
        end
        plot(mu,pkspikeratedetAVG(:,k,1),'k--');hold on;
        ha=errorbar(mu,pkspikeratestoAVG(:,k,j),pkspikeratestoSEM(:,k,j),colors1{thr});
        hb = get(ha,'children');  
        Xdata = get(hb(2),'Xdata');
        temp = 4:3:length(Xdata);
        temp(3:3:end) = [];
        % xleft and xright contain the indices of the left and right
        %  endpoints of the horizontal lines
        xleft = temp; xright = temp+1; 
        % Change line length
        Xdata(xleft) =  Xdata(xleft) + .02;
        Xdata(xright) =  Xdata(xright) - .02;
        set(hb(2),'Xdata',Xdata)
        set(gca,'xdir','Reverse');
        axis([mu(1) mu(end) 0 max(max(max(pkspikeratestoAVG)))]);
        xlabel(ctrlparam);ylabel('Spike rate');
        title(['k = ' num2str(k1(k))]);
    end
    set(hh(j),'Name',['Spike Rate; Noise = ' num2str(noiselevel(j))]);
end

% Coefficient of Dispersion
for j = 1:3
    for k = 1:5
        hh(3+j)=figure(3+j);
        if addplot==1
            an=findobj(gcf,'type','axes');anm=1:length(an);anm=fliplr(anm);
            hold on;subplot(3,2,k,an(anm(k)));hold on;
        else
            subplot(3,2,k);
        end
        plot(mu,CDdetpkAVG(:,k,1),'k--');hold on;
        ha=errorbar(mu,CDstopkAVG(:,k,j),CDstopkSEM(:,k,j),colors1{thr});
        hb = get(ha,'children');  
        Xdata = get(hb(2),'Xdata');
        temp = 4:3:length(Xdata);
        temp(3:3:end) = [];
        % xleft and xright contain the indices of the left and right
        %  endpoints of the horizontal lines
        xleft = temp; xright = temp+1; 
        % Change line length
        Xdata(xleft) =  Xdata(xleft) + .02;
        Xdata(xright) =  Xdata(xright) - .02;
        set(hb(2),'Xdata',Xdata)
        set(gca,'xdir','Reverse');
        set(get(ha,'Parent'),'YScale','log');
        axis([mu(1) mu(end) 0.01 1000]);
        xlabel(ctrlparam);ylabel('Coefficient of dispersion');
        title(['k = ' num2str(k1(k))]);
        hold on;plot(mu,ones(1,length(mu)),'g--');
    end
    set(hh(3+j),'Name',['Coefficient of Dispersion; Noise = ' num2str(noiselevel(j))]);
end


% RMS Magnitude
for j = 1:3
    for k = 1:5
        hh(6+j)=figure(6+j);
        if addplot==1
            an=findobj(gcf,'type','axes');anm=1:length(an);anm=fliplr(anm);hold on;subplot(3,2,k,an(anm(k)));hold on;
        else
            subplot(3,2,k);
        end
        plot(mu,RMSmagdetAVG(:,k,1),'k--');hold on;
        ha=errorbar(mu,RMSmagstoAVG(:,k,j),RMSmagstoSEM(:,k,j),colors1{thr});
        hb = get(ha,'children');  
        Xdata = get(hb(2),'Xdata');
        temp = 4:3:length(Xdata);
        temp(3:3:end) = [];
        % xleft and xright contain the indices of the left and right
        %  endpoints of the horizontal lines
        xleft = temp; xright = temp+1; 
        % Change line length
        Xdata(xleft) =  Xdata(xleft) + .02;
        Xdata(xright) =  Xdata(xright) - .02;
        set(hb(2),'Xdata',Xdata)
        set(gca,'xdir','Reverse');
        axis([mu(1) mu(end) 0 max(max(max(RMSmagstoAVG())))]);
        xlabel(ctrlparam);ylabel('RMS magnitude');
        title(['k = ' num2str(k1(k))]);
    end
    set(hh(6+j),'Name',['RMS magnitude; Noise = ' num2str(noiselevel(j))]);
end


% Peak-to-peak amplitude
for j = 1:3
    for k = stiffind
        hh(9+j)=figure(9+j);
        if addplot==1
            an=findobj(gcf,'type','axes');anm=1:length(an);anm=fliplr(anm);hold on;subplot(3,2,k,an(anm(k)));hold on;
        else
            subplot(3,2,k);
        end
        plot(mu,PKamplmeandetAVG(1,:),'k--');hold on;
        ha=errorbar(mu,PKamplmeanstoAVG(j,:),PKamplmeanstoSEM(j,:),colors1{thr});
        hb = get(ha,'children');  
        Xdata = get(hb(2),'Xdata');
        temp = 4:3:length(Xdata);
        temp(3:3:end) = [];
        % xleft and xright contain the indices of the left and right
        %  endpoints of the horizontal lines
        xleft = temp; xright = temp+1; 
        % Change line length
        Xdata(xleft) =  Xdata(xleft) + .02;
        Xdata(xright) =  Xdata(xright) - .02;
        set(hb(2),'Xdata',Xdata)
        set(gca,'xdir','Reverse');
        axis([mu(1) mu(end) 0 max(max(max(PKamplmeansto())))]);
        xlabel(ctrlparam);ylabel('Peak-to-peak amplitude');
        title(['k = ' num2str(k1(k))]);
    end
    set(hh(9+j),'Name',['Peak-to-peak amplitude; Noise = ' num2str(noiselevel(j))]);
end


% Amplitude from FFT
for j = 1:3
    for k = stiffind
        hh(12+j)=figure(12+j);
        if addplot==1             
            an=findobj(gcf,'type','axes');anm=1:length(an);anm=fliplr(anm);hold on;subplot(3,2,k,an(stiffind));hold on;
        else
            subplot(3,2,k);      
        end
        plot(mu,fftampldetAVG(1,:),'k--');hold on;
        ha=errorbar(mu,fftamplstoAVG(j,:),fftamplstoSEM(j,:),colors1{thr});
        hb = get(ha,'children');  
        Xdata = get(hb(2),'Xdata');
        temp = 4:3:length(Xdata);
        temp(3:3:end) = [];
        % xleft and xright contain the indices of the left and right
        %  endpoints of the horizontal lines
        xleft = temp; xright = temp+1; 
        % Change line length
        Xdata(xleft) =  Xdata(xleft) + .02;
        Xdata(xright) =  Xdata(xright) - .02;
        set(hb(2),'Xdata',Xdata)
        set(gca,'xdir','Reverse');
        axis([mu(1) mu(end) 0 max(max(max(fftamplstoAVG())))]);
        xlabel(ctrlparam);ylabel('Amplitude (FFT)');
        title(['k = ' num2str(k1(k))]);
    end
    set(hh(12+j),'Name',['Amplitude from FFT; Noise = ' num2str(noiselevel(j))]);
end

% Frequency from FFT
for j = 1:3
    for k = stiffind
        hh(15+j)=figure(15+j);
            if addplot==1             
                an=findobj(gcf,'type','axes');anm=1:length(an);anm=fliplr(anm);hold on;subplot(3,2,k,an(stiffind));hold on;         
            else
                subplot(3,2,k);
            end
        plot(mu,fftfreqdetAVG(1,:),'k--');hold on;
        ha=errorbar(mu,fftfreqstoAVG(j,:),fftfreqstoSEM(j,:),colors1{thr});
        hb = get(ha,'children');  
        Xdata = get(hb(2),'Xdata');
        temp = 4:3:length(Xdata);
        temp(3:3:end) = [];
        % xleft and xright contain the indices of the left and right
        %  endpoints of the horizontal lines
        xleft = temp; xright = temp+1; 
        % Change line length
        Xdata(xleft) =  Xdata(xleft) + .02;
        Xdata(xright) =  Xdata(xright) - .02;
        set(hb(2),'Xdata',Xdata)
        set(gca,'xdir','Reverse');
        axis([mu(1) mu(end) 0 max(max(max(fftfreqstoAVG())))]);
        xlabel(ctrlparam);ylabel('Frequency (FFT)');
        title(['k = ' num2str(k1(k))]);
    end
    set(hh(15+j),'Name',['Frequency from FFT; Noise = ' num2str(noiselevel(j))]);
end

%% Save figures

savefig(hh,'/Volumes/Promise Pegasus/manual backup/Simulation Data/toymodel/Force-figures.fig');

%savefig(hh,'/Volumes/Promise Pegasus/manual backup/Simulation Data/toymodel-highnoise/Force-figures.fig');
