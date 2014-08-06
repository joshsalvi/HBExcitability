clear all; close all;
display('Importing...');
%load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp7-kp2-start1end2e5-xfish1.0Noise.mat');
load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/noisysinewave33movingavg.mat');

wins = 1:1:1;                                    % WINDOW RANGE SET TO ZERO (not in use)
nup = 0:0.1:3;
ndown = 0:0.1:3;

display('Find up events...');
up_events(1,1:length(nup),1:length(wins))=0;    %initialize
up_ind(1:length(movingavgX),1:length(nup),1:length(wins))=0;
for l = 1:length(wins)
for j = 1:length(nup)
    Mvec_up = zeros(1,length(Xvec));
for i = 1:length(movingavgX)
    if i <= win
        clear a;
        a = Xvec(1:round(win/2)+i-1) > movingavgX(i) + nup(j)*movingstdX(i);            % left edge
        Mvec_up(1:round(win/2)+i-1) = Mvec_up(1:round(win/2)+i-1)+a(:)';
    else if i < length(movingavgX) - win
            clear a;
            a = Xvec(i-floor(win/2):round(win/2)+i-1) > movingavgX(i) + nup(j)*movingstdX(i);           % middle
            Mvec_up(i-floor(win/2):round(win/2)+i-1) = Mvec_up(i-floor(win/2):round(win/2)+i-1)+a(:)';
        else
            clear a;
            a = Xvec(i-floor(win/2):end) > movingavgX(i) + nup(j)*movingstdX(i);          % right edge
            Mvec_up(i-floor(win/2):end) = Mvec_up(i-floor(win/2):end)+a(:)';
        end
    end
end
Mvec_up(Mvec_up > 0) = 1;   % reduce multiple counts to 1
up_ind(:,j,l) = Mvec_up;
clear updiffend updiffstart updiff updiffind
updiffend = find(diff(up_ind(:,j,l))==-1);           % concatenate nearby points, stepping through window sizes
updiffstart = find(diff(up_ind(:,j,l))==1);

if isempty(updiffend)==0 || isempty(updiffstart)==0
if updiffend(1)<updiffstart(1)
    updiffend(1)=[];
end
for k = 2:length(updiffend)
    if abs(updiffend(k) - updiffstart(k-1)) < wins(l)
        if updiffend(k)>updiffstart(k-1)
            up_ind(updiffstart(k-1):updiffend(k),j,l) = 1;
        else
            up_ind(updiffend(k-1):updiffstart(k),j,l) = 1;
        end
    end
end
end
%}
up_events(:,j,l) = floor(sum(abs(diff(up_ind(:,j,l))))/2);
end
if wins(l) == wins(round(length(wins)/2))
    display('50%');
end
end
events_up_diff = diff(up_events);


display('Find down events...');
down_events(1,1:length(ndown),1:length(wins))=0;
down_ind(1:length(movingavgX),1:length(ndown),1:length(wins))=0;
for l = 1:length(wins)
for j = 1:length(ndown)
    Mvec_down = zeros(1,length(Xvec));
for i = 1:length(movingavgX)
    if i <= win
        clear a;
        a = Xvec(1:round(win/2)+i-1) < movingavgX(i) - ndown(j)*movingstdX(i);            % left edge
        Mvec_down(1:round(win/2)+i-1) = Mvec_down(1:round(win/2)+i-1)+a(:)';
    else if i < length(movingavgX) - win
            clear a;
            a = Xvec(i-floor(win/2):round(win/2)+i-1) < movingavgX(i) - ndown(j)*movingstdX(i);           % middle
            Mvec_down(i-floor(win/2):round(win/2)+i-1) = Mvec_down(i-floor(win/2):round(win/2)+i-1)+a(:)';
        else
            clear a;
            a = Xvec(i-floor(win/2):end) < movingavgX(i) - ndown(j)*movingstdX(i);          % right edge
            Mvec_down(i-floor(win/2):end) = Mvec_down(i-floor(win/2):end)+a(:)';
        end
    end
end
Mvec_down(Mvec_down > 0) = 1;   % reduce multiple counts to 1
down_ind(:,j,l) = Mvec_down;
clear downdiffend downdiffstart downdiff downdiffind
downdiffend = find(diff(down_ind(:,j,l))==-1);           % concatenate nearby points, stepping through window sizes
downdiffstart = find(diff(down_ind(:,j,l))==1);

if isempty(downdiffend)==0 || isempty(downdiffstart)==0
if downdiffend(1)<downdiffstart(1)
    downdiffend(1)=[];
end
for k = 2:length(downdiffend)
    if abs(downdiffend(k) - downdiffstart(k-1)) < wins(l)
        if downdiffend(k)>downdiffstart(k-1)
            down_ind(downdiffstart(k-1):downdiffend(k),j,l) = 1;
        else
            down_ind(downdiffend(k):downdiffstart(k-1),j,l) = 1;
        end
    end
end
end
%}
down_events(:,j,l) = floor(sum(abs(diff(down_ind(:,j,l))))/2);
end
if wins(l) == wins(round(length(wins)/2))
    display('50%');
end
end
events_down_diff = diff(down_events);



for l = 1:length(wins)
    clear nupconv_hold ndownconv_hold
uppeaks{l} = findpeaks(up_events(:,:,l));
downpeaks{l} = findpeaks(down_events(:,:,l));
if length(uppeaks{l}) > 1
    nupconv_hold = find(abs(diff(up_events(1,uppeaks{l}(end-1):uppeaks{l}(end),l)))==min(abs(diff(up_events(1,uppeaks{l}(end-1):uppeaks{l}(end),l)))))+uppeaks{l}(end-1);
else
    nupconv_hold = uppeaks{l}(1);
end
if length(downpeaks{l}) > 1
    ndownconv_hold = find(abs(diff(down_events(1,downpeaks{l}(end-1):downpeaks{l}(end),l)))==min(abs(diff(down_events(1,downpeaks{l}(end-1):downpeaks{l}(end),l)))))+downpeaks{l}(end-1);
else
    ndownconv_hold = downpeaks{l}(1);
end
nupconv(l) = nupconv_hold(1);
ndownconv(l) = ndownconv_hold(1);
end

uppeaks = findpeaks(up_events(:,:,1));
downpeaks = findpeaks(down_events(:,:,1));


for i = uppeaks(end):length(nup)
    clear upindex upinddiff spikeindup_start spikeindup_end  
    upindex = up_ind(:,i);
    upinddiff = diff(upindex);
    spikeindup_start = find(upinddiff==1);
    spikeindup_end = find(upinddiff==-1);
    if isempty(spikeindup_start)
        break;
    end
    if isempty(spikeindup_end)
        break;
    end    
if spikeindup_start(1) > spikeindup_end(1)
    spikeindup_end(1) = [];
end

for j = 1:length(spikeindup_start)-1
    RT_up(j,i) = (time(spikeindup_end(j))-time(spikeindup_start(j))) * 1e3;         % convert to ms
end    
for j = 1:length(spikeindup_start)-1
    IEI_up(j,i) = (time(spikeindup_start(j+1))-time(spikeindup_end(j))) * 1e3;         % convert to ms
end    
end

for i = downpeaks(end):length(ndown)
    clear upindex upinddiff spikeindup_start spikeindup_end  
    downindex = down_ind(:,i);
    downinddiff = diff(downindex);
    spikeinddown_start = find(downinddiff==1);
    spikeinddown_end = find(downinddiff==-1);
    if isempty(spikeinddown_start)
        break;
    end
    if isempty(spikeinddown_end)
        break;
    end     
if spikeinddown_start(1) > spikeinddown_end(1)
    spikeinddown_end(1) = [];
end

for j = 1:length(spikeinddown_start)-1
    RT_down(j,i) = (time(spikeinddown_end(j))-time(spikeinddown_start(j))) * 1e3;         % convert to ms
end    
for j = 1:length(spikeinddown_start)-1
    IEI_down(j,i) = (time(spikeinddown_start(j+1))-time(spikeinddown_end(j))) * 1e3;         % convert to ms
end    
end

% analyze everything
sizeup=size(RT_up);
sizedown=size(RT_down);
for i = uppeaks(end):sizeup(2)
    RT_updist{i} = RT_up(:,i);
    RT_updist{i}(RT_updist{i}==0)=[];
    IEI_updist{i} = IEI_up(:,i);
    IEI_updist{i}(IEI_updist{i}==0)=[];
    if length(RT_updist{i}) <4
        RT_updist{i}(length(RT_updist{i})+1:4) = 0;
    end
    if length(IEI_updist{i}) <4
        IEI_updist{i}(length(IEI_updist{i})+1:4) = 0;
    end 
    [~,RTupgauss_p(i)] = lillietest(RT_updist{i});
    [~,RTupexp_p(i)] = lillietest(RT_updist{i},'Distr','exp');
    [~,IEIupgauss_p(i)] = lillietest(IEI_updist{i});
    [~,IEIupexp_p(i)] = lillietest(IEI_updist{i},'Distr','exp');
end

for i = downpeaks(end):sizedown(2)
    RT_downdist{i} = RT_down(:,i);
    RT_downdist{i}(RT_downdist{i}==0)=[];
    IEI_downdist{i} = IEI_down(:,i);
    IEI_downdist{i}(IEI_downdist{i}==0)=[];
    if length(RT_downdist{i}) <4
        RT_downdist{i}(length(RT_downdist{i})+1:4) = 0;
    end
    if length(IEI_downdist{i}) <4
        IEI_downdist{i}(length(IEI_downdist{i})+1:4) = 0;
    end  
    [~,RTdowngauss_p(i)] = lillietest(RT_downdist{i});
    [~,RTdownexp_p(i)] = lillietest(RT_downdist{i},'Distr','exp');
    [~,IEIdowngauss_p(i)] = lillietest(IEI_downdist{i});
    [~,IEIdownexp_p(i)] = lillietest(IEI_downdist{i},'Distr','exp');
end


figure(1);
subplot(2,2,1);plot(nup,up_events,'k');hold on;scatter(nup(uppeaks),up_events(uppeaks),'ro');plot(nup(uppeaks(end):end),up_events(uppeaks(end):end),'r');
ylabel('Events');xlabel('n up');title('Up Events and Peaks');
subplot(2,2,2);plot(ndown,down_events,'k');hold on;scatter(ndown(downpeaks),down_events(downpeaks),'ro');plot(ndown(downpeaks(end):end),down_events(downpeaks(end):end),'r');
ylabel('Events');xlabel('n down');title('Down Events and Peaks');
subplot(2,2,3);plot(nup(1:length(IEIupgauss_p)),IEIupgauss_p,'r');hold on;plot(nup(1:length(IEIupexp_p)),IEIupexp_p,'b');title('IEI Up p-values (red=gauss,blue=exp)');
subplot(2,2,4);plot(ndown(1:length(IEIdowngauss_p)),IEIdowngauss_p,'r');hold on;plot(ndown(1:length(IEIdownexp_p)),IEIdownexp_p,'b');title('IEI Up p-values (red=gauss,blue=exp)');

% PLOT TIME TRACES AND HISTOGRAMS
figure(2);
subplot(5,3,1);
plot(time,Xvec);hold on; ylabel('displacement');
plot(time(down_ind(:,downpeaks(end))==1),Xvec(down_ind(:,downpeaks(end))==1),'r.');
plot(time(up_ind(:,uppeaks(end))==1),Xvec(up_ind(:,uppeaks(end))==1),'g.');
title(sprintf('%s %s%s %s%s %s%s%s%s','Spikes ','nup= ',num2str(nup(uppeaks(end))),'ndown= ',num2str(ndown(downpeaks(end))),'eventsup/eventsdown= ',num2str(up_events(1,uppeaks(end))),'/',num2str(down_events(1,downpeaks(end)))));
binsizeup = 2*iqr(IEI_updist{uppeaks(end)})*length(IEI_updist{uppeaks(end)})^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(IEI_downdist{downpeaks(end)})*length(IEI_downdist{downpeaks(end)})^(-1/3);
nbinsup = round((max(IEI_updist{uppeaks(end)}) - min(IEI_updist{uppeaks(end)}))/binsizeup);
nbinsdown = round((max(IEI_downdist{downpeaks(end)}) - min(IEI_downdist{downpeaks(end)}))/binsizedown);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 | isnan(nbinsup)
        nbinsup=2;
    end    
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
else if isempty(nbinsdown) ==1 | isnan(nbinsdown)
        nbinsdown=2;
    end
end
subplot(5,3,2);hist(IEI_updist{uppeaks(end)},nbinsup);title('Up Events'); xlabel(num2str(mean(IEI_updist{uppeaks(end)})));
subplot(5,3,3);hist(IEI_downdist{downpeaks(end)},nbinsdown);title('Down Events'); xlabel(num2str(mean(IEI_downdist{downpeaks(end)})));

inputup= input(sprintf('%s%s %s%s %s','Minimum UP= ',num2str(uppeaks(end)),'; Maximum= ',num2str(sizeup(2)),'; Choice for up? #'));
inputdown= input(sprintf('%s%s %s%s %s','Minimum DOWN= ',num2str(downpeaks(end)),'; Maximum= ',num2str(sizedown(2)),'; Choice for down? #'));

subplot(5,3,4);
plot(time,Xvec);hold on; ylabel('displacement');
plot(time(down_ind(:,inputdown)==1),Xvec(down_ind(:,inputdown)==1),'r.');
plot(time(up_ind(:,inputup)==1),Xvec(up_ind(:,inputup)==1),'g.');
title(sprintf('%s %s%s %s%s %s%s%s%s','Spikes ','nup= ',num2str(nup(inputup)),'ndown= ',num2str(ndown(inputdown)),'eventsup/eventsdown= ',num2str(up_events(1,inputup)),'/',num2str(down_events(1,inputdown))));
binsizeup = 2*iqr(IEI_updist{inputup})*length(IEI_updist{inputup})^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(IEI_downdist{inputdown})*length(IEI_downdist{inputdown})^(-1/3);
nbinsup = round((max(IEI_updist{inputup}) - min(IEI_updist{inputup}))/binsizeup);
nbinsdown = round((max(IEI_downdist{inputdown}) - min(IEI_downdist{inputdown}))/binsizedown);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 | isnan(nbinsup)
        nbinsup=2;
    end    
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
else if isempty(nbinsdown) ==1 | isnan(nbinsdown)
        nbinsdown=2;
    end
end
subplot(5,3,5);hist(IEI_updist{inputup},nbinsup);title('Up Events'); xlabel(num2str(mean(IEI_updist{inputup})));
subplot(5,3,6);hist(IEI_downdist{inputdown},nbinsdown);title('Down Events');xlabel(num2str(mean(IEI_downdist{inputdown})));

inputup= input(sprintf('%s%s %s%s %s','Minimum UP= ',num2str(uppeaks(end)),'; Maximum= ',num2str(sizeup(2)),'; Choice for up? #'));
inputdown= input(sprintf('%s%s %s%s %s','Minimum DOWN= ',num2str(downpeaks(end)),'; Maximum= ',num2str(sizedown(2)),'; Choice for down? #'));

subplot(5,3,7);
plot(time,Xvec);hold on; ylabel('displacement');
plot(time(down_ind(:,inputdown)==1),Xvec(down_ind(:,inputdown)==1),'r.');
plot(time(up_ind(:,inputup)==1),Xvec(up_ind(:,inputup)==1),'g.');
title(sprintf('%s %s%s %s%s %s%s%s%s','Spikes ','nup= ',num2str(nup(inputup)),'ndown= ',num2str(ndown(inputdown)),'eventsup/eventsdown= ',num2str(up_events(1,inputup)),'/',num2str(down_events(1,inputdown))));
binsizeup = 2*iqr(IEI_updist{inputup})*length(IEI_updist{inputup})^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(IEI_downdist{inputdown})*length(IEI_downdist{inputdown})^(-1/3);
nbinsup = round((max(IEI_updist{inputup}) - min(IEI_updist{inputup}))/binsizeup);
nbinsdown = round((max(IEI_downdist{inputdown}) - min(IEI_downdist{inputdown}))/binsizedown);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 | isnan(nbinsup)
        nbinsup=2;
    end    
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
else if isempty(nbinsdown) ==1 | isnan(nbinsdown)
        nbinsdown=2;
    end
end
subplot(5,3,8);hist(IEI_updist{inputup},nbinsup);title('Up Events'); xlabel(num2str(mean(IEI_updist{inputup})));
subplot(5,3,9);hist(IEI_downdist{inputdown},nbinsdown);title('Down Events'); xlabel(num2str(mean(IEI_downdist{inputdown})));

inputup= input(sprintf('%s%s %s%s %s','Minimum UP= ',num2str(uppeaks(end)),'; Maximum= ',num2str(sizeup(2)),'; Choice for up? #'));
inputdown= input(sprintf('%s%s %s%s %s','Minimum DOWN= ',num2str(downpeaks(end)),'; Maximum= ',num2str(sizedown(2)),'; Choice for down? #'));

subplot(5,3,10);
plot(time,Xvec);hold on; ylabel('displacement');
plot(time(down_ind(:,inputdown)==1),Xvec(down_ind(:,inputdown)==1),'r.');
plot(time(up_ind(:,inputup)==1),Xvec(up_ind(:,inputup)==1),'g.');
title(sprintf('%s %s%s %s%s %s%s%s%s','Spikes ','nup= ',num2str(nup(inputup)),'ndown= ',num2str(ndown(inputdown)),'eventsup/eventsdown= ',num2str(up_events(1,inputup)),'/',num2str(down_events(1,inputdown))));
binsizeup = 2*iqr(IEI_updist{inputup})*length(IEI_updist{inputup})^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(IEI_downdist{inputdown})*length(IEI_downdist{inputdown})^(-1/3);
nbinsup = round((max(IEI_updist{inputup}) - min(IEI_updist{inputup}))/binsizeup);
nbinsdown = round((max(IEI_downdist{inputdown}) - min(IEI_downdist{inputdown}))/binsizedown);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 | isnan(nbinsup)
        nbinsup=2;
    end    
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
else if isempty(nbinsdown) ==1 | isnan(nbinsdown)
        nbinsdown=2;
    end
end
subplot(5,3,11);hist(IEI_updist{inputup},nbinsup);title('Up Events'); xlabel(num2str(mean(IEI_updist{inputup})));
subplot(5,3,12);hist(IEI_downdist{inputdown},nbinsdown);title('Down Events'); xlabel(num2str(mean(IEI_downdist{inputdown})));

inputup= input(sprintf('%s%s %s%s %s','Minimum UP= ',num2str(uppeaks(end)),'; Maximum= ',num2str(sizeup(2)),'; Choice for up? #'));
inputdown= input(sprintf('%s%s %s%s %s','Minimum DOWN= ',num2str(downpeaks(end)),'; Maximum= ',num2str(sizedown(2)),'; Choice for down? #'));

subplot(5,3,13);
plot(time,Xvec);hold on; ylabel('displacement');
plot(time(down_ind(:,inputdown)==1),Xvec(down_ind(:,inputdown)==1),'r.');
plot(time(up_ind(:,inputup)==1),Xvec(up_ind(:,inputup)==1),'g.');
title(sprintf('%s %s%s %s%s %s%s%s%s','Spikes ','nup= ',num2str(nup(inputup)),'ndown= ',num2str(ndown(inputdown)),'eventsup/eventsdown= ',num2str(up_events(1,inputup)),'/',num2str(down_events(1,inputdown))));
binsizeup = 2*iqr(IEI_updist{inputup})*length(IEI_updist{inputup})^(-1/3);      % freedman-diaconis rule
binsizedown = 2*iqr(IEI_downdist{inputdown})*length(IEI_downdist{inputdown})^(-1/3);
nbinsup = round((max(IEI_updist{inputup}) - min(IEI_updist{inputup}))/binsizeup);
nbinsdown = round((max(IEI_downdist{inputdown}) - min(IEI_downdist{inputdown}))/binsizedown);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 | isnan(nbinsup)
        nbinsup=2;
    end    
end
if nbinsdown > 1e4
    nbinsdown = 1e4;
else if isempty(nbinsdown) ==1 | isnan(nbinsdown)
        nbinsdown=2;
    end
end
subplot(5,3,14);hist(IEI_updist{inputup},nbinsup);title('Up Events'); xlabel(num2str(mean(IEI_updist{inputup})));
subplot(5,3,15);hist(IEI_downdist{inputdown},nbinsdown);title('Down Events'); xlabel(num2str(mean(IEI_downdist{inputdown})));

% PLOT BASED UPON P VALUES
mingaussup = find(IEIupgauss_p == min(IEIupgauss_p));mingaussup=mingaussup(1);
maxgaussup = find(IEIupgauss_p == max(IEIupgauss_p));maxgaussup=maxgaussup(1);
mingaussdown = find(IEIdowngauss_p == min(IEIdowngauss_p));mingaussdown=mingaussdown(1);
maxgaussdown = find(IEIdowngauss_p == max(IEIdowngauss_p));maxgaussdown=maxgaussdown(1);
minexpup = find(IEIupexp_p == min(IEIupexp_p));minexpup=minexpup(1);
maxexpup = find(IEIupexp_p == max(IEIupexp_p));maxexpup=maxexpup(1);
minexpup = find(IEIdownexp_p == min(IEIdownexp_p));minexpup=minexpup(1);
maxexpup = find(IEIdownexp_p == max(IEIdownexp_p));maxexpup=maxexpup(1);

figure(3);
subplot(8,2,1); plot(time,Xvec);hold on; ylabel('displacement');plot(time(up_ind(:,mingaussup)==1),Xvec(up_ind(:,mingaussup)==1),'g.');title(sprintf('%s %s%s %s%s %s%s %s%s','Spikes ','nup= ',num2str(nup(mingaussup)),'upp= ',num2str(IEIupgauss_p(mingaussup)),'eventsup= ',num2str(up_events(1,mingaussup))));
binsizeup = 2*iqr(IEI_updist{mingaussup})*length(IEI_updist{mingaussup})^(-1/3);      % freedman-diaconis rule
nbinsup = round((max(IEI_updist{mingaussup}) - min(IEI_updist{mingaussup}))/binsizeup);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 | isnan(nbinsup)
        nbinsup=2;
    end    
end
subplot(8,2,2);hist(IEI_updist{mingaussup},nbinsup);
subplot(8,2,3); plot(time,Xvec);hold on; ylabel('displacement');plot(time(up_ind(:,maxgaussup)==1),Xvec(up_ind(:,maxgaussup)==1),'g.');title(sprintf('%s %s%s %s%s %s%s %s%s','Spikes ','nup= ',num2str(nup(maxgaussup)),'upp= ',num2str(IEIupgauss_p(maxgaussup)),'eventsup= ',num2str(up_events(1,maxgaussup))));
binsizeup = 2*iqr(IEI_updist{maxgaussup})*length(IEI_updist{maxgaussup})^(-1/3);      % freedman-diaconis rule
nbinsup = round((max(IEI_updist{maxgaussup}) - min(IEI_updist{maxgaussup}))/binsizeup);
if nbinsup > 1e4                    % maximum limit on bins to reduce memory use
    nbinsup = 1e4;
else if isempty(nbinsup) ==1 | isnan(nbinsup)
        nbinsup=2;
    end    
end
subplot(8,2,4);hist(IEI_updist{maxgaussup},nbinsup);

%%%%%%% CONTINUE THIS PATTERN FOR ABOVE

display('Saving...');
%save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp7-kp2-start1end2e5-xfish1.0Noise-movingavg_concwins2.mat');
%save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/noisysinewave22movingavg_concwins2_exp.mat');

display('Finished.');
