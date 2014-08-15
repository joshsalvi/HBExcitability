% Loop through different filters.
clear all; close all;

display('Starting...');
% input MAT file (CHOOSE)
dir_in = '/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-02-21.01/Ear 1/Cell 4/20140221-cell4.mat';
% output directory and prefix (CHOOSE)
dir_out = '/Users/joshsalvi/Downloads/output/20140221-cell4';

% input operating point
Fp = 1;         % force index (CHOOSE)
kp = 1;         % stiffness index (CHOOSE)

winloop = 1;         % loop through windows? (1=yes, 0=no) (CHOOSE)
win1    = 500e-3;    % window (sec) if winloop=0; (CHOOSE)

if winloop==1
    winloopmin = 100e-3;     % minimum window for looping (CHOOSE)
    winloopmax = 5000e-3;   % maximum window for looping (CHOOSE)
    Nwinloop = 10;       % number of window loops (CHOOSE)
    winlooprange = winloopmin:(winloopmax-winloopmin)/(Nwinloop-1):winloopmax;
else
    winlooprange = win1;
end

Fs         = 10e3;     % sampling rate (Hz) (CHOOSE)

freq1min   = 50;     % starting frequency of high-pass filter (CHOOSE)
freq1max   = 1000;    % ending frequency of high-pass filter (CHOOSE)
Nfreq1     = 20;      % number of frequencies to analyze (CHOOSE)
freq1range = freq1min:(freq1max-freq1min)/(Nfreq1-1):freq1max;

mintimeyn = 1;      % loop through minimum times? (CHOOSE)
if mintimeyn == 1
    minTmin   = 1e-3;    % minimum residence time allowed, min (CHOOSE)
    minTmax   = 100e-3;  % minimum residence time allowed, max (CHOOSE)
    NminT     = 15;      % number of minimum times to loop through (CHOOSE)
    minTrange = minTmin:(minTmax-minTmin)/(NminT-1):minTmax;
else
    minT      = 1e-20;   % very small minT if not otherwise described
end


savesig  = 0;       % save the signal? (1=yes; 0=no) (CHOOSE)
outfile = 'xyz';    % initialize outfile    

altsol  = 1;        % run alternative solutions? (1=yes; 0=no) (CHOOSE)
figures = 0;        % generate figures (0=no RECOMMENDED) (CHOOSE)

% Display items
if winloop == 1
    disp('Window Loop: YES');
else
    disp('Window Loop: NO');
end
if mintimeyn == 1
    disp('Use Minimum RT: YES');
else
    disp('Use Minimum RT: NO');
end
if altsol == 1
    disp('Alternative Solution: YES');
else
    disp('Alternative Solution: NO');
end
if savesig == 1
    disp('Save: YES');
else
    disp('Save: NO');
end


display('Initial import...');
load(dir_in);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   LOOPS   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Looping through windows...');
save('a.mat','k_rand','F_rand','Xd');

for a = 1:length(freq1range)
    for b = 1:length(winlooprange)
        for c = 1:length(minTrange)
        clear outfile
        iind = num2str(a);jind = num2str(b);kind = num2str(c);
        outfile = sprintf('%s%s%s%s%s%s%s%s%s%s%s%s',dir_out,'-HPFilt',num2str(freq1range(a)),'-DetrendWin',num2str(winlooprange(b)),'-MinimumResidenceTime',num2str(minTrange(c)),'-Fs',num2str(Fs),'-',[iind jind kind],'.mat');
        [uptimes{a,b,c} downtimes{a,b,c}] = RTfilt('a.mat',Fp,kp,winlooprange(b),freq1range(a),Fs,savesig,outfile,minTrange(c),altsol,figures);
        end
    end
if a == round(length(freq1range)/2)
    display('...50% complete...')
end
end

disp('Saving number of events...');
for a = 1:length(freq1range)
    for b = 1:length(winlooprange)
        for c = 1:length(minTrange)
            upevents(a,b,c) = length(uptimes{a,b,c});
            downevents(a,b,c) = length(downtimes{a,b,c});
        end
    end
end


display('Saving...');
save(sprintf('%s%s',dir_out,'UpDownTimes.mat'),'uptimes','downtimes','dir_in','Fs','winlooprange','freq1range','minTrange','altsol','upevents','downevents');
display('Finished.');

%%  Plot the Data
% If figures==1 above, then other plots will be generated. This is highly
% recommended to avoid a large number of plots from being generated.

figure(1);
subplot(3,3,1); title('High-Pass Filter Frequency');
clear m m1 n n1;m=mean(mean(upevents,3),2);m1(:,:)=m;n=mean(mean(downevents,3),2);n1(:,:)=n;plot(freq1range,m1,'r');hold on;plot(freq1range,n1,'b');ylabel('Mean Events (up=r;do=b)');xlabel('High-Pass Filter Frequency');
[fitupF rupF]=fit(freq1range',m1,'exp1');[fitdownF rdownF]=fit(freq1range',m1,'exp1'); 
subplot(3,3,4);plot(freq1range(2:end),diff(m1),'r');hold on;plot(freq1range(2:end),diff(n1),'b');ylabel('Diff(events)');title('High-Pass Filter Frequency');
subplot(3,3,7);plot(freq1range(3:end),diff(diff(m1)),'r');hold on;plot(freq1range(3:end),diff(diff(n1)),'b');ylabel('Diff(Diff(events))');title('High-Pass Filter Frequency');
clear upFstd downFstd; upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
subplot(3,3,1);hold on;plot(freq1range,upFstd,'r.');plot(freq1range,downFstd,'b.');
subplot(3,3,4);hold on;plot(freq1range(2:end),diff(upFstd),'r.');plot(freq1range(2:end),diff(downFstd),'b.');
subplot(3,3,7);hold on;plot(freq1range(3:end),diff(diff(upFstd)),'r.');plot(freq1range(3:end),diff(diff(downFstd)),'b.');
% Find the appropriate detrending window from the point at which there is
% the least amount of change in the number of events with a change in the
% window size.
u = find(upFstd==min(upFstd)); d = find(upFstd==min(upFstd)); u = u(1); d = d(1);
upF = freq1range(u); downF = freq1range(d);
subplot(3,3,1);hold on;scatter(upF,m1(u),'ro');scatter(downF,n1(d),'ro');

subplot(3,3,2);title('Detrending Window Size');
clear m m1 n n1;m=mean(mean(upevents,3),1);m1(:,:)=m;n=mean(mean(downevents,3),1);n1(:,:)=n;plot(winlooprange,m1,'r');hold on;plot(winlooprange,n1,'b');ylabel('Mean Events (up=r;do=b)');xlabel('Detrending Window Size');
[fitupW rupW]=fit(winlooprange',m1','exp1');[fitdownW rdownW]=fit(winlooprange',m1','exp1');
subplot(3,3,5);plot(winlooprange(2:end),diff(m1),'r');hold on;plot(winlooprange(2:end),diff(n1),'b');ylabel('Diff(events)');xlabel('Detrending Window Size');
subplot(3,3,8);plot(winlooprange(3:end),diff(diff(m1)),'r');hold on;plot(winlooprange(3:end),diff(diff(n1)),'b');ylabel('Diff(Diff(events))');xlabel('Detrending Window Size');
clear upFstd downFstd;upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
subplot(3,3,2);hold on;plot(winlooprange,upFstd,'r.');plot(winlooprange,downFstd,'b.');
subplot(3,3,5);hold on;plot(winlooprange(2:end),diff(upFstd),'r.');plot(winlooprange(2:end),diff(downFstd),'b.');
subplot(3,3,8);hold on;plot(winlooprange(3:end),diff(diff(upFstd)),'r.');plot(winlooprange(3:end),diff(diff(downFstd)),'b.');
% Find the appropriate detrending window from the point at which there is
% the least amount of change in the number of events with a change in the
% window size.
u = findnearest(diff(m1),0)+1; d = findnearest(diff(n1),0)+1; u = u(1); d=d(1);
upW = winlooprange(u); downW = winlooprange(d);
subplot(3,3,2);hold on;scatter(upW,m1(u),'ro');scatter(downW,n1(d),'ro');

subplot(3,3,3);title('Minimum RT Time');
clear m m1 n n1;m=(mean(mean(upevents,2),1));m1(:,:)=m;n=(mean(mean(downevents,2),1));n1(:,:)=n;plot(minTrange,m1,'r');hold on;plot(minTrange,n1,'b');ylabel('Mean Events (up=r;do=b)');xlabel('Minimum Residence Time');
[fitupT rupT]=fit(minTrange',m1,'exp1');[fitdownT rdownT]=fit(minTrange',m1,'exp1');
subplot(3,3,6);plot(minTrange(2:end),diff(m1),'r');hold on;plot(minTrange(2:end),diff(n1),'b');ylabel('Diff(events)');xlabel('Minimum Residence Time');
subplot(3,3,9);plot(minTrange(3:end),diff(diff(m1)),'r');hold on;plot(minTrange(3:end),diff(diff(n1)),'b');ylabel('Diff(Diff(events))');xlabel('Minimum Residence Time');
clear upFstd downFstd;upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
subplot(3,3,3);hold on;plot(minTrange,upFstd,'r.');plot(minTrange,downFstd,'b.');
subplot(3,3,6);hold on;plot(minTrange(2:end),diff(upFstd),'r.');plot(minTrange(2:end),diff(downFstd),'b.');
subplot(3,3,9);hold on;plot(minTrange(3:end),diff(diff(upFstd)),'r.');plot(minTrange(3:end),diff(diff(downFstd)),'b.');
% Find the appropriate minimum residence time from the point with the
% maximum slope.
u = find(diff(diff(m1))==max(diff(diff(m1))))+1; d = find(diff(diff(n1))==max(diff(diff(n1))))+1; u = u(1); d = d(1);
upT = minTrange(u); downT = minTrange(d);
subplot(3,3,3);hold on;scatter(upT,m1(u),'ro');scatter(downT,n1(d),'ro');

% Run RTfilt with the 
outfile = sprintf('%s%s%s%s%s%s%s%s%s%s%s%s',dir_out,'-HPFilt',num2str(mean([upF downF])),'-DetrendWin',num2str(mean([upW downW])),'-MinimumResidenceTime',num2str(mean([upT downT])),'-Fs',num2str(Fs),'-','POSTPLOTS','.mat');
RTfilt('a.mat',Fp,kp,mean([upW downW]),mean([upF downF]),Fs,1,outfile,mean([upT downT]),1,1);

%%  Statistics on Filtered Data 
% If savesig==1 above, then all of this statistics should have been saved
% to corresponding MAT files.

% Pick the index you would like to analyze
i = 1;
j = 1;
k = 1;