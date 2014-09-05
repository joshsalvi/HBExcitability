% Loop through different filters.
clear all; close all;

display('Starting...');

% input MAT file (CHOOSE)
dir_in = '/Users/joshsalvi/Downloads/output/artificialwaveforms/waves.mat';
%dir_in = '/Users/joshsalvi/Documents/Lab/Lab/Data Analysis/State Spaces/20130419/20130419-cell8.mat';
%dir_in = '/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/noisysinewave.mat';
%dir_in = '/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/xfish1.0Noise.mat';

% output directory and prefix (CHOOSE)
dir_out = '/Users/joshsalvi/Downloads/output/artificialwaveforms/waves-Fp7-kp1-';

% input operating point
Fp = 7;        % force index (CHOOSE)
kp = 1;        % stiffness index (CHOOSE)

winloop = 1;         % loop through windows? (1=yes, 0=no) (CHOOSE)
win1    = 500e-3;    % window (sec) if winloop=0; (CHOOSE)

if winloop==1
    winloopmin = 50e-3;     % minimum window for looping (CHOOSE)
    winloopmax = 3000e-3;   % maximum window for looping (CHOOSE)
    Nwinloop = 10;       % number of window loops (CHOOSE)
    winlooprange = winloopmin:(winloopmax-winloopmin)/(Nwinloop-1):winloopmax;
else
    winlooprange = win1;
end

Fs         = 10e3;     % sampling rate (Hz) (CHOOSE)

freq1min   = 50;     % starting frequency of high-pass filter (CHOOSE)
freq1max   = 2000;    % ending frequency of high-pass filter (CHOOSE)
Nfreq1     = 10;      % number of frequencies to analyze (CHOOSE)
freq1range = freq1min:(freq1max-freq1min)/(Nfreq1-1):freq1max;

mintimeyn = 1;      % loop through minimum times? (CHOOSE)
if mintimeyn == 1
    minTmin   = 1e-3;    % minimum residence time allowed, min (CHOOSE)
    minTmax   = 200e-3;  % minimum residence time allowed, max (CHOOSE)
    NminT     = 5;      % number of minimum times to loop through (CHOOSE)
    minTrange = minTmin:(minTmax-minTmin)/(NminT-1):minTmax;
else
    minT      = 1e-20;   % very small minT if not otherwise described
end


savesig = 0;       % save the signal? (1=yes; 0=no) (CHOOSE)
outfile = 'xyz';    % initialize outfile    

altsol  = 1;        % run alternative solutions? (1=yes; 0=no) (CHOOSE; default=1)
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
        [uptimes{a,b,c}, downtimes{a,b,c}, Xvecfilt_detrended{a,b,c},RTupexp_p(a,b,c),RTupgauss_p{a,b,c},RTdownexp_p(a,b,c),RTdowngauss_p{a,b,c},RTupexp_stat(a,b,c),RTdownexp_stat(a,b,c),RTupgauss_stat{a,b,c},RTdowngauss_stat{a,b,c},poisswinmov_up{a,b,c},poisswinmov_down{a,b,c},movwinrange{a,b,c}] = RTfilt('a.mat',Fp,kp,winlooprange(b),freq1range(a),Fs,savesig,outfile,minTrange(c),altsol,figures,0,[]);
        end
    end
if a == round(length(freq1range)/2)
    display('...50% complete...')
end
end

%Initialize
RTdowngauss_p2 = zeros(3,length(freq1range),length(winlooprange),length(minTrange));
RTupgauss_p2 = zeros(3,length(freq1range),length(winlooprange),length(minTrange));
RTupgauss_stat2 = zeros(3,length(freq1range),length(winlooprange),length(minTrange));
RTdowngauss_stat2 = zeros(3,length(freq1range),length(winlooprange),length(minTrange));

disp('Saving number of events...');
for a = 1:length(freq1range)
    for b = 1:length(winlooprange)
        for c = 1:length(minTrange)
            upevents(a,b,c) = length(uptimes{a,b,c});
            downevents(a,b,c) = length(downtimes{a,b,c});
            sigfiltstd(a,b,c) = std(Xvecfilt_detrended{a,b,c});
            RTupgauss_p2(:,a,b,c) = RTupgauss_p{a,b,c};
            RTdowngauss_p2(:,a,b,c) = RTdowngauss_p{a,b,c};
            RTupgauss_stat2(:,a,b,c) = RTupgauss_stat{a,b,c};
            RTdowngauss_stat2(:,a,b,c) = RTdowngauss_stat{a,b,c};
        end
    end
end


display('Saving...');
save(sprintf('%s%s',dir_out,'UpDownTimes.mat'),'uptimes','downtimes','dir_in','Fs','winlooprange','freq1range','minTrange','altsol','upevents','downevents','dir_out','Fp','kp','sigfiltstd','RTupexp_p','RTupgauss_p2','RTdownexp_p','RTdowngauss_p2','RTupexp_stat','RTdownexp_stat','RTupgauss_stat2','RTdowngauss_stat2','poisswinmov_up','poisswinmov_down','movwinrange');
display('Finished.');

%% Pre-load (RUN FIRST - just once per cell)
clear all; close all;
load('/Users/joshsalvi/Documents/Lab/Lab/Data Analysis/State Spaces/20130419/20130419-cell7.mat')
save('a.mat')
clear all;
%% Option(1) Plot the data and run RTfilt based on set parameters
close all;
% (1) Detrending window size: Chosen from the window that reduces the
% standard deviation of the signal.
% (2) High-pass filter frequency: Chosen first from the set from (1) and
% the frequency that yields the smallest change in the number of events
% (found from the std() of the number of events versus filter frequency).
% (3) Minimum RT time: Currently set to find the minimum diff(diff()) of
% the number of events with a change in the minimum window size. However,
% the code also plots the test statistics for each minimum time in case it
% is desired to use these as cutoffs.
runrt = 1;      % Run RTfilt and save?
offsetyn = 1;   % apply an offset?    %% USE IN CASES WHERE COUNTS BEGIN TO INCREASE AGAIN BY MORE THAN TWOFOLD 
offset = 0;

% load the correct files
load('/Users/joshsalvi/Downloads/output/olddata/20130419-cell7-Fp19-kp1-UpDownTimes.mat');

warning off;
% If figures==1 above, then other plots will be generated. This is highly
% recommended to avoid a large number of plots from being generated.
upevents2=upevents;
downevents2=downevents;
RTupexp_pn = RTupexp_p;
RTdownexp_pn = RTdownexp_p;
RTupgauss_p2n = RTupgauss_p2;
RTdowngauss_p2n = RTdowngauss_p2;
clear RTupexp_p RTdownexp_p RTupgauss_p2 RTdowngauss_p2 RTupgauss_p RTdowngauss_p;

figure;
subplot(3,3,1);xlabel('Detrending Window Size');
clear m m1 n n1;m=mean(mean(sigfiltstd,3),1);m1(:,:)=m;n=mean(mean(sigfiltstd,3),1);n1(:,:)=n;plot(winlooprange,m1,'r');hold on;plot(winlooprange,n1,'b');ylabel('RMS of Signal');xlabel('Detrending Window Size');
[~, rupW]=fit(winlooprange',m1','exp1');[~, rdownW]=fit(winlooprange',m1','exp1');
subplot(3,3,4);plot(winlooprange(2:end),diff(m1),'r');hold on;plot(winlooprange(2:end),diff(n1),'b');ylabel('Diff(RMS)');xlabel('Detrending Window Size');
subplot(3,3,7);plot(winlooprange(3:end),diff(diff(m1)),'r');hold on;plot(winlooprange(3:end),diff(diff(n1)),'b');ylabel('Diff(Diff(RMS))');xlabel('Detrending Window Size');
clear upFstd downFstd;upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
subplot(3,3,1);hold on;plot(winlooprange,upFstd,'r.');plot(winlooprange,downFstd,'b.');
subplot(3,3,4);hold on;plot(winlooprange(2:end),diff(upFstd),'r.');plot(winlooprange(2:end),diff(downFstd),'b.');
subplot(3,3,7);hold on;plot(winlooprange(3:end),diff(diff(upFstd)),'r.');plot(winlooprange(3:end),diff(diff(downFstd)),'b.');
% Find the appropriate detrending window from the point at which there is
% the least amount of change in the number of events with a change in the
% window size.
uw = find(m1==min(m1)); dw = find(n1==min(n1)); uw = uw(1); dw=dw(1);
upW = winlooprange(uw); downW = winlooprange(dw);
subplot(3,3,1);hold on;scatter(upW,m1(uw),'ro');scatter(downW,n1(dw),'ro');
% Placeholder for selected part of up/down events
clear upevents downevents;
upevents = upevents2(:,uw,:); downevents = downevents2(:,dw,:);

subplot(3,3,2); xlabel('High-Pass Filter Frequency');
clear m m1 n n1;m=mean(upevents,3);m1(:,:)=m;n=mean(downevents,3);n1(:,:)=n;plot(freq1range,m1,'r');hold on;plot(freq1range,n1,'b');ylabel('Mean Events (up=r;do=b)');xlabel('High-Pass Filter Frequency');
[~, rupF]=fit(freq1range',m1,'exp1');[~, rdownF]=fit(freq1range',m1,'exp1'); 
subplot(3,3,5);plot(freq1range(2:end),diff(m1),'r');hold on;plot(freq1range(2:end),diff(n1),'b');ylabel('Diff(events)');xlabel('High-Pass Filter Frequency');
subplot(3,3,8);plot(freq1range(3:end),diff(diff(m1)),'r');hold on;plot(freq1range(3:end),diff(diff(n1)),'b');ylabel('Diff(Diff(events))');xlabel('High-Pass Filter Frequency');
clear upFstd downFstd; upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
subplot(3,3,2);hold on;plot(freq1range,upFstd,'r.');plot(freq1range,downFstd,'b.');
subplot(3,3,5);hold on;plot(freq1range(2:end),diff(upFstd),'r.');plot(freq1range(2:end),diff(downFstd),'b.');
subplot(3,3,8);hold on;plot(freq1range(3:end),diff(diff(upFstd)),'r.');plot(freq1range(3:end),diff(diff(downFstd)),'b.');
% Find the appropriate detrending window from the point at which there is
% the least amount of change in the number of events with a change in the
% window size.
uf = find(upFstd==min(upFstd)); df = find(upFstd==min(upFstd)); uf = uf(1); df = df(1);
upF = freq1range(uf); downF = freq1range(df);
subplot(3,3,2);hold on;scatter(upF,m1(uf),'ro');scatter(downF,n1(df),'ro');
% Placeholder for selected part of up/down events
clear upevents downevents;
upevents = upevents2(uf,uw,:); downevents = downevents2(df,dw,:);
RTupexp_p = RTupexp_pn(uf,uw,:); RTdownexp_p = RTdownexp_pn(uf,uw,:);
RTupgauss_p = RTupgauss_p2n(:,uf,uw,:); RTdowngauss_p = RTdowngauss_p2n(:,uf,uw,:); 

subplot(3,3,3);xlabel('Minimum RT Time');
clear m m1 n n1;m=upevents;m1(:,:)=m;n=downevents;n1(:,:)=n;plot(minTrange,m1,'r');hold on;plot(minTrange,n1,'b');ylabel('Mean Events (up=r;do=b)');xlabel('Minimum Residence Time');
[~, rupT]=fit(minTrange',m1,'exp1');[~, rdownT]=fit(minTrange',m1,'exp1');
subplot(3,3,6);plot(minTrange(2:end),diff(m1),'r');hold on;plot(minTrange(2:end),diff(n1),'b');ylabel('Diff(events)');xlabel('Minimum Residence Time');
subplot(3,3,9);plot(minTrange(3:end),diff(diff(m1)),'r');hold on;plot(minTrange(3:end),diff(diff(n1)),'b');ylabel('Diff(Diff(events))');xlabel('Minimum Residence Time');
clear upFstd downFstd;upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
subplot(3,3,3);hold on;plot(minTrange,upFstd,'r.');plot(minTrange,downFstd,'b.');
subplot(3,3,6);hold on;plot(minTrange(2:end),diff(upFstd),'r.');plot(minTrange(2:end),diff(downFstd),'b.');
subplot(3,3,9);hold on;plot(minTrange(3:end),diff(diff(upFstd)),'r.');plot(minTrange(3:end),diff(diff(downFstd)),'b.');
% Find the appropriate minimum residence time from the point with the
% maximum slope.
ur = find(diff(diff(m1))==max(diff(diff(m1))))+1; dr = find(diff(diff(n1))==max(diff(diff(n1))))+1; ur = ur(1); dr = dr(1);
upT = minTrange(ur); downT = minTrange(dr);
subplot(3,3,3);hold on;scatter(upT,m1(ur),'ro');scatter(downT,n1(dr),'ro');

% Plot the test statistics for different minimum residence times
figure; subplot(4,2,1);plot(minTrange,squeeze(RTupexp_stat(uf,uw,:)));xlabel('Minimum Residence Time');ylabel('Lillie Test; Exp; UP');
subplot(4,2,2);plot(minTrange,squeeze(RTdownexp_stat(uf,uw,:)));xlabel('Minimum Residence Time');ylabel('Lillie Test; Exp; DOWN');
subplot(4,2,3);plot(minTrange,squeeze(RTupgauss_stat2(1,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('Lillie Test; Gauss; UP');
subplot(4,2,4);plot(minTrange,squeeze(RTdowngauss_stat2(1,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('Lillie Test; Gauss; DOWN');
subplot(4,2,5);plot(minTrange,squeeze(RTupgauss_stat2(2,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('KS Test; Gauss; UP');
subplot(4,2,6);plot(minTrange,squeeze(RTdowngauss_stat2(2,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('KS Test; Gauss; DOWN');
subplot(4,2,7);plot(minTrange,squeeze(RTupgauss_stat2(3,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('JB Test; Gauss; UP');
subplot(4,2,8);plot(minTrange,squeeze(RTdowngauss_stat2(3,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('JB Test; Gauss; DOWN');

% Run RTfilt with the chosen parameters...
if runrt==1
    disp('RTfilt() and Save: ON');
outfile = sprintf('%s%s%s%s%s%s%s%s%s%s%s%s',dir_out,'-HPFilt',num2str(mean([upF downF])),'-DetrendWin',num2str(mean([upW downW])),'-MinimumResidenceTime',num2str(mean([upT downT])),'-Fs',num2str(Fs),'-','POSTPLOTS','.mat');
RTfilt('a.mat',Fp,kp,mean([upW downW]),mean([upF downF]),Fs,1,outfile,mean([upT downT]),1,1,offsetyn,offset);
else
    disp('RTfilt() and Save: OFF');
end
disp('Finished.');

disp(sprintf('%s%s\n%s%s\n%s%s','win= ',num2str(mean([upW downW])),'freq= ',num2str(mean([upF downF])),'minT= ',num2str(mean([upT downT]))))



%% Option(2) Analyze data via LOOPING

clear all; close all;
disp('Loading...');
%load('/Users/joshsalvi/Documents/Lab/Lab/Data Analysis/State Spaces/20130419/20130419-cell7.mat')
load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2014-02-21.01/Ear 1/Cell 5/20140221-cell5.mat')
%load('/Users/joshsalvi/Downloads/output/artificialwaveforms/waves.mat');
save('a.mat')
clear all;

% (1) Detrending window size: Chosen from the window that reduces the
% standard deviation of the signal.
% (2) High-pass filter frequency: Chosen first from the set from (1) and
% the frequency that yields the smallest change in the number of events
% (found from the std() of the number of events versus filter frequency).
% (3) Minimum RT time: Currently set to find the minimum diff(diff()) of
% the number of events with a change in the minimum window size. However,
% the code also plots the test statistics for each minimum time in case it
% is desired to use these as cutoffs.

% initialize
N_F = 15;
N_k = 1;

disp('Looping...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   LOOPS   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:N_F                   
    for j = 1:N_k               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   INPUT   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% load the correct files
%file = '/Users/joshsalvi/Downloads/output/olddata/20130419-cell7-Fp';
file = '/Users/joshsalvi/Downloads/output/20140221-cell5-Fp';
%file = '/Users/joshsalvi/Downloads/output/artificialwaveforms/waves-Fp';
runrt = 1;      % Run RTfilt and save?

% constrain values? (0=no)
constrain1 = 1;
% if you choose to constrain, what values would you like to use?
win_in  = 0.05;
freq_in = 2000;         % taken from individual examples in previous analysis (mean,mode,etc)
minT_in = 0.0008;        % chosen from minimum values seen in individual histograms prior to implementation of a minimum time
% plots? (1=yes)
plotyn = 0;
offsetyn = 1;   % apply an offset?    %% USE IN CASES WHERE COUNTS BEGIN TO INCREASE AGAIN BY MORE THAN TWOFOLD 
offset = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fp = i; kp = j;
load(sprintf('%s%s%s%s%s',file,num2str(Fp),'-kp',num2str(kp),'-UpDownTimes.mat'));
Fp = i; kp = j;
warning off;
% If figures==1 above, then other plots will be generated. This is highly
% recommended to avoid a large number of plots from being generated.
upevents2=upevents;
downevents2=downevents;
RTupexp_pn = RTupexp_p;
RTdownexp_pn = RTdownexp_p;
RTupgauss_p2n = RTupgauss_p2;
RTdowngauss_p2n = RTdowngauss_p2;
clear RTupexp_p RTdownexp_p RTupgauss_p2 RTdowngauss_p2 RTupgauss_p RTdowngauss_p;

if plotyn ==1
figure;
subplot(3,3,1);xlabel('Detrending Window Size');
clear m m1 n n1;m=mean(mean(sigfiltstd,3),1);m1(:,:)=m;n=mean(mean(sigfiltstd,3),1);n1(:,:)=n;plot(winlooprange,m1,'r');hold on;plot(winlooprange,n1,'b');ylabel('RMS of Signal');xlabel('Detrending Window Size');
[~, rupW]=fit(winlooprange',m1','exp1');[~, rdownW]=fit(winlooprange',m1','exp1');
subplot(3,3,4);plot(winlooprange(2:end),diff(m1),'r');hold on;plot(winlooprange(2:end),diff(n1),'b');ylabel('Diff(RMS)');xlabel('Detrending Window Size');
subplot(3,3,7);plot(winlooprange(3:end),diff(diff(m1)),'r');hold on;plot(winlooprange(3:end),diff(diff(n1)),'b');ylabel('Diff(Diff(RMS))');xlabel('Detrending Window Size');
clear upFstd downFstd;upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
subplot(3,3,1);hold on;plot(winlooprange,upFstd,'r.');plot(winlooprange,downFstd,'b.');
subplot(3,3,4);hold on;plot(winlooprange(2:end),diff(upFstd),'r.');plot(winlooprange(2:end),diff(downFstd),'b.');
subplot(3,3,7);hold on;plot(winlooprange(3:end),diff(diff(upFstd)),'r.');plot(winlooprange(3:end),diff(diff(downFstd)),'b.');
% Find the appropriate detrending window from the point at which there is
% the least amount of change in the number of events with a change in the
% window size.
uw = find(m1==min(m1)); dw = find(n1==min(n1)); uw = uw(1); dw=dw(1);
upW = winlooprange(uw); downW = winlooprange(dw);
subplot(3,3,1);hold on;scatter(upW,m1(uw),'ro');scatter(downW,n1(dw),'ro');
% Placeholder for selected part of up/down events
clear upevents downevents;
upevents = upevents2(:,uw,:); downevents = downevents2(:,dw,:);

subplot(3,3,2); xlabel('High-Pass Filter Frequency');
clear m m1 n n1;m=mean(upevents,3);m1(:,:)=m;n=mean(downevents,3);n1(:,:)=n;plot(freq1range,m1,'r');hold on;plot(freq1range,n1,'b');ylabel('Mean Events (up=r;do=b)');xlabel('High-Pass Filter Frequency');
[~, rupF]=fit(freq1range',m1,'exp1');[~, rdownF]=fit(freq1range',m1,'exp1'); 
subplot(3,3,5);plot(freq1range(2:end),diff(m1),'r');hold on;plot(freq1range(2:end),diff(n1),'b');ylabel('Diff(events)');xlabel('High-Pass Filter Frequency');
subplot(3,3,8);plot(freq1range(3:end),diff(diff(m1)),'r');hold on;plot(freq1range(3:end),diff(diff(n1)),'b');ylabel('Diff(Diff(events))');xlabel('High-Pass Filter Frequency');
clear upFstd downFstd; upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
subplot(3,3,2);hold on;plot(freq1range,upFstd,'r.');plot(freq1range,downFstd,'b.');
subplot(3,3,5);hold on;plot(freq1range(2:end),diff(upFstd),'r.');plot(freq1range(2:end),diff(downFstd),'b.');
subplot(3,3,8);hold on;plot(freq1range(3:end),diff(diff(upFstd)),'r.');plot(freq1range(3:end),diff(diff(downFstd)),'b.');
% Find the appropriate detrending window from the point at which there is
% the least amount of change in the number of events with a change in the
% window size.
uf = find(upFstd==min(upFstd)); df = find(upFstd==min(upFstd)); uf = uf(1); df = df(1);
upF = freq1range(uf); downF = freq1range(df);
subplot(3,3,2);hold on;scatter(upF,m1(uf),'ro');scatter(downF,n1(df),'ro');
% Placeholder for selected part of up/down events
clear upevents downevents;
upevents = upevents2(uf,uw,:); downevents = downevents2(df,dw,:);
RTupexp_p = RTupexp_pn(uf,uw,:); RTdownexp_p = RTdownexp_pn(uf,uw,:);
RTupgauss_p = RTupgauss_p2n(:,uf,uw,:); RTdowngauss_p = RTdowngauss_p2n(:,uf,uw,:); 

subplot(3,3,3);xlabel('Minimum RT Time');
clear m m1 n n1;m=upevents;m1(:,:)=m;n=downevents;n1(:,:)=n;plot(minTrange,m1,'r');hold on;plot(minTrange,n1,'b');ylabel('Mean Events (up=r;do=b)');xlabel('Minimum Residence Time');
[~, rupT]=fit(minTrange',m1,'exp1');[~, rdownT]=fit(minTrange',m1,'exp1');
subplot(3,3,6);plot(minTrange(2:end),diff(m1),'r');hold on;plot(minTrange(2:end),diff(n1),'b');ylabel('Diff(events)');xlabel('Minimum Residence Time');
subplot(3,3,9);plot(minTrange(3:end),diff(diff(m1)),'r');hold on;plot(minTrange(3:end),diff(diff(n1)),'b');ylabel('Diff(Diff(events))');xlabel('Minimum Residence Time');
clear upFstd downFstd;upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
subplot(3,3,3);hold on;plot(minTrange,upFstd,'r.');plot(minTrange,downFstd,'b.');
subplot(3,3,6);hold on;plot(minTrange(2:end),diff(upFstd),'r.');plot(minTrange(2:end),diff(downFstd),'b.');
subplot(3,3,9);hold on;plot(minTrange(3:end),diff(diff(upFstd)),'r.');plot(minTrange(3:end),diff(diff(downFstd)),'b.');
% Find the appropriate minimum residence time from the point with the
% maximum slope.
ur = find(diff(diff(m1))==max(diff(diff(m1))))+1; dr = find(diff(diff(n1))==max(diff(diff(n1))))+1; ur = ur(1); dr = dr(1);
upT = minTrange(ur); downT = minTrange(dr);
subplot(3,3,3);hold on;scatter(upT,m1(ur),'ro');scatter(downT,n1(dr),'ro');

% Plot the test statistics for different minimum residence times
figure; subplot(4,2,1);plot(minTrange,squeeze(RTupexp_stat(uf,uw,:)));xlabel('Minimum Residence Time');ylabel('Lillie Test; Exp; UP');
subplot(4,2,2);plot(minTrange,squeeze(RTdownexp_stat(uf,uw,:)));xlabel('Minimum Residence Time');ylabel('Lillie Test; Exp; DOWN');
subplot(4,2,3);plot(minTrange,squeeze(RTupgauss_stat2(1,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('Lillie Test; Gauss; UP');
subplot(4,2,4);plot(minTrange,squeeze(RTdowngauss_stat2(1,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('Lillie Test; Gauss; DOWN');
subplot(4,2,5);plot(minTrange,squeeze(RTupgauss_stat2(2,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('KS Test; Gauss; UP');
subplot(4,2,6);plot(minTrange,squeeze(RTdowngauss_stat2(2,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('KS Test; Gauss; DOWN');
subplot(4,2,7);plot(minTrange,squeeze(RTupgauss_stat2(3,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('JB Test; Gauss; UP');
subplot(4,2,8);plot(minTrange,squeeze(RTdowngauss_stat2(3,uf,uw,:)));xlabel('Minimum Residence Time');ylabel('JB Test; Gauss; DOWN');

else
    clear m m1 n n1;m=mean(mean(sigfiltstd,3),1);m1(:,:)=m;n=mean(mean(sigfiltstd,3),1);n1(:,:)=n;
    [~, rupW]=fit(winlooprange',m1','exp1');[~, rdownW]=fit(winlooprange',m1','exp1');
    uw = find(m1==min(m1)); dw = find(n1==min(n1)); uw = uw(1); dw=dw(1);
    upW = winlooprange(uw); downW = winlooprange(dw);
    % Placeholder for selected part of up/down events
    clear upevents downevents;
    upevents = upevents2(:,uw,:); downevents = downevents2(:,dw,:);
    clear upFstd downFstd;upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
    clear m m1 n n1;m=mean(upevents,3);m1(:,:)=m;n=mean(downevents,3);n1(:,:)=n;
    [~, rupF]=fit(freq1range',m1,'exp1');[~, rdownF]=fit(freq1range',m1,'exp1'); 
    clear upFstd downFstd; upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
    uf = find(upFstd==min(upFstd)); df = find(upFstd==min(upFstd)); uf = uf(1); df = df(1);
    upF = freq1range(uf); downF = freq1range(df);
    % Placeholder for selected part of up/down events
    clear upevents downevents;
    upevents = upevents2(uf,uw,:); downevents = downevents2(df,dw,:);
    RTupexp_p = RTupexp_pn(uf,uw,:); RTdownexp_p = RTdownexp_pn(uf,uw,:);
    RTupgauss_p = RTupgauss_p2n(:,uf,uw,:); RTdowngauss_p = RTdowngauss_p2n(:,uf,uw,:);
    clear m m1 n n1;m=upevents;m1(:,:)=m;n=downevents;n1(:,:)=n;
    clear upFstd downFstd;upFstd=movingstd(m1,2);downFstd=movingstd(n1,2);
    [~, rupT]=fit(minTrange',m1,'exp1');[~, rdownT]=fit(minTrange',m1,'exp1');
    % Find the appropriate minimum residence time from the point with the
    % maximum slope.
    ur = find(diff(diff(m1))==max(diff(diff(m1))))+1; dr = find(diff(diff(n1))==max(diff(diff(n1))))+1; ur = ur(1); dr = dr(1);
    upT = minTrange(ur); downT = minTrange(dr);
end


    
% Run RTfilt with the chosen parameters...

if constrain1 == 0
win_in  = mean([upW downW]);
freq_in = mean([upF downF]);
minT_in = mean([upT downT]);
end

if runrt==1
    %disp('RTfilt() and Save: ON');
outfile = sprintf('%s%s%s%s%s%s%s%s%s%s%s%s',dir_out,'-HPFilt',num2str(freq_in),'-DetrendWin',num2str(win_in),'-MinimumResidenceTime',num2str(minT_in),'-Fs',num2str(Fs),'-','POSTPLOTS','.mat');
[upt{Fp,kp}, downt{Fp,kp}, ~,RTupexp_p5(1,Fp,kp),RTupgauss_p5(1:3,Fp,kp),RTdownexp_p5(1,Fp,kp),RTdowngauss_p5(1:3,Fp,kp),RTupexp_stat5(1,Fp,kp),RTdownexp_stat5(1:3,Fp,kp),RTupgauss_stat5(1:3,Fp,kp),RTdowngauss_stat5(1:3,Fp,kp),poisswinmov_up5{Fp,kp},poisswinmov_down5{Fp,kp},movwinrange5{Fp,kp}]=RTfilt('a.mat',Fp,kp,win_in,freq_in,Fs,1,outfile,minT_in,1,0,offsetyn,offset);
else
    disp('RTfilt() and Save: OFF');
end
%{
if Fp == round(N_F/2)
    disp('...50%...');
end
%}
up_eventsall(Fp,kp)   = length(upt{Fp,kp});
down_eventsall(Fp,kp) = length(downt{Fp,kp});
end
end

disp('Saving...');
save(sprintf('%s%s',file,'-stats_all.mat'),'RTupexp_p5','RTupgauss_p5','RTdownexp_p5','RTdowngauss_p5','RTupexp_stat5','RTdownexp_stat5','RTupgauss_stat5','RTdowngauss_stat5','win_in','freq_in','minT_in','down_eventsall','up_eventsall','poisswinmov_up5','poisswinmov_down5','movwinrange5');
disp('Finished.');


%% Loop again to calculate the following:
% (1) Exponential and Gaussian distribution statistics
% (2) Mean residence times
% (3) Poisson distribution statistics
% (4) Additional work with the up and down times
clear all; close all;

disp('Initializing...');
% LOAD DATA
cd '/Users/joshsalvi/Downloads/output/constrain';
prefix = '20140221-cell5';
%cd '/Users/joshsalvi/Downloads/output/';
%prefix = '20140221-cell4-';
files = dir([prefix '*-POSTPLOTS.mat']);
statfile = dir([prefix '*-stats_all.mat']);
histyn = 0;

disp('Looping through analyzed files...');
for m = 1:length(files)
    load(files(m).name);    % load        
    meanRT_up(m)   = mean(uptimes);     % find mean time spent in the up and down states
    meanRT_down(m) = mean(downtimes);
    upt{m}   = uptimes;
    downt{m} = downtimes;
    % Fit:  normal, exponential, and poisson... RTupdistD{m}(1) will be the
    % distribution RTupdistD{m}(1).DistName is the name of the distribution
    % and RTupdistD{m}(1).Paramci are its fitting parameters
    [RTupdistD{m}, RTupdistPD{m}] = allfitdist(uptimes);
    [RTdowndistD{m}, RTdowndistPD{m}] = allfitdist(downtimes);
    
    % Chi-squared statistic: n.s. if normal distribution
    [~,RTupgausschi2gof_p(m),RTupgausschi2gof_stat{m}]     = chi2gof(uptimes);
    [~,RTdowngausschi2gof_p(m),RTdowngausschi2gof_stat{m}] = chi2gof(downtimes);
    
    clear binsizeup binsizedown
    binsizeup = 2*iqr(uptimes)*length(uptimes)^(-1/3);      % freedman-diaconis rule
    binsizedown = 2*iqr(downtimes)*length(downtimes)^(-1/3);
    nbinsup(m) = round((max(uptimes) - min(uptimes))/binsizeup);
    nbinsdown(m) = round((max(downtimes) - min(downtimes))/binsizedown);
    if isnan(nbinsup(m)) == 1 || isinf(nbinsup(m)) == 1
             nbinsup(m) = 5;
    end
    if isnan(nbinsdown(m)) == 1 || isinf(nbinsdown(m)) == 1
             nbinsdown(m) = 5;
    end
    [upa{m} upb{m}] = hist(uptimes,nbinsup(m));      % make histograms
    [downa{m} downb{m}] = hist(downtimes,nbinsdown(m));
    n1 = sum(upa{m}); n2 = sum(downa{m});
    pdup{m} = fitdist(upb{m}','Exponential','Frequency',upa{m}');
    pddown{m} = fitdist(downb{m}','Exponential','Frequency',downa{m}');
    expcountsup{m} = n1*pdf(pdup{m},upb{m}); expcountsdown{m} = n1*pdf(pddown{m},downb{m});
    % Test for normal and exponential distributions using the chi-squared
    % test
    [~,RTupexpchi2gof_p(m),RTupexpchi2gof_stat{m}] = chi2gof(upb{m},'Ctrs',upb{m}, ...
        'Frequency', upa{m}, ...
        'Expected', expcountsup{m},...
        'Nparams',1);
    [~,RTdownexpchi2gof_p(m),RTdownexpchi2gof_stat{m}] = chi2gof(downb{m},'Ctrs',downb{m}, ...
        'Frequency', downa{m}, ...
        'Expected', expcountsdown{m},...
        'Nparams',1);
    if m == round(length(files)/2)
        disp('...50%...');
    end
end

% Calculate information from stats
disp('Running Statistics...');
load(statfile.name);

spikefreq_up   = up_eventsall./time(end);
spikefreq_down = down_eventsall./time(end);

if histyn ==1
% Preload five histograms
histindex = round(1:(length(files)-1)/(4):length(files));
nbins = 20;

for i = 1:5

figure;
subplot(2,2,1);histfit(upt{histindex(i)},nbins); title('RT up');
subplot(2,2,2);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.5f\nStd. Dev. = %.5f\nMedian = %.5f\nMin = %.5f\nMax = %.5f\nLillie(exp) = %.1e\nLillie(gauss) = %.1e\nKS = %.1e\nJB = %.1e\nChi2(poiss) = %.1e\nNumber of counts = %s',mean(upt{histindex(i)}),std(upt{histindex(i)}),median(upt{histindex(i)}),min(upt{histindex(i)}),max(upt{histindex(i)}),RTupexp_stat5(histindex(i)),RTupgauss_stat5(1,(histindex(i))),RTupgauss_stat5(2,(histindex(i))),RTupgauss_stat5(3,(histindex(i))),(RTuppoiss_pstats.chi2stat),num2str(up_eventsall(histindex(i)))));
subplot(2,2,3);histfit(downt{histindex(i)},nbins); title('RT down');
subplot(2,2,4);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.5f\nStd. Dev. = %.5f\nMedian = %.5f\nMin = %.5f\nMax = %.5f\nLillie(exp) = %.1e\nLillie(gauss) = %.1e\nKS = %.1e\nJB = %.1e\nChi2(poiss) = %.1e\nNumber of counts = %s',mean(downt{histindex(i)}),std(downt{histindex(i)}),median(downt{histindex(i)}),min(downt{histindex(i)}),max(downt{histindex(i)}),RTdownexp_stat5(histindex(i)),RTdowngauss_stat5(1,(histindex(i))),RTdowngauss_stat5(2,(histindex(i))),RTdowngauss_stat5(3,(histindex(i))),(RTdownpoiss_pstats.chi2stat),num2str(down_eventsall(histindex(i)))));
end
end

% Plot spike frequency, meanRT versus a control parameter
% INPUT
forces = 0;
stiffnesses = 0;
force_c = 90e-12;
stiffness_c = 600e-6;

if forces == 1
figure;
plot((force_c-F_rand(1:length(spikefreq_up)))*1e12,spikefreq_up,'r');hold on;plot((force_c-F_rand(1:length(spikefreq_up)))*1e12,spikefreq_down,'b');plot((force_c-F_rand(1:length(spikefreq_up)))*1e12,sqrt((force_c-F_rand(1:length(spikefreq_up)))*1e12),'k--');
title('Red=up, blue=down; dashed=sqrt(Fc-F)');xlabel('Force (pN)');ylabel('Spike Frequency (spikes/sec)');
figure;
plot((force_c-F_rand(1:length(spikefreq_up)))*1e12,meanRT_up,'r');hold on;plot((force_c-F_rand(1:length(spikefreq_up)))*1e12,meanRT_down,'b');
title('Red=up, blue=down');xlabel('Force (pN)');ylabel('Mean Waiting Time (sec)');
elseif stiffnesses == 1
figure;
plot((stiffness_c-k_rand(1:length(spikefreq_up)))*1e6,spikefreq_up,'r');hold on;plot((stiffness_c-k_rand(1:length(spikefreq_up)))*1e6,spikefreq_down,'b');plot((stiffness_c-k_rand(1:length(spikefreq_up)))*1e6,sqrt((stiffness_c-k_rand(1:length(spikefreq_up)))*1e6),'b');
title('Red=up, blue=down');xlabel('Force (pN)');ylabel('Spike Frequency (spikes/sec)');
plot((stiffness_c-k_rand(1:length(spikefreq_up)))*1e6,meanRT_up,'r');hold on;plot((stiffness_c-k_rand(1:length(spikefreq_up)))*1e6,meanRT_down,'b');
title('Red=up, blue=down');xlabel('Force (pN)');ylabel('Mean Waiting Time (sec)');
end

% Plot statistics
if forces == 1
figure;
subplot(2,2,1);
plot(F_rand,RTupexp_stat5,'r');hold on;plot(F_rand,RTdownexp_stat5,'b'); ylabel('Lillie(exp) Stat');title('Red=up, blue=down');
subplot(2,2,2);
plot(F_rand,RTupgauss_stat5(1,:),'r');hold on;plot(F_rand,RTdowngauss_stat5(1,:),'b'); ylabel('Lillie(gauss) Stat');title('Red=up, blue=down');
subplot(2,2,3);
plot(F_rand,RTupgauss_stat5(2,:),'r');hold on;plot(F_rand,RTdowngauss_stat5(2,:),'b'); ylabel('JB Stat');title('Red=up, blue=down');
subplot(2,2,4);
plot(F_rand,RTupgauss_stat5(3,:),'r');hold on;plot(F_rand,RTdowngauss_stat5(3,:),'b'); ylabel('KS Stat');title('Red=up, blue=down');

elseif stiffnesses == 1
figure;
subplot(2,2,1);
plot(k_rand,RTupexp_stat5,'r');hold on;plot(k_rand,RTdownexp_stat5,'b'); ylabel('Lillie(exp) Stat');title('Red=up, blue=down');
subplot(2,2,2);
plot(k_rand,RTupgauss_stat5(1,:),'r');hold on;plot(k_rand,RTdowngauss_stat5(1,:),'b'); ylabel('Lillie(gauss) Stat');title('Red=up, blue=down');
subplot(2,2,3);
plot(k_rand,RTupgauss_stat5(2,:),'r');hold on;plot(k_rand,RTdowngauss_stat5(2,:),'b'); ylabel('JB Stat');title('Red=up, blue=down');
subplot(2,2,4);
plot(k_rand,RTupgauss_stat5(3,:),'r');hold on;plot(k_rand,RTdowngauss_stat5(3,:),'b'); ylabel('KS Stat');title('Red=up, blue=down');

end

% Loop through windows for Poisson distribution
figure;
warning off
for i = 1:length(poisswinmov_down5)
    for j = 1:5
        % Freedman-Diaconis Rule
        poissbinup(i) = 2*iqr(poisswinmov_up5{i}(:,j))*length(poisswinmov_up5{i}(:,j))^(-1/3); 
        poissbindown(i) = 2*iqr(poisswinmov_down5{i}(:,j))*length(poisswinmov_down5{i}(:,j))^(-1/3);
        nbinspoissup(i) = round((max(poisswinmov_up5{i}(:,j)) - min(poisswinmov_up5{i}(:,j)))/poissbinup(i));
        nbinspoissdown(i) = round((max(poisswinmov_down5{i}(:,j)) - min(poisswinmov_down5{i}(:,j)))/poissbindown(i));
        if isnan(nbinspoissup(i)) == 1 || isinf(nbinspoissup(i)) == 1
             nbinspoissup(i) = 5;
        end
        if isnan(nbinspoissdown(i)) == 1 || isinf(nbinspoissdown(i)) == 1
             nbinspoissdown(i) = 5;
        end
        % Create
        [poisshistupa{i,j},poisshistupb{i,j}]=hist(poisswinmov_up5{i}(:,j),nbinspoissup(i));
        [poisshistdowna{i,j},poisshistdownb{i,j}]=hist(poisswinmov_down5{i}(:,j),nbinspoissdown(i));
        % Fit to ALL : poisson, exponential, normal
        [poissfitupD{i,j},poissfitupPD{i,j}] = allfitdist(poisswinmov_up5{i}(:,j));
        [poissfitdownD{i,j},poissfitdownPD{i,j}] = allfitdist(poisswinmov_down5{i}(:,j));
        set(0,'DefaultAxesColorOrder',autumn(5));
        % Use histogram information in a chi-squared goodness-of-fit test
         n1 = sum(poisshistupa{i,j}); n2 = sum(poisshistdowna{i,j});
         pdup_poiss{i,j} = fitdist(poisshistupb{i,j}','Poisson','Frequency',poisshistupa{i,j}');
         pddown_poiss{i,j} = fitdist(poisshistdownb{i,j}','Poisson','Frequency',poisshistdowna{i,j}');
         expcountsup_poiss{i,j} = n1*pdf(pdup_poiss{i,j},poisshistupb{i,j}); expcountsdown_poiss{i,j} = n1*pdf(pddown_poiss{i,j},poisshistdownb{i,j});
         % Test for normal and exponential distributions using the chi-squared
         % test
         [~,Poissupchi2gof_p(i,j),Poissupchi2gof_stat{i,j}] = chi2gof(poisshistupb{i,j},'Ctrs',poisshistupb{i,j}, ...
            'Frequency', poisshistupa{i,j}, ...
            'Expected', expcountsup_poiss{i,j},...
            'Nparams',1);
         [~,Poissdownchi2gof_p(i,j),Poissdownchi2gof_stat{i,j}] = chi2gof(poisshistdownb{i,j},'Ctrs',poisshistdownb{i,j}, ...
            'Frequency', poisshistdowna{i,j}, ...
            'Expected', expcountsdown_poiss{i,j},...
            'Nparams',1);
        
         % Plot histograms
        set(gca,'LooseInset',get(gca,'TightInset'));
        subplot(2,length(poisswinmov_down5),i);set(gca,'LooseInset',get(gca,'TightInset'));hold all;plot(poisshistupb{i,j},poisshistupa{i,j});title('Spike Distribution');ylabel('# Events');xlabel('Spikes/Window');hold all;
        subplot(2,length(poisswinmov_down5),i+length(poisswinmov_down5));set(gca,'LooseInset',get(gca,'TightInset'));hold all;plot(poisshistupb{i,j},poisshistupa{i,j});title('Spike Distribution');ylabel('# Events');xlabel('Spikes/Window');hold all;

    end
end


disp('Saving...');
save([prefix 'analyzeddata.mat'],'poisshistupa','poisshistupb','poisshistdowna','poisshistdownb','poissfitupD','poissfitdownD','RTupdistD','RTdowndistD','meanRT_up','meanRT_down','RTupexp_stat5','RTdownexp_stat5','RTupgauss_stat5','RTdowngauss_stat5','upt','downt','F_rand','k_rand','time','spikefreq_up','spikefreq_down','Poissupchi2gof_p','Poissupchi2gof_stat','Poissdownchi2gof_p','Poissupchi2gof_stat','RTupgausschi2gof_p','RTupgausschi2gof_stat','RTdowngausschi2gof_p','RTdowngausschi2gof_stat','RTupexpchi2gof_p','RTupexpchi2gof_stat','RTdownexpchi2gof_p','RTdownexpchi2gof_stat','upa','upb','downb','downa','poisshistupa','poisshistupb','poisshistdownb','poisshistdowna');

%% DISCONTINUED Method
warning('DISCONTINUED METHOD');
% If figures==1 above, then other plots will be generated. This is highly
% recommended to avoid a large number of plots from being generated.

figure(1);
subplot(3,3,1); xlabel('High-Pass Filter Frequency');
clear m m1 n n1;m=mean(mean(upevents,3),2);m1(:,:)=m;n=mean(mean(downevents,3),2);n1(:,:)=n;plot(freq1range,m1,'r');hold on;plot(freq1range,n1,'b');ylabel('Mean Events (up=r;do=b)');xlabel('High-Pass Filter Frequency');
[fitupF, rupF]=fit(freq1range',m1,'exp1');[fitdownF, rdownF]=fit(freq1range',m1,'exp1'); 
subplot(3,3,4);plot(freq1range(2:end),diff(m1),'r');hold on;plot(freq1range(2:end),diff(n1),'b');ylabel('Diff(events)');xlabel('High-Pass Filter Frequency');
subplot(3,3,7);plot(freq1range(3:end),diff(diff(m1)),'r');hold on;plot(freq1range(3:end),diff(diff(n1)),'b');ylabel('Diff(Diff(events))');xlabel('High-Pass Filter Frequency');
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

subplot(3,3,2);xlabel('Detrending Window Size');
clear m m1 n n1;m=mean(mean(upevents,3),1);m1(:,:)=m;n=mean(mean(downevents,3),1);n1(:,:)=n;plot(winlooprange,m1,'r');hold on;plot(winlooprange,n1,'b');ylabel('Mean Events (up=r;do=b)');xlabel('Detrending Window Size');
[fitupW, rupW]=fit(winlooprange',m1','exp1');[fitdownW, rdownW]=fit(winlooprange',m1','exp1');
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

subplot(3,3,3);xlabel('Minimum RT Time');
clear m m1 n n1;m=(mean(mean(upevents,2),1));m1(:,:)=m;n=(mean(mean(downevents,2),1));n1(:,:)=n;plot(minTrange,m1,'r');hold on;plot(minTrange,n1,'b');ylabel('Mean Events (up=r;do=b)');xlabel('Minimum Residence Time');
[fitupT, rupT]=fit(minTrange',m1,'exp1');[fitdownT, rdownT]=fit(minTrange',m1,'exp1');
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

% Run RTfilt with the chosen parameters...
outfile = sprintf('%s%s%s%s%s%s%s%s%s%s%s%s',dir_out,'-HPFilt',num2str(mean([upF downF])),'-DetrendWin',num2str(mean([upW downW])),'-MinimumResidenceTime',num2str(mean([upT downT])),'-Fs',num2str(Fs),'-','POSTPLOTS','.mat');
RTfilt('a.mat',Fp,kp,mean([upW downW]),mean([upF downF]),Fs,1,outfile,mean([upT downT]),1,1);

