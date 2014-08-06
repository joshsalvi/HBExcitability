clear all; close all;

display('Importing...');
%load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/noisysinewave11movingavg_concwins2.mat');
load('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp4-kp4-start1end4e5-xfish1.0Noise-movingavg_concwins2.mat');
warning off;

N=200;                % choose number of stdev
Nmin=0.5;               % choose minimum value (mean-Nmin)
% Remove outliers
RT_up(RT_up>(mean(RT_up)+N*std(RT_up)) | RT_up<(mean(RT_up)-N*std(RT_up)) | RT_up < mean(RT_up)-Nmin) = [];
RT_down(RT_down>(mean(RT_down)+N*std(RT_down)) | RT_down<(mean(RT_down)-N*std(RT_down)) | RT_down < mean(RT_down)-Nmin) = [];
RT_buffer(RT_buffer>(mean(RT_buffer)+N*std(RT_buffer)) | RT_buffer<(mean(RT_buffer)-N*std(RT_buffer)) | RT_buffer < mean(RT_buffer)-Nmin) = [];
IEI_up(IEI_up>(mean(IEI_up)+N*std(IEI_up)) | IEI_up<(mean(IEI_up)-N*std(IEI_up)) | IEI_up < mean(IEI_up)-Nmin) = [];
IEI_down(IEI_down>(mean(IEI_down)+N*std(IEI_down)) | IEI_down<(mean(IEI_down)-N*std(IEI_down)) | IEI_down < mean(IEI_down)-Nmin) = [];
IEI_buffer(IEI_buffer>(mean(IEI_buffer)+N*std(IEI_buffer)) | IEI_buffer<(mean(IEI_buffer)-N*std(IEI_buffer)) | IEI_buffer < mean(IEI_buffer)-Nmin) = [];


if length(RT_up) < 4                    % If there are fewer than 4 observations, lillietest() will not work. Zero pad
    RT_up(length(RT_up)+1:4) = 0;
end
if length(RT_down) < 4
    RT_down(length(RT_down)+1:4) = 0;
end
if length(RT_buffer) < 4
    RT_buffer(length(RT_buffer)+1:4) = 0;
end
if length(IEI_up) < 4
    IEI_up(length(IEI_up)+1:4) = 0;
end
if length(IEI_down) < 4
    IEI_down(length(IEI_down)+1:4) = 0;
end
if length(IEI_buffer) < 4
    IEI_buffer(length(IEI_buffer)+1:4) = 0;
end

% Lilliefors' composite goodness-of-fit test for normal and exponential
% distributions
[RTupexp_h, RTupexp_p] = lillietest(RT_up,'Distr','exp');
[~, RTupgauss_p] = lillietest(RT_up);
[RTdownexp_h, RTdownexp_p] = lillietest(RT_down,'Distr','exp');
[~, RTdowngauss_p] = lillietest(RT_down);
[RTbufferexp_h, RTbufferexp_p] = lillietest(RT_buffer,'Distr','exp');
[~, RTbuffergauss_p] = lillietest(RT_buffer);
[IEIexp_h, IEIupexp_p] = lillietest(IEI_up,'Distr','exp');
[~, IEIupgauss_p] = lillietest(IEI_up);
[IEIdownexp_h, IEIdownexp_p] = lillietest(IEI_down,'Distr','exp');
[~, IEIdowngauss_p] = lillietest(IEI_down);
[IEIbufferexp_h, IEIbufferexp_p] = lillietest(IEI_buffer,'Distr','exp');
[~, IEIbuffergauss_p] = lillietest(IEI_buffer);
warning on;
% Kolmogorov-Smirnov goodness-of-fit hypothesis test for normal
% distributions
[~, RTupgauss_p(2)] = kstest(RT_up);
[~, RTdowngauss_p(2)] = kstest(RT_down);
[~, RTbuffergauss_p(2)] = kstest(RT_buffer);
[~, IEIupgauss_p(2)] = kstest(IEI_up);
[~, IEIdowngauss_p(2)] = kstest(IEI_down);
[~, IEIbuffergauss_p(2)] = kstest(IEI_buffer);

% Jarque-Bera hypothesis test of composite normality
[RTupgauss_h, RTupgauss_p(3)] = jbtest(RT_up);
[RTdowngauss_h, RTdowngauss_p(3)] = jbtest(RT_down);
[RTbuffergauss_h, RTbuffergauss_p(3)] = jbtest(RT_buffer);
[IEIupgauss_h, IEIupgauss_p(3)] = jbtest(IEI_up);
[IEIdowngauss_h, IEIdowngauss_p(3)] = jbtest(IEI_down);
[IEIbuffergauss_h, IEIbuffergauss_p(3)] = jbtest(IEI_buffer);


%Freedman-diaconis rules
binsizeup = 2*iqr(RT_up)*length(RT_up)^(-1/3);      % freedman-diaconis rule
nbinsup = round((max(RT_up) - min(RT_up))/binsizeup);
if nbinsup > 1e6
    nbinsup = 1e6;
end
binsizedown = 2*iqr(RT_down)*length(RT_down)^(-1/3);      % freedman-diaconis rule
nbinsdown = round((max(RT_down) - min(RT_down))/binsizedown);
if nbinsdown > 1e6
    nbinsdown = 1e6;
end
binsizeupiei = 2*iqr(IEI_up)*length(IEI_up)^(-1/3);      % freedman-diaconis rule
nbinsupiei = round((max(IEI_up) - min(IEI_up))/binsizeupiei);
if nbinsupiei > 1e6
    nbinsupiei = 1e6;
end
binsizedowniei = 2*iqr(IEI_down)*length(IEI_down)^(-1/3);      % freedman-diaconis rule
nbinsdowniei = round((max(IEI_down) - min(IEI_down))/binsizedowniei);
if nbinsdowniei > 1e6
    nbinsdowniei = 1e6;
end
binsizebufferiei = 2*iqr(IEI_buffer)*length(IEI_buffer)^(-1/3);      % freedman-diaconis rule
nbinsbufferiei = round((max(IEI_buffer) - min(IEI_buffer))/binsizebufferiei);
if nbinsbufferiei > 1e6
    nbinsbufferiei = 1e6;
end
binsizebuffer = 2*iqr(RT_buffer)*length(RT_buffer)^(-1/3);      % freedman-diaconis rule
nbinsbuffer = round((max(RT_buffer) - min(RT_buffer))/binsizebuffer);
if nbinsbuffer > 1e6
    nbinsbuffer = 1e6;
end

figure;
subplot(3,4,1);histfit(RT_up,nbinsup); title('RT up');
subplot(3,4,2);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.1f\nStd. Dev. = %.1f\nMedian = %.1f\nMin = %.1f\nMax = %.1f\nLillie(exp) p = %.1e\nLillie(gauss) p = %.1e\nKS p = %.1e\nJB p = %.1e',mean(RT_up),std(RT_up),median(RT_up),min(RT_up),max(RT_up),RTupexp_p,RTupgauss_p(1),RTupgauss_p(2),RTupgauss_p(3)));
subplot(3,4,3);histfit(RT_down,nbinsdown); title('RT down');
subplot(3,4,4);set(gca, 'visible', 'off');text(0,0.6,sprintf('Mean = %.1f\nStd. Dev. = %.1f\nMedian = %.1f\nMin = %.1f\nMax = %.1f\nLillie(exp) p = %.1e\nLillie(gauss) p = %.1e\nKS p = %.1e\nJB p = %.1e',mean(RT_down),std(RT_down),median(RT_down),min(RT_down),max(RT_down),RTdownexp_p,RTdowngauss_p(1),RTdowngauss_p(2),RTdowngauss_p(3)));
subplot(3,4,5);histfit(IEI_up,nbinsupiei); title('IEI up');
subplot(3,4,6);set(gca, 'visible', 'off');text(0,.6,sprintf('Mean = %.1f\nStd. Dev. = %.1f\nMedian = %.1f\nMin = %.1f\nMax = %.1f\nLillie(exp) p = %.1e\nLillie(gauss) p = %.1e\nKS p = %.1e\nJB p = %.1e',mean(IEI_up),std(IEI_up),median(IEI_up),min(IEI_up),max(IEI_up),IEIupexp_p,IEIupgauss_p(1),IEIupgauss_p(2),IEIupgauss_p(3)));
subplot(3,4,7);histfit(IEI_down,nbinsdowniei); title('IEI down');
subplot(3,4,8);set(gca, 'visible', 'off');text(0,.6,sprintf('Mean = %.1f\nStd. Dev. = %.1f\nMedian = %.1f\nMin = %.1f\nMax = %.1f\nLillie(exp) p = %.1e\nLillie(gauss) p = %.1e\nKS p = %.1e\nJB p = %.1e',mean(IEI_down),std(IEI_down),median(IEI_down),min(IEI_down),max(IEI_down),IEIdownexp_p,IEIdowngauss_p(1),IEIdowngauss_p(2),IEIdowngauss_p(3)));
subplot(3,4,9);histfit(RT_buffer,nbinsbuffer); title('RT buffer');
subplot(3,4,10);set(gca, 'visible', 'off');text(0,.6,sprintf('Mean = %.1f\nStd. Dev. = %.1f\nMedian = %.1f\nMin = %.1f\nMax = %.1f\nLillie(exp) p = %.1e\nLillie(gauss) p = %.1e\nKS p = %.1e\nJB p = %.1e',mean(RT_buffer),std(RT_buffer),median(RT_buffer),min(RT_buffer),max(RT_buffer),RTbufferexp_p,RTbuffergauss_p(1),RTbuffergauss_p(2),RTbuffergauss_p(3)));
subplot(3,4,11);histfit(IEI_buffer,nbinsbufferiei); title('IEI buffer');
subplot(3,4,12);set(gca, 'visible', 'off');text(0,.6,sprintf('Mean = %.1f\nStd. Dev. = %.1f\nMedian = %.1f\nMin = %.1f\nMax = %.1f\nLillie(exp) p = %.1e\nLillie(gauss) p = %.1e\nKS p = %.1e\nJB p = %.1e',mean(IEI_buffer),std(IEI_buffer),median(IEI_buffer),min(IEI_buffer),max(IEI_buffer),IEIbufferexp_p,IEIbuffergauss_p(1),IEIbuffergauss_p(2),IEIbuffergauss_p(3)));

figure;
subplot(3,2,1);qqplot(RT_up);title('RT up');
subplot(3,2,2);qqplot(RT_down);title('RT down');
subplot(3,2,3);qqplot(IEI_up);title('IEI up');
subplot(3,2,4);qqplot(IEI_down);title('IEI down');
subplot(3,2,5);qqplot(RT_buffer);title('RT buffer');
subplot(3,2,6);qqplot(IEI_buffer);title('IEI buffer');


display('Saving...');
%save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/Sinusoids/noisysinewave11movingavg_concwins2_stats2_removedoutliersN3.mat','RTupexp_p','RTupgauss_p','RTdownexp_p','RTdowngauss_p','RTbufferexp_p','RTbuffergauss_p','IEIdownexp_p','IEIdowngauss_p','IEIupexp_p','IEIupgauss_p','IEIbufferexp_p','IEIbuffergauss_p');
%save('/Users/joshsalvi/Documents/Lab/Lab/Simulation Data/ONHFishJosh/Analysis/Fp1-kp2-start1end2e5-xfish1.0Noise-movingavg_concwins2_stats_removeoutliersN200Nmin01.mat','RTupexp_p','RTupgauss_p','RTdownexp_p','RTdowngauss_p','RTbufferexp_p','RTbuffergauss_p','IEIdownexp_p','IEIdowngauss_p','IEIupexp_p','IEIupgauss_p','IEIbufferexp_p','IEIbuffergauss_p');

display('Finished.');

