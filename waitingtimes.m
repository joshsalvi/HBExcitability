% INPUT PARAMETERS
close all;
clear file_RT dfit dfit2 file trace RT_up RT_down
 m=10;n=2;
file=Xd_pulse_center_m(10000:29999,m,n);
Timems=time.data(length(file):length(file)*2,4,6);
Fs=Fs;
 
analyze = 1; % 1=up, 2=down
minjump=0;   % minimum time in ms for a jump across 0
 
 
%------------------
trace=smooth(file,150);
 file=diff(trace);
 
figure(1);
plot(file,'k'); hold on; ylabel('nm');xlabel('ms');
file(1)=-0.1;
if file(1)>0
    a=1;
elseif file(1)<0
    a=-1;
elseif file(1)==0
    error('Warning, signal starts with a 0')
end
    
    
for i = 1:length(file)
    if file(i) > 2*std(file)
        file_RT(i)=1;
    elseif file(i) < -2*std(file)
        file_RT(i)=-1;
    else
        file_RT(i)=0;
    end
end
 
up = find(file_RT==1);
down = find(file_RT==0);
ju=1;
jd=1;
i=1;
while i < length(file_RT)
if file_RT(i)==-1
    RT_down(jd)=1/Fs;
    while file_RT(i)<1
        RT_down(jd)=RT_down(jd)+1/Fs;
        i=i+1;
        if i==length(file_RT)
            break;
        end
    end
    jd=jd+1;
end
 
if file_RT(i)==1
    RT_up(ju)=1/Fs;
    while file_RT(i)>-1
        RT_up(ju)=RT_up(ju)+1/Fs;
        i=i+1;
        if i==length(file_RT)
            break;
        end
    end
    ju=ju+1;
end
end

g=fittype('a*exp(-b*x)');
[down_a,down_b]=hist(RT_down,50);
down_a=down_a/trapz(down_b,down_a);     % use trapezoidal integration to normalize histogram 
dfit=fit(down_b(1:20)',down_a(1:20)',g);

[up_a,up_b]=hist(RT_up,50);
up_a=up_a/trapz(up_b,up_a);     % use trapezoidal integration to normalize histogram 
dfit2=fit(up_b(1:20)',up_a(1:20)',g);

RT_both=horzcat(RT_up,RT_down);
[both_a,both_b]=hist(RT_both,50);
both_a=both_a/trapz(both_b,both_a);     % use trapezoidal integration to normalize histogram 
dfit3=fit(both_b(1:20)',both_a(1:20)',g);

figure();
scatter(down_b,down_a,'k+'); hold on; plot(dfit); xlabel('Waiting Time (sec), down'); ylabel('P(t)'); title(['Force = ',num2str(Ord_F(m),3),'Stiffness = ',num2str(Ord_k(n),3),'y = a*exp(-b*x)+c; tau=',num2str(1/dfit.b,3)]);


figure();
scatter(up_b,up_a,'k+'); hold on; plot(dfit2); xlabel('Waiting Time (sec), up'); ylabel('P(t)'); title(['Force = ',num2str(Ord_F(m),3),'Stiffness = ',num2str(Ord_k(n),3),'y = a*exp(-b*x)+c; tau=',num2str(1/dfit2.b,3)]);


figure();
scatter(both_b,both_a,'k+'); hold on; plot(dfit3); xlabel('Waiting Time (sec), both'); ylabel('P(t)'); title(['Force = ',num2str(Ord_F(m),3),'Stiffness = ',num2str(Ord_k(n),3),'y = a*exp(-b*x)+c; tau=',num2str(1/dfit3.b,3)]);
