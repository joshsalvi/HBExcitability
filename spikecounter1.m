Nw = 300;


for j = 1:a
    Xd_pulse_center(:,j) = Xd_pulse_center(:,j)-mean(Xd_pulse_center(:,j));
clear q q1 q2 q3 w1
w1 = smooth(Xd_pulse_center(:,1,j),Nw);
q=find(w1<(mean(w1)-std(w1)));
q1=zeros(1,length(Xd_pulse_center));
q1(q)=1;
q2=diff(q1);
q3 = find(diff(q1)>0.5);

spike_count(j) = length(q3);
spike_freq(j) = spike_count(j)/(length(Xd_pulse_center)/Fs);


end


for j = 1:a
    clear w1
    w1 = smooth(Xd_pulse_center(:,1,j),Nw);
    bw=2*iqr(Xd_pulse_center(:,j))/length(Xd_pulse_center(:,j))^(1/3);
    Nbins=round((max(Xd_pulse_center(:,j))-min(Xd_pulse_center(:,j)))/bw);
    
    subplot(ceil(sqrt(a)),ceil(sqrt(a)),j); hist(w1,Nbins);
end
