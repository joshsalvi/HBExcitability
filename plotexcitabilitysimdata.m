

if biftype == 1
    biftype2='supercritical Hopf';
elseif biftype == 2
    biftype2='SNIC';
elseif biftype == 3
    biftype2='subcritical Hopf';
end

if biftype == 1 || biftype == 3
    I = mu;
end


for i = 1:length(Xdet)
    figure(1);subplot(ceil(length(Xdet)/2),2,i);plot(I,pkspikeratedet(i,:),'k--');hold on;plot(I,pkspikeratesto(i,:),'r');title(sprintf('%s%s %s%s','Peaks, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Spike Rate');
    figure(2);subplot(ceil(length(Xdet)/2),2,i);plot(I,trspikeratedet(i,:),'k--');hold on;plot(I,trspikeratesto(i,:),'r');title(sprintf('%s%s %s%s','Troughs, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Spike Rate');
end;clear i

for i = 1:length(Xdet)
    figure(3);subplot(ceil(length(Xdet)/2),2,i);plot(I,CDdetpk(i,:),'k--');hold on;plot(I,CDstopk(i,:),'r');plot(I,ones(1,length(I)),'g--');title(sprintf('%s%s %s%s','Peaks, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Coeff. Dispersion');
    figure(4);subplot(ceil(length(Xdet)/2),2,i);plot(I,CDdettr(i,:),'k--');hold on;plot(I,CDstotr(i,:),'r');plot(I,ones(1,length(I)),'g--');title(sprintf('%s%s %s%s','Troughs, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Coeff. Dispersion');
end;clear i

for i = 1:length(Xdet)
    figure(5);subplot(ceil(length(Xdet)/2),2,i);plot(I,pkdiffusiondet(i,:),'k--');hold on;plot(I,pkdiffusionsto(i,:),'r');title(sprintf('%s%s %s%s','Peaks, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Diffusion Coeff.');
    figure(6);subplot(ceil(length(Xdet)/2),2,i);plot(I,trdiffusiondet(i,:),'k--');hold on;plot(I,trdiffusionsto(i,:),'r');title(sprintf('%s%s %s%s','Troughs, ',biftype2,' Noise level = ',num2str(noiselevel(i))));xlabel('I');ylabel('Diffusion Coeff.');
end;clear i


if biftype == 1 || biftype == 3
    clear I;
end

