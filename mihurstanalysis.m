clear all; close all; clc

N = 5e3;               % maximum signal length
N2 = 500;              % number of iterations (divisions)
N3 = 20;                % number of averages
gp = 2^12;              % number of grid points
D=0.0001;
dimensions = 1;         % two dimensional simulation
tau = .1;               % time interval in seconds
time = tau * 1:N;       % create a time vector for plotting
Hin = [0.25 0.5 0.75];  % Hurst exponents (input)

for l = 1:N3

    
k = sqrt(D * dimensions * tau);
dy = k * randn(N,1);
%}

%{
for j = 1:length(Hin)
    y{j,l} = generateBM(gp,Hin(j));
    yi{j,l} = imag(hilbert(y{j,l}));
end
N = length(y{1,1});
time = tau * 1:N; 
%}

% Normal diffusion (H ~ 0.5)
y{1,l} = cumsum(dy);yi{1,l}=imag(hilbert(y{1,l}));
disp(['Normal Diffusion; H = ' num2str(genhurst(y{1,l}))]);

% Super diffusion (H > 0.5)
y{2,l} = cumsum(randn(1,N).^3);yi{2,l}=imag(hilbert(y{2,l}));
disp(['Super Diffusion; H = ' num2str(genhurst(y{2,l}))]);

% Sub diffusion (H > 0.5)
y{3,l} = cumsum(randn(1,N))+randn(1,N);yi{3,l}=imag(hilbert(y{3,l}));
disp(['Sub Diffusion; H = ' num2str(genhurst(y{3,l}))]);


%}

end

for l = 1:N3
    % Take middle 60% of signal
for j = 1:3
    y{j,l} = y{j,l}(N/5:end-N/5);
    yi{j,l} = yi{j,l}(N/5:end-N/5);
end

end
for l = 1:N3
setfiguredefaults(4);
if l == 1
time = time(N/5:end-N/5);
N = length(y{j,l});
end
for j = 1:3
    %figure(j);subplot(2,ceil(N3/2),l);plot(time.*1e-2,y{j,l},'k');hold on;plot(time.*1e-2,yi{j,l},'r');
    %xlabel('Time (sec * 10^2)');ylabel('Position');legend([{'X'};{'X_H'}]);
    %xh = num2str(genhurst(y{j,l}));xh=xh(1:4);
    %title(['H = ' xh]);
end
end
%%    
for l = 1:N3
for j = 1:N2
    y1r = y{1,l}(1:j*floor(N/N2));y1i=yi{1,l}(1:j*floor(N/N2));Nb1r(j,l)=freedmandiaconis(y1r);Nb1i(j,l)=Nb1r(j,l);
    y2r = y{2,l}(1:j*floor(N/N2));y2i=yi{2,l}(1:j*floor(N/N2));Nb2r(j,l)=freedmandiaconis(y2r);Nb2i(j,l)=Nb2r(j,l);
    y3r = y{3,l}(1:j*floor(N/N2));y3i=yi{3,l}(1:j*floor(N/N2));Nb3r(j,l)=freedmandiaconis(y3r);Nb3i(j,l)=Nb3r(j,l);
    nn1=j*floor(N/N2);
    [a1r b1r]=hist(y1r,Nb1r(j,l));[a1i b1i]=hist(y1i,Nb1i(j,l));
    [a2r b2r]=hist(y2r,Nb2r(j,l));[a2i b2i]=hist(y2i,Nb2i(j,l));
    [a3r b3r]=hist(y3r,Nb3r(j,l));[a3i b3i]=hist(y3i,Nb3i(j,l));
    E1r(j,l) = h(a1r');E1i(j,l) = h(a1i');
    E2r(j,l) = h(a2r');E2i(j,l) = h(a2i');
    E3r(j,l) = h(a3r');E3i(j,l) = h(a3i');
    MI1(j,l) = mi(a1r',a1i');
    MI2(j,l) = mi(a2r',a2i');
    MI3(j,l) = mi(a3r',a3i');
    H1r(j,l) = genhurst(y1r);H1i(j,l) = genhurst(y1i);
    H2r(j,l) = genhurst(y2r);H2i(j,l) = genhurst(y2i);
    H3r(j,l) = genhurst(y3r);H3i(j,l) = genhurst(y3i);
    E1r2(j,l) = log2( sqrt(2*pi*exp(1)) ) + H1r(j,l)*log2(j*floor(N/N2));
    E1i2(j,l) = log2( sqrt(2*pi*exp(1)) ) + H1i(j,l)*log2(j*floor(N/N2));
    E2r2(j,l) = log2( sqrt(2*pi*exp(1)) ) + H2r(j,l)*log2(j*floor(N/N2));
    E2i2(j,l) = log2( sqrt(2*pi*exp(1)) ) + H2i(j,l)*log2(j*floor(N/N2));
    E3r2(j,l) = log2( sqrt(2*pi*exp(1)) ) + H3r(j,l)*log2(j*floor(N/N2));
    E3i2(j,l) = log2( sqrt(2*pi*exp(1)) ) + H3i(j,l)*log2(j*floor(N/N2));
    
    if mod(j,50)==0
        disp([num2str(j) '/' num2str(N2)]);
    end
    C=cov(y1r,y1i);cov1(j,l)=C(1,2);
    C=cov(y2r,y2i);cov2(j,l)=C(1,2);
    C=cov(y3r,y3i);cov3(j,l)=C(1,2);
    
    MI12(j,l) = 0.5*(H1r(j,l) + H1i(j,l))*log2(nn1) - ((nn1^(2*H1r(j,l))+nn1^(2*H1i(j,l))))/(2*nn1^(H1r(j,l)+H1i(j,l)))+log2(2*pi*exp(1))-log2(sqrt(2*pi));
    MI22(j,l) = 0.5*(H2r(j,l) + H2i(j,l))*log2(nn1) - ((nn1^(2*H2r(j,l))+nn1^(2*H2i(j,l))))/(2*nn1^(H2r(j,l)+H2i(j,l)))+log2(2*pi*exp(1))-log2(sqrt(2*pi));
    MI32(j,l) = 0.5*(H3r(j,l) + H3i(j,l))*log2(nn1) - ((nn1^(2*H3r(j,l))+nn1^(2*H3i(j,l))))/(2*nn1^(H3r(j,l)+H3i(j,l)))+log2(2*pi*exp(1))-log2(sqrt(2*pi));


end
disp([num2str(l) '/' num2str(N3)]);
end
    

E1rmean=mean(E1r,2);E1imean=mean(E1i,2);
E2rmean=mean(E2r,2);E2imean=mean(E2i,2);
E3rmean=mean(E3r,2);E3imean=mean(E3i,2);
H1rmean=mean(H1r,2);H1imean=mean(H1i,2);
H2rmean=mean(H2r,2);H2imean=mean(H2i,2);
H3rmean=mean(H3r,2);H3imean=mean(H3i,2);
MI1mean=mean(MI1,2);MI12mean=mean(MI12,2);
MI2mean=mean(MI2,2);MI22mean=mean(MI22,2);
MI3mean=mean(MI3,2);MI32mean=mean(MI32,2);
E1r2mean=mean(E1r2,2);E1i2mean=mean(E1i2,2);
E2r2mean=mean(E2r2,2);E2i2mean=mean(E2i2,2);
E3r2mean=mean(E3r2,2);E3i2mean=mean(E3i2,2);


    disp('Finished.');
    
    
    
