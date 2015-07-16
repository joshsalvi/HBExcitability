mu = linspace(-0.5,0,500);
t=linspace(0,1e4,5e5);

%% SUPERCRITICAL HOPF
%{
mu = linspace(-1,1,500);
t=linspace(0,1e4,5e5);
rt=mu;
for i = 1:200
    [XHopfdet005{i} XHopfsto005{i}] = hopfstoch(rt(i),1,0.05,0.05,t);
    [XHopfdet01{i} XHopfsto01{i}] = hopfstoch(rt(i),1,0.1,0.1,t);
    [XHopfdet02{i} XHopfsto02{i}] = hopfstoch(rt(i),1,0.2,0.2,t);
    [XHopfdet04{i} XHopfsto04{i}] = hopfstoch(rt(i),1,0.4,0.4,t);
end

Hopfdet{:,1}=XHopfdet005;
Hopfdet{:,2}=XHopfdet01;
Hopfdet{:,3}=XHopfdet02;
Hopfdet{:,4}=XHopfdet04;
Hopfsto{:,1}=XHopfsto005;
Hopfsto{:,2}=XHopfsto01;
Hopfsto{:,3}=XHopfsto02;
Hopfsto{:,4}=XHopfsto04;
Xdet = Hopfdet; Xsto = Hopfsto;
clear XHopfdet005 XHopfdet01 XHopfdet02 XHopfdet04 XHopfsto005 XHopfsto01 XHopfsto02 XHopfsto04
clear Hopfdet Hopfsto
save('/Users/joshsalvi/Desktop/Hopf/Hopfstochoutput-all-1.mat','-v7.3')
%}
disp('Hopf 1');
clear

mu = linspace(-1,1,500);
t=linspace(0,1e4,5e5);
rt=mu;
for i = 201:350
    [XHopfdet005{i} XHopfsto005{i}] = hopfstoch(rt(i),1,0.05,0.05,t);
    [XHopfdet01{i} XHopfsto01{i}] = hopfstoch(rt(i),1,0.1,0.1,t);
    [XHopfdet02{i} XHopfsto02{i}] = hopfstoch(rt(i),1,0.2,0.2,t);
    [XHopfdet04{i} XHopfsto04{i}] = hopfstoch(rt(i),1,0.4,0.4,t);
end

Hopfdet{:,1}=XHopfdet005;
Hopfdet{:,2}=XHopfdet01;
Hopfdet{:,3}=XHopfdet02;
Hopfdet{:,4}=XHopfdet04;
Hopfsto{:,1}=XHopfsto005;
Hopfsto{:,2}=XHopfsto01;
Hopfsto{:,3}=XHopfsto02;
Hopfsto{:,4}=XHopfsto04;
Xdet = Hopfdet; Xsto = Hopfsto;
clear XHopfdet005 XHopfdet01 XHopfdet02 XHopfdet04 XHopfsto005 XHopfsto01 XHopfsto02 XHopfsto04
clear Hopfdet Hopfsto
save('/Users/joshsalvi/Desktop/Hopf/Hopfstochoutput-all-2.mat','-v7.3')
disp('Hopf 2');
clear


mu = linspace(-1,1,500);
t=linspace(0,1e4,5e5);
rt=mu;
for i = 351:500
    [XHopfdet005{i} XHopfsto005{i}] = hopfstoch(rt(i),1,0.05,0.05,t);
    [XHopfdet01{i} XHopfsto01{i}] = hopfstoch(rt(i),1,0.1,0.1,t);
    [XHopfdet02{i} XHopfsto02{i}] = hopfstoch(rt(i),1,0.2,0.2,t);
    [XHopfdet04{i} XHopfsto04{i}] = hopfstoch(rt(i),1,0.4,0.4,t);
end

Hopfdet{:,1}=XHopfdet005;
Hopfdet{:,2}=XHopfdet01;
Hopfdet{:,3}=XHopfdet02;
Hopfdet{:,4}=XHopfdet04;
Hopfsto{:,1}=XHopfsto005;
Hopfsto{:,2}=XHopfsto01;
Hopfsto{:,3}=XHopfsto02;
Hopfsto{:,4}=XHopfsto04;
Xdet = Hopfdet; Xsto = Hopfsto;
clear XHopfdet005 XHopfdet01 XHopfdet02 XHopfdet04 XHopfsto005 XHopfsto01 XHopfsto02 XHopfsto04
clear Hopfdet Hopfsto
save('/Users/joshsalvi/Desktop/Hopf/Hopfstochoutput-all-3.mat','-v7.3')
disp('Hopf 3');
clear

%%
% SUBCRITICAL HOPF

mu = linspace(-0.5,0,500);
t=linspace(0,1e4,5e5);
rt=mu;
for i=137:150
   rt=mu;
[XsubHopfdet01{i}, XsubHopfsto01{i}] = subhopfstoch(rt(i),1,0.1,0.1,t);
[XsubHopfdet02{i}, XsubHopfsto02{i}] = subhopfstoch(rt(i),1,0.2,0.2,t);
[XsubHopfdet04{i}, XsubHopfsto04{i}] = subhopfstoch(rt(i),1,0.4,0.4,t);
[XsubHopfdet005{i}, XsubHopfsto005{i}] = subhopfstoch(rt(i),1,0.05,0.05,t);
end
noiselevel=[0.05 0.1 0.2 0.4 0.6 0.8 1];


subHopfdet{:,1}=XsubHopfdet005;
subHopfdet{:,2}=XsubHopfdet01;
subHopfdet{:,3}=XsubHopfdet02;
subHopfdet{:,4}=XsubHopfdet04;
subHopfsto{:,1}=XsubHopfsto005;
subHopfsto{:,2}=XsubHopfsto01;
subHopfsto{:,3}=XsubHopfsto02;
subHopfsto{:,4}=XsubHopfsto04;
Xdet = subHopfdet; Xsto = subHopfsto;
clear subHopfdet subHopfsto
save('/Users/joshsalvi/Desktop/subHopf/subHopfstochoutput-all-1.mat','-v7.3')
disp('subHopf 1');
clear
%}

mu = linspace(-0.5,0,500);
t=linspace(0,1e4,5e5);
rt=mu;
for i=151:300
   rt=mu;
[XsubHopfdet01{i}, XsubHopfsto01{i}] = subhopfstoch(rt(i),1,0.1,0.1,t);
[XsubHopfdet02{i}, XsubHopfsto02{i}] = subhopfstoch(rt(i),1,0.2,0.2,t);
[XsubHopfdet04{i}, XsubHopfsto04{i}] = subhopfstoch(rt(i),1,0.4,0.4,t);
[XsubHopfdet005{i}, XsubHopfsto005{i}] = subhopfstoch(rt(i),1,0.05,0.05,t);
end
noiselevel=[0.05 0.1 0.2 0.4 0.6 0.8 1];


subHopfdet{:,1}=XsubHopfdet005;
subHopfdet{:,2}=XsubHopfdet01;
subHopfdet{:,3}=XsubHopfdet02;
subHopfdet{:,4}=XsubHopfdet04;
subHopfsto{:,1}=XsubHopfsto005;
subHopfsto{:,2}=XsubHopfsto01;
subHopfsto{:,3}=XsubHopfsto02;
subHopfsto{:,4}=XsubHopfsto04;
Xdet = subHopfdet; Xsto = subHopfsto;
clear subHopfdet subHopfsto
save('/Users/joshsalvi/Desktop/subHopf/subHopfstochoutput-all-2.mat','-v7.3')
disp('subHopf 2');
clear
%}

mu = linspace(-0.5,0,500);
t=linspace(0,1e4,5e5);
rt=mu;
for i=301:400
   rt=mu;
[XsubHopfdet01{i}, XsubHopfsto01{i}] = subhopfstoch(rt(i),1,0.1,0.1,t);
[XsubHopfdet02{i}, XsubHopfsto02{i}] = subhopfstoch(rt(i),1,0.2,0.2,t);
[XsubHopfdet04{i}, XsubHopfsto04{i}] = subhopfstoch(rt(i),1,0.4,0.4,t);
[XsubHopfdet005{i}, XsubHopfsto005{i}] = subhopfstoch(rt(i),1,0.05,0.05,t);
end
noiselevel=[0.05 0.1 0.2 0.4 0.6 0.8 1];


subHopfdet{:,1}=XsubHopfdet005;
subHopfdet{:,2}=XsubHopfdet01;
subHopfdet{:,3}=XsubHopfdet02;
subHopfdet{:,4}=XsubHopfdet04;
subHopfsto{:,1}=XsubHopfsto005;
subHopfsto{:,2}=XsubHopfsto01;
subHopfsto{:,3}=XsubHopfsto02;
subHopfsto{:,4}=XsubHopfsto04;

Xdet = subHopfdet; Xsto = subHopfsto;
clear subHopfdet subHopfsto
save('/Users/joshsalvi/Desktop/subHopf/subHopfstochoutput-all-3.mat','-v7.3')
disp('SubHopf 3');
clear
%}

mu = linspace(-0.5,0,500);
t=linspace(0,1e4,5e5);
rt=mu;
for i=401:500
   rt=mu;
[XsubHopfdet01{i}, XsubHopfsto01{i}] = subhopfstoch(rt(i),1,0.1,0.1,t);
[XsubHopfdet02{i}, XsubHopfsto02{i}] = subhopfstoch(rt(i),1,0.2,0.2,t);
[XsubHopfdet04{i}, XsubHopfsto04{i}] = subhopfstoch(rt(i),1,0.4,0.4,t);
[XsubHopfdet005{i}, XsubHopfsto005{i}] = subhopfstoch(rt(i),1,0.05,0.05,t);
end
noiselevel=[0.05 0.1 0.2 0.4 0.6 0.8 1];


subHopfdet{:,1}=XsubHopfdet005;
subHopfdet{:,2}=XsubHopfdet01;
subHopfdet{:,3}=XsubHopfdet02;
subHopfdet{:,4}=XsubHopfdet04;
subHopfsto{:,1}=XsubHopfsto005;
subHopfsto{:,2}=XsubHopfsto01;
subHopfsto{:,3}=XsubHopfsto02;
subHopfsto{:,4}=XsubHopfsto04;
Xdet = subHopfdet; Xsto = subHopfsto;
clear subHopfdet subHopfsto
save('/Users/joshsalvi/Desktop/subHopf/subHopfstochoutput-all-4.mat','-v7.3')
disp('SubHopf 4');
clear
%% SNIC
%{
mu = linspace(-1,1,500);
t=linspace(0,1e4,5e5);

for i=1:100
    I=mu;rt=I;
    [XSNICdet005{i} XSNICsto005{i}] = thetamodel(rt(i),0.05,t);
[XSNICdet01{i} XSNICsto01{i}] = thetamodel(rt(i),0.1,t);
[XSNICdet02{i} XSNICsto02{i}] = thetamodel(rt(i),0.2,t);
[XSNICdet04{i} XSNICsto04{i}] = thetamodel(rt(i),0.4,t);
end

SNICdet{:,1}=XSNICdet005;
SNICdet{:,2}=XSNICdet01;
SNICdet{:,3}=XSNICdet02;
SNICdet{:,4}=XSNICdet04;
SNICsto{:,1}=XSNICsto005;
SNICsto{:,2}=XSNICsto01;
SNICsto{:,3}=XSNICsto02;
SNICsto{:,4}=XSNICsto04;
Xdet = SNICdet; Xsto = SNICsto;
clear SNICdet SNICsto XSNICdet005 XSNICsto005 XSNICdet01 XSNICdet01 XSNICsto01 XSNICdet02 XSNICsto02 XSNICdet04 XSNICsto04
save('/Volumes/LaCie JDS/Simulation Data/SNIC/SNICcos-finerI-1.mat','-v7.3');
disp('SNIC1');
%}
%{
mu = linspace(-1,1,500);
t=linspace(0,1e4,5e5);
for i=101:200
    I=mu;rt=I;
    [XSNICdet005{i} XSNICsto005{i}] = thetamodel(rt(i),0.05,t);
[XSNICdet01{i} XSNICsto01{i}] = thetamodel(rt(i),0.1,t);
[XSNICdet02{i} XSNICsto02{i}] = thetamodel(rt(i),0.2,t);
[XSNICdet04{i} XSNICsto04{i}] = thetamodel(rt(i),0.4,t);
end

SNICdet{:,1}=XSNICdet005;
SNICdet{:,2}=XSNICdet01;
SNICdet{:,3}=XSNICdet02;
SNICdet{:,4}=XSNICdet04;
SNICsto{:,1}=XSNICsto005;
SNICsto{:,2}=XSNICsto01;
SNICsto{:,3}=XSNICsto02;
SNICsto{:,4}=XSNICsto04;
Xdet = SNICdet; Xsto = SNICsto;
clear SNICdet SNICsto XSNICdet005 XSNICsto005 XSNICdet01 XSNICdet01 XSNICsto01 XSNICdet02 XSNICsto02 XSNICdet04 XSNICsto04
save('/Volumes/LaCie JDS/Simulation Data/SNIC/SNICcos-finerI-2.mat','-v7.3');
disp('SNIC2');
%}
%{
mu = linspace(-1,1,500);
t=linspace(0,1e4,5e5);
for i=201:300
    I=mu;rt=I;
    [XSNICdet005{i} XSNICsto005{i}] = thetamodel(rt(i),0.05,t);
[XSNICdet01{i} XSNICsto01{i}] = thetamodel(rt(i),0.1,t);
[XSNICdet02{i} XSNICsto02{i}] = thetamodel(rt(i),0.2,t);
[XSNICdet04{i} XSNICsto04{i}] = thetamodel(rt(i),0.4,t);
end

SNICdet{:,1}=XSNICdet005;
SNICdet{:,2}=XSNICdet01;
SNICdet{:,3}=XSNICdet02;
SNICdet{:,4}=XSNICdet04;
SNICsto{:,1}=XSNICsto005;
SNICsto{:,2}=XSNICsto01;
SNICsto{:,3}=XSNICsto02;
SNICsto{:,4}=XSNICsto04;
Xdet = SNICdet; Xsto = SNICsto;
clear SNICdet SNICsto XSNICdet005 XSNICsto005 XSNICdet01 XSNICdet01 XSNICsto01 XSNICdet02 XSNICsto02 XSNICdet04 XSNICsto04
save('/Volumes/LaCie JDS/Simulation Data/SNIC/SNICcos-finerI-3.mat','-v7.3');
disp('SNIC3');
%}


mu = linspace(-1,1,500);
t=linspace(0,1e4,5e5);
for i=336:400
    I=mu;rt=I;
    [XSNICdet005{i} XSNICsto005{i}] = thetamodel(rt(i),0.05,t);
[XSNICdet01{i} XSNICsto01{i}] = thetamodel(rt(i),0.1,t);
[XSNICdet02{i} XSNICsto02{i}] = thetamodel(rt(i),0.2,t);
[XSNICdet04{i} XSNICsto04{i}] = thetamodel(rt(i),0.4,t);
end

SNICdet{:,1}=XSNICdet005;
SNICdet{:,2}=XSNICdet01;
SNICdet{:,3}=XSNICdet02;
SNICdet{:,4}=XSNICdet04;
SNICsto{:,1}=XSNICsto005;
SNICsto{:,2}=XSNICsto01;
SNICsto{:,3}=XSNICsto02;
SNICsto{:,4}=XSNICsto04;
Xdet = SNICdet; Xsto = SNICsto;
clear SNICdet SNICsto XSNICdet005 XSNICsto005 XSNICdet01 XSNICdet01 XSNICsto01 XSNICdet02 XSNICsto02 XSNICdet04 XSNICsto04
save('/Volumes/LaCie JDS/Simulation Data/SNIC/SNICcos-finerI-4.mat','-v7.3');
disp('SNIC4');
%}

mu = linspace(-1,1,500);
t=linspace(0,1e4,5e5);
for i=401:500
    I=mu;rt=I;
    [XSNICdet005{i} XSNICsto005{i}] = thetamodel(rt(i),0.05,t);
[XSNICdet01{i} XSNICsto01{i}] = thetamodel(rt(i),0.1,t);
[XSNICdet02{i} XSNICsto02{i}] = thetamodel(rt(i),0.2,t);
[XSNICdet04{i} XSNICsto04{i}] = thetamodel(rt(i),0.4,t);
end

SNICdet{:,1}=XSNICdet005;
SNICdet{:,2}=XSNICdet01;
SNICdet{:,3}=XSNICdet02;
SNICdet{:,4}=XSNICdet04;
SNICsto{:,1}=XSNICsto005;
SNICsto{:,2}=XSNICsto01;
SNICsto{:,3}=XSNICsto02;
SNICsto{:,4}=XSNICsto04;
Xdet = SNICdet; Xsto = SNICsto;
clear SNICdet SNICsto XSNICdet005 XSNICsto005 XSNICdet01 XSNICdet01 XSNICsto01 XSNICdet02 XSNICsto02 XSNICdet04 XSNICsto04
save('/Volumes/LaCie JDS/Simulation Data/SNIC/SNICcos-finerI-5.mat','-v7.3');
disp('SNIC5');
%}
%% HAIR BUNDLE MODEL


k1 = [1.5 1.75 2 2.5];
F1 = 0;
ke = linspace(1.5,4,500);
F = linspace(0,2,500);
noiselevel = [0.05 0.1 0.2];
t=linspace(0,1e6,5e6);

for i = 1:100
    for j = 1:length(k1)
        for k = 1:length(noiselevel)
            [XdetF{i,j,k}, XstoF{i,j,k}] = hbtoymodel(F(i),k1(j),noiselevel(k),0,2,t);    % force scans
        end
    end
end
save('/Users/joshsalvi/Desktop/toymodel/HBmodel-force-5e6-1.mat','-v7.3');
clear
disp('force1 complete');
%}

k1 = [1.5 1.75 2 2.5];
F1 = 0;
ke = linspace(1.5,4,500);
F = linspace(0,2,500);
noiselevel = [0.05 0.1 0.2];
t=linspace(0,1e6,5e6);
for i = 101:200
    for j = 1:length(k1)
        for k = 1:length(noiselevel)
            [XdetF{i,j,k}, XstoF{i,j,k}] = hbtoymodel(F(i),k1(j),noiselevel(k),0,2,t);    % force scans
        end
    end
end
save('/Users/joshsalvi/Desktop/toymodel/HBmodel-force-5e6-2.mat','-v7.3')
clear
disp('force2 complete');
%}

k1 = [1.5 1.75 2 2.5];
F1 = 0;
ke = linspace(1.5,4,500);
F = linspace(0,2,500);
noiselevel = [0.05 0.1 0.2];
t=linspace(0,1e6,5e6);
disp('force2 complete');
for i = 201:300
    for j = 1:length(k1)
        for k = 1:length(noiselevel)
            [XdetF{i,j,k}, XstoF{i,j,k}] = hbtoymodel(F(i),k1(j),noiselevel(k),0,2,t);    % force scans
        end
    end
end
save('/Users/joshsalvi/Desktop/toymodel/HBmodel-force-5e6-3.mat','-v7.3')
clear
%}

k1 = [1.5 1.75 2 2.5];
F1 = 0;
ke = linspace(1.5,4,500);
F = linspace(0,2,500);
noiselevel = [0.05 0.1 0.2];
t=linspace(0,1e6,5e6);
disp('force3 complete');
for i = 301:400
    for j = 1:length(k1)
        for k = 1:length(noiselevel)
            [XdetF{i,j,k}, XstoF{i,j,k}] = hbtoymodel(F(i),k1(j),noiselevel(k),0,2,t);    % force scans
        end
    end
end
save('/Users/joshsalvi/Desktop/toymodel/HBmodel-force-5e6-4.mat','-v7.3')
clear
%}

k1 = [1.5 1.75 2 2.5];
F1 = 0;
ke = linspace(1.5,4,500);
F = linspace(0,2,500);
noiselevel = [0.05 0.1 0.2];
t=linspace(0,1e6,5e6);
disp('force4 complete');
for i = 401:500
    for j = 1:length(k1)
        for k = 1:length(noiselevel)
            [XdetF{i,j,k}, XstoF{i,j,k}] = hbtoymodel(F(i),k1(j),noiselevel(k),0,2,t);    % force scans
        end
    end
end
save('/Users/joshsalvi/Desktop/toymodel/HBmodel-force-5e6-5.mat','-v7.3')
clear
%}

k1 = [1.5 1.75 2 2.5];
F1 = 0;
ke = linspace(1.5,4,500);
F = linspace(0,2,500);
noiselevel = [0.05 0.1 0.2];
t=linspace(0,1e6,5e6);
disp('Finished force scans');
for i = 1:100
    for j = 1:length(F1)
        for k = 1:length(noiselevel)
            [Xdetk{i,j,k}, Xstok{i,j,k}] = hbtoymodel(F1(j),ke(i),noiselevel(k),0,2,t);    % stiffness scans
        end
    end
end
save('/Users/joshsalvi/Desktop/toymodel/HBmodel-stiffness-5e6-1.mat','-v7.3')
clear
%}

k1 = [1.5 1.75 2 2.5];
F1 = 0;
ke = linspace(1.5,4,500);
F = linspace(0,2,500);
noiselevel = [0.05 0.1 0.2];
t=linspace(0,1e6,5e6);
disp('stiffness1 complete');
for i = 101:200
    for j = 1:length(F1)
        for k = 1:length(noiselevel)
            [Xdetk{i,j,k}, Xstok{i,j,k}] = hbtoymodel(F1(j),ke(i),noiselevel(k),0,2,t);    % stiffness scans
        end
    end
end
save('/Users/joshsalvi/Desktop/toymodel/HBmodel-stiffness-5e6-2.mat','-v7.3')
clear
%}

k1 = [1.5 1.75 2 2.5];
F1 = 0;
ke = linspace(1.5,4,500);
F = linspace(0,2,500);
noiselevel = [0.05 0.1 0.2];
t=linspace(0,1e6,5e6);
disp('stiffness2 complete');
for i = 201:300
    for j = 1:length(F1)
        for k = 1:length(noiselevel)
            [Xdetk{i,j,k}, Xstok{i,j,k}] = hbtoymodel(F1(j),ke(i),noiselevel(k),0,2,t);    % stiffness scans
        end
    end
end
save('/Users/joshsalvi/Desktop/toymodel/HBmodel-stiffness-5e6-3.mat','-v7.3')
clear
%}

k1 = [1.5 1.75 2 2.5];
F1 = 0;
ke = linspace(1.5,4,500);
F = linspace(0,2,500);
noiselevel = [0.05 0.1 0.2];
t=linspace(0,1e6,5e6);
disp('stiffness3 complete');
for i = 301:400
    for j = 1:length(F1)
        for k = 1:length(noiselevel)
            [Xdetk{i,j,k}, Xstok{i,j,k}] = hbtoymodel(F1(j),ke(i),noiselevel(k),0,2,t);    % stiffness scans
        end
    end
end
save('/Users/joshsalvi/Desktop/toymodel/HBmodel-stiffness-5e6-4.mat','-v7.3')
clear
%}

k1 = [1.5 1.75 2 2.5];
F1 = 0;
ke = linspace(1.5,4,500);
F = linspace(0,2,500);
noiselevel = [0.05 0.1 0.2];
t=linspace(0,1e6,5e6);
disp('stiffness4 complete');
for i = 401:500
    for j = 1:length(F1)
        for k = 1:length(noiselevel)
            [Xdetk{i,j,k}, Xstok{i,j,k}] = hbtoymodel(F1(j),ke(i),noiselevel(k),0,2,t);    % stiffness scans
        end
    end
end
save('/Users/joshsalvi/Desktop/toymodel/HBmodel-stiffness-5e6-5.mat','-v7.3')
clear
disp('Finished stiffness scans');
%}

%% SUPERCRITICAL HOPF
mu = linspace(-0.5,0,500);
t=linspace(0,1e4,5e5);
rt=mu;
for i = 1:200
    [XHopfdet005{i} XHopfsto005{i}] = hopfstoch(rt(i),1,0.05,0.05,t);
    [XHopfdet01{i} XHopfsto01{i}] = hopfstoch(rt(i),1,0.1,0.1,t);
    [XHopfdet02{i} XHopfsto02{i}] = hopfstoch(rt(i),1,0.2,0.2,t);
    [XHopfdet04{i} XHopfsto04{i}] = hopfstoch(rt(i),1,0.4,0.4,t);
    [XHopfdet06{i} XHopfsto06{i}] = hopfstoch(rt(i),1,0.6,0.6,t);
    [XHopfdet08{i} XHopfsto08{i}] = hopfstoch(rt(i),1,0.8,0.8,t);
    [XHopfdet1{i} XHopfsto1{i}] = hopfstoch(rt(i),1,1,1,t);
end

Hopfdet{:,1}=XHopfdet005;
Hopfdet{:,2}=XHopfdet01;
Hopfdet{:,3}=XHopfdet02;
Hopfdet{:,4}=XHopfdet04;
Hopfdet{:,5}=XHopfdet06;
Hopfdet{:,6}=XHopfdet08;
Hopfdet{:,7}=XHopfdet1;
Hopfsto{:,1}=XHopfsto005;
Hopfsto{:,2}=XHopfsto01;
Hopfsto{:,3}=XHopfsto02;
Hopfsto{:,4}=XHopfsto04;
Hopfsto{:,5}=XHopfsto06;
Hopfsto{:,6}=XHopfsto08;
Hopfsto{:,7}=XHopfsto1;
Xdet = Hopfdet; Xsto = Hopfsto;
clear Hopfdet Hopfsto
save('/Users/joshsalvi/Desktop/Hopf/Hopfstochoutput-all-1.mat','-v7.3')
%}
disp('Hopf 1');
clear

mu = linspace(-0.5,0,500);
t=linspace(0,1e4,5e5);
rt=mu;
for i = 201:350
    [XHopfdet005{i} XHopfsto005{i}] = hopfstoch(rt(i),1,0.05,0.05,t);
    [XHopfdet01{i} XHopfsto01{i}] = hopfstoch(rt(i),1,0.1,0.1,t);
    [XHopfdet02{i} XHopfsto02{i}] = hopfstoch(rt(i),1,0.2,0.0,2,t);
    [XHopfdet04{i} XHopfsto04{i}] = hopfstoch(rt(i),1,0.4,0.4,t);
    [XHopfdet06{i} XHopfsto06{i}] = hopfstoch(rt(i),1,0.6,0.6,t);
    [XHopfdet08{i} XHopfsto08{i}] = hopfstoch(rt(i),1,0.8,0.8,t);
    [XHopfdet1{i} XHopfsto1{i}] = hopfstoch(rt(i),1,1,1,t);
end

Hopfdet{:,1}=XHopfdet005;
Hopfdet{:,2}=XHopfdet01;
Hopfdet{:,3}=XHopfdet02;
Hopfdet{:,4}=XHopfdet04;
Hopfdet{:,5}=XHopfdet06;
Hopfdet{:,6}=XHopfdet08;
Hopfdet{:,7}=XHopfdet1;
Hopfsto{:,1}=XHopfsto005;
Hopfsto{:,2}=XHopfsto01;
Hopfsto{:,3}=XHopfsto02;
Hopfsto{:,4}=XHopfsto04;
Hopfsto{:,5}=XHopfsto06;
Hopfsto{:,6}=XHopfsto08;
Hopfsto{:,7}=XHopfsto1;
Xdet = Hopfdet; Xsto = Hopfsto;
clear Hopfdet Hopfsto
save('/Users/joshsalvi/Desktop/Hopf/Hopfstochoutput-all-2.mat','-v7.3')
disp('Hopf 2');
clear

mu = linspace(-0.5,0,500);
t=linspace(0,1e4,5e5);
rt=mu;
for i = 351:500
    [XHopfdet005{i} XHopfsto005{i}] = hopfstoch(rt(i),1,0.05,0.05,t);
    [XHopfdet01{i} XHopfsto01{i}] = hopfstoch(rt(i),1,0.1,0.1,t);
    [XHopfdet02{i} XHopfsto02{i}] = hopfstoch(rt(i),1,0.2,0.0,2,t);
    [XHopfdet04{i} XHopfsto04{i}] = hopfstoch(rt(i),1,0.4,0.4,t);
    [XHopfdet06{i} XHopfsto06{i}] = hopfstoch(rt(i),1,0.6,0.6,t);
    [XHopfdet08{i} XHopfsto08{i}] = hopfstoch(rt(i),1,0.8,0.8,t);
    [XHopfdet1{i} XHopfsto1{i}] = hopfstoch(rt(i),1,1,1,t);
end

Hopfdet{:,1}=XHopfdet005;
Hopfdet{:,2}=XHopfdet01;
Hopfdet{:,3}=XHopfdet02;
Hopfdet{:,4}=XHopfdet04;
Hopfdet{:,5}=XHopfdet06;
Hopfdet{:,6}=XHopfdet08;
Hopfdet{:,7}=XHopfdet1;
Hopfsto{:,1}=XHopfsto005;
Hopfsto{:,2}=XHopfsto01;
Hopfsto{:,3}=XHopfsto02;
Hopfsto{:,4}=XHopfsto04;
Hopfsto{:,5}=XHopfsto06;
Hopfsto{:,6}=XHopfsto08;
Hopfsto{:,7}=XHopfsto1;
Xdet = Hopfdet; Xsto = Hopfsto;
clear Hopfdet Hopfsto
save('/Users/joshsalvi/Desktop/Hopf/Hopfstochoutput-all-3.mat','-v7.3')
disp('Hopf 3');
clear

%%
% CUSP/FOLD BIFURCATION
% Set -0.5<b2<0 and b2=0.64633 so that fold bifurcations are at b2 = ±0.05

mu = linspace(-0.5,0,500);
b2=0.2565;
t=linspace(0,1e5,5e6);
rt=mu;
for i=1:150
   rt=mu;
[Xcuspdet01{i}, Xcuspsto01{i}] = cuspstoch(rt(i),b2,0.1,t);
[Xcuspdet02{i}, Xcuspsto02{i}] = cuspstoch(rt(i),b2,0.2,t);
[Xcuspdet04{i}, Xcuspsto04{i}] = cuspstoch(rt(i),b2,0.4,t);
[Xcuspdet005{i}, Xcuspsto005{i}] = cuspstoch(rt(i),b2,0.05,t);
end
noiselevel=[0.05 0.1 0.2 0.4 0.6 0.8 1];


cuspdet{:,1}=Xcuspdet005;
cuspdet{:,2}=Xcuspdet01;
cuspdet{:,3}=Xcuspdet02;
cuspdet{:,4}=Xcuspdet04;
cuspsto{:,1}=Xcuspsto005;
cuspsto{:,2}=Xcuspsto01;
cuspsto{:,3}=Xcuspsto02;
cuspsto{:,4}=Xcuspsto04;
Xdet = cuspdet; Xsto = cuspsto;
clear cuspdet cuspsto
clear Xcuspdet01 Xcuspsto01 Xcuspdet02 Xcuspsto02 Xcuspdet04 Xcuspsto04 Xcuspdet005 Xcuspsto005
save('/Users/joshsalvi/Desktop/cusp/cuspstochoutput-all-1.mat','-v7.3')
disp('cusp 1');
clear
%}

mu = linspace(-0.5,0,500);b2=0.2565;
t=linspace(0,1e5,5e6);
rt=mu;
for i=151:300
   rt=mu;
[Xcuspdet01{i}, Xcuspsto01{i}] = cuspstoch(rt(i),b2,0.1,t);
[Xcuspdet02{i}, Xcuspsto02{i}] = cuspstoch(rt(i),b2,0.2,t);
[Xcuspdet04{i}, Xcuspsto04{i}] = cuspstoch(rt(i),b2,0.4,t);
[Xcuspdet005{i}, Xcuspsto005{i}] = cuspstoch(rt(i),b2,0.05,t);
end
noiselevel=[0.05 0.1 0.2 0.4 0.6 0.8 1];


cuspdet{:,1}=Xcuspdet005;
cuspdet{:,2}=Xcuspdet01;
cuspdet{:,3}=Xcuspdet02;
cuspdet{:,4}=Xcuspdet04;
cuspsto{:,1}=Xcuspsto005;
cuspsto{:,2}=Xcuspsto01;
cuspsto{:,3}=Xcuspsto02;
cuspsto{:,4}=Xcuspsto04;
Xdet = cuspdet; Xsto = cuspsto;
clear cuspdet cuspsto
clear Xcuspdet01 Xcuspsto01 Xcuspdet02 Xcuspsto02 Xcuspdet04 Xcuspsto04 Xcuspdet005 Xcuspsto005
save('/Users/joshsalvi/Desktop/cusp/cuspstochoutput-all-2.mat','-v7.3')
disp('cusp 2');
clear
%}

mu = linspace(-0.5,0,500);b2=0.2565;
t=linspace(0,1e5,5e6);
rt=mu;
for i=301:400
   rt=mu;
[Xcuspdet01{i}, Xcuspsto01{i}] = cuspstoch(rt(i),b2,0.1,t);
[Xcuspdet02{i}, Xcuspsto02{i}] = cuspstoch(rt(i),b2,0.2,t);
[Xcuspdet04{i}, Xcuspsto04{i}] = cuspstoch(rt(i),b2,0.4,t);
[Xcuspdet005{i}, Xcuspsto005{i}] = cuspstoch(rt(i),b2,0.05,t);
end
noiselevel=[0.05 0.1 0.2 0.4 0.6 0.8 1];


cuspdet{:,1}=Xcuspdet005;
cuspdet{:,2}=Xcuspdet01;
cuspdet{:,3}=Xcuspdet02;
cuspdet{:,4}=Xcuspdet04;
cuspsto{:,1}=Xcuspsto005;
cuspsto{:,2}=Xcuspsto01;
cuspsto{:,3}=Xcuspsto02;
cuspsto{:,4}=Xcuspsto04;

Xdet = cuspdet; Xsto = cuspsto;
clear cuspdet cuspsto
clear Xcuspdet01 Xcuspsto01 Xcuspdet02 Xcuspsto02 Xcuspdet04 Xcuspsto04 Xcuspdet005 Xcuspsto005
save('/Users/joshsalvi/Desktop/cusp/cuspstochoutput-all-3.mat','-v7.3')
disp('cusp 3');
clear
%}

mu = linspace(-0.5,0,500);b2=0.2565;
t=linspace(0,1e5,5e6);
rt=mu;
for i=401:500
   rt=mu;
[Xcuspdet01{i}, Xcuspsto01{i}] = cuspstoch(rt(i),b2,0.1,t);
[Xcuspdet02{i}, Xcuspsto02{i}] = cuspstoch(rt(i),b2,0.2,t);
[Xcuspdet04{i}, Xcuspsto04{i}] = cuspstoch(rt(i),b2,0.4,t);
[Xcuspdet005{i}, Xcuspsto005{i}] = cuspstoch(rt(i),b2,0.05,t);
end
noiselevel=[0.05 0.1 0.2 0.4 0.6 0.8 1];


cuspdet{:,1}=Xcuspdet005;
cuspdet{:,2}=Xcuspdet01;
cuspdet{:,3}=Xcuspdet02;
cuspdet{:,4}=Xcuspdet04;
cuspsto{:,1}=Xcuspsto005;
cuspsto{:,2}=Xcuspsto01;
cuspsto{:,3}=Xcuspsto02;
cuspsto{:,4}=Xcuspsto04;
Xdet = cuspdet; Xsto = cuspsto;
clear cuspdet cuspsto
clear Xcuspdet01 Xcuspsto01 Xcuspdet02 Xcuspsto02 Xcuspdet04 Xcuspsto04 Xcuspdet005 Xcuspsto005
save('/Users/joshsalvi/Desktop/cusp/cuspstochoutput-all-4.mat','-v7.3')
disp('cusp 4');
clear
%% Partitioning data example

noiselevel = [0.05 0.1 0.2 0.4];
for i = 1:4
for j=205:428
Xdet1{i}{j}=Xdet{i}{j}(1:1e5);
Xdet2{i}{j}=Xdet{i}{j}(1e5+1:2e5);
Xdet3{i}{j}=Xdet{i}{j}(2e5+1:3e5);
Xdet4{i}{j}=Xdet{i}{j}(3e5+1:4e5);
Xdet5{i}{j}=Xdet{i}{j}(4e5+1:end);
Xsto1{i}{j}=Xsto{i}{j}(1:1e5);
Xsto2{i}{j}=Xsto{i}{j}(1e5+1:2e5);
Xsto3{i}{j}=Xsto{i}{j}(2e5+1:3e5);
Xsto4{i}{j}=Xsto{i}{j}(3e5+1:4e5);
Xsto5{i}{j}=Xsto{i}{j}(4e5+1:end);
end
end

save('/Users/joshsalvi/Desktop/SNIC/SNIC-divided/SNICstochoutput-fineI-2-pt1.mat','Xdet1','Xsto1','t','mu','I','-v7.3')
save('/Users/joshsalvi/Desktop/SNIC/SNIC-divided/SNICstochoutput-fineI-2-pt2.mat','Xdet2','Xsto2','t','mu','I','-v7.3')
save('/Users/joshsalvi/Desktop/SNIC/SNIC-divided/SNICstochoutput-fineI-2-pt3.mat','Xdet3','Xsto3','t','mu','I','-v7.3')
save('/Users/joshsalvi/Desktop/SNIC/SNIC-divided/SNICstochoutput-fineI-2-pt4.mat','Xdet4','Xsto4','t','mu','I','-v7.3')
save('/Users/joshsalvi/Desktop/SNIC/SNIC-divided/SNICstochoutput-fineI-2-pt5.mat','Xdet5','Xsto5','t','mu','I','-v7.3')
