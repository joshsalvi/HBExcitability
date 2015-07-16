Nboot=2; 

%load('/Volumes/Promise Pegasus/Manual Backup/Simulation Data/toymodel/HBmodel-force-divided/HBmodel-force-all-pt2-dwnspl1-thresh2.7         2.9         3.1-biftype4-stiffind1.mat')

for j = 1:length(Xdet)
    if mod(j,50)==0
    disp(['j = ' num2str(j) '/' num2str(length(Xdet))]);end

for k = 1:length(Xdet{1})
    for l = 1:length(Xdet{1}{1})
            [dipsto(j,k,l), punisto(j,k,l), Xlowsto(j,k,l), Xupsto(j,k,l)]=HartigansDipSignifTest(Xsto{j}{k}{l},Nboot);
            [dipdet(j,k,l), punidet(j,k,l), Xlowdet(j,k,l), Xupdet(j,k,l)]=HartigansDipSignifTest(Xdet{j}{k}{l},Nboot);
    end
end
end
 
