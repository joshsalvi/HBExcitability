%% subcritical hopf
in2=input('Did you check initialization? (1=yes,0=no).......');
if in2 == 0
    break;
end


for m = 1:5
    if m == 1
               file2=['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-pt' num2str(m) '.mat']; 
mat2=matfile(file2,'Writable',true);
mat2.Xdet=cell(500,4);
mat2.Xsto=cell(500,4);
for n=1:4
    disp(['Seg=' num2str(m) ' iter=' num2str(n)]);
file= ['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-all-' num2str(n) '.mat'];
mat1=matfile(file);
if n==1
     for j = 1:4
         xdet=mat1.Xdet(1,j);
         xsto=mat1.Xsto(1,j);
 for i=1:150
   
            
            xdet2{1}=xdet{1}{i}(1,1:1e5);
            xsto2{1}=xsto{1}{i}(1,1:1e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
    end
    end
elseif n==2
     for j = 1:4
         xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
     for i=151:300
   
            
            xdet2{1}=xdet{1}{i}(1,1:1e5);
            xsto2{1}=xsto{1}{i}(1,1:1e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
elseif n == 3
     for j = 1:4
         xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
    for i=301:400
   
            
            xdet2{1}=xdet{1}{i}(1,1:1e5);
            xsto2{1}=xsto{1}{i}(1,1:1e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n ==4
     for j = 1:4
          xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
      for i=401:500
   
           
            xdet2{1}=xdet{1}{i}(1,1:1e5);
            xsto2{1}=xsto{1}{i}(1,1:1e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
end
end
    elseif m == 2
               file2=['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-pt' num2str(m) '.mat']; 
mat2=matfile(file2,'Writable',true);
mat2.Xdet=cell(500,4);
mat2.Xsto=cell(500,4);
for n=1:4
    disp(['Seg=' num2str(m) ' iter=' num2str(n)]);
file= ['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-all-' num2str(n) '.mat'];
mat1=matfile(file);
if n==1
     for j = 1:4
          xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
 for i=1:150
   
           
            xdet2{1}=xdet{1}{i}(1,1+1e5:2e5);
            xsto2{1}=xsto{1}{i}(1,1+1e5:2e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n==2
     for j = 1:4
          xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
     for i=151:300
   
           
            xdet2{1}=xdet{1}{i}(1,1+1e5:2e5);
            xsto2{1}=xsto{1}{i}(1,1+1e5:2e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
elseif n == 3
      for j = 1:4
           xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
    for i=301:400
  
           
            xdet2{1}=xdet{1}{i}(1,1+1e5:2e5);
            xsto2{1}=xsto{1}{i}(1,1+1e5:2e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n ==4
     for j = 1:4
           xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
      for i=401:500
   
          
            xdet2{1}=xdet{1}{i}(1,1+1e5:2e5);
            xsto2{1}=xsto{1}{i}(1,1+1e5:2e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
end
end

    elseif m == 3
               file2=['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-pt' num2str(m) '.mat']; 
mat2=matfile(file2,'Writable',true);
mat2.Xdet=cell(500,4);
mat2.Xsto=cell(500,4);
for n=1:4
    disp(['Seg=' num2str(m) ' iter=' num2str(n)]);
file= ['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-all-' num2str(n) '.mat'];
mat1=matfile(file);
if n==1
     for j = 1:4
           xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
 for i=1:150
   
          
            xdet2{1}=xdet{1}{i}(1,1+2e5:3e5);
            xsto2{1}=xsto{1}{i}(1,1+2e5:3e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n==2
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
     for i=151:300
   
            xdet2{1}=xdet{1}{i}(1,1+2e5:3e5);
            xsto2{1}=xsto{1}{i}(1,1+2e5:3e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
elseif n == 3
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
    for i=301:400
   
            xdet2{1}=xdet{1}{i}(1,1+2e5:3e5);
            xsto2{1}=xsto{1}{i}(1,1+2e5:3e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n ==4
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
      for i=401:500
   
            xdet2{1}=xdet{1}{i}(1,1+2e5:3e5);
            xsto2{1}=xsto{1}{i}(1,1+2e5:3e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
end
end
    elseif m == 4
        file2=['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-pt' num2str(m) '.mat']; 
mat2=matfile(file2,'Writable',true);
mat2.Xdet=cell(500,4);
mat2.Xsto=cell(500,4);
for n=1:4
    disp(['Seg=' num2str(m) ' iter=' num2str(n)]);
file= ['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-all-' num2str(n) '.mat'];
mat1=matfile(file);
if n==1
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
 for i=1:150
    
            xdet2{1}=xdet{1}{i}(1,1+3e5:4e5);
            xsto2{1}=xsto{1}{i}(1,1+3e5:4e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n==2
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
     for i=151:300
   
            xdet2{1}=xdet{1}{i}(1,1+3e5:4e5);
            xsto2{1}=xsto{1}{i}(1,1+3e5:4e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
elseif n == 3
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
    for i=301:400
   
            xdet2{1}=xdet{1}{i}(1,1+3e5:4e5);
            xsto2{1}=xsto{1}{i}(1,1+3e5:4e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n ==4
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
      for i=401:500
    
            xdet2{1}=xdet{1}{i}(1,1+3e5:4e5);
            xsto2{1}=xsto{1}{i}(1,1+3e5:4e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
end
end
    elseif m == 5
file2=['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-pt' num2str(m) '.mat'];
mat2=matfile(file2,'Writable',true);
mat2.Xdet=cell(500,4);
mat2.Xsto=cell(500,4);
for n=1:4
    disp(['Seg=' num2str(m) ' iter=' num2str(n)]);
file= ['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-all-' num2str(n) '.mat'];
mat1=matfile(file);
if n==1
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
 for i=1:150
   
            xdet2{1}=xdet{1}{i}(1,1+4e5:end);
            xsto2{1}=xsto{1}{i}(1,1+4e5:end);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n==2
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
     for i=151:300
    
            xdet2{1}=xdet{1}{i}(1,1+4e5:end);
            xsto2{1}=xsto{1}{i}(1,1+4e5:end);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
elseif n == 3
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
    for i=301:400
   
            xdet2{1}=xdet{1}{i}(1,1+4e5:end);
            xsto2{1}=xsto{1}{i}(1,1+4e5:end);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n ==4
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
      for i=401:500
    
            xdet2{1}=xdet{1}{i}(1,1+4e5:end);
            xsto2{1}=xsto{1}{i}(1,1+4e5:end);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
end
end
    end
end

disp('subcritical Hopf complete');
% supercritical hopf

for m = 1:5
    if m == 1
               file2=['/Volumes/LaCie JDS/Simulation Data/hopfoffset/hopfstochoffsetoutput-pt' num2str(m) '.mat']; 
mat2=matfile(file2,'Writable',true);
mat2.Xdet=cell(500,4);
mat2.Xsto=cell(500,4);
for n=1:3
    disp(['Seg=' num2str(m) ' iter=' num2str(n)]);
file= ['/Volumes/LaCie JDS/Simulation Data/hopfoffset/hopfstochoffsetoutput-all-' num2str(n) '.mat'];
mat1=matfile(file);
if n==1
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
 for i=1:200
   
            xdet2{1}=xdet{1}{i}(1,1:1e5);
            xsto2{1}=xsto{1}{i}(1,1:1e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n==2
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
     for i=201:350
   
            xdet2{1}=xdet{1}{i}(1,1:1e5);
            xsto2{1}=xsto{1}{i}(1,1:1e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
elseif n == 3
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
    for i=351:500
    
            xdet2{1}=xdet{1}{i}(1,1:1e5);
            xsto2{1}=xsto{1}{i}(1,1:1e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n ==4
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
      for i=401:500
   
            xdet2{1}=xdet{1}{i}(1,1:1e5);
            xsto2{1}=xsto{1}{i}(1,1:1e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
end
end
    elseif m == 2
               file2=['/Volumes/LaCie JDS/Simulation Data/hopfoffset/hopfstochoffsetoutput-pt' num2str(m) '.mat']; 
mat2=matfile(file2,'Writable',true);
mat2.Xdet=cell(500,4);
mat2.Xsto=cell(500,4);
for n=1:3
    disp(['Seg=' num2str(m) ' iter=' num2str(n)]);
file= ['/Volumes/LaCie JDS/Simulation Data/hopfoffset/hopfstochoffsetoutput-all-' num2str(n) '.mat'];
mat1=matfile(file);
if n==1
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
 for i=1:200
    
            xdet2{1}=xdet{1}{i}(1,1+1e5:2e5);
            xsto2{1}=xsto{1}{i}(1,1+1e5:2e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n==2
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
     for i=201:350
    
            xdet2{1}=xdet{1}{i}(1,1+1e5:2e5);
            xsto2{1}=xsto{1}{i}(1,1+1e5:2e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
elseif n == 3
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
    for i=351:500
    
            xdet2{1}=xdet{1}{i}(1,1+1e5:2e5);
            xsto2{1}=xsto{1}{i}(1,1+1e5:2e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n ==4
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
      for i=401:500
   
            xdet2{1}=xdet{1}{i}(1,1+1e5:2e5);
            xsto2{1}=xsto{1}{i}(1,1+1e5:2e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
end
end

    elseif m == 3
               file2=['/Volumes/LaCie JDS/Simulation Data/hopfoffset/hopfstochoffsetoutput-pt' num2str(m) '.mat']; 
mat2=matfile(file2,'Writable',true);
mat2.Xdet=cell(500,4);
mat2.Xsto=cell(500,4);
for n=1:3
    disp(['Seg=' num2str(m) ' iter=' num2str(n)]);
file= ['/Volumes/LaCie JDS/Simulation Data/hopfoffset/hopfstochoffsetoutput-all-' num2str(n) '.mat'];
mat1=matfile(file);
if n==1
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
 for i=1:200
    
            xdet2{1}=xdet{1}{i}(1,1+2e5:3e5);
            xsto2{1}=xsto{1}{i}(1,1+2e5:3e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n==2
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
     for i=201:350
    
            xdet2{1}=xdet{1}{i}(1,1+2e5:3e5);
            xsto2{1}=xsto{1}{i}(1,1+2e5:3e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
elseif n == 3
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
    for i=351:500
    
            xdet2{1}=xdet{1}{i}(1,1+2e5:3e5);
            xsto2{1}=xsto{1}{i}(1,1+2e5:3e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n ==4
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
      for i=401:500
    
            xdet2{1}=xdet{1}{i}(1,1+2e5:3e5);
            xsto2{1}=xsto{1}{i}(1,1+2e5:3e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
end
end
    elseif m == 4
        file2=['/Volumes/LaCie JDS/Simulation Data/hopfoffset/hopfstochoffsetoutput-pt' num2str(m) '.mat']; 
mat2=matfile(file2,'Writable',true);
mat2.Xdet=cell(500,4);
mat2.Xsto=cell(500,4);
for n=1:3
    disp(['Seg=' num2str(m) ' iter=' num2str(n)]);
file= ['/Volumes/LaCie JDS/Simulation Data/hopfoffset/hopfstochoffsetoutput-all-' num2str(n) '.mat'];
mat1=matfile(file);
if n==1
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
 for i=1:200
    
            xdet2{1}=xdet{1}{i}(1,1+3e5:4e5);
            xsto2{1}=xsto{1}{i}(1,1+3e5:4e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n==2
     for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
     for i=201:350
   
            xdet2{1}=xdet{1}{i}(1,1+3e5:4e5);
            xsto2{1}=xsto{1}{i}(1,1+3e5:4e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
elseif n == 3
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
    for i=351:500
    
            xdet2{1}=xdet{1}{i}(1,1+3e5:4e5);
            xsto2{1}=xsto{1}{i}(1,1+3e5:4e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n ==4
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
      for i=401:500
    
            xdet2{1}=xdet{1}{i}(1,1+3e5:4e5);
            xsto2{1}=xsto{1}{i}(1,1+3e5:4e5);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
end
end
    elseif m == 5
file2=['/Volumes/LaCie JDS/Simulation Data/subhopfoffset/subhopfstochoffsetoutput-pt' num2str(m) '.mat'];
mat2=matfile(file2,'Writable',true);
mat2.Xdet=cell(500,4);
mat2.Xsto=cell(500,4);
for n=1:3
    disp(['Seg=' num2str(m) ' iter=' num2str(n)]);
file= ['/Volumes/LaCie JDS/Simulation Data/hopfoffset/hopfstochoffsetoutput-all-' num2str(n) '.mat'];
mat1=matfile(file);
if n==1
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
 for i=1:200
    
            xdet2{1}=xdet{1}{i}(1,1+4e5:end);
            xsto2{1}=xsto{1}{i}(1,1+4e5:end);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n==2
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
     for i=201:350
    
            xdet2{1}=xdet{1}{i}(1,1+4e5:end);
            xsto2{1}=xsto{1}{i}(1,1+4e5:end);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
elseif n == 3
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
    for i=351:500
    
            xdet2{1}=xdet{1}{i}(1,1+4e5:end);
            xsto2{1}=xsto{1}{i}(1,1+4e5:end);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
    end
elseif n ==4
    for j = 1:4
            xdet=mat1.Xdet(1,j);
            xsto=mat1.Xsto(1,j);
      for i=401:500
    
            xdet2{1}=xdet{1}{i}(1,1+4e5:end);
            xsto2{1}=xsto{1}{i}(1,1+4e5:end);
mat2.Xdet(i,j)=xdet2;
mat2.Xsto(i,j)=xsto2;
clear xdet2 xsto2
        end
     end
end
end
    end
end

