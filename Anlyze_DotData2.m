

%% Anlyze_DotData2.m

% Alistair Boettiger                                   Date Begun: 03/14/11
% Levine Lab                                        Last Modified: 03/14/11

%% Description
% Load output of Data files from Unsupervised_DotFinding and compare
% different embryos
%

clear all;

folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
% fname =   'MP10_22C_sna_y_c';  %
fname =  'MP05_22C_sna_y'; 
load([folder,fname,'_slidedata'], 'Data'); 

% load([folder,fname,'_slidedata_3te'], 'Data'); 
% load([folder,'MP05_data_b'], 'Data');

% Quick look at data distributions 

N = length(Data); 
ave = zeros(N,2);
stdev = zeros(N,2); 
ons = cell(N,2);
offs = cell(N,2);

 h = 2048; 

t = .5; 
spread = 1.4;

   % figure(1); clf;
    for chn = 1:2  
            for i=1:N
                %try
                    if i<10
                        emb = ['0',num2str(i)];
                    else
                        emb = num2str(i);
                    end
                    
                    mRNAsadj = Data{i,chn}.mRNAsadj;
                    load([folder,fname,'_',emb,'_nucdata.mat']); 
                    
                       hn = size(NucLabeled,1);  % size of rescaled nuclear image
                       NucLabel = imresize(NucLabeled,.5,'nearest'); % upscale NucLabeled to resolution of mRNA chanel;  
                    
                    Nnucs = max(NucLabel(:)); 
                    Plot_mRNA = NucLabel;  
                     for n=1:Nnucs; % i 
                         Plot_mRNA(NucLabel==n) = mRNAsadj(n);
                     end
                        
                    % subplot(2,N,N*(chn-1)+i);
                     [ons{i,chn},offs{i,chn}] = fxn_regionvar(NucLabel,Plot_mRNA,mRNAsadj,0,spread,Nnucs);
                    
                   
%                      figure(1);
%                          subplot(2,N,N*(chn-1)+i);
%                         hist(mRNAsadj,linspace(0,400,30)); xlim([0,400]);
%                         ave(i,chn) = nanmean(mRNAsadj);
%                         stdev(i,chn) = nanstd(mRNAsadj);
%                         title([ 'mean = ' num2str( ave(i,chn),3), '  std = ',num2str(stdev(i,chn),3) ] ) ; 
%                 catch err
%                     disp(err.message); 
%                     break
%                 end
            end
    end
           
   ave(:,1)./ave(:,2)
   stdev(:,1)./stdev(:,2)
        
   
   
        
   
             saveas(Fig_regvar,[folder,fname,'_',emb,'_chn',num2str(mRNAchn),'rvar.fig']); 
   
        
        