

%% Anlyze_DotData2.m

% Alistair Boettiger                                   Date Begun: 03/15/11
% Levine Lab                                        Last Modified: 03/15/11

%% Description
% Load output of Data files from Unsupervised_DotFinding and compare
% different embryos for Mat-alpha-Gal4 1x vs 2x.  
%

clear all;

folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
fname =  'MGa2x_LacZ_sna' %  'MGa_LacZ'; %  
load([folder,fname,'_slidedata'], 'Data'); 


% Quick look at data distributions 

N = length(Data); 
ave = NaN*zeros(N,2);
stdev = NaN*zeros(N,2); 
ons = cell(N,2);
offs = cell(N,2);

 h = 2048; w = 2048;

t = .5; 
spread = 1.4;

    figure(2); clf;
    for chn = 1:2  
            for i=1:N
                try
                    if i<10
                        emb = ['0',num2str(i)];
                    else
                        emb = num2str(i);
                    end
                    
                    mRNAsadj = Data{i,chn}.mRNAsadj;
                    load([folder,fname,'_',emb,'_nucdata.mat']); 
                    
                 catch err
                     disp(err.message); 
                     break
                end

                      PlotmRNA = imresize(NucLabeled,.5,'nearest');
                      Nnucs =    max( NucLabeled(:) );
                      for n=1:Nnucs;
                          PlotmRNA(PlotmRNA==n) = mRNAsadj(n+1);
                      end
                          
                      
%                      [ons{i,chn},offs{i,chn}] = fxn_regionvar(Data{i,chn}.NucLabeled,Data{i,chn}.PlotmRNA,mRNAsadj,0,spread,Nnucs);
                    
                   
                      figure(2);
                          subplot(2,N,N*(chn-1)+i);
                         hist(mRNAsadj,linspace(0,400,30)); xlim([0,400]);
                         ave(i,chn) = nanmean(mRNAsadj);
                         stdev(i,chn) = nanstd(mRNAsadj);
                         title([ 'mean = ' num2str( ave(i,chn),3), '  std = ',num2str(stdev(i,chn),3) ] ) ; 
                         
                         figure(3);   subplot(2,N,N*(chn-1)+i);
                         imagesc(PlotmRNA); axis off; set(gcf,'color','k');
        %   disp(i);
            end
 
    end
           

        
        