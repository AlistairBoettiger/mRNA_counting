

%% Anlyze_DotData2.m

% Alistair Boettiger                                   Date Begun: 03/14/11
% Levine Lab                                        Last Modified: 03/14/11

%% Description
% Load output of Data files from Unsupervised_DotFinding and compare
% different embryos
%

clear all;

rawfolder = '/Volumes/Data/Lab Data/Raw_Data/02-17-11/YW_ths_sog/';
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
fname = 'YW_ths_sog';
% load([folder,fname,'_slidedata_9-12'], 'Data'); 
% adddata = Data ; clear Data; 

load([folder,fname,'_slidedata'], 'Data'); 

% Data = [Data; cell(2,2)];
% Data(9:12,:) = adddata(9:12,:);
% save([folder,fname,'_slidedata'], 'Data'); 

%%



%load([folder,fname,'_Slidedata'], 'Data');  % fixed mRNA sadj

% load([folder,fname,'_slidedata_3te'], 'Data'); 
% load([folder,'MP05_data_b'], 'Data');

% Quick look at data distributions 

N = length(Data); 
ave = zeros(N,2);
stdev = zeros(N,2); 
ons = cell(N,2);
offs = cell(N,2);

 h = 2048; w = 2048;

t1 = .5; 
spread = 1.4;
minObj = 10000; 
stre = 25; 
cor = 1;


    figure(1); clf;
    for i=1:N   
        for chn = 1:2 
                try
                    if i<10
                        emb = ['0',num2str(i)];
                    else
                        emb = num2str(i);
                    end            
                    
                    mRNAsadj = Data{i,chn}.mRNAsadj;
                catch err
                    % disp(err.message); 
                end
                
                       
                load([rawfolder,fname,'_',emb,'_nucdata.mat']);                                 
              

                try  
                       PlotmRNA = imresize(NucLabeled,.5,'nearest');
                       NucLabel = imresize(NucLabeled,.5,'nearest'); 
                      Nnucs =    max( NucLabel(:) );
                      for n=1:Nnucs;
                          PlotmRNA(PlotmRNA==n) = mRNAsadj(n+cor);
                      end
                      [ons{i,chn},offs{i,chn}] = fxn_regionvar(NucLabel,PlotmRNA,mRNAsadj,t1,spread,Nnucs,minObj,stre);
                catch err 
                     disp(err.message); 
                end
        end
    end
  
  %% 
 
    x = linspace(0,450,50); 
    figure(1); clf; colordef black; set(gcf,'color','k');
    for e=1:12
        figure(1); subplot(6,2,e); 
         hist(ons{e,2},x);
         set(findobj(gca,'Type','patch'),'FaceColor','g');
         hold on; hist(ons{e,1},x); xlim([0,450]);
         title(['sog=',num2str(mean(ons{e,2}),3),...
             ' (',num2str(std(ons{e,2}),3) ')  ths=',num2str(mean(ons{e,1}),3),...
             ' (',num2str(std(ons{e,1}),3), ')']);
  
    end
    
     x = linspace(0,450,50); 
    figure(2); clf; colordef black; set(gcf,'color','k');
    for e=1:12
        figure(2); subplot(6,2,e); 
         hist(offs{e,2},x);
         set(findobj(gca,'Type','patch'),'FaceColor','g');
         hold on; hist(offs{e,1},x); xlim([0,450]);
         title(['sog=',num2str(mean(offs{e,2}),3),...
             ' (',num2str(std(offs{e,2}),3) ')  ths=',num2str(mean(offs{e,1}),3),...
             ' (',num2str(std(offs{e,1}),3), ')']);
  
    end
   
   