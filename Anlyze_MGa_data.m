

%% Anlyze_DotData2.m

% Alistair Boettiger                                   Date Begun: 03/15/11
% Levine Lab                                        Last Modified: 03/15/11

%% Description
% Load output of Data files from Unsupervised_DotFinding and compare
% different embryos for Mat-alpha-Gal4 1x vs 2x.  
%

clear all;



%% combine data from multiple slides
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
%     rawfolder =  '/Volumes/Data/Lab Data/Raw_Data/02-17-11/MGa1x/';
%     fname1 = 'MGa_LacZ';  fname2 = 'MGa1x_LacZ_b' ; combname = 'MGa1x_LacZ_combined';

 rawfolder =  '/Volumes/Data/Lab Data/Raw_Data/02-17-11/MGa2x/';
    fname1 = 'MGa2x_LacZ_sna';  fname2 = 'MGa2x_LacZ_sna_b' ; combname = 'MGa2x_LacZ_combined';

  NDat = []; Dat = []; 
for slide = 1:2
    clear Data; 
    if slide == 1
        fname = fname1;
        load([folder,fname,'_slidedata'], 'Data'); 
    elseif slide == 2
    fname = fname2; 
    load([folder,fname,'_slidedata'], 'Data'); 
    end
  
        NData= cell(15,1);
        k = 0; 
    for i = 1:15
        k = k + 1;
        if i<10
            emb = ['0',num2str(i)];
        else
            emb = num2str(i);
        end

        try
            load([rawfolder,fname,'_',emb,'_nucdata.mat']);
           catch de
               disp(de.message);
            try
            load([folder,fname,'_',emb,'_nucdata.mat']);
            catch de
                disp(de.message);
                disp('giving up');
                break
            end
        end
          NData{i} = NucLabeled;
    end
NDat =  [NDat; NData(1:k-1)];   
Dat =   [Dat; Data(1:k-1,:)] ;    
    
end


% save([folder, combname,'_slidedata'], 'Dat','NDat'); 
%  load([folder, combname,'_slidedata'], 'Dat','NDat');


%%  All embryos

% Quick look at data distributions 

cor = 1; % offset correction factor between NucLabeled and mRNAsadj indices

N = 24; 

 h = 2048; w = 2048;

t = .5; 
spread = 1.4;

for slide = 1:2
    clear NDat Dat; 
    if slide == 1;
        fname =  'MGa1x_LacZ_combined';  load([folder,fname,'_slidedata']); 
    elseif slide == 2; 
        fname = 'MGa2x_LacZ_combined'; load([folder,fname,'_slidedata']); 
    end
    
    
    ave = NaN*zeros(N,2);
    stdev = NaN*zeros(N,2); 
    ons = cell(N,2);
    offs = cell(N,2);


    figure(2); clf;
    for chn = 1:2  
            for i=1:N
                try
                    if i<10
                        emb = ['0',num2str(i)];
                    else
                        emb = num2str(i);
                    end
                    
                    mRNAsadj = Dat{i,chn}.mRNAsadj;
                    
                catch err
                     disp(err.message); 
                     break
                end

                           
                 PlotmRNA = imresize(NDat{i},.5,'nearest');
                 NucLabel = imresize(NDat{i},.5,'nearest'); 
                 Nnucs =    max( NucLabel(:) );
                 for n=1:Nnucs;
                      PlotmRNA(PlotmRNA==n) = mRNAsadj(n+cor);
                 end
                   
                     ave(i,chn) = nanmean(mRNAsadj);
                     stdev(i,chn) = nanstd(mRNAsadj);
                 
%                   figure(2);
%                      subplot(2,N,N*(chn-1)+i);
%                      hist(mRNAsadj,linspace(0,400,30)); xlim([0,400]);
%                      title([ 'mean = ' num2str( ave(i,chn),3), '  std = ',num2str(stdev(i,chn),3) ] ) ; 
% 
%                  figure(3);   
%                     subplot(2,N,N*(chn-1)+i);
%                     imagesc(PlotmRNA); axis off; set(gcf,'color','k');
            end           
    end
    
    if slide == 1
        MGa1x_ave = ave;
        MGa1x_std = stdev;
        clear ave stdev;
    elseif slide == 2
        MGa2x_ave = ave;
        MGa2x_std = stdev;
        clear ave stdev;
    end
end      
     
%%



figure(1); clf; scatter(MGa2x_ave(:,1),MGa2x_std(:,1),'ro');
hold on; scatter(MGa1x_ave(:,1),MGa1x_std(:,1),'g.');



        %% individual embryo

emb = '03';  i = str2double(emb); 




       mRNAsadj1 = Data{i,1}.mRNAsadj;
      mRNAsadj2 = Data{i,2}.mRNAsadj;
      
      ms = linspace(0,300,30);
      
      M1 = hist(mRNAsadj1,ms);
      M2 = hist(mRNAsadj2,ms);
      
      figure(1); clf; colordef black; set(gcf,'color','k');
      set(gca,'FontSize',14);
      bar(ms,M1,'r'); hold on; 
      bar(ms,M2,'g'); xlim([-10,280]);

        