

%% Anlyze_DotData2.m

% Alistair Boettiger                                   Date Begun: 03/14/11
% Levine Lab                                        Last Modified: 03/14/11

%% Description
% Load output of Data files from Unsupervised_DotFinding and compare
% different embryos
%

clear all;

  t = .5; 
spread = 1.4;

folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
% fname =   'MP10_22C_sna_y_c';  load([folder,fname,'_Slidedata'], 'Data');  %
%fname =  'MP05_22C_sna_y'; load([folder,fname,'_Slidedata'], 'Data'); 

cor = 0; 

 fname =  'MP05_22C_sna_y_c';  cor = 1;  emb = '08';  %  emb = '04'; %   emb = '07'; 

 % fname =  'MP05_22C_sna_y_c'; 
 % fname =  'MP10_22C_sna_y_c';  emb = '01'; % emb = '02';% emb = '03';%  emb = '04';%  emb = '05'; % emb = '07';  %** emb = '06';%**  %   emb = '08'; %


chn = 2; i = str2double(emb);

load([folder,fname,'_slidedata'], 'Data'); 




load([folder,fname,'_',emb,'_nucdata.mat']); 
   
              mRNAsadj = Data{i,chn}.mRNAsadj;

  PlotmRNA = imresize(NucLabeled,.5,'nearest');
  NucLabel = imresize(NucLabeled,.5,'nearest'); 
                      Nnucs =    max( NucLabeled(:) );
                      for n=1:Nnucs;
                          PlotmRNA(PlotmRNA==n) = mRNAsadj(n+cor);
                      end
                     
                      chn_names = {'snail';'reporter'};
                      
                      figure(1); clf; colordef black;
                      imagesc(PlotmRNA); colormap hot; colorbar;
                      set(gcf,'color','k'); axis off;
                      title([chn_names{chn},' ',texlabel(fname,'literal')],'FontSize',14);
                      set(gca,'FontSize',14);
          
  [ons,offs] = fxn_regionvar(NucLabel,PlotmRNA,mRNAsadj,0,spread,Nnucs);
   

  

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

t = .5; 
spread = 1.4;

    figure(1); clf;
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
 %                % ------------------------
%                       Nend = max(NucLabeled(:)); % total nuclei 
%                       Nmin = single(NucLabeled); 
%                       Nmin(Nmin==0)=NaN; 
%                       Nstart = min(Nmin(:)); 
%                       Nucs_list = nonzeros(unique(NucLabeled));
%                       Nnucs = length(Nucs_list);
%                    
%                     hn = size(NucLabeled,1);  % size of rescaled nuclear image
%                     NucLabel = imresize(NucLabeled,h/hn,'nearest'); % upscale NucLabeled to resolution of mRNA chanel;  
%                     
%                
%                     inds = floor(Data{i,chn}.dotC(:,2))+floor(Data{i,chn}.dotC(:,1))*h;   
%                     inds(inds>w*h) = w*h;   
%                        
%                     mRNA_cnt = zeros(1,Nnucs); % store counts of mRNA per cell  
%                     mRNA_den = zeros(1,Nnucs);  % store densities of mRNA per cell
%                     nuc_area = zeros(1,Nnucs); 
%                     Plot_mRNA = single(NucLabel);
%                   
%                     for n=1:Nnucs; % i = 4
%                         nn = Nucs_list(n);
%                         imdata.Area(n) = length(find(NucLabel==nn));
%                         imdata.PixelID{n} = find(NucLabel==nn);
%                         mRNA_cnt(n) = length(intersect(imdata.PixelID{n},inds));
%                         mRNA_den(n) = mRNA_cnt(n)/imdata.Area(n); 
%                         nuc_area(n) = length(imdata.PixelID{n});
%                         Plot_mRNA(NucLabel==nn) = single(mRNA_den(n));
%                     end
%                     
%                      speed_scale = 1; 
%                      mRNA_sadj = mRNA_den*mean(imdata.Area);
%                      NucLabeled = imresize(NucLabel,speed_scale,'nearest');  
%                      Plot_mRNA = imresize(Plot_mRNA,speed_scale,'nearest'); 
%                      
%                        Data{i,chn}.mRNAsadj = mRNA_sadj;
%                        Data{i,chn}.NucLabeled = NucLabeled;
%                        Data{i,chn}.PlotmRNA = Plot_mRNA; 
%               % ----------------------------  %
                try
                     Nnucs =    max( Data{i,chn}.NucLabeled(:) );
                catch
                       PlotmRNA = imresize(NucLabeled,.5,'nearest');
                      Nnucs =    max( NucLabeled(:) );
                      for n=1:Nnucs;
                          PlotmRNA(PlotmRNA==n) = mRNAsadj(n+1);
                      end
                end
                
                     [ons{i,chn},offs{i,chn}] = fxn_regionvar(Data{i,chn}.NucLabeled,Data{i,chn}.PlotmRNA,mRNAsadj,0,spread,Nnucs);
                    
                   
%                       figure(1);
%                           subplot(2,N,N*(chn-1)+i);
%                          hist(mRNAsadj,linspace(0,400,30)); xlim([0,400]);
%                          ave(i,chn) = nanmean(mRNAsadj);
%                          stdev(i,chn) = nanstd(mRNAsadj);
%                          title([ 'mean = ' num2str( ave(i,chn),3), '  std = ',num2str(stdev(i,chn),3) ] ) ; 
           disp(i);
            end
 
    end
           
    Ons_mean = cellfun(@mean,ons)
    Ons_std = cellfun(@var,ons)
    
    Ons_mean(:,1)./Ons_mean(:,2)
    Ons_std(:,1)./Ons_std(:,2)
    
    Ons_std(:,2)./Ons_mean(:,2)
    Ons_std(:,1)./Ons_mean(:,1)
    
    
    % save([folder,fname,'_MP10_ons'], 'ons'); 
    
%    ave(:,1)./ave(:,2)
%    stdev(:,1)./stdev(:,2)
        
   
%save([folder,fname,'_Slidedata'], 'Data'); 
        
   
   %          saveas(Fig_regvar,[folder,fname,'_',emb,'_chn',num2str(mRNAchn),'rvar.fig']); 
   
        
        