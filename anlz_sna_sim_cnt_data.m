
%%  anlz_sna_sim_cnt_data

% Alistair Boettiger                                   Date Begun: 07/17/11
% Levine Lab                                        Last Modified: 07/20/11


%% 

clear all;
  folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/';
  rawfolder = '/Volumes/Data/Lab Data/Raw_Data/';


%%

% working off of data from anlz_sna_count_data

date = '2011-04_and_earlier/'; % 
fname ='MP12Hz_snaD_22C_b'; ver = ''; vout = '';
load([folder,date,fname,ver,'_slidedata',vout,'.mat'],'data'); 

figure(1); clf; colordef white; set(gcf,'color','w');
figure(11); clf; colordef white; set(gcf,'color','w');
for e = 1:length(data);
    sna = data{e}.Data_sort(:,2) - min(data{e}.Data_sort(:,2));
    sim = data{e}.Data_sort(:,3);
    figure(1); subplot(2,5,e); 
    scatter(sna,sim); xlim([0,250]);
    
    figure(11); subplot(2,5,e);
    plot(data{e}.Data_sort(:,1),sna,'r.');
    hold on; 
    plot(data{e}.Data_sort(:,1),sim,'b.');
end

%% sna-sim heterogeniety in sna2.8 Hets
  

  slidedate = '2011-04_and_earlier/'; % '2011-05-22/'; %    '2011-06-20/';%    
  subfolder =  'sna2p8het/'; %  's05_MP06/'   ; % 's04_MP10/';%  's06_MP10_sna18/'; %   's07_MP08/' ;%  'sna2.8Hz/';%   'MP07Hz/';% 's21_MP07/';  % '' %    's05_MP06/'; % 's11_G4B/'; %  % s04_MP10/';%    
  fname =  'sna2p8het_sim_sna'; % 's05_MP06Hz_b'; ver = '_v4';%  's04_MP10Hz_b'; ver = '_v2'; % 's06_MP10_sna18_b'; ver = '_v4'; cbar =1;  %     % 's07_MP08Hz_snaD_22C';  % 'sna2.8Hz_snaD_22C';st_channel = 1; % 'MP07Hz_snaD_22C'; % 'MP07het_snaD_22C'; Es=4; %  'MP10_22C_sna_y_d'; ver = '_v3'; % 'MP05_22C_sna_y_c';  ver = ''; %  % 's07_MP05Hz_22C';%  's05_MP06Hz'; ver = '_v2';%   %      's06_MP10_sna18_b'; st_channel = 1;   ver = '_v4'; % ' 's05_MP06Hz_b'; ver = '_v2';  % 's11_G4B_LacZ'; ver = '_v2'; legon =0; %     'MP07het_snaD_22C';%   %  '_v2';  %
 ver = '_v2';

 % work off of raw data counts
 emb = '01';

       load([rawfolder,slidedate,subfolder,fname,'_',emb,'_nucdata.mat']);  


          load([folder,slidedate,fname,'_',emb,'_chn', num2str( st_channel+1),'_data',ver,'.mat']); 
          mRNAsadj = mRNA_sadj; % mRNA_cnt./nuc_area;
          disp(['chn1 thresh: ', num2str(Rpars.min_int)]);
          ipars{1} = Rpars; 

            load([folder,slidedate,fname,'_',emb,'_chn',num2str( st_channel+2),'_data',ver,'.mat']);
             disp(['chn2 thresh: ', num2str(Rpars.min_int)]);
            mRNAsadj2 = mRNA_sadj; % mRNA_cnt./nuc_area;  % 
            ipars{2} = Rpars; 

    
          
    
%%              
       PlotmRNA = imresize(NucLabeled,.5,'nearest');
       PlotmRNA2 = PlotmRNA; 

       
       NucLabeled = imresize(NucLabeled,.5,'nearest'); 
                      Nnucs =    max( NucLabeled(:) );
                      for n=1:Nnucs;
                          PlotmRNA(PlotmRNA==n) = mRNAsadj(n+cor);
                          PlotmRNA2(PlotmRNA2==n) = mRNAsadj2(n+cor);
                            
                      end
         sim_max = 120;   
         GrHot = zeros(sim_max,3); 
     for n=1:sim_max
         GrHot(n,:) = [ (n/sim_max).^3, min(1.5*n/sim_max,1), max((-.25+n/sim_max),0).^.5]; 
     end              
                       
    figure(1); clf; imagesc(PlotmRNA); colormap(GrHot);  colorbar; caxis([40,sim_max]);
    figure(2); clf; imagesc(PlotmRNA2); colormap hot;  colorbar; caxis([40,170]);
    
       
    figure(3); clf; plot(PlotmRNA2(:),PlotmRNA(:),'ko'); 
    xlabel('snail  (mRNA counts)');
    ylabel('sim (mRNA counts)'); 
    xlim([40,200]); ylim([40,160]);
    
    sna = PlotmRNA2-30;
    sna(sna<90) = 0; 
    sim = PlotmRNA-30;
    sim(sna>90) = 0; 
    sim(sim<0) = 0; 
       
    [h,w] = size(PlotmRNA);
    I = zeros(h,w,3,'uint16');
    I(:,:,1) = makeuint(sna,16);
    
     I(:,:,2) = makeuint(3*sim,16);
    figure(4); clf; imagesc(I); set(gcf,'color','k'); 
    axis off;
    
    
    Red = [ (1:190)'./190,zeros(190,1), zeros(190,1)];
    Green = [ zeros(120,1),(1:120)'./120, zeros(120,1)];
    
    figure(1); colordef white;  clf; imagesc(sna); colormap(Red); colorbar; 
    figure(2); colordef white; clf; imagesc(sim); colormap(Green); colorbar; colordef white;
    
     sna = PlotmRNA2-30;
      sim = PlotmRNA-30;
      sna(sna<40 & sim<40) = 0;
      sim(sna<40 & sim<40) = 0;
      sna = nonzeros(sna);
      sim = nonzeros(sim);
     figure(3); clf; scatter((sna(:)),(sim(:)),'ko'); 
    xlabel('snail  (mRNA counts)','FontSize',14);
    ylabel('sim (mRNA counts)','FontSize',14);
    xlim([1,160]); ylim([1,120]);
  lsline; set(gca,'FontSize',14);
 
  [p,S] =  polyfit(sim(:),sna(:),1)
  
    corrcoef(sim(:),sna(:))
    
    figure(5); clf; scatterhist(sna(:),sim(:));
    
    %% cumulative probability of sim given snail
    
   
    figure(5); clf; x = linspace(0,200,20);
    subplot(3,1,1); hist(mRNAsadj(mRNAsadj2>150),x); xlim([0,200]);
    subplot(3,1,2); hist(mRNAsadj(mRNAsadj2<150 & mRNAsadj2>120),x); xlim([0,200]);
    subplot(3,1,3); hist(mRNAsadj(mRNAsadj2<120 & mRNAsadj2>90),x); xlim([0,200]);
    
    p = zeros(1,200);
    for s = 1:200
        p(s) = sum(mRNAsadj(mRNAsadj2>s))/Nnucs;
    end     
    figure(5); clf; plot(p); 
    
    