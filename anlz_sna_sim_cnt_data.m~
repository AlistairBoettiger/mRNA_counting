
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
 fname ='MP12Hz_snaD_22C_b'; 
ver = ''; vout = '';

% fname = 'sna2p8het_30C';
load([folder,date,fname,ver,'_slidedata',vout,'.mat'],'data'); 

figure(1); clf; colordef white; set(gcf,'color','w');
figure(11); clf; colordef white; set(gcf,'color','w');
for e = 1:length(data);
    sna = data{e}.Data_sort(:,2) - min(data{e}.Data_sort(:,2));
    sim = data{e}.Data_sort(:,3);
    figure(1); subplot(2,5,e);
    
    dist = data{e}.Data_sort(:,1);
    dc = [dist/16E4,zeros(length(dist),1),1-dist/16E4];
    scatter(sna,sim,100,dc,'.');
    xlim([0,250]);
    xlabel('sna mRNA counts'); ylabel('sim mRNA counts');
    
    
    figure(3);  subplot(2,5,e);
    x=data{e}.Data_sort(:,1);
    hist(sna,linspace(0,400,50));  xlim([0,400]);
    
    figure(11); subplot(2,5,e);
    plot(data{e}.Data_sort(:,1),sna,'r.');
    hold on; 
    plot(data{e}.Data_sort(:,1),sim,'b.');
end


%%    Bimodal histogram within mesoderm 

date = '2011-04_and_earlier/'; % 
 fname ='MP12Hz_snaD_22C_b'; 
ver = ''; vout = '';

load([folder,date,fname,ver,'_slidedata',vout,'.mat'],'data'); 

       e =5; xmax = 7.5E4; xmin = 1.8E4;
     %  e = 6; xmax = 7.5E4;
    sna = data{e}.Data_sort(:,2) - min(data{e}.Data_sort(:,2));
    sim = data{e}.Data_sort(:,3);
    x=data{e}.Data_sort(:,1);
    b = linspace(0,350,27);
  
    figure(2); 
      plot(data{e}.Data_sort(:,1),sna,'r.');
    hold on; 
    plot(data{e}.Data_sort(:,1),sim,'b.');
    
        figure(3); clf;   
           
 hist(sna(x>xmin & x<xmax),b);  xlim([-5,260]);    
xlabel('mRNA count'); ylabel('frequency');
  
    
           figure(4); clf;
    hist(sna,b);  xlim([0,300]);
    
%%  plot bimodal histogram  

date = '2011-04_and_earlier/'; % 
ver = ''; vout = '';

 fname = 'sna2p8het_30C';
load([folder,date,fname,ver,'_slidedata',vout,'.mat'],'data'); 
       figure(3); clf;
       e =7; xmax = 9.5E4;
     %  e = 6; xmax = 7.5E4;
    sna = data{e}.Data_sort(:,2) - min(data{e}.Data_sort(:,2));
    sim = data{e}.Data_sort(:,3);
    x=data{e}.Data_sort(:,1);
    b = linspace(0,350,27);
 hist(sna(x<xmax),b);  xlim([-5,260]);    
xlabel('mRNA count'); ylabel('frequency');
  
    
           figure(4); clf;
    hist(sna,b);  xlim([0,300]);
    
       figure(2); 
      plot(data{e}.Data_sort(:,1),sna,'r.');
    hold on; 
    plot(data{e}.Data_sort(:,1),sim,'b.');
    
    %%  plot bimodal histogram   B

date = '2011-04_and_earlier/'; % 
ver = ''; vout = '';

 fname = 'sna2p8het_sim_sna';   ver = '_v2'; % 
load([folder,date,fname,ver,'_slidedata',vout,'.mat'],'data'); 
       figure(3); clf;
       e =1; xmax = 11.7E4; xmin = 1.5E4;
     %  e = 6; xmax = 7.5E4;
    sna = data{e}.Data_sort(:,2) - min(data{e}.Data_sort(:,2));
    sim = data{e}.Data_sort(:,3);
    x=data{e}.Data_sort(:,1);
    b = linspace(0,350,27);
 hist(sna(x>xmin & x<xmax),b);  xlim([-5,260]);    
xlabel('mRNA count'); ylabel('frequency');
  
    
           figure(4); clf;
    hist(sna,b);  xlim([0,300]);
    
       figure(2); 
      plot(data{e}.Data_sort(:,1),sna,'r.');
    hold on; 
    plot(data{e}.Data_sort(:,1),sim,'b.');
    
%% sna-sim heterogeniety in sna2.8 Hets
  
st_channel = 0;  cor = 0;
slidedate = '2011-04_and_earlier/'; % 'sna2p8het_sim_sna';   ver = '_v2'; % 
% % ------------------------------------ % % 
%   
%   subfolder =  'sna2p8het/'; % 
%   fname =  'sna2p8het_30C'; % 'sna2p8het_sim_sna';   ver = '_v2'; % 
% 
%  % work off of raw data counts
% emb = '03' ;  ThS = 125; sna_bkg = 25; sim_bkg = 5;
%  emb = '04';  ThS = 80; sna_bkg = 25; sim_bkg = 5;
% % ------------------------------------ % %  
 
 %
 
  subfolder =  'MP12Hz/'; % 
  fname ='MP12Hz_snaD_22C_b'; 

 % work off of raw data counts
emb = '05' ;  ThS = 100; sna_bkg = 30; sim_bkg = 30;

 
%
 

       load([rawfolder,slidedate,subfolder,fname,'_',emb,'_nucdata.mat']);  


          load([folder,slidedate,fname,'_',emb,'_chn', num2str( st_channel+1),'_data',ver,'.mat']); 
          mRNAsadj = mRNA_sadj; % mRNA_cnt./nuc_area;
          disp(['chn1 thresh: ', num2str(Rpars.min_int)]);
          ipars{1} = Rpars; 

            load([folder,slidedate,fname,'_',emb,'_chn',num2str( st_channel+2),'_data',ver,'.mat']);
             disp(['chn2 thresh: ', num2str(Rpars.min_int)]);
            mRNAsadj2 = mRNA_sadj; % mRNA_cnt./nuc_area;  % 
            ipars{2} = Rpars; 

        PlotmRNA = imresize(NucLabeled,.5,'nearest');
       PlotmRNA2 = PlotmRNA; 

       
       NucLabeled = imresize(NucLabeled,.5,'nearest'); 
          
    
%  

   
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
    
    sna = PlotmRNA2-sna_bkg;
    sna(sna<ThS) = 0; 
    sim = PlotmRNA-sim_bkg;
    sim(sna>ThS) = 0; 
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

  %%  
  
  Th = 30;
     sna = PlotmRNA2-sna_bkg;
      sim = PlotmRNA-sim_bkg;
      sna(sna<Th & sim<Th) = 0;
      sim(sna<Th & sim<Th) = 0;
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
    
    
    %%
    
    x = linspace(-2,2,20); y = linspace(0,1,20);
    [X,Y] = meshgrid(x,y);
    Z = exp(-X.^2) + .09*rand(20,20); Z(Z>1) = 1;  
    
    Zn = Z-.3 + .3*rand(20,20);  Zn(Zn<0) = 0;
    
  figure(1); clf; set(gcf,'color','w');
  subplot(2,2,1); surf(X,Y,Z); shading interp; view(22,66); colormap copper; %axis off; 
  Z2 = Z; Z2(Z>.5) = 1; Z2(Z<.5) = 0;
  subplot(2,2,3); surf(X,Y,Z2); shading interp; view(22,76); colormap copper; 
    
   subplot(2,2,2); surf(X,Y,Zn); shading interp; view(22,66); colormap copper;
  Zn2 = Zn; Zn2(Zn>.5) = 1; Zn2(Zn<.5) = 0;
  subplot(2,2,4); surf(X,Y,Zn2); shading interp; view(22,76); colormap copper; 
    
    %%
    
    