%% mRNA dynamics, show different nuclei patterns. 
clear all;

% load all data sorted by age

  folder = 'C:\Users/Alistair/My Documents/Projects/mRNA_counting/Data/';
   rawfolder = 'C:\Users/Alistair/Data/';
   
 %% Load Nuc data
 
 % ----------------------------% exiting cc12 ---------------------------%
% 'wt_sna'  emb = 03, 04, 10
% 'MP06_cflip' emb = 08
% 'MP06Hz_b' emb = 02
Nuc{5} = imread([rawfolder,'2011-12/','max_wt_sna_03.tif']);
Nuc{3} = imread([rawfolder,'2011-12/','max_wt_sna_04.tif']);
Nuc{1} = imread([rawfolder,'2011-12/','max_wt_sna_10.tif']);
Nuc12{4} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_08.tif']);
Nuc{2}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_02.tif']);
%-----------------------------------------------------------------------%

%========================== during cc13====================================
% 'wt_sna'  emb = 09
% 'MP06_cflip' emb = 01, 02, 09, 10 
% 'MP06Hz_b' emb = 04, 05

Nuc{6} = imread([rawfolder,'2011-12/','max_wt_sna_09.tif']);
Nuc{7} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_01.tif']);
Nuc{8} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_02.tif']);
Nuc{12} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_09.tif']);
Nuc{11} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_10.tif']);
Nuc{9}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_04.tif']);
Nuc{10}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_05.tif']);
%========================================================================%


% --------------------- % exiting cc13 -------------------------------------
% 'wt_sna'  emb = 05, 06 
% 'MP06_cflip' emb = 03, 04, 06
% 'MP06Hz_b' emb = 01, 09, 10
Nuc{15} = imread([rawfolder,'2011-12/','max_wt_sna_05.tif']);  % meta
Nuc{14} = imread([rawfolder,'2011-12/','max_wt_sna_06.tif']); % meta
Nuc{13} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_03.tif']); % meta
Nuc{20} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_04.tif']); % ana
Nuc{17} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_06.tif']); %meta
Nuc{18}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_01.tif']); % late meta
Nuc{19}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_09.tif']); % late meta
Nuc{16}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_10.tif']); % ana
%------------------------------------------------------------------%

% =============== % telophase into cc14 ===============================
% 'MP06_cflip' emb = 05 
% 'MP06Hz_b' emb = 03, 07
% 'MP06Hz' emb = 10
Nuc{21}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_07.tif']);
Nuc{22}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_10.tif']);
Nuc{23} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_05.tif']);
Nuc{24}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_03.tif']);
% =====================================================================



% % ---------------- early cc14 ---------------------------%
% 'wt_sna'  emb = 01, 02,  07, 08,
% 'MP06Hz_b' emb = 06, 08, 11
% 'MP06Hz'; emb = 01, 02, 11, 12
Nuc{25} = imread([rawfolder,'2011-12/','max_wt_sna_01.tif']);
Nuc{26} = imread([rawfolder,'2011-12/','max_wt_sna_02.tif']);
Nuc{27} = imread([rawfolder,'2011-12/','max_wt_sna_07.tif']);
Nuc{32} = imread([rawfolder,'2011-12/','max_wt_sna_08.tif']);
Nuc{29}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_06.tif']);
Nuc{30}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_08.tif']);
Nuc{31}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_b_11.tif']);
Nuc{28}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_01.tif']);
Nuc{33}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_02.tif']);
Nuc{34}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_11.tif']);
Nuc{35}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_12.tif']);
%---------------------------------------------------------%


% ========================  steady state cc14 ~170 ==============================% 
% 'MP06Hz' emb = 03, 04, 05, 06, 07, 08
Nuc{36}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_07.tif']);
Nuc{37}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_06.tif']);
Nuc{56}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_03.tif']);
Nuc{57}= imread([rawfolder,'2011-05-22/s05_MP06/','max_s05_MP06Hz_04.tif']);

% =================================================================== % 









   
 %% Load mRNA data  
 
 clear mRNAs mu; 
 
% ----------------------------% exiting cc12 ---------------------------%
% 'wt_sna'  emb = 03, 04, 10
% 'MP06_cflip' emb = 08
% 'MP06Hz_b' emb = 02
slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{5} = data{3}.Data_sort(:,2);  mu{5} = data{3}.mu;
mRNAs{3} = data{4}.Data_sort(:,2);  mu{3} = data{4}.mu;
mRNAs{1} = data{10}.Data_sort(:,2); mu{1} = data{10}.mu;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{4} = data{8}.Data_sort(:,2); mu{4} = data{8}.mu;
%mRNAs{4} = data{8}.Data_sort(:,3); mu{4} = data{8}.mu(:,2); vout = '_o2'

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{2} = data{2}.Data_sort(:,2); mu{2} = data{2}.mu;
%-----------------------------------------------------------------------%

%========================== during cc13====================================
% 'wt_sna'  emb = 09
% 'MP06_cflip' emb = 01, 02, 09, 10 
% 'MP06Hz_b' emb = 04, 05
slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{6} = data{9}.Data_sort(:,2); mu{6} = data{9}.mu;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{7} = data{1}.Data_sort(:,2); mu{7} = data{1}.mu;
mRNAs{8} = data{2}.Data_sort(:,2); mu{8} = data{2}.mu;
mRNAs{12} = data{9}.Data_sort(:,2); mu{12} = data{9}.mu;
mRNAs{11} = data{10}.Data_sort(:,2); mu{11} = data{10}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{9} = data{4}.Data_sort(:,2); mu{9} = data{4}.mu;
mRNAs{10} = data{5}.Data_sort(:,2); mu{10} = data{5}.mu;
%========================================================================%



% --------------------- % exiting cc13 -------------------------------------
% 'wt_sna'  emb = 05, 06 
% 'MP06_cflip' emb = 03, 04, 06
% 'MP06Hz_b' emb = 01, 09, 10

slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{15} = data{5}.Data_sort(:,2); mu{15} = data{5}.mu;
mRNAs{14} = data{6}.Data_sort(:,2); mu{14} = data{6}.mu;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{13} = data{3}.Data_sort(:,2); mu{13} = data{3}.mu;
mRNAs{20} = data{4}.Data_sort(:,2); mu{20} = data{4}.mu;
mRNAs{17} = data{6}.Data_sort(:,2); mu{17} = data{6}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{18} = data{1}.Data_sort(:,2); mu{18} = data{1}.mu;
mRNAs{19} = data{9}.Data_sort(:,2); mu{19} = data{9}.mu;
mRNAs{16} = data{10}.Data_sort(:,2); mu{16} = data{10}.mu;
%------------------------------------------------------------------%

% =============== % telophase into cc14 ===============================
% 'MP06_cflip' emb = 05 
% 'MP06Hz_b' emb = 03, 07
% 'MP06Hz' emb = 10
slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{23} = data{5}.Data_sort(:,2); mu{21} = data{5}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{24} = data{3}.Data_sort(:,2); mu{22} = data{3}.mu;
mRNAs{21} = data{7}.Data_sort(:,2); mu{23} = data{7}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{22} = data{10}.Data_sort(:,2);  mu{24} = data{10}.mu;
% =====================================================================



% % ---------------- early cc14 ---------------------------%
% 'wt_sna'  emb = 01, 02,  07, 08,
% 'MP06Hz_b' emb = 06, 08, 11
% 'MP06Hz'; emb = 01, 02, 11, 12
slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{25} = data{1}.Data_sort(:,2);   mu{25} = data{1}.mu;
mRNAs{26} = data{2}.Data_sort(:,2);   mu{26} = data{2}.mu;
mRNAs{27} = data{7}.Data_sort(:,2);   mu{27} = data{7}.mu;
mRNAs{32} = data{8}.Data_sort(:,2);   mu{32} = data{8}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{29} = data{6}.Data_sort(:,2);   mu{29} = data{6}.mu;
mRNAs{30} = data{8}.Data_sort(:,2);  mu{30} = data{8}.mu;
mRNAs{31} = data{11}.Data_sort(:,2);  mu{31} = data{11}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{28} = data{1}.Data_sort(:,2);  mu{28} = data{1}.mu;
mRNAs{33} = data{2}.Data_sort(:,2);  mu{33} = data{2}.mu;
mRNAs{34} = data{11}.Data_sort(:,2);  mu{34} = data{11}.mu;
mRNAs{35} = data{12}.Data_sort(:,2);  mu{35} = data{12}.mu;
%---------------------------------------------------------%

% ========================  steady state cc14 ~170 ==============================% 
% 'MP06Hz' emb = 03, 04, 05, 06, 07, 08

slidedate = '2011-05-22/'; fname = 's05_MP06Hz' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{36} = data{7}.Data_sort(:,2); mu{36} = data{7}.mu;
mRNAs{37} = data{6}.Data_sort(:,2); mu{37} = data{6}.mu;

mRNAs{48} = data{5}.Data_sort(:,2); mu{48} = data{5}.mu;
mRNAs{49} = data{8}.Data_sort(:,2); mu{49} = data{8}.mu;

mRNAs{56} = data{3}.Data_sort(:,2); mu{56} = data{3}.mu;
mRNAs{57} = data{4}.Data_sort(:,2); mu{57} = data{4}.mu;


slidedate = '2011-02-17/'; fname = 'MP05_22C_sna_y';  ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{38} = data{9}.Data_sort(:,2); mu{38} = data{9}.mu;
mRNAs{39} = data{12}.Data_sort(:,2); mu{39} = data{12}.mu;
mRNAs{40} = data{1}.Data_sort(:,2); mu{40} = data{1}.mu;
mRNAs{41} = data{2}.Data_sort(:,2); mu{41} = data{2}.mu;
mRNAs{42} = data{3}.Data_sort(:,2); mu{42} = data{3}.mu;
mRNAs{43} = data{5}.Data_sort(:,2); mu{43} = data{5}.mu;

mRNAs{50} = data{7}.Data_sort(:,2); mu{50} = data{7}.mu;
mRNAs{51} = data{8}.Data_sort(:,2); mu{51} = data{8}.mu;

% slidedate = '2011-02-17/';   fname ='MP10_22C_sna_y_d'; ver = '_v3';% vout = '';
% load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
% mRNAs{48} = data{7}.Data_sort(:,2); mu{48} = data{7}.mu(:,1);
% mRNAs{49} = data{8}.Data_sort(:,2); mu{49} = data{8}.mu(:,1);

slidedate = '2011-02-17/'; fname = 'MP05_22C_sna_y_c';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs{58} = data{1}.Data_sort(:,2); mu{58} = data{1}.mu(:,1);
mRNAs{59} = data{2}.Data_sort(:,2); mu{59} = data{2}.mu(:,1);
mRNAs{54} = data{4}.Data_sort(:,2); mu{54} = data{4}.mu(:,1);

slidedate = '2011-02-17/';   fname ='MP10_22C_sna_y_d'; ver = '_v3';% vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 

mRNAs{44} = data{4}.Data_sort(:,2); mu{44} = data{4}.mu(:,1);
mRNAs{45} = data{5}.Data_sort(:,2); mu{45} = data{5}.mu(:,1);
mRNAs{46} = data{2}.Data_sort(:,2); mu{46} = data{2}.mu(:,1);
mRNAs{47} = data{3}.Data_sort(:,2); mu{47} = data{3}.mu(:,1);

mRNAs{52} = data{7}.Data_sort(:,2); mu{52} = data{7}.mu(:,1);
mRNAs{53} = data{8}.Data_sort(:,2); mu{53} = data{8}.mu(:,1);
mRNAs{55} = data{1}.Data_sort(:,2); mu{55} = data{1}.mu(:,1);

%%
Es = length(mRNAs);  mRNA = cell(Es,1);bkd=mRNA; CMap = zeros(Es,3);
g = [5,12,20,24,35,Es+1];% group boundaries
group = 1;  % embryo group
gi = 0; % index of embryos within group

cc12telo = cell(g(1),1);
cc13 = cell(g(2)-g(1),1);
cc13meta = cell(g(3)-g(2),1);
cc13telo = cell(g(4)-g(3),1);
cc14early = cell(g(5)-g(4),1);
cc14late = cell(g(6)-g(5),1);
meds_on = zeros(Es,1);
meds_off = zeros(Es,1); 

dyn_fig = figure(2); clf; colordef white; set(gcf,'color','w'); hold on;
    for e=2:Es
        
    if e > g(group)
        group = group + 1;
        gi = 0; 
    end
    CMap(e,:) = [group/length(g),0,1-group/length(g)] ;
    
        
        sna = mRNAs{e} - min(mu{e});
        pk = max(mu{e}-min(mu{e}));
        mRNA{e} =  sna( sna>1/2*pk); 
        bkd{e} = sna( sna<1/2*pk);
        figure(2);  
        L = length(mRNA{e});
        plot(.7*rand(1,L)+e*ones(1,L),mRNA{e},'.','color',CMap(e,:),'MarkerSize',5);
        meds_on(e) = median(mRNA{e});
        
     
        L2 = length(bkd{e});
        plot(.7*rand(1,L2)+e*ones(1,L2),bkd{e},'.','color',CMap(e,:)*.2+[.8,.8,.8],'MarkerSize',5);
        % meds_off(e) = mean(bkd{e}(bkd{e}>0));
         meds_off(e) = sum(bkd{e}>15)/(L2+L);
         
      
    
    end 
    
    ylim([0,300]); set(gca,'FontSize',14);
ylabel('mRNA counts per cell');
    mean(cc13telo)
    
    x = 2:Es;
    
    figure(2); hold on;
    plot(x+.35,meds_on(2:end),'ko'); hold on;
    plot(x+.35,smooth(meds_on(2:end)));
    
    figure(2);  hold on;
    bkd_curve = smooth(meds_off(2:end));
    plot(x+.35,150*bkd_curve); 
    plot(x+.35,150*meds_off(2:end),'o');
    

%%

fout = 'C:\Users/Alistair/My Documents/Projects/Snail Patterning/Results/';
saveas(dyn_fig,[fout,'snail_dyn.eps'],'eps');

%%

fout = 'C:\Users/Alistair/My Documents/Projects/Snail Patterning/Results/';
Nfig = figure(3); clf; 
e = 2; Nchn = Nuc{e}(:,:,3); 
imagesc(Nchn); colormap copper; 
caxis([mean(Nchn(:))/1.1,3*mean(Nchn(:))]); 
set(gcf,'color','k');
saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);

Nfig = figure(3); clf; 
e = 10; Nchn = Nuc{e}(:,:,3); 
imagesc(Nchn); colormap copper; 
caxis([mean(Nchn(:))/2,2*mean(Nchn(:))]); 
set(gcf,'color','k');
saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);

Nfig = figure(3); clf; 
e = 13; Nchn = Nuc{e}(:,:,3); 
imagesc(Nchn); colormap copper; 
caxis([mean(Nchn(:))/2,3*mean(Nchn(:))]); 
set(gcf,'color','k');
saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);

Nfig = figure(3); clf;
e = 15; Nchn = Nuc{e}(:,:,3); 
imagesc(Nchn); colormap copper; 
caxis([mean(Nchn(:))/2,4*mean(Nchn(:))]); 
set(gcf,'color','k');
saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);

Nfig = figure(3); clf;
e = 21; Nchn = Nuc{e}(:,:,3); 
imagesc(Nchn); colormap copper; 
caxis([mean(Nchn(:))/4,4*mean(Nchn(:))]);
set(gcf,'color','k');
saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);

Nfig = figure(3); clf;
e = 23; Nchn = Nuc{e}(:,:,3); 
imagesc(Nchn); colormap copper; 
caxis([mean(Nchn(:))/4,4*mean(Nchn(:))]); 
set(gcf,'color','k');
saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);

Nfig = figure(3); clf;
e = 28; Nchn = Nuc{e}(:,:,3);  % 28
imagesc(Nchn); colormap copper; 
caxis([mean(Nchn(:))/4,4*mean(Nchn(:))]); 
set(gcf,'color','k');
saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);


Nfig = figure(3); clf;
e = 56; Nchn = Nuc{e}(:,:,3);  % 28
imagesc(Nchn); colormap copper; 
caxis([mean(Nchn(:))/1.5,2*mean(Nchn(:))]); 
set(gcf,'color','k');
saveas(Nfig,[fout,'Nucs_emb',num2str(e),'.png']);

%%  half-life

% .5 = exp(-7.5/tau)
% tau = -7.5/log(.5) = 10.82
exp(-2/10.82) 



