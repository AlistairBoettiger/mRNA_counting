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
Nuc{4} = imread([rawfolder,'2011-11/s08_MP06_cflip/','max_MP06_cflip_b_08.tif']);
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
mRNAs12telo{5} = data{3}.Data_sort(:,2);  mu12telo{5} = data{3}.mu;
mRNAs12telo{3} = data{4}.Data_sort(:,2);  mu12telo{3} = data{4}.mu;
mRNAs12telo{1} = data{10}.Data_sort(:,2); mu12telo{1} = data{10}.mu;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs12telo{4} = data{8}.Data_sort(:,2); mu12telo{4} = data{8}.mu;
%mRNAs{4} = data{8}.Data_sort(:,3); mu{4} = data{8}.mu(:,2); vout = '_o2'

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs12telo{2} = data{2}.Data_sort(:,2); mu12telo{2} = data{2}.mu;
%-----------------------------------------------------------------------%

%========================== during cc13====================================
% 'wt_sna'  emb = 09
% 'MP06_cflip' emb = 01, 02, 09, 10 
% 'MP06Hz_b' emb = 04, 05
slidedate = '2011-12/'; fname = 's142_sna';  ver = ''; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{1} = data{1}.Data_sort(:,2); mu13{1} = data{1}.mu;
mRNAs13{2} = data{2}.Data_sort(:,2); mu1{2} = data{2}.mu;

slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{1} = data{9}.Data_sort(:,2); mu13{1} = data{9}.mu;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{2} = data{1}.Data_sort(:,2); mu13{2} = data{1}.mu;
mRNAs13{3} = data{2}.Data_sort(:,2); mu13{3} = data{2}.mu;
mRNAs13{7} = data{9}.Data_sort(:,2); mu13{7} = data{9}.mu;
mRNAs13{6} = data{10}.Data_sort(:,2); mu13{6} = data{10}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13{4} = data{4}.Data_sort(:,2); mu13{4} = data{4}.mu;
mRNAs13{5} = data{5}.Data_sort(:,2); mu13{5} = data{5}.mu;
%========================================================================%



% --------------------- % metaphase cc13 ----------------------------------
% 'wt_sna'  emb = 05, 06 
% 'MP06_cflip' emb = 03, 04, 06
% 'MP06Hz_b' emb = 01, 09, 10

slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{3} = data{5}.Data_sort(:,2); mu13meta{3} = data{5}.mu;
mRNAs13meta{2} = data{6}.Data_sort(:,2); mu13meta{2} = data{6}.mu;

slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{1} = data{3}.Data_sort(:,2); mu13meta{1} = data{3}.mu;
mRNAs13meta{8} = data{4}.Data_sort(:,2); mu13meta{8} = data{4}.mu;
mRNAs13meta{5} = data{6}.Data_sort(:,2); mu13meta{5} = data{6}.mu;

slidedate = '2011-12/';   fname ='s140_sna'; ver = '_v2';% vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data');
mRNAs13meta{5} = data{7}.Data_sort(:,2); mu13meta{5} = data{7}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13meta{6} = data{1}.Data_sort(:,2); mu13meta{6} = data{1}.mu;
mRNAs13meta{7} = data{9}.Data_sort(:,2); mu13meta{7} = data{9}.mu;
mRNAs13meta{4} = data{10}.Data_sort(:,2); mu13meta{4} = data{10}.mu;
%------------------------------------------------------------------%

% =============== % telophase into cc14 ===============================
% 'MP06_cflip' emb = 05 
% 'MP06Hz_b' emb = 03, 07
% 'MP06Hz' emb = 10
slidedate = '2011-11/'; fname = 'MP06_cflip_b' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13telo{3} = data{5}.Data_sort(:,2); mu13telo{3} = data{5}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13telo{4} = data{3}.Data_sort(:,2); mu13telo{4} = data{3}.mu;
mRNAs13telo{1} = data{7}.Data_sort(:,2); mu13telo{1} = data{7}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs13telo{2} = data{10}.Data_sort(:,2);  mu13telo{2} = data{10}.mu;
% =====================================================================



% % ---------------- early cc14 ---------------------------%
% 'wt_sna'  emb = 01, 02,  07, 08,
% 'MP06Hz_b' emb = 06, 08, 11
% 'MP06Hz'; emb = 01, 02, 11, 12
slidedate = '2011-12/'; fname = 'wt_sna' ; ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14e{1} = data{1}.Data_sort(:,2);   mu14e{1} = data{1}.mu;
mRNAs14e{2} = data{2}.Data_sort(:,2);   mu14e{2} = data{2}.mu;
mRNAs14e{11} = data{7}.Data_sort(:,2);   mu14e{11} = data{7}.mu;
mRNAs14e{14} = data{8}.Data_sort(:,2);   mu14e{14} = data{8}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz_b' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14e{3} = data{6}.Data_sort(:,2);   mu14e{3} = data{6}.mu;
mRNAs14e{12} = data{8}.Data_sort(:,2);  mu14e{12} = data{8}.mu;
mRNAs14e{13} = data{11}.Data_sort(:,2);  mu14e{13} = data{11}.mu;

slidedate = '2011-05-22/'; fname = 's05_MP06Hz' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14e{4} = data{1}.Data_sort(:,2);  mu14e{4} = data{1}.mu;
mRNAs14e{15} = data{2}.Data_sort(:,2);  mu14e{15} = data{2}.mu;
mRNAs14e{16} = data{11}.Data_sort(:,2);  mu14e{16} = data{11}.mu;
mRNAs14e{17} = data{12}.Data_sort(:,2);  mu14e{17} = data{12}.mu;

slidedate = '2011-12/';   fname ='s140_sna'; ver = '_v2';% vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data');
mRNAs14e{7} = data{1}.Data_sort(:,2); mu14e{7} = data{1}.mu(:,1);
mRNAs14e{5} = data{3}.Data_sort(:,2); mu14e{5} = data{3}.mu(:,1);
mRNAs14e{9} = data{4}.Data_sort(:,2); mu14e{9} = data{4}.mu(:,1);
mRNAs14e{6} = data{5}.Data_sort(:,2); mu14e{6} = data{5}.mu(:,1);
mRNAs14e{8} = data{8}.Data_sort(:,2); mu14e{8} = data{8}.mu(:,1);
mRNAs14e{10} = data{9}.Data_sort(:,2); mu14e{10} = data{9}.mu(:,1);


%---------------------------------------------------------%

% ========================  steady state cc14 ~170 ==============================% 
% 'MP06Hz' emb = 03, 04, 05, 06, 07, 08

slidedate = '2011-05-22/'; fname = 's05_MP06Hz' ; ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14s{1} = data{7}.Data_sort(:,2); mu14s{1} = data{7}.mu;
mRNAs14s{2} = data{6}.Data_sort(:,2); mu14s{2} = data{6}.mu;

mRNAs14s{13} = data{5}.Data_sort(:,2); mu14s{13} = data{5}.mu;
mRNAs14s{14} = data{8}.Data_sort(:,2); mu14s{14} = data{8}.mu;

mRNAs14s{21} = data{3}.Data_sort(:,2); mu14s{21} = data{3}.mu;
mRNAs14s{22} = data{4}.Data_sort(:,2); mu14s{22} = data{4}.mu;


slidedate = '2011-02-17/'; fname = 'MP05_22C_sna_y';  ver = '_vN'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14s{3} = data{9}.Data_sort(:,2); mu14s{3} = data{9}.mu;
mRNAs14s{4} = data{12}.Data_sort(:,2); mu14s{4} = data{12}.mu;
mRNAs14s{5} = data{1}.Data_sort(:,2); mu14s{5} = data{1}.mu;
mRNAs14s{6} = data{2}.Data_sort(:,2); mu14s{6} = data{2}.mu;
mRNAs14s{7} = data{3}.Data_sort(:,2); mu14s{7} = data{3}.mu;
mRNAs14s{8} = data{5}.Data_sort(:,2); mu14s{8} = data{5}.mu;

mRNAs14s{15} = data{7}.Data_sort(:,2); mu14s{15} = data{7}.mu;
mRNAs14s{16} = data{8}.Data_sort(:,2); mu14s{16} = data{8}.mu;


slidedate = '2011-02-17/';   fname ='MP10_22C_sna_y_d'; ver = '_vN';% vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 

mRNAs14s{9} = data{4}.Data_sort(:,2); mu14s{9} = data{4}.mu(:,1);
mRNAs14s{10} = data{5}.Data_sort(:,2); mu14s{10} = data{5}.mu(:,1);
mRNAs14s{11} = data{2}.Data_sort(:,2); mu14s{11} = data{2}.mu(:,1);
mRNAs14s{12} = data{3}.Data_sort(:,2); mu14s{12} = data{3}.mu(:,1);

mRNAs14s{17} = data{7}.Data_sort(:,2); mu14s{17} = data{7}.mu(:,1);
mRNAs14s{18} = data{8}.Data_sort(:,2); mu14s{18} = data{8}.mu(:,1);
mRNAs14s{20} = data{1}.Data_sort(:,2); mu14s{20} = data{1}.mu(:,1);

slidedate = '2011-02-17/'; fname = 'MP05_22C_sna_y_c';  ver = '_v2'; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAs14s{23} = data{1}.Data_sort(:,2); mu14s{23} = data{1}.mu(:,1);
mRNAs14s{24} = data{2}.Data_sort(:,2); mu14s{24} = data{2}.mu(:,1);
mRNAs14s{19} = data{4}.Data_sort(:,2); mu14s{19} = data{4}.mu(:,1);

%% new data
mRNAsN{1} = NaN; muN{1} = NaN; 

slidedate = '2011-12/'; fname = 's142_sna';  ver = ''; vout = '';
load([folder,slidedate,fname,ver,'_slidedata',vout,'.mat'],'data'); 
mRNAsN{1} = data{1}.Data_sort(:,2); muN{1} = data{1}.mu;
mRNAsN{2} = data{2}.Data_sort(:,2); muN{2} = data{2}.mu;
mRNAsN{3} = data{3}.Data_sort(:,2); muN{3} = data{3}.mu;
mRNAsN{4} = data{4}.Data_sort(:,2); muN{4} = data{4}.mu;
mRNAsN{5} = data{5}.Data_sort(:,2); muN{5} = data{5}.mu;
mRNAsN{6} = data{6}.Data_sort(:,2); muN{6} = data{6}.mu;
mRNAsN{7} = data{7}.Data_sort(:,2); muN{7} = data{7}.mu;
mRNAsN{8} = data{8}.Data_sort(:,2); muN{8} = data{8}.mu;
mRNAsN{9} = data{9}.Data_sort(:,2); muN{9} = data{9}.mu;
mRNAsN{10} = data{10}.Data_sort(:,2); muN{10} = data{10}.mu;
mRNAsN{11} = data{11}.Data_sort(:,2); muN{11} = data{11}.mu;
mRNAsN{12} = data{12}.Data_sort(:,2); muN{12} = data{12}.mu;
mRNAsN{13} = data{13}.Data_sort(:,2); muN{13} = data{13}.mu;
mRNAsN{14} = data{14}.Data_sort(:,2); muN{14} = data{14}.mu;
mRNAsN{15} = data{15}.Data_sort(:,2); muN{15} = data{15}.mu;
mRNAsN{16} = data{16}.Data_sort(:,2); muN{16} = data{16}.mu;

mRNAs = [mRNAs12telo,mRNAs13,mRNAs13meta,mRNAs13telo,mRNAs14e,mRNAs14s,mRNAsN];
mu = [mu12telo,mu13,mu13meta,mu13telo,mu14e,mu14s,muN];

all_mRNA = {mRNAs12telo(2:end),mRNAs13,[mRNAs13meta(1:4),mRNAs13meta(6:7)],mRNAs13telo,mRNAs14e,mRNAs14s,mRNAsN};
all_mu = {mu12telo(2:end),mu13,[mu13meta(1:4),mu13meta(6:7)],mu13telo,mu14e,mu14s,muN};
G = length(all_mRNA);


%%
Es = length(mRNAs);  mRNA = cell(Es,1);bkd=mRNA; CMap = zeros(Es,3);

Ls = [length(mRNAs12telo),length(mRNAs13),length(mRNAs13meta),length(mRNAs13telo),length(mRNAs14e),length(mRNAs14s),length(mRNAsN)];
g = [Ls(1),sum(Ls(1:2)),sum(Ls(1:3)),sum(Ls(1:4)),sum(Ls(1:5)),sum(Ls(1:6)),sum(Ls(1:7))];

group = 1;  % embryo group
gi = 0; % index of embryos within group

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
    
    
    x = 2:Es;
    
    figure(2); hold on;
    plot(x+.35,meds_on(2:end),'ko'); hold on;
    plot(x+.35,smooth(meds_on(2:end)));
    
    figure(2);  hold on;
    bkd_curve = smooth(meds_off(2:end));
    plot(x+.35,150*bkd_curve); 
    plot(x+.35,150*meds_off(2:end),'o');

    
    
  
%% new section

mRNAon = cell(6,1);

med_cnt = zeros(G,1);
upp_cnt = zeros(G,1);
low_cnt = zeros(G,1);
std_cnt = zeros(G,1); 
for k=1:G
    for e=1:length(all_mRNA{k})
            sna = all_mRNA{k}{e} - min(all_mu{k}{e});
            pk = max(all_mu{k}{e}-min(all_mu{k}{e}));
            mRNAon{k}{e} =  sna( sna>1/2*pk); 
    end  
     med_cnt(k) = median(cellfun(@mean,mRNAon{k}));
     
     sort_mean_cnt = sort(cellfun(@mean,mRNAon{k}));
     upp_cnt(k) = sort_mean_cnt(round(3/4*length(sort_mean_cnt)));
     low_cnt(k) = sort_mean_cnt(round(1/4*length(sort_mean_cnt)));
     std_cnt(k) = std(sort_mean_cnt);
end    

%%  half-life

% .5 = exp(-7.5/tau)
% tau = -7.5/log(.5) = 10.82
% exp(-2/10.82) 

prev_tau = 10.82;
T = 2; % estimate of mitosis length

my_tau = -T/log( 2*med_cnt(4)/med_cnt(3) );
max_tau =  max([-T/log( 2*upp_cnt(4)/med_cnt(3) ),-T/log( 2*med_cnt(4)/low_cnt(3))]);
min_tau = min([-T/log( 2*low_cnt(4)/med_cnt(3) ),-T/log( 2*med_cnt(4)/upp_cnt(3))]);
% tau = T/log(2m1/m2)
% 
% Dtau = sqrt( (dtau/dT*DT)^2 + dtau/dm1*Dm1)^2 + (dtau/dm2*Dm2)^2 )

m1 = med_cnt(4); m2 = med_cnt(3); 
 
 DT = .25; % 15 seconds uncertainty 
Dtau = sqrt( (-1/log( 2*m1/m2)*DT)^2+(T/(log( 2*m1/m2)^2*m1)*std_cnt(4))^2 + (T/(log( 2*m1/m2)^2*m2)*std_cnt(3))^2)
% however error in tau is not semetric, so  

t = linspace(0,10,25);
degred = med_cnt(3)*exp(-t/my_tau); 
low_d = med_cnt(3)*exp(-t/max_tau);
high_d = med_cnt(3)*exp(-t/min_tau);
figure(4); clf; plot(t,degred,'k-','linewidth',3); hold on;
plot(t,low_d,'k--');

plot(0,med_cnt(3),'.','Color',[.5,0,.5],'MarkerSize',30);
plot(2,2*med_cnt(4),'.','Color',[.7,0,.3],'MarkerSize',30);
plot(7.5,.5*med_cnt(3),'.','Color',[.6,.5,.1],'MarkerSize',30);
plot(-10,-10, 'k+','MarkerSize',10);
plot(0*ones(1,length(mRNAon{3})),cellfun(@mean,mRNAon{3}),'+','Color',[.5,0,.5],'MarkerSize',10);
plot([2*ones(1,length(mRNAon{4}))],[2*cellfun(@mean,mRNAon{4})],'+','MarkerSize',10,'Color',[.7,0,.3]);
 plot(t,high_d,'k--');
 
xlim([min(t),max(t)]);
legend(['measured halflife =',num2str(my_tau*log(2),3),' min ',...
    '(',num2str(min_tau*log(2),2),',',num2str(max_tau*log(2),2),')'],'error bounds',...
    'prometaphase','late telophase','measured ftz halflife = 7.5 min',...
    'mean counts from indiv. embryos');
set(gcf,'color','w');
ylabel('mRNA count'); xlabel('time (min)');
title('mRNA degredation during cell division');
ylim([0,300]);

%%  mRNA synthesis
t = linspace(0,20,25);
lifetime = my_tau/log(2);

Tdiv = 16; % time of division / prophase
m1 = med_cnt(1); m2 = med_cnt(3);  
synthesis_rate = (m2 - m1*exp(-Tdiv/lifetime))/(lifetime*(1-exp(-Tdiv/lifetime)))  
min_syn = (low_cnt(3) - upp_cnt(1)*exp(-Tdiv/lifetime))/(lifetime*(1-exp(-Tdiv/lifetime)))  
max_syn = (upp_cnt(3) - low_cnt(1)*exp(-Tdiv/lifetime))/(lifetime*(1-exp(-Tdiv/lifetime)))  
 
% % with my tau (min tau really blows this up)
% min_syn = (low_cnt(3) - upp_cnt(1)*exp(-Tdiv/max_tau))/(max_tau*(1-exp(-Tdiv/max_tau)))  
% max_syn = (upp_cnt(3) - low_cnt(1)*exp(-Tdiv/min_tau))/(min_tau*(1-exp(-Tdiv/min_tau)))  
 

prod = med_cnt(1)*exp(-t/lifetime) + synthesis_rate*lifetime*(1-exp(-t/lifetime));
upp_prod = med_cnt(1)*exp(-t/lifetime) + min_syn*lifetime*(1-exp(-t/lifetime));
low_prod = med_cnt(1)*exp(-t/lifetime) + max_syn*lifetime*(1-exp(-t/lifetime));

figure(5); clf; plot(t,prod,'k--','linewidth',3); hold on;
plot(t,upp_prod,'k--');
plot(0,med_cnt(1),'b.','MarkerSize',30);
plot(6.5,med_cnt(2),'.','Color',[2/8,0,1-2/8],'MarkerSize',30);
plot(16,med_cnt(3),'.','Color',[4/8,0,1-4/8],'MarkerSize',30);
plot(-10,-10, 'k+','MarkerSize',10);

plot([0*ones(1,length(mRNAon{1}))],[cellfun(@mean,mRNAon{1})],'b+','MarkerSize',10);
plot([16*ones(1,length(mRNAon{3}))],[cellfun(@mean,mRNAon{3})],'+','Color',[4/8,0,1-4/8],'MarkerSize',10);
plot(2*(1:length(mRNAon{2})),[cellfun(@mean,mRNAon{2})],'+','Color',[2/8,0,1-2/8],'MarkerSize',10);
 plot(t,low_prod,'k--');
set(gcf,'color','w');
ylabel('mRNA count'); xlabel('time (min)');

xlim([0,20]); ylim([0,270]);

legend(['tx rate =', num2str(synthesis_rate/2,2),'mRNA/min (', num2str(min_syn/2,2),',',num2str(max_syn/2,2), ')'],...
    'uncertainty bound','ave # mRNA at telophase of cc12',...
    'ave # mRNA of interphase embryos','ave # mRNA at prometaphase cc13',...
     'mean counts from indiv. embryos','Location','South')
 title('mRNA synthesis rate');
%%  mRNA synthesis
t = linspace(0,20,25);
synthesis_rate = 21; lifetime = my_tau/log(2);
prod = med_cnt(1)*exp(-t/lifetime) + synthesis_rate*lifetime*(1-exp(-t/lifetime));
figure(5); clf; plot(t,prod); hold on;

   %  plot(.7*rand(1,L)+e*ones(1,L),mRNA{e},'.','color',CMap(e,:),'MarkerSize',5);
   
cc12s = cell2mat(mRNAon{1}'); 
plot(.7*rand(1,length(cc12s)), cc12s,'.'); hold on;
cc13t = cell2mat(mRNAon{3}'); 
plot(16+.7*rand(1,length(cc13t)), cc13t,'.','color','r');
for e=1:length(mRNAon{2})
    L = length(mRNAon{2}{e});
    plot(.7*rand(1,L)+2*e*ones(1,L),mRNAon{2}{e},'.','color','m','MarkerSize',5);
end


 plot(t,prod,'k-','linewidth',3);

set(gcf,'color','w');
ylabel('mRNA count'); xlabel('time (min)');

legend(['tx rate =', num2str(synthesis_rate/2,2),'mRNA/min'],...
    '# mRNA at telophase of cc12','# mRNA at prometaphase cc13',...
    '# mRNA of interphase embryos')

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





