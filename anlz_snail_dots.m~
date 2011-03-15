
%% anlz_2chn_dots.m

% Alistair Boettiger                                Date Begun: 02/02/11
% Levine Lab                                    Last Modified: 02/11/11

% Analyze snail and yellow data

%% Load data, histogram cell_counts
clear all;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
%fname = 'MP05_sna_y_22C_04';
% fname = 'MP10_22C_sna_y_06';
% fname = 'MP05_sna_y_22C_06_bleachdata';
% fname = 'MP05_sna_y_22C_06_bleach_lowT';


 %fname = 'MP05_sna_y_22C_04';
 fname = 'MP05_22C_sna_y_01_p1p05';
 
 
load([folder,fname,'.mat']);  

    % Notes:
    % MP05_22C_sna_y from 02-17-11 (slide is 02-06-11?), threshs: .1, .05
    %        miss some possible dots at .1.  pretty clean at .05.  
    % MP05_22C_sna_y_01_p1p05
    % MP05_22C_sna_y_02_p1p05
    %
    % at .07, .03
    % MP05_22C_sna_y_01
    % MP05_22C_sna_y_02
    % 
    
    
    % 
    %   with thresholds: .07 .03
    % 'MP05_sna_y_22C_01' is sna mean 600
    % 'MP05_sna_y_22C_02' is young 300
     % 'MP05_sna_y_22C_03' is sna mean 600
      % 'MP05_sna_y_22C_04' is sna mean 600
      % _06 bleached
      % 07 sna mean 600 800
      % 08 is m = 700 but off peak is 150
      % 08_snap7 is 600
      % 09_snap15 is 500.  Not same number of cells!
    
 % Plotting for trouble shooting
%     figure(1); clf; 
%     subplot(2,1,1); hist(mRNA_sadj1,30);  title('snail');
%     subplot(2,1,2);  hist(mRNA_sadj2,30); title('y');
   


    MP05_sna = mRNA_sadj1;
    MP05_y = mRNA_sadj2;
    MP05_Nucs = NucLabeled;  Nnucs05 = max(NucLabeled(:)); 
    MP05_conn =  connectivity(NucLabeled);
    
    mRNA = linspace(0,1000,50); 
    
    figure(2); clf; set(gcf,'color','k'); colordef black;
    subplot(2,2,1); hist(MP05_sna,mRNA);  title('snail'); xlim([0,max(mRNA)]);
    subplot(2,2,2);  hist(MP05_y,mRNA); title('MP05 y'); xlim([0,max(mRNA)]);
    
    
    % fname = 'MP10_sna_y_22C_01';  
    fname = 'MP10_22C_sna_y_c_06'; 

    % Notes:
    % Data from 02-08-11   (MOSTLY STARTS TOO DEEP, EXCEPT IM 01)
    % 'MP10_22C_sna_y_01' starts correctly.  
    % 'MP10_22C_sna_y_02'  misses some apical endogenous transcripts.  probably cc13
    % 'MP10_22C_sna_y_04' very young
    % 'MP10_22C_sna_y_05' misses a few apical endogenous transcripts
    % 'MP10_22C_sna_y_06' starts mid nuclei misses MANY apical endogenous transcripts
    % 'MP10_22C_sna_y_07' starts w/ nucs, misses many apical endogenous transcripts
    % 'MP10_22C_sna_y_08' starts w/ nucs, misses many apical endogenous transcripts  
    % 'MP10_22C_sna_y_09' starts w/ nucs, misses many apical endogenous transcripts  
    % 'MP10_22C_sna_y_10' starts w/ nucs, misses many apical endogenous transcripts  
    
    % Data from 01-31-11   (MOSTLY STARTS TOO DEEP, EXCEPT IM 01)
    % 'MP10_sna_y_22C_01'; % good image, 
    % 'MP10_sna_y_22C_01_0503'; poor threshold choice
    %  'MP10_sna_y_22C_02'; -- starts Much Much two deep, still shows nice  correlation of y to sna, though misses half of transcripts.
    % 03 - 06 all WAY too deep, miss the top half of every cell
    % 07 -- young embryo.  missing some?
    % 08,09 too deep.  10 + no data. 
    % 
    % Series b 01-31-11   (MOSTLY STARTS TOO DEEP, EXCEPT IM 01) 
    % 'MP10_sna_y_22C_b_01'; % this one is pretty good
    % 'MP10_sna_y_22C_b_02' % starts just a bit low, misses some apical transcripts   
    % 'MP10_sna_y_22_b_03' starts w/ nucs, misses many apical endogenous transcripts     
    %
    % Data from 02-06-11 folder (taken 02-12-11).  (GOOD DATA)
    % MP10_22C_sna_y_c_01 mean sna =400, bit young
    % _02 very younger mean sna =200,
    %_03  younger mean sna =300,
    %_04  younger mean sna =300,
    % _05 mean sna = 600, but y less than half... :(
    % _06 mean sna =700. 
    % _07 mean sna =500. 
    % 08
    %09   incomplete i,age
    
    
load([folder,fname,'.mat']);
    
    MP10_sna = mRNA_sadj1;
    MP10_y = mRNA_sadj2;
    MP10_Nucs = NucLabeled;  Nnucs10 = max(NucLabeled(:)); 
     MP10_conn =  connectivity(NucLabeled);
 
    subplot(2,2,3); hist(MP10_sna,mRNA);  title('snail'); xlim([0,max(mRNA)]);
    subplot(2,2,4);  hist(MP10_y,mRNA); title('MP10 y'); xlim([0,max(mRNA)]);
 
    
    %%
    
[MP05_sna_plot,MP05_y_plot,MP05_sna_var,MP05_y_var] = fxn_compdotvar(MP05_Nucs,MP05_conn,MP05_sna,MP05_y,Nnucs05);
[MP10_sna_plot,MP10_y_plot,MP10_sna_var,MP10_y_var] = fxn_compdotvar(MP10_Nucs,MP10_conn,MP10_sna,MP10_y,Nnucs10);
            
            
MP05_sna_on = MP05_sna > mean(MP05_sna)*.6;   
MP05_y_on = MP05_y > mean(MP05_y)*.6; 
            
MP05_mean_sna_var = mean(MP05_sna_var(MP05_sna_on));
MP05_mean_y_var = mean(MP05_y_var(MP05_y_on));

  
MP10_sna_on = MP10_sna > mean(MP10_sna)*.6;   
MP10_y_on = MP10_y > mean(MP10_y)*.6; 
            
MP10_mean_sna_var = mean(MP10_sna_var(MP10_sna_on));
MP10_mean_y_var = mean(MP10_y_var(MP10_y_on));
 


% figure(4); clf; 
% subplot(2,2,1); scatter(MP05_sna,MP05_sna_var);
% subplot(2,2,2); scatter(MP05_y,MP05_y_var);
% subplot(2,2,3); scatter(MP10_sna,MP10_sna_var);
% subplot(2,2,4); scatter(MP10_y,MP10_y_var);

 %          
       
  
figure(3); clf; cmax = 800;
 colordef black;
 subplot(2,2,1); imagesc(MP05_sna_plot); colormap('hot'); colorbar; 
 set(gcf,'color','k');% caxis([0,cmax]); 
 title(['MP05 sna, local var = ',num2str(MP05_mean_sna_var,2)]  );
 
 subplot(2,2,2); imagesc(MP05_y_plot); colormap('hot'); colorbar;
 set(gcf,'color','k'); 
 title(['MP05 y, local var = ',num2str(MP05_mean_y_var,2)]  );
 
  subplot(2,2,3); imagesc(MP10_sna_plot); colormap('hot'); colorbar; 
 set(gcf,'color','k');% caxis([0,cmax]); 
 title(['MP10 sna, local var = ',num2str(MP10_mean_sna_var,2)]  );
 
 subplot(2,2,4); imagesc(MP10_y_plot); colormap('hot'); colorbar;
 set(gcf,'color','k'); 
 title(['MP10 y, local var = ',num2str(MP10_mean_y_var,2)]  );
  
    
 %% Define expression region
spread = 1.35;
t1 = .25;

figure(5); clf; subplot(2,2,1);
[MP05_sna_on_cnts,MP05_sna_off_cnts]= fxn_regionvar(MP05_Nucs,MP05_sna_plot,MP05_sna,t1,spread,Nnucs05);
subplot(2,2,2);
[MP05_y_on_cnts,MP05_y_off_cnts]= fxn_regionvar(MP05_Nucs,MP05_y_plot,MP05_y,t1,spread,Nnucs05);
subplot(2,2,3);
[MP10_sna_on_cnts,MP10_sna_off_cnts]= fxn_regionvar(MP10_Nucs,MP10_sna_plot,MP10_sna,t1,spread,Nnucs10);
subplot(2,2,4);
[MP10_y_on_cnts,MP10_y_off_cnts]= fxn_regionvar(MP10_Nucs,MP10_y_plot,MP10_y,t1,spread,Nnucs10);































 %%
 
 figure(4); clf;
 
 DepthDots(In,Cell_bnd,inds_Z,h,w);

  
%% Define expression region

spread = 1.5;
t1 = .25;

% Automatic Threshold
C1 = uint8( cell_sadj1/max(cell_sadj1(:))*255);
%t1 = graythresh(C1);
bw1 = im2bw(C1,t1); % 
bw1 = imclose(bw1,strel('disk',100)); 
bw1 = imfill(bw1,'holes');
bw1 = bwareaopen(bw1,2E4);

bndry1 = bwboundaries(bw1)
figure(3); clf; subplot(2,2,1); imshow(C1);  hold on;
plot(bndry1{1}(:,2),bndry1{1}(:,1),'c');
subplot(2,2,2); imshow(bw1); hold on;



onMask = ismember(NucLabeled,inReg);

outReg = setdiff(1:Nnucs,inReg);
on_cnts = mRNA_sadj1(inReg);
off_cnts = mRNA_sadj1(outReg);


on_mean = mean(on_cnts);
on_std = std(on_cnts);
off_mean = mean(off_cnts);


m = linspace(0,cmax,30);
figure(2); clf; 
subplot(2,2,1); hist(on_cnts,m); 
title(['sna mean=',num2str(on_mean,4),' std=',num2str(on_std,4),' cov=',num2str(on_std/on_mean,3)]);
subplot(2,2,2); hist(off_cnts,m); title(['mean=',num2str(off_mean,4)]);


over1 = find(mRNA_sadj1>on_mean*spread);
under1 = find(mRNA_sadj1<on_mean/spread); 

Over1 = ismember(NucLabeled,over1).*onMask; 
Under1 = ismember(NucLabeled,under1).*onMask; 

I = uint8(zeros(h,w,3));
I(:,:,1) =  C1;
I(:,:,2) = C1 - uint8(255*Over1);
I(:,:,3) = uint8(255*Under1)+C1 -uint8(255*Over1);
figure(4); clf; subplot(1,2,1); imshow(I);
title(['sna  mean=',num2str(on_mean,4),' std=',num2str(on_std,4),' cov=',num2str(on_std/on_mean,3)]);
set(gcf,'color','k');
 hold on;
plot(bndry1{1}(:,2),bndry1{1}(:,1),'c');
%

% 
% 
% Automatic Threshold
C1 = uint8( cell_sadj2/max(cell_sadj2(:))*255);
%t1 = graythresh(C1);
bw1 = im2bw(C1,t1); % 
bw1 = imclose(bw1,strel('disk',100)); 
bw1 = imfill(bw1,'holes');
bw1 = bwareaopen(bw1,2E4);

bndry1 = bwboundaries(bw1);
figure(3);  subplot(2,2,3); imshow(C1); hold on;
plot(bndry1{1}(:,2),bndry1{1}(:,1),'c');
subplot(2,2,4); imshow(bw1); hold on;

inReg = unique(bw1.*NucLabeled);
inReg(inReg==0) = [];
outReg = setdiff(1:Nnucs,inReg);

onMask = ismember(NucLabeled,inReg);

on_cnts = mRNA_sadj2(inReg);
off_cnts = mRNA_sadj2(outReg);

on_mean = mean(on_cnts);
on_std = std(on_cnts);
off_mean = mean(off_cnts);


over1 = find(mRNA_sadj2>on_mean*spread);
under1 = find(mRNA_sadj2<on_mean/spread); 

Over1 = ismember(NucLabeled,over1).*onMask; 
Under1 = ismember(NucLabeled,under1).*onMask; 

I = uint8(zeros(h,w,3));
I(:,:,1) =  C1;
I(:,:,2) = C1 - uint8(255*Over1);
I(:,:,3) = uint8(255*Under1)+C1 -uint8(255*Over1);
figure(4); subplot(1,2,2); imshow(I);
title(['y  mean=',num2str(on_mean,4),' std=',num2str(on_std,4),' cov=',num2str(on_std/on_mean,3)]);
 hold on;
plot(bndry1{1}(:,2),bndry1{1}(:,1),'c');



m = linspace(0,cmax,30);
figure(2); 
subplot(2,2,3); hist(on_cnts,m);
title(['y  mean=',num2str(on_mean,4),' std=',num2str(on_std,4),' cov=',num2str(on_std/on_mean,3)]);
subplot(2,2,4); hist(off_cnts,m); title(['mean=',num2str(off_mean,4)]);

 
 
 

 
 
%% Orient image along major axis      
    bw = imresize(uint16(cell_sadj),.25); 
    thresh = graythresh(bw);
    bw = im2bw(bw,thresh); 
    figure(8); clf; imshow(bw);

    L = bwlabel(bw);
    rprops = regionprops(L,'MajorAxis','Orientation');
    NucLabeled = imrotate(NucLabeled,360-rprops.Orientation,'nearest'); 
    figure(8); clf; imagesc(NucLabeled);
    
    
%% convert nucleus centroids to indexed postions
S = regionprops(NucLabeled,'Centroid');
nuc_cents = reshape([S.Centroid],2,length(S));
[h,w] = size(NucLabeled); 
c_inds = sub2ind([h,w],floor(nuc_cents(2,:)),floor(nuc_cents(1,:)));
% compute distances

% c_inds = sub2ind([h,w],floor(nuc_cents(2,:)*h/hn),floor(nuc_cents(1,:)*w/wn));
% compute distances
%d = sqrt( (nuc_cents(1,:)*w/wn).^2 + (nuc_cents(2,:)*h/hn).^2);
d = nuc_cents(2,:)*h/hn;



%  % Plotting for troubleshooting
%     C = false(h,w);
%     C(c_inds) = 1; 
%     figure(1); clf; imshow(C);

% Sort by distance from upper left corner (max bcd).  
nuc_order = NucLabeled(c_inds);
[b,m,n] = unique(nuc_order);
dists = d(m);
figure(1); clf; plot(dists,mRNA_sadj,'k.');  % check results
Data = [dists; mRNA_sadj]';
Data_sort = sortrows(Data); 


% Convert distances to bcd concentrations.  
bcd = log(dists); % just logorithmic.  Don't need arbitrary scaling coeff.
%figure(1); clf; plot(bcd,mRNA_sadj,'k.');



Data2 = [bcd; mRNA_sadj; dists]';
Data2 = sortrows(Data2);

%%
Sects = round(sqrt(Nnucs));
mu = zeros(1,Sects);
sigma = zeros(1,Sects);

Q = cell(1,Sects); 

dbnd = zeros(Sects,1);
for j=1:Sects
    Q{j} = Data2( floor((j-1)*Nnucs/Sects) + 1: floor(j*Nnucs/Sects), 2);
    mu(j) = nanmean(Q{j});
    sigma(j) = nanstd(Q{j});
    fano = sigma(j)^2/mu(j); 
end

x = linspace(0,max(dists)*50/1000,Sects);


figure(1); clf; colordef white; set(gcf,'color','w');

plot(Data_sort(:,1)*50/1000,Data_sort(:,2),'k'); % check results 
hold on; errorbar(x,mu,sigma,'linestyle','none','linewidth',3,'color','r');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');


figure(2); clf; 
colordef white; set(gcf,'color','w');
plot(mu,sigma./mu,'k'); ylim([0,1]);
ylabel('CoV'); xlabel('mean count');

figure(3); clf; 
colordef white; set(gcf,'color','w');
plot(x,sigma./mu,'k'); ylim([0,1]);
ylabel('CoV'); xlabel('distance (\mum)');



% hold on;
%     scatter(nuc_cents(1,:),nuc_cents(2,:)); 
%     
%     M = NucLabeled - 150*C;
%     figure(2); clf; imagesc(M);hold on;
%         scatter(nuc_cents(1,:),nuc_cents(2,:)); 
