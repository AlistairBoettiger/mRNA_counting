
%% anlz_hb_gradient_data.m

% Alistair Boettiger                                Date Begun: 02/02/11
% Levine Lab                                    Last Modified: 03/30/11

% uses radial distance.  should use planar distance.
% Will work better if original images are oriented along the AP axis.  

clear all;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
%fname = 'hbCent_hb_LacZ_v2';
fname = 'gtX_hb_grad_01';
% fname = MP09_22C_hb_y_01
% fname = MP09_22C_hb_y_d_01


load([folder,fname,'.mat']);
[h,w] = size(NucLabeled);
[hn,wn] = size(In); % size associated with nuc centroid positions. 
Nnucs =  max(NucLabeled(:)); 

 % Plotting for trouble shooting
    figure(1); clf; hist(mRNA_sadj,30);
    figure(2); clf; imagesc(NucLabeled); 
    hold on; plot(nuc_cents(1,:)*w/wn,nuc_cents(2,:)*h/hn,'wo');

  
    
    
    % Plot mRNA count per cell, size normalized, and cell area.  
        cell_cnt = zeros(h,w); 
        cell_sadj = zeros(h,w); 
        cell_area = zeros(h,w);      
        reg_data = regionprops(NucLabeled,'PixelIdxList');
            for k=1:Nnucs
                pixes = reg_data(k).PixelIdxList;             
                cell_cnt(pixes) = mRNA_cnt(k);
                cell_sadj(pixes) = mRNA_sadj(k);
                cell_area(pixes) = nuc_area(k);
            end
        colordef black;
        figure(2); clf; imagesc(cell_cnt); colormap('jet'); colorbar; set(gcf,'color','k');
        figure(3); clf; imagesc(cell_sadj); colormap('jet'); colorbar; set(gcf,'color','k');
        figure(4); clf; imagesc(cell_area); colormap('jet'); colorbar; set(gcf,'color','k');%
  %% New data-method
  clear all;
   
  folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
  fname = 'MP09_22C_hb_y_d';
  emb = '05';  i = str2double(emb); cor = 1; 

  chn = 1;
  
  load([folder,fname,'_',emb,'_nucdata.mat']); 
  load([folder,fname,'_slidedata'], 'Data'); 
   
              mRNAsadj = Data{i,chn}.mRNAsadj;

  PlotmRNA = imresize(NucLabeled,.5,'nearest');
  NucLabel = imresize(NucLabeled,.5,'nearest'); 
                      Nnucs =    max( NucLabeled(:) );
                      for n=1:Nnucs;
                          PlotmRNA(PlotmRNA==n) = mRNAsadj(n+cor);
                      end
    figure(1); clf; imagesc(PlotmRNA); colormap hot; 
  
  
        
%% Orient image along major axis     
    bw = imresize(uint16(PlotmRNA),.2); 
    
    thresh = graythresh(bw);
    bw = im2bw(bw,thresh); 
    bw = bwareaopen(bw,100);
    figure(8); clf; imshow(bw);

    L = bwlabel(bw);
    rprops = regionprops(L,'MajorAxis','Orientation');
    NucLabel = imrotate(NucLabel,180-rprops.Orientation,'nearest'); 
    figure(8); clf; imagesc(NucLabel);
    
    
%% convert nucleus centroids to indexed postions
S = regionprops(NucLabel,'Centroid');
nuc_cents = reshape([S.Centroid],2,length(S));
[h,w] = size(NucLabel); 
[hn,wn] = size(NucLabel); 
c_inds = sub2ind([h,w],floor(nuc_cents(2,:)),floor(nuc_cents(1,:)));
% compute distances

% c_inds = sub2ind([h,w],floor(nuc_cents(2,:)*h/hn),floor(nuc_cents(1,:)*w/wn));
% compute distances
%d = sqrt( (nuc_cents(1,:)*w/wn).^2 + (nuc_cents(2,:)*h/hn).^2);
d = nuc_cents(2,:)*h/hn;

if length(c_inds) < length(mRNAsadj);
    mRNAsadj = mRNAsadj(2:end);
end


maternal = min(mRNAsadj);

%  % Plotting for troubleshooting
%     C = false(h,w);
%     C(c_inds) = 1; 
%     figure(1); clf; imshow(C);

% Sort by distance from upper left corner (max bcd).  
nuc_order = NucLabel(c_inds);
[b,m,n] = unique(nuc_order);
dists = d(m);
figure(1); clf; plot(dists,mRNAsadj - maternal,'g.');  % check results
Data = [dists; mRNAsadj - maternal]';
Data_sort = sortrows(Data); 


% Convert distances to bcd concentrations.  
bcd = log(dists); % just logorithmic.  Don't need arbitrary scaling coeff.
%figure(1); clf; plot(bcd,mRNA_sadj,'k.');



Data2 = [bcd; mRNAsadj-maternal; dists]';
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


figure(1); clf; colordef black; set(gcf,'color','k');

plot(Data_sort(:,1)*50/1000,flipud(Data_sort(:,2)),'w.'); % check results 
hold on; errorbar(x,fliplr(mu),fliplr(sigma),'linestyle','none','linewidth',3,'color','r');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');
set(gca,'FontSize',14);


figure(1); clf; colordef white; set(gcf,'color','w');

plot(Data_sort(:,1)*50/1000,Data_sort(:,2),'k.'); % check results 
hold on; errorbar(x,mu,sigma,'linestyle','none','linewidth',3,'color','r');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');


figure(5); clf; 
colordef white; set(gcf,'color','w');
plot(mu,sigma,'k.'); hold on; plot(mu,sigma.*sqrt([0,diff(mu)]),'g.'); 
 xlabel('mean count');

figure(2); clf; 
colordef white; set(gcf,'color','w');
plot(x,sigma,'r.'); 
hold on; plot(x,mu,'m.'); legend('\sigma','\mu','Location','Best');
xlabel('distance (\mum)');

figure(3); clf; 
colordef white; set(gcf,'color','w');
plot(diff(sigma),'r-'); 
hold on; plot(diff(mu),'m-'); legend('\sigma','\mu');
xlabel('distance (\mum)');




figure(4); clf; 
colordef white; set(gcf,'color','w');
plot(x,sigma./mu,'k.'); ylim([0,1]);
ylabel('CoV'); xlabel('distance (\mum)');

%%

fimage = 'gtX_hb_grad_01_01_max.tif';
dataf = '/Volumes/Data/Lab Data/Raw_Data/02-03-11/';

I = imread([dataf,fimage]);
[h,w] = size(I(:,:,1));



I2 = uint16(zeros(h,w,3));
I2(:,:,1) = 2*I(:,:,1);
I2(:,:,2) = I(:,:,3);
I2(:,:,3) = I(:,:,3);

 I2 = imrotate(I2,90);

figure(4); clf; image(I2(500:1500,:,:));
hold on;

ysc = 1;
xsc = 4.65;

plot(50+xsc*Data_sort(:,1)*50/1000,ysc*flipud(Data_sort(:,2)),'w.','MarkerSize',10); % check results 
hold on; errorbar(50+xsc*x,ysc*fliplr(mu),ysc*fliplr(sigma),'linestyle','none','linewidth',3,'color','r');
axis on; axis xy;
ylabel('number of mRNA transcripts per cell','FontSize',14); 
xlabel('distance (\mum)','FontSize',14);
set(gca,'FontSize',14);






