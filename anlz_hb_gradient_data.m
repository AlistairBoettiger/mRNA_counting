
%% anlz_hb_gradient_data.m

% Alistair Boettiger                                Date Begun: 02/02/11
% Levine Lab                                    Last Modified: 02/02/11

% uses radial distance.  should use planar distance.
% Will work better if original images are oriented along the AP axis.  

clear all;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
%fname = 'hbCent_hb_LacZ_v2';
%fname = 'gtX_hb_grad_01';
fname = 'MP05_sna_y_22C_theshp1';

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
