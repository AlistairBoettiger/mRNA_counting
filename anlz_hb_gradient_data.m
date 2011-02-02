
%% anlz_hb_gradient_data.m


clear all;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
fname = 'hbCent_hb_LacZ_v2';

load([folder,fname,'.mat']);
[h,w] = size(NucLabeled);
[hn,wn] = size(In); % size associated with nuc centroid positions. 
Nnucs =  max(NucLabeled(:)); 

%  % Plotting for trouble shooting
%     figure(1); clf; hist(mRNA_sadj,30);
%     figure(2); clf; imagesc(NucLabeled); 
%     hold on; plot(nuc_cents(1,:)*w/wn,nuc_cents(2,:)*h/hn,'wo');

    
% convert nucleus centroids to indexed postions
c_inds = sub2ind([h,w],floor(nuc_cents(2,:)*h/hn),floor(nuc_cents(1,:)*w/wn));
% compute distances
d = sqrt( (nuc_cents(1,:)*w/wn).^2 + (nuc_cents(2,:)*h/hn).^2);




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
figure(4); clf; imagesc(cell_area); colormap('jet'); colorbar; set(gcf,'color','k');
      



%  % Plotting for troubleshooting
%     C = false(h,w);
%     C(c_inds) = 1; 
%     figure(1); clf; imshow(C);

% Sort by distance from upper left corner (max bcd).  
nuc_order = NucLabeled(c_inds);
[b,m,n] = unique(nuc_order);
dists = d(m);
figure(1); clf; plot(dists,mRNA_sadj,'y.');  % check results
Data = [dists; mRNA_sadj]';
Data_sort = sortrows(Data); 
figure(1); clf; plot(Data_sort(:,1),Data_sort(:,2)); % check results 


% Convert distances to bcd concentrations.  
bcd = log(dists); % just logorithmic.  Don't need arbitrary scaling coeff.
figure(1); clf; plot(bcd,mRNA_sadj,'y.');

Data2 = [bcd; mRNA_sadj; dists]';
Data2 = sortrows(Data2);

t = linspace(0,2*pi,100);
Sects = 6;
xc = zeros(Sects,100);
yc = zeros(Sects,100);
Q = cell(1,Sects); 
dbnd = zeros(Sects,1);
for j=1:Sects
    Q{j} = Data2( floor((j-1)*Nnucs/Sects) + 1: floor(j*Nnucs/Sects), 2);
    r = Data2(floor(j*Nnucs/Sects), 3);
    xc(j,:) = r*cos(t);
    yc(j,:) = r*sin(t); 
    figure(3); hold on; plot(xc(j,:),yc(j,:),'w','linewidth',3);
end

figure(1); clf;
mRNA = linspace(0,1400,15);
for j=1:Sects
    subplot(Sects,1,j);
    hist(Q{j},mRNA); xlim([0,1500]);
end








% hold on;
%     scatter(nuc_cents(1,:),nuc_cents(2,:)); 
%     
%     M = NucLabeled - 150*C;
%     figure(2); clf; imagesc(M);hold on;
%         scatter(nuc_cents(1,:),nuc_cents(2,:)); 
