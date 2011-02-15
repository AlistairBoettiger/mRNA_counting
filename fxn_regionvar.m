
%                           fxn_regionvar.m
% Alistair Boettiger                                Date Begun: 02/14/11
%                                                   Last Modified: 02/14/11

%% Description
% NucLabeled -- nuclei label matrix
% cell_sadj --  mRNA counts for all (nuclei) objects in Label matrix
% t1 --         bw threshold to define on region
% spread --     fraction of mean cutoff to be defined as over / under expression  
%
%% Called by:
% anlz_snail_dots.m link established 02/14/11



function [on_cnts,off_cnts]= fxn_regionvar(NucLabeled,cell_sadj1,mRNA_sadj1,t1,spread,Nnucs)
%%

% Parameters that might want to be user controlled. 
minObjSize = 2E4;
strel_close = 100;


[h,w] = size(NucLabeled); 
% Automatic Threshold
C1 = uint8(cell_sadj1/max(cell_sadj1(:))*255);
%t1 = graythresh(C1);
bw1 = im2bw(C1,t1); % 
bw1 = imclose(bw1,strel('disk',strel_close)); 
bw1 = imfill(bw1,'holes');
bw1 = bwareaopen(bw1,minObjSize);

bndry1 = bwboundaries(bw1);
% figure(3); clf; subplot(2,2,1); imshow(C1);  hold on;
% plot(bndry1{1}(:,2),bndry1{1}(:,1),'c');
% subplot(2,2,2); imshow(bw1); hold on;


% 
% S = regionprops(NucLabeled,'Centroid');
% cents = reshape([S.Centroid],2,length(S));
% inpolygon(cents(:,1),cents(:,2),bndry1{1}(:,2),bndry1{1}(:,1));


temp = bw1.*NucLabeled;
temp2 = hist(temp(:),1:max(NucLabeled(:)));
%figure(2); clf; hist(temp(:),1:max(NucLabeled(:)));


inReg = find(temp2>1E4); 
% [b,m] = unique(bw1.*NucLabeled);
% inReg  = b(m>2E6);
% inReg(inReg==0) = [];

onMask = ismember(NucLabeled,inReg);

outReg = setdiff(1:Nnucs,inReg);
on_cnts = mRNA_sadj1(inReg);
off_cnts = mRNA_sadj1(outReg);


on_mean = mean(on_cnts);
on_std = std(on_cnts);
off_mean = mean(off_cnts);


% figure(2); clf; 
% subplot(2,2,1); hist(on_cnts,mRNA); 
% title(['sna mean=',num2str(on_mean,4),' std=',num2str(on_std,4),'
% cov=',num2str(on_std/on_mean,3)]);
% subplot(2,2,2); hist(off_cnts,mRNA); title(['mean=',num2str(off_mean,4)]);


% internal plotting
over1 = find(mRNA_sadj1>on_mean*spread);
under1 = find(mRNA_sadj1<on_mean/spread); 

Over1 = ismember(NucLabeled,over1).*onMask; 
Under1 = ismember(NucLabeled,under1).*onMask; 

I = uint8(zeros(h,w,3));
I(:,:,1) =  C1;
I(:,:,2) = C1 - uint8(255*Over1);
I(:,:,3) = uint8(255*Under1)+C1 -uint8(255*Over1);

C = gray(max(mRNA_sadj1));
imshow(I); colormap(C); colorbar; cbfreeze; 
title(['sna  mean=',num2str(on_mean,4),' std=',num2str(on_std,4),...
    ' cov=',num2str(on_std/on_mean,3)]);
set(gcf,'color','k');
 hold on;  plot(bndry1{1}(:,2),bndry1{1}(:,1),'c');