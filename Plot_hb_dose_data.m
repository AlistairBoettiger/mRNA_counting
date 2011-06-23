
%%                  
% Alistair Boettiger                                   Date Begun: 06/10/11
% Levine Lab                                        Last Modified: 06/10/11

clear all;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/2011-05-22/'; 

fname = 's10_bcd1x';
load([folder,fname,'_graddata'],'hbdata'); 
bcd1x = hbdata; 
e=5;% 4; 
e1=e;
o1 = 1.75E4;

fname = 's11_bcd6x';
load([folder,fname,'_graddata'],'hbdata'); 
bcd6x = hbdata; 
e=1; % 8;
e2=e;
o2 = 0;

%%
figure(3); clf;
ofs6 = [7E4,0,3E5,10E4,10E4,3E5];
ofs1 = [3E5,6E4,7E4,9E4,7E4,7E4];
for e=1:6
    corr_mRNA = bcd6x{e}.Nnucs/bcd6x{1}.Nnucs;
   plot(bcd6x{e}.Data_sort(:,1)+ofs6(e),corr_mRNA*bcd6x{e}.Data_sort(:,2),'b.'); hold on;
   
   corr_mRNA = bcd1x{e}.Nnucs/bcd6x{1}.Nnucs;
   plot(bcd1x{e}.Data_sort(:,1)+ofs1(e),corr_mRNA*bcd1x{e}.Data_sort(:,2),'r.'); hold on;
end
 corr_mRNA = bcd6x{3}.Nnucs/bcd6x{1}.Nnucs;
 plot(bcd6x{3}.Data_sort(:,1),corr_mRNA*bcd6x{3}.Data_sort(:,2),'c.'); hold on;
xlim([0,2.1E5]); 
legend('6x bcd, N=4 embryos','1x bcd, N=6 embryos');
ylabel('number of mRNA transcripts per cell','Fontsize',F); 
xlabel('distance (\mum)','Fontsize',F);
%% Correct alignment AP


e =e1;   
     bcd1x{e}.x = bcd1x{e}.x + o1; 
     bcd1x{e}.Data_sort(:,1)= bcd1x{e}.Data_sort(:,1)+ o1; 

e =e2;
        bcd6x{e}.x = bcd6x{e}.x + o2;
        bcd6x{e}.Data_sort(:,1)= bcd6x{e}.Data_sort(:,1)+ o2; 


 %%  Combined plot of gradients
figure(1); clf; F = 14;
colordef white; set(gcf,'color','w');


e = e1; plot(bcd1x{e}.Data_sort(:,1),bcd1x{e}.Data_sort(:,2),'r.'); 
hold on;
e = e2; plot(bcd6x{e}.Data_sort(:,1),bcd6x{e}.Data_sort(:,2),'b.');
ylabel('number of mRNA transcripts per cell','Fontsize',F); 
xlabel('distance (\mum)','Fontsize',F);

legend(['bcd1x, N=',num2str(bcd1x{e1}.Nnucs)],...
    ['bcd6x, N=',num2str(bcd6x{e2}.Nnucs)]);
set(gca,'Fontsize',F);
title(['mid cc13 mRNA counts']);

        
%% Plot CoV

figure(2); clf; 
colordef white; set(gcf,'color','w');  F = 14;

% The endogenous variances 
e=e1;  errorbar(bcd1x{e}.x,bcd1x{e}.sigma(:,1)./bcd1x{e}.mu(:,1),bcd1x{e}.bssigma(:,1)./bcd1x{e}.mu(:,1),  'r.','MarkerSize',10); hold on;
e=e2;  errorbar(bcd6x{e}.x,bcd6x{e}.sigma(:,1)./bcd6x{e}.mu(:,1),bcd6x{e}.bssigma(:,1)./bcd6x{e}.mu(:,1),  'b.','MarkerSize',10); 


 
 % plot(MP09{e}.x,sqrt(MP09{e}.mu(:,1))./MP09{e}.mu(:,1),'k.','MarkerSize',10);
    ylabel('CoV','FontSize',F);
    xlabel('distance (\mum)','Fontsize',F); hold on;
    set(gca,'Fontsize',F);
    legend(['bcd1x, N=',num2str(bcd1x{e1}.Nnucs)],...
    ['bcd6x, N=',num2str(bcd6x{e2}.Nnucs)] );

%% Separate plot of gradients
figure(3); clf; set(gcf,'color','w');

subplot(1,3,1); 
e=e1;
plot(bcd1x{e}.Data_sort(:,1),bcd1x{e}.Data_sort(:,2),'r.'); % check results  
title(['Nuclei = ',num2str(bcd1x{e}.Nnucs)]);
hold on;
errorbar(Mbcd1x{e}.x,bcd1x{e}.mu(:,1),bcd1x{e}.sigma(:,1),'linestyle','none','linewidth',3,'color','r');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');
xlim([0,2E5]);  ylim([0,400]);




%% 

figure(4); clf; colordef black; set(gcf,'color','k');
subplot(2,3,1);  imagesc(MP01{e1}.PlotmRNA2); colormap hot; colorbar; axis off; caxis([50,250]);
title('no proximal');
subplot(2,3,2);  imagesc(MP02{e2}.PlotmRNA2); colormap hot; colorbar; axis off; caxis([50,150]);
title('no distal');
subplot(2,3,3);  imagesc(MP09{e9}.PlotmRNA2); colormap hot; colorbar; axis off; caxis([50,300]);
title('control');
subplot(2,3,4);  imagesc(MP01{e1}.PlotmRNA); colormap hot; colorbar; axis off; caxis([100,275]);
subplot(2,3,5);  imagesc(MP02{e2}.PlotmRNA); colormap hot; colorbar; axis off; caxis([100,350]);
subplot(2,3,6);  imagesc(MP09{e9}.PlotmRNA); colormap hot; colorbar; axis off;caxis([100,350]);




%%
clear Iin hb_map;
Iin = imread('/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_modeling/Results/hb_expression.tif');
figure(5); clf; imagesc(Iin(:,:,1)); colormap hot; axis off; set(gcf,'color','k');
mask = im2bw(Iin(:,:,1),.03);
figure(5); clf; imagesc(mask);

hb_map = mask.*double(100+.5*Iin(:,:,1) -  .5*Iin(:,:,3));
C1 = colormap(hot); C2 = colormap(1-hot);
C = [0,0,0; C2; C1(1:end,:)];

figure(6); clf; imagesc(hb_map); colormap(C); %colorbar;
axis off; colordef black; set(gcf,'color','k');

