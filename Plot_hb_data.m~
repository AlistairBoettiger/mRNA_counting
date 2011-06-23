
%%                  
% Alistair Boettiger                                   Date Begun: 06/10/11
% Levine Lab                                        Last Modified: 06/10/11

clear all;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/2011-05-22/'; 



fname = 's02_MP01_Hz_22C'; 
load([folder,fname,'_graddata'],'hbdata'); 
MP01 = hbdata;
e=4;% 4; 
e1=e;
o1 = 4E4;

fname = 's02_MP01_Hz_22C_b'; 
load([folder,fname,'_graddata'],'hbdata'); 
MP01_b = hbdata;

fname = 's03_MP02_Hz_22C'; 
load([folder,fname,'_graddata'],'hbdata'); 
MP02 = hbdata;
e=1; % 8;
e2=e;
o2 = 0;

fname = 's03_MP02_Hz_22C_b'; 
load([folder,fname,'_graddata'],'hbdata'); 
MP02_b = hbdata;


fname = 's01_MP09_Hz_22C_b'; 
load([folder,fname,'_graddata'],'hbdata'); 
MP09 = hbdata;
    e = 5; % 5 9; 
    e9 = e;
    o9 = 4E4;

fname = 's01_MP09_Hz_22C'; 
load([folder,fname,'_graddata'],'hbdata'); 
MP09_a = hbdata;

fname = 's01_MP09_Hz_22C_c'; 
load([folder,fname,'_graddata'],'hbdata'); 
MP09_c = hbdata;
    
%%  
figure(3); clf; figure(4); clf; figure(5); clf;

F = 14;
Es = 12;
%       1  2  3  4  5  6   7   8  9  10  11  12
ofs9 = [0, 0, 0, 1.2E4, 0,2.4E4,2.4E4,0,0,0, 2E4,0];
ofs9_a = [0, 0, 0, 0E4, 0,0E4,0E4,0,0,0, 0E4,0];

i=0; ii =0; ia = 0;
for e=1:Es
    try
        corr_mRNA =  1% MP09{e}.Nnucs/MP09{1}.Nnucs;
        if MP09{e}.Nnucs > 200 && MP09{e}.Nnucs < 330;
            i=i+1;
            figure(3);
            plot(MP09{e}.Data_sort(:,1)+ofs9(e),corr_mRNA*MP09{e}.Data_sort(:,2),'.','color',[i/Es,0,1-i/Es]); hold on;
            Nucs_14{1+2*(i-1)} =['wt',num2str(e),' N=', num2str( MP09{e}.Nnucs) ];
            plot(MP09{e}.Data_sort(:,1)+ofs9(e),corr_mRNA*MP09{e}.Data_sort(:,3),'.','color',[0,1-i/Es,i/Es]);
            Nucs_14{2+2*(i-1)} =['y',num2str(e),' N=', num2str( MP09{e}.Nnucs) ];
        end
     
               
        if MP09{e}.Nnucs > 100 && MP09{e}.Nnucs < 200;
            ii=ii+1;
            figure(4);
            plot(MP09{e}.Data_sort(:,1)+ofs9(e),corr_mRNA*MP09{e}.Data_sort(:,2),'.','color',[ii/Es,0,1-ii/Es]); hold on;
            Nucs_13{1+2*(ii-1)} =['wt',num2str(e),' N=', num2str( MP09{e}.Nnucs) ];
             plot(MP09{e}.Data_sort(:,1)+ofs9(e),corr_mRNA*MP09{e}.Data_sort(:,3),'.','color',[0,1-ii/Es,ii/Es]);
            Nucs_13{2+2*(ii-1)} =['y',num2str(e),' N=', num2str( MP09{e}.Nnucs) ];
        end
        
        
    catch
        % Nucs{e} = 'miss';
        continue
    end
    
    
    
    try
        if MP09_c{e}.Nnucs > 100 && MP09_c{e}.Nnucs < 330;
            ia=ia+1;
            figure(5);
            plot(MP09_c{e}.Data_sort(:,1)+ofs9_a(e),corr_mRNA*MP09_c{e}.Data_sort(:,2),'.','color',[ia/Es,0,1-ia/Es]); hold on;
            Nucs_14a{1+2*(ia-1)} =['wt',num2str(e),' N=', num2str( MP09_c{e}.Nnucs) ];
            plot(MP09_c{e}.Data_sort(:,1)+ofs9(e),corr_mRNA*MP09_c{e}.Data_sort(:,3),'.','color',[0,1-ia/Es,ia/Es]);
            Nucs_14a{2+2*(ia-1)} =['y',num2str(e),' N=', num2str( MP09_c{e}.Nnucs) ];
        end
    catch
        continue     
    end
    
        
end
figure(3); 
legend(Nucs_14{:});
ylabel('number of mRNA transcripts per cell','Fontsize',F); 
xlabel('distance (\mum)','Fontsize',F);

figure(4); 
legend(Nucs_13{:});
ylabel('number of mRNA transcripts per cell','Fontsize',F); 
xlabel('distance (\mum)','Fontsize',F);

figure(5); 
legend(Nucs_14a{:});
ylabel('number of mRNA transcripts per cell','Fontsize',F); 
xlabel('distance (\mum)','Fontsize',F);
    
    
%% Correct alignment AP


e =e1;   
     MP01{e}.x = MP01{e}.x + o1; 
     MP01{e}.Data_sort(:,1)= MP01{e}.Data_sort(:,1)+ o1; 

e =e2;
        MP02{e}.x = MP02{e}.x + o2;
          MP02{e}.Data_sort(:,1)= MP02{e}.Data_sort(:,1)+ o2; 

e =e9;
        MP09{e}.x = MP09{e}.x +o9;
         MP09{e}.Data_sort(:,1) = MP09{e}.Data_sort(:,1)+ o9 ;

%% Plot CoV

figure(1); clf; 
colordef white; set(gcf,'color','w');  F = 14;

e=e1; 
% plot(MP01{e}.x,sqrt(MP01{e}.mu(:,1))./MP01{e}.mu(:,1),'k.','MarkerSize',10);
errorbar(MP01{e}.x,MP01{e}.sigma(:,2)./MP01{e}.mu(:,2),MP01{e}.bssigma(:,2)./MP01{e}.mu(:,2),  'g*','MarkerSize',10); ylim([0,1]);
hold on;
e = e2;
% plot(MP02{e}.x,sqrt(MP02{e}.mu(:,1))./MP02{e}.mu(:,1),'k.','MarkerSize',10);
errorbar(MP02{e}.x,MP02{e}.sigma(:,2)./MP02{e}.mu(:,2),MP02{e}.bssigma(:,2)./MP02{e}.mu(:,2),  'm*','MarkerSize',10); ylim([0,1]);
hold on;
    e = e9; 
errorbar(MP09{e}.x,MP09{e}.sigma(:,2)./MP09{e}.mu(:,2),MP09{e}.bssigma(:,2)./MP09{e}.mu(:,2),  'b.','MarkerSize',10); ylim([0,1]);

% The endogenous variances 
e=e1; errorbar(MP01{e}.x,MP01{e}.sigma(:,1)./MP01{e}.mu(:,1),MP01{e}.bssigma(:,1)./MP01{e}.mu(:,1),  'k.','MarkerSize',10); 
e = e2; errorbar(MP02{e}.x,MP02{e}.sigma(:,1)./MP02{e}.mu(:,1),MP02{e}.bssigma(:,1)./MP02{e}.mu(:,1),  'k.','MarkerSize',10); 
e = e9; errorbar(MP09{e}.x,MP09{e}.sigma(:,1)./MP09{e}.mu(:,1),MP09{e}.bssigma(:,1)./MP09{e}.mu(:,1),  'k.','MarkerSize',10); 

 
 % plot(MP09{e}.x,sqrt(MP09{e}.mu(:,1))./MP09{e}.mu(:,1),'k.','MarkerSize',10);
    ylabel('CoV','FontSize',F);
    xlabel('distance (\mum)','Fontsize',F); hold on;
    set(gca,'Fontsize',F);
    legend(['no proximal, N=',num2str(MP01{e1}.Nnucs)],...
    ['no distal, N=',num2str(MP02{e2}.Nnucs)],...
    ['cntrl N=',num2str(MP09{e9}.Nnucs)],...
    'Combined endogenous');

%% Separate plot of gradients
figure(2); clf; set(gcf,'color','w');

subplot(1,3,1); 
e=e1;
plot(MP01{e}.Data_sort(:,1),MP01{e}.Data_sort(:,2),'k.'); % check results  
title(['Nuclei = ',num2str(MP01{e}.Nnucs)]);
hold on;
errorbar(MP01{e}.x,MP01{e}.mu(:,1),MP01{e}.sigma(:,1),'linestyle','none','linewidth',3,'color','r');
plot(MP01{e}.Data_sort(:,1),MP01{e}.Data_sort(:,3),'g.'); 
errorbar(MP01{e}.x,MP01{e}.mu(:,2),MP01{e}.sigma(:,2),'linestyle','none','linewidth',3,'color','c');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');
xlim([0,2E5]);  ylim([0,400]);

subplot(1,3,2); 
e = e2;
plot(MP02{e}.Data_sort(:,1),MP02{e}.Data_sort(:,2),'k.'); % check results  
title(['Nuclei = ',num2str(MP02{e}.Nnucs)]);
hold on;
errorbar(MP02{e}.x,MP02{e}.mu(:,1),MP02{e}.sigma(:,1),'linestyle','none','linewidth',3,'color','r');
plot(MP02{e}.Data_sort(:,1),MP02{e}.Data_sort(:,3),'m.'); 
errorbar(MP02{e}.x,MP02{e}.mu(:,2),MP02{e}.sigma(:,2),'linestyle','none','linewidth',3,'color','c');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');
xlim([0,2E5]); ylim([0,400]);

subplot(1,3,3); 
e =  e9;
plot(MP09{e}.Data_sort(:,1),MP09{e}.Data_sort(:,2),'k.'); % check results  
title(['Nuclei = ',num2str(MP09{e}.Nnucs)]);
hold on;
errorbar(MP09{e}.x,MP09{e}.mu(:,1),MP09{e}.sigma(:,1),'linestyle','none','linewidth',3,'color','r');
plot(MP09{e}.Data_sort(:,1),MP09{e}.Data_sort(:,3),'b.'); 
errorbar(MP09{e}.x,MP09{e}.mu(:,2),MP09{e}.sigma(:,2),'linestyle','none','linewidth',3,'color','c');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');
xlim([0,2E5]); ylim([0,400]);

%%  Combined plot of gradients
figure(3); clf; F = 14;
colordef white; set(gcf,'color','w');

e = e1; plot(MP01{e}.Data_sort(:,1),MP01{e}.Data_sort(:,3),'g.');  hold on;
e = e2; plot(MP02{e}.Data_sort(:,1),MP02{e}.Data_sort(:,3),'m.');
e = e9; plot(MP09{e}.Data_sort(:,1),MP09{e}.Data_sort(:,3),'b.'); 

e = e1; plot(MP01{e}.Data_sort(:,1),MP01{e}.Data_sort(:,2),'k.'); % check results  
e = e2; plot(MP02{e}.Data_sort(:,1),MP02{e}.Data_sort(:,2),'k.'); % check results  
e = e9; plot(MP09{e}.Data_sort(:,1),MP09{e}.Data_sort(:,2),'k.'); % check results 
ylabel('number of mRNA transcripts per cell','Fontsize',F); 
xlabel('distance (\mum)','Fontsize',F);

legend(['no proximal, N=',num2str(MP01{e1}.Nnucs)],...
    ['no distal, N=',num2str(MP02{e2}.Nnucs)],['cntrl N=',num2str(MP09{e9}.Nnucs)],'Endogenous, combined');
set(gca,'Fontsize',F);
title(['mid cc13 mRNA counts']);


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
% clear Iin hb_map;
% Iin = imread('/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Results/hb_expression.tif');
% figure(5); clf; imagesc(Iin(:,:,1)); colormap hot; axis off; set(gcf,'color','k');
% mask = im2bw(Iin(:,:,1),.03);
% figure(5); clf; imagesc(mask);
% 
% hb_map = mask.*double(100+.5*Iin(:,:,1) -  .5*Iin(:,:,3));
% C1 = colormap(hot); C2 = colormap(1-hot);
% C = [0,0,0; C2; C1(1:end,:)];
% 
% figure(6); clf; imagesc(hb_map); colormap(C); %colorbar;
% axis off; colordef black; set(gcf,'color','k');

