
%%                  
% Alistair Boettiger                                   Date Begun: 06/10/11
% Levine Lab                                        Last Modified: 06/10/11

clear all;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/'; 

fname = 's01_MP09_Hz_22C_b'; 
load([folder,fname,'_graddata'],'hbdata'); 
MP09 = hbdata;


fname = 's02_MP01_Hz_22C'; 
load([folder,fname,'_graddata'],'hbdata'); 
MP01 = hbdata;


fname = 's03_MP02_Hz_22C'; 
load([folder,fname,'_graddata'],'hbdata'); 
MP02 = hbdata;

%% Correct orientations and Choose embryos to graph 

e=4;% 4; 
e1=e;
o1 = 4E4;

% Orient Anterior Left
     or1 = MP01{e}.mu(1,1) -  MP01{e}.mu(end,1);
     MP01{e}.x = MP01{e}.x + o1; 
     MP01{e}.Data_sort(:,1)= MP01{e}.Data_sort(:,1)+ o1; 
        if or1<0
           MP01{e}.Data_sort(:,2) = flipud( MP01{e}.Data_sort(:,2));
           MP01{e}.Data_sort(:,3) = flipud( MP01{e}.Data_sort(:,3));
           MP01{e}.mu(:,1) = flipud(MP01{e}.mu(:,1));
           MP01{e}.sigma(:,1) = flipud(MP01{e}.sigma(:,1));
           MP01{e}.mu(:,2) = flipud(MP01{e}.mu(:,2));
           MP01{e}.sigma(:,2) = flipud(MP01{e}.sigma(:,2));
           MP01{e}.bssigma(:,1) = flipud(MP01{e}.bssigma(:,1));
           MP01{e}.bssigma(:,2) = flipud(MP01{e}.bssigma(:,2));
        end

e=1; % 8;
e2=e;
o2 = 0;

    % Orient Anterior Left
        or2 = MP02{e}.mu(1,1) -  MP02{e}.mu(end,1);
        MP02{e}.x = MP02{e}.x + o2;
          MP02{e}.Data_sort(:,1)= MP02{e}.Data_sort(:,1)+ o2; 
        if or2<0
           MP02{e}.Data_sort(:,2) = flipud( MP02{e}.Data_sort(:,2));
           MP02{e}.Data_sort(:,3) = flipud( MP02{e}.Data_sort(:,3));
           MP02{e}.mu(:,1) = flipud(MP02{e}.mu(:,1));
           MP02{e}.sigma(:,1) = flipud(MP02{e}.sigma(:,1));
           MP02{e}.mu(:,2) = flipud(MP02{e}.mu(:,2));
           MP02{e}.sigma(:,2) = flipud(MP02{e}.sigma(:,2));
           MP02{e}.bssigma(:,1) = flipud(MP02{e}.bssigma(:,1));
           MP02{e}.bssigma(:,2) = flipud(MP02{e}.bssigma(:,2));
        end

% MP09 Data: control embryos
    e = 5; % 5 9; 
    e9 = e;
    o9 = 4E4;
    % Orient Anterior Left
        or9 = MP09{e}.mu(1,1) -  MP09{e}.mu(end,1);
        MP09{e}.x = MP09{e}.x +o9;
         MP09{e}.Data_sort(:,1) = MP09{e}.Data_sort(:,1)+ o9 ;
        if or9<0
           MP09{e}.Data_sort(:,2) = flipud( MP09{e}.Data_sort(:,2));
           MP09{e}.Data_sort(:,3) = flipud( MP09{e}.Data_sort(:,3));
           MP09{e}.mu(:,1) = flipud(MP09{e}.mu(:,1));
           MP09{e}.sigma(:,1) = flipud(MP09{e}.sigma(:,1));
           MP09{e}.mu(:,2) = flipud(MP09{e}.mu(:,2));
           MP09{e}.sigma(:,2) = flipud(MP09{e}.sigma(:,2));
           MP09{e}.bssigma(:,1) = flipud(MP09{e}.bssigma(:,1));
           MP09{e}.bssigma(:,2) = flipud(MP09{e}.bssigma(:,2));
        end
    % Plot
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
errorbar(MP01{e}.x,MP01{e}.sigma(:,1)./MP01{e}.mu(:,1),MP01{e}.bssigma(:,1)./MP01{e}.mu(:,1),  'k.','MarkerSize',10); 
errorbar(MP01{e}.x,MP01{e}.sigma(:,1)./MP01{e}.mu(:,1),MP01{e}.bssigma(:,1)./MP01{e}.mu(:,1),  'k.','MarkerSize',10); 
errorbar(MP09{e}.x,MP09{e}.sigma(:,1)./MP09{e}.mu(:,1),MP09{e}.bssigma(:,1)./MP09{e}.mu(:,1),  'k.','MarkerSize',10); 

 
 % plot(MP09{e}.x,sqrt(MP09{e}.mu(:,1))./MP09{e}.mu(:,1),'k.','MarkerSize',10);
    ylabel('CoV','FontSize',F);
    xlabel('distance (\mum)','Fontsize',F); hold on;
    set(gca,'Fontsize',F);
    legend(['no proximal, N=',num2str(MP01{e1}.Nnucs)],...
    ['no distal, N=',num2str(MP02{e2}.Nnucs)],['cntrl N=',num2str(MP09{e9}.Nnucs)]);

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

