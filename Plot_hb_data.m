
%%                  
% Alistair Boettiger                                   Date Begun: 06/10/11
% Levine Lab                                        Last Modified: 06/27/11

clear all;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/2011-05-22/'; 
slides = {'MP01b','MP01','MP02b','MP02','MP09b','MP09c'}; 

 slides = {'MP09b','MP09c'}; 
 slides = {'MP01b','MP01'}

chns = 2; ver = '';

figure(1); clf; figure(21); clf; figure(22); clf; set(gcf,'color','w');
j = 0; NucsT = {};

S = length(slides);

pars = cell(1,S);


for s = 1:S

    switch slides{s}
        case 'MP01'
        fname = 's02_MP01_Hz_22C';  ver = '_v2';
        % MP01 = hbdata;
        e=4;% 4; 
        e1=e;
      %  o1 = 4E4;

        case 'MP01b'
        fname = 's02_MP01_Hz_22C_b';  ver = '';
        % MP01_b = hbdata;

        case 'MP02'
        fname = 's03_MP02_Hz_22C';  ver = '';
        % MP02 = hbdata;
        e=1; % 8;
        e2=e;
       %  o2 = 0;

        case 'MP02b'
        fname = 's03_MP02_Hz_22C_b';  ver = '_v2'; % little difference
        % MP02_b = hbdata;

        case 'MP09b'
        fname = 's01_MP09_Hz_22C_b'; ver = '_v2';
        % MP09 = hbdata;
            e = 5; % 5 9; 
            e9 = e;
         %   o9 = 4E4;

        case 'MP09'
        fname = 's01_MP09_Hz_22C';  ver = '_v2';
        % MP09_a = hbdata;

        case 'MP09c'
        fname = 's01_MP09_Hz_22C_c'; ver = '_v2';
        % MP09_c = hbdata;
    end
    
    
    load([folder,fname,'_graddata',ver],'hbdata'); 
    % Show what were working on
    disp(fname);
    disp(['# embryos in dataset = ',num2str(sum(1-cellfun('isempty',hbdata)))])
    try
        disp(['chn1 thresh: ', num2str(hbdata{e}.ipars{1}.min_int)])
        disp(['chn2 thresh: ', num2str(hbdata{e}.ipars{2}.min_int)])
        disp(['missG: ', num2str(hbdata{e}.ipars{3})])
    catch er
        disp('processing parameter data not available'); 
    end
    
    %
    Es = length(hbdata);
    offsets = zeros(1,Es); 
    pars{s} = zeros(Es,4); % store curve fit parameters for alignment.  

    % Get boundary point of first curve as a reference
    if s == 1;
        e=1;
        gx = hbdata{e}.Data_sort(:,1);  
        grad = hbdata{e}.Data_sort(:,2);  

        [pars{s}(e,:),fit] = fxn_fit_sigmoid(gx',grad',[4,mean(gx),max(grad),min(grad)],'r');
        offsets(e) = pars{s}(e,2); 
        p1 = pars{s}(e,2);

        figure(20); clf; 
        plot(gx', grad,'b');
        hold on;  
        plot(gx,fit,'r');
        plot(p1,max(grad)/4,'r*','MarkerSize',20);
        k = 2;
    else
        k=1;
    end  

    for e= k:Es;  % Get curve fit parameters for the rest of the curves
        try
             gx = hbdata{e}.Data_sort(:,1);
             grad = hbdata{e}.Data_sort(:,2);  

            figure(20); clf; plot(gx,grad);
            [pars{s}(e,:),fit] = fxn_fit_sigmoid(gx',grad',[4,mean(gx),max(grad),min(grad)],'r');
            offsets(e) = pars{s}(e,2); 

            figure(20); clf; plot(gx', grad,'b');
            hold on; 
            plot(gx,fit,'r');
            plot(offsets(e),max(grad)/4,'r*','MarkerSize',20);
            
      if  pars{s}(e,3)+ pars{s}(e,4) < 370 &&  pars{s}(e,3)+ pars{s}(e,4) > 210 && pars{s}(e,1) > 3.5 && pars{s}(e,1) < 6  % 
          % pars{s}(e,3)+ pars{s}(e,4) < 250 &&  pars{s}(e,3)+ pars{s}(e,4) > 150 && pars{s}(e,1) > 3.5  && pars{s}(e,1) < 6
          %
          figure(21); hold on;
               plot(gx+ p1-offsets(e),fit,'r');
               disp(pars{s}(e,3))           
      end
            
        catch er
            disp(er.message); 
            continue
        end
    end


   %   figure(s);  clf;

     % figure(1);
    colordef white; set(gcf,'color','w');
   % i = 1; Nucs = {};
    for e = 1:Es 
       % if  isempty(hbdata{e}) ~= 1;
            if  pars{s}(e,3)+ pars{s}(e,4) < 370 &&  pars{s}(e,3)+ pars{s}(e,4) > 210 && pars{s}(e,1) > 3.5 && pars{s}(e,1) < 6  % pars{s}(e,3)+ pars{s}(e,4) < 390 &&  pars{s}(e,3)+ pars{s}(e,4) > 310 && pars{s}(e,1) > 3.5  && pars{s}(e,1) < 6

                % pars{s}(e,3)+ pars{s}(e,4) < 250 &&  pars{s}(e,3)+ pars{s}(e,4) > 150 && pars{s}(e,1) > 3.5  && pars{s}(e,1) < 6

                j = j+1;     
               figure(22); hold on; 
               mcorr =  1; %   hbdata{e}.Nnucs./hbdata{1}.Nnucs;
               plot(hbdata{e}.Data_sort(:,1) + p1-offsets(e),mcorr*hbdata{e}.Data_sort(:,2),'.','color',[1-e/Es,0,e/(2*Es)],'MarkerSize',5); % check results  
               hold on; 
               ylim([0,550]); 
               % errorbar(hbdata{e}.x+ p1-offsets(e),mcorr*hbdata{e}.mu(:,1),mcorr*hbdata{e}.sigma(:,1),'linestyle','none','linewidth',1,'color',[e/Es,0,1-e/Es],'MarkerSize',1);
              % Nucs{i} =['wt',num2str(e),' N=', num2str( hbdata{e}.Nnucs) ];
                 NucsT{j} = ['wt',' emb',num2str(e),' N=', num2str( hbdata{e}.Nnucs) ];

                % Nucs{1+2*(i-1)} =['wt',num2str(e),' N=', num2str( hbdata{e}.Nnucs) ];
                % i = i + 1; 
                if chns == 2;
                     figure(22); hold on; 
                    plot(hbdata{e}.Data_sort(:,1)+  p1-offsets(e),mcorr*hbdata{e}.Data_sort(:,3),'.','color',[s/S,1-e/Es,e/(s*Es)],'MarkerSize',5); 
                    hold on; 
               %     errorbar(hbdata{e}.x+ p1-offsets(e),mcorr*hbdata{e}.mu(:,2),mcorr*hbdata{e}.sigma(:,2),'linestyle','none','linewidth',1,'color',[s/S,1-e/Es,e/(s*Es)],'MarkerSize',1);
                    j = j+1;
                    NucsT{j} = [slides{s},' emb',num2str(e),' N=', num2str( hbdata{e}.Nnucs) ];
                end
            end
      %  end
    end
%             ylabel('number of mRNA transcripts per cell'); xlabel('distance (nm)');
%            legend(Nucs{:});
%             title(texlabel(fname,'literal')); 

end
            figure(22); 
            ylabel('number of mRNA transcripts per cell'); xlabel('distance (nm)');
           legend(NucsT{:});
            title(texlabel(fname,'literal')); 
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
xlabel('distance (nm)','Fontsize',F);

figure(4); 
legend(Nucs_13{:});
ylabel('number of mRNA transcripts per cell','Fontsize',F); 
xlabel('distance (nm)','Fontsize',F);

figure(5); 
legend(Nucs_14a{:});
ylabel('number of mRNA transcripts per cell','Fontsize',F); 
xlabel('distance (nm)','Fontsize',F);
    
    
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
    xlabel('distance (nm)','Fontsize',F); hold on;
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
ylabel('number of mRNA transcripts per cell'); xlabel('distance (nm)');
xlim([0,2E5]);  ylim([0,400]);

subplot(1,3,2); 
e = e2;
plot(MP02{e}.Data_sort(:,1),MP02{e}.Data_sort(:,2),'k.'); % check results  
title(['Nuclei = ',num2str(MP02{e}.Nnucs)]);
hold on;
errorbar(MP02{e}.x,MP02{e}.mu(:,1),MP02{e}.sigma(:,1),'linestyle','none','linewidth',3,'color','r');
plot(MP02{e}.Data_sort(:,1),MP02{e}.Data_sort(:,3),'m.'); 
errorbar(MP02{e}.x,MP02{e}.mu(:,2),MP02{e}.sigma(:,2),'linestyle','none','linewidth',3,'color','c');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (nm)');
xlim([0,2E5]); ylim([0,400]);

subplot(1,3,3); 
e =  e9;
plot(MP09{e}.Data_sort(:,1),MP09{e}.Data_sort(:,2),'k.'); % check results  
title(['Nuclei = ',num2str(MP09{e}.Nnucs)]);
hold on;
errorbar(MP09{e}.x,MP09{e}.mu(:,1),MP09{e}.sigma(:,1),'linestyle','none','linewidth',3,'color','r');
plot(MP09{e}.Data_sort(:,1),MP09{e}.Data_sort(:,3),'b.'); 
errorbar(MP09{e}.x,MP09{e}.mu(:,2),MP09{e}.sigma(:,2),'linestyle','none','linewidth',3,'color','c');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (nm)');
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
xlabel('distance (nm)','Fontsize',F);

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

