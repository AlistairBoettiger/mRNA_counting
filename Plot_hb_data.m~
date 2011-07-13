
%%                  
% Alistair Boettiger                                   Date Begun: 06/10/11
% Levine Lab                                        Last Modified: 07/12/11

clear all;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/2011-05-22/'; 
slides = {'No prox A'; 'No prox B';'No distal A'; 'No distal B';  'cntrl A'; 'cntrl B'};% ; 'cntrl C' %  {'MP01b';'MP01';'MP02b';'MP02';'MP09';'MP09b';'MP09c'}; 

%slides = {'cntrl A'; 'cntrl B';'No distal A'; 'No distal B';'No prox A'; 'No prox B'  };

 % slides = {'No distal A'; 'No distal B'}

chns = 2; 
ver = ''; vout = '';
oldset = 0; 


    indiv_endog = 1; 

% OPTIONS
Ndatasets = 5;  % N datasets
min_hb = .001; % 35;
max_hb = .025;
max_curves = 50;
subtract_maternal = 1; 
plot_grads = 1; 
plot_fits = 0; % plot curve fitting results
other_plots = 0; 

% ~~~~~~~~~~ INITIALIZE DATA STRUCTURES  ~~~~~~~~~~~~~~~~~ %
S = length(slides);
Tembs = 40; 
j = 0; 

NucsT = cell(max_curves,1);
pcurve = cell(Tembs,2,Ndatasets);
pxdata = cell(Tembs,Ndatasets); 
pars = cell(1,S);

y_med_cov = NaN*zeros(Tembs,S);
hb_med_cov = NaN*zeros(Tembs,S);
yhb_ave = NaN*zeros(Tembs,S); 
yhb_std = NaN*zeros(Tembs,S);

if plot_grads == 1
    figure(21); clf; set(gcf,'color','w'); colordef white;
    figure(22); clf; set(gcf,'color','w'); colordef white;
    figure(23); clf; set(gcf,'color','w'); colordef white;
    figure(24); clf; set(gcf,'color','w'); colordef white;
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

for s = 1:S
    switch slides{s}
        case 'No prox A'  % 'MP01'
        fname = 's02_MP01_Hz_22C';  ver = '_v2';
        skip = 9;
        dataset = 1;

        case 'No prox B' % 'MP01b'
        fname = 's02_MP01_Hz_22C_b';  ver = '_v2'; % 
        skip = 25;
        dataset = 1;

        case 'No distal A' % 'MP02'
        fname = 's03_MP02_Hz_22C';  ver = '';
        skip = 25;
        dataset = 2;

        case 'No distal B' % 'MP02b'
        fname = 's03_MP02_Hz_22C_b';  ver = ''; % _v2 little difference
        skip = 25;
        dataset = 2;

        case 'cntrl A' % 'MP09'
        fname = 's01_MP09_Hz_22C';  ver = '_v3';
        skip = 10;
        dataset = 3;
        
        case 'cntrl B' %'MP09b'
        fname = 's01_MP09_Hz_22C_b'; ver = '_v2';
        skip = 7;
        dataset = 3;

        case 'cntrl C' % 'MP09c'
        fname = 's01_MP09_Hz_22C_c'; ver = '_v2';
        skip = 25;
        dataset = 3;
        
        case '6x bcd'
            
        case '1x bcd'
        
    end
    
    load([folder,fname,ver,'_slidedata',vout],'hbdata'); 
    % load([folder,fname,'_graddata',ver],'hbdata'); 
    
    
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
        e = 1;     
        gx = hbdata{e}.Data_sort(:,1)'; %  + 1E6;  % Need to keep x positive for data fitting funciton to work right
        mcorr =   1; %  hbdata{e}.Nnucs./hbdata{1}.Nnucs;
        
        if subtract_maternal == 1; 
            hb_cnt = mcorr*hbdata{e}.Data_sort(:,2) - min(mcorr*hbdata{e}.Data_sort(:,2)) ;  
        else
             hb_cnt = mcorr*hbdata{e}.Data_sort(:,2);
        end

        [pars{s}(e,:),fit] = fxn_fit_sigmoid(gx,hb_cnt',[3,mean(gx),max(hb_cnt),min(hb_cnt)],'r');
        offsets(e) = pars{s}(e,2); 
        p1 = pars{s}(e,2);

        if plot_fits == 1
            figure(20); clf; 
            plot(gx', hb_cnt,'b');
            hold on;  
            plot(gx,fit,'r');
            plot(p1,max(hb_cnt)/4,'r*','MarkerSize',20);
        end
    end
    
    
    if dataset - oldset ~= 0; % restart embryo counter if we've switched datasets  
        emb = 0;
        disp('starting new dataset'); 
    end
    
    for e= 1:Es;  % Get curve fit parameters for the rest of the curves
        try
             mcorr =   1; %  hbdata{e}.Nnucs./hbdata{1}.Nnucs;
             if subtract_maternal == 1; 
                hb_cnt = mcorr*hbdata{e}.Data_sort(:,2) - min(mcorr*hbdata{e}.Data_sort(:,2)) ;  
                hb_m = mcorr*hbdata{e}.mu(:,1)  - min(mcorr*hbdata{e}.Data_sort(:,2));
            else
                 hb_cnt = mcorr*hbdata{e}.Data_sort(:,2);
                 hb_m = mcorr*hbdata{e}.mu(:,1);
            end   
             gx = hbdata{e}.Data_sort(:,1) ;%   + 1E6; % Need to keep x positive for data fitting funciton to work right
            [pars{s}(e,:),fit] = fxn_fit_sigmoid(gx',hb_cnt',[4,mean(gx),max(hb_cnt),min(hb_cnt)],'r');
            if s ==1 && e == 1
            else
                offsets(e) = pars{s}(e,2); % s = 1; e=7
 
              % % % % Manual work around for bug at the moment.    
                if pars{s}(e,2) > 2E5
                   offsets(e) = pars{s}(e,2) -4.2E5;
                   disp('offset error');
                end
            end
            
            if plot_fits ==1
                figure(20); clf; plot(gx',hb_cnt,'b');
                hold on; 
                plot(gx,fit,'r');
                plot(offsets(e),max(hb_cnt)/4,'r*','MarkerSize',20);
            end
            
            cond =hbdata{e}.Nnucs<250 && hbdata{e}.Nnucs>150 &&  pars{s}(e,1) > 2 && pars{s}(e,3)+ pars{s}(e,4) <11E-3 && pars{s}(e,3)+ pars{s}(e,4) >3E-3 ;  
            %cond =hbdata{e}.Nnucs<250 && hbdata{e}.Nnucs>150 &&  pars{s}(e,1) > 2 && pars{s}(e,3)+ pars{s}(e,4) <20E-3 && pars{s}(e,3)+ pars{s}(e,4) >10E-3 ;  
            
            
        catch er
            disp(er.message); 
            continue
        end
   
            if  isempty(find(e==skip,1))  && cond;%  
                j = j+1; 
                emb = emb+1;
                NucsT{j} = ['wt',' emb',num2str(e),' N=', num2str( hbdata{e}.Nnucs) ];
                xdata = hbdata{e}.Data_sort(:,1) + p1-offsets(e)  ; % + 1E6;
               pxdata{emb,dataset} = xdata;
               pcurve{emb,1,dataset} = hb_cnt;
                           
                       if plot_grads == 1
                           figure(21); hold on;
                           plot(gx+ p1-offsets(e),fit,'r');
                           disp(pars{s}(e,3))   

                           figure(22); hold on; 
                           plot(xdata,hb_cnt ,'.','color','k','MarkerSize',5); % check results    % Alt color: [1-s*e/(S*Es),0,s*e/(S*Es)]
                           ylim([0,max_hb]); 


                           figure(24); hold on;
                           hb_e = mcorr*hbdata{e}.sigma(:,1);
                           errorbar(hbdata{e}.x+ p1-offsets(e),hb_m,hb_e,'linestyle','--','linewidth',1,'color','k','MarkerSize',1); 
                           errorbar(hbdata{e}.x+ p1-offsets(e),hb_m,hb_e,'linestyle','none','linewidth',1,'color','k','MarkerSize',1); 
                           hold on; plot(xdata,fit,'--','color',[s/S,0,1-s/S]);
                           ylim([0,max_hb]);
                       end
   
                if chns == 2;
                    j = j+1;
                    NucsT{j} = [slides{s},' emb',num2str(e),' N=', num2str( hbdata{e}.Nnucs)];
                    
                    y_cnt =  mcorr*hbdata{e}.Data_sort(:,3) - min(mcorr*hbdata{e}.Data_sort(:,3));
                    y_m = mcorr*hbdata{e}.mu(:,2)  - min(mcorr*hbdata{e}.Data_sort(:,3));
                    y_e =  mcorr*hbdata{e}.sigma(:,2); 
                    yhb = y_cnt./hb_cnt;
                    pcurve{emb,2,dataset} = y_cnt;        

                            if plot_grads == 1
                                figure(22); hold on; 
                                plot(xdata,y_cnt,'o','color',[s/S,0,1-s/S],'MarkerSize',5); 
                                hold on; 

                                figure(23); hold on; 
                                plot(xdata,yhb,'o','color',[s/S,0,1-s/S],'MarkerSize',5);  % Alt color
                                ylim([0,2]); 

                                figure(24); hold on;
                                  %   errorbar(hbdata{e}.x+ p1-offsets(e),y_m,y_e,'linestyle','-','linewidth',1,'color',[s/S,0,1-s/S],'MarkerSize',1); 
                                errorbar(hbdata{e}.x+ p1-offsets(e),y_m,y_e,'linestyle','none','linewidth',1,'color',[s/S,0,1-s/S],'MarkerSize',1); 
                                hold on; plot(xdata,fit,'color',[s/S,0,1-s/S]);
                            end

                    % Summary statistics
                    yhb_ave(e,s) = nanmean(y_cnt(hb_cnt>min_hb)./hb_cnt(hb_cnt>min_hb));
                    yhb_std(e,s) = nanstd(y_cnt(hb_cnt>min_hb)./hb_cnt(hb_cnt>min_hb));
                    
                    cov = hbdata{e}.sigma(:,1)./hbdata{e}.mu(:,1); % hb
                    hb_med_cov(e,s) = median(cov(hb_m>min_hb));
                    
                    cov = hbdata{e}.sigma(:,2)./hbdata{e}.mu(:,2); % y
                    y_med_cov(e,s) = median(cov(hb_m>min_hb));
                    
                end
            end          
            
    end % loop over embryos
    oldset = dataset; 
end % loop of slides

% get rid of unused, pre-allocated cells;
NucsT = NucsT(1:j);


          if plot_grads == 1
            figure(22); 
             ylabel('number of mRNA transcripts per cell');
             xlabel('distance (nm)');
            legend(NucsT{:});
            % title(texlabel(fname,'literal')); 
          end
            
            

figure(1); clf; boxplot(yhb_ave,'labels',slides);
ylabel('y:hb ratio'); set(gcf,'color','w');

% figure(2); clf; boxplot(yhb_std./yhb_ave,'labels',slides);
% ylabel('cov of y:hb ratio');  set(gcf,'color','w');


figure(1); clf; boxplot([[yhb_ave(:,1);yhb_ave(:,2)], [yhb_ave(:,3);yhb_ave(:,4)], [yhb_ave(:,5);yhb_ave(:,6)]],'labels',slides(1:2:end));
ylabel('y:hb ratio'); set(gcf,'color','w');



figure(3); clf; boxplot([y_med_cov,hb_med_cov]); 
ylabel('cov of mRNA counts'); set(gcf,'color','w');

hb_med_com = reshape(hb_med_cov,S*Tembs,1);
hb_med_com = hb_med_com(logical(1-isnan(hb_med_com)));
hb_med_com = [hb_med_com; NaN*ones(Tembs-length(hb_med_com),1)];

figure(3); clf; boxplot([y_med_cov,hb_med_com],'labels',[slides(:); 'wt']); 
ylabel('cov of mRNA counts'); ylim([0,.3]); set(gcf,'color','w');


figure(2); clf; boxplot([[y_med_cov(:,1);y_med_cov(:,2)],[y_med_cov(:,3);y_med_cov(:,4)],[y_med_cov(:,5);y_med_cov(:,6)],[hb_med_com;hb_med_com]],'labels',[slides(1:2:end); 'wt']); 
ylabel('median cov of mRNA counts'); ylim([0,.3]); set(gcf,'color','w');

figure(4); clf; boxplot(y_med_cov./hb_med_cov,'labels',slides);
ylabel('ratio of cov of mRNA counts');  set(gcf,'color','w');




    %%

    
pts = 500;
solid = .2*pts:.9*pts;
xmin = min(cat(1,pxdata{:}));
xmax = max(cat(1,pxdata{:}));
xs = linspace(xmin,xmax,pts)+abs(xmin);

N = sum(1-cellfun('isempty', pxdata(1,:)));
Es = zeros(1,N);

figure(21); clf;
   figure(22); clf;

mcurve = zeros(N,pts);
stdcurve = zeros(N,pts); 
fitcurve = zeros(N,length(solid));
fit_hb = zeros(N,length(solid));
mhb = zeros(N,pts);
shb = zeros(N,pts); 

    hbpool = {}; hbx = {};
    
 for s = 1:N   
    c = zeros(50,pts); 
    hb = zeros(50,pts); 
    Es(s) = sum(1-cellfun('isempty',pxdata(:,s)));
    for emb = 1:Es(s)
            x = pxdata{emb,s}+abs(xmin);
            y = pcurve{emb,2,s};
            c(emb,:) = interp1(x',y',xs);% linspace(min(x),max(x),100));
            figure(21); plot(x,y,'.','color',[s/N,0,1-s/N],'MarkerSize',3); hold on;
            
            if indiv_endog == 1  
                hb(emb,:) = interp1(x',pcurve{emb,1,s}',xs);
                plot(x,pcurve{emb,1,s},'.','color','k','MarkerSize',3); hold on;
            end
           % plot(xs,c(emb,:),'bo','MarkerSize',2); hold on;
    end

    
    c(c==0) = NaN;    
    mcurve(s,:) = nanmean(c); % yellow
    stdcurve(s,:) = nanstd(c);  % yellow variation
    y = smooth(mcurve(s,:),.3,'rloess');
    [jnk, fitcurve(s,:)] = fxn_fit_sigmoid(xs(solid),y(solid)',[5,mean(xs),max(mcurve(s,:)),0],'r',[NaN,NaN,NaN,0]);
    
    figure(21); 
    plot(xs,mcurve(s,:),'color',[s/N,0,1-s/N]);
    plot(xs(solid),fitcurve(s,:),'color',[s/N,0,1-s/N]); hold on;
    
    figure(22); 
    plot(xs(solid),fitcurve(s,:),'color',[s/N,0,1-s/N]); hold on;
     plot(xs,mcurve(s,:),'o','color',[s/N,0,1-s/N],'MarkerSize',3); hold on;
    errorbar(xs(2*s:12:end),mcurve(s,2*s:12:end),stdcurve(s,2*s:12:end),'linestyle','none','color',[s/N,0,1-s/N]);
    
    if indiv_endog == 1 ;%  && s~=3
        hb(hb==0) = NaN;
        mhb(s,:) = nanmean(hb);
        shb(s,:) = nanstd(hb); 
        h = smooth(mhb(s,:),.3,'rloess');
        [jnk2, fit_hb(s,:)] = fxn_fit_sigmoid(xs(solid),h(solid)',[4,mean(xs),max(mhb(s,:)),0],'r',[NaN,NaN,NaN,0]);

        figure(21); 
         plot(xs, mhb(s,:) ,'k'); 
         plot(xs(solid),fit_hb(s,:),'color','k'); hold on;
         
         figure(22); 
         plot(xs(solid),fit_hb(s,:),'k--'); hold on;
     plot(xs,mhb(s,:),'o','color','k','MarkerSize',3); hold on;
   % errorbar(xs(2*s:12:end),mhb(s,2*s:12:end),shb(s,2*s:12:end),'linestyle','none','color','k');
   
         
     else
        hbpool = [hbpool; pcurve(:,1,s)];
        hbx = [hbx;pxdata(:,s)]; 
    end
 end

 if indiv_endog ~=1
     has_data = logical(1-cellfun('isempty',hbx));
     hbpool2 = hbpool(has_data);
     hbx2 = hbx(has_data);
   
     hb = zeros(length(hbx2),pts); 
     for e = 1:length(hbx2)
            x = cell2mat(hbx2)'+abs(xmin);
            hb(e,:)=interp1(x,cell2mat(hbpool2)',xs);% linspace(min(x),max(x),100));
            figure(21); plot(x,cell2mat(hbpool2)','o','color','k','MarkerSize',1); hold on;
     end         
        hb(hb==0) = NaN;
        mhb(1,:) = nanmean(hb);
        shb(1,:) = nanstd(hb); 
        h = smooth(mhb(1,:),.3,'rloess');
        [jnk2, fit_hb(1,:)] = fxn_fit_sigmoid(xs(solid),h(solid)',[4,mean(xs),max(mhb(s,:)),0],'r',[NaN,NaN,NaN,0]);
        figure(21); 
         plot(xs, mhb(1,:) ,'k--'); 
         plot(xs(solid),fit_hb(1,:),'color','k'); hold on;
            
 end
 
 figure(21); 
 xlim([4E4,17E4]);
 ylim([0,0.012]); 
 xlabel('distance (nm)'); ylabel('mRNA density (molecules/50nm voxel)');
 legend(['no proximal N=', num2str(Es(1))],['no distal N=', num2str(Es(2))],...
     ['control N=', num2str(Es(3))],['wt N=', num2str(sum(Es))]);
 
  figure(22); 
 xlim([4E4,17E4]);
 ylim([0,0.012]); 
 xlabel('distance (nm)'); ylabel('mRNA density (molecules/50nm voxel)');
 legend(['no proximal N=', num2str(Es(1))],['no distal N=', num2str(Es(2))],...
     ['control N=', num2str(Es(3))],['wt N=', num2str(sum(Es))]);

%%

cmin = .00; cmax = .008;

 fname = 's02_MP01_Hz_22C';  ver = '_v2';  e = 6;
    load([folder,fname,ver,'_slidedata',vout],'hbdata'); 
    figure(2); clf; set(gcf,'color','k');
    
    subplot(1,3,1);
    imagesc(hbdata{e}.PlotmRNA2-min(hbdata{e}.Data_sort(:,3))  ); colormap hot;  
    colorbar;
    axis off; 
   caxis([0,.95*max(hbdata{e}.mu(:,2))]);
      
      
      fname = 's03_MP02_Hz_22C';  ver = '';e = 7;
      load([folder,fname,ver,'_slidedata',vout],'hbdata'); 

      subplot(1,3,2);
    imagesc(hbdata{e}.PlotmRNA2-min(hbdata{e}.Data_sort(:,3))  ); colormap hot;   colorbar;
    axis off; 
      caxis([0,.95*max(hbdata{e}.mu(:,2))]);
    
     fname = 's01_MP09_Hz_22C_b';  ver = '_v2'; e=5;
           load([folder,fname,ver,'_slidedata',vout],'hbdata'); 
subplot(1,3,3); 
    imagesc(hbdata{e}.PlotmRNA2-min(hbdata{e}.Data_sort(:,3))  ); colormap hot;  colorbar; 
    axis off; 
      caxis([0,.95*max(hbdata{e}.mu(:,2))]);
    


%%

if other_plots == 1;




   
    
figure(4); clf; colordef black; set(gcf,'color','k');
subplot(2,3,1);  
title('no proximal');
subplot(2,3,2);  imagesc(MP02{e2}.PlotmRNA2); colormap hot; colorbar; axis off; caxis([50,150]);
title('no distal');
subplot(2,3,3);  imagesc(MP09{e9}.PlotmRNA2); colormap hot; colorbar; axis off; caxis([50,300]);
title('control');
subplot(2,3,4);  imagesc(MP01{e1}.PlotmRNA); colormap hot; colorbar; axis off; caxis([100,275]);
subplot(2,3,5);  imagesc(MP02{e2}.PlotmRNA); colormap hot; colorbar; axis off; caxis([100,350]);
subplot(2,3,6);  imagesc(MP09{e9}.PlotmRNA); colormap hot; colorbar; axis off;caxis([100,350]);

    
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
title('mid cc13 mRNA counts');


%% 






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

%%
     fname = 's03_MP02_Hz_22C';  ver = '';
         
    load([folder,fname,'_graddata',ver],'hbdata'); 
    figure(1); clf; set(gcf,'color','k');
    imagesc(hbdata{3}.PlotmRNA2); colormap hot; colorbar; axis off; caxis([20,275]);
    
    Hb = makeuint(hbdata{3}.PlotmRNA2,8);
    mask = im2bw(Hb,.0001); 
    Kr = 1-im2bw(Hb,.55);
    
    HbKr = (double(Hb) + 400*Kr).*mask;
    
    figure(2); clf; set(gcf,'color','k');
    imagesc(HbKr);  colormap hot; colorbar; axis off; caxis([20,275]);
    
    %%
      figure(1); clf; set(gcf,'color','k');
      
       Hb = makeuint(hbdata{3}.PlotmRNA,8);
      
    imagesc(Hb); colormap hot; colorbar; axis off; caxis([20,275]);
    
   
    mask = im2bw(Hb,.0001); 
    Kr = 1-im2bw(Hb,.6);
    
    HbKr = (double(Hb) + 400*Kr).*mask;
    
    figure(2); clf; set(gcf,'color','k');
    imagesc(HbKr);  colormap hot; colorbar; axis off; caxis([20,275]);
    
end     

