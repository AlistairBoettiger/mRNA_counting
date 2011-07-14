
%%                          Plot_sna_data.m
% Alistair Boettiger                                Date Begun: 06/27/11
% Levine Lab         

%% Notes
%  this script is written to work on the ouptut of anlz_sna_counting data.
%  It takes combines multiple slides into a common plotting and statistical
%  analysis framework (previous script combines multiple embryos from the
%  same slide into a common dataset).  

clear all;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/';
 
% slides = {'y control'; 'y no distal'; 'y no proximal';'y control het'; 'y no proximal het';'cntrl y, sna[IIG]'};
% slides = {'cntrl y, sna[IIG]'};
 slides = {'primary alone','shadow removed', 'primary removed', 'wt'}; % { 'sna2.8Hz', 'MP08Hz','MP07Hz', 'wt'};
%slides = {'sna2.8Hz','wt'};   % slides = {'MP10Hz', 'MP06Hz'
%'MP10het','MP05het'};
ver = '';
vout= '';

% record stats
S = length(slides);
Ndatasets = 1; 
oldset = 0;
dset = 0; 

% set up figures for legend
figure(1); clf;  figure(3); clf;  figure(4); clf;
colordef white; set(gcf,'color','w');
 for s = 1:S
        figure(1); hold on; plot(0,500,'-','color',[s/S,0,1-s/S]); 
        colordef white; set(gcf,'color','w');
        figure(3); hold on; plot(0,500,'-','color',[s/S,0,1-s/S]);
        colordef white; set(gcf,'color','w');
        figure(4);hold on;  plot(0,500,'-','color',[s/S,0,1-s/S]);
        colordef white; set(gcf,'color','w');
 end
figure(1); legend(slides);
 figure(3); legend(slides);
 figure(4); legend(slides); 

ave_sna = cell(S,1);
std_sna = cell(S,1); 
cov_sna = cell(S,1);
ave_y = cell(S,1); 
std_y = cell(S,1); 
cov_y = cell(S,1);

Tembs  = 40; 


pcurve = cell(Tembs,2,Ndatasets);
pxdata = cell(Tembs,Ndatasets); 


for s = 1:S
     ave_sna{s} = NaN*zeros(Tembs,1); 
     std_sna{s} = NaN*zeros(Tembs,1);
     cov_sna{s} = NaN*zeros(Tembs,1);
     ave_y{s} = NaN*zeros(Tembs,1); 
     std_y{s} = NaN*zeros(Tembs,1); 
     cov_y{s} = NaN*zeros(Tembs,1);
    
     % Slide label information 
switch slides{s}
    case 'y control' % 'MP10Hz'
       date = '2011-05-22/';  
      fname = 's04_MP10Hz';% 
      chns = 2; ver = '_v2';
      skip = [1,3];   
      dataset = 1;
      
    case 'y no distal' % 'MP06Hz'
        date = '2011-05-22/'; 
      fname =  's05_MP06Hz' ;%  
      chns = 2;
      ver = ''; %'_v2';
      skip = 10;% 
      dataset = 2;
      
    case 'y no distal B' % 'MP06Hz'
        date = '2011-05-22/'; 
      fname =  's05_MP06Hz_b' ;%  
      chns = 2;
      ver = '_v4';
      skip = 20;% 
      dataset = 2;
 
    case 'y no proximal' % 'MP05Hz'
        date ='2011-06-20/';%  
        fname ='s07_MP05Hz_22C';
          ver = ''; % '_v2';
        chns = 2;
        skip = 7;
        dataset = 3;
      
      
    case 'y control het' % 'MP10het'
        date = '2011-02-17/'; 
        fname ='MP10_22C_sna_y_d'; 
        ver = '_v3';%
        chns = 2; 
       skip = 9;   
       dataset = 4;
      
    case 'y no proximal het' % 'MP05het'
        date = '2011-02-17/'; 
        fname = 'MP05_22C_sna_y_c'; 
        ver = '_v2';
        chns = 2;
        skip = 3;
        dataset = 5;
        
   
    case 'cntrl y, sna[IIG]'
      date = '2011-05-22/'; 
      fname =  's06_MP10_sna18_b'; %  
      chns = 2;
      ver = '_v4';
      skip = 20;% 
      dataset = 6;
        
        
        
        
        
   % % Rescues, used separately, restart dataset indexing from 1 
    case 'wt'
        date = '2011-02-17/'; 
        fname ='MP10_22C_sna_y_d'; 
        ver = '_v3';%
        chns = 1; 
       skip = 9;    
        dataset = 1;       
        
   case 'primary alone' % 'sna2.8Hz'
        date = '2011-06-20/';
        fname = 'sna2.8Hz_snaD_22C';
        ver = '';
        chns = 1;
        skip = 15; 
        dataset = 2;
    
   case 'shadow removed' %  'MP08Hz'
        date = '2011-05-22/'; 
        fname =  's07_MP08Hz_snaD_22C';
        ver = '';
        chns = 1; 
        skip = [2,3,7]; % 3 misoriented
        dataset = 3;
        
    case 'primary removed'% 'MP07Hz'
        date = '2011-06-20/'; 
        fname = 'MP07Hz_snaD_22C';%
         ver = '';
        chns =1;
        skip =[1,3,5];
        dataset = 4;
        
    case '1x primary removed' % 'MP07het'
        date = '2011-05-22/'; 
        fname =  'MP07het_snaD_22C';%
         ver = '';
        chns =1;
        skip =15;
        dataset = 5;
  
end

% Actually grab the data from the chosen slide ('case').  
   load([folder,fname,ver,'_slidedata',vout,'.mat'],'data'); 
   % load([folder,date,fname,'_graddata',ver,'.mat'],'data');  
    Es = length(data); 
    offsets = zeros(1,Es); 
    
    % Show what were working on
    disp(fname);
    disp(['# embryos in dataset = ',num2str(sum(1-cellfun('isempty',data)))])
    try
        disp(['chn1 thresh: ', num2str(data{e}.ipars{1}.min_int)])
        disp(['chn2 thresh: ', num2str(data{e}.ipars{2}.min_int)])
        disp(['missG: ', num2str(data{e}.ipars{3})])
    catch er
        disp('processing parameter data not available'); 
    end


 % Get boundary point of first curve (from embryo 1 on slide 1)
 % Boundary alignment is done by fitting Hill functions and aligning their
 % threshold values.  
 if s == 1;
     e=1; % use embryo 1 as baseline
 % shorthand relabeling of data
  gx = data{e}.Data_sort(:,1);  
  grad = data{e}.Data_sort(:,2);  
  L = length(grad); st = floor(L/5);
  grad = grad(st:end); 
  gx = gx(st:end);
 
  % First guess parameters to fit curve for aligment
  n = 4; 
  theta = mean(gx);  
  A = max(grad);  
  b=min(grad); 
  [p,fit] = fxn_fit_sigmoid(gx',grad',[n,theta,A,b],'r');
  p1 = p(2); 
   offsets(e) = p1; 
   
    figure(20); clf; plot(gx', grad,'b');
    hold on; 
    plot(gx,fit,'r');
    plot(p(2),A/4,'r*','MarkerSize',20);
 
    k = 2;  % for slide 1 we'll start indexing curves to shift from k = 2;
 else
     k = 1;    % for all other slides we'll start indexing from 1
 end
  
 % fit curves to all other data to get different thresholds needed for
 % curve alignment.  
for e= k:Es;
    try
     gx = data{e}.Data_sort(:,1);
     grad = data{e}.Data_sort(:,2);  
     L = length(grad);st = floor(L/5);
     grad = grad(st:end); 
     gx = gx(st:end);

    figure(20); clf; plot(gx,grad);
    n = 3;
    theta = 2*mean(gx);  
    A = max(grad); 
    b=min(grad); 
    [p,fit] = fxn_fit_sigmoid(gx',grad',[n,theta,A,b],'r');

    offsets(e) = abs(p(2)); 

    figure(20); clf; plot(gx', grad,'b');
    hold on; 
    plot(gx,fit,'r');
    plot(p(2),A/4,'r*','MarkerSize',20);
    catch er
        disp(er.message); 
        continue
    end
end


figure(1);
colordef white; set(gcf,'color','w');



% restart embryo counter if we've switched datasets  
    if dataset - oldset ~= 0; 
        emb = 0;
        disp('starting new dataset'); 
        dset = dset + 1; 
    end


for e = 1:Es 
        if isempty(find(e==skip,1))  && isempty(data{e}) ~= 1;
            emb =emb+1;
            
             
            
            mcorr =   1;%   data{e}.Nnucs./data{1}.Nnucs; 
            sna_cnt = mcorr*data{e}.Data_sort(:,2) - min(mcorr*data{e}.Data_sort(:,2));
            sna_m = mcorr*data{e}.mu(:,1)  - min(mcorr*data{e}.Data_sort(:,2));
                      
            xdata = data{e}.Data_sort(:,1) + p1-offsets(e)  ; 
            
            % store whole curves
            pxdata{emb,dset} = xdata;  % dataset
            pcurve{emb,1,dset} = sna_cnt; %  dataset
            
            
            % compute median coefficient of variation
            med_cov  = data{e}.sigma(:,1);
            cov_sna{s}(e) = median(med_cov(sna_m>max(sna_m)/2)./sna_m(sna_m>max(sna_m)/2)); 
            
     %     cnt_dist{s}{e} = sna_cnt(sna_cnt>max(sna_m)/2);
            ave_sna{s}(e) = mean(sna_cnt(sna_cnt>max(sna_m)/2)) ;
            std_sna{s}(e) = std(sna_cnt(sna_cnt>max(sna_m)/2)) ;
            
            if  chns == 1
                figure(1);
                plot(xdata,sna_cnt ,'o','color',[s/S,0,1-s/S],'MarkerSize',2); % check results    [e/Es,0,1-e/Es]
                hold on; 
                plot(data{e}.x+ p1-offsets(e), sna_m,'.-','color',[s/S,0,1-s/S],'MarkerSize',5); % [e/Es,0,1-e/Es]
                ylim([0,250]); 
                
                sna_hist  =  hist(sna_cnt(sna_cnt>max(sna_m)/2),0:10:250); 
                figure(4); hold on;
                plot(0:10:250,sna_hist,'color',[s/S,0,1-s/S],'linewidth',2);
                 ylim([0,200]); xlim([0,250]); 
                 xlabel('mRNA counts');
                 ylabel('frequency (# cells)');
            end
          

            
      %  errorbar(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,1) - min(mcorr*data{e}.mu(:,1)) ,mcorr*data{e}.sigma(:,1),'linestyle','none','linewidth',3,'color',[e/Es,0,1-e/Es],'MarkerSize',1);
           if chns == 2;
                figure(1);
                 plot(data{e}.Data_sort(:,1) + p1-offsets(e),sna_cnt ,'.','color','k','MarkerSize',1); % check results    [e/Es,0,1-e/Es]
                 hold on; 
                 plot(data{e}.x+ p1-offsets(e), sna_m,'.','color','k','MarkerSize',5); % [e/Es,0,1-e/Es]
               
                y_cnt = mcorr*data{e}.Data_sort(:,3)- min(mcorr*data{e}.Data_sort(:,3));
                y_m = mcorr*data{e}.mu(:,2)- min(mcorr*data{e}.Data_sort(:,3));
                
                % whole curve for later plotting
                pcurve{emb,2,dset} = y_cnt;  % dataset
                
                figure(1); 
                plot(data{e}.Data_sort(:,1)+  p1-offsets(e),y_cnt,'.','color',[s/S,0,1-s/S],'MarkerSize',1); %[(s-1)/2,1-e/Es,e/Es]
                hold on; 
                plot(data{e}.x+ p1-offsets(e),y_m,'o','color',[s/S,0,1-s/S],'MarkerSize',5); % [0,1-e/Es,e/Es]

                ave_y{s}(e) = mean(y_cnt(sna_cnt>max(sna_m)/2)) ;
                std_y{s}(e) = std(y_cnt(sna_cnt>max(sna_m)/2)) ;
                
                med_cov  = data{e}.sigma(:,2);
                cov_y{s}(e) = median(med_cov(sna_m>max(sna_m)/2)./y_m(sna_m>max(sna_m)/2)); 

                figure(3); hold on; 
            % plot(data{e}.Data_sort(:,1)+  p1-offsets(e), data{e}.Data_sort(:,3)./data{e}.Data_sort(:,2) ,'.','color',[e/Es,0,1-e/Es],'MarkerSize',5); ylim([0,2]); 
            % plot(mean(data{e}.Data_sort(:,3)./data{e}.Data_sort(:,2)),'.','color',[e/Es,0,1-e/Es]);

              errorbar(data{e}.x+  p1-offsets(e),data{e}.sigma(:,2)./data{e}.mu(:,2),data{e}.bssigma(:,2)./data{e}.mu(:,2),'o','color',[s/S,0,1-s/S],'MarkerSize',5);
               ylabel('CoV of mRNA count'); xlabel('distance (nm)');
              ylim([0,.50]);
            
           if e == 2;% 2
            y_hist  =  hist(y_cnt(sna_cnt>max(sna_m)/2),0:5:250); 
            figure(4); hold on; plot(0:5:250,y_hist,'color',[s/S,0,1-s/S],'linewidth',4);
           end
           
              % errorbar(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,2)-min(mcorr*data{e}.mu(:,2)),mcorr*data{e}.sigma(:,2),'linestyle','none','linewidth',3,'color',[0,1-e/Es,e/Es],'MarkerSize',1);          
          end
        end
end

    
end % end loop over slides

figure(1);  ylabel('number of mRNA transcripts per cell'); xlabel('distance (nm)');


%%

if chns == 1;
    
    sna_lev = zeros(Tembs,S);
    sna_std = zeros(Tembs,S);
    sna_cov = zeros(Tembs,S); 
    
    for s = 1:S
        sna_lev(:,s) = ave_sna{s};
        sna_std(:,s) = std_sna{s};
        sna_cov(:,s) = cov_sna{s}; 
    end

    
    figure(2); clf; boxplot(sna_lev,'width',.8,'labels',slides);
    ylabel('mRNA counts'); ylim([0,225]);
     figure(3); clf; boxplot(sna_cov,'width',.8,'labels',slides);
    ylabel('CoV for mRNA counts'); ylim([0,0.25]);
end

if chns == 2
    wt_sna = [];

    ybox = zeros(Tembs,2*S);
    yhb = zeros(Tembs,S);
    yvar = zeros(Tembs,2*S); 
    ycov = zeros(Tembs,2*S);
    for s = 1:S
        ybox(:,2*s) = ave_y{s};
        ybox(:,2*s-1) = ave_sna{s};

        yhb(:,s) = ave_y{s}./ave_sna{s};
        
        yvar(:,2*s) = std_y{s}./ave_y{s};
        yvar(:,2*s-1) = std_sna{s}./ave_sna{s};  
        
        ycov(:,2*s) = cov_y{s};
        ycov(:,2*s-1) = cov_sna{s};

        wt_sna = [wt_sna; ave_sna{s}];
    end


    figure(2); clf; boxplot(ybox,'colors',[1,0,0;0,1,0],'width',.8);
        ylim([0,250]);
        set(gcf,'color','w'); ylabel('mRNA counts');

    wt_sna = wt_sna(logical(1-isnan(wt_sna)));
    wt_sna = [wt_sna; NaN*zeros(Tembs-length(wt_sna),1)];

    labs = ['wt'; slides(:)];
    
%     figure(2); clf; boxplot( [wt_sna,ybox(:,2:2:end)],'whisker',0,...
%     'labels',labs);
%      ylim([0,250]);
%  
%      figure(5); clf; boxplot(ycov,'colors',[1,0,0;0,1,0],'width',.8,...
%         'labels',{'wt','2x y-cntrl','wt','2x no-shadow','wt','2x no-prox','wt','1x y-cntrl','wt','1x no-prox'});
%         ylim([0,0.3]);
%         set(gcf,'color','w'); ylabel('CoV for mRNA counts');
    
        
   figure(1); clf; boxplot(ybox); 
        
   figure(2); boxplot(yhb,'width',.8,'labels',slides);
     ylabel('ratio, y / hb');
     
    
end
 

%% Playing with different plots of saved curve data

indiv_endog = 1; chns = 1;



pts = 500;
solid = .3*pts:.9*pts;
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
fit_dat = zeros(N,length(solid));
mdat = zeros(N,pts);
sdat = zeros(N,pts); 

    datpool = {}; 
    datx = {};
    
 for s = 1:N   % s =1
    c = zeros(50,pts); 
    dat = zeros(50,pts); 
    Es(s) = sum(1-cellfun('isempty',pxdata(:,s)));
    for emb = 1:Es(s)
        
        if chns == 2 
            x = pxdata{emb,s}+abs(xmin);
            y = pcurve{emb,2,s};
            [x,v] = unique(x);
            y = y(v);
            
            c(emb,:) = interp1(x',y',xs);% linspace(min(x),max(x),100));
            figure(21); plot(x,y,'.','color',[s/N,0,1-s/N],'MarkerSize',1); hold on;
        end
        
            if indiv_endog == 1  
                x = pxdata{emb,s}+abs(xmin);
                [x,v] = unique(x);
                y = pcurve{emb,1,s}';
                y = y(v);
                dat(emb,:) = interp1(x',y,xs);
                if chns == 2
                    plot(x,y,'.','color','k','MarkerSize',1); hold on;
                else
                    plot(x,y,'.','color',[s/N,0,1-s/N],'MarkerSize',1); hold on;
                end
            end
           % plot(xs,c(emb,:),'bo','MarkerSize',2); hold on;
    end

    if chns ==2;
        c(c==0) = NaN;    
        mcurve(s,:) = nanmean(c); % yellow
        stdcurve(s,:) = nanstd(c);  % yellow variation

        y = smooth(mcurve(s,:),.1,'rloess');
        [jnk, fitcurve(s,:)] = fxn_fit_sigmoid(xs(solid),y(solid)',[5,mean(xs),max(mcurve(s,:)),0],'r',[NaN,NaN,NaN,0]);

        figure(21); 
            plot(xs,mcurve(s,:),'color',[s/N,0,1-s/N]);
            plot(xs(solid),fitcurve(s,:),'color',[s/N,0,1-s/N],'linewidth',2); hold on;   
        figure(22); 
            plot(xs(solid),fitcurve(s,:),'color',[s/N,0,1-s/N],'linewidth',2); hold on;
             plot(xs,mcurve(s,:),'o','color',[s/N,0,1-s/N],'MarkerSize',3); hold on;
            errorbar(xs(2*s:12:end),mcurve(s,2*s:12:end),stdcurve(s,2*s:12:end),'linestyle','none','color',[s/N,0,1-s/N]);
    end
    
    if indiv_endog == 1 ;%  && s=3
        dat(dat==0) = NaN;
 
        dat_peak = nanmedian(dat,2);

        low = nanmean(  dat(dat_peak<60,:)  );
        mid = nanmean( dat(dat_peak>60,:)   );
        
        

       %  figure(1); clf; plot(xs,low); hold on; plot(xs,mid);
        
        mdat(s,:) = nanmean(dat);
        sdat(s,:) = nanstd(dat); 
        h = smooth(mdat(s,:),.1,'rloess');
        [jnk2, fit_dat(s,:)] = fxn_fit_sigmoid(xs(solid),h(solid)',[4,mean(xs),max(mdat(s,:)),0],'r',[NaN,NaN,NaN,0]);

        if chns == 2;
             figure(21); 
                 plot(xs, mdat(s,:) ,'k'); 
                 plot(xs(solid),fit_dat(s,:),'color','k','linewidth',2); hold on; 
             figure(22); 
                 plot(xs(solid),fit_dat(s,:),'k--','linewidth',2); hold on;
                plot(xs,mdat(s,:),'o','color','k','MarkerSize',3); hold on;
                errorbar(xs(7:12:end),mdat(s,7:12:end),sdat(s,7:12:end),'linestyle','none','color','k');
        elseif chns ==1;
             figure(21); 
                 plot(xs, mdat(s,:),'color' ,[s/N,0,1-s/N]); 
                 plot(xs(solid),fit_dat(s,:),'color',[s/N,0,1-s/N],'linewidth',2); hold on; 
             figure(22); 
                 plot(xs(solid),fit_dat(s,:),'--','color',[s/N,0,1-s/N],'linewidth',2); hold on;
                plot(xs,mdat(s,:),'o','color',[s/N,0,1-s/N],'MarkerSize',3); hold on;
                errorbar(xs(7:12:end),mdat(s,7:12:end),sdat(s,7:12:end),'linestyle','none','color',[s/N,0,1-s/N]);
        end
     else
        datpool = [datpool; pcurve(:,1,s)];
        datx = [datx;pxdata(:,s)]; 
    end
 end
 
 if indiv_endog ~=1
     has_data = logical(1-cellfun('isempty',datx));
     datpool2 = datpool(has_data);
     datx2 = datx(has_data);
   
     dat = zeros(length(datx2),pts); 
     for e = 1:length(datx2)
            x = cell2mat(datx2)'+abs(xmin);
            y = cell2mat(datpool2)';
            [x,v] = unique(x);
            y = y(v);     
            
            dat(e,:)=interp1(x,y,xs);% linspace(min(x),max(x),100));
            if chns == 2
                 figure(21); plot(x,y,'o','color','k','MarkerSize',1); hold on;
            else
                figure(21); plot(x,y,'o','color',[e/N,0,1-e/N],'MarkerSize',1); hold on;
            end
     end         
        dat(dat==0) = NaN;
        mdat(1,:) = nanmean(dat);
        sdat(1,:) = nanstd(dat); 
        h = smooth(mdat(1,:),.1,'rloess');
        [jnk2, fit_dat(1,:)] = fxn_fit_sigmoid(xs(solid),h(solid)',[5,mean(xs),max(h(solid)),0],'r',[NaN,NaN,NaN,0]);
        figure(21); 
             % plot(xs, mdat(1,:) ,'k--'); figure(1); clf;
             plot(xs, h ,'k--');
             plot(xs(solid),fit_dat(1,:),'color','k','linewidth',2); hold on;      
        figure(22); 
            plot(xs(solid),fit_dat(1,:),'k--'); hold on;
            plot(xs,mdat(1,:),'o','color','k','MarkerSize',3); hold on;
            errorbar(xs(1:12:end),mdat(1,1:12:end),sdat(1,1:12:end),'linestyle','none','color','k');     
 end
 
 figure(21); 
 ylim([0,250]);
 xlim([.1E5,1.6E5]); 
 xlabel('distance (nm)'); ylabel('mRNA density (molecules/50nm voxel)');
  if chns == 2; % label for yellow analysis only  
      legend(['no proximal N=', num2str(Es(1))],['no distal N=', num2str(Es(2))],...
    ['control N=', num2str(Es(3))],['wt N=', num2str(sum(Es))]);
  end

  figure(22); 
 ylim([0,250]);
 xlim([.1E5,1.6E5]); 
 xlabel('distance (nm)'); ylabel('mRNA density (molecules/50nm voxel)');

   if chns ==2; % label for yellow analysis only    
       legend(['no proximal N=', num2str(Es(1))],['no distal N=', num2str(Es(2))],...
      ['control N=', num2str(Es(3))],['wt N=', num2str(sum(Es))]);
   end




%%

folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/';
date = '2011-05-22/';
fname = 's05_MP06Hz_b'; ver = '_v4';
load([folder,date,fname,'_graddata',ver],'data');  


e = 10;

 figure(8); clf; colordef white; set(gcf,'color','w'); 
        plot(data{e}.Data_sort(:,1)  ,data{e}.Data_sort(:,2) - min(data{e}.Data_sort(:,2)) ,'.','color','k','MarkerSize',5); % check results    [e/Es,0,1-e/Es]
        hold on; 
        errorbar(data{e}.x,data{e}.mu(:,1)  - min(data{e}.Data_sort(:,2)), data{e}.sigma(:,1),'linestyle','none','linewidth',3,'color','k','MarkerSize',5); % [e/Es,0,1-e/Es]
        
         plot(data{e}.Data_sort(:,1),data{e}.Data_sort(:,3)- min(data{e}.Data_sort(:,3)),'.','color','g','MarkerSize',5); %[(s-1)/2,1-e/Es,e/Es]
            hold on; 
         errorbar(data{e}.x,data{e}.mu(:,2)- min(data{e}.Data_sort(:,3)), data{e}.sigma(:,2),'linestyle','none','linewidth',3,'color','g','MarkerSize',5); % [0,1-e/Es,e/Es]

           ylabel('number of mRNA transcripts per cell','Fontsize',14); xlabel('distance (nm)','Fontsize',14);
         set(gcf,'color','w'); set(gca,'Fontsize',14);
         
         
%%  Loss of snail in snaN


%   fin = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Snail Patterning/Results/';
% 
%   I = imread([fin,'snaN_sna_sim.tif']);
%   I = imread([fin,'snaN_snail_sog.tif']);
%   
%   Iout = imadjust(I(:,:,1),[0,1]);
%   figure(1); clf; imagesc(Iout); colormap hot;
%  axis off; set(gcf,'color','k'); caxis([20,350]);
%  
%  I2 = imread(['/Volumes/Data/Lab Data/Shadow_data/MP10_22C_y_sna_09.tif']);
%   I2 = imflip(I2,2);
%   
%   
%   Iout = imadjust(I2(:,:,2),[.08,.7]);
%   figure(1); clf; imagesc(Iout); colormap hot;
%  axis off; set(gcf,'color','k');
%  
 