
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
 
 slides = {'MP10Hz', 'MP06Hz' 'MP10het','MP05het'};
% slides = {'primary alone','shadow removed', 'primary removed', 'wt'}; % { 'sna2.8Hz', 'MP08Hz','MP07Hz', 'wt'};
%slides = {'sna2.8Hz','wt'};
ver = '';

% record stats
S = length(slides);


% set up figures for legend
figure(1); clf;  figure(3); clf;  figure(4); clf;
 for s = 1:S
        figure(1); hold on; plot(0,500,'-','color',[s/S,0,1-s/S]);
        figure(3); hold on; plot(0,500,'-','color',[s/S,0,1-s/S]);
        figure(4);hold on;  plot(0,500,'-','color',[s/S,0,1-s/S]);
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

for s = 1:S 
     ave_sna{s} = NaN*zeros(Tembs,1); 
     std_sna{s} = NaN*zeros(Tembs,1);
     cov_sna{s} = NaN*zeros(Tembs,1);
     ave_y{s} = NaN*zeros(Tembs,1); 
     std_y{s} = NaN*zeros(Tembs,1); 
     cov_y{s} = NaN*zeros(Tembs,1);
    
     % Slide label information 
switch slides{s}
case 'MP10Hz'
       date = '2011-05-22/';  
      fname = 's04_MP10Hz';% 's05_MP06Hz' ;%  's07_MP08Hz_snaD_22C';
      chns = 2; ver = '_v3';
      skip = [1,3];   
      
    case 'MP06Hz'
        date = '2011-05-22/'; 
      fname =  's05_MP06Hz' ;%  's07_MP08Hz_snaD_22C';
      chns = 2;
      ver = '_v2';
      skip = 10;% 
               
    case 'MP10het'
        date = '2011-02-17/'; 
        fname ='MP10_22C_sna_y_d'; 
        ver = '_v3';%
        chns = 2; 
       skip = 9;    
      
    case 'MP05het'
        date = '2011-02-17/'; 
        fname = 'MP05_22C_sna_y_c'; 
        ver = '_v2';
        chns = 2;
        skip = 3;
        
        
        
   
    case 'wt'
        date = '2011-02-17/'; 
        fname ='MP10_22C_sna_y_d'; 
        ver = '_v3';%
        chns = 1; 
       skip = 9;    
        
%        date = '2011-05-22/';  
%       fname = 's04_MP10Hz';% 's05_MP06Hz' ;%  's07_MP08Hz_snaD_22C';
%       chns = 1; ver = '_v3';
%       skip = [1,3];   
        
        
   case 'primary alone' % 'sna2.8Hz'
        date = '2011-06-20/';
        fname = 'sna2.8Hz_snaD_22C';
        ver = '';
        chns = 1;
        skip = 15; 
    
   case 'shadow removed' %  'MP08Hz'
        date = '2011-05-22/'; 
        fname =  's07_MP08Hz_snaD_22C';
        ver = '';
        chns = 1; 
        skip = [2,3,7]; % 3 misoriented
        
    case 'primary removed'% 'MP07Hz'
        date = '2011-06-20/'; 
        fname = 'MP07Hz_snaD_22C';%
         ver = '';
        chns =1;
        skip =[1,3,5];
        
    case '1x primary removed' % 'MP07het'
        date = '2011-05-22/'; 
        fname =  'MP07het_snaD_22C';%
         ver = '';
        chns =1;
        skip =15;
        
  
end

% Actually grab the data from the chosen slide ('case').  
    load([folder,date,fname,'_graddata',ver,'.mat'],'data');  
    Es = length(data); 
    offsets = zeros(1,Es); 
    
    % Show what were working on
    disp(fname);
    disp(['# embryos in dataset = ',num2str(sum(1-cellfun('isempty',data)))])



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

for e = 1:Es 
        if isempty(find(e==skip,1))  && isempty(data{e}) ~= 1;
            
            mcorr =   1;%   data{e}.Nnucs./data{1}.Nnucs; 
            sna_cnt = mcorr*data{e}.Data_sort(:,2) - min(mcorr*data{e}.Data_sort(:,2));
            sna_m = mcorr*data{e}.mu(:,1)  - min(mcorr*data{e}.Data_sort(:,2));
            
            % compute median coefficient of variation
            med_cov  = data{e}.sigma(:,1);
            cov_sna{s}(e) = median(med_cov(sna_m>max(sna_m)/2)./sna_m(sna_m>max(sna_m)/2)); 
            
     %     cnt_dist{s}{e} = sna_cnt(sna_cnt>max(sna_m)/2);
            ave_sna{s}(e) = mean(sna_cnt(sna_cnt>max(sna_m)/2)) ;
            std_sna{s}(e) = std(sna_cnt(sna_cnt>max(sna_m)/2)) ;
            
            if  chns == 1
                figure(1);
                plot(data{e}.Data_sort(:,1) + p1-offsets(e),sna_cnt ,'o','color',[s/S,0,1-s/S],'MarkerSize',2); % check results    [e/Es,0,1-e/Es]
                hold on; 
                plot(data{e}.x+ p1-offsets(e), sna_m,'.-','color',[s/S,0,1-s/S],'MarkerSize',5); % [e/Es,0,1-e/Es]
                ylim([0,250]); 
                
                sna_hist  =  hist(sna_cnt(sna_cnt>max(sna_m)/2),0:10:250); 
                figure(4); hold on;
                plot(0:10:250,sna_hist,'color',[s/S,0,1-s/S],'linewidth',2);
                 ylim([0,200]); xlim([0,250]); 
                 xlabel('mRNA counts'); ylabel('frequency (# cells)');
            end
          



            
      %  errorbar(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,1) - min(mcorr*data{e}.mu(:,1)) ,mcorr*data{e}.sigma(:,1),'linestyle','none','linewidth',3,'color',[e/Es,0,1-e/Es],'MarkerSize',1);
           if chns == 2;
                figure(1);
                 plot(data{e}.Data_sort(:,1) + p1-offsets(e),sna_cnt ,'.','color','k','MarkerSize',1); % check results    [e/Es,0,1-e/Es]
                 hold on; 
                 plot(data{e}.x+ p1-offsets(e), sna_m,'.','color','k','MarkerSize',5); % [e/Es,0,1-e/Es]
               
                y_cnt = mcorr*data{e}.Data_sort(:,3)- min(mcorr*data{e}.Data_sort(:,3));
                y_m = mcorr*data{e}.mu(:,2)- min(mcorr*data{e}.Data_sort(:,3));
                
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

    D = zeros(Tembs,S);
    ybox = zeros(Tembs,2*S);
    yvar = zeros(Tembs,2*S); 
    ycov = zeros(Tembs,2*S);
    for s = 1:S
        D(1:Tembs,s) = ave_y{s}./ave_sna{s};
       % ave_y{s}(ave_y{s}<1) = NaN; 
       % ave_sna{s}(ave_sna{s}<1) = NaN; 
        ybox(:,2*s) = ave_y{s};
        ybox(:,2*s-1) = ave_sna{s};

       % std_y{s}(std_y{s}<1) = NaN; 
       % std_sna{s}(std_sna{s}<1) = NaN; 
        yvar(:,2*s) = std_y{s}./ave_y{s};
        yvar(:,2*s-1) = std_sna{s}./ave_sna{s};  
        
        ycov(:,2*s) = cov_y{s};
        ycov(:,2*s-1) = cov_sna{s};

        wt_sna = [wt_sna; ave_sna{s}];
    end
    figure(2); clf; boxplot(D); ylim([0,1]);

    figure(2); clf; boxplot(ybox,'colors',[1,0,0;0,1,0],'width',.8,...
        'labels',{'wt','2x y-cntrl','wt','2x no-shadow','wt','1x y-cntrl','wt','1x no-prox'});
        ylim([0,250]);
        set(gcf,'color','w'); ylabel('mRNA counts');

    figure(6); clf; boxplot(yvar,'colors',[0,0,1;0,1,0],'width',.8,'whisker',0,...
        'labels',{'wt','2x y-cntrl','wt','2x no-shadow','wt','1x y-cntrl','wt','1x no-prox'});
        set(gcf,'color','w'); ylabel('CoV for mRNA counts'); 
        ylim([0,.3]);

    wt_sna = wt_sna(logical(1-isnan(wt_sna)));
    wt_sna = [wt_sna; NaN*zeros(Tembs-length(wt_sna),1)];

    figure(2); clf; boxplot( [wt_sna,1.15*ybox(:,2:2:end)],'whisker',0,...
    'labels',{'wt','2x y-cntrl','2x no-shadow','1x y-cntrl','1x no-prox'});
     ylim([0,200]);
 
     figure(5); clf; boxplot(ycov,'colors',[1,0,0;0,1,0],'width',.8,...
        'labels',{'wt','2x y-cntrl','wt','2x no-shadow','wt','1x y-cntrl','wt','1x no-prox'});
        ylim([0,0.3]);
        set(gcf,'color','w'); ylabel('CoV for mRNA counts');
     
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
