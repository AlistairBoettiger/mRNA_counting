
%%                          Plot_sna_data.m
% Alistair Boettiger                                Date Begun: 06/27/11
% Levine Lab         

%% Notes
%  this script is written to work on the ouptut of anlz_sna_counting data.
%  It takes combines multiple slides into a common plotting and statistical
%  analysis framework (previous script combines multiple embryos from the
%  same slide into a common dataset).  

clear all;
 figure(1); clf;  figure(3); clf;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/';
 
slides = {'MP10Hz', 'MP06Hz' 'MP10het','MP05het'};
%slides = {'MP10Hz', 'MP08Hz','MP07Hz', 'MP07het'};
ver = '';

% record stats
S = length(slides);
ave_sna = cell(S,1);
std_sna = cell(S,1); 
ave_y = cell(S,1); 
std_y = cell(S,1); 

for s = 1:S 
     ave_sna{s} = zeros(20,1); 
     std_sna{s} = zeros(20,1); 
     ave_y{s} = zeros(20,1); 
     std_y{s} = zeros(20,1); 
    
     % Slide label information 
switch slides{s}
case 'MP10Hz'
       date = '2011-05-22/';  
      fname = 's04_MP10Hz';% 's05_MP06Hz' ;%  's07_MP08Hz_snaD_22C';
      chns = 2; ver = '_v3';
      skip = 3;   
      
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
        
  
end

% Actually grab the data from the chosen slide ('case').  
    load([folder,date,fname,'_graddata',ver],'data');  
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
            
            figure(1);
            plot(data{e}.Data_sort(:,1) + p1-offsets(e),sna_cnt ,'.','color','k','MarkerSize',1); % check results    [e/Es,0,1-e/Es]
            hold on; 
            plot(data{e}.x+ p1-offsets(e), sna_m,'.','color','k','MarkerSize',5); % [e/Es,0,1-e/Es]

            ave_sna{s}(e) = mean(sna_cnt(sna_cnt>max(sna_m)/2)) ;
            std_sna{s}(e) = std(sna_cnt(sna_cnt>max(sna_m)/2)) ;

        
      %  errorbar(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,1) - min(mcorr*data{e}.mu(:,1)) ,mcorr*data{e}.sigma(:,1),'linestyle','none','linewidth',3,'color',[e/Es,0,1-e/Es],'MarkerSize',1);
           if chns == 2;
                y_cnt = mcorr*data{e}.Data_sort(:,3)- min(mcorr*data{e}.Data_sort(:,3));
                y_m = mcorr*data{e}.mu(:,2)- min(mcorr*data{e}.Data_sort(:,3));
                
                figure(1); 
                plot(data{e}.Data_sort(:,1)+  p1-offsets(e),y_cnt,'.','color',[s/4,0,1-s/4],'MarkerSize',1); %[(s-1)/2,1-e/Es,e/Es]
                hold on; 
                plot(data{e}.x+ p1-offsets(e),y_m,'o','color',[s/4,0,1-s/4],'MarkerSize',5); % [0,1-e/Es,e/Es]

                ave_y{s}(e) = mean(y_cnt(y_cnt>max(y_m)/2)) ;
                std_y{s}(e) = std(y_cnt(y_cnt>max(y_m)/2)) ;

                figure(3); hold on; 
            % plot(data{e}.Data_sort(:,1)+  p1-offsets(e), data{e}.Data_sort(:,3)./data{e}.Data_sort(:,2) ,'.','color',[e/Es,0,1-e/Es],'MarkerSize',5); ylim([0,2]); 
            % plot(mean(data{e}.Data_sort(:,3)./data{e}.Data_sort(:,2)),'.','color',[e/Es,0,1-e/Es]);

              errorbar(data{e}.x+  p1-offsets(e),data{e}.sigma(:,2)./data{e}.mu(:,2),data{e}.bssigma(:,2)./data{e}.mu(:,2),'o','color',[s/4,0,1-s/4],'MarkerSize',5);
                % errorbar(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,2)-min(mcorr*data{e}.mu(:,2)),mcorr*data{e}.sigma(:,2),'linestyle','none','linewidth',3,'color',[0,1-e/Es,e/Es],'MarkerSize',1);          
          end
        end
end
        ylabel('number of mRNA transcripts per cell'); xlabel('distance (nm)');
        ylim([0,250]);

    
end % end loop over slides

%%
D = zeros(20,S);
ybox = zeros(20,2*S);
yvar = zeros(20,2*S); 
for s = 1:S
    D(1:20,s) = ave_y{s}./ave_sna{s};
    ave_y{s}(ave_y{s}<1) = NaN; 
    ave_sna{s}(ave_sna{s}<1) = NaN; 
    ybox(1:20,2*s) = ave_y{s};
    ybox(1:20,2*s-1) = ave_sna{s};
    
    std_y{s}(std_y{s}<1) = NaN; 
    std_sna{s}(std_sna{s}<1) = NaN; 
    yvar(1:20,2*s) = std_y{s}./ave_y{s};
    yvar(1:20,2*s-1) = std_sna{s}./ave_sna{s};  
end
figure(2); clf; boxplot(D); ylim([0,1]);

figure(2); clf; boxplot(ybox,'colors',[1,0,0;0,1,0],'width',.8,...
    'labels',{'wt','2x y-cntrl','wt','2x no-shadow','wt','1x y-cntrl','wt','1x no-prox'});
ylim([0,250]);
set(gcf,'color','w'); ylabel('mRNA counts');

figure(3); clf; boxplot(yvar,'colors',[1,0,0;0,1,0],'width',.8,...
    'labels',{'wt','2x y-cntrl','wt','2x no-shadow','wt','1x y-cntrl','wt','1x no-prox'});
set(gcf,'color','w'); ylabel('CoV for mRNA counts'); 
ylim([0,.25]);


%% For rescues... (not developed yet -- 7/5/11).  
 e=1;
switch 'MP07het'; % 'MP07Hz'; % 'MP06Hz'; % 'MP10het' % e =2;  % 'MP05het'; % e =7;

        
   case 'MP08Hz'
        date = '2011-05-22/'; 
        fname =  's07_MP08Hz_snaD_22C';
        ver = '';
        chns = 1; 
        skip = 15;
        
    case 'MP07Hz'
        date = '2011-06-20/'; 
        fname = 'MP07Hz_snaD_22C';%
         ver = '';
        chns =1;
        skip =15;
        
    case 'MP07het'
        date = '2011-05-22/'; 
        fname =  'MP07het_snaD_22C';%
         ver = '';
        chns =1;
        skip =15;
end

disp(fname); 

    load([folder,date,fname,'_graddata',ver],'data');  

if chns == 2; 
 figure(3); clf; imagesc(data{e}.PlotmRNA2); colordef black; set(gcf,'color','k'); colormap hot; colorbar; caxis([25,200]); axis off;
 end
 figure(4); clf; imagesc(data{e}.PlotmRNA); colordef black; set(gcf,'color','k'); colormap hot; colorbar; caxis([25,200]); axis off;


%%

folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/';
date = '2011-05-22/';
fname = 's05_MP06Hz_b'; ver = '_v4';
load([folder,date,fname,'_graddata',ver],'data');  


e = 10;

 figure(3); clf; colordef white; set(gcf,'color','w'); 
        plot(data{e}.Data_sort(:,1)  ,data{e}.Data_sort(:,2) - min(data{e}.Data_sort(:,2)) ,'.','color','k','MarkerSize',5); % check results    [e/Es,0,1-e/Es]
        hold on; 
        errorbar(data{e}.x,data{e}.mu(:,1)  - min(data{e}.Data_sort(:,2)), data{e}.sigma(:,1),'linestyle','none','linewidth',3,'color','k','MarkerSize',5); % [e/Es,0,1-e/Es]
        
         plot(data{e}.Data_sort(:,1),data{e}.Data_sort(:,3)- min(data{e}.Data_sort(:,3)),'.','color','g','MarkerSize',5); %[(s-1)/2,1-e/Es,e/Es]
            hold on; 
         errorbar(data{e}.x,data{e}.mu(:,2)- min(data{e}.Data_sort(:,3)), data{e}.sigma(:,2),'linestyle','none','linewidth',3,'color','g','MarkerSize',5); % [0,1-e/Es,e/Es]

           ylabel('number of mRNA transcripts per cell','Fontsize',14); xlabel('distance (nm)','Fontsize',14);
         set(gcf,'color','w'); set(gca,'Fontsize',14);
