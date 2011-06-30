

clear all;
 figure(1); clf;  figure(3); clf;
folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/';
 
slides = {'MP10Hz', 'MP06Hz' 'MP10het','MP05het'};
ver = '';

for s = 1:length(slides) 
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
      skip = 10;% [1,2,3,4,5,6,8,9,10]; 
               
    case 'MP05het'
        date = '2011-02-17/'; 
        fname = 'MP05_22C_sna_y_c'; 
        ver = '_v2';
        skip = 3;
        
    case 'MP10het'
        date = '2011-02-17/'; 
        fname ='MP10_22C_sna_y_d'; 
        ver = '_v3';% 
       skip = 9;
      
   case 'MP08Hz'
        fname =  's07_MP08Hz_snaD_22C';
        chns = 1; 
        skip = 15;

end

disp(fname); 

    load([folder,date,fname,'_graddata',ver],'data');  
Es = length(data); 
    offsets = zeros(1,Es); 


 %  
 
 % Get boundary point of first curve as a reference
 if s == 1;
 e=1;
  gx = data{e}.Data_sort(:,1);  
  grad = data{e}.Data_sort(:,2);  
  L = length(grad); st = floor(L/5);
  grad = grad(st:end); 
  gx = gx(st:end);
 
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
 
    k = 2;
 else
     k = 1;
    
 end
 
 
 
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
    catch
        continue
    end
end

%
% offsets(3) = 5E5;

figure(1);
% figure(s); clf; 
colordef white; set(gcf,'color','w');

for e = 1:Es 
    try
        if isempty(find(e==skip))
      %   subplot(4,3,e);
      figure(1); % subplot(2,2,s); title(texlabel(fname,'literal'));
        mcorr =     data{e}.Nnucs./data{1}.Nnucs;
        plot(data{e}.Data_sort(:,1) + p1-offsets(e),mcorr*data{e}.Data_sort(:,2) - min(mcorr*data{e}.Data_sort(:,2)) ,'.','color','k','MarkerSize',1); % check results    [e/Es,0,1-e/Es]
        hold on; 
        plot(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,1)  - min(mcorr*data{e}.Data_sort(:,2) ),'.','color','k','MarkerSize',5); % [e/Es,0,1-e/Es]


        
      %  errorbar(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,1) - min(mcorr*data{e}.mu(:,1)) ,mcorr*data{e}.sigma(:,1),'linestyle','none','linewidth',3,'color',[e/Es,0,1-e/Es],'MarkerSize',1);
        if chns == 2;
            plot(data{e}.Data_sort(:,1)+  p1-offsets(e),mcorr*data{e}.Data_sort(:,3)- min(mcorr*data{e}.Data_sort(:,3)),'.','color',[s/4,0,1-s/4],'MarkerSize',1); %[(s-1)/2,1-e/Es,e/Es]
            hold on; 
         plot(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,2)- min(mcorr*data{e}.Data_sort(:,3)),'o','color',[s/4,0,1-s/4],'MarkerSize',5); % [0,1-e/Es,e/Es]

        figure(3); hold on; % subplot(2,2,s); 
        % plot(data{e}.Data_sort(:,1)+  p1-offsets(e), data{e}.Data_sort(:,3)./data{e}.Data_sort(:,2) ,'.','color',[e/Es,0,1-e/Es],'MarkerSize',5); ylim([0,2]); 
         
         % plot(mean(data{e}.Data_sort(:,3)./data{e}.Data_sort(:,2)),'.','color',[e/Es,0,1-e/Es]);
         
          errorbar(data{e}.x+  p1-offsets(e),data{e}.sigma(:,2)./data{e}.mu(:,2),data{e}.bssigma(:,2)./data{e}.mu(:,2),'o','color',[s/4,0,1-s/4],'MarkerSize',5);
            % errorbar(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,2)-min(mcorr*data{e}.mu(:,2)),mcorr*data{e}.sigma(:,2),'linestyle','none','linewidth',3,'color',[0,1-e/Es,e/Es],'MarkerSize',1);
        
        end
        end
        
    catch er
        disp(er.message); 
        continue
    end
end
        ylabel('number of mRNA transcripts per cell'); xlabel('distance (nm)');
        ylim([0,250]);

    
end % end loop over slides

%%
 e=2;
switch  'MP06het'; % 'MP10het' % e =2;  % 'MP05het'; % e =7;

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
      skip = 10;% [1,2,3,4,5,6,8,9,10]; 
               
    case 'MP05het'
        date = '2011-02-17/'; 
        fname = 'MP05_22C_sna_y_c'; 
        ver = '_v2';
        skip = 3;
        
    case 'MP10het'
        date = '2011-02-17/'; 
        fname ='MP10_22C_sna_y_d'; 
        ver = '_v3';% 
       skip = 9;
      
   case 'MP08Hz'
        fname =  's07_MP08Hz_snaD_22C';
        chns = 1; 
        skip = 15;

end

disp(fname); 

    load([folder,date,fname,'_graddata',ver],'data');  


 figure(3); clf; imagesc(data{e}.PlotmRNA2); colordef black; set(gcf,'color','k'); colormap hot; colorbar; caxis([25,100]); axis off;
 figure(4); clf; imagesc(data{e}.PlotmRNA); colordef black; set(gcf,'color','k'); colormap hot; colorbar; caxis([25,200]); axis off;


