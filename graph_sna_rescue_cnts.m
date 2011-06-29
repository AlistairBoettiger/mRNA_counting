

  folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/2011-05-22/'; 
  fname = 'MP07het_snaD_22C';
  ver = '';
  load([folder,fname,'_graddata',ver],'data');  
  
  c1 = [.7,0,.3];
  c2 = [.2,0,.8];
  c3 = [.4,0,.6];
  
  p1 = 0; 
  MP07 = data;
  
    mcorr = 1; e = 1; st =50;
      n = 4; 
      theta = mean(data{e}.Data_sort(st:end,1));  
      A = max(data{e}.Data_sort(st:end,2));  
      b=min(data{e}.Data_sort(st:end,2)); 
      [p,fit] = fxn_fit_sigmoid(data{e}.Data_sort(st:end,1)',data{e}.Data_sort(st:end,2)',[n,theta,A,b],'r');
  p1 = p(2); 
      
      figure(20); clf; plot(data{e}.Data_sort(st:end,1), data{e}.Data_sort(st:end,2),'b');
    hold on; 
    plot(data{e}.Data_sort(st:end,1),fit,'r');
    plot(p(2),A/4,'r*','MarkerSize',20); 
      
  
  
  figure(1); clf; 

     plot(data{e}.Data_sort(:,1) ,mcorr*data{e}.Data_sort(:,2) - min(data{e}.Data_sort(:,2))  ,'.','color',c1,'MarkerSize',5); % check results  
        hold on; 
        errorbar(data{e}.x,  mcorr*data{e}.mu(:,1)- min(data{e}.Data_sort(:,2)),mcorr*data{e}.sigma(:,1),'linestyle','none','linewidth',2,'color',c1,'MarkerSize',1);

        e=2
          plot(data{e}.Data_sort(:,1) - 6E4 ,mcorr*data{e}.Data_sort(:,2) - min(data{e}.Data_sort(:,2))  ,'.','color','r','MarkerSize',5); % check results  
        hold on; 
        errorbar(data{e}.x - 6E4,  mcorr*data{e}.mu(:,1)- min(data{e}.Data_sort(:,2)),mcorr*data{e}.sigma(:,1),'linestyle','none','linewidth',2,'color',c1,'MarkerSize',1);

        
        
        
fname = 's07_MP08Hz_snaD_22C';
  ver = '';
  load([folder,fname,'_graddata',ver],'data');  
  
  
  e = 1; o1= -2.3E4;
    plot(data{e}.Data_sort(:,1) + o1,  mcorr*data{e}.Data_sort(:,2),'.','color',c2,'MarkerSize',5); % check results  
        hold on; 
        errorbar(data{e}.x + o1  ,mcorr*data{e}.mu(:,1),mcorr*data{e}.sigma(:,1),'linestyle','none','linewidth',2,'color',c2,'MarkerSize',1);

        
          e = 9; o2 = 0; 
    plot(data{e}.Data_sort(:,1) + o2,mcorr*data{e}.Data_sort(:,2),'.','color',c3,'MarkerSize',5); % check results  
        hold on; 
        errorbar(data{e}.x+ o2 ,mcorr*data{e}.mu(:,1),mcorr*data{e}.sigma(:,1),'linestyle','none','linewidth',2,'color',c3,'MarkerSize',1);
        
        ylim([0,200]);
        
        
        
 fname = 's04_MP10Hz';
   ver = '_v3';
  load([folder,fname,'_graddata',ver],'data');  
  % figure(3); clf; hold on;
  
  figure(2); clf; hold on;
  
  offsets = zeros(1,10); 

  
  for e = 4:5
      n = 4; 
      theta = mean(data{e}.Data_sort(st:end,1));  
      A = max(data{e}.Data_sort(st:end,2));  
      b=min(data{e}.Data_sort(st:end,2)); 
      [p,fit] = fxn_fit_sigmoid(data{e}.Data_sort(st:end,1)',data{e}.Data_sort(st:end,2)',[n,theta,A,b],'r'); % pause
      offsets(e) = p(2); 
      
    
      figure(1); 
    plot(data{e}.Data_sort(:,1) + p1 - offsets(e),  mcorr*data{e}.Data_sort(:,2)- min(data{e}.Data_sort(:,2)),'.','color','k','MarkerSize',5); % check results  
     errorbar(data{e}.x+ p1 - offsets(e),   mcorr*data{e}.mu(:,1)- min(data{e}.Data_sort(:,2)),mcorr*data{e}.sigma(:,1),'linestyle','none','linewidth',1,'color','k','MarkerSize',1);

    hold on;
    
    figure(2); 
    plot(data{e}.Data_sort(st:end,1)' + p1 - offsets(e),fit);
      
    figure(20); clf; plot(data{e}.Data_sort(st:end,1), data{e}.Data_sort(st:end,2),'b');
    hold on; 
    plot(data{e}.Data_sort(st:end,1),fit,'r');
    plot(p(2),A/4,'r*','MarkerSize',20); 
      
  end
  
  
  figure(1); 
set(gcf,'color','w'); set(gca,'FontSize',14);
ylabel('mRNA per cell'); xlim([-1.5E4,12E4]);
legend('no primary, hets','', 'no primary hets','', 'no shadow, Hz','', 'no shadow Hz','','wt','','wt',''s);
 
