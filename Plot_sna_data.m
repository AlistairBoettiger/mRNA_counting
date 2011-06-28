

clear all;


slides = {'MP10Hz', 'MP06Hz'}


for s = 1:length(slides) 
switch slides{s}

    case 'MP10Hz'

      rawfolder = '/Volumes/Data/Lab Data/Raw_Data/2011-05-22/s04_MP10/';% s05_MP06/'   ; %   s07_MP08/'
      folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/2011-05-22/'; 
      fname = 's04_MP10Hz';% 's05_MP06Hz' ;%  's07_MP08Hz_snaD_22C';
      chns = 2;

    case 'MP06Hz'
      rawfolder = '/Volumes/Data/Lab Data/Raw_Data/2011-05-22/s05_MP06/'   ; %   s07_MP08/'
      folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/2011-05-22/'; 
      fname =  's05_MP06Hz' ;%  's07_MP08Hz_snaD_22C';
      chns = 2;

end

disp(fname); 

    load([folder,fname,'_graddata','ver'],'data');  
    Es =  length(data);
    offsets = zeros(1,Es); 


 %%  
 
 % Get boundary point of first curve as a reference
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
  
 
for e= 2:Es;
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

    offsets(e) = p(2); 

    figure(20); clf; plot(gx', grad,'b');
    hold on; 
    plot(gx,fit,'r');
    plot(p(2),A/4,'r*','MarkerSize',20);
    catch
        continue
    end
end

%%
% offsets(3) = 5E5;

figure(s); clf; 
colordef white; set(gcf,'color','w');

for e = 1:Es 
    try
        
        mcorr =   data{e}.Nnucs./data{1}.Nnucs;
        plot(data{e}.Data_sort(:,1) + p1-offsets(e),mcorr*data{e}.Data_sort(:,2),'.','color',[e/Es,0,1-e/Es],'MarkerSize',5); % check results  
        hold on; 
        errorbar(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,1),mcorr*data{e}.sigma(:,1),'linestyle','none','linewidth',1,'color',[e/Es,0,1-e/Es],'MarkerSize',1);
        if chns == 2;
            plot(data{e}.Data_sort(:,1)+  p1-offsets(e),mcorr*data{e}.Data_sort(:,3),'.','color',[0,1-e/Es,e/Es],'MarkerSize',1); 
            hold on; 
            errorbar(data{e}.x+ p1-offsets(e),mcorr*data{e}.mu(:,2),mcorr*data{e}.sigma(:,2),'linestyle','none','linewidth',1,'color',[0,1-e/Es,e/Es],'MarkerSize',1);
        end
    catch
        continue
    end
end
        ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');
        ylim([0,350]);

    
end % end loop over slides


