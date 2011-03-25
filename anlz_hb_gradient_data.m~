
%% anlz_hb_gradient_data.m

% Alistair Boettiger                                Date Begun: 02/02/11
% Levine Lab                                    Last Modified: 03/30/11


  %% New data-method
  clear all;
    rawfolder =  '/Volumes/Data/Lab Data/Raw_Data/02-17-11/MP01_22C/';
  folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
   fname = 'MP09_22C_hb_y_e';
 % fname = 'MP01_22C_hb_y';
  load([folder,fname,'_slidedata'], 'Data'); 
  

  
  figure(10); clf; figure(11); clf; figure(12); clf;

  
  hbdata = cell(12,1); 
for e = 1:12
  %emb = '03'; % '08'; e = 8
    if e<10
        emb = ['0',num2str(e)];
    else
        emb = num2str(e);
    end
  
  i = str2double(emb); cor = 1; 

  
  try
  load([folder,fname,'_',emb,'_nucdata.mat']); 
  catch er
      try 
           load([rawfolder,fname,'_',emb,'_nucdata.mat']);  
      catch er
      break
      end
  end
      
              mRNAsadj = Data{i,1}.mRNAsadj;
              mRNAsadj2= Data{i,2}.mRNAsadj;

      PlotmRNA = imresize(NucLabeled,.5,'nearest');
       NucLabeled = imresize(NucLabeled,.5,'nearest'); 
                      Nnucs =    max( NucLabeled(:) );
                      for n=1:Nnucs;
                          PlotmRNA(PlotmRNA==n) = mRNAsadj(n+cor);
                      end
    figure(1); clf; imagesc(PlotmRNA); colormap hot; 
        
%% Orient image along major axis     
    bw = imresize(makeuint(PlotmRNA,16),.5,'nearest');  
    thresh = graythresh(bw);
    
    bw = im2bw(bw,thresh); 
    bw = bwareaopen(bw,500);
   % figure(1); clf; imagesc(bw);

    rprops = regionprops(bwlabel(bw),'MajorAxis','Orientation');
    
 meanvar = 100;
  rotes = 0; 
  try
  
 while meanvar > 75  && rotes<4;
    
    NucLabel = imrotate(NucLabeled,(0+90*rotes)-rprops.Orientation,'nearest'); 
    figure(1); clf; imagesc(NucLabel);
    rotes = rotes + 1; 
    
% convert nucleus centroids to indexed postions


    S = regionprops(NucLabel,'Centroid');
    nuc_cents = reshape([S.Centroid],2,length(S));
    [h,w] = size(NucLabel); 
 %    [hn,wn] = size(NucLabel); 
    c_inds = sub2ind([h,w],floor(nuc_cents(2,:)),floor(nuc_cents(1,:)));
    % compute distances
    d = nuc_cents(2,:)*h; % /hn;

    if length(c_inds) < length(mRNAsadj);
        mRNAsadj = mRNAsadj(2:end);
        mRNAsadj2 = mRNAsadj2(2:end);
    end


    maternal = min(mRNAsadj);

    %  % Plotting for troubleshooting
    %     C = false(h,w);
    %     C(c_inds) = 1; 
    %     figure(1); clf; imshow(C);

    % Sort by distance from upper left corner (max bcd).  
    nuc_order = NucLabel(c_inds);
    [b,m,n] = unique(nuc_order);
    dists = d(m);
    figure(1); clf; plot(dists,mRNAsadj - maternal,'g.');  % check results
    Data_notsort = [dists; mRNAsadj - maternal; mRNAsadj2]';
    Data_sort = sortrows(Data_notsort); 

    meanvar = std(Data_sort(10:20,2));
 end
  catch err
      disp(err.message);
      continue 
  end
 
% % Convert distances to bcd concentrations.  
% bcd = log(dists); % just logorithmic.  Don't need arbitrary scaling coeff.
% %figure(1); clf; plot(bcd,mRNA_sadj,'k.');
% Data2 = [bcd; mRNAsadj-maternal; dists]';
% Data2 = sortrows(Data2);

%%
Sects = round(sqrt(Nnucs));
Q = cell(1,Sects); 
mu = zeros(Sects,2);
sigma = zeros(Sects,2);
x = zeros(Sects,1); 
Q2 = cell(1,Sects); 


dbnd = zeros(Sects,1);
for j=1:Sects
    Q{j} = Data_sort( floor((j-1)*Nnucs/Sects) + 1: floor(j*Nnucs/Sects), 2);
    Q2{j} = Data_sort( floor((j-1)*Nnucs/Sects) + 1: floor(j*Nnucs/Sects), 3);
      
    mu(j,:) = [nanmean(Q{j}),nanmean(Q2{j})]   ;
    sigma(j,:) = [nanstd(Q{j}),nanstd(Q2{j})]   ;

    x(j) = mean(Data_sort( floor((j-1)*Nnucs/Sects) + 1: floor(j*Nnucs/Sects), 1));
    fano = sigma(j)^2/mu(j); 
end

% x = linspace(mean(dists(1:Sects)),mean(dists(end-Sects:1:end)),Sects)*50/1000;


figure(1); clf; colordef black; set(gcf,'color','k');
plot(Data_sort(:,1),flipud(Data_sort(:,2)),'w.'); % check results 
hold on; errorbar(x,flipud(mu(:,1)),flipud(sigma(:,1)),'linestyle','none','linewidth',3,'color','r');
ylabel('number of mRNA transcripts per cell','FontSize',16);
xlabel('distance (\mum)','FontSize',16);
set(gca,'FontSize',16);

figure(1); clf; colordef white; set(gcf,'color','w');
plot(Data_sort(:,1),flipud(Data_sort(:,2)),'k.'); % check results 
hold on; errorbar(x,flipud(mu(:,1)),flipud(sigma(:,1)),'linestyle','none','linewidth',3,'color','r');
ylabel('number of mRNA transcripts per cell','FontSize',16);
xlabel('distance (\mum)','FontSize',16);
set(gca,'FontSize',16);

figure(2); clf; colordef white; set(gcf,'color','w');
plot(x,flipud(sigma(:,1))./flipud(mu(:,1)),'ro','MarkerSize',10); ylim([0,1]);
hold on; plot(x,sqrt(flipud(mu(:,1)))./flipud(mu(:,1)),'ko','MarkerSize',10);
ylabel('CoV','FontSize',16); xlabel('distance (\mum)','FontSize',16);
legend('Actual CoV','Poisson CoV');
set(gca,'FontSize',16);


figure(2); clf; colordef white; set(gcf,'color','w');
plot(x,flipud(sigma(:,1).^2)./flipud(mu(:,1)),'ro','MarkerSize',10);
hold on; plot(x,flipud(mu(:,1))./flipud(mu(:,1)),'ko','MarkerSize',10);
ylabel('Fano Factor','FontSize',16); xlabel('distance (\mum)','FontSize',16);
legend('Actual','Poisson');
set(gca,'FontSize',16);


figure(10); subplot(3,4,e);  colordef white; set(gcf,'color','w');
    plot(Data_sort(:,1),Data_sort(:,2),'k.'); % check results  

%plot(Data_sort(:,1)*50/1000,Data_sort(:,2),'k.'); % check results 
%figure(1); clf; set(gcf,'color','w'); 
hold on; errorbar(x,mu(:,1),sigma(:,1),'linestyle','none','linewidth',3,'color','r');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');
 
plot(Data_sort(:,1),Data_sort(:,3),'b.'); 
hold on; errorbar(x,mu(:,2),sigma(:,2),'linestyle','none','linewidth',3,'color','c');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');



figure(11); subplot(3,4,e);  
colordef white; set(gcf,'color','w');
plot(x,sigma(:,1)./mu(:,1),'ro','MarkerSize',10); ylim([0,1]);
hold on; plot(x,sqrt(mu(:,1))./mu(:,1),'ko','MarkerSize',10);
ylabel('CoV'); xlabel('distance (\mum)');
legend('Actual CoV','Poisson CoV');


figure(12); subplot(3,4,e);  
colordef white; set(gcf,'color','w');
plot(x,sigma(:,2)./mu(:,2),'bo','MarkerSize',10); ylim([0,1]);
hold on; plot(x,sqrt(mu(:,2))./mu(:,2),'ko','MarkerSize',10);
ylabel('CoV'); xlabel('distance (\mum)');
legend('Actual CoV','Poisson CoV');

hbdata{e}.sortData = Data_sort;
hbdata{e}.mu = mu;
hbdata{e}.sigma = sigma;
end

save([folder,fname,'grad_data'], 'hbdata');
