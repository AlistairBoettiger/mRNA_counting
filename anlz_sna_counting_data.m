
%% anlz_sna_gradient_data.m

% Alistair Boettiger                                Date Begun: 02/02/11
% Levine Lab                                    Last Modified: 06/23/11


%% Modifications
% Modified 06/12/11 to scale green channel to correct for missed detection
% rate
% Modified 06/14/11 to save oriented mRNA plots with corrected counts.  
% 6/23/11 converted from hb to sna. 


  %% New data-method
  clear all;
  maternal =0;  ver = '';  Es = 14; cor = 0; % defaults
  rawfolder = '/Volumes/Data/Lab Data/Raw_Data/2011-05-22/s07_MP08/';
  folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/mRNA_counting/Data/2011-05-22/'; 
  fname =  's07_MP08Hz_snaD_22C'; missG = 1.3; %1.3; 
  
  chns = 1; 
 
  
  figure(10); clf; figure(11); clf;
  
  data = cell(Es,1); 
for e = 1:Es %  
    if e<10
        emb = ['0',num2str(e)];
    else
        emb = num2str(e);
    end
  
  i = str2double(emb);


      try 
           load([rawfolder,fname,'_',emb,'_nucdata.mat']);  
      catch er
          disp(er.message);
         continue
      end

          
          try
              load([folder,fname,'_',emb,'_chn1','_data',ver,'.mat']); 
              mRNAsadj = mRNA_sadj; % mRNA_cnt./nuc_area;
              if chns ==2
                load([folder,fname,'_',emb,'_chn2','_data',ver,'.mat']); 
                mRNAsadj2 = mRNA_sadj; % mRNA_cnt./nuc_area;  % 
              end
          catch er
              disp(er.message);
              disp('unable to load image');
              break
          end
    
              
       PlotmRNA = imresize(NucLabeled,.5,'nearest');
       
       if chns == 2; 
            PlotmRNA2 = PlotmRNA; 
       end
       NucLabeled = imresize(NucLabeled,.5,'nearest'); 
                      Nnucs =    max( NucLabeled(:) );
                      for n=1:Nnucs;
                          PlotmRNA(PlotmRNA==n) = mRNAsadj(n+cor);
                          if chns == 2
                            PlotmRNA2(PlotmRNA2==n) = mRNAsadj2(n+cor);
                          end
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
    
    NucLabel = imrotate(NucLabeled,(0+90*rotes)-rprops(1).Orientation,'nearest'); 
    figure(1); clf; imagesc(NucLabel);
    rotes = rotes + 1; 
    
    

    
% convert nucleus centroids to indexed postions


    S = regionprops(NucLabel,'Centroid');
    nuc_cents = reshape([S.Centroid],2,length(S));
    [h,w] = size(NucLabel); 
     c_inds = sub2ind([h,w],floor(nuc_cents(2,:)),floor(nuc_cents(1,:)));
    % compute distances
    d = nuc_cents(2,:)*h; % /hn;

    if length(c_inds) < length(mRNAsadj);
        mRNAsadj = mRNAsadj(2:end);
        if chns == 2;
            mRNAsadj2 = mRNAsadj2(2:end);
        end
    end


   % maternal =0; % min(mRNAsadj); %0; % 

     
    %  % Plotting for troubleshooting
    %     C = false(h,w);
    %     C(c_inds) = 1; 
    %     figure(1); clf; imshow(C);

    % Sort by distance from upper left corner (max bcd).  
    nuc_order = NucLabel(c_inds);
    [b,m,n] = unique(nuc_order);
    dists = d(m);
    figure(1); clf; plot(dists,mRNAsadj ,'g.');  % check results
    if chns == 2
        Data_notsort = [dists; mRNAsadj ; missG*(mRNAsadj2)]';
    else
        Data_notsort = [dists; mRNAsadj]';
    end
    Data_sort = sortrows(Data_notsort); 
    meanvar = std(Data_sort(10:20,2));
 end
 
 
 % show rotated image
    % figure(2); clf; 
        PlotmRNA_r = imrotate(PlotmRNA,(0+90*rotes)-rprops(1).Orientation,'nearest'); 
    if chns == 2
        PlotmRNA2_r = imrotate(missG*PlotmRNA2,(0+90*rotes)-rprops(1).Orientation,'nearest'); 
        imagesc(PlotmRNA2_r); colormap hot;
    end
 
  catch err
      disp(err.message);
      continue 
 end
 

%% Cluster cells based on distance from anterior pole
Sects = round(sqrt(Nnucs));
Q = cell(1,Sects); 
if chns == 2;
    mu = zeros(Sects,2);
    sigma = zeros(Sects,2);
    bsmu = zeros(Sects,2);
    bssigma = zeros(Sects,2);
    x = zeros(Sects,1); 
    Q2 = cell(1,Sects); 
else
    mu = zeros(Sects,1);
    sigma = zeros(Sects,1);
    bsmu = zeros(Sects,1);
    bssigma = zeros(Sects,1);
    x = zeros(Sects,1); 
end

dbnd = zeros(Sects,1);
for j=1:Sects
    Q{j} = Data_sort( floor((j-1)*Nnucs/Sects) + 1: floor(j*Nnucs/Sects), 2);
    
    if chns == 2 
    Q2{j} = Data_sort( floor((j-1)*Nnucs/Sects) + 1: floor(j*Nnucs/Sects), 3);
    mu(j,:) = [nanmean(Q{j}),nanmean(Q2{j})]   ;
    bsmu(j,:) =  [bserr(Q{j},'nanmean'),bserr(Q2{j},'nanmean')] ;
    sigma(j,:) = [nanstd(Q{j}),nanstd(Q2{j})]   ;
    bssigma(j,:) =  [bserr(Q{j},'nanstd'),bserr(Q2{j},'nanstd')] ;
    else
    mu(j,:) = [nanmean(Q{j})]   ;
    bsmu(j,:) =  [bserr(Q{j},'nanmean')] ;
    sigma(j,:) = [nanstd(Q{j})]   ;
    bssigma(j,:) =  [bserr(Q{j},'nanstd')] ;
    end

    x(j) = mean(Data_sort( floor((j-1)*Nnucs/Sects) + 1: floor(j*Nnucs/Sects), 1));
    fano = sigma(j)^2/mu(j); 
end

% x = linspace(mean(dists(1:Sects)),mean(dists(end-Sects:1:end)),Sects)*50/1000;

% fix orientation
or = mu(1,1) -  mu(end,1);
if or<0
   Data_sort(:,2) = flipud( Data_sort(:,2));
   mu(:,1) = flipud(mu(:,1));
   sigma(:,1) = flipud(sigma(:,1));
   bssigma(:,1) = flipud(bssigma(:,1)); 
   PlotmRNA_r = fliplr(PlotmRNA_r);

   if chns == 2
        Data_sort(:,3) = flipud( Data_sort(:,3));
        mu(:,2) = flipud(mu(:,2));
        sigma(:,2) = flipud(sigma(:,2));
        bssigma(:,2) = flipud(bssigma(:,2));
        PlotmRNA2_r = fliplr(PlotmRNA2_r); 
   end   
end

%% Plotting



figure(10); subplot(4,4,e);  colordef white; set(gcf,'color','w');
    plot(Data_sort(:,1),Data_sort(:,2),'k.'); % check results  
    hold on; 
    errorbar(x,mu(:,1),sigma(:,1),'linestyle','none','linewidth',3,'color','r');
    if chns == 2;
        plot(Data_sort(:,1),Data_sort(:,3),'b.'); 
        hold on; 
        errorbar(x,mu(:,2),sigma(:,2),'linestyle','none','linewidth',3,'color','c');
    end
    ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');
    title(['Nuclei = ',num2str(Nnucs)]);


figure(11); subplot(4,4,e);  
colordef white; set(gcf,'color','w');
errorbar(x,sigma(:,1)./mu(:,1),bssigma(:,1)./mu(:,1),  'r.','MarkerSize',10); 
hold on;
if chns == 2;
    errorbar(x,sigma(:,2)./mu(:,2),bssigma(:,2)./mu(:,2),  'b.','MarkerSize',10); ylim([0,1]);
end
plot(x,sqrt(mu(:,1))./mu(:,1),'k.','MarkerSize',10);
ylabel('CoV'); xlabel('distance (\mum)');
if chns == 2
    legend('sna CoV','y CoV','Poisson CoV');
else
  legend('sna CoV','Poisson CoV');  
end
title(['Nuclei = ',num2str(Nnucs)]);
 

    
data{e}.Data_sort = Data_sort;
data{e}.mu = mu;
data{e}.bsmu = bsmu;
data{e}.sigma = sigma;
data{e}.bssigma = bssigma;
data{e}.x = x; 
data{e}.Nnucs = Nnucs; 

data{e}.PlotmRNA = PlotmRNA_r; 
if chns == 2
    data{e}.PlotmRNA2 = PlotmRNA2_r;
end

end


save([folder,fname,'_graddata','ver'],'data'); 

%%

e =1;
figure(1); clf; imagesc(data{e}.PlotmRNA); colordef black; set(gcf,'color','k'); colormap hot; colorbar; caxis([15,110]); axis off;
figure(2); clf; imagesc(data{2}.PlotmRNA); colordef black; set(gcf,'color','k'); colormap hot; colorbar; caxis([15,300]); axis off;

