
%% anlz_hb_gradient_data.m

% Alistair Boettiger                                Date Begun: 02/02/11
% Levine Lab                                    Last Modified: 03/30/11


  %% New data-method
  clear all;
   
  folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
  fname = 'MP09_22C_hb_y_e';
  emb = '03'; % '05'; 
  i = str2double(emb); cor = 1; 

  chn = 1;
  
  load([folder,fname,'_',emb,'_nucdata.mat']); 
  load([folder,fname,'_slidedata'], 'Data'); 
   
              mRNAsadj = Data{i,chn}.mRNAsadj;

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
    figure(8); clf; imagesc(bw);

    L = bwlabel(bw);
    rprops = regionprops(L,'MajorAxis','Orientation');
    
 meanvar = 100;
  rotes = 0;  
 while meanvar > 50  && rotes<4;
    
    NucLabel = imrotate(NucLabeled,(0+90*rotes)-rprops.Orientation,'nearest'); 
    figure(8); clf; imagesc(NucLabel);
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
    Data = [dists; mRNAsadj - maternal]';
    Data_sort = sortrows(Data); 

    meanvar = std(Data_sort(10:20,2));
 end
 
% Convert distances to bcd concentrations.  
bcd = log(dists); % just logorithmic.  Don't need arbitrary scaling coeff.
%figure(1); clf; plot(bcd,mRNA_sadj,'k.');
Data2 = [bcd; mRNAsadj-maternal; dists]';
Data2 = sortrows(Data2);

%%
Sects = round(sqrt(Nnucs));
mu = zeros(1,Sects);
sigma = zeros(1,Sects);
Q = cell(1,Sects); 

dbnd = zeros(Sects,1);
for j=1:Sects
    Q{j} = Data2( floor((j-1)*Nnucs/Sects) + 1: floor(j*Nnucs/Sects), 2);
    mu(j) = nanmean(Q{j});
    sigma(j) = nanstd(Q{j});
    fano = sigma(j)^2/mu(j); 
end

x = linspace(0,max(dists)*50/1000,Sects);


% figure(1); clf; colordef black; set(gcf,'color','k');
% plot(Data_sort(:,1)*50/1000,flipud(Data_sort(:,2)),'w.'); % check results 
% hold on; errorbar(x,fliplr(mu),fliplr(sigma),'linestyle','none','linewidth',3,'color','r');
% ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');
% set(gca,'FontSize',14);


figure(1); clf; colordef white; set(gcf,'color','w');
plot(Data_sort(:,1)*50/1000,Data_sort(:,2),'k.'); % check results 
hold on; errorbar(x,mu,sigma,'linestyle','none','linewidth',3,'color','r');
ylabel('number of mRNA transcripts per cell'); xlabel('distance (\mum)');


figure(5); clf; 
colordef white; set(gcf,'color','w');
plot(mu,sigma,'k.'); hold on; plot(mu,sigma.*sqrt([0,diff(mu)]),'g.'); 
 xlabel('mean count');

figure(2); clf; 
colordef white; set(gcf,'color','w');
plot(x,sigma,'r.'); 
hold on; plot(x,mu,'m.'); legend('\sigma','\mu','Location','Best');
xlabel('distance (\mum)');

% figure(3); clf; 
% colordef white; set(gcf,'color','w');
% plot(diff(sigma),'r-'); 
% hold on; plot(diff(mu),'m-'); legend('\sigma','\mu');
% xlabel('distance (\mum)');


figure(3); clf;
plot(sigma.^2./mu,'r');
hold on; plot(mu/10,'m');



figure(4); clf; 
colordef white; set(gcf,'color','w');
plot(x,sigma./mu,'ro','MarkerSize',10); ylim([0,1]);
hold on; plot(x,sqrt(mu)./mu,'ko','MarkerSize',10);
ylabel('CoV'); xlabel('distance (\mum)');
legend('Actual CoV','Poisson CoV');

