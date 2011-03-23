%%                          Anlz_mRNA_deg.m
%
% Alistair Boettiger                                   Date Begun: 03/22/11
% Levine Lab                                        Last Modified: 03/22/11

%% Description
% Load output of Data files from Unsupervised_DotFinding and compare
% different embryos
%

clear all; 


folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
rawfolder = '/Volumes/Data/Lab Data/Raw_Data/02-06-11/MP10_22C/';

 emb = '01'; 

fname =   'MP10_22C_sna_y_d1';  load([folder,fname,'_Slidedata'], 'Data');  %
M{1}  = Data;
load([rawfolder,fname,'_',emb,'_nucdata.mat']); 
NucData{1} = NucLabeled;
fname =   'MP10_22C_sna_y_d2';  load([folder,fname,'_Slidedata'], 'Data');  %
M{2} = Data;
load([rawfolder,fname,'_',emb,'_nucdata.mat']); 
NucData{2} = NucLabeled;
fname =   'MP10_22C_sna_y_d3';  load([folder,fname,'_Slidedata'], 'Data');  %
M{3} = Data; 
load([rawfolder,fname,'_',emb,'_nucdata.mat']); 
NucData{3} = NucLabeled;

clear Data;

cor = 1;

Nnucs = zeros(3,1); 
on_cnts = cell(3,2); 
off_cnts = cell(3,2); 
for m = 1:3
       Nnucs(m) =  max(NucData{m}(:));

  mRNAsadj1 = M{m}{1,1}.mRNAsadj;
  mRNAsadj2 = M{m}{1,2}.mRNAsadj;

  % build mRNA count matrices
  PlotmRNA1 = NucData{m};
  PlotmRNA2 = NucData{m};
  for n=1:Nnucs(m)
    PlotmRNA1(PlotmRNA1 == n) = mRNAsadj1(n+cor);
    PlotmRNA2(PlotmRNA2 == n) = mRNAsadj2(n+cor);  
  end
  
 
  spread = 1.5; x = linspace(0,450,50); 
   [on_cnts{m,1},~]= fxn_regionvar(NucData{m},PlotmRNA1,mRNAsadj1,.6,spread,Nnucs(m));
   [~,off_cnts{m,1}]= fxn_regionvar(NucData{m},PlotmRNA1,mRNAsadj1,.6,spread,Nnucs(m));
    figure(10+m); clf; subplot(2,2,1);
    hist(on_cnts{m,1},x); 
    title(['mean=',num2str(mean(on_cnts{m,1})*length(on_cnts{m,1}))]);
      xlim([min(x),max(x)]);
    subplot(2,2,2); 
    hist(off_cnts{m,1},x); xlim([min(x),max(x)/2]); 
    title(['mean=',num2str(mean(off_cnts{m,1})*length(off_cnts{m,1}))]);
    
    [on_cnts{m,2},~]= fxn_regionvar(NucData{m},PlotmRNA2,mRNAsadj2,.6,spread,Nnucs(m));
    [~,off_cnts{m,2}]= fxn_regionvar(NucData{m},PlotmRNA2,mRNAsadj2,.4,spread,Nnucs(m));
    figure(10+m);  subplot(2,2,3); 
    hist(on_cnts{m,2},x); 
    title(['mean=',num2str(mean(on_cnts{m,2})*length(on_cnts{m,2}))]);
     xlim([min(x),max(x)/2]);
    subplot(2,2,4); 
    hist(off_cnts{m,2},x); xlim([min(x),max(x)/2]);
    title(['mean=',num2str(mean(off_cnts{m,2})*length(off_cnts{m,2}))]);
    
end   
   % Needs to relative to the population of on nuclei.
   
   %%
   figure(1); clf; 
   x = linspace(0,1,20); 
   C = hsv(4); 
   for m=1:3   
       % subplot(3,1,m); 
       histdata{m} = hist(on_cnts{m,2}*length(on_cnts{m,2})/2E4,x); hold on;
       data{m} = on_cnts{m,2}*length(on_cnts{m,2})/2E4;
      % bar(x,histdata{m}); 
       %bar(x,histdata{m},'FaceColor',C(m,:));      
        plot(x,histdata{m},'Marker','s','Color',C(m,:),'LineStyle','none');     
       xlim([min(x),max(x)]);
         title(['mean=',num2str(mean(on_cnts{m,2})*length(on_cnts{m,2}))]);
   end
   alpha = .5;
   
   figure(2); clf; colordef black;
   names = {'early';'mid';'late'}
   xlab = 'mRNA density';
   F = 15; 
   cityscape(data,names,xlab,F,[25,26,27]);
   set(gcf,'color','k');
   
   
   
   
