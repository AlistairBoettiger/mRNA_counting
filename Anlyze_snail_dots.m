

%% Anlyze_snail_dots.m

% Alistair Boettiger                                   Date Begun: 03/14/11
% Levine Lab                                        Last Modified: 03/23/11

%% Description
% Load output of Data files from Unsupervised_DotFinding and compare
% different embryos from MP05 and MP10 slides.  
%

 clear all;

clear ons offs


folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
 %fname =   'MP10_22C_sna_y_c'; cor = 0;  load([folder,fname,'_Slidedata'], 'Data');  %
%  fname =  'MP05_22C_sna_y'; cor = 1; load([folder,fname,'_slidedata'], 'Data'); 
 
%  fname =  'MP10_22C_sna_y_d'; load([folder,fname,'_slidedata'], 'Data'); 

 % fname =  'MP05_22C_sna_y_c';  cor = 1;  emb = '08';  %  emb = '04'; %   emb = '07'; 
 % fname =  'MP05_22C_sna_y_c';  cor = 1;  emb = '01';
 % fname = 'MP10_22C_sna_y_d'; cor = 1; emb = '05';
% 
 % fname =  'MP10_22C_sna_y_c';  emb = '01'; % emb = '02';% emb = '03';%
 % emb = '04';%  emb = '05'; % emb = '07';  %** emb = '06';%**  %   emb = '08'; %
% load([folder,fname,'_slidedata'], 'Data'); 


%%



%load([folder,fname,'_Slidedata'], 'Data');  % fixed mRNA sadj

% load([folder,fname,'_slidedata_3te'], 'Data'); 
% load([folder,'MP05_data_b'], 'Data');

% Quick look at data distributions 

N = 12;
ave = zeros(N,2);
stdev = zeros(N,2); 
ons = cell(N,2);
offs = cell(N,2);

 h = 2048; w = 2048;

t1 = .6; 
spread = 1.4;
minObj = 10000; 
stre = 25; 

for slide = 1:2
    if slide == 1
      % fname =   'MP10_22C_sna_y_c'; cor = 0;  load([folder,fname,'_Slidedata'], 'Data');  % 
       fname =   'MP10_22C_sna_y_d'; cor = 1;  load([folder,fname,'_slidedata_v4']);  % 
    elseif slide == 2
        fname =  'MP05_22C_sna_y_c'; cor = 1; load([folder,fname,'_slidedata_v2'], 'Data'); 
    end

    figure(1); clf;
    for i=1:N   
        for chn = 1:2 
                try
                    if i<10
                        emb = ['0',num2str(i)];
                    else
                        emb = num2str(i);
                    end            
                    
                    mRNAsadj = Data{i,chn}.mRNAsadj;
                    load([folder,fname,'_',emb,'_nucdata.mat']); 
                                                   
                 catch err
                     disp(err.message); 
                     break
                end

                try  
                       PlotmRNA = imresize(NucLabeled,.5,'nearest');
                       NucLabel = imresize(NucLabeled,.5,'nearest'); 
                      Nnucs =    max( NucLabel(:) );
                      for n=1:Nnucs;
                          PlotmRNA(PlotmRNA==n) = mRNAsadj(n+cor);
                      end
                      [ons{i,chn},offs{i,chn}] = fxn_regionvar(NucLabel,PlotmRNA,mRNAsadj,t1,spread,Nnucs,minObj,stre);
                catch err 
                     disp(err.message); 
%                       Nnucs =    max( Data{i,chn}.NucLabeled(:) );
%                      [ons{i,chn},offs{i,chn}] = fxn_regionvar(Data{i,chn}.NucLabeled,Data{i,chn}.PlotmRNA,mRNAsadj,t1,spread,Nnucs,minObj,stre); 
                end
        end
    end
   
    if slide == 1
        MP10_ons = ons;
        MP10_offs = offs;
        clear ons offs
    elseif slide == 2
        MP05_ons = ons;
        MP05_offs = offs;
        clear ons offs;
    end
end    
    
 %%   
    x = linspace(0,450,50); 
    figure(1); clf; colordef black; set(gcf,'color','k');
    for e=1:5
        figure(1); subplot(4,2,e); 
         hist(MP05_ons{e,2},x);
         set(findobj(gca,'Type','patch'),'FaceColor','y');
         hold on; hist(MP05_ons{e,1},x); xlim([0,450]);
         title(['y05=',num2str(mean(MP05_ons{e,2}),3),...
             ' (',num2str(std(MP05_ons{e,2}),3) ')  sna=',num2str(mean(MP05_ons{e,1}),3),...
             ' (',num2str(std(MP05_ons{e,1}),3), ')']);
    end
   
       figure(22); clf;  colordef black; set(gcf,'color','k');
    for e=1:10
        figure(22); subplot(5,2,e); 
         hist(MP10_ons{5,2},x);
         set(findobj(gca,'Type','patch'),'FaceColor','y');
         hold on; hist(MP10_ons{e,1},x);  xlim([0,450]);
          title(['y10=',num2str(mean(MP10_ons{e,2}),3),...
              ' (',num2str(std(MP10_ons{e,2}),3), ')  sna=',num2str(mean(MP10_ons{e,1}),3),...
              ' (',num2str(std(MP10_ons{e,1}),3), ')' ]);

    end
    
    
        %%

% 
% chn = 2; i = str2double(emb);
% load([folder,fname,'_',emb,'_nucdata.mat']); 
%    
%               mRNAsadj = Data{i,chn}.mRNAsadj;
% 
%   PlotmRNA = imresize(NucLabeled,.5,'nearest');
%   NucLabel = imresize(NucLabeled,.5,'nearest'); 
%                       Nnucs =    max( NucLabeled(:) );
%                       for n=1:Nnucs;
%                           PlotmRNA(PlotmRNA==n) = mRNAsadj(n+cor);
%                       end
%                      
%                       chn_names = {'snail';'reporter'};
%                       
%                       figure(1); clf; colordef black;
%                       imagesc(PlotmRNA); colormap hot; colorbar;
%                       set(gcf,'color','k'); axis off;
%                       title([chn_names{chn},' ',texlabel(fname,'literal')],'FontSize',14);
%                       set(gca,'FontSize',14);
%           
%   [ons,offs] = fxn_regionvar(NucLabel,PlotmRNA,mRNAsadj,0,spread,Nnucs);
%    
% 
%  %% TOT Histograms 
%   
% 
%        mRNAsadj1 = Data{i,1}.mRNAsadj;
%       mRNAsadj2 = Data{i,2}.mRNAsadj;
%       
%       ms = linspace(0,400,50);
%       
%       M1 = hist(mRNAsadj1,ms);
%       M2 = hist(mRNAsadj2,ms);
%       
%       figure(1); clf; colordef black; set(gcf,'color','k');
%       set(gca,'FontSize',14);
%       bar(ms,M1,'r','EdgeColor','k');alpha .7; hold on; 
%       bar(ms,M2,'g','EdgeColor','k'); 
%       xlim([-5,max(ms)]);
%   
%   
%   
%         