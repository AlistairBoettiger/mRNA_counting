

 clear all;
    rawfolder =  '/Volumes/Data/Lab Data/Raw_Data/02-17-11/MP09_22C/';
  folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/'; 
   fname = 'MP09_22C_hb_y_e';
 % fname = 'MP01_22C_hb_y';
  % load([folder,fname,'_slidedata'], 'Data'); 
  
  
  
  emb = '05';
  Imax = imread([rawfolder,'max_',fname,'_',emb,'.tif']);
  
  [h,w] = size(Imax(:,:,1)); 
  load([folder,fname,'_',emb,'_nucdata.mat']); 
  
  figure(4); clf; imagesc(Cell_bnd);
  CellBnd = imresize(Cell_bnd,4,'nearest'); 
  CellBnd = imerode(CellBnd,strel('disk',1));
  figure(4); clf; imagesc(CellBnd);

  % % Max project
%   Idisp = Imax(:,:,1);
%   Idisp(:,:,1) = Idisp(:,:,1) + makeuint(CellBnd,16);
% 
%   figure(1); clf; imagesc(Idisp); 
%   figure(2); clf; colordef black; set(gcf,'color','k'); axis square; 
%   imagesc(Idisp(1:300,100:400,:)); colormap hot; colorbar; caxis([0,5E4]);
%     figure(3); clf; colordef black; set(gcf,'color','k'); axis square; 
%     imagesc(Idisp(h-300:h,w-400:w-100,:)); colormap hot; caxis([0,5E4]); colorbar;
    
    z = 22; % 16;
      Ilayer = imread([rawfolder,fname,'_',emb,'_z',num2str(z),'.tif']);
      
        Idisp = Ilayer(:,:,1);
  Idisp = Idisp +  makeuint(CellBnd,16);

  
   figure(1); clf; imagesc(Idisp); 
  
  % 
  figure(2); clf; colordef black; set(gcf,'color','k'); 
  imagesc(Idisp(41:340,101:400,:)); colormap hot; caxis([0,5E4]); axis square; axis off; %  colorbar; 
    figure(3); clf; colordef black; set(gcf,'color','k');
    imagesc(Idisp(h-400+1:h-100,w-500+1:w-200,:)); colormap hot; caxis([0,5E4]);  axis square;  axis off; %  colorbar;
  
%%
  
  emb = '08';
  Imax = imread([rawfolder,'max_',fname,'_',emb,'.tif']);
  
  figure(1); clf; imagesc(Imax);
  
  [h,w] = size(Imax(:,:,1)); 
  load([folder,fname,'_',emb,'_nucdata.mat']); 
  
  figure(4); clf; imagesc(Cell_bnd);
  CellBnd = imresize(Cell_bnd,4,'nearest'); 
  CellBnd = imerode(CellBnd,strel('disk',1));
  figure(4); clf; imagesc(CellBnd);

  % % Max project
%   Idisp = Imax(:,:,1);
%   Idisp(:,:,1) = Idisp(:,:,1) + makeuint(CellBnd,16);
% 
%   figure(1); clf; imagesc(Idisp); 
%   figure(2); clf; colordef black; set(gcf,'color','k'); axis square; 
%   imagesc(Idisp(1:300,100:400,:)); colormap hot; colorbar; caxis([0,5E4]);
%     figure(3); clf; colordef black; set(gcf,'color','k'); axis square; 
%     imagesc(Idisp(h-300:h,w-400:w-100,:)); colormap hot; caxis([0,5E4]); colorbar;
    
    z =20; % 22; % 16;
      Ilayer = imread([rawfolder,fname,'_',emb,'_z',num2str(z),'.tif']);
      
        Idisp = Ilayer(:,:,1);
  Idisp = Idisp +  makeuint(CellBnd,16);

  
  % figure(1); clf; imagesc(Idisp); 
  
  y1 = 200; x1 = 600; y2 = h-600; x2 = 500; 
  figure(2); clf; colordef black; set(gcf,'color','k'); 
  imagesc(Idisp(y1:y1+299,x1:x1+299,:)); colormap hot; caxis([0,5E4]); axis square; axis off; %  colorbar; 
    figure(3); clf; colordef black; set(gcf,'color','k');
    imagesc(Idisp(y2:y2+299,x2:x2+299,:)); colormap hot; caxis([0,5E4]);  axis square;  axis off; %  colorbar;
  
