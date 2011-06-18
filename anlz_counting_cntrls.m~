%%                          anlz_counting_cntls.m
% Alistair Boettiger                                   Date Begun: 06/07/11
% Levine Lab 



clear all;


Ttot = tic; 


disp('loading data...');
folder = '/Volumes/Data/Lab Data/Raw_Data/2011-05-22/';
% stackfolder =  's14_comp_cntrl/'; % 's12_cntrl_2label/'; % 

fname = 's12_cntrl_2label'; ver = '_2'; emb ='02';  xp1 = 50; yp1 = 900; ws = 200; hs = 200; nuc_chn = 3; nt = 90;   stackfolder = 's12_cntrl_2label/'; %      

% fname = 's14_comp_cntrl'; ver = '_1'; emb ='01';   xp1 = 130; yp1 = 1350; ws = 200; hs = 200; nuc_chn = 4; nt = 66;  stackfolder =  's14_comp_cntrl/';    

 disp_nascent = 0; shownucs = 0;
 
 filename = [folder,'/',fname];
% xp1 =1; yp1 = 1;  ws = 206; hs = 206;  % ws = 512; hs = 512;  % ws = 2048; hs = 2048; 

%---- Dot Finding Parameters ----- %
    sigmaE = 3;%  IMPORTANT
    sigmaI = 4; % IMPORTANT
  %  min_int  = 0.04;  %5  ;% .05 % not necessary Fix at Zero
    FiltSize = 30;% 
    min_size = 30;% 
    
   consec_layers = 2; 
   ovlap = 4; 
%---------------------------------%


%------------ Display Option 1 or 0 / ON or OFF -------------%
   show_projected = 1; % show max-project with all dots and linked dots.  
   plotZdata = 1 ;% show z-map of data
   getpreciseZ = 1;%
%------------------------------------------------------------%

load([filename,ver,'.mat'])  
Zs = Datas.LSM_info.DimensionZ; 
w = Datas.Stack1.Image1.IMG.width;
h = Datas.Stack1.Image1.IMG.height;


%% Set-up

xp2 = xp1+ws-1;
yp2 = yp1+hs-1;

    Imax = imreadfast([folder,stackfolder,'max_',fname,'_',emb,'.tif']); 
    Imax_dots = Imax(xp1:xp2,yp1:yp2,1:3);  
    figure(1); clf; imagesc(Imax(:,:,1:3));
    figure(2); clf; imagesc(Imax_dots);
 
    
disp(['Coordinates:  ', num2str(xp1), ' : ', num2str(xp2), ',   ' num2str(yp1), ' : ', num2str(yp2) ] );

   
    % Build the Gaussian Filter   
    Ex = fspecial('gaussian',FiltSize,sigmaE); % excitatory gaussian
    Ix = fspecial('gaussian',FiltSize,sigmaI); % inhibitory gaussian
    
%% Find all dots


% Initialize variables
 im_folder = cell(1,Zs);
 DotLabels1 = cell(1,Zs);    
 DotData1 = cell(1,Zs);    
 Inds1 = cell(1,Zs); 
 Ints1 = cell(1,Zs); 
    
 DotLabels2 = cell(1,Zs);    
 DotData2 = cell(1,Zs);    
 Inds2 = cell(1,Zs); 
 Ints2 = cell(1,Zs); 

 Inuc = zeros(ws,hs,Zs,'uint8'); 
 
 tic; disp('finding dots...'); 
for z = 1:Zs % z = 20        
      % Dot finding for channel 1
       im_folder{z} = [folder,stackfolder,fname,'_',emb,'_z',num2str(z),'.tif'];
        Iin_z = imreadfast(im_folder{z});
       
       mRNAchn = 1;   min_int =.06;  % .05; % 
       [DotLabels1{z},DotData1{z},Inds1{z},Ints1{z}] = dotfinder(Iin_z(xp1:xp2,yp1:yp2,mRNAchn),Ex,Ix,min_int,min_size);
         
       mRNAchn = 2;    min_int = .05;   % .05;    %
       [DotLabels2{z},DotData2{z},Inds2{z},Ints2{z}] = dotfinder(Iin_z(xp1:xp2,yp1:yp2,mRNAchn),Ex,Ix,min_int,min_size);   
       
       
    % Nuclear finding            
    In =  Iin_z(xp1:xp2,yp1:yp2,nuc_chn); % Im{1,z}{3}( xp1:xp2,yp1:yp2 ); 
    nucnoise = double(im2bw(In,.15)); 
    nucnoise(nucnoise==0) = NaN; 
    Inuc(:,:,z) = makeuint(In,8);
   % figure(2); clf; imagesc(Inuc(:,:,z))   
       
end
toc;
%%
 
intype = class(Iin_z);
ck_dots = tic;
        NewDotC1 = CheckDotUpDown(DotLabels1,DotData1,Inds1,Ints1,plotZdata,getpreciseZ,consec_layers,ovlap,xp1,xp2,yp1,yp2,intype);
        NewDotC2 = CheckDotUpDown(DotLabels2,DotData2,Inds2,Ints2,plotZdata,getpreciseZ,consec_layers,ovlap,xp1,xp2,yp1,yp2,intype);
 toc(ck_dots) 

 
 %% Plots! 

 % Convert cell-arrays to matrices for easy reference plotting by stacking.
D1 = cell2mat(DotData1');
D2 = cell2mat(DotData2');


    Imax = imreadfast([folder,stackfolder,'max_',fname,'_',emb,'.tif']); 
    Imax_dots = 1.5*Imax(xp1:xp2,yp1:yp2,1:3);  
    Imax_dots(:,:,3) = .1*Imax_dots(:,:,3);
   
    figure(2);  clf;  
    imagesc(Imax_dots);
    colordef black; set(gcf,'color','k'); 
    hold on;
    plot(  NewDotC1(:,1),NewDotC1(:,2),'mo','MarkerSize',5 );
    plot(  NewDotC2(:,1),NewDotC2(:,2),'co','MarkerSize',5 ); 
    plot(  D1(:,1),D1(:,2),'m+','MarkerSize',4);
    plot(  D2(:,1),D2(:,2),'c+','MarkerSize',4);

   figure(4); clf; 
   imagesc(Imax_dots(:,:,1)+Imax_dots(:,:,2));
    hold on; set(gcf,'color','k');
    plot(  NewDotC1(:,1),NewDotC1(:,2),'b+','MarkerSize',5);
    plot(  NewDotC2(:,1),NewDotC2(:,2),'go','MarkerSize',5);
    colormap hot; colorbar;

    
    
%%
   minS =10;% 8; % 10;    
   figure(3); clf; 
        subplot(2,1,1); set(gcf,'color','k');
        Imax_r = Imax_dots; Imax_r(:,:,2) = 0*Imax_r(:,:,2);
        imagesc(Imax_r); hold on;           
        plot(  NewDotC1(:,1),NewDotC1(:,2),'m+','MarkerSize',5 );

        subplot(2,1,2); 
        Imax_g = Imax_dots; Imax_g(:,:,1) = 0*Imax_g(:,:,1);
        imagesc(Imax_g); hold on;
        plot(  NewDotC2(:,1),NewDotC2(:,2),'go','MarkerSize',5 ); 





    d2 = 20*zeros(1,length(NewDotC2));
    jn = zeros(1,length(NewDotC2));
    for n = 1:length(NewDotC2)
        [d2(n),jn(n)] = min(sqrt( (NewDotC2(n,1) - NewDotC1(:,1)).^2 + (NewDotC2(n,2) - NewDotC1(:,2)).^2 + (NewDotC2(n,3) - NewDotC1(:,3)).^2) );
    end
    figure(10); clf; hist(d2);
    d2miss = sum(d2< minS)/length(d2);

    d1 = 10*zeros(1,length(NewDotC1));
    for n = 1:length(NewDotC1)
        d1(n) = min(sqrt( (NewDotC1(n,1) - NewDotC2(:,1)).^2 + (NewDotC1(n,2) - NewDotC2(:,2)).^2 + (NewDotC1(n,3) - NewDotC2(:,3)).^2)  );
    end
    figure(10); clf; hist(d1);
    d1miss = sum(d1< minS)/length(d1);

     figure(4); hold on; plot(NewDotC1((d1< minS),1),NewDotC1((d1< minS),2),'c*')

  figure(3);
      subplot(2,1,1); hold on; plot(NewDotC2((d2> minS),1),NewDotC2((d2> minS),2),'y*')
      title(['miss rate = ',num2str(d2miss,3)]);
      subplot(2,1,2); hold on; plot(NewDotC1((d1> minS),1),NewDotC1((d1> minS),2),'y*')
      title(['miss rate = ',num2str(d1miss,3)]);
  
      
%%  3D dot plotting


zspace = 330;


BothDots = [NewDotC1((d1<minS),1), NewDotC1((d1<minS),2), NewDotC1((d1<minS),3) ];
justRed = [NewDotC1((d1>minS),1), NewDotC1((d1>minS),2), NewDotC1((d1>minS),3) ];
justGreen = [NewDotC2((d2>minS),1), NewDotC2((d2>minS),2), NewDotC2((d2>minS),3) ];

%%  3D-isosurfaces view, Prep
st = 1; stp = 1;
Zmax = Zs-st+1; % nt = 66;

x = (1:stp:ws)*50;
y = (1:stp:hs)*50;
 z = Zs*zspace-zspace*(st:Zmax);

[X,Y,Z] = meshgrid(x,y,z);

data = Inuc(1:stp:end,1:stp:end,st:Zmax);
data(isnan(data)) = 0;
data = smooth3(data,'box',15);

% 3D-isosurfaces view, Draw
tic
disp('plotting data...'); 



% adust intensity threshold for nascent transcript detection
lev1 = .1;
lev2 = .1; 

figure(5); clf;
if shownucs == 1;
    patch(isosurface(Y,X,Z,data,nt,'material','dull'),'FaceColor','blue','EdgeColor','none'); hold on;  alpha .3; 
end

if disp_nascent == 1
patch(isosurface(Y,X,Z,Isect1(1:stp:end,1:stp:end,1:Zmax),lev1*25000,'material','shiny'),'FaceColor','red','EdgeColor','none'); alpha .3
patch(isosurface(Y,X,Z,Isect2(1:stp:end,1:stp:end,1:Zmax),lev2*15100,'material','shiny'),'FaceColor','green','EdgeColor','none');  
end

 camlight(5,30);
light('Position',[1000,1.5E4 8000],'Style','local');
 lighting phong; % view(112,34); view(-111,38);
view(-116,70);
set(gcf,'color','k'); set(gca,'color','k');




%%  3D-isosurfaces Draw dots


pointsR = [justRed(:,1)*50-yp1*50,justRed(:,2)*50-xp1*50,zspace*Zs-zspace*justRed(:,3)]; 
pointsG = [justGreen(:,1)*50-yp1*50,justGreen(:,2)*50-xp1*50,zspace*Zs-zspace*justGreen(:,3)]; 
pointsY = [BothDots(:,1)*50-yp1*50,BothDots(:,2)*50-xp1*50,zspace*Zs-zspace*BothDots(:,3)]; 

V1 = zeros(size(X));  
V2 = V1; 
VY = V1;

ind1 = cell(Zs,1);
ind2 = cell(Zs,1);
ind3 = cell(Zs,1);
for zs =1:Zs; 
    layer = justRed(:,3)<zs+1 & justRed(:,3)>zs -1; 
    xs = pointsR(layer,1)/(50*ws)*length(x);
    ys = pointsR(layer,2)/(50*hs)*length(y); 
    ind1{zs} = floor(xs) + length(y)*(floor(ys)) + length(x)*length(y)*(zs-1);
    
    layer = justGreen(:,3)<zs+1 & justGreen(:,3)>zs -1; 
    xs = pointsG(layer,1)/(50*ws)*length(x);
    ys = pointsG(layer,2)/(50*hs)*length(y); 
    ind2{zs} = floor(xs) + length(y)*(floor(ys)) + length(x)*length(y)*(zs-1);
    
    layer = BothDots(:,3)<zs+1 & BothDots(:,3)>zs -1; 
    xs = pointsY(layer,1)/(50*ws)*length(x);
    ys = pointsY(layer,2)/(50*hs)*length(y); 
    ind3{zs} = floor(xs) + length(y)*(floor(ys)) + length(x)*length(y)*(zs-1);
    
end
inds1 = cell2mat(ind1); inds1 = inds1(inds1>0);
inds2 = cell2mat(ind2); inds2 = inds2(inds2>0);
inds3 = cell2mat(ind3); inds3 = inds3(inds3>0);
V1(inds1) = 2;
V2(inds2) = 2;
VY(inds3) = 2;

figure(5); hold on;
 patch(isosurface(X,Y,Z,V1,1,'material','dull'),'FaceColor','red','EdgeColor','none');
patch(isosurface(X,Y,Z,V2,1,'material','dull'),'FaceColor','green','EdgeColor','none');
patch(isosurface(X,Y,Z,VY,1,'material','dull'),'FaceColor','yellow','EdgeColor','none');
camlight left; 
light('Position',[1000,1.5E4 8000],'Style','local');
 lighting phong; % view(112,34); 
view(194,80);
%  

toc








      
      
toc(Ttot); 