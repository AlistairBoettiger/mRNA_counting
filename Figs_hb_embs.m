

%%
fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Results/';
folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/';
edat =  'MP09_22C_y_hb34_data.mat';
 edat = 'BAC09_30C_y_hb01_data.mat'; % cc12

load([folder,edat]); 

Nnucs = length(cent);
Hreg = ismember(H,1:250);

I = handles.It;  
[h,w] = size(I(:,:,1));
 Iz = uint8(zeros(h,w,3));
 Iz(:,:,1) = 255*uint8(Nuc_cntr-Hreg);
 Iz(:,:,2) = uint8(255*Nuc_cntr);
 Iz(:,:,3) = uint8(255*Nuc_cntr);

figure(4); clf; imshow(Iz);
Io = imresize(Iz,.2);
imwrite(Io,[fout,'hb_reg.tif']);


inR = 1:250;
Roff = rand(1,250)>.9;
inR_off = inR(Roff); 

Hreg_off = ismember(H,inR_off);

I = handles.It;  
[h,w] = size(I(:,:,1));
 Iz = uint8(zeros(h,w,3));
 Iz(:,:,1) = 255*uint8(Nuc_cntr-Hreg+Hreg_off);
 Iz(:,:,2) = 255*uint8(Nuc_cntr-Hreg_off);
 Iz(:,:,3) = 255*uint8(Nuc_cntr-Hreg_off);

figure(4); clf; imshow(Iz);
Io = imresize(Iz,.2);
imwrite(Io,[fout,'hb_reg_sim2.tif']);



%%

% 
fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Results/';
folder = '/Volumes/Data/Lab Data/Shadow_data';
im1 = 'MP09_22C_y_hb_09.tif'; % no transcription during mitotic division.
% im1 ='BAC09_22C_y_hb_12.tif'; 
im1 = 'BAC01_30C_y_hb_08.tif';

I = imread([folder,'/',im1]);  

C = [0,0,0;
    0,0,0;
    1,0,0];

T = [.1,1;
    .12,1;
    .23,1];

f = [0,0];

I1 = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I1);

Io = imresize(I1,.4);

imwrite(Io,[fout,'cc10_mitotic.tif']);
%%


fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Results/';
folder = '/Volumes/Data/Lab Data/Shadow_data';
im1 = 'MP09_22C_y_hb_09.tif'; % no transcription during mitotic division.
% im1 ='BAC09_22C_y_hb_12.tif'; 

I = imread([folder,'/',im1]);  

C = [0,0,0;
    0,0,0;
    1,0,0];

T = [.1,1;
    .12,1;
    .23,1];

f = [0,0];

I1 = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I1);

Io = imresize(I1,.4);

imwrite(Io,[fout,'cc11_mitotic.tif']);


%%
 edat = 'BAC09_30C_y_hb01_data.mat'; % cc12
 edat =  'MP09_22C_y_hb34_data.mat';
load([folder,edat]); 

       age = getage(H,cent);
       disp(age);

I = handles.It;  

C = [1.5,0,0;
    0,0,0;
    0,0,.8];

T = [.05,1;
    .1,1;
    .0,1];

f = [0,0];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);


L1n2 = L1 - L2 > 0 ;

8/length(pts1)
11/(length(cent)-length(pts1)+11)


[h,w] = size(I(:,:,1));
 Iz = uint8(zeros(h,w,3));
    Iz(:,:,1) = imadd(uint8(0*Reg1.*L1n2.*Cell_bnd.*Reg1),I(:,:,1)) ;
    Iz(:,:,2) =  imadd(uint8(0*Reg1.*L2n1a.*Cell_bnd),1*I(:,:,2)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Iz(:,:,3) =    imadd(uint8(0*Reg1.*L2n1a.*Cell_bnd),I(:,:,3)) -handles.It(:,:,1)- .2*handles.It(:,:,2) ;
    % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));

    h = figure(2);  clf; imshow(Iz); hold on;
% hold on; plot(bndrys2{1}(:,1),bndrys2{1}(:,2),'g');
Iz = imflip(Iz,1);
Io = imresize(Iz,.5);
figure(3); clf; imshow(Io);



fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Results/';
%imwrite(Io,[fout,'cc12_hb.tif']);
imwrite(Io,[fout,'cc11_hb_yonly.tif']);

%%
folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/';
edat =  'MP09_22C_y_hb34_data.mat'; % cc11
 % edat = 'MP09_22C_y_hb05_data.mat'; % cc10
 
 % edat = 'MP09_22C_y_hb49_data.mat'; % cc11
% edat = 'BAC09_30C_y_hb01_data.mat'; % cc
load([folder,edat]); 

       age = getage(H,cent);
       disp(age);

I = handles.It;  

C = [1.25,0,0;
    0,1.25,0;
    0,0,.8];

T = [.05,1;
    .1,1;
    .0,1];

f = [0,0];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);


L1n2 = L1 - L2 > 0 ;

8/length(pts1)
11/(length(cent)-length(pts1)+11)


[h,w] = size(I(:,:,1));
 Iz = uint8(zeros(h,w,3));
    Iz(:,:,1) = imadd(uint8(255*Reg1.*L1n2.*Cell_bnd.*Reg1),I(:,:,1)) ;
    Iz(:,:,2) =  imadd(uint8(255*Reg1.*L2n1a.*Cell_bnd),1*I(:,:,2)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Iz(:,:,3) =    imadd(uint8(255*Reg1.*L2n1a.*Cell_bnd),I(:,:,3)) -handles.It(:,:,1)- .2*handles.It(:,:,2) ;
    % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
Iz = imflip(Iz,1);
    h = figure(2);  clf; imshow(Iz); hold on;
% hold on; plot(bndrys2{1}(:,1),bndrys2{1}(:,2),'g');

Io = imresize(Iz,.5);
figure(3); clf; imshow(Io);



fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Results/';
% imwrite(Io,[fout,'cc11_hb_7min.tif']);

%%

folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/';
edat =  'MP09_22C_y_hb34_data.mat'; % cc11
 % edat = 'MP09_22C_y_hb05_data.mat'; % cc10
 
 % edat = 'MP09_22C_y_hb49_data.mat'; % cc11
% edat = 'BAC09_30C_y_hb01_data.mat'; % cc
load([folder,edat]); 

       age = getage(H,cent);
       disp(age);

I = handles.It;  

C = [1.25,0,0;
    0,1.25,0;
    0,0,.8];

T = [.05,.6;
    .1,1;
    .0,1];

f = [0,0];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);





y = intersect(pts1,ptr_nucin1); 
notbord = union(ptr_nucin1,ptr_nucin2); 

err1 = 2;  err = 1.5; 
Nn = all_nucs(end);

err = zeros(1,Nn);
fin = zeros(1,Nn);
fou = zeros(1,Nn);

for t=1:Nn
err1 = err;
onreg = all_nucs(1:t);
offreg = all_nucs(t+1:end);



    fin(t) = length(intersect(onreg,y))/(t);
    fou(t) = length( intersect( intersect(offreg,notbord),setdiff(all_nucs,pts1)  )  )/ length(y);% (all_nucs(end)-t);
%                               (in the off region, not a border,  * nuclei not active)  / nuclei active    
%                                       number not active nuclei in region / number total active nuclei 
err(t) = fin(t)+fou(t);
end
[er,t] = min(err);


figure(1); clf; plot(err,'k'); hold on; plot(fou,'c'); plot(fin,'r');

figure(1); clf; plot(fliplr(err),'k','linewidth',3); hold on; 
plot(all_nucs,fliplr(fou),'c','linewidth',3); plot(fliplr(fin),'r','linewidth',3);
set(gcf,'color','w'); plot([Nn-t,Nn-t],[0,1],'k--','linewidth',1); ylim([0,max(fou)]);
xlabel('nuclei number'); ylabel('error fraction'); 

% Black background version of figure 1; 
figure(11);  clf; colordef black; plot(fliplr(err),'w','linewidth',3); hold on; 
 set(gca,'color','k');
plot(all_nucs,fliplr(fou),'c','linewidth',3); plot(fliplr(fin),'r','linewidth',3);
set(gcf,'color','k'); plot([Nn-t,Nn-t],[0,1],'w--','linewidth',1); ylim([0,max(fou)]);
xlabel('nuclei number'); ylabel('error fraction'); 


Roff = ismember(H,(1:t));
Ron = 1 - Roff; 

notBord = Reg1 | Reg2; 
[h,w] = size(I(:,:,1));
 Iz = uint8(zeros(h,w,3));
    Iz(:,:,1) = imadd(uint8(255*Roff.*L1.*Cell_bnd.*notBord),I(:,:,1)) ;
    Iz(:,:,2) =  imadd(uint8(255*Ron.*(1-L1).*Cell_bnd.*notBord),0*I(:,:,3)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Iz(:,:,3) =    imadd(uint8(255*Ron.*Cell_bnd.*notBord),I(:,:,3)) -handles.It(:,:,1)- .2*handles.It(:,:,2) ;
 
    Iz = imflip(Iz,1);
    
    figure(2);  clf; imshow(Iz); hold on;

Io = imresize(Iz,.5);
figure(3); clf; imshow(Io);
fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Results/';
% imwrite(Io,[fout,'cc11_hb_7min.tif']);


%%
