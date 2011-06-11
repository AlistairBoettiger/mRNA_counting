


    
clear all; 
N = 30; 

[X2,Y2,Z2] = meshgrid(1:N,1:N,1:N);

DNA = zeros(N,N,N);

dx = zeros(1,4*N); dy=zeros(1,4*N); 
dna2 = zeros(N); 

k = 2; 
dx(1) = .5*N; dy(1) = 1; 
while dx(k-1) < N && dy(k-1) < N

    dx(k) = round(dx(k-1) + .52*rand);
    dy(k) = round(dy(k-1) + .9*rand);
    dna2(dx(k),dy(k)) = 2; 
    dna2(dx(k)-1,dy(k)) = 1;
    dna2(dx(k)+1,dy(k)) = 1;
     k = k+1;
end
DNA(:,:,9) = dna2;
dx = nonzeros(dx); 
dy = nonzeros(dy); % collapse

   set(gcf,'color','w');
    axis off; 
    
    
    
    
   
   % DNA 
  figure(2); clf; axis off;
    patch(isosurface(Y2,X2,Z2,DNA,.5),'FaceColor','blue','EdgeColor','none'); hold on; 
    lighting phong; camlight left;
    xlim([0,N]); ylim([0,N]); zlim([0,N]);  axis off;
    
  
    % Promoter
        pr = [N/2,1,9]; % promoter location 
            [Xp,Yp,Zp] = sphere(10); Cp = ones(11); 
          surf(Xp+pr(1),Yp+pr(2),Zp+pr(3),Cp); 

    
% Proteins    
    [X,Y,Z] = sphere(20);  X = 2*X; Y=2*Y; Z = 2*Z;  
    Xo = X; Yo = Y; Zo=Z;
    
    C1 = 2*ones(21); 
    C2 = 2*ones(21); 
   
    
p1x = pr(1) + 3;  p1y = pr(2) + 3;  p1z =  pr(3) + 1; % starting centroid for molecule 1

P1x =  X + p1x; P1y = Y + p1y; P1z = Z + p1z;  % move 3D molecule 1 over centroid
P2x = X + 4; P2y = Y; P2z = Z; % move 3D-version of molecule 2 over the centroid. 

loaded = 0; % molecule 1 is not loaded.
tx_time = 100;  

for t = 1:1000;

    if loaded1 == 0; 
        % random diffusion for molecule 1
        s1x = round(2*(.5 - rand)); 
        s1y = round(2*(.5 - rand)); 
        s1z = round(2*(.5 - rand)); 

        p1x = p1x +s1x;  % keep track of centroid x 
        p1y = p1y +s1y;  % keep track of centroid
        p1z = p1z +s1z;  % keep track of centroid

        P1x = 1*(P1x + s1x);
        P1y = 1*(P1y + s1y);
        P1z = 1*(P1z + s1z);
     end 
     
    % If molecule gets close to the promoter, bind. 
    d = sqrt( (pr(1) - p1x).^2 + (pr(2) - p1y).^2  + (pr(3) - p1z).^2 );
    if  d < 9  && loaded1 == 0
      tx_time = 1;
      %  bind
      loaded1 = 1;  % record as loaded
      P1x = Xo + pr(1); % move to bound position
      P1y = Yo + pr(2);
      P1z = Zo + pr(3); 

      % plot bound molecule
        figure(1); clf; axis off;
        b1 = surf(P1x,P1y,P1z,C1); hold on;  set(b1,'facecolor','y'); 
        b2 = surf(P2x,P2y,P1z,C2); shading interp; set(b2,'facecolor','g'); 
        patch(isosurface(Y2,X2,Z2,DNA,.5),'FaceColor','blue','EdgeColor','none'); hold on; 
        lighting phong; camlight left;
        xlim([0,N]); ylim([0,N]); zlim([0,N]);  axis off;
    end    
      
    % Move forward; 
      if tx_time < length(dx); % still transcribing
          tx_time = tx_time + 1; % step forward in tx time 
          p1x = dx(tx_time) - dx(tx_time-1);
          p1y = dy(tx_time)- dy(tx_time-1);
                   
        P1x = 1*(P1x + p1x);
        P1y = 1*(P1y + p1y);
      else 
          loaded1 = 0; % finished transcription, no longer bound
      end
    
      
     if loaded2 == 0; 
        % random diffusion for molecule 1
        s2x = round(2*(.5 - rand)); 
        s2y = round(2*(.5 - rand)); 
        s2z = round(2*(.5 - rand)); 

        p2x = p2x +s2x;  % keep track of centroid x 
        p2y = p2y +s2y;  % keep track of centroid
        p2z = p2z +s2z;  % keep track of centroid

        P2x = 1*(P1x + s2x);
        P2y = 1*(P1y + s2y);
        P2z = 1*(P1z + s2z);
     end 
     
    % If molecule gets close to the promoter, bind. 
    d2 = sqrt( (pr(1) - p2x).^2 + (pr(2) - p2y).^2  + (pr(3) - p2z).^2 );
    if  d2 < 9  && loaded2 == 0
      tx_time = 1;
      %  bind
      loaded2 = 1;  % record as loaded
      P2x = Xo + pr(1); % move to bound position
      P2y = Yo + pr(2);
      P2z = Zo + pr(3); 
   end    
      
      
      
      
    figure(1); clf; axis off;
    b1 = surf(P1x,P1y,P1z,C1); hold on; set(b1,'facecolor','g'); 
    b2 = surf(P2x,P2y,P1z,C2); shading interp;  set(b1,'facecolor','r'); 
    patch(isosurface(Y2,X2,Z2,DNA,.5),'FaceColor','blue','EdgeColor','none'); hold on; 
    lighting phong; camlight left;
    xlim([0,N*.7]); ylim([0,N*.5]); zlim([0,N*.7]);  axis off;
    
 pause(.01);
end
