
%%  'Gaussian' Sphere


     set(gcf,'color','w');
    
    [X,Y,Z] = sphere(20);
    C = 2*ones(21); 
    
    figure(1); clf; 
    r = linspace(1,5,5);
        N = length(r); 
    for j= 1:N
        figure(1); 
         h = surf(r(j)*X,r(j)*Y,r(j)*Z,C);
         alpha .4; % shading interp; 
         hold on;
    end
    
    [X1,Y1] = meshgrid(linspace(-8,8,21));
    Z1 = 3*ones(21); 
    surf(X1,Y1,Z1,.9*C); % shading interp;  
    light('Position',[0,0,0],'Style','local'); 
    %    camlight left;
    lighting phong; axis off;