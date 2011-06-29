  
% Needs work, scrap 0

x = 0:10;
y = sin(x);
xx = 0:.25:10;
yy = spline(x,y,xx);
    figure(5); clf; plot(x,y,'o',xx,yy)
    %%
    figure(4); clf; colordef white; set(gcf,'color','w');
 
    
 x = linspace(0,3,100);
 a = 5*sin(x) + 2*randn(1,100);
 b = 5*sin(x-1) + 2*randn(1,100);
 [jnk,shift] = max(xcorr(a,b))
 figure(5); clf; plot(x,a); hold on; plot(x,b,'r');
 hold on; plot(x,a,'b.'); 
 plot(x,[b(end-shift+1:end),b(1:end-shift)],'r.');
    
 
 
 
    for e = 1:14 ; % e =5;
        try
 
            mcor = MP09{e}.Nnucs./MP09{1}.Nnucs;
            grad = mcor*MP09{e}.Data_sort(:,3);
            stnd = MP09{1}.Data_sort(:,3);
            gx = MP09{e}.Data_sort(:,1);
            sx = MP09{1}.Data_sort(:,1);            
            
            L = min(length(grad),length(stnd));
            [jnk,shift] = max(xcorr(stnd(1:L),grad(1:L)))
            
            
            
%grad = flipud(grad);
figure(5); clf; plot(gx,flipud(grad),'r.'); hold on;
figure(5); plot(sx,flipud(stnd),'b.');

xdata = linspace(min(gx(1),sx(1)),max(gx(end),sx(end)),1000)';

n = 3; theta = mean(gx);  A = max(grad);  b=min(grad); 
[p,fit] = fxn_fit_sigmoid(gx',flipud(grad)',[n,theta,A,b]);
grad = p(3).*(xdata).^p(1)./(p(2).^p(1)+(xdata).^p(1)) + p(4) ; 

n = 3; theta = mean(sx);  A = max(stnd);  b=min(stnd); 
[p,fit] = fxn_fit_sigmoid(sx',flipud(stnd)',[n,theta,A,b]);
stnd = p(3).*(xdata).^p(1)./(p(2).^p(1)+(xdata).^p(1)) + p(4) ; 

figure(5); clf; plot(xdata,flipud(grad),'r.'); hold on;
figure(5); plot(xdata,flipud(stnd),'b.');


[jnk,shift] = max(xcorr(fliplr(grad),fliplr(stnd)))


figure(5); hold on; plot(xdata,stnd,'b'); hold on; plot(xdata,grad,'r');
            
            stnd = [zeros(3E2,1); stnd; zeros(3E2,1)]; % room to move around
            L = length(stnd);
            grad = [grad; zeros(L-length(grad),1)];
                      
              err2 = 0;  
              err1 = sum(abs(stnd-grad));
            while err2 <= err1 
                    err1 = sum(abs(stnd-grad));
                    grad = [grad(L-1:L);grad(1:L-2)];
                    err2 = sum(abs(stnd-grad));
            end
            
            figure(5); clf; plot(grad,'b.'); hold on; plot(stnd,'r.'); 
       
            
            
             plot(MP09{e}.Data_sort(:,1),MP09{e}.Data_sort(:,3),'g.'); hold on;
        catch
            continue
        end
    end
    
   %%
%               mcor = MP09{e}.Nnucs./MP09{1}.Nnucs;
%             grad = mcor*MP09{e}.Data_sort(:,3);
%             stnd = MP09{1}.Data_sort(:,3);
%             gx = MP09{e}.Data_sort(:,1);
%             sx = MP09{1}.Data_sort(:,1);
%             
%             st = min(gx(1),sx(1)); %   gx(1) ~= sx(1)~ = 0;
%             en = max(gx(end),sx(end)); 
%             
%             figure(5); clf; hold on; plot(gx,grad,'m.');
%             plot(sx,stnd,'c.');
%             xi = linspace(st,en,100);
% 
%             grad = smooth(grad,.3,'rloess');
%             stnd = smooth(stnd,.3,'rloess');
%            % figure(5); plot(gx,grad,'b.');
%             grad = interp1(gx,grad,xi);
%             stnd = interp1(sx,stnd,xi);
%             figure(5); hold on; plot(xi,grad,'r.');
%           figure(5); hold on; plot(xi,stnd,'b.');
%           
%             
%            
%             offset = -1E5;
%             gx = gxi
%         
%       
%             err1 = 1; err2 = 0; ostp = 100;
%             while err2<=err1
%                 err1 = sum((grad - stnd).^2 + (gxi-sx).^2) ;
%                 offset = offset +ostp;
%                 gx = gx + ostp;
%                 err2 = sum((grad - stnd).^2 + (gxi-sx).^2) ;
%                 
%                 figure(5); clf; plot(sx,stnd,'r');
%                 hold on; plot(gx,grad,'b');
%                 
%             end
%             figure(5); clf; 
%             plot(MP09{1}.Data_sort(:,1),MP09{1}.Data_sort(:,3),'b'); hold on;
%             plot(MP09{e}.Data_sort(:,1)+offset,MP09{e}.Data_sort(:,3),'r'); 
%     
    