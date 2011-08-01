
%% bistability model for snail

% Alistair Boettiger
clear all;

R = linspace(0,100,5000);
S = linspace(0,400,5000);


nS = 10; KmS = 60;   VS = 100;  Ts = 1.8;
nR = 7; KmR = 55;   VR = 100;
Sm = VS*KmS^nS./(R.^nS + KmS^nS);
Rm = VR*KmR^nR./(S.^nR + KmR^nR);

figure(1); clf; plot(R,Ts*Sm,'linewidth',3);
xlabel('[Sna-Repressors]','FontSize',14);
ylabel('steady state [sna-mRNA]','FontSize',14);
set(gca,'FontSize',14); 
xlim([0,100]);
ylim([0,250]);
figure(2); clf; plot(S,Rm,'linewidth',3);
xlabel('[Sna]','FontSize',14); 
ylabel('steady state [repressor]','FontSize',14); 
set(gca,'FontSize',14);
xlim([0,250]);



%%

figure(1); clf; set(gcf,'color','w');
plot(S,Rm,'g','linewidth',3); hold on; 
plot(Ts*Sm,R,'r','linewidth',3); 
xlabel('[sna mRNA]','FontSize',14); 
ylabel('repressors of sna','FontSize',14);
set(gca,'FontSize',14); 


figure(1); 
plot(Ts*2*Sm,R,'r--','linewidth',2); 
plot(Ts*.5*Sm,R,'r-.','linewidth',2); 

[X,Y] = meshgrid(linspace(2,400,20),linspace(.02,125,20));
V = VR*KmR^nR./(X.^nR + KmR^nR) - Y;
U = (VS*KmS^nS./(Y.^nS + KmS^nS) - X/Ts);
quiver(X,Y,U,V,'MaxHeadSize',.1);
xlim([0,400]); ylim([0,125]);


Sm = VS*KmS^nS./(Rm.^nS + KmS^nS);
figure(2); clf;  set(gcf,'color','w');
plot(S,Sm,'b','linewidth',3); hold on;
plot(S,S/Ts,'r','linewidth',3');
plot(S,2*Sm,'b--','linewidth',3); hold on;
plot(S,.5*Sm,'b-.','linewidth',3); hold on;
set(gca,'FontSize',14);
xlabel('[Sna]','FontSize',14); 
ylabel('dS/dt','FontSize',14); 
xlim([0,400]);

%%
 Sm = VS*KmS^nS./(Rm.^nS + KmS^nS);
figure(3); clf; plot(S,Ts*Sm,'linewidth',3);  
xlabel('[Sna]','FontSize',14); 
ylabel('steady state [sna-mRNA]','FontSize',14);
set(gca,'FontSize',14); 
ylim([0,250]); xlim([0,250]);

%%
Ss = Sm-S/Ts;
Ss_p = Ss;  Ss_p(Ss<-.02) = NaN;
Ss_n = Ss;  Ss_n(Ss>.02) = NaN;

figure(3); clf; set(gcf,'color','w'); 
plot(S,Ss_p ,'r-','linewidth',3); hold on;
plot(S,Ss_n ,'b-','linewidth',3);
xlim([0,max(S)]); 
ylim([-30,100]);
xlabel('[Sna]','FontSize',14); 
ylabel('d[Sna]/dt','FontSize',14);
set(gca,'FontSize',14);

[Sg,Yg] = meshgrid( S(1:200:end),linspace(-20,100,25));
[u,v] = meshgrid(Ss(1:200:end),zeros(1,25));

quiver(Sg,Yg,u,v);

% quiver(S(1:250:end),0*ones(1,20),Ss(1:250:end),zeros(1,20));

% % Nullclines
% figure(1); clf; set(gcf,'color','w'); 
% plot(S,Sm ,'r-','linewidth',3); hold on;
% plot(S,S/Ts,'b-','linewidth',3);
% xlim([0,max(S)]); 
% ylim([0,150]);
% xlabel('[Sna]','FontSize',14); 
% ylabel('d[Sna]/dt','FontSize',14);
% set(gca,'FontSize',14);



%% Represor Model

 KmSS = 150; nSS = 20; % strong
  % KmSS = 110; nSS = 2;

 Ts = 5; VS = 100;

Sm = KmSS.^nSS*(VS*KmS^nS./(Rm.^nS + KmS^nS))./(KmSS.^nSS  + S.^nSS); 

S2_syn = KmSS.^nSS*(VS*KmS^nS./(R.^nS + KmS^nS))./(KmSS.^nSS  + S.^nSS); 
figure(3); clf; plot(R,S2_syn,'linewidth',3);
set(gca,'FontSize',14);
xlabel('[sna Repressors]','FontSize',14); 
ylabel('dS/dt','FontSize',14); 

figure(3); clf; plot(S,Sm,'linewidth',3); hold on;
plot(S,S/Ts,'r','linewidth',3);


figure(3);  
Sm = KmSS.^nSS*(2*VS*KmS^nS./(Rm.^nS + KmS^nS))./(KmSS.^nSS  + S.^nSS); 
plot(S,Sm,'b--','linewidth',3); hold on;
Sm = KmSS.^nSS*(.5*VS*KmS^nS./(Rm.^nS + KmS^nS))./(KmSS.^nSS  + S.^nSS); 
plot(S,Sm,'b--','linewidth',3); hold on;
set(gca,'FontSize',14);
xlabel('[Sna]','FontSize',14); 
ylabel('dS/dt','FontSize',14); 
xlim([0,max(S)]); 
legend('sna synthesis','sna degredation', '2x sna dose','1/2x sna dose');


% 
% Kiss = zeros(length(R)/10,nSS+1);
% for r = 1:max(R);
%  Kiss(r,:)  =  roots([1,zeros(1,nSS-1),KmSS.^nSS,-Ts*VS*KmSS.^nSS*KmS.^nS./(KmS.^nS + (r-1).^nS)] )';
% end
% figure(5); clf;  set(gcf,'color','w'); 
% plot(real(Kiss(:,end)),'linewidth',3); 
% xlabel('[R]','FontSize',14);  
% ylabel('steady state [Sna]','FontSize',14); 
% set(gca,'FontSize',14);


Kiss = zeros(length(R)/10,nSS+1);
for r = 1:length(R)/10;
 Kiss(r,:)  =  roots([1,zeros(1,nSS-1),KmSS.^nSS,-Ts*VS*KmSS.^nSS*KmS.^nS/(KmS.^nS + Rm(10*r).^nS)] )';
end
figure(5); clf;  set(gcf,'color','w'); 
plot(real(Kiss(:,end)),fliplr(R(1:10:end)),'r','linewidth',3); hold on;
plot(S,Rm,'g','linewidth',3); 

[X,Y] = meshgrid(linspace(2,250,20),linspace(.02,125,20));
V = VR*KmR^nR./(X.^nR + KmR^nR) - Y;
U = (KmSS^nSS*(VS*KmS^nS./(Y.^nS + KmS^nS))./(KmSS^nSS+X.^nSS) - X/Ts);
quiver(X,Y,U,V,'MaxHeadSize',.1);
xlim([0,250]); ylim([0,125]);




VS = 2*VS;
Kiss = zeros(length(R)/10,nSS+1);
for r = 1:length(R)/10;
 Kiss(r,:)  =  roots([1,zeros(1,nSS-1),KmSS.^nSS,-Ts*VS*KmSS.^nSS*KmS.^nS/(KmS.^nS + Rm(10*r).^nS)] )';
end
figure(5); hold on;
plot(real(Kiss(:,end)),fliplr(R(1:10:end)),'r--','linewidth',2); hold on;
% legend('sna nullcline','sna-rep nullcline', '2x sna copy','Location','SouthOutside');
xlabel('[S]','FontSize',14);  
ylabel('[R]','FontSize',14); 
set(gca,'FontSize',14);


VS = VS/4;
Kiss = zeros(length(R)/10,nSS+1);
for r = 1:length(R)/10;
 Kiss(r,:)  =  roots([1,zeros(1,nSS-1),KmSS.^nSS,-Ts*VS*KmSS.^nSS*KmS.^nS/(KmS.^nS + Rm(10*r).^nS)] )';
end
figure(5); hold on;
plot(real(Kiss(:,end)),fliplr(R(1:10:end)),'r-.','linewidth',2); hold on;
% legend('sna nullcline','sna-rep nullcline', '2x sna copy','Location','SouthOutside');
xlabel('[S]','FontSize',14);  
ylabel('[R]','FontSize',14); 
set(gca,'FontSize',14);





%%

Ss = Sm-S/Ts;
Ss_p = Ss;  Ss_p(Ss<-.02) = NaN;
Ss_n = Ss;  Ss_n(Ss>.02) = NaN;

figure(4); clf; set(gcf,'color','w'); 
plot(S,Ss_p ,'r-','linewidth',3); hold on;
plot(S,Ss_n ,'b-','linewidth',3);
xlim([0,max(S)]); 
ylim([-30,100]);
xlabel('[Sna]','FontSize',14); 
ylabel('d[Sna]/dt','FontSize',14);
set(gca,'FontSize',14);

[Sg,Yg] = meshgrid( S(1:200:end),linspace(-20,100,25));
[u,v] = meshgrid(Ss(1:200:end),zeros(1,25));
quiver(Sg,Yg,u,v);




