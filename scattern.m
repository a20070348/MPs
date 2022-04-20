clear;
clc;

p2pz10n = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10p2pn.dat');
p3pz10n = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10p3pn.dat');
p4pz10n = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10p4pn.dat');
p5pz10n = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10p5pn.dat');
wpz10n = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10wpn.dat');

subfig1 = subplot(2,2,1);
binscatter(p2pz10n,wpz10n,[100,100]);
colormap(subfig1,'jet');
set(gca,'FontSize',24);
title('scatter plot for horizontally aligned fiber','fontsize',18,'position',[100,6.5])
xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
xlim([-100 300])

subfig2 = subplot(2,2,2);
binscatter(p3pz10n,wpz10n,[100,100]);
colormap(subfig2,'jet');
set(gca,'FontSize',24);
title('scatter plot for randomly aligned fiber','fontsize',18,'position',[100,6.5])
xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
xlim([-100 300])

subfig3 = subplot(2,2,3);
binscatter(p4pz10n,wpz10n,[100,100]);
colormap(subfig3,'jet');
set(gca,'FontSize',24);
title('scatter plot for fiber settling dependent on local averaged dissipate rate','fontsize',18,'position',[100,6.5])
xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
xlim([-100 300])

subfig4 = subplot(2,2,4);
binscatter(p5pz10n,wpz10n,[100,100]);
colormap(subfig4,'jet');
set(gca,'FontSize',24);
title('scatter plot for fiber settling dependent on horizontally averaged dissipate rate','fontsize',18,'position',[100,6.5])
xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
xlim([-100 300])
% 
% p2pz05cc = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz05p2pcc.dat');
% p3pz05cc = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz05p3pcc.dat');
% p4pz05cc = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz05p4pcc.dat');
% p5pz05cc = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz05p5pcc.dat');
% wpz05cc = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz05wpcc.dat');
% 
% subfig1 = subplot(2,2,1);
% binscatter(p2pz05cc,wpz05cc,[100,100]);
% colormap(subfig1,'jet');
% set(gca,'FontSize',24);
% title('scatter plot for horizontally aligned fiber','fontsize',18,'position',[100,4.5])
% xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
% ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
% xlim([-100 300])
% 
% subfig2 = subplot(2,2,2);
% binscatter(p3pz05cc,wpz05cc,[100,100]);
% colormap(subfig2,'jet');
% set(gca,'FontSize',24);
% title('scatter plot for randomly aligned fiber','fontsize',18,'position',[100,4.5])
% xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
% ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
% xlim([-100 300])
% 
% subfig3 = subplot(2,2,3);
% binscatter(p4pz05cc,wpz05cc,[100,100]);
% colormap(subfig3,'jet');
% set(gca,'FontSize',24);
% title('scatter plot for fiber settling dependent on local averaged dissipate rate','fontsize',18,'position',[100,4.5])
% xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
% ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
% xlim([-100 300])
% 
% subfig4 = subplot(2,2,4);
% binscatter(p5pz05cc,wpz05cc,[100,100]);
% colormap(subfig4,'jet');
% set(gca,'FontSize',24);
% title('scatter plot for fiber settling dependent on horizontally averaged dissipate rate','fontsize',18,'position',[100,4.5])
% xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
% ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
% xlim([-100 300])

% p2pz30cc = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz30p2pcc.dat');
% p3pz30cc = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz30p3pcc.dat');
% p4pz30cc = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz30p4pcc.dat');
% p5pz30cc = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz30p5pcc.dat');
% wpz30cc = load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz30wpcc.dat');
% 
% subfig1 = subplot(2,2,1);
% binscatter(p2pz30cc,wpz30cc,[100,100]);
% colormap(subfig1,'jet');
% set(gca,'FontSize',24);
% title('scatter plot for horizontally aligned fiber','fontsize',18,'position',[100,4.5])
% xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
% ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
% xlim([-100 300])
% 
% subfig2 = subplot(2,2,2);
% binscatter(p3pz30cc,wpz30cc,[100,100]);
% colormap(subfig2,'jet');
% set(gca,'FontSize',24);
% title('scatter plot for randomly aligned fiber','fontsize',18,'position',[100,4.5])
% xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
% ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
% xlim([-100 300])
% 
% subfig3 = subplot(2,2,3);
% binscatter(p4pz30cc,wpz30cc,[100,100]);
% colormap(subfig3,'jet');
% set(gca,'FontSize',24);
% title('scatter plot for fiber settling dependent on local averaged dissipate rate','fontsize',18,'position',[100,4.5])
% xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
% ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
% xlim([-100 300])
% 
% subfig4 = subplot(2,2,4);
% binscatter(p5pz30cc,wpz30cc,[100,100]);
% colormap(subfig4,'jet');
% set(gca,'FontSize',24);
% title('scatter plot for fiber settling dependent on horizontally averaged dissipate rate','fontsize',18,'position',[100,4.5])
% xlabel('$p^\prime$','fontsize',36,'Interpreter','latex') 
% ylabel('$w^\prime$','fontsize',36,'Interpreter','latex')
% xlim([-100 300])

