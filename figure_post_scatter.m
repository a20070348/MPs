open('ws_scatter.fig');
lh = findobj(gca,'type','scatter');% ??????????lh?????
xc = get(lh, 'xdata');            % ??x????xc???????
yc = get(lh, 'ydata');            % ??y????yc???????
zc = get(lh, 'zdata');
gscatter(xc,yc,zc,24);
title(' Curved Histogram of Fiber Settling Velocity','fontsize',18)
xlabel('${w_s (cm/s)}$','fontsize',36,'Interpreter','latex') 
ylabel('count','fontsize',36,'Interpreter','latex')
legend({'flat fiber','cylindrical fiber'},'Location','northeast','fontsize',18)