open('Fiber_WS_Field_Data.fig');
lh = findall(gca, 'type', 'line');% ??????????lh?????
xc = get(lh, 'xdata');            % ??x????xc???????
yc = get(lh, 'ydata');            % ??y????yc???????
title(' Curved Histogram of Fiber Settling Velocity','fontsize',18)
xlabel('${w_s (cm/s)}$','fontsize',36,'Interpreter','latex') 
ylabel('count','fontsize',36,'Interpreter','latex')
legend({'flat fiber','cylindrical fiber'},'Location','northeast','fontsize',18)