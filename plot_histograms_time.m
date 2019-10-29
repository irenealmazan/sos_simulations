% this is a crappy and lazy script to plot the summary of the statistical
% analysis of the speckle

delta_III = (III-Ibar)./Ibar;

figure(101);
clf



time_frame = 100;

col = [20 65]; % col(1) = column number, col(2) = time frame

subplot(3,3,1);plot(squeeze(III(:,col(1),col(2))),'b');
hold on;
plot(Ibar(:,col(1)),'r');
title(['col = ' num2str(col)]);
legend('III','Ibar');

subplot(3,3,4);
plot(squeeze(delta_III(:,col(1),col(2))),'k');
title('III - Ibar');

subplot(3,3,7);
histogram(delta_III(:,col(1),col(2)));
title('III - Ibar');


col = [44 time_frame]; % col(1) = column number, col(2) = time frame



subplot(3,3,2);plot(squeeze(III(:,col(1),col(2))),'b');
hold on;
plot(Ibar(:,col(1)),'r');
title(['col = ' num2str(col(1)) 'time frame =' num2str(col(2))]);

subplot(3,3,5);
plot(squeeze(delta_III(:,col(1),col(2))),'k');
title('III - Ibar');

subplot(3,3,8);
histogram(delta_III(:,col(1),col(2)));
title('III - Ibar');

col = [64 time_frame]; % col(1) = column number, col(2) = time frame



subplot(3,3,3);plot(squeeze(III(:,col(1),col(2))),'b');
hold on;
plot(Ibar(:,col(1)),'r');
title(['col = ' num2str(col(1)) 'time frame =' num2str(col(2))]);

subplot(3,3,6);
plot(squeeze(delta_III(:,col(1),col(2))),'k');
title('III - Ibar');

subplot(3,3,9);
histogram(delta_III(:,col(1),col(2)));
title('III - Ibar');

set(gcf,'Position',[69    47   890   704]);
