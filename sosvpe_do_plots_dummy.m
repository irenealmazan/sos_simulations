%pixels = [20 64];
%pixels = [64 20];
%pixels = [80 80];
pixels = [20 100];

figure;

subplot(141);
plot(ddam,squeeze(III_resize(pixels(1),pixels(2),:)),'ob'); 
hold on;
%plot(ddam,ones(size(ddam,1),1).*Ibar(pixels(1),pixels(2)),'LineWidth',3.0,'Color','r');
plot(ddam,squeeze(Ibar(pixels(1),pixels(2),:)),'LineWidth',3.0,'Color','r');
legend(['III @ (' num2str(pixels(1)) ',' num2str(pixels(2)) ')'],'Ibar');
xlabel('time (sim. units)');
ylabel('III');
title(['III vs time @ (' num2str(pixels(1)) ',' num2str(pixels(2)) ')']);
set(gca,'FontSize',15);

subplot(142);
if ~exist('dnlI','var')
    dlnI = III_resize./Ibar - 1;
end
plot(ddam,squeeze(dlnI(pixels(1),pixels(2),:)),'ob');
xlabel('time (sim. units)');
ylabel('dlnI');
title(['dnlI vs time at pixels = (' num2str(pixels(1)) ',' num2str(pixels(2)) ')']);
set(gca,'FontSize',15);

subplot(143);
plot(ddam,squeeze(Cdt(pixels(1),pixels(2),:)),'ob');
xlabel('time (sim. units)');
ylabel('Cdti');
title(['Cdti vs time at pixels = (' num2str(pixels(1)) ',' num2str(pixels(2)) ')']);
set(gca,'FontSize',15);

subplot(144);
III_integrated = squeeze(sum(sum(III_resize,1),2));
plot(ddam,III_integrated,'ob');

xlabel('time (sim. units)');
ylabel('III\_integrated');
title(['Integrated intensity']);
set(gca,'FontSize',15);

set(gcf,'Position',[50 433 1441 303]);


return;

%%%% plot the intensities



if ~exist('number_of_runs')
    number_of_runs = 1;
end

figure;
hold on;
for ss = 1:number_of_runs
   
    subplot(1,number_of_runs,ss)
    III_integrated = squeeze(sum(sum(damHW_struct(ss).III,1),2));
    plot(III_integrated,'o','Color',color_rates(ss));
    title([runname_titles(ss,:)]);
    
    if ss == 3
         ylim([0 0.11]);
    elseif ss == 4
           ylim([0 0.024]);
    else
        ylim([0 0.25]);
    
    end
    
    set(gca,'Yscale','log');
    
end

xlabel('time (sim. units)');
ylabel('III\_integrated');
title(['Integrated intensity']);
set(gca,'FontSize',15);

set(gcf,'Position',[50 433 1441 303]);