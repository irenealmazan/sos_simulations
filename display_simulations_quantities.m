% with this script we display the quantites which govern the temporal
% evolution of the s-o-s simulations

figure(1);
clf;

subplot(211);
plot(dtime);
title('dtime where Delta time = 1/sumtot');

subplot(212);
plot(1./diff(dtime));
set(gca,'Yscale','log');
title('1/diff(dtime) = sumtot');

et(gcf,'Name','gasconc = 0');


%{
figure(1);
clf;

subplot(6,2,1);
plot(damono);
title('damono where amono=natoms/(nrow*ncol)');


subplot(6,2,2);
plot(diff(damono));
title('diff(damono)');


subplot(6,2,3);
plot(dtime);
title('dtime where time = 1/sumtot');

subplot(6,2,4);
plot(1./diff(dtime));
title('1/diff(dtime) = sumtot');

subplot(6,2,5);
plot(dnmove);
title('total moves dnmove');

% subplot(7,2,7);
% plot(dneven);
% title('neven ???');

subplot(6,2,7);
plot(dsurfc);
title('surf conc (atoms/a^2) (dsurfc)');

subplot(6,2,6);
plot((dsum3 + max(dsurfc,0).*dsumd3));
title('dsum3 + dsurfc.*dsumd3');

subplot(6,2,9);
plot(dsum3);
title('sum of jump rates (1-5) for all the sites (dsum3)');

subplot(6,2,11);
plot(dsumd3);
title('evaportation rates for all the sites (dsumd3)');
%}