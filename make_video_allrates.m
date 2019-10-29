% this scripts makes the video of the simulations at different growth rates


% growth rate 3.3e-9
path = 'miscut_32_29/';%'miscut_32_28/';
runname = [path 'miscut_32_29'];%[path 'miscut_32_28'];
make_vid(runname,[0 35],'jet');


% growth rate 3.3e-9/2
path = 'miscut_32_30/';%'miscut_32_28/';
runname = [path 'miscut_32_30'];%[path 'miscut_32_28'];
make_vid(runname,[0 35],'jet');

% growth rate [3.3e-9 for 2 ML 0]
path = 'miscut_32_31/';%'miscut_32_28/';
runname = [path 'miscut_32_31'];%[path 'miscut_32_28'];
make_vid(runname,[0 35],'jet');