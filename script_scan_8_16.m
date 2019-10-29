 TCV = [800]; % Equilibrium Thermocouple (new heater)
                    
% power in Watts
PWV = [25];

% Filenames
specfilenameM = ['miscut_'];

% scans number (see lab book)
SCNstrM = ['32_18'];

% Center pixels in between the CTRS
%XCENV = [230;236;240;240];% del
%YCENV = [120;118;120;120]; % nu

XCENV = [64];% del
YCENV = [64]; % nu



% Width of the larger ROI
XWIDV = [60];
% Pilatus detector on sevchex arm (X is del and Y is nu)
YWIDV = [60];


% Offset of two of the ROIS around the Positions of CTRs
ymax = [8];

% Min and max on time range for delta-time average
MLdata = 0.005;
tminv = [1];
tmaxv = [6000];
%tminv = [720; 2500;1110;200];
%tmaxv = [1400;3500;2000;600];