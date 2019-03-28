% Matlab code created 09-Aug-14 from fortran code of 1993
%     Also followed isflipNM2.m
%     Runs SOS growth simulation
%     Miscut surfaces
%     Standard kinetic barrier models for input of rate table:
%     Absolute, Relative, Eyring
%     Time-dependent growth rate / gas concentration
%     Brian Stephenson
%
% Passes many variables to change.m function as globals
%
% Number of directions for moves:
%     Directions 1-4 are jumps to the in-plane neighbors;
%     For rate matrix, direction 5 is "evaporation",
%     and direction 6 is "deposition".
%
% Transition rate matrix: rate(idir,irow,icol)
%     lists rate values for all possible moves
%     idir is direction, irow icol give (initial) location
%
% Number of states of an atom:
%     States 1-5 are surface atoms
%     at sites with 0-4 bonds to in-plane neighbors,
%     plus 1 bond to the atom beneath;
%     State 6 is an atom after evaporation or before deposition.
%
% Energy of each state
%     ewell(state) gives energy of each state, 
%     relative to unbound atom (state 6)
%
% Transition rate table: rtable(istate,fstate) 
%     gives transition rate for each type of move, 
%     from initial to final state.
%     For deposition, rtable(6,fstate) gives rate for unity surfconc;
%     actual deposition rates scaled by surfconc.
%     35 types total; rtable(6,6) not used.
%
% Rates are consistent with absolute rate theory description
%     where ewell(state) gives the energy of each state,
%     separated by barriers.
%
% Assumes a surface layer above crystal
%     with concentration surfconc, determined by mass balance 
%     between deposition, evaporation between layer and crystal
%     and addition to layer from vapor phase
%
% Equilibrium (fixed surfconc) growth rate (ML/s) = gasconc
% Equilibrium surfconc = (gasconc + rtable(ii,6)) / rtable(6,ii)

addpath(genpath('../'));

global kgrp nrow ncol nsteps neven ih ibi ibf ibfd sum3 sum2g sum2 sum1g sum1 sumd3 sumd2g sumd2 sumd1g rate rtable;
kstate = 6; % number of states


% Number of members in each group of rows or columns:
% This grouping added to increase speed of choosing a move.
kgrp = 16;

nrow = 128;      % Simulation size
ncol = 128;

%nrow = 32;       % Simulation size (test with small)
%ncol = 32;

nsteps = 8;      % number of miscut steps on initial surface
                 % should divide evenly into ncol/2

% Values to set rtable:
rmodel = 'R';     % IN-PLANE KINETIC MODEL TYPE (A, E, R), default E
xybond =  7;   % IN-PLANE BOND ENERGY / kT
zbond =   7;   % SURFACE NORMAL BOND ENERGY / kT
xybarr =  7;   % IN-PLANE BARRIER ENERGY / kT
zbarr =   7*[1 1 1 1 1];  % SURFACE NORMAL BARRIER ENERGIES / kT

amsurf = 1;      % NUMBER OF MONOLAYERS IN SURFACE PHASE

runname = 'miscut_8_16';

% Parameters that change at each interval:
gasconc =  [3.333e-9]/2;   % GAS CONCENTRATION
MLdata =   [0.005];       % ML INTERVAL BETWEEN DATA POINTS
dtdata =   MLdata./gasconc;   % TIME INTERVAL BETWEEN DATA POINTS
MLend =    [5];     % END ML FOR INTERVALS
tend =     MLend./gasconc;   % END TIMES FOR INTERVALS
saveflag = [1   ];   % FLAG: SAVE RESULTS
dispflag = [0   ];   % FLAG: DISPLAY RESULTS
%plotflag = [0   ];   % FLAG: DISABLES FIGURES TO RUN CODE WITH MATLAB OPTION -NODISPLAY
endimage = [0   ];   % FLAG: SAVE LAST IMAGE OF INTERVAL
freq_saving = 10; % FREQUENCY TO SAVE PARTIAL RESULTS

tstart = 0;

% Parameters for plotted images
clim = [0 2+nsteps+sum(gasconc.*tend)]; 
textclr = 'k'; 

% Calculate rate table from model type and input parameters:
% rtable(ii,jj) gives rate constant for transition from state ii to jj
% rtable(ii,6) are evaporations
% rtable(6,ii) are depositions

ewell = NaN*ones(kstate,1); % energies of each state
for ii = 1:kstate-1
    ewell(ii) = - (ii-1)*xybond - zbond;
end
ewell(kstate) = 0.;

rtable = NaN*ones(kstate,kstate);

% Set values with either index = kstate
for ii = 1:kstate-1    
    etrans = zbarr(ii) + ewell(kstate); % Always relative to max
    rtable(kstate,ii) = exp(ewell(kstate) - etrans);
    rtable(ii,kstate) = exp(ewell(ii) - etrans);
end
%rtable(kstate,kstate) = 0.; % Not used

% Set remaining values
for jj = 1:kstate-1
    for ii = 1:kstate-1
        if rmodel(1) == 'A'
            etrans = xybarr + ewell(1);
        elseif rmodel(1) == 'R'
            etrans = xybarr + max(ewell(ii),ewell(jj));  
        else
            etrans = xybarr + (ewell(ii)+ewell(jj))/2.;
        end
        rtable(ii,jj) = exp(ewell(ii) - etrans);
    end
end       

% Fix rtable to debug: 
if 0
rtable = ...
[1 1 1 1 1 0; ...
 1 1 1 1 1 0; ...
 1 1 1 1 1 0; ...
 1 1 1 1 1 0; ...
 1 1 1 1 1 0; ...
 0 0 0 0 0 0];

for ii = 1:kstate-1    % no evap or dep
    rtable(kstate,ii) = 0;
    rtable(ii,kstate) = 0;
end
end

disp('**  CALCULATED RATE TABLE:');
disp(rtable);

%% DO THE SIMULATION

sosvpe_mainloop;


