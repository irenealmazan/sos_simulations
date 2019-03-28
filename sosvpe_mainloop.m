% Matlab code created 09-Aug-14 from fortran code of 1993
% Modifid by I. Calvo-Almazan, 28-Mar-19
%     It has to be called through sosvpe.m which sets all the initialize
%     flags

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








% Initialize state variables
amono = 0.;
time = tstart;
nmove = 0;
neven = 0;
natoms = 0;
ndep = 0;
nevap = 0;
%cai = dmimage*nrow*ncol;
%cti = dtimage;
%ifm = 0;
%ift = 0;
%idat = 0;
nerr = zeros(1,3);
ntype = zeros(kstate,kstate);

ipt = 1;        % Image frame number
idt = 1;        % Saved data frame number 

intvl = 1;      % interval counter
ctd = tstart + dtdata(intvl);
surfconc = (gasconc(intvl)+rtable(5,6))./rtable(6,5); % initialize close to equil value
surfcdel = 1./(amsurf*nrow*ncol); % Change in surfconc from addition of 1 atom

%  Initialize height matrix:
%  Can have integer values
%ih = ones(nrow,ncol,'uint8');
ih = ones(nrow,ncol); % height of crystal surface at each site

%  Initialize bond matrices
%  These can have values from 1 to 5 based on number of bonds at each site
ibi = 5*ones(nrow,ncol); % initial number of bonds (state)
ibf = ones(4,nrow,ncol); % final number of bonds (state) after move in each direction
ibfd = ones(nrow,ncol); % final number of bonds after deposition

%  Modify for miscut:

if nsteps > 0
   width = ncol/nsteps;
   for istep = 1:nsteps
       istcol = ncol - (istep - 0.5)*width;
       ibi(:,istcol-1) = 4;
       ibf(2,:,istcol) = 2;
       ibf(4,:,istcol) = 2;
       ibfd(:,istcol) = 2;
       ibf(3,:,istcol+1) = 2;
       ih(:,1:istcol-1) = ih(:,1:istcol-1) + 1;
       if mod(ih(1,1),2) == 0
           neven = neven + nrow*(istcol-1);
       else
           neven = neven - nrow*(istcol-1);
       end
   end
end

%  Initialize rate matrix:
%  Copy from rate table for each possible move, 
%  based on initial and final states (move type)

for icol = 1:ncol
    for irow = 1:nrow
        for idir = 1:4
            rate(idir,irow,icol) = rtable(ibi(irow,icol),ibf(idir,irow,icol));
        end
        rate(5,irow,icol) = rtable(ibi(irow,icol),6);
        rate(6,irow,icol) = rtable(6,ibfd(irow,icol));
    end
end

%  Initialize sums and partial sums of rates used to locate next move
%
% sum3: sum of rates of jumps and evaporations (directions 1-4 and 5) at all sites
% sum2g: sums of groups of kgrp columns of sites
% sum2: sums of columns of sites
% sum1g: sums of groups of kgrp rows in a column of sites
% sum1: sum of directions 1-5 at a site

% sumd3 sum of rates of depositions (direction 6) at all sites
% sumd2g: sums of groups of kgrp columns of sites
% sumd2: sums of columns of sites
% sumd1g: sums of groups of kgrp rows in a column of sites

sum3 = 0.;  
sumd3 = 0.;
for icolg = 1:ncol/kgrp
    sum2g(icolg) = 0.;
    sumd2g(icolg) = 0.;
    for icol = icolg*kgrp-kgrp+1 : icolg*kgrp
        sum2(icol) = 0.;
        sumd2(icol) = 0.;
        for irowg = 1:nrow/kgrp
            sum1g(irowg,icol) = 0.;
            sumd1g(irowg,icol) = 0.;
            for irow = irowg*kgrp-kgrp+1 : irowg*kgrp
                sum1(irow,icol) = 0.;
                for idir = 1:5
                    sum1(irow,icol)=sum1(irow,icol)+rate(idir,irow,icol);
                end
                sum1g(irowg,icol)=sum1g(irowg,icol)+sum1(irow,icol);
                sumd1g(irowg,icol)=sumd1g(irowg,icol)+rate(6,irow,icol);
            end
            sum2(icol) = sum2(icol) + sum1g(irowg,icol);
            sumd2(icol) = sumd2(icol) + sumd1g(irowg,icol);
        end
        sum2g(icolg) = sum2g(icolg) + sum2(icol);
        sumd2g(icolg) = sumd2g(icolg) + sumd2(icol);
    end
    sum3 = sum3 + sum2g(icolg);
    sumd3 = sumd3 + sumd2g(icolg);
end


%ihs = NaN*uint8(ones(nrow,ncol,length(endimage))); % To save end images
ihs = NaN*ones(nrow,ncol,length(endimage)); % To save end images

%        Initialize image for real-time display

if sum(dispflag) > 0

  
    figure;
    set(gcf,'Position',[600 200 500 400]);
    set(gcf,'PaperPosition',[1 1 5 4]);
    axes('Box','on');
    imagesc(ih,clim);
    axis image
    shading flat;
    set(gca,'Visible','off');
    
    xp = ncol*[0 0.6 0.6 0];
    yp = nrow*[0 0 0.1 0.1];
    patch(xp,yp,'w');
    xtx = ncol/40;
    ytx = nrow/20;
    ht = text(xtx,ytx,[num2str(amono,'%5.2f'),' ML, time = ',int2str(time)]);
    set(ht,'Color',textclr,'FontSize',12);
    colorbar;
    pause(0.3);
    hc = get(gca,'Children');
    
else
    %{
    figure;
    set(gcf,'Position',[600 200 500 400]);
    set(gcf,'PaperPosition',[1 1 5 4]);
    axes('Box','on');
    imagesc(ih,clim);
    axis image
    shading flat;
    set(gca,'Visible','off');
    %}
    xp = ncol*[0 0.6 0.6 0];
    yp = nrow*[0 0 0.1 0.1];
%    patch(xp,yp,'w');
    xtx = ncol/40;
    ytx = nrow/20;
    ht = text(xtx,ytx,[num2str(amono,'%5.2f'),' ML, time = ',int2str(time)]);
%    set(ht,'Color',textclr,'FontSize',12);
%    colorbar;
    pause(0.3);
%    hc = get(gca,'Children');
    
end

%if dispflag(intvl) > 0

    % Save image matrix
    ihm(:,:,ipt) = ih;
    ipt = ipt + 1;
%end

if sum(saveflag) > 0
    amono=natoms/(nrow*ncol);
    damono = amono;
    dtime = time;
    dnmove = nmove;
    dneven = neven;
    dsurfc = surfconc;
    dnerr = nerr;
    dntype(:,:,idt) = ntype;
    dsum3 = sum3;
    dsumd3 = sumd3;
    idt = idt + 1;
end

%*******START LOOP FOR EACH INTERVAL ***********************

while(1)

%  For each gas conc interval, set dep rate multiplier:
%  rate per unit time of increase of surfconc
%  at steady-state, gasconc will be growth rate (ML/s)
surfcrate = gasconc(intvl)*surfcdel*nrow*ncol; 

%*******START LOOP FOR EACH MOVE****************************

while(1)

%  Calculate sumtot, sum of rates of all moves
%  If surfconc is negative, do not include deposition moves in sumtot:

sumtot = sum3 + sumd3*max(surfconc,0);

%  Pick move at random
%  Get position of move in rate matrix

pos = rand * sumtot;

%  If pos exceeds sum3, move is a deposition (see else below):

if pos <= sum3

%       Move is an in-plane jump or evaporation.

%       Find column group where sum2g first exceeds or equals pos,
%       move lies in that column group:

    psum = -pos;
    for icolg = 1:ncol/kgrp
        psum = psum + sum2g(icolg);
        if (psum >= 0.); break; end
    end
    if (psum < 0.)
        nerr(3) = nerr(3) + 1;
        disp(['Sum3 error: psum = ' num2str(psum)]);
        return
    end

%       Subtract from group sum to get position in group:

    pos = sum2g(icolg) - psum;

%       Find column where sum2 first exceeds or equals pos,
%       move lies in that column:

    psum = -pos;
    for icol = icolg*kgrp-kgrp+1 : icolg*kgrp
        psum = psum + sum2(icol);
        if (psum >= 0.); break; end
    end
    if (psum < 0)
        nerr(3) = nerr(3) + 1;
        disp(['Sum2g error: psum = ' num2str(psum)]);
        return
    end

%       Subtract from column sum to get position in column:

    pos = sum2(icol) - psum;

%       Likewise, locate site in column:

    psum = -pos;
    for irowg = 1 : nrow/kgrp
        psum = psum + sum1g(irowg,icol);
        if (psum >= 0.); break; end
    end
    if (psum < 0)
        nerr(2) = nerr(2) + 1;
        disp(['Sum2 error: icol, psum = ' num2str(icol) ' ' num2str(psum)]);
        return
    end

    pos = sum1g(irowg,icol) - psum;

    psum = -pos;
    for irow = irowg*kgrp-kgrp+1 : irowg*kgrp
        psum = psum + sum1(irow,icol);
        if (psum >= 0.); break; end
    end
    if (psum < 0)
        nerr(2) = nerr(2) + 1;
        disp(['Sum1g error: icol, psum = ' num2str(icol) ' ' num2str(psum)]);
        return
    end

    pos = sum1(irow,icol) - psum;

%       Choose move at site

    psum = -pos;
    for idir = 1:5
        psum = psum + rate(idir,irow,icol);
        if (psum >= 0.); break; end
    end
    if (psum < 0)
        nerr(1) = nerr(1) + 1;
        disp(['Sum1 error: irow, icol, psum = ' num2str(irow) ' ' num2str(icol) ' ' num2str(psum)]);
        return
    end

    if idir == 5 % Move is evaporation:
        natoms = natoms - 1; % Keep track of total atoms in crystal
        nevap = nevap + 1; % Keep track of evaporations
        surfconc = surfconc + surfcdel; % Evaporation adds an atom to the surface layer
        ntype(ibi(irow,icol),6) = ntype(ibi(irow,icol),6) + 1; % Keep track of total moves of each type
        change(irow,icol,-1); % Call routine to change height, bond, rate, sum matrices

    else  % Move is a jump

% Keep track of total moves of each type        
        ntype(ibi(irow,icol),ibf(idir,irow,icol)) = ntype(ibi(irow,icol),ibf(idir,irow,icol)) + 1;

% To debug, report chosen move type
% disp([nmove idir, irow, icol, ibi(irow,icol), ibf(idir,irow,icol)]);

% Remove atom from initial site
        change(irow,icol,-1); % Call routine to change height, bond, rate, sum matrices

%       Calculate final position of move
        if idir == 1
            icol = icol + 1;
        elseif idir == 2 
            irow = irow + 1;
        elseif idir == 3 
            icol = icol - 1;
        elseif idir == 4
            irow = irow - 1;
        end
%       Apply periodic boundary conditions
        if icol <= 0 
            icol = ncol;
        elseif icol > ncol 
            icol = 1;
        end
        if irow <= 0 
            irow = nrow;
        elseif irow > nrow 
            irow = 1;
        end

% Add atom to final site
        change(irow,icol,1); % Call routine to change height, bond, rate, sum matrices

    end

else    % Move is a deposition

%       Find position for deposition:

    pos = (pos - sum3)/surfconc;

%       Find column group where sumd2g first exceeds or equals pos,
%       move lies in that column group:

    psum = -pos;
    for icolg = 1:ncol/kgrp
        psum = psum + sumd2g(icolg);
        if (psum >= 0.); break; end
    end
    if (psum < 0.)
        nerr(3) = nerr(3) + 1;
        disp(['Sum3 error: psum = ' num2str(psum)]);
        return
    end

%       Subtract from group sum to get position in group:

    pos = sumd2g(icolg) - psum;

%       Find column where sumd2 first exceeds or equals pos,
%       move lies in that column:

    psum = -pos;
    for icol = icolg*kgrp-kgrp+1 : icolg*kgrp
        psum = psum + sumd2(icol);
        if (psum >= 0.); break; end
    end
    if (psum < 0)
        nerr(3) = nerr(3) + 1;
        disp(['Sumd2g error: psum = ' num2str(psum)]);
        return
    end

%       Subtract from column sum to get position in column:

    pos = sumd2(icol) - psum;

%       Likewise, locate site in column:

    psum = -pos;
    for irowg = 1 : nrow/kgrp
        psum = psum + sumd1g(irowg,icol);
        if (psum >= 0.); break; end
    end
    if (psum < 0)
        nerr(2) = nerr(2) + 1;
        disp(['Sumd2 error: icol, psum = ' num2str(icol) ' ' num2str(psum)]);
        return
    end

    pos = sumd1g(irowg,icol) - psum;

    psum = -pos;
    for irow = irowg*kgrp-kgrp+1 : irowg*kgrp
        psum = psum + rate(6,irow,icol);
        if (psum >= 0.); break; end
    end
    if (psum < 0)
        nerr(2) = nerr(2) + 1;
        disp(['Sumd1g error: icol, psum = ' num2str(icol) ' ' num2str(psum)]);
        return
    end

%       Make deposition move:

    natoms = natoms + 1; % Keep track of total atoms in crystal
    ndep = ndep + 1; % Keep track of depositions
    surfconc = surfconc - surfcdel; % Deposition removes atom from surface layer
    ntype(6,ibfd(irow,icol)) = ntype(6,ibfd(irow,icol)) + 1; % Keep track of total moves of each type
    change(irow,icol,1); % Call routine to change height, bond, rate, sum matrices

end

% Finished with whatever type of move, now increment time step

time = time + 1./sumtot;
surfconc = surfconc + surfcrate*(1./sumtot); % Add atoms to layer from gas
nmove = nmove + 1; % Keep track of total moves

% if nmove > 25; break; end

%       Write to output if it is time:

if time > ctd
    
    ctd = ctd + dtdata(intvl);
    
    if saveflag(intvl) > 0
        amono=natoms/(nrow*ncol);
        damono = [damono; amono];
        dtime = [dtime; time];
        dnmove = [dnmove; nmove];
        dneven = [dneven; neven];
        dsurfc = [dsurfc; surfconc];
        dnerr = [dnerr; nerr];
        dntype(:,:,idt) = ntype;
        dsum3 = [dsum3; sum3];
        dsumd3 = [dsumd3; sumd3];
        idt = idt + 1;
    end

%       Update screen
    if dispflag(intvl) > 0
        
        set(hc(end),'CData',ih);
        %imagesc(ih,clim);
        %axis image;
        %set(gca,'Visible','off');
        patch(xp,yp,'w');
        ht = text(xtx,ytx,[num2str(amono,'%5.2f'),' ML, time = ',int2str(time)]);
        set(ht,'Color',textclr,'FontSize',12);
        pause(0.05);
    
        % Save image matrix
        ihm(:,:,ipt) = ih;
        ipt = ipt + 1;
        
    else
        %set(hc(end),'CData',ih);
        %imagesc(ih,clim);
        %axis image;
        %set(gca,'Visible','off');
        %patch(xp,yp,'w');
        ht = text(xtx,ytx,[num2str(amono,'%5.2f'),' ML, time = ',int2str(time)]);
        %set(ht,'Color',textclr,'FontSize',12);
        pause(0.05);
    
        % Save image matrix
        ihm(:,:,ipt) = ih;
        ipt = ipt + 1;
        
    end

%       Display stats

    disp(['nmove = ' num2str(nmove) ', time = ' num2str(time) ', neven = ' num2str(neven) ', amono = ' num2str(amono) ', surfc = ' num2str(surfconc)]);
end

if (time >= tend(intvl)); break; end
%*******END LOOP FOR EACH MOVE****************************
end

if endimage(intvl)
   ihs(:,:,intvl) = ih;
   eamono(intvl) = natoms/(nrow*ncol);
   etime(intvl) = time;
end
    
    
intvl = intvl + 1;
if (intvl > length(tend)); break; end
ctd = tend(intvl-1) + dtdata(intvl);
%*******END LOOP FOR EACH INTERVAL****************************
end

if sum(endimage) > 0
    disp('Printing end images.');
    for ii = 1:length(endimage)
        if endimage(ii)
           
            if sum(dispflag)
            figure;
            set(gcf,'Position',[600 200 400 400]);
            set(gcf,'PaperPosition',[1 1 4 4]);
            axes('Box','on');
            %colormap('hot');
            imagesc(ihs(:,:,ii),clim);
            axis image
            shading flat;
            set(gca,'Visible','off');
            end
            
            if 1
                xp = ncol*[0 0.6 0.6 0];
                yp = nrow*[0 0 0.1 0.1];
                patch(xp,yp,'w');
                xtx = ncol/40;
                ytx = nrow/20;
                ht = text(xtx,ytx,[num2str(eamono(ii),'%5.2f'),' ML, time = ',int2str(etime(ii))]);
                if sum(dispflag)
                    set(ht,'Color',textclr,'FontSize',12);
                end
            end
            
%            eval(['print -depsc ' runname '_endimage_' num2str(ii)]);
%            eval(['print -djpeg90 ' runname '_endimage_' num2str(ii)]);
        end
    end
end


if sum(saveflag) > 0
    disp('Writing statistics file.');
    eval(['save ' runname '_stats nrow ncol nsteps rmodel xybond zbond xybarr zbarr amsurf gasconc dtdata tstart tend ewell rtable damono dtime dnmove dneven dsurfc dnerr dntype dsum3 dsumd3 ihm']);
end

if sum(dispflag)
    figure
    plot(damono,(2*dneven/(nrow*ncol)-1).^2);
    xlabel('Monolayers');
    ylabel('Anti-Bragg Intensity');
end