function change(iy, ix, iflag)
% function change(iy, ix, iflag)
% add or remove atom at ix, iy, depending on sign of flag:
%   add atom if iflag >= 0
%   remove atom if iflag < 0
% and recalculate rate and bond tables
% called by sosvpe.m
% edited from subroutine in sosvpe.for
% GBS 09-Aug-14

global kgrp nrow ncol nsteps neven ih ibi ibf ibfd sum3 sum2g sum2 sum1g sum1 sumd3 sumd2g sumd2 sumd1g rate rtable;

%  Change surface height,
%  initialize height variables:

if iflag >= 0
    idh = 1;
else
    idh = -1;
end

ih(iy,ix) = ih(iy,ix) + idh;
ihc = ih(iy,ix);
ihcp = ihc + 1;

if iflag >= 0
    ihr1 = ihc - 1;
    ihr2 = ihc - 2;
    ihr3 = ihc;
else
    ihr1 = ihc;
    ihr2 = ihc - 1;
    ihr3 = ihc + 1;
end

%  Keep track of number of atoms at even heights:

if mod(ihc,2) == 0
    neven = neven + 1;
else
    neven = neven - 1;
end

%  Initialize pointers to neighbors,
%  taking into account periodic boundaries:

iyp = iy + 1;
ixp = ix + 1;
iypp = iy + 2;
ixpp = ix + 2;
iym = iy - 1;
ixm = ix - 1;
iymm = iy - 2;
ixmm = ix - 2;
if iyp > nrow; iyp = iyp - nrow; end
if ixp > ncol; ixp = ixp - ncol; end
if iypp > nrow; iypp = iypp - nrow; end
if ixpp > ncol; ixpp = ixpp - ncol; end
if iym <= 0; iym = iym + nrow; end
if ixm <= 0; ixm = ixm + ncol; end
if iymm <= 0; iymm = iymm + nrow; end
if ixmm <= 0; ixmm = ixmm + ncol; end

%  Calculate groups of pointers:

iyg = floor((iy-1)/kgrp) + 1;
ixg = floor((ix-1)/kgrp) + 1;
iypg = floor((iyp-1)/kgrp) + 1;
ixpg = floor((ixp-1)/kgrp) + 1;
iyppg = floor((iypp-1)/kgrp) + 1;
ixppg = floor((ixpp-1)/kgrp) + 1;
iymg = floor((iym-1)/kgrp) + 1;
ixmg = floor((ixm-1)/kgrp) + 1;
iymmg = floor((iymm-1)/kgrp) + 1;
ixmmg = floor((ixmm-1)/kgrp) + 1;

%  Initialize flags for changes to second-neighbor rates:

nnnfl = zeros(4,1);

%  Update bonding arrays ibi and ibf:

%  Count new number of bonds of center,
%  number of bonds of center with additional atom:

nbr = 1;
nbrp = 1;

ihn = ih(iy,ixp);
if ixp == 1; ihn = ihn - nsteps; end
if ihc <= ihn
    nbr = nbr + 1;
    if ihc < ihn; nbrp = nbrp + 1; end
end

ihn = ih(iyp,ix);
if ihc <= ihn
    nbr = nbr + 1;
    if ihc < ihn; nbrp = nbrp + 1; end
end

ihn = ih(iy,ixm);
if ixm == ncol; ihn = ihn + nsteps; end
if ihc <= ihn
    nbr = nbr + 1;
    if ihc < ihn; nbrp = nbrp + 1; end
end

ihn = ih(iym,ix);
if ihc <= ihn
    nbr = nbr + 1;
    if ihc < ihn; nbrp = nbrp + 1; end
end

%  For moves from center,
%  change initial bond number:

ibi(iy,ix) = nbr;

%  For moves to center,
%  change final bond numbers:

ibf(3,iy,ixp) = nbrp;
ibf(4,iyp,ix) = nbrp;
ibf(1,iy,ixm) = nbrp;
ibf(2,iym,ix) = nbrp;

%  For deposition to center,
%  change final bond number:

ibfd(iy,ix) = nbrp;

%  Make adjustments based on relative height of neighbor:

%  For moves to neighbor from second neighbors
%  and deposition into neighbor,
%  adjust final bond numbers:

ihn = ih(iy,ixp);
if ixp == 1; ihn = ihn - nsteps; end
if ihr1 == ihn
    nnnfl(1) = 1;
    ibfd(iy,ixp) = ibfd(iy,ixp) + idh;
    ibf(2,iym,ixp) = ibf(2,iym,ixp) + idh;
    ibf(3,iy,ixpp) = ibf(3,iy,ixpp) + idh;
    ibf(4,iyp,ixp) = ibf(4,iyp,ixp) + idh;

%   For move from center to neighbor,
%   adjust final bond number:

elseif ihr2 == ihn
    ibf(1,iy,ix) = ibf(1,iy,ix) + idh;

%   For move to center from neighbor,
%   adjust initial bond number:

elseif ihr3 == ihn
    ibi(iy,ixp) = ibi(iy,ixp) + idh;
end

%   For move to center from neighbor,
%   adjust final bond number:

if ihcp == ihn; ibf(3,iy,ixp) = ibf(3,iy,ixp) - 1; end

%   Repeat adjustments based on heights of other neighbors:

ihn = ih(iyp,ix);
if ihr1 == ihn
    nnnfl(2) = 1;
    ibfd(iyp,ix) = ibfd(iyp,ix) + idh;
    ibf(3,iyp,ixp) = ibf(3,iyp,ixp) + idh;
    ibf(4,iypp,ix) = ibf(4,iypp,ix) + idh;
    ibf(1,iyp,ixm) = ibf(1,iyp,ixm) + idh;
elseif ihr2 == ihn
    ibf(2,iy,ix) = ibf(2,iy,ix) + idh;
elseif ihr3 == ihn
    ibi(iyp,ix) = ibi(iyp,ix) + idh;
end
if ihcp == ihn; ibf(4,iyp,ix) = ibf(4,iyp,ix) - 1; end

ihn = ih(iy,ixm);
if ixm == ncol; ihn = ihn + nsteps; end
if ihr1 == ihn
    nnnfl(3) = 1;
    ibfd(iy,ixm) = ibfd(iy,ixm) + idh;
    ibf(4,iyp,ixm) = ibf(4,iyp,ixm) + idh;
    ibf(1,iy,ixmm) = ibf(1,iy,ixmm) + idh;
    ibf(2,iym,ixm) = ibf(2,iym,ixm) + idh;
elseif ihr2 == ihn
    ibf(3,iy,ix) = ibf(3,iy,ix) + idh;
elseif ihr3 == ihn
    ibi(iy,ixm) = ibi(iy,ixm) + idh;
end
if ihcp == ihn; ibf(1,iy,ixm) = ibf(1,iy,ixm) - 1; end

ihn = ih(iym,ix);
if ihr1 == ihn
    nnnfl(4) = 1;
    ibfd(iym,ix) = ibfd(iym,ix) + idh;
    ibf(1,iym,ixm) = ibf(1,iym,ixm) + idh;
    ibf(2,iymm,ix) = ibf(2,iymm,ix) + idh;
    ibf(3,iym,ixp) = ibf(3,iym,ixp) + idh;
elseif ihr2 == ihn
    ibf(4,iy,ix) = ibf(4,iy,ix) + idh;
elseif ihr3 == ihn
    ibi(iym,ix) = ibi(iym,ix) + idh;
end
if ihcp == ihn; ibf(2,iym,ix) = ibf(2,iym,ix) - 1; end

%   Update rate array and sum vectors:

sum3         = sum3         - sum2(ix);
sum2g(ixg)   = sum2g(ixg)   - sum2(ix);
sumd3        = sumd3        - sumd2(ix);
sumd2g(ixg)  = sumd2g(ixg)  - sumd2(ix);
sum3         = sum3         - sum2(ixp);
sum2g(ixpg)  = sum2g(ixpg)  - sum2(ixp);
sumd3        = sumd3        - sumd2(ixp);
sumd2g(ixpg) = sumd2g(ixpg) - sumd2(ixp);
sum3         = sum3         - sum2(ixm);
sum2g(ixmg)  = sum2g(ixmg)  - sum2(ixm);
sumd3        = sumd3        - sumd2(ixm);
sumd2g(ixmg) = sumd2g(ixmg) - sumd2(ixm);

%   Moves from center and neighbors:

sum2(ix)        = sum2(ix)        - sum1(iy,ix);
sum1g(iyg,ix)   = sum1g(iyg,ix)   - sum1(iy,ix);
sumd2(ix)       = sumd2(ix)       - rate(6,iy,ix);
sumd1g(iyg,ix)  = sumd1g(iyg,ix)  - rate(6,iy,ix);
sum2(ix)        = sum2(ix)        - sum1(iyp,ix);
sum1g(iypg,ix)  = sum1g(iypg,ix)  - sum1(iyp,ix);
sumd2(ix)       = sumd2(ix)       - rate(6,iyp,ix);
sumd1g(iypg,ix) = sumd1g(iypg,ix) - rate(6,iyp,ix);
sum2(ix)        = sum2(ix)        - sum1(iym,ix);
sum1g(iymg,ix)  = sum1g(iymg,ix)  - sum1(iym,ix);
sumd2(ix)       = sumd2(ix)       - rate(6,iym,ix);
sumd1g(iymg,ix) = sumd1g(iymg,ix) - rate(6,iym,ix);
sum2(ixp)       = sum2(ixp)       - sum1(iy,ixp);
sum1g(iyg,ixp)  = sum1g(iyg,ixp)  - sum1(iy,ixp);
sumd2(ixp)      = sumd2(ixp)      - rate(6,iy,ixp);
sumd1g(iyg,ixp) = sumd1g(iyg,ixp) - rate(6,iy,ixp);
sum2(ixm)       = sum2(ixm)       - sum1(iy,ixm);
sum1g(iyg,ixm)  = sum1g(iyg,ixm)  - sum1(iy,ixm);
sumd2(ixm)      = sumd2(ixm)      - rate(6,iy,ixm);
sumd1g(iyg,ixm) = sumd1g(iyg,ixm) - rate(6,iy,ixm);

sum1(iy,ix)  = 0.;
sum1(iy,ixp) = 0.;
sum1(iyp,ix) = 0.;
sum1(iy,ixm) = 0.;
sum1(iym,ix) = 0.;

for idir = 1:4

rate(idir,iy,ix)  = rtable(ibi(iy,ix), ibf(idir,iy,ix));
rate(idir,iy,ixp) = rtable(ibi(iy,ixp),ibf(idir,iy,ixp));
rate(idir,iyp,ix) = rtable(ibi(iyp,ix),ibf(idir,iyp,ix));
rate(idir,iy,ixm) = rtable(ibi(iy,ixm),ibf(idir,iy,ixm));
rate(idir,iym,ix) = rtable(ibi(iym,ix),ibf(idir,iym,ix));

sum1(iy,ix)  = sum1(iy,ix)  + rate(idir,iy,ix);
sum1(iy,ixp) = sum1(iy,ixp) + rate(idir,iy,ixp);
sum1(iyp,ix) = sum1(iyp,ix) + rate(idir,iyp,ix);
sum1(iy,ixm) = sum1(iy,ixm) + rate(idir,iy,ixm);
sum1(iym,ix) = sum1(iym,ix) + rate(idir,iym,ix);

end

rate(5,iy,ix) = rtable(ibi(iy,ix),6);
rate(6,iy,ix) = rtable(6,ibfd(iy,ix));
sum1(iy,ix) = sum1(iy,ix) + rate(5,iy,ix);

rate(5,iy,ixp) = rtable(ibi(iy,ixp),6);
rate(6,iy,ixp) = rtable(6,ibfd(iy,ixp));
sum1(iy,ixp) = sum1(iy,ixp) + rate(5,iy,ixp);

rate(5,iyp,ix) = rtable(ibi(iyp,ix),6);
rate(6,iyp,ix) = rtable(6,ibfd(iyp,ix));
sum1(iyp,ix) = sum1(iyp,ix) + rate(5,iyp,ix);

rate(5,iy,ixm) = rtable(ibi(iy,ixm),6);
rate(6,iy,ixm) = rtable(6,ibfd(iy,ixm));
sum1(iy,ixm) = sum1(iy,ixm) + rate(5,iy,ixm);

rate(5,iym,ix) = rtable(ibi(iym,ix),6);
rate(6,iym,ix) = rtable(6,ibfd(iym,ix));
sum1(iym,ix) = sum1(iym,ix) + rate(5,iym,ix);

sum2(ix)        = sum2(ix)        + sum1(iy,ix);
sum1g(iyg,ix)   = sum1g(iyg,ix)   + sum1(iy,ix);
sumd2(ix)       = sumd2(ix)       + rate(6,iy,ix);
sumd1g(iyg,ix)  = sumd1g(iyg,ix)  + rate(6,iy,ix);
sum2(ix)        = sum2(ix)        + sum1(iyp,ix);
sum1g(iypg,ix)  = sum1g(iypg,ix)  + sum1(iyp,ix);
sumd2(ix)       = sumd2(ix)       + rate(6,iyp,ix);
sumd1g(iypg,ix) = sumd1g(iypg,ix) + rate(6,iyp,ix);
sum2(ix)        = sum2(ix)        + sum1(iym,ix);
sum1g(iymg,ix)  = sum1g(iymg,ix)  + sum1(iym,ix);
sumd2(ix)       = sumd2(ix)       + rate(6,iym,ix);
sumd1g(iymg,ix) = sumd1g(iymg,ix) + rate(6,iym,ix);
sum2(ixp)       = sum2(ixp)       + sum1(iy,ixp);
sum1g(iyg,ixp)  = sum1g(iyg,ixp)  + sum1(iy,ixp);
sumd2(ixp)      = sumd2(ixp)      + rate(6,iy,ixp);
sumd1g(iyg,ixp) = sumd1g(iyg,ixp) + rate(6,iy,ixp);
sum2(ixm)       = sum2(ixm)       + sum1(iy,ixm);
sum1g(iyg,ixm)  = sum1g(iyg,ixm)  + sum1(iy,ixm);
sumd2(ixm)      = sumd2(ixm)      + rate(6,iy,ixm);
sumd1g(iyg,ixm) = sumd1g(iyg,ixm) + rate(6,iy,ixm);

%  Moves from second-neighbors to neighbors:

if nnnfl(1) ~= 0

    sum1(iym,ixp)   = sum1(iym,ixp)   - rate(2,iym,ixp);
    sum1g(iymg,ixp) = sum1g(iymg,ixp) - rate(2,iym,ixp);
    sum2(ixp)       = sum2(ixp)       - rate(2,iym,ixp);
    rate(2,iym,ixp) = rtable(ibi(iym,ixp),ibf(2,iym,ixp));
    sum1(iym,ixp)   = sum1(iym,ixp)   + rate(2,iym,ixp);
    sum1g(iymg,ixp) = sum1g(iymg,ixp) + rate(2,iym,ixp);
    sum2(ixp)       = sum2(ixp)       + rate(2,iym,ixp);

    sum1(iy,ixpp)   = sum1(iy,ixpp)   - rate(3,iy,ixpp);
    sum1g(iyg,ixpp) = sum1g(iyg,ixpp) - rate(3,iy,ixpp);
    sum2(ixpp)      = sum2(ixpp)      - rate(3,iy,ixpp);
    sum2g(ixppg)    = sum2g(ixppg)    - rate(3,iy,ixpp);
    sum3            = sum3            - rate(3,iy,ixpp);
    rate(3,iy,ixpp) = rtable(ibi(iy,ixpp),ibf(3,iy,ixpp));
    sum1(iy,ixpp)   = sum1(iy,ixpp)   + rate(3,iy,ixpp);
    sum1g(iyg,ixpp) = sum1g(iyg,ixpp) + rate(3,iy,ixpp);
    sum2(ixpp)      = sum2(ixpp)      + rate(3,iy,ixpp);
    sum2g(ixppg)    = sum2g(ixppg)    + rate(3,iy,ixpp);
    sum3            = sum3            + rate(3,iy,ixpp);

    sum1(iyp,ixp)   = sum1(iyp,ixp)   - rate(4,iyp,ixp);
    sum1g(iypg,ixp) = sum1g(iypg,ixp) - rate(4,iyp,ixp);
    sum2(ixp)       = sum2(ixp)       - rate(4,iyp,ixp);
    rate(4,iyp,ixp) = rtable(ibi(iyp,ixp),ibf(4,iyp,ixp));
    sum1(iyp,ixp)   = sum1(iyp,ixp)   + rate(4,iyp,ixp);
    sum1g(iypg,ixp) = sum1g(iypg,ixp) + rate(4,iyp,ixp);
    sum2(ixp)       = sum2(ixp)       + rate(4,iyp,ixp);

end

if nnnfl(2) ~= 0

    sum1(iyp,ixp)   = sum1(iyp,ixp)   - rate(3,iyp,ixp);
    sum1g(iypg,ixp) = sum1g(iypg,ixp) - rate(3,iyp,ixp);
    sum2(ixp)       = sum2(ixp)       - rate(3,iyp,ixp);
    rate(3,iyp,ixp) = rtable(ibi(iyp,ixp),ibf(3,iyp,ixp));
    sum1(iyp,ixp)   = sum1(iyp,ixp)   + rate(3,iyp,ixp);
    sum1g(iypg,ixp) = sum1g(iypg,ixp) + rate(3,iyp,ixp);
    sum2(ixp)       = sum2(ixp)       + rate(3,iyp,ixp);

    sum1(iypp,ix)   = sum1(iypp,ix)   - rate(4,iypp,ix);
    sum1g(iyppg,ix) = sum1g(iyppg,ix) - rate(4,iypp,ix);
    sum2(ix)        = sum2(ix)        - rate(4,iypp,ix);
    rate(4,iypp,ix) = rtable(ibi(iypp,ix),ibf(4,iypp,ix));
    sum1(iypp,ix)   = sum1(iypp,ix)   + rate(4,iypp,ix);
    sum1g(iyppg,ix) = sum1g(iyppg,ix) + rate(4,iypp,ix);
    sum2(ix)        = sum2(ix)        + rate(4,iypp,ix);

    sum1(iyp,ixm)   = sum1(iyp,ixm)   - rate(1,iyp,ixm);
    sum1g(iypg,ixm) = sum1g(iypg,ixm) - rate(1,iyp,ixm);
    sum2(ixm)       = sum2(ixm)       - rate(1,iyp,ixm);
    rate(1,iyp,ixm) = rtable(ibi(iyp,ixm),ibf(1,iyp,ixm));
    sum1(iyp,ixm)   = sum1(iyp,ixm)   + rate(1,iyp,ixm);
    sum1g(iypg,ixm) = sum1g(iypg,ixm) + rate(1,iyp,ixm);
    sum2(ixm)       = sum2(ixm)       + rate(1,iyp,ixm);

end

if nnnfl(3) ~= 0

    sum1(iyp,ixm)   = sum1(iyp,ixm)   - rate(4,iyp,ixm);
    sum1g(iypg,ixm) = sum1g(iypg,ixm) - rate(4,iyp,ixm);
    sum2(ixm)       = sum2(ixm)       - rate(4,iyp,ixm);
    rate(4,iyp,ixm) = rtable(ibi(iyp,ixm),ibf(4,iyp,ixm));
    sum1(iyp,ixm)   = sum1(iyp,ixm)   + rate(4,iyp,ixm);
    sum1g(iypg,ixm) = sum1g(iypg,ixm) + rate(4,iyp,ixm);
    sum2(ixm)       = sum2(ixm)       + rate(4,iyp,ixm);

    sum1(iy,ixmm)   = sum1(iy,ixmm)   - rate(1,iy,ixmm);
    sum1g(iyg,ixmm) = sum1g(iyg,ixmm) - rate(1,iy,ixmm);
    sum2(ixmm)      = sum2(ixmm)      - rate(1,iy,ixmm);
    sum2g(ixmmg)    = sum2g(ixmmg)    - rate(1,iy,ixmm);
    sum3            = sum3            - rate(1,iy,ixmm);
    rate(1,iy,ixmm) = rtable(ibi(iy,ixmm),ibf(1,iy,ixmm));
    sum1(iy,ixmm)   = sum1(iy,ixmm)   + rate(1,iy,ixmm);
    sum1g(iyg,ixmm) = sum1g(iyg,ixmm) + rate(1,iy,ixmm);
    sum2(ixmm)      = sum2(ixmm)      + rate(1,iy,ixmm);
    sum2g(ixmmg)    = sum2g(ixmmg)    + rate(1,iy,ixmm);
    sum3            = sum3            + rate(1,iy,ixmm);

    sum1(iym,ixm)   = sum1(iym,ixm)   - rate(2,iym,ixm);
    sum1g(iymg,ixm) = sum1g(iymg,ixm) - rate(2,iym,ixm);
    sum2(ixm)       = sum2(ixm)       - rate(2,iym,ixm);
    rate(2,iym,ixm) = rtable(ibi(iym,ixm),ibf(2,iym,ixm));
    sum1(iym,ixm)   = sum1(iym,ixm)   + rate(2,iym,ixm);
    sum1g(iymg,ixm) = sum1g(iymg,ixm) + rate(2,iym,ixm);
    sum2(ixm)       = sum2(ixm)       + rate(2,iym,ixm);

end

if nnnfl(4) ~= 0

    sum1(iym,ixm)   = sum1(iym,ixm)   - rate(1,iym,ixm);
    sum1g(iymg,ixm) = sum1g(iymg,ixm) - rate(1,iym,ixm);
    sum2(ixm)       = sum2(ixm)       - rate(1,iym,ixm);
    rate(1,iym,ixm) = rtable(ibi(iym,ixm),ibf(1,iym,ixm));
    sum1(iym,ixm)   = sum1(iym,ixm)   + rate(1,iym,ixm);
    sum1g(iymg,ixm) = sum1g(iymg,ixm) + rate(1,iym,ixm);
    sum2(ixm)       = sum2(ixm)       + rate(1,iym,ixm);

    sum1(iymm,ix)   = sum1(iymm,ix)   - rate(2,iymm,ix);
    sum1g(iymmg,ix) = sum1g(iymmg,ix) - rate(2,iymm,ix);
    sum2(ix)        = sum2(ix)        - rate(2,iymm,ix);
    rate(2,iymm,ix) = rtable(ibi(iymm,ix),ibf(2,iymm,ix));
    sum1(iymm,ix)   = sum1(iymm,ix)   + rate(2,iymm,ix);
    sum1g(iymmg,ix) = sum1g(iymmg,ix) + rate(2,iymm,ix);
    sum2(ix)        = sum2(ix)        + rate(2,iymm,ix);

    sum1(iym,ixp)   = sum1(iym,ixp)   - rate(3,iym,ixp);
    sum1g(iymg,ixp) = sum1g(iymg,ixp) - rate(3,iym,ixp);
    sum2(ixp)       = sum2(ixp)       - rate(3,iym,ixp);
    rate(3,iym,ixp) = rtable(ibi(iym,ixp),ibf(3,iym,ixp));
    sum1(iym,ixp)   = sum1(iym,ixp)   + rate(3,iym,ixp);
    sum1g(iymg,ixp) = sum1g(iymg,ixp) + rate(3,iym,ixp);
    sum2(ixp)       = sum2(ixp)       + rate(3,iym,ixp);

end

sum3         = sum3         + sum2(ix);
sum2g(ixg)   = sum2g(ixg)   + sum2(ix);
sumd3        = sumd3        + sumd2(ix);
sumd2g(ixg)  = sumd2g(ixg)  + sumd2(ix);
sum3         = sum3         + sum2(ixp);
sum2g(ixpg)  = sum2g(ixpg)  + sum2(ixp);
sumd3        = sumd3        + sumd2(ixp);
sumd2g(ixpg) = sumd2g(ixpg) + sumd2(ixp);
sum3         = sum3         + sum2(ixm);
sum2g(ixmg)  = sum2g(ixmg)  + sum2(ixm);
sumd3        = sumd3        + sumd2(ixm);
sumd2g(ixmg) = sumd2g(ixmg) + sumd2(ixm);

return

end

