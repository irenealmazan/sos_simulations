%This script calculates the 3D data set obtained by slicing the reciprocal
%space along the L direction 

close all;clear all;

skip = 1;
path = 'miscut_8_14/';
runname = [path 'miscut_8_14'];
runname_title = 'miscut\_8\_14';
namefile = [runname '_sf_allL.mat'];

addpath(genpath(['/Users/ialmazn/Box' ' Sync/XPCS_sputtering_ZnO_TiO2_2017_2019/XPCS_analysis_on_thefly']));

if ~skip
    
    load([runname '_stats.mat']);

    if ~exist('nsteps','var'); nsteps = 0; end

    if ~exist('damono','var')
        havg = squeeze(mean(mean(ihm))); % Average height, gives growth amount
        damono = havg - 1 - nsteps/2;
    end

    
    % Calculate for an array of Ls
    L = [0.1:.01:1.5];
    %L = 0.6;
    %L = 0.7;
    %L = 0.8;
    %L = 0.9;
    %L = 1.0;
    %L = 1.5;

% This works for zero miscut:
% Anti-Bragg scattering will have opposite phase for odd and even heights
%    phm = 2*mod(ihm,2) - 1;

% For miscut, each column has different phase

% Calculate for arbitrary L
    zsub = nsteps*ones(nrow,1)*(1 - [1:ncol]/ncol); % height ramp from miscut
    ihsub = floor(zsub); % height of top atoms in substrate
    
    ftm = NaN*ones(nrow,ncol,numel(L));

    
    for lll = 1:numel(L)
        xLsub = exp(2i*pi*L(lll)*(ihsub - zsub));
    xL = exp(2i*pi*L(lll));

% Calculate CTRs of substrate
    mm = ncol/nsteps; % number of unit cells per terrace of substrate (even)

    % Use formulas from Trainor, surface coord
    % See miscut_m_CTR3_LK.m
    
    offsetz = -1 + 1/(mm); % substrate offset in z, unit cells; works for nsteps = 8
    offsetx = -1; % substrate offset in x, unit cells
    ftms = zeros(nrow,ncol);
    
    % Center of ftm(:,:,ii) is nrow/2+1, ncol/2+1 (e.g. 65 65 for 128)

    iy = nrow/2 + 1;
    ix = ncol/2 + 1;
    %hsctr = [0:mm-1];
    for cc = -mm/2:mm/2-1 % Loop over all CTRs
        % CTRs are spaced by ncol/mm = nsteps pixels
        ixc = ix + cc*nsteps;
        hhh = cc/mm + L(lll)/mm; % true H value along CTR
        x = exp(2i*pi*hhh);
        num = exp(2i*pi*offsetz*L(lll))*exp(2i*pi*offsetx*hhh)*x/mm; % to match shift of substrate top from sim by one unit cell
        denom = x - 1;
        if abs(denom) > 0.001 % H, L value not on Bragg position
            ftms(iy,ixc) = num/denom; % comparable to ifft2
        else % H, L value on Bragg position
            ftms(iy,ixc) = 1e4; % comparable to ifft2
        end
        
        
        
    end
    
  
    % calculate slice of reciprocal space:
    
    % Loop over all time steps
    nt = 1;%size(ihm,3); % Number of time steps
    
   % Do not use first dt_minML of growth for time correlations
    dt_minML = 2;%if in ML
    idt = damono > dt_minML;
    
    
    ii =  find(idt>0,1);
    
    
    ph = xLsub.*(1 - xL.^(ihm(:,:,ii) - ihsub))./(1 - xL);
    ftm(:,:,lll) = fftshift(ifft2(ph)) +ftms;    
    
    
    
    
    end
    
    
     III = conj(ftm).*ftm;
    
    save([namefile]);

else
    load([namefile]);
end

% central CTR is at nrow/2+1, ncol/2+1
ixcen = ncol/2 + 1;
iycen = nrow/2 + 1;
% DQ of CTRs is nsteps (pixels) in x (col)

QX = [1:ncol]-ixcen;
QY = [1:nrow]-iycen;


L_index = 7; 

figure(1000);
clf;

subplot(131);
imagesc(QX,QY,log10(III(:,:,L_index)));
title(['L = ' num2str(L(L_index))]);
xlabel('Q_x (pixels)');
ylabel('Q_y (pixels)')

subplot(132);
imagesc(QX,L,log10(squeeze(III(:,64,:))));
set(gca,'YDir','normal');
title(['L = ' num2str(L(L_index))]);
xlabel('Q_x');
ylabel('L');
title('Q_y = 64');

subplot(133);
imagesc(QY,L,log10(squeeze(III(64,:,:))));
title(['L = ' num2str(L(L_index))]);
set(gca,'YDir','normal');
xlabel('Q_y');
ylabel('L');
title('Q_x = 64');

