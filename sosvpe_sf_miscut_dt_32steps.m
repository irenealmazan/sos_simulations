% Script to plot structure factor from sosvpe simulations
% From sosvpe_sf_2t_new2
% Dropped azimuthal averages
% Added Cdt calc
% 22-FEB-19 GBS



path = 'miscut_32_29/';%'miscut_32_40/';%'miscut_32_28/';
runname = [path 'miscut_32_29'];%[path 'miscut_32_40'];%[path 'miscut_32_28'];
runname_title = 'miscut\_32\_29';%'miscut\_32\_40';%'miscut\_32\_28';
namefile = [runname '_corr_dt_L0p5.mat'];%[runname 'half_corr_dt_L0p5.mat'];
damHWfilename = [path 'damHW.mat'];

isequilibrium = 1;
skip = 1;


color_rates = 'b';
growth_rates = [0];%[1.65e-9/32];%[1.65e-9]/2;%[0.0];%2*[1.65e-9];


addpath(genpath(['/Users/ialmazn/Documents/data_analysis_ALL/XPCS_sputtering_ZnO_TiO2_2017_2019/XPCS_analysis_on_thefly']));


%make_vid(runname,[0 10],'jet');
%make_vid_sf(runname,[0 100],'jet');

if ~skip
    
   
    
    %
     load([runname '_stats.mat']); 

     if strcmp(runname_title ,[ 'miscut\_32\_26'])
         load([runname '_ihm_last2000.mat']);
         ihm = ihm_resize;
         clear ihm_resize;
     elseif strcmp(runname_title ,[ 'miscut\_32\_32'])
         load([runname '_ihm.mat']);
     elseif strcmp(runname_title ,[ 'miscut\_32\_33'])
         load([runname '_ihm.mat']);
    elseif strcmp(runname_title ,[ 'miscut\_32\_37'])
        load([runname '_ihm_half.mat']);
         ihm = ihm_half;
         dtime_orig = load([runname '_stats.mat'],'dtime'); 
         dtime = dtime_orig.dtime(end-round(size(ihm_half,3)/2):end);
         clear ihm_half
         clear dtime_orig
    elseif strcmp(runname_title ,[ 'miscut\_32\_38'])
        load([runname '_ihm_half.mat']);
         ihm = ihm_half;
         dtime_orig = load([runname '_stats.mat'],'dtime'); 
         dtime = dtime_orig.dtime(end-size(ihm_half,3):end);
         clear ihm_half
         clear dtime_orig
         
    elseif strcmp(runname_title ,[ 'miscut\_32\_39'])
        load([runname '_ihm_half.mat']);
         ihm = ihm_half;
         dtime_orig = load([runname '_stats.mat'],'dtime'); 
         dtime = dtime_orig.dtime(end-size(ihm_half,3):end);
         clear ihm_half
         clear dtime_orig
        
    elseif strcmp(runname_title ,[ 'miscut\_32\_40'])
        load([runname '_ihm_half.mat']);
         ihm = ihm_half;
         dtime_orig = load([runname '_stats.mat'],'dtime'); 
         dtime = dtime_orig.dtime(end-size(ihm_half,3):end);
         clear ihm_half
         clear dtime_orig
     else
        load([runname '_ihm.mat']);
    end
    
    
    if ~exist('nsteps','var'); nsteps = 0; end

    if ~exist('damono','var')
        havg = squeeze(mean(mean(ihm))); % Average height, gives growth amount
        damono = havg - 1 - nsteps/2;
    end

    dt_minML = dtime(end-2.50e4);%dtime(30);%dtime(end-round(2.5e4/8));%dtime(end-6.3e3);% % dtime(end-1e4);
    dt_maxML = dtime(end);%dtime(250);%dtime(end-1.25e4);%
    %idt = dtime >dt_minML;
    idt = (dtime>dt_minML & dtime<dt_maxML);
    ddam = dtime(idt);%damono(idt);
    %{
    if ~isequilibrium
        dt_minML = 1.99;%if in ML
        idt = damono > dt_minML;
        ddam = dtime(idt);%damono(idt);
    else
        %dt_minML = dtime(end-2000);%;dtime(end-1000);%1.4160e+09;%2.02e10;%6.0606e+08;% if in seconds which is 2/3.3e-9 to preserve the same parameters than before%2;if in ML
        %dt_maxML = dtime(end-10e3);
        dt_maxML = dtime(end-1e3);
        %idt = dtime >dt_minML;
        idt = dtime>dt_maxML;
        ddam = dtime(idt);%damono(idt);
    end
    %}
    
    ihm_orig = ihm(:,:,idt);
    clear ihm;
    
    % store the original number of pixels and redefine the window over which we want to calculate the fft:
    ncol_orig = ncol;
    ncol = ncol_orig/4;
    ncol_range = [-ncol/2+1:ncol/2]+ncol_orig/2;
    
    % redefine the number of steps:
    nsteps_orig = nsteps;
    nsteps  = nsteps_orig/4;
    
    ihm = ihm_orig(:,ncol_range,:);
    
    
   
    % Calculate for fixed L
    L = 0.5;

% This works for zero miscut:
% Anti-Bragg scattering will have opposite phase for odd and even heights
%    phm = 2*mod(ihm,2) - 1;

% For miscut, each column has different phase

% select a narrow window of the surface to calculate the fft:
    zsub_orig = nsteps_orig*ones(nrow,1)*(1 - [1:ncol_orig]/ncol_orig); % height ramp from miscut
    zsub = zsub_orig(:,ncol_range);%nsteps*ones(nrow,1)*(1 - [1:ncol]/ncol); % height ramp from miscut
    ihsub = floor(zsub); % height of top atoms in substrate
    
    
    
    % Calculate for the fixed L defined above
    xLsub = exp(2i*pi*L*(ihsub - zsub));
    xL = exp(2i*pi*L);

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
        hhh = cc/mm + L/mm; % true H value along CTR
        x = exp(2i*pi*hhh);
        num = exp(2i*pi*offsetz*L)*exp(2i*pi*offsetx*hhh)*x/mm; % to match shift of substrate top from sim by one unit cell
        denom = x - 1;
        if abs(denom) > 0.001 % H, L value not on Bragg position
            ftms(iy,ixc) = num/denom; % comparable to ifft2
        else % H, L value on Bragg position
            ftms(iy,ixc) = 1e4; % comparable to ifft2
        end
        
        
        
    end
    
  
    
    % Loop over all time steps
    nt = size(ihm,3); % Number of time steps
    ftm = NaN*ones(nrow,ncol,nt);
    
    
    %figure(20);
    %clf;
    for ii = 1:nt
        %step_probe = probe.*(ihm(:,:,ii) - ihsub);
        %ph = xLsub.*(1 - xL.^(step_probe))./(1 - xL);
        ph = xLsub.*(1 - xL.^(ihm(:,:,ii) - ihsub))./(1 - xL);
        %ph = probe.*(xLsub.*(1 - xL.^(ihm(:,:,ii) - ihsub))./(1 - xL));
        
        % Calculate Fourier transform of phased height
        % It matters here whether to do fft2 or ifft2; 
        % ifft2 does exp(+iqr), and divides by nrow * ncol
        % also add substrate contribution to CTR positions here 
        ftm(:,:,ii) = fftshift(ifft2(ph)) + ftms;
      
   
    end
    
    
   
    

    % Calc delta-t correlations for each pixel
    
    % Do not use first dt_minML of growth for time correlations
    
    
    
    III = abs(ftm).^2;
    ndt = size(III,3);
    ddam = ddam - ddam(1); % delta time coord (ML) (assumes evenly spaced)
    
    % Use 2-D Savitzky-Golay for getting Ibar with imputs:
    use_2Dsg_smooth = 1;
    if use_2Dsg_smooth
        % remove the CTRs:
        CTR_pixels_col = 65+nsteps*[-2 -1 0 1];
        CTR_pixels_row = 65;
        
        delta_pixel_row = [-1 0 0 1];
        delta_pixel_col = [0 -1 1 0];   
        
        for cc = 1:numel(CTR_pixels_col)
            for tt = 1:size(III,3)
                III_mean  = 0;
                for ll = 1:numel(delta_pixel_row)
                   III_mean = III_mean + III(CTR_pixels_row+delta_pixel_row(ll),CTR_pixels_col(cc)+delta_pixel_col(ll),tt);
                end
                III_mean = III_mean/numel(delta_pixel_row);
                III(CTR_pixels_row,CTR_pixels_col(cc),tt) = III_mean;
            end
        end
        
        
        maxpd = 2; % determines maximum degree of smoothing polynomial
        iii = 1;%2;%3; %5;% half-width of pixel range in columns
        jjj = 1;%3;%2;%3; %5;% half-width of pixel range in rows
        
        [Ibar,filt] = Functions_sos.smooth_sg(III,maxpd,iii,jjj);%mean(III,3);
        
        Ibar_norm = squeeze(sum(sum(Ibar,1),2));
        
        % redefinition of nrow and ncol
        nrow_crop_orig = nrow;
        ncol_crop_orig = ncol;
        
        nrow = size(Ibar,1);
        ncol = size(Ibar,2);
        
        
        III_resize = III(jjj+1:jjj+nrow,iii+1:iii+ncol,:);
        
        III_norm = squeeze(sum(sum(III_resize,1),2));
       
        for tt = 1:size(III_resize,3)          
            dlnI(:,:,tt) = (Ibar_norm(tt).*III_resize(:,:,tt))./(III_norm(tt).*Ibar(:,:,tt)) - 1;           
        end
        figNum = 100;
        row_to_plot = 50;
        frames_to_plot = [1:500:1000];%[1:500:2000];
        %Function_display_sos.show_filt_Ibar(III,Ibar,dlnI,filt,maxpd,row_to_plot,frames_to_plot,figNum);
        
    else
        Ibar = mean(III,3);
        III_resize = III;
        dlnI = III./Ibar - 1;
    end
    
    Cdt = NaN*ones(size(dlnI));
    I1 = ones(ndt,1);
    Idt0 = conv(I1,I1);
    
    for ii = 1:nrow%nrow
        for jj = 1:ncol%ncol
            It = squeeze(dlnI(ii,jj,:));
            Idt = conv(It,flip(It))./Idt0;
            Cdt(ii,jj,:) = Idt(end-ndt+1:end);
        end
    end
    
    % calculate the contrast:
    rowcen = size(Cdt,1)/2;colcen = size(Cdt,2)/2;
    Cdt_mean = mean(mean(Cdt(rowcen-7:rowcen-2,colcen-5:colcen+5,1)));
    
    save([namefile],'III_resize','Ibar','Cdt','ddam','nsteps','ncol_orig','-v7.3');

else
    load([namefile]);
end


%% Initialize important parameters:

nrow = size(III_resize,1);
ncol = size(III_resize,2);
ndt = size(III_resize,3);

if ~exist('nsteps','var')
    nsteps = 8;
end

if ~exist('ncol_orig','var')
    ncol_orig = ncol*4;
end

%% Get time HW of Cdt for each pixel, units of ML



damHW = ddam(end)*ones(nrow,ncol);

skip_HWcalc = 1;
fit_Cdti = 0;
if  ~skip_HWcalc
    
    for ii = 1:nrow
        for jj = 1:ncol
            Cdti = squeeze(Cdt(ii,jj,:));
            Cdti = Cdti/Cdti(1);
            err1 = 1;
            
           fitfunc_str = 'FittingFunctions.CCN2single_fit';

            if fit_Cdti                
                %%{
                fit_range = [1:1:round(1/2.5*length(Cdti))];
                pin = [0 1 1e4 0];%pin_iiT(iT,:);
                dp =  [[0 1 1 0]*0.0001];%dp_iiT(iT,:);
                w = ones(length(fit_range),1);
                
                CCfunc.time_1D = ddam';
                CCfunc.CCNdtV = Cdti;
                
                [ fitres] = FittingFunctions.fit_2time_corr_with_leasqr(CCfunc,fitfunc_str,fit_range,pin,dp, w);
                
                damHW(ii,jj) = fitres.pout(3);
                
                figure(100);
                clf;
                plot(CCfunc.time_1D(fit_range),CCfunc.CCNdtV(fit_range),'ob');
                hold on;
                plot(CCfunc.time_1D(fit_range),fitres.fitfunc,'r','LineWidth',3.0);
                title(['jj = ' num2str(jj) 'ii = ' num2str(ii) ' tau = ' num2str(damHW(ii,jj),'%e')]);
                pause(.1);
                %}
            else
                %%{
                for kk = 1:ndt
                    if (Cdti(kk) < 0.5);err1 = 0; break; end
                end
                if (err1==0)
                    isp = [(kk-1):kk];
                    dindex = [1:numel(ddam)];
                    %damHW(ii,jj) = interp1(Cdti(isp),dindex(isp),0.5);
                    damHW(ii,jj) = interp1(Cdti(isp),ddam(isp),0.5);
                    %damHW(ii,jj) = interp1(Cdti(isp),dtime(isp),0.5);
                end
                
                %%{
                if mod(ii,40) == 0
                    if mod(jj,40) == 0
                        figure(100);
                        clf;
                        plot(ddam,Cdti,'ob');
                        hold on;
                        plot(ddam,feval(fitfunc_str,ddam,[0 1.0 damHW(ii,jj)/log(2) 0]),'r','LineWidth',3.0);
                        title(['jj = ' num2str(jj) 'ii = ' num2str(ii) ' tau = ' num2str(damHW(ii,jj),'%e')]);
                        ylim([-1 1]);
                        pause();
                    end
                end
                %}
            end
            
           
        end
    end
   save([damHWfilename],'damHW');
else
    load(damHWfilename);
end
    
%% Plot section

% Plot correlation time in ML
% Scale from HW assuming exponential
tauML = damHW/log(2);

% Plot correlatin time in s
%tauML = tauML;

for kk = 1:size(runname,1)
    
    % Plot correlation time in ML
    % Scale from HW assuming exponential
     damHW_struct(kk).tauML = tauML;
     damHW_struct(kk).Cdt = Cdt;
     damHW_struct(kk).III = III_resize;
     
     
end


% central CTR is at nrow/2+1, ncol/2+1
ixcen = ncol/2 + 1;
iycen = nrow/2 + 1;
% DQ of CTRs is nsteps (pixels) in x (col)

QX = [1:ncol]-ixcen;
QY = [1:nrow]-iycen;

% QX original:
QX_orig = [1:ncol_orig]-ncol_orig/2;

% average ranges offsets and half-widths
ixahw = 56;
ixaoff = -nsteps/2;%5*nsteps;%3*nsteps;%nsteps;%nsteps/2;%0;%-nsteps/2;
iyahw = 15;
iyaoff = 0;%12;%4;%24;%16;%0;%
%iyaoff = 16;

yaoff_array = iyaoff;%[0];%[0 4 12 16 24];
xaoff_array = ixaoff;%[-nsteps/2];%[-nsteps/2 0 nsteps/2 nsteps 3*nsteps 5*nsteps];%[0 4 12 16 24];

color_multiple_offset = ['k' 'g' 'r' 'm' 'b'];


POSITION = [100 100 400 300];
PAPERPOSITION = [1 1 4 3];



% construct structure to pass all this information quickly to plotting
% functions:
struct_ranges.iycen = iycen;
struct_ranges.ixcen = ixcen;
struct_ranges.QX = QX;
struct_ranges.QY = QY;
struct_ranges.yaoff_array = yaoff_array;
struct_ranges.iyahw = iyahw;
struct_ranges.iyaoff = iyaoff;
struct_ranges.xaoff_array = xaoff_array;
struct_ranges.ixahw = ixahw;
struct_ranges.ixaoff = ixaoff;
struct_ranges.nsteps = nsteps;

 % Mark CTR positions
iCTR = [-5:4];


% Simple theory for tau(Q)
%tauMLth = nsteps./(2*pi*abs(QX-ixaoff));
% if time axis is in frames:
%{

tauMLth_Qx = .3e4*nsteps./(pi*abs(QX-ixaoff).^2); % Unclear factor
tauMLth_Qx_1 = .25e3*nsteps./(pi*abs(QX-ixaoff)); % Unclear factor
tauMLth_Qy = .08e4*nsteps./(pi*abs(QY).^2); % Unclear factor
tauMLth_Qy_3 = .1e5*nsteps./(pi*abs(QY).^3); % Unclear factor
tauMLth_Qy_2p5 = .03e5*nsteps./(pi*abs(QY).^2.5); % Unclear factor

%}

% if time axis is in seconds:
tauMLth_Qx = .5e10*nsteps./(pi*abs(QX-ixaoff).^2); % Unclear factor
tauMLth_Qx_1 = .25e9*nsteps./(pi*abs(QX-ixaoff)); % Unclear factor
tauMLth_Qy = .2e10*nsteps./(pi*abs(QY).^2); % Unclear factor
tauMLth_Qy_3 = .1e5*nsteps./(pi*abs(QY).^3); % Unclear factor
tauMLth_Qy_2p5 = .03e5*nsteps./(pi*abs(QY).^2.5); % Unclear factor


% Figure 1: tau
clear map2D;

figNum = 1;
map2D(1).tauML = tauML;
Function_display_sos.make_figure_2D([runname_title ' Correlation Time (sim. units) in log scale'],...
    map2D,'tauML',struct_ranges,iCTR,figNum,...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,'CLIM',[6 10],'ZSCALE','linear');

% Figure 1000: Contrast
figNum = 1000;
map2D(1).contrast = 10.^(Cdt(:,:,1));
Function_display_sos.make_figure_2D([runname_title ' Contrast '],...
    map2D,'contrast',struct_ranges,iCTR,figNum,...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,'CLIM',[.5 1.5],'ZSCALE','linear');


% Figure 12: diffracted intensity
figNum = 12;
% time frame to represent:
time_frame = 100;
map2D(1).III = III_resize(:,:,time_frame);
Function_display_sos.make_figure_2D([runname_title ' Diffracted intensity, time frame = ' num2str(time_frame)],...
    map2D,'III',struct_ranges,iCTR,figNum,...
    'CLIM',[-10 -2],'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,'plot_averages',0,'ZSCALE','linear');


% Figure 13: Ibar
figNum = 13;
% time frame to represent:
time_frame = 100;
map2D(1).Ibar = abs(Ibar(:,:,time_frame));
Function_display_sos.make_figure_2D([runname_title ' Ibar, time frame = ' num2str(time_frame)],...
    map2D,'Ibar',struct_ranges,iCTR,figNum,...
    'CLIM',[-10 -2],'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,'plot_averages',0,'ZSCALE','linear');


% average over rows
range_rows = iycen+iyaoff+[-iyahw:iyahw];

for kk = 1:size(runname,1)
    tauML_1D = mean(damHW_struct(kk).tauML(range_rows,:),1);
    damHW1D_struct(kk).tauML_1D_Qx = tauML_1D;

    Qx_offset = QX-ixaoff;
    Qx_positive = Qx_offset(Qx_offset>0);
    tauML_Qx_positive = tauML_1D(Qx_offset>0);
    tauML_th_QX_positive = tauMLth_Qx(Qx_offset>0);
    
    damHW1D_struct(kk).tauML_1D_Qx_positive = tauML_Qx_positive;
    damHW1D_struct(kk).tauML_th_QX_positive = tauML_th_QX_positive;
    
    Qx_negative = abs(Qx_offset(Qx_offset<=0));
    tauML_Qx_negative = tauML_1D(Qx_offset<=0);
    tauML_th_QX_negative = tauMLth_Qx(Qx_offset<=0);
    
    damHW1D_struct(kk).tauML_1D_Qx_negative = tauML_Qx_negative;
    damHW1D_struct(kk).tauML_th_QX_negative = tauML_th_QX_negative;
    
    damHW1D_struct(kk).Color = color_rates(kk);
    legend_str{kk} = [runname ' growth rate (ML/s) = ' num2str(growth_rates(kk))];
end


figNum = 2;

Function_display_sos.make_figures_tauML1D(['offset in Q_y = ' num2str(iyaoff)],...
    QX,damHW1D_struct,'tauML_1D_Qx',figNum,...
    'iCTR',iCTR,'nsteps',nsteps,...
    'XLABEL','Q_x (pixels)','YLIM',[1e5 1e10]);

legend(legend_str);


figNum = 20;

Function_display_sos.make_figures_tauML1D(['QX is positive,  offset in Qy = ' num2str(iyaoff)],...
    Qx_positive,damHW1D_struct,'tauML_1D_Qx_positive',...
    figNum,'tauMLth',tauML_th_QX_positive,'iCTR',iCTR,'nsteps',nsteps,'ixaoff',ixaoff,...
    'XLABEL','Q_x (pixels)',...
    'XSCALE','log','YSCALE',...
    'log','YLIM',[1e5 1e10]);

legend(legend_str);

figNum = 21;

Function_display_sos.make_figures_tauML1D(['QX is negative,  offset in Qy = ' num2str(iyaoff)],Qx_negative,damHW1D_struct,'tauML_1D_Qx_negative',...
    figNum,'tauMLth',tauML_th_QX_negative,'iCTR',iCTR,'nsteps',nsteps,'ixaoff',ixaoff,'XLABEL','Q_x (pixels)',...
    'XSCALE','log','YSCALE',...
    'log','YLIM',[1e5 1e10]);

legend(legend_str);


% average over columns
range_cols = ixcen+ixaoff+[-ixahw:ixahw];

for kk = 1:size(runname,1)
    tauML_1D = mean(damHW_struct(kk).tauML(:,range_cols),2);

    damHW1D_struct(kk).tauML_1D_Qy = tauML_1D;

    Qy_offset = QY-iyaoff;
    Qy_positive = Qy_offset(Qy_offset>0);
    tauML_Qy_positive = tauML_1D(Qy_offset>0);
    tauML_th_Qy_positive = tauMLth_Qy(Qy_offset>0);
    
    damHW1D_struct(kk).tauML_1D_Qy_positive = tauML_Qy_positive;
    damHW1D_struct(kk).tauML_th_QY_positive = tauML_th_Qy_positive;
    
    Qy_negative = abs(Qy_offset(Qy_offset<=0));
    tauML_Qy_negative = tauML_1D(Qy_offset<=0);
    tauML_th_Qy_negative = tauMLth_Qy(Qy_offset<=0);
    
    damHW1D_struct(kk).tauML_1D_Qy_negative = tauML_Qy_negative;
    damHW1D_struct(kk).tauML_th_QY_negative = tauML_th_Qy_negative;
    
    damHW1D_struct(kk).Color = color_rates(kk);
    legend_str{kk} = [runname_title ' growth rate (ML/s) = ' num2str(growth_rates(kk))];
end


figNum = 4;

Function_display_sos.make_figures_tauML1D(['offset in Q_x = ' num2str(ixaoff)],...
    QY,damHW1D_struct,'tauML_1D_Qy',figNum,...
    'XLABEL','Q_y (pixels)','YLIM',[1e5 1e10]);

legend(legend_str);


figNum = 40;

Function_display_sos.make_figures_tauML1D(['QY is positive,  offset in Qx = ' num2str(ixaoff)],Qy_positive,damHW1D_struct,'tauML_1D_Qy_positive',...
    figNum,'tauMLth',tauML_th_Qy_positive,'iyaoff',iyaoff,'XLABEL','Q_y (pixels)',...
    'XSCALE','log','YSCALE',...
    'log','YLIM',[1e5 1e10]);

legend(legend_str);

figNum = 41;

Function_display_sos.make_figures_tauML1D(['QY is negative,  offset in Qx = ' num2str(ixaoff)],Qy_negative,damHW1D_struct,'tauML_1D_Qy_negative',...
    figNum,'tauMLth',tauML_th_Qy_negative,'iyaoff',iyaoff,'XLABEL','Q_y (pixels)',...
    'XSCALE','log','YSCALE',...
    'log','YLIM',[1e5 1e10]);
legend(legend_str);


return;

%% Average over diagonals or other directions:


% average over rows

alpha = 45;

% points through which the line has to pass
y_star = nrow/2;
x_star = ncol/2;

a = tand(alpha);
b = y_star - a*x_star;
Delta = 10;
b_array = b + [Delta  -Delta]./cosd(alpha);

alpha2 = 90 - alpha;
a2 = -tand(alpha2);
b2 = y_star-a2*x_star;
Delta2 = 60;
b2_array = b2+ [ Delta2 -Delta2]./cosd(alpha2);%[-60 60];
% 
% for ooo = 1:numel(b_array)
%     b2(ooo) = (a-a2)*(ncol/2)+b_array(ooo)+b2_array(ooo)/cosd(alpha);
% end

upper_limit1 = round(a).*[1:ncol]+b_array(1);
lower_limit1 = round(a).*[1:ncol]+b_array(2);

upper_limit2 = round(a2).*[1:ncol]+ b2_array(1);
lower_limit2 = round(a2).*[1:ncol]+ b2_array(2);

[X_index_grid, Y_index_grid_flip] = meshgrid([1:nrow],[1:ncol]);

Y_index_grid = flipud(Y_index_grid_flip);

mask1 = (Y_index_grid >lower_limit1 &Y_index_grid < upper_limit1 );
mask2 = ( Y_index_grid >lower_limit2 &Y_index_grid < upper_limit2 );

total_mask = mask1.*mask2;

figure(1000);clf;
subplot(131);imagesc(mask1);axis image;title('mask1');
subplot(132);imagesc(mask2);axis image;title('mask2');
subplot(133);imagesc(mask1.*mask2);axis image;

figure(1001);clf;
plot(lower_limit1,'r');
hold on;
plot(lower_limit2,'b');
legend('lower_limit1','lower_limit2');
axis image;

figure(700);clf;
imagesc(log10(damHW_struct(1).tauML.*total_mask));colorbar;

for kk = 1:size(runname,1)
    % select a strip of 2D map
    tauML_2Dstrip = damHW_struct(kk).tauML.*mask1.*mask2;
    
    % calculate manually mean:
    
    tau_mean = zeros(ncol,1);
    offset_range = zeros(ncol,1);
    for jjj = 1:ncol
        if sum(total_mask(:,jjj)) ~= 0
           offset_range(jjj) = (a-a2)*jjj+mean(b_array);
           for ooo = 1:numel(b_array)
              x_lim(ooo) = (b_array(ooo)-offset_range(jjj))/(a2-a);
              y_lim(ooo) = a*x_lim(ooo)+b_array(ooo);
           end
           
           x_range = [x_lim(1):x_lim(2)];
           y_range = [y_lim(2):y_lim(1)];
           
           %{
           if x_lim(2)>x_lim(1)
               x_range = [x_lim(1):x_lim(2)];
               y_range = [y_lim(2):y_lim(1)];
           else
               x_range = [x_lim(2):x_lim(1)];
               y_range = [y_lim(1):y_lim(2)];
           end
           %}
           
           for llll =1:numel(x_range)
               tau_mean(jjj) = tau_mean(jjj) + tauML_2Dstrip(y_range(llll),x_range(llll))/numel(x_range);
           end
            
        end
        
           
           Delta_Q_diag(jjj) = QY(jjj)/sind(alpha);
    end
    
    non_zero_index = find(mask1.*mask2>0);
    
    
    
    tauML_1D =  tau_mean;%mean(damHW_struct(kk).tauML.*mask1.*mask2,2);

    damHW1D_struct(kk).tauML_1D_Qy = tauML_1D;
    Qy_offset = Delta_Q_diag;%QY-iyaoff;
    
    Qy_positive = Qy_offset(Qy_offset>0);
    tauML_Qy_positive = tauML_1D(Qy_offset>0);
    tauML_th_Qy_positive = tauMLth_Qy(Qy_offset>0);
    
    damHW1D_struct(kk).tauML_1D_Qy_positive = tauML_Qy_positive;
    damHW1D_struct(kk).tauML_th_QY_positive = tauML_th_Qy_positive;
    
    Qy_negative = abs(Qy_offset(Qy_offset<=0));
    tauML_Qy_negative = tauML_1D(Qy_offset<=0);
    tauML_th_Qy_negative = tauMLth_Qy(Qy_offset<=0);
    
    damHW1D_struct(kk).tauML_1D_Qy_negative = tauML_Qy_negative;
    damHW1D_struct(kk).tauML_th_QY_negative = tauML_th_Qy_negative;
    
    damHW1D_struct(kk).Color = color_rates(kk);
    legend_str{kk} = [runname_title ' growth rate (ML/s) = ' num2str(growth_rates(kk))];
end


figNum = 5;

Function_display_sos.make_figures_tauML1D(['offset in Q_x = ' num2str(ixaoff)],...
    Qy_offset,damHW1D_struct,'tauML_1D_Qy',figNum,...
    'XLABEL','Q_y (pixels)','YLIM',[1 1e10]);

legend(legend_str);


figNum = 50;

Function_display_sos.make_figures_tauML1D(['QY is positive,  offset in Qx = ' num2str(ixaoff)],Qy_positive,damHW1D_struct,'tauML_1D_Qy_positive',...
    figNum,'tauMLth',tauML_th_Qy_positive,'iyaoff',iyaoff,'XLABEL','Q_y (pixels)',...
    'XSCALE','log','YSCALE',...
    'log','YLIM',[1e5 1e10]);

legend(legend_str);

figNum = 51;

Function_display_sos.make_figures_tauML1D(['QY is negative,  offset in Qx = ' num2str(ixaoff)],Qy_negative,damHW1D_struct,'tauML_1D_Qy_negative',...
    figNum,'tauMLth',tauML_th_Qy_negative,'iyaoff',iyaoff,'XLABEL','Q_y (pixels)',...
    'XSCALE','log','YSCALE',...
    'log','YLIM',[1e5 1e10]);
legend(legend_str);


return;

figNum = 5;

Function_display_sos.make_figures_tauML1D(['offset in Q_x = ' num2str(ixaoff)],...
     Q_diag ,damHW1D_struct,'tauML_1D_diag',figNum,...
    'XLABEL','Q\_diagonal (pixels)','YLIM',[1e5 1e9]);



return;


% in log scale:

figNum = 20;
Qx_offset = QX-ixaoff;
Qx_positive = Qx_offset(Qx_offset>0);
tauML_Qx_positive = tauML_1D(Qx_offset>0);
tauML_th_QX_positive = tauMLth_Qx(Qx_offset>0);
Function_display_sos.make_figures_tauML1D([runname_title ' QX is positive, offset in y = ' num2str(iyaoff)],Qx_positive,tauML_Qx_positive,figNum,...
    'tauMLth',tauML_th_QX_positive,'iCTR',iCTR,'nsteps',nsteps,'ixaoff',ixaoff,'XLABEL','Q_x (pixels)','XSCALE','log','YSCALE',...
    'log','YLIM',[1e6 3e9]);

figNum = 21;
Qx_negative = abs(Qx_offset(Qx_offset<=0));
tauML_Qx_negative = tauML_1D(Qx_offset<=0);
tauML_th_QX_negative = tauMLth_Qx(Qx_offset<=0);
Function_display_sos.make_figures_tauML1D([runname_title ' QX is negative, offset in y = ' num2str(iyaoff)],Qx_negative,tauML_Qx_negative,figNum,...
    'tauMLth',tauML_th_QX_negative,'iCTR',iCTR,'nsteps',nsteps,'ixaoff',ixaoff,'XLABEL','Q_x (pixels)','XSCALE','log','YSCALE',...
    'log','YLIM',[1e6 3e9]);



figNum = 4;

% average over columns
range_cols = ixcen+ixaoff+[-ixahw:ixahw];
tauML_1D = mean(tauML(:,range_cols),2);
Function_display_sos.make_figures_tauML1D([runname_title],QY,tauML_1D,figNum,...
    'XLABEL','Q_y (pixels)','YLIM',[1e6 3e9]);

% in log scale:

figNum = 40;
Qy_offset = QY-iyaoff;
Qy_positive = Qy_offset(Qy_offset>0);
tauML_Qy_positive = tauML_1D(Qy_offset>0);
tauML_th_Qy_positive = tauMLth_Qy(Qy_offset>0);
Function_display_sos.make_figures_tauML1D([runname_title ' QY is positive, log scale'],...
    Qy_positive,tauML_Qy_positive,figNum,...
    'tauMLth',tauML_th_Qy_positive,'XLABEL','Q_y (pixels)','XSCALE','log','YSCALE',...
    'log','YLIM',[1e6 3e9]);

figNum = 41;
Qy_negative = abs(Qy_offset(Qy_offset<=0));
tauML_Qy_negative = tauML_1D(Qy_offset<=0);
tauML_th_Qy_negative = tauMLth_Qy(Qy_offset<=0);
Function_display_sos.make_figures_tauML1D([runname_title ' QY is negative, log scale'],Qy_negative,tauML_Qy_negative,figNum,...
    'tauMLth',tauML_th_Qy_negative,'XLABEL','Q_y (pixels)','XSCALE','log','YSCALE',...
    'log','YLIM',[1e6 3e9]);








return;


figure(2);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(QX,mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o','Color','r');
pa = axis;
hl = line(QX,tauMLth_Qx);
set(hl,'LineStyle','-','Color','m','LineWidth',2.0);
axis(pa);
% Show CTR positions
for ii = iCTR
    hl = line(nsteps*ii*[1 1],pa(3:4));
    set(hl,'LineStyle','--','Color','r');
end
xlabel('Q_X (pixels)');
%ylabel('Mean Correlation Time (ML)');
ylabel('Mean Correlation Time (s)');
title([runname_title]);
legend(['Offset Q_y = ' num2str(iyaoff)]);





% plot log(tau) vs log(Q_x-Q_1)

figure(3);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
%hl = line(abs(QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
hl = line((QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
%[legend_str] = Function_display_sos.show_multiple_offsets(QX,QY,yaoff_array,iycen,iyahw,ixaoff,iycen,ixahw,tauML,color_multiple_offset)
set(hl,'LineStyle','none','Marker','o');
set(gca,'Xscale','log','Yscale','log');
set(gca,'Ylim',[1e6 1e10]);
set(gca,'Xlim',[1 200]);
pa = axis;
hl = line(abs(QX-ixaoff),tauMLth_Qx);
set(hl,'LineStyle','-','Color','m','LineWidth',2.0);
axis(pa);
hl = line(abs(QX-ixaoff),tauMLth_Qx_1);
set(hl,'LineStyle','-','Color','k','LineWidth',2.0);
axis(pa);
% Show CTR positions
for ii = iCTR
    hl = line(abs(nsteps*ii*[1 1]-ixaoff),pa(3:4));
    set(hl,'LineStyle','--','Color','r');
end
xlabel('|Q_X - Q_1| (pixels)');
ylabel('Mean Correlation Time (ML)');
%ylabel('Mean Correlation Time (s)');
title([runname_title ' offset in Q_y = ' num2str(iyaoff) ]);
legend('128 x 512','comparison with 1/Q_x^2','comparison with 1/Q_x ');
%legend(legend_str);


figure(4);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(QY,mean(tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o');
hold on;
hl2 = line(QY,tauMLth_Qy);
set(hl2,'LineStyle','-','Color','m');
xlabel('Q_Y (pixels)');
%ylabel('Mean Correlation Time (ML)');
ylabel('Mean Correlation Time (s)');
title([runname_title]);

figure(5);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(abs(QY),mean(tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o');
hold on;
hl2 = line(abs(QY),tauMLth_Qy);
set(hl2,'LineStyle','-','Color','m','LineWidth',2.0);
set(gca,'Xscale','log','Yscale','log');
xlabel('Q_Y (pixels)');
%ylabel('Mean Correlation Time (ML)');
ylabel('Mean Correlation Time (s)');
title([runname_title ' Offset Q_x = ' num2str(ixaoff)]);
legend(['128 x 512'],'comparison with 1/Q_y^2');
%legend(['Offset Q_y = ' num2str(iyaoff)]);


figure(100);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
pcolor(QX,QY,ihm(:,:,time_frame));
shading flat;
xlabel('Q_X (pixels)');
ylabel('Q_Y (pixels)');
colorbar;
%title([runname_title ' Correlation Time (ML) in log scale']);
%title([runname_title ' Correlation Time (sim. unit time) in log scale']);
set(gca,'FontSize',15);
axis image;
title(['cropped surface at time frame = ' num2str(time_frame)]);

figure(101);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
pcolor(QX_orig,QY,ihm_orig(:,:,time_frame));
shading flat;
xlabel('Q_X (pixels)');
ylabel('Q_Y (pixels)');
colorbar;
%title([runname_title ' Correlation Time (ML) in log scale']);
%title([runname_title ' Correlation Time (sim. unit time) in log scale']);
set(gca,'FontSize',15);
axis image;
title(['total surface at time frame = ' num2str(time_frame)]);

% Irene's modification plot log(tau) vs log(Q_x-Q_1)

%{
figure(6);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(abs(QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
hold on;
%QX_recentered = QX-ixaoff;
%index_pos = find(QX_recentered>0);
%hl = line((QX_recentered(index_pos)),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],index_pos),1));
set(hl,'LineStyle','none','Marker','o');
set(gca,'Xscale','log','Yscale','log');
pa = axis;
hl = line(abs(QX-ixaoff),tauMLth_Qx);
%hl = line(QX_recentered(index_pos),tauMLth((index_pos)));
set(hl,'LineStyle','-','Color','m','LineWidth',2.0);
axis(pa);
% Show CTR positions
for ii = iCTR
    hl = line(abs(nsteps*ii*[1 1]-ixaoff),pa(3:4));
    set(hl,'LineStyle','--','Color','r');
end
xlabel('Q_X (pixels)');
%xlabel('|Q_X - Q_1| (pixels)');
%ylabel('Mean Correlation Time (ML)');
ylabel('Mean Correlation Time (s)');
title([runname_title]);

%}


%%%%%%%%% for the crystal growth conference presentation
%{
path_save_fig = '/Users/ialmazn/Box Sync/XPCS_sputtering_ZnO_TiO2_2017_2019/fig_growth_conf/';

figure(100);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
pcolor(QX,QY,log10(tauML));
shading flat;
xlabel('Q_X (pixels)');
ylabel('Q_Y (pixels)');
colorbar;
%title([runname_title ' Correlation Time (ML) in log scale']);
%title([runname_title ' Correlation Time (sim. unit time) in log scale']);
set(gca,'FontSize',15);
axis image;

print('-r600',figure(100),[path_save_fig 'timedecay_coh_simulation'],'-dpng');

figure(7);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
pcolor((ihm(:,:,1)));
shading flat;
colormap jet;
axis off;

print('-r600',figure(7),[path_save_fig 'steps_simulation_1'],'-dpng');





figure(10);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
pcolor(QX,QY,log10(III(:,:,370)));
%pcolor(QX,QY,(III(:,:,370)));
hold on;
set(gca,'Clim',[2 6])
% Mark CTR positions
iCTR = [-5:4];
hl = line(nsteps*iCTR,zeros(size(iCTR)));
set(hl,'Linestyle','none','Marker','x','MarkerSize',0.4,'Color','y');
shading flat;
axis off;
colorbar;
title('diffracted intensity');

print('-r600',figure(10),[path_save_fig 'det_coh_simulation_370'],'-dpng');

%{
figure(13);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
%hl = line(abs(QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
hl = line((QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o','Color','r');
set(gca,'Xscale','log','Yscale','log');
pa = axis;
hl = line(abs(QX-ixaoff),tauMLth_Qx);
set(hl,'LineStyle','-','Color','r','LineWidth',2.0);
axis(pa);
ylim([ 1    5e2]);
xlabel('Q_x (pixels)')
ylabel('\tau (frames)');
set(gca,'FontSize',15);

print('-r600',figure(13),[path_save_fig 'tauQx_simulation_Q2'],'-dpng');

figure(131);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
%hl = line(abs(QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
hl = line((QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o','Color','r');
set(gca,'Xscale','log','Yscale','log');
pa = axis;
%hl = line(abs(QX-ixaoff),tauMLth_Qx);
%set(hl,'LineStyle','-','Color','r','LineWidth',2.0);
axis(pa);
ylim([ 1    5e2]);
xlabel('Q_x (pixels)')
ylabel('\tau (frames)');
set(gca,'FontSize',15);

print('-r600',figure(13),[path_save_fig 'tauQx_simulation_3regions'],'-dpng');
%}




figure(15);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(abs(QY),mean(tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o','Color','k');
hold on;
hl2 = line(abs(QY),tauMLth_Qy);
set(hl2,'LineStyle','-','Color','k','LineWidth',2.0);
hold on;
hl3 = line(abs(QY),tauMLth_Qy_3)
set(hl3,'LineStyle','-','Color','b','LineWidth',2.0);
set(gca,'Xscale','log','Yscale','log');
xlabel('Q_Y (pixels)');
%ylabel('Mean Correlation Time (ML)');
ylabel('\tau (sim. unit) ');
ylim([ 1    5e2]);
set(gca,'FontSize',15);
%title([runname_title]);
%legend(['offset Q_x =' num2str(ixaoff)]);

print('-r600',figure(15),[path_save_fig 'tauQy_simulation_Q2'],'-dpng');



figure(16);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
%hl = line(abs(QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
hl = line((QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o','Color','r');
set(gca,'Xscale','log','Yscale','log');
pa = axis;
hl = line(abs(QX-ixaoff),tauMLth_Qx);
set(hl,'LineStyle','-','Color','r','LineWidth',2.0);
hold on;
hl2 = line(abs(QY),mean(tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl2,'LineStyle','none','Marker','o','Color','k');
hold on;
hl3 = line(abs(QY),tauMLth_Qy_2p5)
set(hl3,'LineStyle','-','Color','k','LineWidth',2.0);
axis(pa);
ylim([ 1    5e2]);
xlabel('Q_x and Qy (pixels)')
ylabel('\tau (frames)');
set(gca,'FontSize',15);

print('-r600',figure(16),[path_save_fig 'tauQxandQy_simulation'],'-dpng');





figure(1000);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
pcolor(log10(tauML));
shading flat;
xlabel('Q_X (pixels)');
ylabel('Q_Y (pixels)');
%colorbar;
%title([runname_title ' Correlation Time (ML) in log scale']);
title([runname_title ' map to choose interesting pixels']);

ii_row = 90;
jj_col = 90;

Cdti_toshow = squeeze(Cdt(ii_row,jj_col,:));
Cdti_toshow = Cdti_toshow/Cdti_toshow(1);

%figure(16);
%clf;
%plot(ddam,Cdti_toshow);
%xlabel('Time Delta (sim. unit)');
%ylabel('Correlation');

fit_range = [1:1:round(1/2.5*length(Cdti_toshow))];
fitfunc_str = 'FittingFunctions.CCN2single_fit';
pin = [0 1 1 0];%pin_iiT(iT,:);
dp =  [[0 1 1 0]*0.0001];%dp_iiT(iT,:);
w = ones(length(fit_range),1);

CCfunc.time_1D = ddam';
CCfunc.CCNdtV = Cdti_toshow;

figure;
[ fitres] = FittingFunctions.fit_2time_corr_with_leasqr(CCfunc,fitfunc_str,fit_range,pin,dp, w);

damHW_fit = fitres.pout(3);

figure(18);
clf;
%plot(CCfunc.time_1D(fit_range),CCfunc.CCNdtV(fit_range),'ob');
plot(CCfunc.CCNdtV(fit_range),'or','LineWidth',3.0);
hold on;
plot(fitres.fitfunc,'r','LineWidth',3.0);
%plot(CCfunc.time_1D(fit_range),fitres.fitfunc,'r','LineWidth',3.0);
%title(['col jj = ' num2str(jj_col) 'row ii = ' num2str(ii_row) 'tau = ' num2str(damHW_fit)] );
ylim([-1e-1 1]);
xlim([0 100])
xlabel('Time Delta (frames)');
ylabel('Correlation');
set(gca,'FontSize',20);
set(gca,'Position',[ 0.1622    0.1590    0.7428    0.7660]);


print('-r600',figure(18),[path_save_fig 'exp_decay_fast'],'-dpng');


ii_row = 76;
jj_col = 88;

Cdti_toshow = squeeze(Cdt(ii_row,jj_col,:));
Cdti_toshow = Cdti_toshow/Cdti_toshow(1);

%figure(16);
%clf;
%plot(ddam,Cdti_toshow);
%xlabel('Time Delta (sim. unit)');
%ylabel('Correlation');

fit_range = [1:1:round(1/2.5*length(Cdti_toshow))];
fitfunc_str = 'FittingFunctions.CCN2single_fit';
pin = [0 1 1 0];%pin_iiT(iT,:);
dp =  [[0 1 1 0]*0.0001];%dp_iiT(iT,:);
w = ones(length(fit_range),1);

CCfunc.time_1D = ddam';
CCfunc.CCNdtV = Cdti_toshow;

figure;
[ fitres] = FittingFunctions.fit_2time_corr_with_leasqr(CCfunc,fitfunc_str,fit_range,pin,dp, w);

damHW_fit = fitres.pout(3);

figure(18);
hold on;
plot(CCfunc.CCNdtV(fit_range),'og','LineWidth',3.0);
hold on;
plot(fitres.fitfunc,'g','LineWidth',3.0);
set(gca,'Position',[ 0.1622    0.1590    0.7428    0.7660]);

print('-r600',figure(18),[path_save_fig 'exp_decay_fast_and_slow'],'-dpng');


savefig(figure(18),'./figures_conference/exp_decay_107_90.fig');


figure(2000);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
pcolor(log10(III(:,:,round(end/2))));
shading flat;
xlabel('Q_X (pixels)');
ylabel('Q_Y (pixels)');
%colorbar;
%title([runname_title ' Correlation Time (ML) in log scale']);
title([runname_title ' Intensity']);
%}