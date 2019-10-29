% Script to plot structure factor from sosvpe simulations
% From sosvpe_sf_2t_new2
% Dropped azimuthal averages
% Added Cdt calc
% 22-FEB-19 GBS

skip = 1;
path = 'miscut_8_14/';
runname = [path 'miscut_8_14'];
runname_title = 'miscut\_8\_14';
namefile = [runname '_corr_dt_L0p5.mat'];

addpath(genpath(['/Users/ialmazn/Box' ' Sync/XPCS_sputtering_ZnO_TiO2_2017_2019/XPCS_analysis_on_thefly']));


%make_vid(runname,[0 10],'jet');
%make_vid_sf(runname,[0 100],'jet');

if ~skip
    
    load([runname '_stats.mat']);

    if ~exist('nsteps','var'); nsteps = 0; end

    if ~exist('damono','var')
        havg = squeeze(mean(mean(ihm))); % Average height, gives growth amount
        damono = havg - 1 - nsteps/2;
    end

    
     % simulation of probe:
    [Xgrid,Ygrid] = meshgrid([1:ncol],[1:nrow]);
    FWHM_X = 2e-4;
    probe = exp(-FWHM_X.*(Xgrid-round(ncol/2)).^2).*ones(nrow,ncol);
    
    probe_fft = fftshift(fftn(fftshift(probe)));

    % Calculate for fixed L
    L = 0.3;
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
        %ftm_noprobe(:,:,ii) = fftshift(ifft2(ph)) +ftms;
        
         %ftm(:,:,ii) = conv2(ftm_noprobe(:,:,ii),probe_fft,'same');
        
        %{
        subplot(211);
        imagesc( abs(ph));
        axis image;
        colorbar;
        
%         subplot(222);
%         imagesc( angle(ph));
%         axis image;
%         colorbar;
        
        subplot(212);
        imagesc(log10(abs(ftm(:,:,ii))));
        axis image;
        colorbar;
        
%         subplot(224);
%         imagesc(angle(ftm(:,:,ii)));
%         axis image;
%         colorbar;
        title(['time = ' num2str(ii)]);
        
        pause(.1)
        %}
    end
    
    
   
    

    % Calc delta-t correlations for each pixel
    
    % Do not use first dt_minML of growth for time correlations
    dt_minML = 2;%if in ML
    idt = damono > dt_minML;
    
    III = abs(ftm(:,:,idt)).^2;
    ndt = size(III,3);
    ddam = damono(idt);
    ddam = ddam - ddam(1); % delta time coord (ML) (assumes evenly spaced)
    
    % Use time average as Ibar
    
    Ibar = mean(III,3);
    
    dlnI = III./Ibar - 1;

    Cdt = NaN*ones(size(dlnI));
    I1 = ones(ndt,1);
    Idt0 = conv(I1,I1);
    for ii = 1:nrow
        for jj = 1:ncol
            It = squeeze(dlnI(ii,jj,:));
            Idt = conv(It,flip(It))./Idt0;
            Cdt(ii,jj,:) = Idt(end-ndt+1:end);
        end
    end
    
    save([namefile]);

else
    load([namefile]);
end

xp = ncol*[0 0.2 0.2 0];
yp = nrow*[0 0 0.1 0.1];
xtx = ncol/40;
ytx = nrow/20;
textclr = 'k';

imovie = 0;
if imovie % Write movie of speckle
    writerObj = VideoWriter([runname '.avi']);
    open(writerObj);

    figure
    set(gcf,'Position',[600 200 400 400]);
    set(gcf,'PaperPosition',[1 1 4 4]);
    axes('Box','on');
    imagesc(log(abs(ftm(:,:,1)).^2),[-16 -6]);
    axis image
    shading flat;
    set(gca,'Visible','off');
    patch(xp,yp,'w');
    ht = text(xtx,ytx,[num2str(damono(ii),'%5.2f'),' ML']);
    set(ht,'Color',textclr,'FontSize',12);
    pause(0.05);
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
    hc = get(gca,'Children');
    for ii = 2:nt
        %imagesc(log(abs(ftm(:,:,ii)).^2),[2 12]);
        set(hc(end),'CData',log(abs(ftm(:,:,ii)).^2));
        %axis image
        %shading flat;
        %set(gca,'Visible','off');
        patch(xp,yp,'w');
        ht = text(xtx,ytx,[num2str(damono(ii),'%5.2f'),' ML']);
        set(ht,'Color',textclr,'FontSize',12);
        pause(0.05);
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
    close(writerObj);

end

% get time HW of Cdt for each pixel, units of ML

damHW = ddam(end)*ones(nrow,ncol);

skip_HWcalc = 0;
if  ~skip_HWcalc  
    for ii = 1:nrow
        for jj = 1:ncol
            Cdti = squeeze(Cdt(ii,jj,:));
            Cdti = Cdti/Cdti(1);
            err1 = 1;
            
            %{
            fit_range = [1:1:round(1/2.5*length(Cdti))];
            fitfunc_str = 'FittingFunctions.CCN2single_fit';
            pin = [0 1 1 0];%pin_iiT(iT,:);
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
            title(['jj = ' num2str(jj) 'ii = ' num2str(ii)]);
            pause(.1);
            %}
            
            for kk = 1:ndt
               if (Cdti(kk) < 0.5);err1 = 0; break; end
            end
            if (err1==0)
               isp = [(kk-1):kk];
               dindex = [1:numel(ddam)];
               %damHW(ii,jj) = interp1(Cdti(isp),dindex(isp),0.5);
               %damHW(ii,jj) = interp1(Cdti(isp),ddam(isp),0.5);
               damHW(ii,jj) = interp1(Cdti(isp),dtime(isp),0.5);
            end
        end
    end
else
    load('damHW_fit.mat');
end

% central CTR is at nrow/2+1, ncol/2+1
ixcen = ncol/2 + 1;
iycen = nrow/2 + 1;
% DQ of CTRs is nsteps (pixels) in x (col)

QX = [1:ncol]-ixcen;
QY = [1:nrow]-iycen;

% average ranges offsets and half-widths
ixahw = 0;
ixaoff = -nsteps/2;
iyahw = 0;
iyaoff = 0;
%iyaoff = 16;

POSITION = [100 100 400 300];
PAPERPOSITION = [1 1 4 3];

% Plot correlation time in ML
% Scale from HW assuming exponential
tauML = damHW/log(2);

% Plot correlation time in s
%tauML = tauML;
yaoff_array = [0];%[0 4 12 16 24];%
xaoff_array = [-nsteps/2 ];%[-nsteps/2 0 nsteps/2 nsteps 3*nsteps 5*nsteps];%[-nsteps/2];%[0 4 12 16 24];


% Plot correlatin time in s
%tauML = tauML;

figure(1);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
pcolor(QX,QY,log10(tauML));
set(gca,'Clim',[6.1971    9.3297]);
shading flat;
% Show ranges averaged
hl = line(QX,(iyaoff+iyahw)*ones(size(QY)));
set(hl,'Linestyle','--','Color','w');
hl = line(QX,(iyaoff-iyahw)*ones(size(QY)));
set(hl,'Linestyle','--','Color','w');
hl = line((ixaoff+ixahw)*ones(size(QX)),QY);
set(hl,'Linestyle','--','Color','w');
hl = line((ixaoff-ixahw)*ones(size(QX)),QY);
set(hl,'Linestyle','--','Color','w');
% Mark CTR positions
iCTR = [-5:4];
hl = line(nsteps*iCTR,zeros(size(iCTR)));
set(hl,'Linestyle','none','Marker','x','Color','r');
xlabel('Q_X (pixels)');
ylabel('Q_Y (pixels)');
colorbar;
%title([runname_title ' Correlation Time (ML) in log scale']);
title([runname_title ' Correlation Time (s) in log scale at L = ' num2str(L)]);

figure(10);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
pcolor(QX,QY,log10(III(:,:,1)));
%set(gca,'Clim',[2 6]);
set(gca,'Clim',[-6 -2]);
shading flat;
% Show ranges averaged
hl = line(QX,(iyaoff+iyahw)*ones(size(QY)));
set(hl,'Linestyle','--','Color','w');
hl = line(QX,(iyaoff-iyahw)*ones(size(QY)));
set(hl,'Linestyle','--','Color','w');
hl = line((ixaoff+ixahw)*ones(size(QX)),QY);
set(hl,'Linestyle','--','Color','w');
hl = line((ixaoff-ixahw)*ones(size(QX)),QY);
set(hl,'Linestyle','--','Color','w');
% Mark CTR positions
iCTR = [-5:4];
hl = line(nsteps*iCTR,zeros(size(iCTR)));
set(hl,'Linestyle','none','Marker','x','Color','r');
xlabel('Q_X (pixels)');
ylabel('Q_Y (pixels)');
colorbar;
%title([runname_title ' Correlation Time (ML)']);
title([runname_title ' Diffracted intensity, time frame = 600, at L = ' num2str(L)]);

% Simple theory for tau(Q)
%tauMLth = nsteps./(2*pi*abs(QX-ixaoff));

% if time axis is in frames:
%{
tauMLth_Qx = .3e4*nsteps./(pi*abs(QX-ixaoff).^2); % Unclear factor
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
set(hl,'LineStyle','none','Marker','o');
set(gca,'Xscale','log','Yscale','log');
pa = axis;
hl = line(abs(QX-ixaoff),tauMLth_Qx);
set(hl,'LineStyle','-','Color','m','LineWidth',2.0);
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
title([runname_title ' at L = ' num2str(L)]);
legend('128 x 128','comparison with 1/Q_x^2','comparison with 1/Q_x');



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
title([runname_title ' at L = ' num2str(L)]);
%legend(['offset Q_x =' num2str(ixaoff)],'comparison with 1/Q_y^2');
legend(['128 x 128' ],'comparison with 1/Q_y^2');



return;

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







% Figure 1:
figNum = 1;
tauML_struct.tauML = tauML;
Function_display_sos.make_figure_2D([runname_title ' Correlation Time (sim. units) in log scale'],...
    tauML_struct,'tauML',struct_ranges,iCTR,figNum,'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION);

% Figure 10:
figNum = 10;
tauML_struct.III = III(:,:,time_frame);

Function_display_sos.make_figure_2D([runname_title ' Diffracted intensity, time frame = ' num2str(time_frame)],...
  tauML_struct,'III',struct_ranges,iCTR,figNum,'CLIM',[-6 -2],'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION);

% Figures 2 and 4 (averaged tau_ML over colums or rows):
figNum = 2;

% average over rows
range_rows = iycen+iyaoff+[-iyahw:iyahw];
tauML_1D = mean(tauML(range_rows,:),1);
Function_display_sos.make_figures_tauML1D([runname_title],QX,tauML_1D,figNum,'iCTR',iCTR,'nsteps',nsteps,...
    'XLABEL','Q_x (pixels)','YLIM',[1e6 3e9]);

% in log scale:

figNum = 20;
Qx_offset = QX-ixaoff;
Qx_positive = Qx_offset(Qx_offset>0);
tauML_Qx_positive = tauML_1D(Qx_offset>0);
tauML_th_QX_positive = tauMLth_Qx(Qx_offset>0);
Function_display_sos.make_figures_tauML1D([runname_title ' QX is positive, log scale'],Qx_positive,tauML_Qx_positive,figNum,...
    'tauMLth',tauML_th_QX_positive,'iCTR',iCTR,'nsteps',nsteps,'ixaoff',ixaoff,'XLABEL','Q_x (pixels)','XSCALE','log','YSCALE',...
    'log','YLIM',[1e6 3e9]);

figNum = 21;
Qx_negative = abs(Qx_offset(Qx_offset<=0));
tauML_Qx_negative = tauML_1D(Qx_offset<=0);
tauML_th_QX_negative = tauMLth_Qx(Qx_offset<=0);
Function_display_sos.make_figures_tauML1D([runname_title ' QX is negative, log scale'],Qx_negative,tauML_Qx_negative,figNum,...
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
set(gca,'Clim',[-6 -2])
% Mark CTR positions
iCTR = [-5:4];
hl = line(nsteps*iCTR,zeros(size(iCTR)));
set(hl,'Linestyle','none','Marker','x','MarkerSize',0.4,'Color','y');
shading flat;
axis off;
colorbar;

print('-r600',figure(10),[path_save_fig 'det_coh_simulation_370'],'-dpng');

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