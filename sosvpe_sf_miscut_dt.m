% Script to plot structure factor from sosvpe simulations
% From sosvpe_sf_2t_new2
% Dropped azimuthal averages
% Added Cdt calc
% 22-FEB-19 GBS

skip = 0;
runname = 'miscut_8_13';
runname_title = 'miscut\_8\_13';

if ~skip
    
    load([runname '_stats.mat']);

    if ~exist('nsteps','var'); nsteps = 0; end

    if ~exist('damono','var')
        havg = squeeze(mean(mean(ihm))); % Average height, gives growth amount
        damono = havg - 1 - nsteps/2;
    end

% Calculate for fixed L
    L = 0.5;

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
    
    for ii = 1:nt
        ph = xLsub.*(1 - xL.^(ihm(:,:,ii) - ihsub))./(1 - xL);

        % Calculate Fourier transform of phased height
        % It matters here whether to do fft2 or ifft2; 
        % ifft2 does exp(+iqr), and divides by nrow * ncol
        % also add substrate contribution to CTR positions here 
        ftm(:,:,ii) = fftshift(ifft2(ph)) + ftms;
    end

    % Calc delta-t correlations for each pixel
    
    % Do not use first dt_minML of growth for time correlations
    dt_minML = 2;
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
    
    save([runname '_corr_dt.mat']);

else
    load([runname '_corr_dt.mat']);
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
for ii = 1:nrow
    for jj = 1:ncol
        Cdti = squeeze(Cdt(ii,jj,:));
        Cdti = Cdti/Cdti(1);
        err1 = 1;
        for kk = 1:ndt
            if (Cdti(kk) < 0.5);err1 = 0; break; end
        end
        if (err1==0) 
            isp = [(kk-1):kk];
            damHW(ii,jj) = interp1(Cdti(isp),ddam(isp),0.5);
        end
    end
end

% central CTR is at nrow/2+1, ncol/2+1
ixcen = ncol/2 + 1;
iycen = nrow/2 + 1;
% DQ of CTRs is nsteps (pixels) in x (col)

QX = [1:ncol]-ixcen;
QY = [1:nrow]-iycen;

% average ranges offsets and half-widths
ixahw = 4;
ixaoff = -nsteps/2;
iyahw = 4;
iyaoff = 0;
%iyaoff = 8;

POSITION = [100 100 400 300];
PAPERPOSITION = [1 1 4 3];

% Plot correlation time in ML
% Scale from HW assuming exponential
tauML = damHW/log(2);

figure
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
pcolor(QX,QY,tauML);
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
title([runname_title ' Correlation Time (ML)']);

% Simple theory for tau(Q)
%tauMLth = nsteps./(2*pi*abs(QX-ixaoff));
tauMLth = nsteps./(pi*abs(QX-ixaoff)); % Unclear factor


figure
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(QX,mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o');
pa = axis;
hl = line(QX,tauMLth);
set(hl,'LineStyle','-','Color','m');
axis(pa);
% Show CTR positions
for ii = iCTR
    hl = line(nsteps*ii*[1 1],pa(3:4));
    set(hl,'LineStyle','--','Color','r');
end
xlabel('Q_X (pixels)');
ylabel('Mean Correlation Time (ML)');
title([runname_title]);

figure
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(abs(QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o');
set(gca,'Xscale','log','Yscale','log');
pa = axis;
hl = line(abs(QX-ixaoff),tauMLth);
set(hl,'LineStyle','-','Color','m');
axis(pa);
% Show CTR positions
for ii = iCTR
    hl = line(abs(nsteps*ii*[1 1]-ixaoff),pa(3:4));
    set(hl,'LineStyle','--','Color','r');
end
xlabel('Q_X (pixels)');
ylabel('Mean Correlation Time (ML)');
title([runname_title]);

figure
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(QY,mean(tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o');
xlabel('Q_Y (pixels)');
ylabel('Mean Correlation Time (ML)');
title([runname_title]);

figure
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(abs(QY),mean(tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o');
set(gca,'Xscale','log','Yscale','log');
xlabel('Q_Y (pixels)');
ylabel('Mean Correlation Time (ML)');
title([runname_title]);
