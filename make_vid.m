function make_vid(runname,clim,cmap)
% function make_vid(runname)
% to make video file of images created by sosvpe
% 15-Aug-14 GBS

load([runname '_stats.mat']);
load([runname '_ihm.mat']);
%load([runname '_corr_dt.mat']);

clear MOV1
% Prepare the new file.
vidObj = VideoWriter([runname '_vid.avi']);
open(vidObj);

textclr = 'k'; % Parameters for plotted images
if nargin < 2
    clim = [1 6.5]; 
end
if nargin < 3
    cmap = 'jet'; 
end

figure;
set(gcf,'Position',[600 200 400 400]);
set(gcf,'PaperPosition',[1 1 4 4]);
axes('Box','on');
imagesc(squeeze(ihm(:,:,1)),clim);
axis image
shading flat;
colormap(cmap);
set(gca,'Visible','off');
xp = ncol*[0 0.6 0.6 0];
yp = nrow*[0 0 0.1 0.1];
patch(xp,yp,'w');
xtx = ncol/40;
ytx = nrow/20;
%ht = text(xtx,ytx,[num2str(damono(1),'%5.2f'),' ML, time = ',int2str(dtime(1))]);
ht = text(xtx,ytx,[num2str(damono(1),'%5.2f'),' ML, time frame = 1']);
set(ht,'Color',textclr,'FontSize',12);
pause(0.1);

hc = get(gca,'Children');

MOV1(1) = getframe;
    % Write frame to the video file.
writeVideo(vidObj,MOV1(1));
    
for ipt = 2:size(ihm,3)
    set(hc(end),'CData',squeeze(ihm(:,:,ipt)));
    patch(xp,yp,'w');
    ht = text(xtx,ytx,[num2str(damono(ipt),'%5.2f'),' ML, time frame= ',int2str(ipt)]);
    %ht = text(xtx,ytx,[num2str(damono(ipt),'%5.2f'),' ML, time = ',int2str(dtime(ipt))]);
    set(ht,'Color',textclr,'FontSize',12);
    pause(0.1);
    
    MOV1(ipt) = getframe;
        % Write frame to the video file.
    writeVideo(vidObj,MOV1(ipt));        
end

close(vidObj);
return;
end
