function make_vid_sf(runname,clim,cmap)
% function make_vid(runname)
% to make video file of images created by sosvpe
% 15-Aug-14 GBS

load([runname '_corr_dt.mat']);

clear MOV1
% Prepare the new file.
vidObj = VideoWriter([runname 'sf_vid.avi']);
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
imagesc(squeeze(log(abs(ftm(:,:,1)).^2)),clim);
axis image
shading flat;
colormap(cmap);
set(gca,'Visible','off');
xp = ncol*[0 0.6 0.6 0];
yp = nrow*[0 0 0.1 0.1];
patch(xp,yp,'w');
xtx = ncol/40;
ytx = nrow/20;
ht = text(xtx,ytx,[num2str(damono(1),'%5.2f'),' ML, time frame = ',int2str(1)]);
set(ht,'Color',textclr,'FontSize',12);
pause(0.1);

hc = get(gca,'Children');

MOV1(1) = getframe;
    % Write frame to the video file.
writeVideo(vidObj,MOV1(1));
    
for ipt = 2:size(ftm,3)
    set(hc(end),'CData',squeeze(log(abs(ftm(:,:,ipt)).^2)));
    patch(xp,yp,'w');
    ht = text(xtx,ytx,[num2str(damono(ipt),'%5.2f'),' ML, timeframe = ',int2str(ipt)]);
    set(ht,'Color',textclr,'FontSize',12);
    pause(0.1);
    
    MOV1(ipt) = getframe;
        % Write frame to the video file.
    writeVideo(vidObj,MOV1(ipt));        
end

close(vidObj);
return;
end
