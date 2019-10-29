% Script to plot structure factor from sosvpe simulations at different
% growth rates
% From sosvpe_sf_miscut_dt by GBS
% 31-JULY-19 ICA

close all;clear all;

addpath(genpath(['/Users/ialmazn/Documents/data_analysis_ALL/XPCS_sputtering_ZnO_TiO2_2017_2019/XPCS_analysis_on_thefly']));


 %paths = ['miscut_8_14/';'miscut_8_17/';'miscut_8_21/'];
 %runnames = ['miscut_8_14';'miscut_8_17';'miscut_8_21'];
 %runname_titles = ['miscut\_8\_14';'miscut\_8\_17'; 'miscut\_8\_21'];

paths = {'miscut_32_29/';'miscut_32_30/';'miscut_32_35/';'miscut_32_37/';'miscut_32_38/';'miscut_32_39/';'miscut_32_40/';'miscut_32_36/'};%{'miscut_32_36secondhalf/';'miscut_32_31half/';}%{'miscut_32_36transient/';'miscut_32_36firsthalf/';'miscut_32_36secondhalf/'};%{'miscut_32_31firsthalf/';'miscut_32_31half/'};%{'miscut_32_31/';'miscut_32_31half/';'miscut_32_31quarter/';'miscut_32_31eight/'};%['miscut_32_31/';'miscut_32_31half/'];%['miscut_32_33/';'miscut_32_34/'];%['miscut_32_30/';'miscut_32_33/'];%['miscut_32_23/';'miscut_32_30/'];%['miscut_32_23/';'miscut_32_30/'];%['miscut_32_22/';'miscut_32_23/';'miscut_32_28/';'miscut_32_26/'];%['miscut_32_27/'];%['miscut_8_14/'];%['miscut_32_22/'];%
runnames = {'miscut_32_29';'miscut_32_30';'miscut_32_35';'miscut_32_37';'miscut_32_38';'miscut_32_39';'miscut_32_40';'miscut_32_36'};%{'miscut_32_36secondhalf';'miscut_32_31half';}% {'miscut_32_36transient';'miscut_32_36firsthalf';'miscut_32_36secondhalf'}; %{'miscut_32_29';'miscut_32_30';'miscut_32_35';'miscut_32_31'};%{'miscut_32_31firsthalf';'miscut_32_31half'};%{'miscut_32_31';'miscut_32_31half';'miscut_32_31quarter' ;'miscut_32_31eight'};%%['miscut_32_33';'miscut_32_34'];%['miscut_32_30';'miscut_32_33'];%['miscut_32_29';'miscut_32_30';'miscut_32_31'];%['miscut_32_22';'miscut_32_23';'miscut_32_28';'miscut_32_26'];%['miscut_32_27'];%['miscut_8_14'];%['miscut_32_22'];%
runname_titles = {'miscut\_32\_29';'miscut\_32\_30';'miscut\_32\_35';'miscut\_32\_37';'miscut\_32\_38';'miscut\_32\_39';'miscut\_32\_40';'miscut\_32\_36'};%{'miscut\_32\_36secondhalf';'miscut\_32\_31half';}%{'miscut\_32\_36transient';'miscut\_32\_36firsthalf';'miscut\_32\_36secondhalf'};%{'miscut\_32\_29';'miscut\_32\_30';'miscut\_32\_35';'miscut\_32\_31'};%{'miscut\_32\_31firsthalf';'miscut\_32\_31half'};%{'miscut\_32\_31';'miscut\_32\_31half';'miscut\_32\_31quarter';'miscut\_32\_31eight'};%['miscut\_32\_33';'miscut\_32\_34'];%['miscut\_32\_31';'miscut\_32\_32'];%['miscut\_32\_29';'miscut\_32\_30';'miscut\_32\_31'];%['miscut\_32\_22';'miscut\_32\_23'; 'miscut\_32\_28'; 'miscut\_32\_26'];%['miscut\_32\_22'];%['miscut\_32\_27'];%['miscut\_8\_14'];%

growth_rates = [3.3e-9 3.3e-9/2 3.3e-9/4 3.3e-9/8 3.3e-9/16 3.3e-9/32 3.3e-9/64 0.0];%[0.0 0.0 0.0];%[0.0 0.0 0.0 0.0];%[0.0 0.0];%[3.3e-9 3.3e-9]./2;%[0.0 0.0];%[ 3.0e-9 1.65e-9 1.0e-15];% in ML/s%[ 3.0e-9];%
color_rates = ['r' 'k' 'm' 'b' 'g' 'y' 'c' 'b'];%['k'];%
path_save_fig = ['./figures_presentation_v6'];
skip = 1;

addpath(genpath(['/Users/ialmazn/Box' ' Sync/XPCS_sputtering_ZnO_TiO2_2017_2019/XPCS_analysis_on_thefly']));


%make_vid(runname,[0 10],'jet');
%make_vid_sf(runname,[0 100],'jet');
if ~ skip
    for ss = 1:size(runnames,1)
        path = paths{ss,:};
        runname = [paths{ss,:} runnames{ss,:}];
        runname_title = runname_titles{ss,:};
        
        load([runname '_corr_dt_L0p5.mat']);
        
        % get time HW of Cdt for each pixel, units of ML
        
        skip_HWcalc = 1;
        if  ~skip_HWcalc
            %{
            if kk~=3
                dt_minML = 2;%if in ML
                idt = damono > dt_minML;
                ddam = dtime(idt);
            else
                dt_minML = 2.02e10;%6.0606e+08;% if in seconds which is 2/3.3e-9 to preserve the same parameters than before%2;if in ML
                idt = dtime > dt_minML;%damono > dt_minML;
                ddam = dtime(idt);%damono(idt);
            end
            %}
            ddam = ddam - ddam(1); % delta time coord (ML) (assumes evenly spaced)

            nrow = size(III_resize,1);
            ncol = size(III_resize,2);
            ndt = size(III_resize,3);
            
            damHW = Functions_sos.calc_damHW(nrow,ncol,Cdt,ndt,ddam);
            filename = char([path 'damHW_' runnames(ss,:) '.mat']);

            % Do not use first dt_minML of growth for time correlations
            %{
            
            damHW = Functions_sos.calc_damHW_byfit(nrow,ncol,Cdt,dtime(idt));
           filename = char([path 'damHW_fit_' runnames(kk,:) '.mat']);

            %}
            save(filename,'damHW');
            %save(filename,'damHW','Cdt','III','Ibar');
        else
            load([path 'damHW_' runnames{ss,:} '.mat']);
        end
        
    end
    
end

skip_do_damHWstruct = 1;



for ss = 1:size(runnames,1)
    
    path = paths{ss,:};
    runname = [paths{ss,:} runnames{ss,:}];
    runname_title = runname_titles{ss,:};
    
    if ~skip_do_damHWstruct
        load([runname '_corr_dt_L0p5.mat'],'Cdt','III_resize','Ibar');
        
        %damHW_struct(ss).Cdt = Cdt;
        damHW_struct(ss).III = III_resize;
        damHW_struct(ss).Contrast = Cdt(:,:,1);
        damHW_struct(ss).Ibar = Ibar(:,:,1);
        
        damHW_struct_tosave = damHW_struct(ss);
        
        save([path 'damHW_struct.mat'],'damHW_struct_tosave','-v7.3');
    else
        load([path 'damHW_struct.mat']);
        damHW_struct(ss) = damHW_struct_tosave;
    end
    
end


%% load ihm


for ss = 5:size(runnames,1)
    
     path = paths{ss,:};
     runname = [paths{ss,:} runnames{ss,:}];
     runname_title = runname_titles{ss,:};
 
    
     
     
%      if strcmp(runname_title ,[ 'miscut\_32\_26'])
%          load([runname '_ihm_last2000.mat']);
%          ihm = ihm_resize;
%          clear ihm_resize;
%          dtime_orig = load([runname '_stats.mat'],'dtime');
%          dtime = dtime_orig.dtime(end-size(ihm,3):end);
%          clear dtime_orig
%      elseif strcmp(runname_title ,[ 'miscut\_32\_32'])
%          load([runname '_ihm.mat']);
%          dtime_orig = load([runname '_stats.mat'],'dtime');
%          dtime = dtime_orig.dtime;
%          clear dtime_orig
%      elseif strcmp(runname_title ,[ 'miscut\_32\_33'])
%          load([runname '_ihm.mat']);
%          dtime_orig = load([runname '_stats.mat'],'dtime');
%          dtime = dtime_orig.dtime;
%          clear dtime_orig

    if strcmp(runname_title ,[ 'miscut\_32\_29'])
        load([runname '_ihm.mat']);
        load([runname '_stats.mat'],'dtime');
    elseif strcmp(runname_title ,[ 'miscut\_32\_30'])
        load([runname '_ihm.mat']);
        load([runname '_stats.mat'],'dtime');
    elseif strcmp(runname_title ,[ 'miscut\_32\_35'])
        load([runname '_ihm.mat']);
        load([runname '_stats.mat'],'dtime');
    elseif strcmp(runname_title ,[ 'miscut\_32\_36'])
        load([runname '_ihm.mat']);
        load([runname '_stats.mat'],'dtime');

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

    end

    dt_minML = dtime(end-2.50e4);%dtime(30);%dtime(end-round(2.5e4/8));%dtime(end-6.3e3);% % dtime(end-1e4);
    dt_maxML = dtime(end);%dtime(250);%dtime(end-1.25e4);%
    %idt = dtime >dt_minML;
    idt = (dtime>dt_minML & dtime<dt_maxML);
    
    idt_half =round(size(idt,1)/2);
    
    ihm_struct(ss).ihm = squeeze(ihm(:,:,idt_half));
    

    
    
end

nrow_orig = size(ihm_struct(1).ihm,1);
ncol_orig = size(ihm_struct(1).ihm,2);

%% load damHW

for ss = 1:size(runnames,1)
    
     path = paths{ss,:};
    
    
    load([path 'damHW.mat']);
    
    if ss == 1
        damHW_denom = damHW;
    end
    
    damHW_struct(ss).tauML = damHW;%./damHW_denom;
    
    
    
end

%% Initialize parameters section

number_of_runs = size(runnames,1);

ncol = size(damHW_struct(1).Ibar,2);%128;
nrow = size(damHW_struct(1).Ibar,1);%128;
nsteps = 8;

%xp = ncol*[0 0.2 0.2 0];
%yp = nrow*[0 0 0.1 0.1];
%xtx = ncol/40;
%ytx = nrow/20;
%textclr = 'k';

% central CTR is at nrow/2+1, ncol/2+1
ixcen = ncol/2 + 1;
iycen = nrow/2 + 1;
% DQ of CTRs is nsteps (pixels) in x (col)

QX = [1:ncol]-ixcen;
QY = [1:nrow]-iycen;

% average ranges offsets and half-widths
ixahw = 56;
ixaoff = -nsteps/2;%0;%5*nsteps;%3*nsteps;%nsteps;%nsteps/2;%0;%-nsteps/2;
iyahw = 15;
iyaoff = 0;%12;%4;%24;%16;%0;%
%iyaoff = 16;

POSITION = [ 98 457 2410 275];%[1 415 1233 287];%[1 12 473 690];%[51 47 387 692];
PAPERPOSITION = [1 1 4 3];
PAPERSIZE = [ 8.5000   11.0000];



% Plot correlation time in s
%tauML = tauML;
yaoff_array = [0];%[16 24];%[24];%[0 4 12 16 24];%[0];%
xaoff_array =  [-nsteps/2];%[nsteps 3*nsteps 5*nsteps];%[-nsteps/2 0 nsteps/2 ];%[0 4 12 16 24];%[-nsteps/2 ];%


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

struct_ranges_ihm = struct_ranges;
struct_ranges_ihm.QX = [1:ncol];
struct_ranges_ihm.QY = [1:nrow];
struct_ranges_ihm.row_range = [1:nrow];
struct_ranges_ihm.col_range = round(ncol_orig/2)+[-round(ncol/2):round(ncol/2)-1]; 

 % Mark CTR positions
iCTR = [-5:4];

color_multiple_offset = ['k' 'g' 'r' 'm' 'b'];


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

%% Theoretical guidelines section

% if time axis is in seconds:
tauMLth_Qx = 7.0e10*nsteps./(pi*abs(QX-ixaoff).^2); % Unclear factor
tauMLth_Qx_1 = .25e9*nsteps./(pi*abs(QX-ixaoff)); % Unclear factor
tauMLth_Qy = .5e10*nsteps./(pi*abs(QY).^2); % Unclear factor
tauMLth_Qy_3 = .1e5*nsteps./(pi*abs(QY).^3); % Unclear factor
tauMLth_Qy_2p5 = .03e5*nsteps./(pi*abs(QY).^2.5); % Unclear factor


%% 2D maps plot section

% Figure 1: time constants
figNum = 1;
clear title_figure;
for ss = 1:size(runname_titles,1)
    %title_figure(ss,:) = {[runname_titles(ss,:)],[ ' tau (sim. units)'],['in log scale']};
    title_figure(ss,:) = {[runname_titles{ss,:} ' tau (sim. units)']};
end
Function_display_sos.make_figure_2D(title_figure,damHW_struct,'tauML',struct_ranges,...
    iCTR,figNum,'CLIM',[5 10],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,'PAPERSIZE',PAPERSIZE,...
    'plot_averages',1);

%set(figure(figNum),'PaperOrientation','landscape');
%print(figure(figNum),[path_save_fig '/tau'],'-dpdf','-fillpage');

% Figure 10: diffraction pattern
figNum = 10;
% time frame to represent:
%time_frame = [1e4 1e4 5e3 1e3];
clear title_figure; clear time_frame;
for ss = 1:size(runname_titles,1)
    iii = 1;
    jjj = 1;
    time_frame(ss) = round(size(damHW_struct(ss).III,3)/2);
    title_figure(ss,:) = {[runname_titles{ss,:} ' III, @ ' num2str(time_frame(ss),'%.1e')]};%[runname_titles(ss,:) ' Diffracted intensity, time frame = ' num2str(time_frame)];
    map2D(ss).III = damHW_struct(ss).III(:,:,time_frame(ss));
end


Function_display_sos.make_figure_2D( title_figure,map2D,'III',struct_ranges,iCTR,figNum,...
    'CLIM',[-10 -2],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);

%print('-r600',figure(figNum),'III','-dpng');


% Figure 11: IBar
figNum = 11;
% time frame to represent:
time_frame = [1e4 1e4 5e3];
clear title_figure;
for ss = 1:size(runname_titles,1)
    title_figure(ss,:) = {[runname_titles{ss,:} ' Ibar']};
    map2D(ss).Ibar = abs(damHW_struct(ss).Ibar);
end


Function_display_sos.make_figure_2D( title_figure,map2D,'Ibar',struct_ranges,iCTR,figNum,...
    'CLIM',[-10 -2],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);

% Figure 12: Contrast
figNum = 12;
% time frame to represent:
%time_frame = 130;
clear title_figure;
for ss = 1:size(runname_titles,1)
    title_figure(ss,:) = {[runname_titles{ss,:} ' Contrast']};
    map2D(ss).Contrast = 10.^(damHW_struct(ss).Contrast);
end


Function_display_sos.make_figure_2D( title_figure,map2D,'Contrast',struct_ranges,iCTR,figNum,...
    'CLIM',[0 2],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);




% Figure 13: ihm
figNum = 13;
% time frame to represent:
%time_frame = 130;
clear title_figure;
for ss = 1:size(runname_titles,1)
    title_figure(ss,:) = {[runname_titles{ss,:} ' ihm']};
    ihm_struct(ss).ihm_square = ihm_struct(ss).ihm(struct_ranges_ihm.row_range,struct_ranges_ihm.col_range,:);
end


Function_display_sos.make_figure_2D_realspace( title_figure,ihm_struct,'ihm_square',struct_ranges_ihm,figNum,...
    'CLIM',[0 35],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);




%% 1 D averages

% Figures 2 and 4 (averaged tau_ML over colums or rows):

% average over rows
range_rows = iycen+iyaoff+[-iyahw:iyahw];

for ss = 1:size(runnames,1)
    tauML_1D = mean(damHW_struct(ss).tauML(range_rows,:),1);
    damHW1D_struct(ss).tauML_1D_Qx = tauML_1D;

   %{ 
   if ss== 1
        tauML_1D_Qx_denom = tauML_1D;
        damHW1D_struct(ss).tauML_1D_Qx = tauML_1D./damHW1D_struct(1).tauML_1D_Qx;
    else
        damHW1D_struct(ss).tauML_1D_Qx = tauML_1D./tauML_1D_Qx_denom;
    end
    %}
   
    Qx_offset = QX-ixaoff;
    Qx_positive = Qx_offset(Qx_offset>0);
    tauML_Qx_positive = damHW1D_struct(ss).tauML_1D_Qx(Qx_offset>0);%tauML_1D(Qx_offset>0);
    tauML_th_QX_positive = tauMLth_Qx(Qx_offset>0);%tauMLth_Qx_1(Qx_offset>0);%
    
    damHW1D_struct(ss).tauML_1D_Qx_positive = tauML_Qx_positive;
    damHW1D_struct(ss).tauML_th_QX_positive = tauML_th_QX_positive;
    
    Qx_negative = abs(Qx_offset(Qx_offset<=0));
    tauML_Qx_negative = damHW1D_struct(ss).tauML_1D_Qx(Qx_offset<=0);%tauML_1D(Qx_offset<=0);
    tauML_th_QX_negative =  tauMLth_Qx(Qx_offset<=0);%tauMLth_Qx_1(Qx_offset<=0);%
    
    damHW1D_struct(ss).tauML_1D_Qx_negative = tauML_Qx_negative;
    damHW1D_struct(ss).tauML_th_QX_negative = tauML_th_QX_negative;
    
    damHW1D_struct(ss).Color = color_rates(ss);
    legend_str{ss} = [runname_titles{ss,:} ' growth rate (ML/s) = ' num2str(growth_rates(ss))];
end


figNum = 2;

Function_display_sos.make_figures_tauML1D(['offset in Q_y = ' num2str(iyaoff)],...
    QX,damHW1D_struct,...
    'tauML_1D_Qx',figNum,...
    'iCTR',iCTR,'nsteps',nsteps,...
    'XLABEL','Q_x (pixels)','YLIM',[1e5 1e10]);%[1 10]);%

legend(legend_str);


figNum = 20;

Function_display_sos.make_figures_tauML1D(['QX is positive,  offset in Qy = ' num2str(iyaoff)],Qx_positive,damHW1D_struct,'tauML_1D_Qx_positive',...
    figNum,'tauMLth',tauML_th_QX_positive,'iCTR',iCTR,'nsteps',nsteps,'ixaoff',ixaoff,'XLABEL','Q_x (pixels)',...
     'POSITION',[534    64   769   538],'FONTSIZE',15,'XSCALE','log','YSCALE',...
    'log','YLIM',[1e5 1e10]);%[1 100]);%,

legend(legend_str);

figNum = 21;

Function_display_sos.make_figures_tauML1D(['QX is negative,  offset in Qy = ' num2str(iyaoff)],Qx_negative,damHW1D_struct,'tauML_1D_Qx_negative',...
    figNum,'tauMLth',tauML_th_QX_negative,'iCTR',iCTR,'nsteps',nsteps,'ixaoff',ixaoff,'XLABEL','Q_x (pixels)',...
     'POSITION',[534    64   769   538],'FONTSIZE',15,'XSCALE','log','YSCALE',...
    'log','YLIM',[1e5 1e10]);%[1e6 1e10]);%[1 100]);%,

legend(legend_str);


% average over columns
range_cols = ixcen+ixaoff+[-ixahw:ixahw];

for ss = 1:size(runnames,1)
    tauML_1D = mean(damHW_struct(ss).tauML(:,range_cols),2);
    damHW1D_struct(ss).tauML_1D_Qy = tauML_1D;
    
    %{ 
    if ss== 1
        tauML_1D_Qy_denom = tauML_1D;
  
     end
     damHW1D_struct(ss).tauML_1D_Qy = tauML_1D./tauML_1D_Qy_denom;
    %}
    
    Qy_offset = QY-iyaoff;
    Qy_positive = Qy_offset(Qy_offset>0);
    tauML_Qy_positive = tauML_1D(Qy_offset>0);
    tauML_th_Qy_positive = tauMLth_Qy(Qy_offset>0);
    
    damHW1D_struct(ss).tauML_1D_Qy_positive = tauML_Qy_positive;
    damHW1D_struct(ss).tauML_th_QY_positive = tauML_th_Qy_positive;
    
    Qy_negative = abs(Qy_offset(Qy_offset<=0));
    tauML_Qy_negative = tauML_1D(Qy_offset<=0);
    tauML_th_Qy_negative = tauMLth_Qy(Qy_offset<=0);
    
    damHW1D_struct(ss).tauML_1D_Qy_negative = tauML_Qy_negative;
    damHW1D_struct(ss).tauML_th_QY_negative = tauML_th_Qy_negative;
    
    damHW1D_struct(ss).Color = color_rates(ss);
    legend_str{ss} = [runname_titles{ss,:} ' growth rate (ML/s) = ' num2str(growth_rates(ss))];
end


figNum = 4;

Function_display_sos.make_figures_tauML1D(['offset in Q_x = ' num2str(ixaoff)],...
    QY,damHW1D_struct,'tauML_1D_Qy',figNum,...
    'XLABEL','Q_y (pixels)','YLIM',[1e5 1e10]);%,[1 10]);%

legend(legend_str);


figNum = 40;

Function_display_sos.make_figures_tauML1D(['QY is positive,  offset in Qx = ' num2str(ixaoff)],Qy_positive,damHW1D_struct,'tauML_1D_Qy_positive',...
    figNum,'tauMLth',tauML_th_Qy_positive,'iyaoff',iyaoff,'XLABEL','Q_y (pixels)',...
    'POSITION',[534    64   769   538],'FONTSIZE',15,'XSCALE','log','YSCALE',...
    'log','YLIM',[1e5 1e10]);%,[1 10]);%

legend(legend_str,'1/Q^2');

figNum = 41;

Function_display_sos.make_figures_tauML1D(['QY is negative,  offset in Qx = ' num2str(ixaoff)],Qy_negative,damHW1D_struct,'tauML_1D_Qy_negative',...
    figNum,'tauMLth',tauML_th_Qy_negative,'iyaoff',iyaoff,'XLABEL','Q_y (pixels)',...
    'POSITION',[534    64   769   538],'FONTSIZE',15,'XSCALE','log','YSCALE',...
    'log','YLIM',[1e5 1e10]);%,[1 10]);%
legend(legend_str,'1/Q^2');


% for one single Q vs growth rate

QY_val = 30;
figure(5);
clf;
hold on;
for ss = 1:size(runnames,1)
    tauML_1D = mean(damHW_struct(ss).tauML(:,range_cols),2);
    damHW1D_struct(ss).tauML_1D_oneQ = tauML_1D(QY_val);
    
    plot(growth_rates(ss),damHW1D_struct(ss).tauML_1D_oneQ,'Marker','o','MarkerSize',15,'Color',color_rates(ss));
    
end
xlabel('Growth rates (ML/sim units)');
ylabel('\tau (sim units)');
set(gca,'FontSize',15);
box on;
title(['Time constant vs growth rate at QY = ' num2str(QY_val)]);
set(gcf,'Position',[534    64   769   538]);
set(gca,'YScale','log');
ylim([1e5 1e10]);

QX_val = 30;
figure(6);
clf;
hold on;
for ss = 1:size(runnames,1)
    tauML_1D = mean(damHW_struct(ss).tauML(range_rows,:),1);
    damHW1D_struct(ss).tauML_1D_oneQ = tauML_1D(QX_val);
    
    plot(growth_rates(ss),damHW1D_struct(ss).tauML_1D_oneQ,'Marker','o','MarkerSize',15,'Color',color_rates(ss));
    
end
xlabel('Growth rates (ML/sim units)');
ylabel('\tau (sim units)');
set(gca,'FontSize',15);
box on;
title(['Time constant vs growth rate at QX = ' num2str(QY_val)]);
set(gcf,'Position',[534    64   769   538]);
set(gca,'YScale','log');
ylim([1e5 1e10]);





return;


figure;
plot(dtime,damono,'LineWidth',3.0);hold on;
hl = line(dtime(end-2.5e4).*ones(numel(dtime),1),[1:numel(dtime)].*2/numel(dtime)); 
set(hl,'Color','r','LineWidth',3.0);
hl = line(dtime(end-1e4).*ones(numel(dtime),1),[1:numel(dtime)].*2/numel(dtime)); 
set(hl,'Color','k','LineWidth',3.0);
hl = line(dtime(end-round(2.5e4/4)).*ones(numel(dtime),1),[1:numel(dtime)].*2/numel(dtime)); 
set(hl,'Color','m','LineWidth',3.0);
hl = line(dtime(end-round(2.5e4/8)).*ones(numel(dtime),1),[1:numel(dtime)].*2/numel(dtime)); 
set(hl,'Color','b','LineWidth',3.0);
xlabel('dtime (unit sim.)');ylabel('damono (ML)');
set(gca,'FontSize',15);    


%{
figure(1);
clf;
set(gcf,'Position',[90         245        1031         387]);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
subplot(1,number_of_runs,1);
pcolor(QX,QY,log10(damHW_struct(1).tauML));
set(gca,'Clim',[6.1971    9.3297]);
shading flat;
axis image;
% Show ranges averaged
Function_display_sos.show_ranges_averaged(QX,QY,yaoff_array,iyahw,xaoff_array,ixahw)
% Mark CTR positions
iCTR = [-5:4];
hl = line(nsteps*iCTR,zeros(size(iCTR)));
set(hl,'Linestyle','none','Marker','x','Color','r');
xlabel('Q_X (pixels)');
ylabel('Q_Y (pixels)');
colorbar;
%title([runname_title ' Correlation Time (ML) in log scale']);
title([runname_titles(1,:) ' Tau (sim. unit time) in log scale']);
set(gca,'FontSize',10)

subplot(1,number_of_runs,2);
pcolor(QX,QY,log10(damHW_struct(2).tauML));
set(gca,'Clim',[6.1971    9.3297]);
shading flat;
axis image;
% plot averages
Function_display_sos.show_ranges_averaged(QX,QY,yaoff_array,iyahw,xaoff_array,ixahw)
% Mark CTR positions
iCTR = [-5:4];
hl = line(nsteps*iCTR,zeros(size(iCTR)));
set(hl,'Linestyle','none','Marker','x','Color','r');
xlabel('Q_X (pixels)');
ylabel('Q_Y (pixels)');
colorbar;
title([runname_titles(2,:) ' Tau (sim. unit time) in log scale']);
set(gca,'FontSize',10);

subplot(1,number_of_runs,3);
pcolor(QX,QY,log10(damHW_struct(3).tauML));
set(gca,'Clim',[6.1971    9.3297]);
shading flat;
axis image;
% plot averages
Function_display_sos.show_ranges_averaged(QX,QY,yaoff_array,iyahw,xaoff_array,ixahw)
% Mark CTR positions
iCTR = [-5:4];
hl = line(nsteps*iCTR,zeros(size(iCTR)));
set(hl,'Linestyle','none','Marker','x','Color','r');
xlabel('Q_X (pixels)');
ylabel('Q_Y (pixels)');
colorbar;
title([runname_titles(2,:) ' Tau (sim. unit time) in log scale']);
set(gca,'FontSize',10)

%print('-r600',figure(1),[path_save_fig '/tau2Maps_allgrowthrates'],'-dpng');

% Simple theory for tau(Q)
%tauMLth = nsteps./(2*pi*abs(QX-ixaoff));
tauMLth_Qx = .5e10*nsteps./(pi*(QX-ixaoff).^2); % Unclear factor
tauMLth_Qy = .08e4*nsteps./(pi*(QY).^2); % Unclear factor
tauMLth_Qy_3 = .1e5*nsteps./(pi*(QY).^3); % Unclear factor
tauMLth_Qy_2p5 = .03e5*nsteps./(pi*(QY).^2.5); % Unclear factor



figure(2);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(QX,mean(damHW_struct(1).tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o','Color','r');
pa = axis;
%hl = line(QX,tauMLth_Qx);
%set(hl,'LineStyle','-','Color','m','LineWidth',2.0);
axis(pa);
hl = line(QX,mean(damHW_struct(2).tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o','Color','k');
axis(pa);
hl = line(QX,mean(damHW_struct(3).tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o','Color','b');
axis(pa);
% Show CTR positions
for ii = iCTR
    hl = line(nsteps*ii*[1 1],pa(3:4));
    set(hl,'LineStyle','--','Color','r');
end
xlabel('Q_X (pixels)');
%ylabel('Mean Correlation Time (ML)');
ylabel('Mean Correlation Time (frames)');
ylim([1e-3 3e9]);
title([runname_title]);
legend(runname_titles(1,:),runname_titles(2,:),runname_titles(3,:));





% plot log(tau) vs log(Q_x-Q_1)

figure(3);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
%hl = line(abs(QX-ixaoff),mean(tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
hl = line((QX-ixaoff),mean(damHW_struct(1).tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o','Color','r');
set(gca,'Xscale','log','Yscale','log');
set(gca,'Ylim',[1e6 5e10]);
pa = axis;
hl = line((QX-ixaoff),mean(damHW_struct(2).tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o','Color','k');
axis(pa);
hl = line((QX-ixaoff),mean(damHW_struct(3).tauML(iycen+iyaoff+[-iyahw:iyahw],:),1));
set(hl,'LineStyle','none','Marker','o','Color','b');
axis(pa);
hl = line(abs(QX-ixaoff),tauMLth_Qx);
set(hl,'LineStyle','-','Color','m','LineWidth',2.0);
axis(pa);
% Show CTR positions
for ii = iCTR
    hl = line(abs(nsteps*ii*[1 1]-ixaoff),pa(3:4));
    set(hl,'LineStyle','--','Color','r');
end
xlabel('|Q_X - Q_1| (pixels)');
ylabel('Mean Correlation Time (frames)');
%ylabel('Mean Correlation Time (s)');
%title([runname_title]);
legend(runname_titles(1,:),runname_titles(2,:),runname_titles(3,:),'fit to 1/Q^2');
title(['offset in y =' num2str(iyaoff)]);


figure(4);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(QY,mean(damHW_struct(1).tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o','Color','r');
hold on;
%hl2 = line(QY,tauMLth_Qy);
%set(hl2,'LineStyle','-','Color','m');
hl = line(QY,mean(damHW_struct(2).tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o','Color','k');
hl = line(QY,mean(damHW_struct(3).tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o','Color','b');
xlabel('Q_Y (pixels)');
%ylabel('Mean Correlation Time (ML)');
ylabel('Mean Correlation Time (frames)');
%title([runname_title]);
legend(runname_titles(1,:),runname_titles(2,:),runname_titles(3,:));



figure(5);
clf;
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);
axes('Box','on');
hl = line(abs(QY),mean(damHW_struct(1).tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o','Color','r');
hold on;
%hl2 = line(abs(QY),tauMLth_Qy);
%set(hl2,'LineStyle','-','Color','m','LineWidth',2.0);
hl = line(abs(QY),mean(damHW_struct(2).tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o','Color','k');
hl = line(abs(QY),mean(damHW_struct(3).tauML(:,ixcen+ixaoff+[-ixahw:ixahw]),2));
set(hl,'LineStyle','none','Marker','o','Color','b');
set(gca,'Xscale','log','Yscale','log');
set(gca,'Ylim',[1e6 1e10]);
xlabel('Q_Y (pixels)');
%ylabel('Mean Correlation Time (ML)');
ylabel('Mean Correlation Time (frames)');
%title([runname_title]);
legend(runname_titles(1,:),runname_titles(2,:),runname_titles(3,:));
title(['offset in Q_x = ' num2str(ixaoff)]);

% plot the depencece of the tau at the first CTR for all the growth rates

% arrhenius law?
tau_CTR = 5e-9;
damHW_CTR_arrhenius = 1e11*exp(-tau_CTR./growth_rates);

figure(6);
clf;
hold on;
for kk = 1:size(runnames,1)
   damHW_vsQx = mean(damHW_struct(kk).tauML(iycen+iyaoff+[-iyahw:iyahw],:),1);
   
   damHW_CTR(kk) = damHW_vsQx(find(QX == 0)); 
   x_arrhenius(kk) = growth_rates(kk);
   plot(x_arrhenius(kk),damHW_CTR(kk),'o','Color',color_rates(kk),'MarkerSize',10);
   
end
%set(gca,'YScale','log');
%set(gca,'XScale','log');
box 'on';
xlabel('growth rates (ML/s)');
ylabel('Mean correlation time (frames)');
set(gca,'FontSize',15);
set(gcf,'Position',POSITION);
set(gcf,'PaperPosition',PAPERPOSITION);

%figure(6);
%hold on;
%plot( x_arrhenius,damHW_CTR_arrhenius,'r','LineWidth',3.0);


%%% plot the intensity correlation function for a single pixel:
growth_rate_index = 3;

pixel_row = 40;
pixel_col = 20;

Cdti = squeeze(damHW_struct(3).Cdt(pixel_row,pixel_col,:));
Cdti = Cdti/Cdti(1);

figure(1000);
plot(dtime(idt),Cdti,'o');
xlabel('time (sim. unit time)');
ylabel('Cdti');
title([runname_titles(growth_rate_index,:) ' at pixels = ' num2str(pixel_row) ',' num2str(pixel_col)]);
%}