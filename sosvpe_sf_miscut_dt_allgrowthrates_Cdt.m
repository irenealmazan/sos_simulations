close all;clear all;

computer_script; % set the value of computer_flag

if strcmp(computer_flag,'blueshift')
    addpath(genpath(['./XPCS_analysis_on_thefly']));
    
else
    addpath(genpath(['/Users/ialmazn/Documents/data_analysis_ALL/XPCS_sputtering_ZnO_TiO2_2017_2019/XPCS_analysis_on_thefly']));
end

 %paths = ['miscut_8_14/';'miscut_8_17/';'miscut_8_21/'];
 %runnames = ['miscut_8_14';'miscut_8_17';'miscut_8_21'];
 %runname_titles = ['miscut\_8\_14';'miscut\_8\_17'; 'miscut\_8\_21'];

paths = {'miscut_32_36/';'miscut_32_40/';'miscut_32_39/';'miscut_32_38/';'miscut_32_37/';'miscut_32_35/';'miscut_32_30/';'miscut_32_29/'};%
runnames = {'miscut_32_36';'miscut_32_40';'miscut_32_39';'miscut_32_38';'miscut_32_37';'miscut_32_35';'miscut_32_30';'miscut_32_29'};%{'miscut_32_36secondhalf';'miscut_32_31half';}% {'miscut_32_36transient';'miscut_32_36firsthalf';'miscut_32_36secondhalf'}; %{'miscut_32_29';'miscut_32_30';'miscut_32_35';'miscut_32_31'};%{'miscut_32_31firsthalf';'miscut_32_31half'};%{'miscut_32_31';'miscut_32_31half';'miscut_32_31quarter' ;'miscut_32_31eight'};%%['miscut_32_33';'miscut_32_34'];%['miscut_32_30';'miscut_32_33'];%['miscut_32_29';'miscut_32_30';'miscut_32_31'];%['miscut_32_22';'miscut_32_23';'miscut_32_28';'miscut_32_26'];%['miscut_32_27'];%['miscut_8_14'];%['miscut_32_22'];%
runname_titles = {'miscut\_32\_36';'miscut\_32\_40';'miscut\_32\_39';'miscut\_32\_38';'miscut\_32\_37';'miscut\_32\_35';'miscut\_32\_30';'miscut\_32\_29'};%{'miscut\_32\_36secondhalf';'miscut\_32\_31half';}%{'miscut\_32\_36transient';'miscut\_32\_36firsthalf';'miscut\_32\_36secondhalf'};%{'miscut\_32\_29';'miscut\_32\_30';'miscut\_32\_35';'miscut\_32\_31'};%{'miscut\_32\_31firsthalf';'miscut\_32\_31half'};%{'miscut\_32\_31';'miscut\_32\_31half';'miscut\_32\_31quarter';'miscut\_32\_31eight'};%['miscut\_32\_33';'miscut\_32\_34'];%['miscut\_32\_31';'miscut\_32\_32'];%['miscut\_32\_29';'miscut\_32\_30';'miscut\_32\_31'];%['miscut\_32\_22';'miscut\_32\_23'; 'miscut\_32\_28'; 'miscut\_32\_26'];%['miscut\_32\_22'];%['miscut\_32\_27'];%['miscut\_8\_14'];%


growth_rates_all = [0.0 3.3e-9/64 3.3e-9/32 3.3e-9/16 3.3e-9/8 3.3e-9/4 3.3e-9/2 3.3e-9];
color_rates_allgrowth = ['or';'ok';'om';'ob';'og';'oy';'oc';'xb'];%['k'];%
path_save_fig = ['./figures_presentation_v6'];



for ss = 1:size(runnames,1)
    
    path = paths{ss,:};
    runname = [paths{ss,:} runnames{ss,:}];
    runname_title = runname_titles{ss,:};
    
    
    load([runname '_corr_dt_L0p5.mat'],'Cdt','ddam');
    
    Cdt_struct(ss).Cdt = Cdt;
    Cdt_struct(ss).ddam = ddam;
    
    %load([path 'damHW.mat']);
    %Cdt_struct(ss).damHW = damHW;
    
end

nrow = size(Cdt,1);
ncol = size(Cdt,2);
nsteps = 8; 

%% fit the CdtI with two single exponential decays (taking as the fast decay time constant, 
%  the time constant obtained from the fit of the 0 growth rate)


if ~skip_fit
    
    %{
    ss = 1;
    path = paths{ss,:};
    damHW_file = [path 'damHW_fit_singleexp_fixcontr.mat'];%%[path 'damHW_fit_singleexp.mat'];%[path 'damHW_fit.mat'];%[path 'damHW.mat'];%
    load(damHW_file);
    damHW0 = fit_res.damHW;
    
    nrow = size(damHW0,1);
    ncol = size(damHW0,2);
    nsteps = 8;
    
    
    Cdt_struct(ss).damHW_fast = fit_res.damHW;
    Cdt_struct(ss).damHW_slow = zeros(nrow,ncol);
    Cdt_struct(ss).back = fit_res.back;
    Cdt_struct(ss).contrast_fast = fit_res.contrast;
    Cdt_struct(ss).contrast_slow = zeros(nrow,ncol);
    Cdt_struct(ss).slope =  zeros(nrow,ncol);
    
    
    Cdt_struct(ss).damHW_slow_sigma = zeros(nrow,ncol);
    Cdt_struct(ss).damHW_fast_sigma = fit_res.sigma_damHW;
    Cdt_struct(ss).back_sigma = fit_res.sigma_back;
    Cdt_struct(ss).contrast_slow_sigma = zeros(nrow,ncol);
    Cdt_struct(ss).contrast_fast_sigma = fit_res.sigma_contrast;
    Cdt_struct(ss).slope_sigma = fit_res.sigma_slope;
    
    %}
    
    
    
    for ss = 1:size(runnames,1)
        path = paths{ss,:};
        
        display(path)
        
        damHW_file = [path 'damHW_fit_singleexp_fixcontr.mat'];%[path 'damHW_fit_doubleexp.mat'];%
        
        
        [fit_res] = Functions_sos.calc_damHW_byfit(nrow,ncol,Cdt_struct(ss).Cdt,Cdt_struct(ss).ddam,path);
        %Functions_sos.calc_damHW_byfit_double_exp_contrast(nrow,ncol,Cdt_struct(ss).Cdt,damHW0,Cdt_struct(ss).ddam);
        %
        
        
        save([damHW_file],'fit_res');
        
        switch damHW_file
            case [path 'damHW_fit_doubleexp.mat']
                Cdt_struct(ss).damHW_slow = fit_res.damHW_slow;
                Cdt_struct(ss).damHW_fast = fit_res.damHW_fast;
                Cdt_struct(ss).back = fit_res.back;
                Cdt_struct(ss).contrast_slow = fit_res.contrast_slow;
                Cdt_struct(ss).contrast_fast = fit_res.contrast_fast;
                Cdt_struct(ss).slope = fit_res.slope;
                
                Cdt_struct(ss).damHW_slow_sigma = fit_res.sigma_damHW_slow;
                Cdt_struct(ss).damHW_fast_sigma = fit_res.sigma_damHW_fast;
                Cdt_struct(ss).back_sigma = fit_res.sigma_back;
                Cdt_struct(ss).contrast_slow_sigma = fit_res.sigma_contrast_slow;
                Cdt_struct(ss).contrast_fast_sigma = fit_res.sigma_contrast_fast;
                Cdt_struct(ss).slope_sigma = fit_res.sigma_slope;
                
            case [path 'damHW_fit_singleexp.mat']
                Cdt_struct(ss).damHW_fast = zeros(nrow,ncol);
                Cdt_struct(ss).damHW_slow = fit_res.damHW;
                Cdt_struct(ss).back = fit_res.back;
                Cdt_struct(ss).contrast_fast = zeros(nrow,ncol);
                Cdt_struct(ss).contrast_slow = fit_res.contrast;
                Cdt_struct(ss).slope =  zeros(nrow,ncol);
                
                Cdt_struct(ss).damHW_slow_sigma = fit_res.sigma_damHW;
                Cdt_struct(ss).damHW_fast_sigma = zeros(nrow,ncol);
                Cdt_struct(ss).back_sigma = fit_res.sigma_back;
                Cdt_struct(ss).contrast_slow_sigma =  fit_res.sigma_contrast;
                Cdt_struct(ss).contrast_fast_sigma =zeros(nrow,ncol);
                Cdt_struct(ss).slope_sigma = fit_res.sigma_slope;
                
           case [path 'damHW_fit_singleexp_fixcontr.mat']
                Cdt_struct(ss).damHW_fast = zeros(nrow,ncol);
                Cdt_struct(ss).damHW_slow = fit_res.damHW;
                Cdt_struct(ss).back = fit_res.back;
                Cdt_struct(ss).contrast_fast = zeros(nrow,ncol);
                Cdt_struct(ss).contrast_slow = fit_res.contrast;
                Cdt_struct(ss).slope =  zeros(nrow,ncol);
                
                Cdt_struct(ss).damHW_slow_sigma = fit_res.sigma_damHW;
                Cdt_struct(ss).damHW_fast_sigma = zeros(nrow,ncol);
                Cdt_struct(ss).back_sigma = fit_res.sigma_back;
                Cdt_struct(ss).contrast_slow_sigma =  fit_res.sigma_contrast;
                Cdt_struct(ss).contrast_fast_sigma =zeros(nrow,ncol);
                Cdt_struct(ss).slope_sigma = fit_res.sigma_slope;
                
        end    
            
    end
    
else
    
    ss = 1;
    
    path = paths{ss,:};
    
    damHW_file = [path 'damHW_fit_singleexp.mat'];%[path 'damHW_fit.mat'];%[path 'damHW.mat'];%
    load(damHW_file);
    
    nrow = size(fit_res.damHW,1);
    ncol = size(fit_res.damHW,2);
    nsteps = 8;
    
    
    Cdt_struct(ss).damHW_fast = fit_res.damHW;
    Cdt_struct(ss).damHW_slow = zeros(nrow,ncol);
    Cdt_struct(ss).back = fit_res.back;
    Cdt_struct(ss).contrast_fast = fit_res.contrast;
    Cdt_struct(ss).contrast_slow = zeros(nrow,ncol);
    Cdt_struct(ss).slope =  zeros(nrow,ncol);
    
    
    Cdt_struct(ss).damHW_slow_sigma = zeros(nrow,ncol);
    Cdt_struct(ss).damHW_fast_sigma = fit_res.sigma_damHW;
    Cdt_struct(ss).back_sigma = fit_res.sigma_back;
    Cdt_struct(ss).contrast_slow_sigma = zeros(nrow,ncol);
    Cdt_struct(ss).contrast_fast_sigma = fit_res.sigma_contrast;
    Cdt_struct(ss).slope_sigma = fit_res.sigma_slope;
    
    %%{
    for ss = 2:size(runnames,1)
        path = paths{ss,:};
        
        damHW_file = [path 'damHW_fit_singleexp.mat'];%%[path 'damHW_fit_doubleexp.mat'];%[path 'damHW.mat'];%[path 'damHW_fit_doubleexp.mat'];%
        
        load([damHW_file]);
        
        
        switch damHW_file
            case [path 'damHW_fit_doubleexp.mat']
                Cdt_struct(ss).damHW_slow = fit_res.damHW_slow;
                Cdt_struct(ss).damHW_fast = fit_res.damHW_fast;
                Cdt_struct(ss).back = fit_res.back;
                Cdt_struct(ss).contrast_slow = fit_res.contrast_slow;
                Cdt_struct(ss).contrast_fast = fit_res.contrast_fast;
                Cdt_struct(ss).slope = fit_res.slope;
                
                Cdt_struct(ss).damHW_slow_sigma = fit_res.sigma_damHW_slow;
                Cdt_struct(ss).damHW_fast_sigma = fit_res.sigma_damHW_fast;
                Cdt_struct(ss).back_sigma = fit_res.sigma_back;
                Cdt_struct(ss).contrast_slow_sigma = fit_res.sigma_contrast_slow;
                Cdt_struct(ss).contrast_fast_sigma = fit_res.sigma_contrast_fast;
                Cdt_struct(ss).slope_sigma = fit_res.sigma_slope;
                
            case [path 'damHW_fit_singleexp.mat']
                Cdt_struct(ss).damHW_fast = fit_res.damHW;
                Cdt_struct(ss).damHW_slow = zeros(nrow,ncol);
                Cdt_struct(ss).back = fit_res.back;
                Cdt_struct(ss).contrast_fast = fit_res.contrast;
                Cdt_struct(ss).contrast_slow = zeros(nrow,ncol);
                Cdt_struct(ss).slope =  zeros(nrow,ncol);
                
                
                Cdt_struct(ss).damHW_slow_sigma = zeros(nrow,ncol);
                Cdt_struct(ss).damHW_fast_sigma = fit_res.sigma_damHW;
                Cdt_struct(ss).back_sigma = fit_res.sigma_back;
                Cdt_struct(ss).contrast_slow_sigma = zeros(nrow,ncol);
                Cdt_struct(ss).contrast_fast_sigma = fit_res.sigma_contrast;
                Cdt_struct(ss).slope_sigma = fit_res.sigma_slope;
                
            case [path 'damHW_fit_singleexp_fixcontr.mat']
                Cdt_struct(ss).damHW_fast = fit_res.damHW;
                Cdt_struct(ss).damHW_slow = zeros(nrow,ncol);
                Cdt_struct(ss).back = fit_res.back;
                Cdt_struct(ss).contrast_fast = fit_res.contrast;
                Cdt_struct(ss).contrast_slow = zeros(nrow,ncol);
                Cdt_struct(ss).slope =  zeros(nrow,ncol);
                
                
                Cdt_struct(ss).damHW_slow_sigma = zeros(nrow,ncol);
                Cdt_struct(ss).damHW_fast_sigma = fit_res.sigma_damHW;
                Cdt_struct(ss).back_sigma = fit_res.sigma_back;
                Cdt_struct(ss).contrast_slow_sigma = zeros(nrow,ncol);
                Cdt_struct(ss).contrast_fast_sigma = fit_res.sigma_contrast;
                Cdt_struct(ss).slope_sigma = fit_res.sigma_slope;
        end
    end
    
    %}
end


%%

return;

%% Initialize parameters to plot

% central CTR is at nrow/2+1, ncol/2+1
ixcen = ncol/2 + 1;
iycen = nrow/2 + 1;
% DQ of CTRs is nsteps (pixels) in x (col)

QX = [1:ncol]-ixcen;
QY = [1:nrow]-iycen;

% QX original:
%QX_orig = [1:ncol_orig]-ncol_orig/2;

% average ranges offsets and half-widths
ixahw = 56;
ixaoff = -nsteps/2;%5*nsteps;%3*nsteps;%nsteps;%nsteps/2;%0;%-nsteps/2;
iyahw = 15;
iyaoff = 0;%12;%4;%24;%16;%0;%
%iyaoff = 16;


% Plot correlation time in s
%tauML = tauML;
yaoff_array = [0];%[16 24];%[24];%[0 4 12 16 24];%[0];%
xaoff_array =  [-nsteps/2];%[nsteps 3*nsteps 5*nsteps];%[-nsteps/2 0 nsteps/2 ];%[0 4 12 16 24];%[-nsteps/2 ];%



iCTR = [-5:4];


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


POSITION = [55 212 723 588];%[67 298 1236 497];%[ 98 457 2410 275];%[1 415 1233 287];%[1 12 473 690];%[51 47 387 692];
PAPERPOSITION = [1 1 4 3];
PAPERSIZE = [ 8.5000   11.0000];

%% Plot fits

fitfunc_str = 'FittingFunctions.CCN2single_fit';
%fitfunc_str = 'FittingFunctions.CCN2single_fit_double_exp';
%fitfunc_str = 'FittingFunctions.CCN2single_fit_double_exp_contrast';
pixel_to_plot_row = [[1:30:nrow]-nrow/2]' ;
pixel_to_plot_col = [[1:30:ncol]-ncol/2]';%[-13 6];
%figNum = 12;

counter = 0;
for pp = 1:size(pixel_to_plot_row,1)
    for pp2 = 1:size(pixel_to_plot_col,1)
        pixel_to_plot_index = [pixel_to_plot_row(pp) pixel_to_plot_col(pp2)];
        
        figNum = 12 + counter;
        Function_display_sos.display_fits(Cdt_struct,runname_titles,fitfunc_str,pixel_to_plot_index,growth_rates_all,figNum);
        counter = counter + 1;
    end
end

%% Plot 2d maps

% Figure 1: time constants
figNum = 100;
clear title_figure;
for ss = 1:size(runname_titles,1)
    %title_figure(ss,:) = {[runname_titles(ss,:)],[ ' tau (sim. units)'],['in log scale']};
    title_figure(ss,:) = {[runname_titles{ss,:} ' tau fast  (sim. units)']};
end
Function_display_sos.make_figure_2D(title_figure,Cdt_struct,'damHW_fast',struct_ranges,...
    iCTR,figNum,'CLIM',[5 10],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,'PAPERSIZE',PAPERSIZE,...
    'plot_averages',1);



% Figure 1: time constants
figNum = 1000;
clear title_figure;
for ss = 1:size(runname_titles,1)
    %title_figure(ss,:) = {[runname_titles(ss,:)],[ ' tau (sim. units)'],['in log scale']};
    title_figure(ss,:) = {[runname_titles{ss,:} 'sigma tau fast (sim. units)']};

end
Function_display_sos.make_figure_2D(title_figure,Cdt_struct,'damHW_fast_sigma',struct_ranges,...
    iCTR,figNum,'CLIM',[5 10],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,'PAPERSIZE',PAPERSIZE,...
    'plot_averages',1);


% Figure 150: time constants
figNum = 150;
clear title_figure;
for ss = 1:size(runname_titles,1)
    %title_figure(ss,:) = {[runname_titles(ss,:)],[ ' tau (sim. units)'],['in log scale']};
    title_figure(ss,:) = {[runname_titles{ss,:} ' tau slow  (sim. units)']};
end
Function_display_sos.make_figure_2D(title_figure,Cdt_struct,'damHW_slow',struct_ranges,...
    iCTR,figNum,'CLIM',[5 10],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,'PAPERSIZE',PAPERSIZE,...
    'plot_averages',1);

% Figure 1: time constants
figNum = 1500;
clear title_figure;
for ss = 1:size(runname_titles,1)
    %title_figure(ss,:) = {[runname_titles(ss,:)],[ ' tau (sim. units)'],['in log scale']};
    title_figure(ss,:) = {[runname_titles{ss,:} 'sigma tau slow (sim. units)']};

end
Function_display_sos.make_figure_2D(title_figure,Cdt_struct,'damHW_slow_sigma',struct_ranges,...
    iCTR,figNum,'CLIM',[5 10],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,'PAPERSIZE',PAPERSIZE,...
    'plot_averages',1);




% Figure 200: Contrast fast
figNum = 200;
% time frame to represent:
%time_frame = 130;
clear title_figure;
for ss = 1:size(runname_titles,1)
    title_figure(ss,:) = {[runname_titles{ss,:} ' Contrast fast']};
    map2D(ss).Contrast = 10.^(Cdt_struct(ss).contrast_fast);
end


Function_display_sos.make_figure_2D( title_figure,map2D,'Contrast',struct_ranges,iCTR,figNum,...
    'CLIM',[0 1],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);


% Figure 250: Contrast slow
figNum = 250;
% time frame to represent:
%time_frame = 130;
clear title_figure;

ss = 1;
title_figure(ss,:) = {[runname_titles{ss,:} ' Contrast']};
map2D(ss).Contrast = 10.^(Cdt_struct(ss).contrast_fast);
for ss = 2:size(runname_titles,1)
    title_figure(ss,:) = {[runname_titles{ss,:} ' Contrast slow']};
    map2D(ss).Contrast = 10.^(Cdt_struct(ss).contrast_slow);
end


Function_display_sos.make_figure_2D( title_figure,map2D,'Contrast',struct_ranges,iCTR,figNum,...
    'CLIM',[-.3 .3],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);



% Figure 2000: Contrast manual
%{
figNum = 2000;
% time frame to represent:
%time_frame = 130;
clear title_figure;
for ss = 1%1:size(runname_titles,1)
    title_figure(ss,:) = {[runname_titles{ss,:} ' Contrast manual']};
    
    for jj = 1:ncol
        for ii = 1:nrow
            map2D(ss).Contrast_manual(ii,jj) = Cdt_struct(ss).Cdt(ii,jj,1);
        end
    end
end


Function_display_sos.make_figure_2D( title_figure,map2D,'Contrast_manual',struct_ranges,iCTR,figNum,...
    'CLIM',[0.0 1.1],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);
%}

% Figure 300: background
%{

figNum = 300;
% time frame to represent:
%time_frame = 130;
clear title_figure;
for ss = 1:size(runname_titles,1)
    title_figure(ss,:) = {[runname_titles{ss,:} ' Background']};
    map2D(ss).back = Cdt_struct(ss).back;
end


Function_display_sos.make_figure_2D( title_figure,map2D,'back',struct_ranges,iCTR,figNum,...
    'CLIM',[0 0.2],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);
%}


return;


%% Plot averages


% theoretical profiles
% if time axis is in seconds:
tauMLth_Qx = .5e10*nsteps./(pi*abs(QX-ixaoff).^2); % Unclear factor
tauMLth_Qx_1 = .25e9*nsteps./(pi*abs(QX-ixaoff)); % Unclear factor
tauMLth_Qy = .2e10*nsteps./(pi*abs(QY).^2); % Unclear factor
tauMLth_Qy_3 = .1e5*nsteps./(pi*abs(QY).^3); % Unclear factor
tauMLth_Qy_2p5 = .03e5*nsteps./(pi*abs(QY).^2.5); % Unclear factor


% average over rows
range_rows = iycen+iyaoff+[-iyahw:iyahw];

for kk = 1:size(runnames,1)
    tauML_1D = mean(Cdt_struct(kk).damHW_fast(range_rows,:),1);%mean(Cdt_struct(ss).damHW_slow(range_rows,:),1);
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
    
    damHW1D_struct(kk).Color = color_rates_allgrowth(kk,:);
    legend_str{kk} = [' growth rate (ML/s) = ' num2str(growth_rates_all(kk))];%[runname ' growth rate (ML/s) = ' num2str(growth_rates_all(kk))];
end


figNum = 2;

Function_display_sos.make_figures_tauML1D(['\tau_s offset in Q_y = ' num2str(iyaoff)],...
    QX,damHW1D_struct,'tauML_1D_Qx',figNum,...
    'iCTR',iCTR,'nsteps',nsteps,'POSITION',POSITION,...
    'XLABEL','Q_x (pixels)','YLIM',[1e6 1e10]);

legend(legend_str);


figNum = 20;

Function_display_sos.make_figures_tauML1D(['QX is positive,  offset in Qy = ' num2str(iyaoff)],...
    Qx_positive,damHW1D_struct,'tauML_1D_Qx_positive',...
    figNum,'tauMLth',tauML_th_QX_positive,'iCTR',iCTR,'nsteps',nsteps,'ixaoff',ixaoff,...
    'XLABEL','Q_x (pixels)','POSITION',POSITION,...
    'XSCALE','log','YSCALE',...
    'log','YLIM',[1e6 1e10]);

legend(legend_str);

figNum = 21;

Function_display_sos.make_figures_tauML1D(['QX is negative,  offset in Qy = ' num2str(iyaoff)],...
    Qx_negative,damHW1D_struct,'tauML_1D_Qx_negative',...
    figNum,'tauMLth',tauML_th_QX_negative,'iCTR',iCTR,'nsteps',nsteps,'ixaoff',ixaoff,...
    'XLABEL','Q_x (pixels)','POSITION',POSITION,...
    'XSCALE','log','YSCALE',...
    'log','YLIM',[1e6 1e10]);

legend(legend_str);


% average over columns
range_cols = ixcen+ixaoff+[-ixahw:ixahw];

for kk = 1:size(runnames,1)
    %kk = ss-1;
    tauML_1D = mean(Cdt_struct(kk).damHW_fast(:,range_cols),2);%mean(Cdt_struct(ss).damHW_slow(:,range_cols),2);

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
    
    damHW1D_struct(kk).Color = color_rates_allgrowth(kk,:);
    legend_str{kk} = [' growth rate (ML/s) = ' num2str(growth_rates_all(kk))];%[runname_title ' growth rate (ML/s) = ' num2str(growth_rates_all(kk))];
end


figNum = 4;

Function_display_sos.make_figures_tauML1D(['\tau_s offset in Q_x = ' num2str(ixaoff)],...
    QY,damHW1D_struct,'tauML_1D_Qy',figNum,'POSITION',POSITION,...
    'XLABEL','Q_y (pixels)','YLIM',[1e6 1e10]);

legend(legend_str);


figNum = 40;

Function_display_sos.make_figures_tauML1D(['QY is positive,  offset in Qx = ' num2str(ixaoff)],...
    Qy_positive,damHW1D_struct,'tauML_1D_Qy_positive',...
    figNum,'tauMLth',tauML_th_Qy_positive,'iyaoff',iyaoff,'XLABEL','Q_y (pixels)',...
    'XSCALE','log','YSCALE','log','POSITION',POSITION,...
    'YLIM',[1e6 1e10]);

legend(legend_str);

figNum = 41;

Function_display_sos.make_figures_tauML1D(['QY is negative,  offset in Qx = ' num2str(ixaoff)],...
    Qy_negative,damHW1D_struct,'tauML_1D_Qy_negative',...
    figNum,'tauMLth',tauML_th_Qy_negative,'iyaoff',iyaoff,'XLABEL','Q_y (pixels)',...
    'XSCALE','log','YSCALE','log','POSITION',POSITION,...
    'YLIM',[1e6 1e10]);
legend(legend_str);


%% Single Q values:

QY_val = 30;
figure(5);
clf;
hold on;
for ss = 1:size(runnames,1)
    tauML_1D = mean(Cdt_struct(ss).damHW_fast(:,range_cols),2);
    damHW1D_struct(ss).tauML_1D_oneQ = tauML_1D(QY_val);
    
    plot(growth_rates_all(ss),damHW1D_struct(ss).tauML_1D_oneQ,color_rates_allgrowth(ss,:),'LineWidth',3.0,'MarkerSize',20);
    
end
xlabel('Growth rates (ML/sim units)');
ylabel('\tau (sim units)');
set(gca,'FontSize',30);
box on;
title(['Time constant vs growth rate at QY = ' num2str(QY_val)]);
set(gcf,'Position',[534    64   769   538]);
set(gca,'YScale','log');
ylim([4e6 3e7]);

QX_val = 30;
figure(6);
clf;
hold on;
for ss = 1:size(runnames,1)
    tauML_1D = mean(Cdt_struct(ss).damHW_fast(range_rows,:),1);
    damHW1D_struct(ss).tauML_1D_oneQ = tauML_1D(QX_val);
    
    plot(growth_rates_all(ss),damHW1D_struct(ss).tauML_1D_oneQ,color_rates_allgrowth(ss,:),'LineWidth',3.0,'MarkerSize',20);
    
end
xlabel('Growth rates (ML/sim units)');
ylabel('\tau (sim units)');
set(gca,'FontSize',30);
box on;
title(['Time constant vs growth rate at QX = ' num2str(QY_val)]);
set(gcf,'Position',[534    64   769   538]);
set(gca,'YScale','log');
ylim([2e7 2e9]);

return;



%% load ihm

computer_script;

if ihm_flag == 0
    
    for ss = 1:size(runnames,1)
        
        path = paths{ss,:};
        runname = [paths{ss,:} runnames{ss,:}];
        runname_title = runname_titles{ss,:};
                
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
        
        if strcmp(runname_title ,[ 'miscut\_32\_37'])
            dt_minML = dtime(end-round(2.50e4/2));%dtime(30);%dtime(end-round(2.5e4/8));%dtime(end-6.3e3);% % dtime(end-1e4);
            
        else
            dt_minML = dtime(end-2.50e4);%dtime(30);%dtime(end-round(2.5e4/8));%dtime(end-6.3e3);% % dtime(end-1e4);
        end
        dt_maxML = dtime(end);%dtime(250);%dtime(end-1.25e4);%
        %idt = dtime >dt_minML;
        idt = (dtime>dt_minML & dtime<dt_maxML);
        
        idt_half =round(size(idt,1)/2);
        
        ihm_struct(ss).ihm = squeeze(ihm(:,:,idt_half));
        
        
        
        
    end
    
    save('./allmiscut_mat_files/ihm_struct_file.mat','ihm_struct')
else
    load('./allmiscut_mat_files/ihm_struct_file.mat','ihm_struct');
end

nrow_orig = size(ihm_struct(1).ihm,1);
ncol_orig = size(ihm_struct(1).ihm,2);



%% Plot ihm


struct_ranges_ihm = struct_ranges;
struct_ranges_ihm.QX = [1:ncol];
struct_ranges_ihm.QY = [1:nrow];
struct_ranges_ihm.row_range = [1:nrow];
struct_ranges_ihm.col_range = round(ncol_orig/2)+[-round(ncol/2):round(ncol/2)-1]; 

% Figure 13: ihm
figNum = 13;
% time frame to represent:
%time_frame = 130;
clear title_figure;
for ss = 1:size(runname_titles,1)
    title_figure(ss,:) = {[runname_titles{ss,:}];[' ihm; g.r = ' num2str(growth_rates_all(ss))]};
    ihm_struct(ss).ihm_square = ihm_struct(ss).ihm(struct_ranges_ihm.row_range,struct_ranges_ihm.col_range,:);
end


Function_display_sos.make_figure_2D_realspace( title_figure,ihm_struct,'ihm_square',struct_ranges_ihm,figNum,...
    'CLIM',[0 35],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);


%% Load and plot the diffraction patterns

 clear time_frame;
computer_script;
 
 if III_flag == 0
     for ss = 1:size(runnames,1)
         path = paths{ss,:};
         runname = [paths{ss,:} runnames{ss,:}];
         runname_title = runname_titles{ss,:};
         
         load([runname '_corr_dt_L0p5.mat'],'III_resize');
         
         time_frame_index = round(size(III_resize,3)/2);
         
         III_struct(ss).III = III_resize(:,:,time_frame_index);
     end
     save('./allmiscut_mat_files/III_struct_file.mat','III_struct')
     
     
 else
     load('./allmiscut_mat_files/III_struct_file.mat','III_struct')
     
 end
 
 
 
 % Figure 10: diffraction pattern
figNum = 14;
% time frame to represent:
%time_frame = [1e4 1e4 5e3 1e3];
clear title_figure;
for ss = 1:size(runname_titles,1)
    iii = 1;
    jjj = 1;
    time_frame(ss) = time_frame_index;
    title_figure(ss,:) = {[runname_titles{ss,:}];[ ' III, @ ' num2str(time_frame(ss),'%.1e') ' g.r = ' num2str(growth_rates_all(ss))]};%[runname_titles(ss,:) ' Diffracted intensity, time frame = ' num2str(time_frame)];
    %map2D(ss).III = damHW_struct(ss).III(:,:,time_frame(ss));
end


Function_display_sos.make_figure_2D( title_figure,III_struct,'III',struct_ranges,iCTR,figNum,...
    'CLIM',[-10 -2],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);



 
 %% Load and plot the Ibar
 clear time_frame;
commputer_script;
 
 if Ibar_flag == 0
     for ss = 1:size(runnames,1)
         path = paths{ss,:};
         runname = [paths{ss,:} runnames{ss,:}];
         runname_title = runname_titles{ss,:};
         
         load([runname '_corr_dt_L0p5.mat'],'Ibar');
         
         time_frame_index = round(size(Ibar,3)/2)
         
         Ibar_struct(ss).Ibar = Ibar(:,:,time_frame_index);
     end
     save('./allmiscut_mat_files/Ibar_struct_file.mat','Ibar_struct')
     
     
 else
     load('./allmiscut_mat_files/Ibar_struct_file.mat','Ibar_struct')
     
 end
 

 
 % Figure 10: diffraction pattern
figNum = 15;
% time frame to represent:
%time_frame = [1e4 1e4 5e3 1e3];
clear title_figure;
for ss = 1:size(runname_titles,1)
    iii = 1;
    jjj = 1;
    time_frame(ss) = time_frame_index;
    title_figure(ss,:) = {[runname_titles{ss,:} ' Ibar, @ ' num2str(time_frame(ss),'%.1e')]};%[runname_titles(ss,:) ' Diffracted intensity, time frame = ' num2str(time_frame)];
    %map2D(ss).III = damHW_struct(ss).III(:,:,time_frame(ss));
end


Function_display_sos.make_figure_2D( title_figure,Ibar_struct,'Ibar',struct_ranges,iCTR,figNum,...
    'CLIM',[-10 -2],...
    'POSITION',POSITION,'PAPERPOSITION',PAPERPOSITION,...
    'plot_averages',0);

