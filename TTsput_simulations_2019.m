% Two-time analysis of March 2018 TiO2 data
% Read and analyze dataset(s)
% March 17, 2018 GBS
% just use equil 2-time, average over a time range

 %startup;

addpath(genpath(['/Users/ialmazn/Box' ' Sync/XPCS_sputtering_ZnO_TiO2_2017_2019/XPCS_analysis_on_thefly']));
%addpath('../XPCS_analysis_on_thefly/common/');
% General flag to decide what to do:
%0 means read again the data; 
%1 means read previously saved data;

skip = 0; 


% General flag to initialize the parameters for the 'equilibrium' or the
% 'growth' or 'power_series' or 'only_data' (for a test) data:
flag_equil_or_growth = 'only_data';

% simulation flag:
sim_flag = 1;

[plotsmooth,plotall,plotorig,plotnew,ixplotmin,ixplotmax,...
    pcolorCC2,plothalf,plotfull,plotcorrorig,plotcorrnew,...
    plotcontrast,pcolorextra,posflag,fitorig,fitnew, altnorm,plotdeltat]...
    = XPCS_initialize_parameters.TTsput_flags();

[DOCU0,DOCU1,LOGFLAG,delay,XFRAC] = XPCS_initialize_parameters.TTparameters_allruns();

[TCV,nT,iiT,specfilenameM,SCNstrM,XCENV,YCENV,XWIDV,YWIDV,...
    ymax,tminv,tmaxv,clrs, DXImin, ...
    DXImax,NRY,PWV] = XPCS_initialize_parameters.TTparameters_singlerun(flag_equil_or_growth);

[fitrange_time_iiT,pin_iiT,dp_iiT,qfitrange ] ...
    = XPCS_initialize_parameters.TTparameters_fit_singlerun(flag_equil_or_growth);


[hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,...
    tbin_allT,CWID_allT,N_degree_allT,mask] =  XPCS_initialize_parameters.TTparameters_2timecorr_calc(flag_equil_or_growth);

[POSITION,PAPERPOSITION,FONTSIZE,CMAX,CLIM,XCOLlabel,YROWlabel,...
                AXISdet,DOCUclim,INFOstr] = XPCS_initialize_parameters.TTplot_parameters();


 [ImageJ,Xstepcol,SINGLE,BKG,...
 scanflag,imname,p_image,ending,POINTSUMS] = XPCS_initialize_parameters.TTsput_read_ini();

[D_ds,kvector,lambda,pixel_size,th_Bragg,SplitQ_CTR,t_osc] = XPCS_initialize_parameters.TT_experiment();

% set the scans that you want to read
iTV = [1:numel(iiT)];

% in which direction we will do the analysis:

 flag_row_col = 'row'; % 'row' = columns are fixed and rows are variable; 'col' = rows are fixed and columns are variable
 num_col_or_row_V = [1]; %[1:9] % the number of columns if 'row' is selected or the number of rows if "column" is selected
    


%_._._._._._._._._._._._._._._._._._._._._._._._._.._._._._._._._._._.._.._._._______________________
if ~skip
    disp('Recalculating 2-time from experimental datasets.');
    %clear Read_Allscans;

    for iT = iTV
        % read data
        
         warning('off','all');
        [Read_Allscans(iT).IIstruct] = XPCS_read_data_2018_08.TTsimulation_read(iiT(iT),TCV,specfilenameM,SCNstrM,DOCU0,DOCU1,ImageJ,POINTSUMS); % reads a dataset, makes IInormb
        warning('on','all');
        
        
        
        
        %%{
        %Save the data
        Singlescan_struct = Read_Allscans(iT);
        
        Singlescan_namefile = ['Scan' SCNstrM(iT,:) '.mat'];
        save(Singlescan_namefile,'Singlescan_struct','-v7.3');
        %}
        
        
    end
    
else
    
    disp('Loading previously stored scans.');
    for iT = iTV
        
        %SCNstrM_numer = str2num(SCNstrM(iT,:)) ;  
        Singlescan_namefile = ['Scan' SCNstrM(iT,:) '.mat'];
        load(Singlescan_namefile);
        %load('/Users/ialmazn/Box Sync/2017_12_huber_sput_TiO2/Scan059.mat')
        Read_Allscans(iT).IIstruct = Singlescan_struct.IIstruct;
    end
    
end
      

for iT = iTV
    % Initialize ROIs:
    [Allscans(iT).ROIS_struct] = XPCS_read_data_2018_08.TTsput_prepare_ROIS(iiT(iT),XCENV,YCENV,XWIDV,YWIDV,ymax);
    
    
    % Calculate Q-range
     [Allscans(iT).Qval_struct] = XPCS_analysis.calculate_qval(XCENV(iT),YCENV(iT),[1:Read_Allscans(iT).IIstruct.Nc],[1:Read_Allscans(iT).IIstruct.Nr],sim_flag);

    
    
    %Initialize the range of the 2 times correlation functions:
    [Allscans(iT).itt_range_struct] = XPCS_analysis.prepare_subrange_for2corr(iT,Read_Allscans(iT).IIstruct.timeX,XCENV,YCENV,XWIDV,YWIDV,tminv,tmaxv);
    
    % Bin IInorm
    [Allscans(iT).IIbin_struct] = XPCS_analysis.bin_scans_in_time(Read_Allscans(iT),Allscans(iT),iT,tbin_allT,ImageJ,POINTSUMS);
    
    % Calculate 2 time correlaiont functions
    flag_mean_or_poly = 'mean';
    [Allscans(iT).CCN2_struct,Allscans(iT).IIbin_struct] = XPCS_analysis.calc_2time_corr(Allscans(iT).IIbin_struct,iT,N_degree_allT,flag_mean_or_poly);
    
    [Allscans(iT).CCN2avg_struct,Allscans(iT).qvector_struct] = XPCS_analysis.from_CCN2V_to_CCN2avg(Allscans(iT),iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,sim_flag);
    
    % calculate the IInormbb_ref averaged per box
    %[Allscans(iT).IInormbb_ref_avg] = XPCS_analysis.calc_IInormb_avg_inpixels(Allscans(iT).IIbin_struct,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,D_ds,kvector,pixel_size,th_Bragg,mask);
  
    
    %DisplayFunctions_XPCS.display_IInormbbref(Allscans(iT).IIbin_struct,Allscans(iT).CCN2avg_struct.boxcenterrc,50+iT,AXISdet,D_ds,kvector,pixel_size,th_Bragg);

     % plot the CTR images in pixels:
    %DisplayFunctions_XPCS.display_CCN2avg(Allscans(iT).CCN2avg_struct,num_col_or_row,flag_row_col,40+iT+num_col_or_row,AXISdet);
    QvalFlag = 0;
    fignum = 90;
    DisplayTTM_plot.display_IInormb(Allscans(iT).IIbin_struct,Allscans(iT).IIbin_struct.IInormbb,'Log 10 binned IInormb',QvalFlag,fignum+iT,ImageJ,SINGLE,XCOLlabel,YROWlabel,AXISdet,INFOstr,sim_flag,CLIM);
    DisplayFunctions_XPCS.display_grid_CCN2avg(Allscans(iT).CCN2avg_struct.ittccen,Allscans(iT).CCN2avg_struct.ittrcen,Allscans(iT).CCN2avg_struct.Ncq_Nrq(1),Allscans(iT).CCN2avg_struct.Ncq_Nrq(2),QvalFlag,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,fignum+iT,ImageJ,sim_flag);
    %DisplayFunctions_XPCS.display_mask_CCN2avg(Allscans(iT).CCN2avg_struct.ittccen,Allscans(iT).CCN2avg_struct.ittrcen,Allscans(iT).CCN2avg_struct.Ncq_Nrq(1),Allscans(iT).CCN2avg_struct.Ncq_Nrq(2),QvalFlag,iT,mask,fignum+iT,ImageJ,D_ds,kvector,pixel_size,th_Bragg);
    
    % plot the CTR images in reciprocal space units
    QvalFlag = 1;
    fignum = 80;
    DisplayTTM_plot.display_IInormb(Allscans(iT).IIbin_struct,Allscans(iT).IIbin_struct.IInormbb,'Log 10 binned IInormb',QvalFlag,fignum+iT+100,ImageJ,SINGLE,XCOLlabel,YROWlabel,AXISdet,INFOstr,sim_flag,CLIM);
    DisplayFunctions_XPCS.display_grid_CCN2avg(Allscans(iT).CCN2avg_struct.ittccen,Allscans(iT).CCN2avg_struct.ittrcen,Allscans(iT).CCN2avg_struct.Ncq_Nrq(1),Allscans(iT).CCN2avg_struct.Ncq_Nrq(2),QvalFlag,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,fignum+iT+100,ImageJ,sim_flag);
    %DisplayFunctions_XPCS.display_mask_CCN2avg(Allscans(iT).CCN2avg_struct.ittccen,Allscans(iT).CCN2avg_struct.ittrcen,Allscans(iT).CCN2avg_struct.Ncq_Nrq(1),Allscans(iT).CCN2avg_struct.Ncq_Nrq(2),QvalFlag,iT,mask,fignum+iT,ImageJ,D_ds,kvector,pixel_size,th_Bragg);

    
    [Allscans(iT).CCN2S_struct] = XPCS_analysis.from_CCN2avg_to_CCN2S(Allscans(iT).CCN2avg_struct);
    
    
end


for iT = iTV
    
    for num_col_or_row = num_col_or_row_V
        [Allscans(iT).CCN2S_struct] = XPCS_analysis.fit_CCN2S_with_leasqr(Allscans(iT).CCN2S_struct,iT,pin_iiT,dp_iiT,'FittingFunctions.CCN2single_fit',[1:1:round(fitrange_time_iiT(iT)*Allscans(iT).IIbin_struct.Ntb)],num_col_or_row,flag_row_col,200);
        
        figh = 1000*iT+100*num_col_or_row;
        
        % display the fitted parameters vs the qvector
        %[Allscans(iT).fitres_CCN2S(num_col_or_row)] = DisplayFunctions_XPCS.display_fit_result(Allscans(iT).CCN2S_struct,num_col_or_row,flag_row_col,[1:numel(Allscans(iT).CCN2S_struct.scancq(num_col_or_row).scanrq(1).CCNdtV_fit.pout)-1],figh);
          
        
        %%{
        switch flag_row_col
            
            case 'row'
                [Allscans_fit(iT).pout_all(num_col_or_row).pout_struct] = DisplayResults_Fit.display_fit_result(Allscans(iT).CCN2S_struct,num_col_or_row,flag_row_col,[1:numel(Allscans(iT).CCN2S_struct.scancq(num_col_or_row).scanrq(1).CCNdtV_fit.pout)-1],figh);
            case 'col'
                [Allscans_fit(iT).pout_all(num_col_or_row).pout_struct] = DisplayResults_Fit.display_fit_result(Allscans(iT).CCN2S_struct,num_col_or_row,flag_row_col,[1:numel(Allscans(iT).CCN2S_struct.scancq(1).scanrq(num_col_or_row).CCNdtV_fit.pout)-1],figh);
        end
        %}
        
       % Display the 2time correlation function and the fit
        DisplayFunctions_XPCS.display_IIbinstruct_CCN2S_CCN2avg(Allscans(iT),num_col_or_row,flag_row_col,101+iT,[2:round(Allscans(iT).IIbin_struct.Ntb/2)],AXISdet);
        
        
        
        
    end
    
    
end

% Display the parameters resulting from the fit

% Display the fit of the averaged two time correlation function:
ppvector = [3];
for iT = iTV
    Allscans(iT).pcolor_struct = DisplayResults_Fit.display_fit_result_allcol_row_pcolor(Allscans(iT),flag_row_col,[ppvector],1000*iT);
    
    %Allscans(iT).pout_avg = XPCS_analysis.calc_average_row_col( Allscans(iT).pcolor_struct,ppvector);
    
    %DisplayFunctions_XPCS.display_fit_average_results(Allscans(iT).pout_avg,ppvector,2000+iT);
    % Allscans_fit(iT).pout_all_int = DisplayFunctions_XPCS.display_fit_result_allcol_row_pcolor_weight(Allscans_fit(iT),flag_row_col,[ppvector],1000);
    
    if plotall
        
        for num_col_or_row = 1
            DisplayFunctions_XPCS.display_IIbinstruct_CCN2S_CCN2avg(Allscans(iT),num_col_or_row,flag_row_col,101+iT,[2:round(Allscans(iT).IIbin_struct.Ntb/3)],AXISdet);
           % DisplayFunctions_XPCS.display_IInormb_avg(Allscans(iT).IInormbb_ref_avg,num_col_or_row,flag_row_col,300,AXISdet)

        end
        
    end
end


%%{


%disp('Fit the data');

 % Fit the taus at different temperatures 
pp_index = 3;

for jj = num_col_or_row_V
    num_col_or_row = num_col_or_row_V(jj);
    [pout_struct] = DisplayResults_Fit.display_fit_result_log(Allscans(iTV),pp_index,num_col_or_row,flag_row_col,700+jj);
    [pout_struct] = DisplayResults_Fit.append_q_law(700+jj,pout_struct,SplitQ_CTR,60,'r');
  
   % [pout_struct] = DisplayFunctions_XPCS.append_q_law(600+jj,pout_struct,SplitQ_CTR,600,'g');

end

%DisplayFunctions_XPCS.display_fit_result_allcol_row_inonefig(pout_struct,flag_row_col,pp_index,'rgkbm',600+jj+1);

% fit pout(3,:)= tau to a power law 1/x
%{
for iT = iTV
    
    figh = 400+iT+num_col_or_row_V;
    [Allscans_fit(iT).pout,Allscans_fit(iT).sigma] = DisplayFunctions_XPCS.display_fit_result(Allscans(iT).CCN2S_struct,num_col_or_row_V,flag_row_col,[3],figh);
   
    
    figure;
   
    pout_struct(iT).pout_fit = XPCS_analysis.fit_tau_with_leasqr(pout_struct(iT),'FittingFunctions.powerlaw_tau',qfitrange.nu);
      
    
    figure(601);
    subplot(1,numel(iiT),iT);
    %title(['#' SCNstrM(iT,:)]);
    hold on;
    plot(pout_struct(iT).pout_fit.pout.x,pout_struct(iT).pout_fit.pout.fitfunc,'r','LineWidth',3.0);
    
    %figure(602);
    %hold on;
    %plot(pout_struct(iT).pout_fit.pout.x,pout_struct(iT).pout_fit.pout.fitfunc,'r','LineWidth',3.0);

end

% display the velocity vs temperature:
island_size = 2e3; % in Angstroms
[vel_struct] = DisplayFunctions_XPCS.display_vel_vs_temp(TCV,pout_struct,'Tau vs temperature for equilibrium',island_size,700);


island_size = 2e3; % in Angstroms
%[vel_struct] = DisplayFunctions_XPCS.display_vel_vs_temp(PWV,pout_struct,'Tau vs power',island_size,800);

DisplayFunctions_XPCS.display_growth_vs_power(PWV,vel_struct,[1:numel(PWV)], '800 C',1000)


% display the velocity vs the power
%power_vect = [10 20 50];
%DisplayFunctions_XPCS.display_growth_vs_power(power_vect,vel_struct,[5 6 7],'900 C',800);
%}
%{
 for iT = iTV
 % Initialize ROIs:
    [Allscans(iT).ROIS_struct] = XPCS_read_data.TTsput_prepare_ROIS(iiT(iT),XCENV,YCENV,XWIDV,YWIDV,ymax);
    
    
    % Calculate Q-range
     [Allscans(iT).Qval_struct] = XPCS_analysis.calculate_qval(XCENV(iT),YCENV(iT),[1:Read_Allscans(iT).IIstruct.Nc],[1:Read_Allscans(iT).IIstruct.Nr],D_ds,kvector,pixel_size,th_Bragg);

    
    
    %Initialize the range of the 2 times correlation functions:
    [Allscans(iT).itt_range_struct] = XPCS_analysis.prepare_subrange_for2corr(iT,Read_Allscans(iT).IIstruct.Xamount,XCENV,YCENV,XWIDV,YWIDV,tminv,tmaxv);
    
    TTM_plot1; 
 end
%}



