function [NO3] = calc_FLOAT_NO3(spec, ncal, Pcor_flag)
% 
% PURPOSE:
%   Calculate nitrate concentration for each sample in a profile(msg file)
%   for APEX or NAVIS floats. This function calculates nitrate
%   concentration from data supplied as a structure retuned from the
%   function parse_NO3msg.m which includes all needed data.
%   See parse_NO3msg.m for details on output structure
%
% USAGE:
%	NO3 = calc_FLOAT_NO3(spec, ncal,Pcor_flag)
%
% INPUTS:
%	spec      =  a structure built by parse_NO3msg.m
%   ncal      =  a structure built by get_float_cals.m or parseNO3cal.m with NO3 cal data
%   Pcor_flag = 1 to pressure correct ESW, otherwise don't
%
% OUTPUTS:
%   NO3 =   a matrix of nitrate concentrations (umol/L) and other info
%       [SDN, DarkCur,Pres,Temp,Sal,NO3,BL_B,BL_M,RMS ERROR,WL~240,ABS~240]

% HISTORY
%   04/06/2017 Added check for Intensites <= dark current and set to NaN
%   01/18/2022 JP - Added code to check and exclude saturated pixels from calc
%   03/30/2022 JP - Switched over to dLN[ESW]dT ESW temperature correction
%   05/22/2024 ECR (Sea-Bird) - calculate TS correction using adjusted pressure
%                   coordinates using cal.Nitrate_Sensor_Offset as per V1.2+

% *************************************************************************
% TESTING
% cal_file  = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\' ...
%                'DATA\CAL\NO3_CAL\wn1357.cal'];
% cal_file  = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\' ...
%                'DATA\CAL\NO3_CAL\ua20150.cal'];
% cal_file  = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\' ...
%                'DATA\CAL\NO3_CAL\ua20603.cal'];
% ncal      = parseNO3cal(cal_file);
% % spec      = parse_NO3msg('C:\temp\20150.002.isus');
% spec      = parse_NO3msg('C:\temp\20603.006.isus');
% Pcor_flag = 1;

% cal_file  = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\' ...
%                'DATA\CAL\NO3_CAL\ua20620.cal'];
% ncal      = parseNO3cal(cal_file);
% spec      = parse_NO3msg('C:\temp\20620.008.isus');
% Pcor_flag = 1;
% ncal.pres_coef = 0.026;


% cal_file  = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\' ...
%                'DATA\CAL\NO3_CAL\wn1347.cal'];
% ncal      = parseNO3cal(cal_file);
% spec      = parse_NO3msg('C:\temp\1347.010.isus');
% Pcor_flag = 1;
% ncal.pres_coef = 0.026;

% cal_file  = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\' ...
%                'DATA\CAL\NO3_CAL\wn1349.cal'];
% ncal      = parseNO3cal(cal_file);
% spec      = parse_NO3msg('C:\temp\1349.007.isus');
% Pcor_flag = 1;
% ncal.pres_coef = 0.026;

% cal_file  = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\' ...
%                'DATA\CAL\NO3_CAL\ua18796.cal'];
% ncal      = parseNO3cal(cal_file);
% spec      = parse_NO3msg('C:\temp\18796.001.isus');
% Pcor_flag = 1;
% ncal.pres_coef = 0.026;

% ************************************************************************
% CONSTANTS
% ************************************************************************
fig_flag  = 1; % 1 = show nitrate fit
WL_offset = ncal.WL_offset;  % Optical wavelength offset
sat_value = 64500; % Pixel intensity saturation limit count

% TESTING
%WL_offset = 208; %TESTING
% ncal.min_fit_WL = 218.5;
% ncal.max_fit_WL = 250;
%ncal.ESW = ncal.ESW + 0.03*ncal.ENO3;
%ncal.ENO3 = ncal.ENO3 * 0.95;
%ncal.ESW = ncal.ESW*1.05; %TESTING!!!
%ncal.pres_coef = 0.08; % 0.0265 = default
%ncal.CalTemp = ncal.CalTemp - 2.5;

% ************************************************************************
% CALCULATE RAW ABSORBANCE 
% CALIBRATION DATA MEASURED FROM ~190 to ~400 nm IN THE LAB. FLOAT CAN'T
% STORE ALL WAVELENGTHS SPECTROPHOTOMETER MEASURES SO ONLY A SUBSET IS
% RETURNED. 
% ************************************************************************
d = spec; clear spec

% MAKE SOME CALIBRATION ADJUSTMENTS BASED ON DATA FLAGS IN CAL STRUCTURE
% 0 based count - pixel registration to wavelength is off by one for some
% instruments. if 217 was pixel 30 it is now pix 31. DOUBLE CHECK THIS!!!

if ncal.pixel_base == 0 
    d.spectra_pix_range = d.spectra_pix_range +1;
    d.pix_fit_win = d.pix_fit_win + 1;
    disp('Pixel registration offset by +1')
end

% CHOOSE REGULAR or SEAWATER DARK CURRENT FOR CALCULATION
DC = d.DC; % default
if ncal.DC_flag == 0
    DC = d.SWDC;
end

% GET PIXEL FIT WINDOW FOR SUNA FLOATS BASED ON WL FIT RANGE
% SUNA FLOATS DON'T HAVE THIS AS A HEADER LINE
if isnan(sum(d.pix_fit_win)) && isfield(d,'WL_fit_win')
    pfit_low = find(ncal.WL >= d.WL_fit_win(1),1,'first');
    pfit_hi  = find(ncal.WL <= d.WL_fit_win(2),1,'last');
    d.pix_fit_win = [pfit_low pfit_hi];
    clear pfit_low pfit_hi
end

% ************************************************************************
%  DO NOT USE PIXEL FIT RANGE IN CAL FILE! USE 1st >= 217 and
%  LAST <= 240 unless cal file has updated fit window
% ************************************************************************
% ADJUST FIT WINDOW IF IT HAS BEEN UPDATED IN THE CAL STRUCTURE
% OTHERWISE USE DEFAULT WINDOW >=217 & <=240
if isfield(ncal, 'min_fit_WL') && ~isempty(ncal.min_fit_WL)
    d.pix_fit_win(1) = find(ncal.WL >= ncal.min_fit_WL,1,'first');
else
    d.pix_fit_win(1) = find(ncal.WL >= 217,1,'first');
end

if isfield(ncal, 'max_fit_WL') && ~isempty(ncal.max_fit_WL)
    d.pix_fit_win(2)  = find(ncal.WL <= ncal.max_fit_WL,1,'last');
else
    d.pix_fit_win(2)  = find(ncal.WL <= 240 ,1,'last');
end

% SAMPLE PIXEL RANGE INDEX IDENTIFIED FROM SAMPLE LINE
saved_pixels = d.spectra_pix_range(1): d.spectra_pix_range(2); % index

% SIZE OF SAMPLE INTENSITY SPECTRA MATRIX(samples X wavelengths)
[rows,cols]  = size(d.UV_INTEN);

% SUBSET CALIBRATION DATA OVER PIXEL RANGE IDENTIFIED FROM SAMPLE LINE
% THEN TRANSPOSE TO ROWS & REPEAT ROWS SO SAME SIZE AS UV_INTEN
REF  = ones(rows,1) * ncal.Ref(saved_pixels)';  % Reference intensities
WL   = ones(rows,1) * ncal.WL(saved_pixels)';   % Wavelengths

ESW  = ncal.ESW(saved_pixels)';    % Extinction coefs SW
TSWA = ncal.TSWA(saved_pixels)';   % Extinction coefs TS
ENO3 = ncal.ENO3(saved_pixels)';   % Extinction coefs NO3 (array)

if isfield(ncal,'EHS') % NOT USED YET BUT CARRY ALONG
    EHS  = ncal.EHS(saved_pixels)';  % Extinction coefs HS  (array)
end

% CHECK FOR SATURATED SAMPLE INTENSITY VALUES JP 01/18/2022
tPIX_SAT = d.UV_INTEN > sat_value; % logical matrix
if sum(tPIX_SAT(:)) > 0
    disp('WARNING: Saturated sample pixel intensities detected in profile');
    disp('Saturated values will be excluded. Nitrate estimates may be compromised');
    d.UV_INTEN(tPIX_SAT) = NaN;
end
    
% SUBTRACT DARK CURRENT INTENSITIES
UV_INTEN_SW = d.UV_INTEN - DC * ones(1,cols);

% CHECK FOR UV_INTEN_SW <= 0 (DC >= UV_INTEN_SW) & SET TO NaN
tINTEN_SW = UV_INTEN_SW <= 0;
if sum(tINTEN_SW(:)) > 0
    UV_INTEN_SW(tINTEN_SW) = NaN;
    disp(['UV Intensities <= DC found for ',d.file_name, ...
        '. Setting these low intensities to NaN']);
    %clear tINTEN_SW
end
    
% ABSORBANCE SPECTRA FOR ALL SAMPLES IN PROFILE
ABS_SW = -log10(UV_INTEN_SW ./ REF);

% ************************************************************************
% CALCULATE TEMPERATURE CORRECTED ABSORBANCES
% ************************************************************************

% Interpolate T,S to pressure offset by ncal.Nitrate_Sensor_Offset
Pint = d.P + ncal.Nitrate_Sensor_Offset; 
Tint = interp1( d.P, d.T, Pint); 
Sint = interp1( d.P, d.S, Pint); 

% Build matrices for calculations all size of spectra - avoid some loops 
% calculate corrected absorbances for all samples all at once
ctd_temp  = Tint * ones(1,cols);           % measured sw temperature, C
ctd_sal   = Sint * ones(1,cols);           % measured sw salinity
cal_temp  = ones(rows,cols) * ncal.CalTemp;  % temp isus calibrated at, C

% Temperature / Salinity Correction Algorithm Processing
switch (ncal.TSalgorithm) 
    case 'None'
        fprintf(1, 'No TS Correction\n');

    case 'TCSS_Sakamoto2009'
        % SAKAMOTO 2009 TCOR
        % Refer to Computing DAC Nitrate v1.1

        A = 1.1500276; % Coefficients from Sakamoto et al., 2009
        B = 0.02840;
        C = -0.3101349;
        D = 0.001222;

        % Sakamoto temperature correction "f function" at ctd temperature
        F_sw = (A + B .* ctd_temp) .* exp((C + D .* ctd_temp) .* ...     % Eqn. 3
               (WL - WL_offset));
        
        % Sakamoto temperature correction "f function" at calibration temperature
        % WL_offset
        F_cal = (A + B .* cal_temp) .* exp((C + D .* cal_temp) .* ...    % Eqn. 3
                (WL - WL_offset));
        
        % Calculate in situ sea water extinction coefficients (mostly due to Br)
        ESW_in_situ = ESW .* F_sw ./ F_cal;      % Eqn. 2


    case 'TCSS_Plant2023'
        % UPDATED ESW TEMPERATURE CORRECTION ALGORITHM TO dLN[ESW]/dT 03/30/2022 JP
        % (03/29/2022 JP)EQUITECH PROBES ONLY with 3 anamolously low float experiments
        % removed 05/20/2016, 11/21/2016, 10/09/2017 (P CHAMBER EXPERIMENTS?)
        Tcorr_coef  = [1.27353e-07 -7.56395e-06 2.91898e-05 1.67660e-03 1.46380e-02];
        Tcorr       = polyval(Tcorr_coef,(WL - WL_offset)) .* (ctd_temp - cal_temp);
        %Tcorr       = polyval(Tcorr_coef,(WL - WL_offset)) .* (ctd_temp - cal_temp)*1.03; % TESTING        
        ESW_in_situ = ESW .* exp(Tcorr);

    otherwise
        error('Unknown TS correction algorithn %s', cal.TSalgorithm);
end


if (~strcmp(ncal.TSalgorithm, 'None'))
    % Add pressure corrrection if flag is on
    if Pcor_flag == 1
        pres_term   = (1- Pint ./ 1000 .* ncal.pres_coef) * ones(1,cols);
        ESW_in_situ = ESW_in_situ .* pres_term;
    end
    
    % Estimated SW (bromide) absorbances from salinity
    ABS_Br_tcor      = ESW_in_situ .* ctd_sal; % in situ bromide absorbances
    
    % Remove seawater (bromide) absorbance from sample absorbance (TCSS method)
    % This leaves nitrate + baseline absorbance by difference
    ABS_cor = ABS_SW - ABS_Br_tcor  ;    % Eqn. 4
end
clear ctd_temp cal_temp ctd_sal REF pres_term

% ************************************************************************
% CALCULATE THE NITRATE CONCENTRATION, BASELINE SLOPE AND INTERCEPT OF THE
% BASELINE ABSORBANCE
% This calulation occurs over the "fit window" which is usually a subset
% of the "spectra window" which is a subset of the calibration spectra
% window. The calibration window is the full wavelength range of the
% spectraphotometer
% ************************************************************************

% SUBSET OF FIT WINDOW PIXELS
t_fit     = saved_pixels >= d.pix_fit_win(1) & ...
            saved_pixels <= d.pix_fit_win(2); % fit win pixel indicies
      
Fit_WL      = WL(1,t_fit)'; %n x 1
Fit_ESW     = ESW(t_fit)';  %n x 1   % Used for ncal.TSalgorithm = 'None';
Fit_TSWA    = TSWA(t_fit)'; %n x 1   % Used for ncal.TSalgorithm = 'None';
Fit_ENO3    = ENO3(t_fit)'; % n x 1
Ones        = ones(size(Fit_ENO3));

% /1000 & /100 for numerical stability or equalizing weights of unknowns
% Need to /100 & /1000 at the end to get back to actual baseline and intercept
if (strcmp(ncal.TSalgorithm, 'None'))
    M     = [Fit_ENO3 Ones/100 Fit_WL/1000 Fit_ESW Fit_TSWA]; % WL x 5
    Fit_ABS_cor = (ABS_SW(:,t_fit))';
else
    M     = [Fit_ENO3 Ones/100 Fit_WL/1000]; % WL x 3
    Fit_ABS_cor = (ABS_cor(:,t_fit))'; % NO3 + BL. Transposed, now WL x sample
end

% Matrix pseudo-inverse
M_INV = pinv(M);

colsM = size(M,2); % Columns may need to be added if more unknowns
NO3 = ones(rows, colsM +3)* NaN; % # samples x (# fit parameters + QC values)
S = ones(rows, 1)* NaN; % # samples x 1
TS = ones(rows, 1)* NaN; % # samples x 1

for i = 1:rows
    tg = ~tPIX_SAT(i,t_fit)'; % non saturated intensities for sample i, n x 1
    m_inv = M_INV; % default unless saturated pixels
    if sum(~tg,1) > 0 % saturated intensities exist
        m_inv = pinv(M(tg,:)); % subset to non saturated WL's in M
    end


    if (strcmp(ncal.TSalgorithm, 'None'))
        % NO3 = [ baseline intercept,  slope, ESW, TSWA]
        NO3(i,1:colsM) = m_inv * Fit_ABS_cor(tg,i);
        NO3(i,2) = NO3(i,2)/100; % baseline intercept
        NO3(i,3) = NO3(i,3)/1000; % baseline slope

        S(i,1)   = NO3(i,4);  % Salinity
        TS(i, 1) = NO3(i,5); % TSWA
        
        ABS_BL  = WL(1,:) .* NO3(i,3) + NO3(i,2); % Baseline absorbance
        ABS_S_EXP = ESW .* S(i,1); 
        ABS_TS_EXP = ESW .* TS(i,1); 
        ABS_NO3 = ABS_SW(i,:) - ABS_BL - ABS_S_EXP - ABS_TS_EXP; % Calculated nitrate absorbance
        ABS_NO3_EXP = ENO3 .* NO3(i,1);
        
        
        FIT_DIF = ABS_SW(i,:) - ABS_BL - ABS_NO3_EXP - ABS_S_EXP - ABS_TS_EXP;
    else
        % NO3 = [ baseline intercept, and slope]
        NO3(i,1:3) = m_inv * Fit_ABS_cor(tg,i);  % (#fit params x #WL) X (#WL x #sammples) = (#fit params x #samples)
        NO3(i,2) = NO3(i,2)/100; % baseline intercept
        NO3(i,3) = NO3(i,3)/1000; % baseline slope
        
        ABS_BL  = WL(1,:) .* NO3(i,3) + NO3(i,2); % Baseline absorbance
        ABS_NO3 = ABS_cor(i,:) - ABS_BL; % Calculated nitrate absorbance
        ABS_NO3_EXP = ENO3 .* NO3(i,1);
        
        FIT_DIF = ABS_cor(i,:) - ABS_BL -ABS_NO3_EXP;
    end
    FIT_DIF_tfit = FIT_DIF(t_fit);
    RMS_ERROR = sqrt((sum(FIT_DIF_tfit(tg').^2)./sum(t_fit(tg))));
    ind_240  = find(Fit_WL <= 240,1,'last');
    ABS_240  = [Fit_WL(ind_240), Fit_ABS_cor(ind_240)]; % WL ~240 & ABS
    
    NO3(i,colsM+1:colsM+3) = [RMS_ERROR ABS_240]; % RMS WL~240 ABS~240
    
    
    % ********************************************************************
    % TEST -- PLOT SPECTRA
    % ********************************************************************
    %     if RMS_ERROR > 0.003 || ABS_240(2) > 0.8 % ONLY PLOT IF NOT GOOD
    if fig_flag == 1 && ~all(isnan(ABS_NO3))
        % clf
        figure(1);clf;
        set(gcf, 'Units', 'normalized');
        set(gcf,'Position',[0.17 0.23 0.52 0.55])
        a1 = axes('Position', [0.13 0.12 0.775 0.59]);
        plot(WL(1,:), ABS_SW(i,:), 'k-','LineWidth', 2) % sw spectrum
        hold on
        x_lim = [min(WL(1,:)) max(WL(1,:))];
        xlim(x_lim)
        y_lim = ylim;
        set(gca, 'Linewidth',3)
        xlabel('Wavelength, nm', 'FontSize', 14)
        ylabel('Absorbance','FontSize', 14)
        set(gca, 'FontSize', 14)
        
        if (strcmp(ncal.TSalgorithm, 'None'))
            plot(WL(1,:), ABS_BL, 'g-','LineWidth', 2) % BL
            plot(WL(1,:), ABS_S_EXP, 'c-','LineWidth', 2) % BL
            plot(WL(1,:), ABS_TS_EXP, 'm-','LineWidth', 2) % BL
            plot(WL(1,:), ABS_SW(i,:), 'ro','LineWidth', 1) % NO3+BL observed
            plot(WL(1,:), ABS_NO3_EXP + ABS_S_EXP + ABS_TS_EXP + ABS_BL, 'r-','LineWidth', 2) % NO3+BL exp
        else
            plot(WL(1,:), ABS_BL, 'g-','LineWidth', 2) % BL
            plot(WL(1,:), ABS_Br_tcor(i,:)+ABS_BL, 'b-','LineWidth', 2) % Br spectrum
            plot(WL(1,:), ABS_cor(i,:), 'ro','LineWidth', 1) % NO3+BL observed
            plot(WL(1,:), ABS_NO3_EXP +ABS_BL, 'r-','LineWidth', 2) % NO3+BL exp
        end

        
        y_lim = ylim;
        ylim(y_lim)
        
        fit_patch = fill([Fit_WL(1),Fit_WL(1), Fit_WL(end),Fit_WL(end)], ...
            [y_lim(1), y_lim(2), y_lim(2), y_lim(1)],'y');
        set(fit_patch,'FaceAlpha',0.20, 'LineStyle','--','EdgeColor', ...
            [160 160 160]/256, 'LineWidth',2)
        
        s1 = {['Float:', d.float, '   Cast:',d.cast, '   Depth:', ...
            num2str(d.P(i),'%0.1f')],['FIT WINDOW (',num2str(Fit_WL(1), ...
            '%4.1f'),' - ', num2str(Fit_WL(end),'%4.1f'),')'], ...
            ['NO3 = ', num2str(NO3(i,1),'%2.1f'), ' ?M'], ...
            ['RMS error = ', num2str(RMS_ERROR,'%2.3e')], ...
            ['CTD T = ', num2str(d.T(i),'%2.1f')], ...
            ['CAN T = ', num2str(d.CAN_T(i),'%2.1f')]};
        
        text (WL(1,1) + (WL(1,end) - WL(1,1))* 0.4, ...
            y_lim(1) + (y_lim(2) - y_lim(1))* 0.95,...
            s1,'FontSize', 16, 'HorizontalAlignment','center', ...
            'VerticalAlignment', 'Top')
    
        set(gca, 'box', 'on')

        set(get(get(fit_patch,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off'); % Don't show patch
        legend('Sample','Baseline', 'Br+BL', 'NO3+BL obs', 'NO3+BL fit', ...
            'LOcation','NorthEast')
        legend('boxoff')
        
        hold off
        
        % PLOT RESIDUALS
        a2 = axes('Position', [0.13 0.75 0.775 0.20]);
        FIT_DIF_E3 = FIT_DIF*1000;
        max_y = max(abs(FIT_DIF_E3),[],'omitnan');
        a2_ylim = [-max_y - 0.1*max_y, max_y + 0.1*max_y];
        ylim(a2_ylim)
        xlim(x_lim)
        ylabel({'Obs - Fit','x 10^{3}'},'FontSize', 14)
        set(a2,'Xticklabel','','box', 'on', 'FontSize', 14, 'Linewidth',3);
        hold on
        
        dif_patch = fill([Fit_WL(1),Fit_WL(1), Fit_WL(end),Fit_WL(end)], ...
            [a2_ylim(1), a2_ylim(2), a2_ylim(2), a2_ylim(1)],'y');
        set(dif_patch,'FaceAlpha',0.20, 'LineStyle','--','EdgeColor', ...
            [160 160 160]/256, 'LineWidth',2)
        
        plot(x_lim,[0 0],'--','color','k' ,'LineWidth', 2)
        plot(WL(1,:), FIT_DIF_E3, 'ko','LineWidth', 1, ...
            'MarkerFaceColor','r') % Resid * 1000
        
        hold off
        pause
    end
end

% ************************************************************************
% Create an output matrix and plot the data
% [SDN, DarkCur,Pres,Temp,Sal,NO3,BL_B,BL_M,RMS ERROR,WL~240,ABS~240]







