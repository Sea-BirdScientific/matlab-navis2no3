%% Process ISUS data from a NAVIS float
%
% Process nitrate spectra from reduced binary spectra in SBS NAVIS ISUS 
% message files.
% 
% Requires:
% 
% TEOS-10 GSW Toolkit,
% findextension.m
% parse_NO3msg.m
% parseNO3cal.m
% calc_FLOAT_NO3.m
% plotSpectrum3D.m
%
% 05/22/23 - ECR (Sea-Bird)   
% 
 
%% Pathnames
srcdir = '/Users/ericrehm/OneDrive - Danaher/Sea-Bird/SUNA/src/NAVIS2NO3_MBARI/matlab-navis2no3/';
datapath = '/Volumes/nasshare/common/dev/D1046-NAVIS/data/';
% calpath = '/Volumes/nasshare/common/calibration/SUNA-V2/';
calpath = fullfile(srcdir, 'caltemp/SUNA-V2/');

doProfilePlots = true;
doRawSpectralPlots = true;

%% load the data from .isus message files: t, s, p, spectra, etc.
% floatid = '1115'; calVariant = 'Cal_01 SW';
% floatid = '1114'; calVariant = 'Cal_01 SW';
% floatid = '0060'; calVariant = 'Cal_02SW';
% floatid = '0061'; calVariant = 'Cal_01 SW';
% floatid = '0062'; calVariant = 'Cal_01 SW';
floatid = '0065'; calVariant = 'Cal_01 SW';

fprintf(1, 'Float ID: %s, calVariant: %s, ', floatid,  calVariant);

% DEV DEV Is time-consuming issus message data retrieval already cached?
% cacheName = "isus_cache_"+floatid+".mat";
% if (~exist(cacheName,'file'))    
%     %% load the isus file from Navis float;
%     % define the folder
%     data_target_dir = fullfile( datapath, floatid ); 
%     % find all the files with isus
%     isuses = findextension( data_target_dir, '.isus' );
%     clear isus
%     % cycle through the filenames
%     for ii = 1:length( isuses )
%         data_target = fullfile( data_target_dir, isuses(ii).name );
%         % load the data from message file
%         [isus(ii)] = parse_NO3msg( data_target);
%     end
%     save(cacheName, 'isus')
% else
%     load(cacheName)
% end

%% load the isus file from Navis float;
% define the folder
data_target_dir = fullfile( datapath, floatid ); 
% find all the files with isus
isuses = findextension( data_target_dir, '.isus' );
clear isus
% cycle through the filenames
for ii = 1:length( isuses )
    data_target = fullfile( data_target_dir, isuses(ii).name );
    % load the data from message file
    [isus(ii)] = parse_NO3msg( data_target);
end


%% Load the SUNA calibration file
% grab the SUNA serial number   
sunaSN = regexp(isus(1).SN, '(?<=:)\d*', 'match');
sunaSN = sunaSN{1};
fprintf(1, ' SUNA s/n: %s\n', sunaSN); 
% define the folder
cal_target_dir = [calpath, '/SUNA', sunaSN];
% define the calibration file
cal_target_name = strtrim(isus(1).CalibrationFile); 
% full name length
cal_target = fullfile( cal_target_dir, calVariant, cal_target_name); % there is not a consistent folder where the file is held, so you may have to adjust the folder name inside the folder cal_target_name 
% load the calibration data
cal = parseNO3cal(cal_target);

%% Choose processing option
cal.ProcessingOption = 'FW';           
% cal.ProcessingOption = 'NoTS';           
% cal.ProcessingOption = 'ArgoV1.1'    
% cal.ProcessingOption = 'ArgoV1.2';
fprintf(1, 'Processing Option: %s\n', cal.ProcessingOption)

%% Processing options that best reproduce implementations
switch (cal.ProcessingOption)
    case 'NoTS'
        % Current SUNA V2 Firmware (non-APF mode)
        cal.WL_offset = 210.0; % [nm] (unused_
        cal.pres_coef = 0.026; % (unused)
        Pcorr_flag = false;  
        cal.TSalgorithm = 'None';
        cal.Nitrate_Sensor_Offset = 0.0;  % No pressure offset in FW
    
    case 'FW'
        % Current SUNA V2 Firmware (damn close)
        cal.WL_offset = 210.0; % [nm] default
        cal.pres_coef = 0.026; % (unused)
        Pcorr_flag = false;
        cal.TSalgorithm = 'TCSS_Sakamoto2009';
        cal.Nitrate_Sensor_Offset = 0.0;  % No pressure offset in FW

    case 'ArgoV1.1'
        % ARGO Processing Bio-Argo concentration at the DAC Level V1.1
        %
        % "The operator may choose to treat OPTICAL_WAVELENGTH_OFFSET as a
        % tunable parameter for each float to eliminate these biases in near
        % surface nitrate. Thus, the value of OPTICAL_WAVELENGTH_OFFSET should
        % be reported."
        cal.WL_offset = 208.5; % [nm] de Fommervault et al. 2015
        cal.pres_coef = 0.026; % de Fommervault et al. 2015
        Pcorr_flag = true;
        cal.TSalgorithm = 'TCSS_Sakamoto2009';
        cal.Nitrate_Sensor_Offset = 0.8;  %  [m] PROVOR=1.5, NAVIS=0.8, SOLO=1.1

    case 'ArgoV1.2'
        % ARGO Processing Bio-Argo concentration at the DAC Level V1.2
        cal.WL_offset = 210.0;  % [nm]  Plant et al. 2023 (submitted)
        cal.pres_coef = 0.0265; % de Fommervault et al. 2015 (increased resolution)
        Pcorr_flag = true;
        cal.TSalgorithm = 'TCSS_Plant2023';
        cal.Nitrate_Sensor_Offset = 0.8;  %  [m] PROVOR=1.5, NAVIS=0.8, SOLO=1.1
end

%% cycle through each profile and calculate NO3, RMS error, etc.

Nprofiles = length(isus)
for ii = 1:Nprofiles
    disp(isus(ii).MessageFile)

    % calculate nitrate
    [NO3] = calc_FLOAT_NO3(isus(ii), cal, Pcorr_flag);
    isus(ii).Nmolar = NO3(:,1);
    isus(ii).NO3 = NO3;
end


%% Plot last two profiles

if (doProfilePlots)
    % figure(2);  cla
    
    nplots = 5;
    
    for ii = Nprofiles-1:Nprofiles
    
        figure;

        % Temp
        subplot( 1,nplots,1);
        plot( isus(ii).T, isus(ii).P, '.-'); axis ij
        title( ['Navis float ', floatid, ' profile ', num2str(ii), ' ', isus(ii).CurrentTimeUTC]);
        ylim( [0, 2000]);
        
        % Salinity
        subplot( 1,nplots,2);
        plot( isus(ii).S, isus(ii).P, '.-'); axis ij;
        ylim( [0, 2000]);
            
        % Firmware and MATLAB-estimated nitrate
        subplot( 1,nplots,3);
        plot( isus(ii).NO3fw, isus(ii).P, '-', isus(ii).Nmolar, isus(ii).P);
        axis ij;
        xlim( [0, 60]);
        ylim( [0, 2000])
        legend( 'SUNA FW', regexprep(cal.TSalgorithm, '_', '\\_'));
        title('NO3-');  
        
        % RMS Error (FW and MATLAB fit)
        subplot( 1,nplots,4 )
        plot( isus(ii).AbsorbanceFitResiduals, isus(ii).P, isus(ii).NO3(:,4), isus(ii).P, '.-')
        axis ij;
        % xlim( [0, 1e-3]);
        ylim( [0, 2000]);
        title('RMS Error');
        
        % subplot( 1,nplots,4 )
        % plot( isus(ii).SupplyVoltage, isus(ii).CTD.prSM, '.-')
        % axis ij;
        % xlim( [9.2, 10.2]);
        % ylim( [0, 2000]);
        % title('SupplyVoltage');
        
        % subplot( 1,nplots,5 )
        % plot( isus(ii).InternalRelativeHumidity, isus(ii).CTD.prSM, '.-');
        % axis ij;
        % xlim( [0.2, 1.0]);
        % ylim( [0, 2000]);
        % title('InternalRelativeHumidity');
    
        % Potential density anomoly    
        SA = gsw_SA_from_SP(isus(ii).S, isus(ii).P, -149.9926, -37.9816);
        CT = gsw_CT_from_t(SA, isus(ii).T, isus(ii).P);
        sigma0 = gsw_rho(SA, CT, isus(ii).P) - 1000.0;
    
        subplot( 1,nplots,5 )
        plot( sigma0, isus(ii).P, '.-');
        axis ij;
        ylim( [0, 2000]);
        title('\sigma_\theta');
    
        % pause
    end
end

%% Plot raw spectra vs. depth in 3D
if (doRawSpectralPlots)
    for ii = Nprofiles-1:Nprofiles
      dat = isus(ii);
      figure(ii*10);
      plotSpectrum3D(dat, cal);
      title( ['Navis float ', floatid, ' profile ', num2str(ii), ' ', isus(ii).CurrentTimeUTC]);
    end
end

