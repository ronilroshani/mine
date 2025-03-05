
%%
global p_path N_GPS_path N_GLONASS_path N_BEIDOU_path I_path S_path F_path N_MIX_path
p_path = [pwd '\'];                                       % Program path
% R_path = [p_path 'RINEX\' num2str(YEAR) '\' DOY '\'];     % RINEX path
N_GPS_path = [p_path 'NAV\GPS'];                          % Navigation path GPS
N_GLONASS_path = [p_path 'NAV\GLONASS'];                  % Navigation path GLONASS
N_BEIDOU_path = [p_path 'NAV\BEIDOU'];                    % Navigation path BEIDOU
N_MIX_path = [p_path 'NAV'];% num2str(YEAR) '\' DOY '\'];     % Navigation path BEIDOU

I_path = [p_path 'IONEX\'];                               % IONEX path
% S_path = [p_path 'Results\' num2str(YEAR) '\' DOY '\'];   % Save Results
F_path = [p_path 'Rinex_read_functions\'];

%% frequencies
global  c f1G f2G f1C f2C order processing_interval Re h Cutoff A PG PC PR fig
global GPS_flag GLO_flag GAL_flag BDS_flag QZS_flag SBS_flag


c=299792458;            %speed of light
%--number of ionosphere model parameter groups
fig=12;

% weight of each system observation
PG=1; PC=1; PR=1;

%frequencies GPS [MHz]
f1G = 1575.42e6;   %hz      
f2G = 1227.60e6;   

%frequencies BeiDou [MHz]
f1C  = 1561.098e6;       
f2C  = 1207.140e6;     

%% % Compute STEC & VTEC
Re = 6371e3;                                     % Radius of Earth
h = 450e3;                                      % Raidus of Single Layer
A=40.3;
order=0;                                        %order
processing_interval = 30;

%% 
% load ('rcf');
Cutoff =[cut(1);cut(2);cut(3)]*pi/180;                      % GRC Cutoff Angle 
GPS_flag = flag(1); GLO_flag = flag(2); GAL_flag = 0; BDS_flag = flag(3); QZS_flag = 0;SBS_flag =0;

