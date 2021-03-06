%Author: Torfinn Johnsrud
%Created: Fall 2017  UPDATED FOR SPECIFIC ANALYSIS 8/23/2018 by HLH - see
%Torfinn backop code for original file
% This program blends the Scale_height verticalwinds.m with the orginal
% TIEGCM_contour.m program.

% Version 2 - HLH
%   instead of having fixed z0 value, we make this an array of bottom
%   intergration values to understand how the pure diffusive profile
%   changes depending on where you set the bottom value.

% Version 3 - HLH, 1/17/19
%   Their is an issue with comparing the N2 model output at 400 km with the
%   integrated value at 400 km because of interpolation. The issue arises
%   because the pressure levels are not exactly at 400 km and when you
%   interpolated the values get screwed. This is seen more drastically with 
%   N2 because it has a much smaller scale height at these altitudes as
%   compared to HE. Therefore, everything is outputed both using pressure
%   coordinates and geometric interpolated points. For N2, the interpolated
%   output is not benefical to research. THINGS ARE IN GEOMETRIC HEIGHT.


% Version 4 - HLH 1/21/19
%   Do in pressure levels 


clear all;
close all;
clc;

%----------------
ut_want = 1;    %
alt_want = 400;
alpha = 0;              % thermal diffusion factor
alpha = -0.38; 
pdrag = 1;
species = 'HE';         % which species do you want? 'HE', 'N2', 'O', or 'O2'


aa1 = '~/Documents/MATLAB/TIEGCM/TIEGCM_output/';

% where to output things on pressure alt. levels
aa2_press = ['~/Documents/MATLAB/TIEGCM/Contour_textfiles/press_coord/', species, '_alpha_', num2str(alpha), '/'];

% where to output things that are on geometric alt. levels (interpolated)
aa2_geom = ['~/Documents/MATLAB/TIEGCM/Contour_textfiles/geom_coord/', species, '_alpha_', num2str(alpha), '/'];
linear=1;

%% -----Loading Viki's tiegcm simulation-----
% Follows (lon,lat,ilev,UT) format
if pdrag == 1
    filename = [aa1, 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc'];
    id = 'pdrag';
end
if pdrag == 0
    filename = [aa1, 'HSUVW.tiegcm2.0_dres.nodragtest_ctrSS_f107_180_001.nc'];
    id = 'ctrSS';
end

den = ncread(filename,'DEN');           % [g/cm^3]
zg = ncread(filename,'ZG')/1e5;         % geometric height [km]

he = ncread(filename,'HE');
n2 = ncread(filename,'N2');
o1 = ncread(filename,'O1');
o2 = ncread(filename,'O2');
tn = ncread(filename,'TN');
mbar = 1./(he/4+n2/28+o1/16+o2/32); %Getting mean molecular mass from mmr
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
wn = ncread(filename,'WN');

% fixed UT time 
den=squeeze(den(:,:,:,ut_want+1));      % Total neutral density [g/cm^3]
zg_tp1=squeeze(zg(:,:,:,ut_want+1));    % Geometric height [km]
% z_tp1=squeeze(z(:,:,:,ut_want+1));    % Geopotential height
he=squeeze(he(:,:,:,ut_want+1));        % Helium mass mixing ratio
n2=squeeze(n2(:,:,:,ut_want+1));        % N2 mass mixing ratio
o1=squeeze(o1(:,:,:,ut_want+1));        % O mass mixing ratio
o2=squeeze(o2(:,:,:,ut_want+1));        % O2 mass mixing ratio
tn=squeeze(tn(:,:,:,ut_want+1));        % Neutral temperature [K]
mbar=squeeze(mbar(:,:,:,ut_want+1));    % Mean molecular mass [kg/kmol]

%%
% now each matrix above has a value for every lat, lon, and pressure level. We
% care about helium, temperature, and mean moleulcar mass. 
%-------------------------------------------------------------------------

r=6372;                                 % Earth Radius [km]
g_tp = 9.81 .* (r./(r+zg_tp1)).^2;      % Find g as altitude changes. g(90km) = 8.7 m/s^2 ??
k=1.38e-23;                             % Boltzman's Constant
atom_unit=1.67e-27;                     % [kg/unit]
Av = 6.022141*10^23;                    % [#/mol]
kmol=6.02214*10^26;
mmw_he=0.004;                           % Helium atomic mass [kg/mol]
mmw_N2=0.02801;                         % N2 molecular mass
mmw_O1=0.016;                           % O1 molecular mass
meanmass = mbar./1000;                  % mean molecular mass [kg/mol]

if strcmp(species, 'HE') 
    mass_mix = he;
    mmw = mmw_he;
elseif strcmp(species, 'N2')
    mass_mix = n2;
    mmw = mmw_N2;
elseif strcmp(species, 'O2')
    mass_mix = o2;
    mmw = 2*mmw_O1;
elseif strcmp(species, 'O')
    mass_mix = o1;
    mmw = mmw_O1;
else
    ERROR = 'BAD SPECIES. PLEASE REDEFINE.' 
end

% ---- Species Number Density  
nhe_bf = (mass_mix.*den).*Av./(mmw.*1000);      % actual species number density [#/cm^3]

%------ Species Mass Density
mhe = mass_mix.*den;                            % species Mass Density [g/cm^3] 


%   (I know everthing following this says helium... please ignore that. It is 
%   specific for the desired species)

%-----Finding Pressure Scale Heights for gas constituents from kT/mg-----
Hp_he=k.*tn./(mmw/Av.*g_tp)/1000;               %in Km
Hp_mean=k.*tn./(meanmass/Av.*g_tp)/1000;


%----- Finding temperature scale height for Helium -----------------------
%-------------------------------------------------------------------------
% Calculate using 3 Point Differentiation Technique

% 144 longitudes (360/2.5), 72 latitudes (180/2.5), 57 altitudes
lon_num = 144;
lat_num = 72;
alt_num = 57;
H_he_star = zeros(size(mhe));
H_tn = zeros(size(tn));                             % 1/temperature scale height?
nhe_diff = zeros(lon_num, lat_num, alt_num);        % number density of helium at z if pure diff. profile.

          
% ******* VALUES CALCULATED HERE ARE equal to 1/the  corrresponding scale heights!!!!!!

for i = 1:lon_num
    for j = 1:lat_num
        for k = 1:alt_num          
            if k == 1 % First Point   
                coeff1 = (2*zg_tp1(i, j, 1)-zg_tp1(i, j, 2)-zg_tp1(i, j, 3))/((zg_tp1(i, j, 1)-...
                    zg_tp1(i, j, 2))*(zg_tp1(i, j, 1)-zg_tp1(i, j, 3)));

                coeff2 = (2*zg_tp1(i, j, 1)-zg_tp1(i, j, 1)-zg_tp1(i, j, 3))/((zg_tp1(i, j, 2)-...
                    zg_tp1(i, j, 1))*(zg_tp1(i, j, 2)-zg_tp1(i, j, 3)));

                coeff3 = (2*zg_tp1(i, j, 1)-zg_tp1(i, j, 1)-zg_tp1(i, j, 2))/((zg_tp1(i, j, 3)-...
                    zg_tp1(i, j, 1))*(zg_tp1(i, j, 3)-zg_tp1(i, j, 2)));

                H_he_star(i,j, 1) = -1/mhe(i, j, 1)*(mhe(i, j, 1)*coeff1+mhe(i, j, 2)*coeff2+...
                    mhe(i, j, 3)*coeff3);
                H_tn(i, j, 1) = 1/tn(i, j, 1)*((tn(i, j, 1)*coeff1+tn(i, j, 2)*coeff2+...
                    tn(i, j, 3)*coeff3));


            elseif k == alt_num %Last point
                coeff1 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k-1)-zg_tp1(i, j, k))/((zg_tp1(i, j, k-2)-...
                    zg_tp1(i, j, k-1))*(zg_tp1(i, j, k-2)-zg_tp1(i, j, k)));

                coeff2 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k-2)-zg_tp1(i, j, k))/((zg_tp1(i, j, k-1)-...
                    zg_tp1(i, j, k-2))*(zg_tp1(i, j, k-1)-zg_tp1(i, j, k)));

                coeff3 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k-2)-zg_tp1(i, j, k-1))/((zg_tp1(i, j, k)-...
                    zg_tp1(i, j, k-2))*(zg_tp1(i, j, k)-zg_tp1(i, j, k-1)));

                H_he_star(i, j, k) = -1/mhe(i, j, k)*(mhe(i, j, k-2)*coeff1+mhe(i, j, k-1)*coeff2+...
                    mhe(i, j, k)*coeff3);
                H_tn(i, j, k) = 1/tn(i, j, k)*((tn(i, j, k-2)*coeff1+tn(i, j, k-1)*coeff2+...
                    tn(i, j, k)*coeff3));


            else %Middle Points
                coeff1 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k)-zg_tp1(i, j, k+1))/((zg_tp1(i, j, k-1)-...
                    zg_tp1(i, j, k))*(zg_tp1(i, j, k-1)-zg_tp1(i, j, k+1)));

                coeff2 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k-1)-zg_tp1(i, j, k+1))/((zg_tp1(i, j, k)-...
                    zg_tp1(i, j, k-1))*(zg_tp1(i, j, k)-zg_tp1(i, j, k+1)));

                coeff3 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k-1)-zg_tp1(i, j, k))/((zg_tp1(i, j, k+1)-...
                    zg_tp1(i, j, k-1))*(zg_tp1(i, j, k+1)-zg_tp1(i, j, k)));

                H_he_star(i, j, k) = -1/mhe(i, j, k)*(mhe(i, j, k-1)*coeff1+mhe(i, j, k)*coeff2+...
                    mhe(i, j, k+1)*coeff3);
                H_tn(i, j, k) = 1/tn(i, j, k)*((tn(i, j, k-1)*coeff1+tn(i, j, k)*coeff2+...
                    tn(i, j, k+1)*coeff3));

            end
        end
    end
end


% ****** Now the ACTUAL scale heights are found....

%-----Get Scale Height from Inverse-----
% The values below are the scale heights for every lat, long and altitude 

H_he_star = 1./H_he_star;                % actual scale height from model density fit

if strcmp(species, 'HE')
    H_temp = 1./(H_tn .* (1 + alpha));   % helium temperature scale height includes thermal diffusion 
    H_diff = 1./(1./H_temp + 1./Hp_he);  % helium diffusive profile (i.e. the diffusive density scale height)
else
    H_temp = 1./H_tn;                               % other species temperature scale height
    H_diff = 1./(1./H_temp + 1./Hp_he);           %  the diffusive density scale height)
end


%% -----Select Altitude for Interpolation and Diffusive Profile Integration -----
dz_geom = 10;
z0 = (110:dz_geom:alt_want);          % bottom integration values [km]
num_z0 = length(z0);

% use these if you want to stay in pressure coordinates for N2
nhe_model_z0_press = zeros(lon_num, lat_num);       % global model number density closest to z0 
nhe_model_400_press = zeros(lon_num, lat_num);
nhe_int_400_press = zeros(lon_num, lat_num);        % global integrated number density at ~400 km
nhe_perc_diff_400_press = zeros(lon_num, lat_num);  % global percent difference between model and integated density at ~400 km

% use these for interpolating values for HE
nhe_model_z0_geom = zeros(lon_num, lat_num);
nhe_model_400_geom = zeros(lon_num, lat_num);      % global model number density interpolated to 400 km
nhe_int_400_geom = zeros(lon_num, lat_num);        % global integrated number density interpoltated to 400 km for every z0
nhe_perc_diff_400_geom = zeros(lon_num, lat_num);  % global percent difference between model and integated density for every z0, intepolated to 400 km

H_diff_geom_z0 = zeros(lon_num, lat_num);          % global interpolated diffusive scale height to use for integration
geom_index_2 = zg_tp1-alt_want;

% Find pressure altitude index to closest to 400 km end point 
for i = 1:num_z0
    for m=1:lon_num
        for n=1:lat_num

            % ----------------------------------------------------------------
            %                   VALUES AT ~400 km
            % ----------------------------------------------------------------

            % Find pressure altitude index to closest to 400 km end point             
            geom_index_end=squeeze(geom_index_2(m,n,:));  
            [val_2, p_index_400] = min(abs(geom_index_end)); 

            % ---------------- PRESSURE COORDINATES ---------------------------
            % Model values at ~400km  (using pressure surfaces)      
            nhe_model_400_press(m, n) = nhe_bf(m, n, p_index_400);

            % ---------------- GEOMETRIC COORDINATES --------------------------
            % Interpolate the model values to 400 km
                if zg_tp1(m,n,p_index_400) >= alt_want      % Linearly interpolate backwards
                    x0 = zg_tp1(m,n,p_index_400-1);
                    x1 = zg_tp1(m,n,p_index_400);

                    y0 = nhe_bf(m,n,p_index_400-1);
                    y1 = nhe_bf(m,n,p_index_400);  
                    nhe_model_400_geom(m, n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);  

                elseif zg_tp1(m,n,p_index_400) < alt_want     % Linearly interpolate forwards
                    x0 = zg_tp1(m,n,p_index_400);
                    x1 = zg_tp1(m,n,p_index_400+1);

                    y0 = nhe_bf(m,n,p_index_400);
                    y1 = nhe_bf(m,n,p_index_400+1);
                    nhe_model_400_geom(m, n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);                   
                    
                end

            % ----------------------------------------------------------------
            %                   VALUES AT z0 km
            % ----------------------------------------------------------------


            % ---------------- PRESSURE COORDINATES ---------------------------        
            % Find pressure altitude index to closest to z0 starting integration point

                geom_index_1 = zg_tp1-z0(i);            
                geom_index_start=squeeze(geom_index_1(m,n,:));
                [val_1, p_index_z0] = min(abs(geom_index_start));

                % Model Values at ~z0 km  (using pressure surfaces)
                nhe_model_z0_press(m, n) = nhe_bf(m, n, p_index_z0);

                % Integrate the model density at z0 up to the 400 km
                sum = 0;
                for j = (p_index_z0) : (p_index_400-1)

                    dz = zg_tp1(m, n, j + 1) - zg_tp1(m, n, j);
                    H_avg = (H_diff(m, n, j + 1) + H_diff(m, n, j)) / 2;
                    sum = sum + dz/H_avg; 
                end

                % Integrated Diffusive Density at ~400 km 
                nhe_int_400_press(m, n) = nhe_model_z0_press(m, n) * exp(-sum); 

                % Percent difference at ~400 km           
                nhe_perc_diff_400_press(m, n) = ((nhe_bf(m, n, p_index_400) ./ nhe_int_400_press(m, n) ) - 1) .* 100;           

        end     % for 1:lat_num
    end         % for 1:lon_num
    
        
    % ---------------------------------------------------------------------
    %                      OUTPUTNG VALUES 
    % ---------------------------------------------------------------------

    % ----------------- PRESSURE COORDINATES ------------------------

    % make folder for specific z0 starting altitude
    new_fold = ['z0_alt_', num2str(z0(i)), '/'];
    status = mkdir([aa2_press, new_fold]);
    fprintf('Outputting %s\n', aa2_press + new_fold);      

    % find percent difference of model vs. integrated densities @ every lat/lon
    % turns lon/lat ---> lat/lon
    dlmwrite([aa2_press, new_fold, species, '_dens_RealVsDiff_z0_', num2str(z0(i)), '_', id, '.txt'], nhe_perc_diff_400_press.');

    % to look at the calculated diffusive density
    dlmwrite([aa2_press, new_fold, species, '_dens_Calculated_Diff_z0_' num2str(z0(i)), '_', id, '.txt'], nhe_int_400_press.');

    % to look at the real model density at z0 
    dlmwrite([aa2_press, new_fold, species, '_dens_model_', num2str(z0(i)), 'km_', id, '.txt'], nhe_model_z0_press.');  
    
end   % for 1:num_z0



%% to look at the model density at 400 km (either ~400km or interpolated)
% ----------------- PRESSURE COORDINATES ------------------------
dlmwrite([aa2_press, species, '_dens_model_400km_', id, '.txt'], nhe_model_400_press.');

% ----------------- GEOMETRIC COORDINATES ------------------------
dlmwrite([aa2_geom, species, '_dens_model_400km_', id, '.txt'], nhe_model_400_geom.');

% output z0 altitude points
dlmwrite([aa2_press 'z0_altitude_values.txt'], z0)
dlmwrite([aa2_geom 'z0_altitude_values.txt'], z0)


fprintf('\nDONE!\n')



