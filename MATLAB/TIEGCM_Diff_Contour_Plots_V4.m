%Author: Hannah Holt (from Torfinn Johnsrud)
%Created: Fall 2017  UPDATED FOR SPECIFIC ANALYSIS 8/23/2018 by HLH - see
%Torfinn backup code for original file
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
%   output is not benefical to research.


% Version 4 - HLH 1/21/19
%   Do in pressure levels 

% OUTPUT CAN BE USED WITH THE FOLLOWING PYTHON FILES
%   1) TIEGCM_contour_percDiff_V2.py
%   2) TIEGCM_10perc_Change_V2.py
%   3) TIEGCM_10perc_Change_thermONLY_V2.py
%   4) TIEGCM_contour_percDiff_therm_V0.py
%   5) ZpressVZpot.py

clear all;
close all;
clc;

%----------------
ut_want = 1;    %

z0_low = 110E3;         % lowest value we want to integrate from [m]
z0_high = 400E3;        % highest value we want to integrate too ~ 400 m

alpha = 0;              % thermal diffusion factor
alpha = -0.38; 
pdrag = 1;
species = 'HE';         % which species do you want? 'HE', 'N2', 'O', or 'O2'
%----------------

aa1 = '~/Documents/MATLAB/TIEGCM/TIEGCM_output/';

% where to output things on pressure levels
aa2_press = ['~/Documents/MATLAB/TIEGCM/Contour_textfiles/pressure/', species, '_alpha_', num2str(alpha), '/'];

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

den = ncread(filename,'DEN');               % [g/cm^3] calculated on ilev
zp = ncread(filename, 'Z')/100;             % geopotential height [m] on ilev
p0 = ncread(filename, 'p0_model');          % p0 reference pressure used by model [microbars]
g0 = ncread(filename, 'grav')/100;          % const. gravitational acceleration [m/s]

plvl = ncread(filename, 'lev');             % midpoint pressure level
iplvl = ncread(filename, 'ilev');           % interface pressure level

n2 = ncread(filename,'N2');                 % on lev
he = ncread(filename,'HE');                 % on lev
tn = ncread(filename,'TN');                 % on lev
o1 = ncread(filename,'O1');                 % on lev               
o2 = ncread(filename,'O2');                 % on lev
mbar = 1./(he/4 + n2/28 + o1/16 + o2/32);   % mean molecular mass from mmr
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');

% condense to fixed UT time 
den = squeeze(den(:,:,:,ut_want+1));      % Total neutral density [g/cm^3] on ilev
zp = squeeze(zp(:,:,:,ut_want+1));        % Geopotential height [m], on ilev
he = squeeze(he(:,:,:,ut_want+1));        % Helium mass mixing ratio
n2 = squeeze(n2(:,:,:,ut_want+1));        % N2 mass mixing ratio
tn = squeeze(tn(:,:,:,ut_want+1));        % Neutral temperature [K] on lev
mbar = squeeze(mbar(:,:,:,ut_want+1));    % Mean molecular mass [g/mol]

%o1=squeeze(o1(:,:,:,ut_want+1));        % O mass mixing ratio
%o2=squeeze(o2(:,:,:,ut_want+1));        % O2 mass mixing ratio


%%
% now each matrix above has a value for every lat, lon, and pressure level. We
% care about helium, temperature, and mean moleulcar mass. 
%-------------------------------------------------------------------------
R = 8.314;                              % R = [J/mol K]
Av = 6.022141*10^23;                    % [#/mol]

mmw_he = 4;                             % Helium atomic mass [g/mol]
mmw_N2 = 28.01;                         % N2 molecular mass  [g/mol]
mmw_O1 = 16.;                           % O1 molecular mass  [g/mol]

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
n_dens = mass_mix .* den ./ mmw .* Av ;      % actual species number density [#/cm^3]

%------ Species Mass Density
m_dens = mass_mix .* den;                    % species Mass Density [g/cm^3] 

%-----Finding Pressure Scale Heights for gas constituents from RT/M g0 -----

Hp = R .* tn ./ (g0 .* mmw ./ 1000);            % species pressure scale height in m              
Hp_mean = R .* tn ./ ( g0 .* mbar ./ 1000);     % mean scale height [m]


%% ----- Finding temperature and pressure scale height  -------------------
% -------------------------------------------------------------------------
% Calculate using 3 Point Differentiation Technique

% 144 longitudes (360/2.5), 72 latitudes (180/2.5), 57 pressure altitudes
lon_num = 144;
lat_num = 72;
alt_num = 57;                       % number of pressure level coordinates

% p index closest to 110 km generally
[M, I] = min(abs(zp - z0_low), [], 3);
p_indx_low = mode(I(:));

% p index cloest to 400 km generally
[M, I] = min(abs(zp - z0_high), [], 3);
p_indx_high = mode(I(:));

dzp = zeros(lon_num, lat_num, alt_num-1);              % find the geopotential altitude differences btw pressure levels (different values for each lat/lon)                                                                                                           
H_tn = zeros(lon_num, lat_num, alt_num);               % 1/temperature scale height
n_diffuse = zeros(lon_num, lat_num, alt_num);          % number density of helium at z if pure diffusive profile.

% ******* VALUES CALCULATED HERE ARE equal to 1/the  corrresponding scale heights!!!!!!
% find 1 / temperature scale height by 3-pt differentiation method

% NOTE - This gives bad values for top 2 and bottom 2 model pressure levels
for k = 1:alt_num 
    
    % first, find global p index closest to ~110 km.
    
    if k == 1 % First Point   
        coeff1 = (2.*zp(:, :, 1)-zp(:, :, 2)-zp(:, :, 3))./((zp(:, :, 1)-...
            zp(:, :, 2)).*(zp(:, :, 1)-zp(:, :, 3)));

        coeff2 = (2.*zp(:, :, 1)-zp(:, :, 1)-zp(:, :, 3))./((zp(:, :, 2)-...
            zp(:, :, 1)).*(zp(:, :, 2)-zp(:, :, 3)));

        coeff3 = (2.*zp(:, :, 1)-zp(:, :, 1)-zp(:, :, 2))./((zp(:, :, 3)-...
            zp(:, :, 1)).*(zp(:, :, 3)-zp(:, :, 2)));

        H_tn(:, :, 1) = 1./tn(:, :, 1) .* ((tn(:, :, 1).*coeff1+tn(:, :, 2).*coeff2+...
            tn(:, :, 3).*coeff3));
        
        dzp(:,:,k) = zp(:,:,k+1) - zp(:,:,k); 

        
    elseif k == alt_num %Last point
        coeff1 = (2.*zp(:, :, k)-zp(:, :, k-1)-zp(:, :, k))./((zp(:, :, k-2)-...
            zp(:, :, k-1)).*(zp(:, :, k-2)-zp(:, :, k)));

        coeff2 = (2.*zp(:, :, k)-zp(:, :, k-2)-zp(:, :, k))./((zp(:, :, k-1)-...
            zp(:, :, k-2)).*(zp(:, :, k-1)-zp(:, :, k)));

        coeff3 = (2.*zp(:, :, k)-zp(:, :, k-2)-zp(:, :, k-1))./((zp(:, :, k)-...
            zp(:, :, k-2)).*(zp(:, :, k)-zp(:, :, k-1)));
        
        H_tn(:, :, k) = 1./tn(:, :, k) .* ((tn(:, :, k-2).*coeff1+tn(:, :, k-1).*coeff2+...
            tn(:, :, k).*coeff3));

    else %Middle Points
        coeff1 = (2.*zp(:, :, k)-zp(:, :, k)-zp(:, :, k+1))./((zp(:, :, k-1)-...
            zp(:, :, k)).*(zp(:, :, k-1)-zp(:, :, k+1)));

        coeff2 = (2.*zp(:, :, k)-zp(:, :, k-1)-zp(:, :, k+1))./((zp(:, :, k)-...
            zp(:, :, k-1)).*(zp(:, :, k)-zp(:, :, k+1)));

        coeff3 = (2.*zp(:, :, k)-zp(:, :, k-1)-zp(:, :, k))./((zp(:, :, k+1)-...
            zp(:, :, k-1)).*(zp(:, :, k+1)-zp(:, :, k)));
                   
        H_tn(:, :, k) = 1./tn(:, :, k) .* ((tn(:, :, k-1).*coeff1+tn(:, :, k).*coeff2+...
            tn(:, :, k+1).*coeff3));
        
        dzp(:,:,k) = zp(:,:,k+1) - zp(:,:,k);
         
    end
    
end


% ****** Now the ACTUAL scale heights are found....

%-----Get Scale Height from Inverse-----
% The values below are the scale heights for every lat, long and altitude 

if strcmp(species, 'HE')
    H_temp = 1./(H_tn .* (1 + alpha));      % temperature scale height includes thermal diffusion [m]
    H_rho = 1./H_temp + 1./Hp;
    H_diffuse = 1./(1./H_temp + 1./Hp);     % diffusive density scale height
else
    H_temp = 1./H_tn;                       % other species temperature scale height (no thermal diffusion)
    H_rho = 1./H_temp + 1./Hp;              % 1/H_rho
    H_diffuse = 1./(1./H_temp + 1./Hp);     % the diffusive density scale height
end


%% ------------------ Diffusive Profile Integration -----------------------

% at this point we have the scale height for the species at every pressure
% level (technically every geopotential coordinate that equals the pressure lvl) 
% and we want to integrate this scale height up to the constant global pressure
% level closest to our desired point.

num_p = p_indx_high - p_indx_low;           % go from p index low to p_index_high includeing the top value (42-10) + 1

n_dens_model_z0 = zeros(lon_num, lat_num);
n_dens_diffuse_400 = zeros(lon_num, lat_num);

% get model output at p_indx_high
n_dens_model_400 = n_dens(:, :, p_indx_high);               % number density [#/cm^3]
%%
for k = 0:num_p % 0 to 32 (i.e. 32 values)
 
    % get model output at the z0 pressure level in order to integrate (we want to start at p_indx_low)
    n_dens_model_z0 = n_dens(:, :, p_indx_low + k);     % number density [#/cm^3]
 
    % now integrate from the z0 value to p_indx_high   
    sum = zeros(lon_num, lat_num);
    for i = (p_indx_low + k):p_indx_high-1 
        sum = sum + dzp(:,:,i)./H_diffuse(:,:,i);          
    end

    n_dens_diffuse_400 = n_dens(:, :, p_indx_low + k) .* exp(-sum);
    n_dens_perc_diff_400 = ( (n_dens_model_400 ./  n_dens_diffuse_400) - 1) .* 100;

    % --------------------- IMPORTANT NOTE -----------------------------------
    % >>> n_dens_diffuse_400(:, :, 1) = the density at ~400km if we integrated from
    % z0 = p_indx_low to p_indx_high

    % >>> n_dens_diffuse_400(:, :, 33) = the density at ~400 km if we integrated
    % from z0 = p_indx_high to p_indx_high ---> i.e. NOTHING HAPPENED. The values are
    % exactly the same as the model output at the 400 km p lvl. Really only
    % indices 1-32 have been integrated.


    % ------------------ Outputting Values to Text Files ----------------------

    % make folder for specific p level starting integration altitudes
    new_fold = ['plvl_', num2str(plvl(p_indx_low+k)), '/'];
    status = mkdir([aa2_press, new_fold]);
    fprintf('Outputting %s\n', aa2_press + new_fold);      
    
    % to look at the percent different btw integrated and model density at ~400km
    % turn lon/lat ---> lat/lon
    dlmwrite([aa2_press, new_fold, species, '_dens_perc_diff_400_', id, '.txt'], n_dens_perc_diff_400.');

    % to look at the calculated diffusive density at ~400 km
    dlmwrite([aa2_press, new_fold, species, '_dens_Diffusive_400_', id, '.txt'], n_dens_diffuse_400.');

    % to look at the model density at z0 
    dlmwrite([aa2_press, new_fold, species, '_dens_model_z0_', id, '.txt'], n_dens_model_z0.');  

    % to look at geopotential altitude values [km] for this specific z0 pressure level 
    dlmwrite([aa2_press, new_fold, 'zp_altitudes_km_', id, '.txt'], (zp(:,:,p_indx_low+k)./1000).');  
    
end

% the model density at ~400km
dlmwrite([aa2_press, species, '_dens_model_400_', id, '.txt'], n_dens_model_400.');

% output of the pressure levels 
dlmwrite([aa2_press 'p_lvls.txt'], plvl(p_indx_low:p_indx_high))

fprintf('\nDONE!\n')


%% ------ THE FOLLOWING IS FOR COMPARING THE SCALE HEIGHTS. ------
% 
% 
% x = lon;    % longitude values [deg]
% y = lat;    % latitude values [deg]
% [X, Y] = meshgrid(x, y);
% 
% start = p_indx_low;
% stop = p_indx_high;
% 
% lshift = -0.1;
% wshift = 0.07;
% 
% len = stop - start;
% A = 1./H_temp;
% B = 1./Hp;
% C = H_rho;
% perc_Ht = A./C .* 100;                  % percent of 1/H_diffusive that is 1/H_temp
% perc_Hp = B./C .* 100;
% num_cont = 25;
% 
% mag_eq = importdata('Magnetic_equator_lat_lon.txt'); % first column is lon, second column is lat
% 
% for i=start:stop   
%     thisfig = figure('units','normalized','outerposition',[1 1 1 .6]);
%     
%     c = jet(100);
%     colormap(c);  
%     
%     h = subplot(1,3,1);
% %     contourf(X, Y, (H_diffuse(:,:,i)./1000)', num_cont)
% %     title('Density Scale Height')    
% %     cbar.Label.String = 'H_{\rho} [km]';
% 
%     contourf(X, Y, (H_rho(:,:,i).*1000)', num_cont)
%     title('1/H_{\rho}')
%     cbar = colorbar();
%     cbar.Label.String = '1/H_{\rho} [1/km]';
%     
%     hold on;   
%     plot(mag_eq(:,1), mag_eq(:,2), 'r', 'Linewidth', 2.5)
%     
%     xlabel('LT [Hr]')
%     ylabel('Lat [Deg]')
%     xticks(linspace(-180, 180, 9))
%     xticklabels({'12', '15', '18', '21', '0', '3', '6', '9', '12'})
%     
%     p = get(h, 'pos');
%     p(1) = p(1) + lshift;       % adjust left starting point      
%     p(3) = p(3) + wshift;     % adjust width
%     set(h, 'pos', p);
%        
%     h = subplot(1,3,2);
% %     contourf(X, Y, (H_temp(:,:,i)./1000)', num_cont)
% %     title('Temperature Scale Height')
% %     cbar.Label.String = 'H_T';
%  
%     contourf(X, Y, perc_Ht(:,:,i)', num_cont)
%     caxis([0 100])
%     title('1/H_T percentage of 1/H_{\rho}')
%     cbar = colorbar();
%     cbar.Label.String = 'Percent 1/H_T';
%     
%     hold on; 
%     plot(mag_eq(:,1), mag_eq(:,2), 'r', 'Linewidth', 2.5)  
%     xlabel('LT [Hr]')
%     ylabel('Lat [Deg]')
%     
%     xticks(linspace(-180, 180, 9))
%     xticklabels({'12', '15', '18', '21', '0', '3', '6', '9', '12'})
%     
%     
%     p = get(h, 'pos');
%     p(1) = p(1) + lshift + .05;      % adjust left starting point      
%     p(3) = p(3) + wshift;            % adjust width
%     set(h, 'pos', p);
%     
%     h = subplot(1,3,3);
% %     contourf(X, Y, (Hp(:,:,i)./1000)', num_cont)       
% %     cbar.Label.String = 'H_p [km]';
% %     title('Pressure Scale Height')
% %     
%     contourf(X, Y, perc_Hp(:,:,i)', num_cont)
%     caxis([0 100]); 
%     title('1/H_P percentage of 1/H_{\rho}')
%     cbar = colorbar(); 
%     cbar.Label.String = 'Percent 1/H_P';
%     
%     hold on;   
%     plot(mag_eq(:,1), mag_eq(:,2), 'r', 'Linewidth', 2.5)   
%     xlabel('LT [Hr]')
%     ylabel('Lat [Deg]')
%    
%     xticks(linspace(-180, 180, 9))
%     xticklabels({'12', '15', '18', '21', '0', '3', '6', '9', '12'})
%     
%     
%     p = get(h, 'pos');
%     p(1) = p(1) + lshift + .1;       % adjust left starting point      
%     p(3) = p(3) + wshift;     % adjust width
%     set(h, 'pos', p);      
%     
%     t1 = ['TIEGCM ', species, ', UT = 0, plvl = ', num2str(plvl(i))];
%     t2 = ['Approximate Height [km] = ', num2str(zp(1,1,i)/1000)];
%     sgtitle({t1; t2; ' '})
%     
% %     if i > 9    
% %         name = ['~/Desktop/', species, '/ScaleHeightComp_', num2str(i), '.png'];
% %     else
% %         name = ['~/Desktop/', species, '/ScaleHeightComp_0', num2str(i), '.png'];
% %     end
% 
%     if i > 9    
%         name = ['~/Desktop/', species, '/ProfileComp_', num2str(i), '.png'];
%     else
%         name = ['~/Desktop/', species, '/ProfileComp_0', num2str(i), '.png'];
%     end
% 
%         
%     saveas(thisfig, name)    
% end
% 
% 

%% ---- The Following is to find the actual pressure value [mb] for each pressure level for a given lat/lon

feat = 1;

% ----- Global Features -------
if feat == 1
lon_want = 102.5;       % North He nighttime maximum feature
lat_want = 63.75;       
end
if feat == 2
lon_want = 40;          % South He nighttime maximum feature
lat_want = 58.75;  
end
if feat == 3
lat_want = 58.75;       % Intermediate He density area
lon_want = 40;      
end
if feat == 4 
lon_want = -82.5;       % North Daytime He minimum feature
lat_want = 21.25;
end
if feat == 5
lon_want = -82.5;       % South Daytime He minimum feature
lat_want = -53.75;  
end

% condense geopotential height to fixed lat and long
lat_ind = find(lat == lat_want);
lon_ind = find(lon == lon_want);

% condense geopotential height and total mass densitity
% geopotential height of every pressure level at specific lat/lon pt [m]
zp_cond = squeeze(zp(lon_ind, lat_ind, :));     
den_cond = squeeze(den(lon_ind, lat_ind, :)) * 1000;   % total mass density condensed [ kg/m^3]

%
p0_Pa = p0 / 10;  % TIEGCM gives pressure in terms of microbars, convert to Pascals
p_Pa = p0_Pa .* exp(-iplvl);

% now plot geopotential altitude versus height for specific feat

thisfig = figure();

semilogx(p_Pa, zp_cond/1000);
xlabel('Pressure [Pa]');
ylabel('Geopotential Altitude [km]');
ylim([95 450]);
grid on;
title(['TIEGCM Pressure for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)])
% saveas(thisfig, ['./Figures/TIEGCM_press_latlon_', num2str(lat_want), '_', num2str(lon_want), '.pdf'])
close(thisfig);

% Now calculate the dp/dz gradient using Torfinn's threepointgradient.m
% function

linear = 0;
dpdz = threepointgradient(zp_cond, p_Pa, linear);
zdenom = (-den_cond.*g0)';
perc_dif = (dpdz ./ (-den_cond.*g0)' - 1 ) .* 100;


this_fig = figure('units','normalized','outerposition',[1 1 1 .6]);
subplot(121);
semilogx(dpdz, zp_cond/1000);
hold on;
semilogx(-den_cond*g0, zp_cond/1000);
legend('dp/dz', '-\rho * g0', 'location', 'best', 'FontSize', 16);
grid on;
xlabel('Value [Pa/m]')
ylabel('Geopotential Altitude [km]')
ylim([95 450]);

subplot(122)
plot(perc_dif, zp_cond/1000)
grid on;
ylabel('Geopotential Altitude [km]')
xlabel('Percent Difference')
ylim([95 450]);
sgtitle(['TIEGCM Pressure Comparison (ilev) for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)], 'FontWeight', 'bold');
% saveas(this_fig, ['./Figures/TIEGCM_hydrostatic_ilev_latlon_', num2str(lat_want), '_', num2str(lon_want), '.png'])
close(this_fig)



%% COMPARING THE MIDPOINTS TO INTERFACES - 7/16/19

% first, we want to put our Hp_mean (which is currently on midpoint lvls)
% onto interface levels.
Hp_lev = squeeze(Hp_mean(lon_ind, lat_ind, :))/1000;     % [km]
Hp_ilev = zeros(length(Hp_lev), 1, 'single');       % [km]

% find average between elements
for i = 1:length(Hp_lev)
    if i ~= 1
        Hp_ilev(i) = 0.5*(Hp_lev(i) + Hp_lev(i-1));
    else
        Hp_ilev(i) = Hp_lev(i); 
    end
end

% now we want to compute the Hp_ilev using the TIEGCM model documentation way,
% i.e. d(Zp)/d(iplvl)
linear = 0;
Hp_ilev_3pt = (threepointgradient(iplvl, zp_cond, 0)/1000)';  % [km]

dz = iplvl(2) - iplvl(1);       % pressure level spacing
Hp_ilev_cendiff = gradient(zp_cond, dz)/1000;              % [km]

% compare the values of Hp using different method
thisfig = figure('units','normalized','outerposition',[1 1 1 .6]);
subplot(121);
plot(Hp_lev, zp_cond/1000, 'Linewidth', 2, 'Color', '#77AC30')
hold on;
plot(Hp_ilev, zp_cond/1000, 'b', 'Linewidth', 2);
hold on;
plot(Hp_ilev_3pt, zp_cond/1000, '.r', 'Markersize', 10)
hold on;
plot(Hp_ilev_cendiff, zp_cond/1000, '--', 'Color', '#EDB120', 'Linewidth', 2)
ylim([100 400])
lgd = legend('H_p^{lev}', 'H_p^{ilev}, Method 1: Hp_{mid} raw avg', 'H_p^{ilev}, Method 2: dZ/dz (3pt)', 'H_p^{ilev}, Method 3: dZ/dz (Cent. Diff.)');
lgd.Location = 'best';
lgd.FontSize = 12;
grid on;
title('TIEGCM Pressure Scale Height Comparison')
xlabel('Mean Pressure Scale Height, Hp [km]');
ylabel('Geopotential Altitude [km]');

% zoomed plot
subplot(122);
plot(Hp_lev, zp_cond/1000, 'Linewidth', 2, 'Color', '#77AC30')
hold on;
plot(Hp_ilev, zp_cond/1000, 'b', 'Linewidth', 2 );
hold on;
plot(Hp_ilev_3pt, zp_cond/1000, '-r', 'Linewidth', 2)
hold on;
plot(Hp_ilev_cendiff, zp_cond/1000, '--', 'Color', '#EDB120', 'Linewidth', 2)
ylim([300 305])
lgd = legend('H_p^{lev}', 'H_p^{ilev}, Method 1: Hp_{mid} raw avg', 'H_p^{ilev}, Method 2: dZ/dz (3pt)', 'H_p^{ilev}, Method 3: dZ/dz (Cent. Diff.)');
lgd.Location = 'best';
lgd.FontSize = 12;
grid on;
title('TIEGCM Pressure Scale Height Comparison ZOOMED')
xlabel('Mean Pressure Scale Height, Hp [km]');
ylabel('Geopotential Altitude [km]');
 



















