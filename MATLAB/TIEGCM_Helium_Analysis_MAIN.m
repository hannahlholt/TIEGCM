
%Author: Torfinn Johnsrud
% Created: Feb. 6 2019

% Version 2 - Updated 7/16/19 by Hannah Holt. All the altitudes are
% switched from geometric to geopotential. 

addpath('./Utilities');

clear all;
close all;

aa1 = '~/Documents/MATLAB/TIEGCM/TIEGCM_output/';

%----------------
ut_want = 1;        % what time segment desired from simulation
feat = 1;           % Select latitude and longitude desired
pdrag = 1;          % 1 if pdrag file used, 0 if not

% ----- Global Features -------
if feat == 1
lon_want = 95;      % North nighttime maximum feature
lat_want = 61.25;
name = 'N_He_max';
name2 = 'HELIUM ENHANCEMENT AT 400 KM';
end
if feat == 2 
lon_want = 85;      % South nighttime maximum feature
lat_want = -58.75;
name = 'S_He_max';
name2 = 'HELIUM ENHANCEMENT AT 400 KM';
end
if feat == 3 
lon_want = -77.5;       % North Daytime minimum feature
lat_want = 18.75;
name = 'N_He_min';
name2 = 'HELIUM DEPLETION AT 400 KM';
end
if feat == 4
lon_want = -77.5;       % South Daytime minimum feature
lat_want = -46.25;      % true minimum is at -13.75 in the Southern hemi!!
% lat_want = -13.75;
name = 'S_He_min';
name2 = 'HELIUM DEPLETION AT 400 KM';

end
% ------ ETA ------ 
if feat == 6
lon_want = -150;   %  afternoon (15 LT) ETA north crest
lat_want = 11.25;  
end
if feat == 7
lon_want = -150;   %  afternoon (15 LT) ETA south crest
lat_want = -11.25;  
end
if feat == 8
lon_want = -150;   %  afternoon (15 LT) ETA trough
lat_want = -3.75;  
end
% -----------------
if feat == 9
lon_want = 95;      % Select your own 
lat_want = 1.25;
name = 'Eq_He_max';
name2 = 'HELIUM ENHANCEMENT AT 400 KM';
end
%------------------
atom_unit = 1.6605e-27; % kg/unit
k = 1.38064852e-23;% Boltzman's Constant
Av = 6.022141*10^23;
kmol = 6.02214*10^26;
mmw_he = 0.004; %Helium atomic mass kg/mol
mmw_N2 = 0.02801; %N2 molecular mass
mmw_O1 = 0.016; % O1 molecular mass
mmw_O2 = 0.032; %O2 molecular mass

%-----------------

%-----Loading Viki's tiegcm simulation-----
% format is (lon,lat,ilev,UT) format
if pdrag == 1
    filename = [aa1, 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc'];
    id = 'pdrag';
end
if pdrag == 0
    filename = [aa1, 'HSUVW.tiegcm2.0_dres.nodragtest_ctrSS_f107_180_001.nc'];
    id = 'no Ion Drag';
end

% ----- Parse Data from file ----- %
den = ncread(filename,'DEN');
zp = ncread(filename, 'Z');             % geopotential height [cm] on ilev
p0 = ncread(filename, 'p0_model');      % p0 reference pressure used by model [microbars]
g0 = ncread(filename, 'grav')/100;      % const. gravitational acceleration [m/s]
iplvl = ncread(filename, 'ilev');       % interface pressure level

he = ncread(filename,'HE'); % Units of mass mixing ratio
n2 = ncread(filename,'N2');
o1 = ncread(filename,'O1');
o2 = ncread(filename,'O2');
tn = ncread(filename,'TN');
mbar = 1./(he/4+n2/28+o1/16+o2/32); % Getting mean molecular mass from mmr
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
wn = ncread(filename,'WN');         % neutral vertical winds [cm/s] ilev surfaces


%fixed latitude and longitude  
i_index=find(lat==lat_want);

% Condense to selected latitude and UT time
den_alt_1 = squeeze(den(:,i_index,:,ut_want+1));        % Total neutral density [g/cm^2]
he_alt_1 = squeeze(he(:,i_index,:,ut_want+1));          % TIEGCM output is mass mixing ratio
N2_alt_1 = squeeze(n2(:,i_index,:,ut_want+1));
O1_alt_1 = squeeze(o1(:,i_index,:,ut_want+1));
O2_alt_1 = squeeze(o2(:,i_index,:,ut_want+1));
tn_alt_1 = squeeze(tn(:,i_index,:,ut_want+1));          % Neutral temperature [K]
mbar_alt_1 = squeeze(mbar(:,i_index,:,ut_want+1));      % Mean molecular mass
wn_alt_1 = squeeze(wn(:,i_index,:,ut_want+1));          % model winds [cm/s]

% -----Condense to selected Longitude-----
i_index = find(lon==lon_want);
den_alt = squeeze(den_alt_1(i_index,:));
He_alt = squeeze(he_alt_1(i_index,:));
N2_alt = squeeze(N2_alt_1(i_index,:));
O1_alt = squeeze(O1_alt_1(i_index,:));
O2_alt = squeeze(O2_alt_1(i_index,:));
tn_alt = squeeze(tn_alt_1(i_index,:));
mbar_alt = squeeze(mbar_alt_1(i_index,:));  %kg/kmol
wn_alt = squeeze(wn_alt_1(i_index,:));

p0_Pa = p0 / 10;  % TIEGCM gives pressure in terms of microbars, convert to Pascals
p_Pa = p0_Pa .* exp(-iplvl);        % [Pa]

% Number Density
nhe = He_alt.*den_alt*1e3/4/atom_unit;      % He number density per m^3
nN2 = N2_alt.*den_alt*1e3/28.01/atom_unit;  % N2 number density per m^3
nO1 = O1_alt.*den_alt*1e3/16/atom_unit;     % O1 number density per m^3
nO2 = O2_alt.*den_alt*1e3/32/atom_unit;     % O2 number density per m^3

% Mass Density
mhe = He_alt.*den_alt;  % Helium Mass Density g/cm^3 
mN2 = N2_alt.*den_alt;  % N2 Mass Density
mO1 = O1_alt.*den_alt;  % Oxygen Mass Density
mO2 = O2_alt.*den_alt;



%% Scale Height Calculations units of [km]
% scale height calculations for altitude slices! i.e. arrays.
[H_He_star, H_He_diff, H_tot_star, H_T, H_N2_star, H_N2_diff,...
    H_O1_star, H_O1_diff, geop_alt, meanmass, H_He_diff_no_alpha, H_tot_diff, H_P, H_m] = Scaleheight_calc_V2HLH(lon_want, lat_want, ut_want);

% all outputs on are ilev
% scale height calculation for latitude slice (i.e. 144x57 matrices) ON
% ILEVS!!!
[N2_mmr_1, O2_mmr_1, O1_mmr_1, He_mmr_1, H_He_star_ilev_1, H_He_diff_ilev_1, H_tot_star_ilev_1, H_T_ilev_1, H_N2_star_ilev_1, H_N2_diff_ilev_1,...
    H_O1_star_ilev_1, H_O1_diff_ilev_1, H_O2_star_ilev_1, H_O2_diff_ilev_1, geop_alt_1, meanmass_ilev_1, H_tot_diff_ilev_1, H_P_ilev_1, H_m_1] = Scaleheight_calc_V3HLH(lat_want, ut_want);
% ----- TESTing____


%%-----Diffusion Calculation-----%%
% Calculating the Di coefficient using TIEGCM 
D_N2 = Di_TIEGCM('N2', tn_alt, p_Pa');      %[m^2/s]
D_O2 = Di_TIEGCM('O2', tn_alt, p_Pa');
D_O1 = Di_TIEGCM('O1', tn_alt, p_Pa');
D_He = Di_TIEGCM('He', tn_alt, p_Pa');

D_Ti = -0.36*D_He; % use new version

winds_model = wn_alt/100;       %[m/s]   


%% EXTRAS BEFORE PLOTTING....

% PUT EVERYTHING ON ILEVS
% first get mass mixing ratios onto ilevs 
num_levs = length(He_alt);
sizeofMat = size(N2_mmr_1);
N2_mmr_ilev = CONVERT2ILEV(N2_alt, 0, num_levs);        % 1x57 arrays!!
O2_mmr_ilev = CONVERT2ILEV(O2_alt, 0, num_levs);
O1_mmr_ilev = CONVERT2ILEV(O1_alt, 0, num_levs);
He_mmr_ilev = CONVERT2ILEV(He_alt, 0, num_levs);
meanmass_ilev = CONVERT2ILEV(meanmass, 0, num_levs);

% find gradients of MMR's with respect to the pressure coord
N2_mmr_grad = threepointgradient(iplvl, N2_mmr_ilev, 0);
O2_mmr_grad = threepointgradient(iplvl, O2_mmr_ilev, 0);
O1_mmr_grad = threepointgradient(iplvl, O1_mmr_ilev, 0);
He_mmr_grad = threepointgradient(iplvl, He_mmr_ilev, 0);

N2_mmr_grad_1 = zeros(sizeofMat);       % 144x57 matrices of mmr gradients
O2_mmr_grad_1 = zeros(sizeofMat);
O1_mmr_grad_1 = zeros(sizeofMat);
He_mmr_grad_1 = zeros(sizeofMat);

omega_ez_grad_1 = zeros(sizeofMat); % gradient of the omega x e^-z "winds" (used below!)
omega_grad_1 = zeros(sizeofMat);    % real gradient of omega

omega = winds_model./(H_P_ilev_1(i_index, :).*1000);         % omega [1/s]
omega_1 = (wn_alt_1/100)./(H_P_ilev_1.*1000);                % omega for every lon/alt [1/s]

% make a 144x57 matrix with each whole column being one of the iplvls. (i.e. rows are all the same for a given column)
iplvl_1 = repmat(iplvl', size(He_mmr_1,1), 1 ); 

% find gradients for MMR matrices.
for l = 1:144
    N2_mmr_grad_1(l,:) = threepointgradient(iplvl_1(l,:), N2_mmr_1(l,:), 0);
    O2_mmr_grad_1(l,:) = threepointgradient(iplvl_1(l,:), O2_mmr_1(l,:), 0);
    O1_mmr_grad_1(l,:) = threepointgradient(iplvl_1(l,:), O1_mmr_1(l,:), 0);
    He_mmr_grad_1(l,:) = threepointgradient(iplvl_1(l,:), He_mmr_1(l,:), 0);  
    
    omega_ez_grad_1(l,:) = threepointgradient(iplvl_1(l,:), ( exp(-iplvl_1(l,:)).*omega_1(l,:) ), 0);   
    omega_grad_1(l,:) = threepointgradient(iplvl_1(l,:), omega_1(l,:), 0);   
end

% find mass densities of each species on the ilevs
% mN2_ilev = N2_mmr_ilev'.*den_alt.*1E3;    % N2 Mass Density [kg/m^3]
% mO2_ilev = O2_mmr_ilev'.*den_alt.*1E3;
% mO1_ilev = O1_mmr_ilev'.*den_alt.*1E3;    % Oxygen Mass Density
% mHe_ilev = He_mmr_ilev'.*den_alt.*1E3;    % Helium Mass Density [kg/m^3] 

% calculate gradient of vertical wind with height
% positive gradient = convergence of horizontal wind field
% negatie gradient = divergence of horizontal wind field

% percent difference of diffusive scale height versus actual.
H_He_perc_1 = ((H_He_star_ilev_1 ./ H_He_diff_ilev_1) - 1) * 100;
H_N2_perc_1 = ((H_N2_star_ilev_1 ./ H_N2_diff_ilev_1) - 1) * 100;

vert_adv = p0_Pa / g0 .* exp(-iplvl').* omega .* He_mmr_grad;   
vert_adv_1 = p0_Pa / g0 .* exp(-iplvl_1) .* omega_1 .* He_mmr_grad_1;
N2_vert_adv_1 = p0_Pa / g0 .* exp(-iplvl_1) .* omega_1 .* N2_mmr_grad_1;

he_mass_flux_div = -p0_Pa / g0 .* He_mmr_ilev .* threepointgradient(iplvl, ( exp(-iplvl').*omega ), 0);
he_mass_flux_div_1 = -p0_Pa / g0 .* He_mmr_1 .* omega_ez_grad_1;
N2_mass_flux_div_1 = -p0_Pa / g0 .* N2_mmr_1 .* omega_ez_grad_1;

total_mass_div_1 = -exp(iplvl_1) .* omega_ez_grad_1;     % total mass divergence 


%% -------- Plotting --------

%% 1. Diffusive & Actual Scale Heights vs Alt
% figure
% hold on;
% plot(H_N2_diff_ilev(5:50), geop_alt(5:50), '--r' , 'Linewidth', 1.5);
% plot(H_N2_star_ilev(5:50), geop_alt(5:50), 'r' , 'Linewidth', 1.5);
% plot(H_O2_diff_ilev(5:50), geop_alt(5:50), '--b' , 'Linewidth', 1.5);
% plot(H_O2_star_ilev(5:50), geop_alt(5:50), 'b' , 'Linewidth', 1.5);
% plot(H_O1_diff_ilev(5:50), geop_alt(5:50), '--k' , 'Linewidth', 1.5);
% plot(H_O1_star_ilev(5:50), geop_alt(5:50), 'k' , 'Linewidth', 1.5);
% plot(H_tot_star_ilev(5:50), geop_alt(5:50), 'g' , 'Linewidth', 1.5); % already on interfaces
% xlabel('Scale Height [km]');
% ylabel('Altitude [km]');
% % xlim([0 350])
% ylim([geop_alt(5) 550])
% grid on;
% title('Neutral Species Scale Heights (ilev surfaces)');
% legend('N2 diff','N2 actual','O2 diff','O2 actual','O1 diff','O1 actual','Total Gas', 'location', 'best');

% figure
% hold on
% plot(H_He_diff_ilev(5:50), geop_alt(5:50), '--k' , 'Linewidth', 1.5);
% plot(H_He_star_ilev(5:50), geop_alt(5:50), 'k' , 'Linewidth', 1.5);
% plot(H_tot_star(5:50), geop_alt(5:50), 'g' , 'Linewidth', 1.5);
% xlabel('Scale Height [km]');
% ylabel('Altitude [km]');
% ylim([geop_alt(5) 550])
% title('Helium Scale Heights (ilev surfaces)');
% legend('HE diff','HE actual','Total Gas', 'location', 'best');
% grid on;

%% 2. Diffusion Coeff Vs Alt
% figure
% top = 150;
% plot(D_He,geop_alt,'Linewidth',2);
% hold on
% plot(D_O2,geop_alt,'Linewidth',2);
% hold on
% plot(D_N2,geop_alt,'Linewidth',2);
% hold on
% plot(D_O1,geop_alt,'Linewidth',2);
% hold on
% % plot(K_e,geop_alt,'Linewidth',2);
% set(gca, 'XScale', 'log')
% xlabel('Diffusion Coefficients [m^2/s]');
% ylabel('Altitude (km)');
% title('Molecular Diffusion Coefficient - NEW TIEGCM VERSION');
% ylim([95 top]);
% grid on;
% legend('He','O2','N2','O1','Eddy', 'location', 'best');
% 
 
%% 3. TIEGCM output winds vs Alt
% plot(winds_model(5:35), geop_alt(5:35),'Linewidth', 1.5);
% xlabel('Vertical wind [m/s]');
% ylabel('Altitude [km]');
% title('Vertical Wind for Helium Diffusion');
% % set(gca, 'XScale', 'log');
% legend('He Diffusive Balance Winds', 'TIEGCM Output Winds');

%% 5. He Scale Height w and w/o thermal diffusion 
% figure
% plot(H_He_diff(5:50), geop_alt(5:50),'Linewidth', 1.5);
% hold on
% plot(H_He_diff_no_alpha(5:50), geop_alt(5:50),'Linewidth', 1.5);
% title('Helium Diffusive Equilibrium Scale Height');
% ylabel('Altitude [km]');
% xlabel('Scale Height [km]');
% legend('With Thermal Diffusion', 'No Thermal Diffusion');

%% 6. MMR vs Alt
% figure()
% subplot(121);
% plot(N2_mmr_ilev(1:50), geop_alt(1:50), 'Linewidth', 1.5);
% hold on
% plot(O2_mmr_ilev(1:50), geop_alt(1:50), 'Linewidth', 1.5);
% hold on
% plot(O1_mmr_ilev(1:50), geop_alt(1:50), 'Linewidth', 1.5);
% hold on
% plot(He_mmr_ilev, geop_alt, 'Linewidth', 1.5);
% 
% legend('N2','O2','O1','He', 'location', 'best');
% title(['TIEGCM Mass Mixing Ratios for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([120 700])
% xlabel('Mass Mixing Ratio');
% ylabel('Geopotential Altitude [km]');
% grid on;

%% 7. grad(MMR) vs Alt
% subplot(122)
% plot(N2_mmr_grad(1:50), geop_alt(1:50), 'Linewidth', 1.5);
% hold on
% plot(O2_mmr_grad(1:50), geop_alt(1:50), 'Linewidth', 1.5);
% hold on
% plot(O1_mmr_grad(1:50), geop_alt(1:50), 'Linewidth', 1.5);
% hold on
% plot(He_mmr_grad(1:50), geop_alt(1:50), 'Linewidth', 1.5);
% 
% legend('N2','O2','O1','He', 'location', 'best');
% title(['TIEGCM Mass Mixing Ratio Gradients for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([100 550])
% xlabel('Mass Mixing Ratio Gradient [1/km]');
% ylabel('Geopotential Altitude [km]');
% grid on;

%% 8. Mean Molecular Weight vs Alt
% figure
% plot(meanmass_ilev*1000, geop_alt)
% title(['TIEGCM Mean Molecular Mass for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([100 550])
% xlabel('Mean Molecular Mass [g/mol]');
% ylabel('Geopotential Altitude [km]');
% grid on;

%% 9. MMR scale heights vs Alt
% figure
% hold on
% plot(H_He_mmr(1:45), geop_alt(1:45), 'Linewidth', 1.5);
% plot(H_N2_mmr(15:45), geop_alt(15:45), 'Linewidth', 1.5);
% plot(H_O2_mmr(1:45), geop_alt(1:45), 'Linewidth', 1.5);
% legend('He','N2','O2');
% title('Mass Mixing Ratio Scale Heights');
% ylabel('Altitude [km]');
% 

%% 10. Compare Wind Term to Diffusion Terms using two methods 
% % from Eq. 5.11 in Torfinn's thesis.
% % METHOD 1 -  see 8/5/2019 entry in my gray/black notebook
% top = 200;
% bot = 95;
% figure()
% subplot(221)
% plot(abs(diff_term_N2), geop_alt);
% hold on;
% plot(abs(wind_term_N2), geop_alt);
% set(gca, 'XScale', 'log');
% legend('N2 diff', 'N2 wind',  'location', 'best');
% title(['TIEGCM N2 Wind and Diffusion Term Magnitutdes for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([bot top])
% ylabel('Geopotential Altitude [km]');
% grid on;
% 
% subplot(222)
% plot(abs(diff_term_He), geop_alt);
% hold on;
% plot(abs(wind_term_He), geop_alt);
% set(gca, 'XScale', 'log');
% legend('He diff', 'He wind', 'location', 'best');
% title(['TIEGCM He Wind and Diffusion Term Magnitutdes for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([bot top])
% ylabel('Geopotential Altitude [km]');
% grid on;
% 
% subplot(223)
% plot(abs(diff_term_O2), geop_alt);
% hold on;
% plot(abs(wind_term_O2), geop_alt);
% set(gca, 'XScale', 'log');
% legend('O2 diff', 'O2 wind',  'location', 'best');
% title(['TIEGCM O2 Wind and Diffusion Term Magnitutdes for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([bot top])
% ylabel('Geopotential Altitude [km]');
% grid on;
% 
% subplot(224)
% plot(abs(diff_term_O1), geop_alt);
% hold on;
% plot(abs(wind_term_O1), geop_alt);
% set(gca, 'XScale', 'log');
% legend('O1 diff', 'O1 wind', 'location', 'best');
% title(['TIEGCM O1 Wind and Diffusion Term Magnitutdes for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([bot top])
% ylabel('Geopotential Altitude [km]');
% grid on;
% sgtitle('METHOD 1')
% 
% % METHOD 2 -  see 8/5/2019 entry in notebook
% 
% figure()
% subplot(221)
% plot(abs(diff_term_N2x), geop_alt);
% hold on;
% plot(abs(wind_term_N2x), geop_alt);
% set(gca, 'XScale', 'log');
% legend('N2 diff', 'N2 wind',  'location', 'best');
% title(['TIEGCM N2 Wind and Diffusion Term Magnitutdes for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([bot top])
% ylabel('Geopotential Altitude [km]');
% grid on;
% 
% subplot(222)
% plot(abs(diff_term_Hex), geop_alt);
% hold on;
% plot(abs(wind_term_Hex), geop_alt);
% set(gca, 'XScale', 'log');
% legend('He diff', 'He wind', 'location', 'best');
% title(['TIEGCM He Wind and Diffusion Term Magnitutdes for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([bot top])
% ylabel('Geopotential Altitude [km]');
% grid on;
% 
% subplot(223)
% plot(abs(diff_term_O2x), geop_alt);
% hold on;
% plot(abs(wind_term_O2x), geop_alt);
% set(gca, 'XScale', 'log');
% legend('O2 diff', 'O2 wind',  'location', 'best');
% title(['TIEGCM O2 Wind and Diffusion Term Magnitutdes for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([bot top])
% ylabel('Geopotential Altitude [km]');
% grid on;
% 
% subplot(224)
% plot(abs(diff_term_O1x), geop_alt);
% hold on;
% plot(abs(wind_term_O1x), geop_alt);
% set(gca, 'XScale', 'log');
% legend('O1 diff', 'O1 wind', 'location', 'best');
% title(['TIEGCM O1 Wind and Diffusion Term Magnitutdes for Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want)]);
% ylim([bot top])
% ylabel('Geopotential Altitude [km]');
% grid on;
% sgtitle('METHOD 2')
% 

%% ---- Comparing Terms in Momentum Eq ------

%% 11. 1D Altitude Profiles
% ytop = geop_alt(end-2);
% ybot = 120;
% logic = geop_alt > ybot & geop_alt < ytop;
% 
% thisfig = figure('units','normalized','outerposition',[0 1 .7 1]);     % left start, bottom start, side length, tall length
% 
% % ----------------
% subplot(221)
% bounds = omega_1(:,logic);
% ztop = (max(abs(bounds(:))))*1;  % make the max 110% of the largest value 
% zbot = -ztop;
% 
% x = lon;
% y = geop_alt;
% [X, Y] = meshgrid(x, y);
% 
% num_cont = 300;
% c = redblue(300);
% colormap(c); 
% Q = omega_1';      % omega for every longitude and alt for a specific latitute [m/s]
% 
% contourf(X, Y, Q, num_cont, 'linecolor', 'none')
% hold on;   
% plot(lon_want*ones(length(geop_alt), 1), geop_alt, 'k--', 'Linewidth', 2.5)   
% cbar = colorbar();
% caxis([zbot ztop]); 
% ylim([ybot ytop]);
% cbar.Label.String = 'Omega [1/s]';
% xlabel('Longtitude [Deg]');
% ylabel('Geopotential Altitude [km]');
% title('TIEGCM Omega Output (ilev surfaces)')
% 
% % ----------------
% subplot(222)
% hold on
% plot(H_He_diff_ilev_1(i_index, :), geop_alt, '--k' , 'Linewidth', 1.5);
% plot(H_He_star_ilev_1(i_index, :), geop_alt, 'k' , 'Linewidth', 1.5);
% xlabel('Scale Height [km]');
% ylabel('Geopotential Altitude [km]');
% ylim([ybot ytop]);
% title('Helium Scale Heights (ilev surfaces)');
% legend('HE diff','HE actual', 'location', 'best');
% grid on;
% box on;
% ax = gca;
% ax.XRuler.Axle.LineWidth = 2;
% ax.YRuler.Axle.LineWidth = 2;
% 
% % ----------------
% subplot(223)
% xtop = (max(abs(vert_adv_1(i_index,logic))))*1.10;  % make the max 110% of the largest value 
% xbot = -xtop;
% xdis = (xtop - xbot)/2;
% hold on;
% r1 = rectangle('Position',[xbot ybot xdis ytop]);
% r1.FaceColor = [.85 .85 1];     % very light blue
% r2 = rectangle('Position',[0 ybot xdis ytop]);
% r2.FaceColor = [1 .9 .9];       % very light red
% plot(NaN, 'color', [.85 .85 1], 'Linewidth', 10);
% plot(NaN, 'color', [1 .9 .9], 'Linewidth', 10);
% plot(vert_adv_1(i_index, :), geop_alt, 'k-', 'Linewidth', 2);
% hold off;
% grid on;
% set(gca, 'Layer', 'top')
% legend({'- Advection','+ Advection'},'Fontsize', 12, 'location', 'best')
% ylim([ybot ytop]);
% xlim([xbot xtop]);
% ylabel('Geopotential Altitude [km]');
% xlabel('Vertical Advection [kg/(m s) * m^2]');
% title('Vertical Advection (ilev surfaces)');
% box on;
% ax = gca;
% ax.XRuler.Axle.LineWidth = 2;
% ax.YRuler.Axle.LineWidth = 2;
% 
% % ----------------
% subplot(224)
% xtop = (max(abs(he_mass_flux_div_1(i_index,logic))))*1.10;
% xbot = -xtop;
% xdis = (xtop - xbot)/2;
% hold on;
% r1 = rectangle('Position',[xbot ybot xdis ytop]);
% r1.FaceColor = [.85 .85 1];     % very light blue
% r2 = rectangle('Position',[0 ybot xdis ytop]);
% r2.FaceColor = [1 .9 .9];       % very light red
% plot(NaN, 'color', [.85 .85 1], 'Linewidth', 10);
% plot(NaN, 'color', [1 .9 .9], 'Linewidth', 10);
% plot(he_mass_flux_div_1(i_index, :), geop_alt, 'k-', 'Linewidth', 2)
% hold on;
% grid on;
% set(gca, 'Layer', 'top')
% legend({'Convergence','Divergence'},'Fontsize', 12, 'location', 'best')
% ylim([ybot ytop]);
% xlim([xbot xtop]);
% ylabel('Geopotential Altitude [km]');
% xlabel('Mass Flux Divergence [kg/(m s) * m^2]')
% title('Helium Mass Flux Divergence (ilev surfaces)');
% box on;
% ax = gca;
% ax.XRuler.Axle.LineWidth = 2;
% ax.YRuler.Axle.LineWidth = 2;
% 
% sgtitle(['TIEGCM Model Run at UT=0, Lat = ', num2str(lat_want) , ', Lon = ', num2str(lon_want),', ', name2]);
% 
% % % saveas(thisfig, ['./Figures/He_mass_flux/TIEGCM_vert_mass_flux_', name, '.png'])
% % 

%% 12. 2D Longitude vs. Alt profiles
% % this is the same as figure 14 above, except everything has moved from a
% % altitude profile to a latitude slice - i.e. contour plots with lon, alt,
% % and variable on the x, y and z axis respectively.
% 
% % ###################################################################################
% % ###################################################################################
% 
% % since the geopotential altitude is different for every longitude value,
% % let your y value be the average of each column 

y_avg = mean(geop_alt_1, 1);        % row vector with average of each column.
x = lon;
[X, Y] = meshgrid(x, y_avg);        % mesh grid used for every subplot

ytop = y_avg(end-2);                % For N2, do end-5??
ybot = 100;
logic = y_avg >= ybot & y_avg <= ytop;

X = X(logic,:);              % needed to do this b/c had trouble plotting
Y = Y(logic,:);

num_cont = 300;
c = redblue(300);

thisfig = figure('units','normalized','outerposition',[0 1 .7 1]);     % left start, bottom start, side length, tall length

% ----------------
subplot(221)
Q = omega_1(:,logic)';      % omega for every longitude and alt for a specific latitute [m/s]
ztop = (max(abs(Q(:))));  % make the max 100% of the largest value 
zbot = -ztop;

colormap(c); 
contourf(X, Y, Q, num_cont, 'linecolor', 'none')
hold on;   
plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5)   
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = '\omega [1/s]';
xlabel('Longtitude [Deg]');
ylabel('Geopotential Altitude [km]');
title('TIEGCM Omega Output (ilev surfaces)')
grid on;

% ----------------
subplot(222)
% Q = H_N2_perc_1(:,logic)';
% Q = H_He_perc_1(:,logic)';
Q = omega_grad_1(:,logic)';
ztop = (max(abs(Q(:))))*1;  
zbot = -ztop;

colormap(c);
contourf(X, Y, Q, num_cont, 'linecolor', 'none')
hold on;   
plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5) 
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = 'Scale Height Perc. Difference';
% cbar.Label.String = '\partial \omega / \partial z [1/s]';
xlabel('Longtitude [Deg]');
ylabel('Geopotential Altitude [km]');
% title('N2 Scale Height Perc. Difference from Diffusive Eq. (ilev surfaces)');
title('Gradient of Omega (ilev surfaces)');
grid on;

% ----------------
subplot(223)
% Q = N2_vert_adv_1(:,logic)';
% Q = vert_adv_1(:,logic)';
Q = omega_1(:,logic)' - omega_grad_1(:,logic)';
ztop = (max(abs(Q(:))))*1;  
zbot = -ztop;

colormap(c);
contourf(X, Y, Q, num_cont, 'linecolor', 'none')
hold on;   
plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5) 
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = 'Vertical Advection [kg/(m s) * m^2]';
% cbar.Label.String = '\omega - \partial \omega/ \partial z';

xlabel('Longtitude [Deg]');
ylabel('Geopotential Altitude [km]');
% title('N2 Vertical Advection (+ Up) (ilev surfaces)');
title('\omega - \partial \omega/ \partial z')
grid on;

% % ----------------
subplot(224)
% Q = N2_mass_flux_div_1(:,logic)';
% Q = he_mass_flux_div_1(:,logic)';
Q = total_mass_div_1(:,logic)';      % postive divergence = negative sign
ztop = (max(abs(Q(:))))*1;   
zbot = -ztop;

colormap(c);
contourf(X, Y, Q, num_cont, 'linecolor', 'none')
hold on;   
plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5) 
hold on;
for i =(14:4:50)
    plot(x, geop_alt_1(:,i));
    
    %label the lines
    xpt = 155;
    ypt = double(geop_alt_1(end, i));
    lbl = ['z = ', num2str(iplvl(i))];
    text(xpt, ypt, lbl, 'FontSize', 8, 'BackgroundColor', 'white', 'HorizontalAlignment', 'Center')
    hold on;   
end
hold off;
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
% cbar.Label.String = 'N2 Mass Flux Divergence [kg/(m s) * m^2]';
cbar.Label.String = 'Total Divergence [1/s]';
xlabel('Longtitude [Deg]');
ylabel('Geopotential Altitude [km]');
% title('N2 Mass Flux Divergence (ilev surfaces)');
title('Total Divergence (ilev surfaces)');
grid on;

sgtitle(['TIEGCM Model Run at UT=0, Lat = ', num2str(lat_want), ', ', name2]);

% % saveas(thisfig, ['./Figures/Total_mass_flux/TIEGCM_lon_alt_', name, '.png'])
% saveas(thisfig, ['./Figures/N2_mass_flux/TIEGCM_lon_alt_', name, '.png'])

%% 13. Percentages of Total Wind Field Div

% A = abs(omega_1(:,logic)');
% B = abs(omega_grad_1(:,logic)');
% C = A + B;
% 
% perc_A = abs(A./C .* 100);
% perc_B = abs(B./C .* 100);
% 
% thisfig = figure('units','normalized','outerposition',[0 1 .7 .5]);     % left start, bottom start, side length, tall length
% 
% subplot(121);
% c = jet(200);
% colormap(c);
% contourf(X, Y, perc_A, 300, 'linecolor', 'none');
% hold on;   
% plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5)   
% cbar = colorbar();
% title('Percentage \omega')
% caxis([0 100]); 
% ylim([ybot ytop]);
% 
% subplot(122);
% c = jet(200);
% colormap(c);
% contourf(X, Y, perc_B, 300, 'linecolor', 'none');
% hold on;   
% plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5)
% cbar = colorbar();
% ylim([ybot ytop]);
% title('Percentage \partial \omega / \partial z')
% caxis([0 100]);
% 
% sgtitle(['TIEGCM Model Run at UT=0, Lat = ', num2str(lat_want), ', ', name2]);
% 
% saveas(thisfig, ['./Figures/TIEGCM_Percentages_', name, '.png'])

%% 14. Long/Alt Slice of H and p0 e^-z / g0

% ytop = geop_alt(end) + 10;
% ybot = 100;
% logic = geop_alt > ybot & geop_alt < ytop;
% 
% x = lon;
% y = geop_alt;
% 
% thisfig = figure('units','normalized','outerposition',[0 1 .7 1]);
% for i =(10:4:57)
%     plot(x, geop_alt_1(:,i));
%     
%     % label the lines
%     xpt = 160;
%     ypt = double(geop_alt_1(end, i));
%     lbl = ['z = ', num2str(iplvl(i))];
%     text(xpt, ypt, lbl, 'FontSize', 10, 'BackgroundColor', 'white', 'HorizontalAlignment', 'Center')
%     hold on;   
% end
% 
% plot(lon_want*ones(length(y), 1), y, 'k--', 'Linewidth', 2.5) 
% ylim([ybot ytop]);
% xlim([-180 180])
% xlabel('Longtitude [Deg]');
% ylabel('Geopotential Altitude [km]');
% grid on;
% 
% sgtitle(['TIEGCM Model Run Pressure Levels at UT = 0, Lat = ', num2str(lat_want), ', ', name2]);

% saveas(thisfig, ['./Figures/PressureLvl_plot_', name, '.png'])

%% 15. rho * H and p0/g0 * e^(-iplvl) vs. Alt

plot1 = (H_P_ilev_1*1000) .* den_alt_1 * (100^3)/1000;       % [kg / m^2]
plot2 = p0_Pa / g0 .* exp(-iplvl_1);

Q1 = ((plot2 ./ plot1) - 1) * 100 ;

y_avg = mean(geop_alt_1, 1);        % row vector with average of each column.
x = lon;
[X, Y] = meshgrid(x, y_avg);        % mesh grid used for every subplot

ytop = y_avg(end-2);                % For N2, do end-5??
ybot = 100;
logic = y_avg >= ybot & y_avg <= ytop;

X = X(logic,:);              % needed to do this b/c had trouble plotting
Y = Y(logic,:);

num_cont = 300;
c = redblue(300);

thisfig = figure('units','normalized','outerposition',[0 1 1 1]);     % left start, bottom start, side length, tall length

% ----------------
Q = Q1(:,logic)';          
ztop = 15;  
zbot = -15;

colormap(c); 
contourf(X, Y, Q, num_cont, 'linecolor', 'none')
hold on;   
plot(lon_want*ones(length(y_avg), 1), y_avg, 'k--', 'Linewidth', 2.5)   
cbar = colorbar();
caxis([zbot ztop]); 
ylim([ybot ytop]);
cbar.Label.String = 'Perc Diff. btw \rho* H and p0/g0 e^{-z}';
xlabel('Longtitude [Deg]');
ylabel('Geopotential Altitude [km]');
title('TIEGCM Percent Difference Output (ilev surfaces)')
grid on;








 
