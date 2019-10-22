function [N2_lz, O2_lz, O1_lz, He_lz, H_He_star, H_He_diff, H_tot_star, H_temp, H_N2_star, H_N2_diff,...
    H_O1_star, H_O1_diff, H_O2_star, H_O2_diff, zp_lz, meanmass, H_tot_diff, Hp_mean, H_mass] = Scaleheight_calc_V3HLH(lat_want, ut_want)
%Calculates scale heights and parses .nc file data. Can be used to condense
%main program code.
% CALCULATES EVERYTHING in terms of a MATRIX
% Calculates all the scale heights on ilevs !!!

filename = 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc';

k = 1.38e-23;         % Boltzman's Constant
Av = 6.022141*10^23;
mmw_he = 0.004;      % Helium atomic mass [kg]
mmw_N2 = 0.028;     % N2 molecular mass
mmw_O1 = 0.016;     % O1 molecular mass
mmw_O2 = 0.032;     % O2 molecular mass 

den = ncread(filename,'DEN');
zp = ncread(filename, 'Z'); 
g0 = ncread(filename, 'grav')/100;          % const. gravitational acceleration [m/s]
he = ncread(filename,'HE');                 % Units of mass mixing ratio
n2 = ncread(filename,'N2');
o1 = ncread(filename,'O1');
o2 = ncread(filename,'O2');
tn = ncread(filename,'TN');
mbar = 1./(he/4+n2/28+o1/16+o2/32); % Getting mean molecular mass from mmr
lat = ncread(filename,'lat');

%fixed latitude and longitude  
i_index = find(lat==lat_want);

% fixed UT time and latitude 
den_lz = squeeze(den(:,i_index,:,ut_want+1));       % Total neutral density ilev
zp_lz = squeeze(zp(:,i_index,:,ut_want+1))/1e5;     % Geopotential height ilev [km]
He_lz_lev = squeeze(he(:,i_index,:,ut_want+1));     % Helium mass mixing ratio
N2_lz_lev = squeeze(n2(:,i_index,:,ut_want+1));
O1_lz_lev = squeeze(o1(:,i_index,:,ut_want+1));
O2_lz_lev = squeeze(o2(:,i_index,:,ut_want+1));
tn_lz = squeeze(tn(:,i_index,:,ut_want+1));         % Neutral temperature [K]
mbar_lz_lev = squeeze(mbar(:,i_index,:,ut_want+1)); % Mean molecular mass

meanmass_lev = mbar_lz_lev/1000;  % [kg/mol]


% the only thing not on the ilevs is the mmrs... CONVERT THEM ALL
lonPts = 144;
altPts = 57;
He_lz = CONVERT2ILEV(He_lz_lev, lonPts, altPts);
N2_lz = CONVERT2ILEV(N2_lz_lev, lonPts, altPts);
O1_lz = CONVERT2ILEV(O1_lz_lev, lonPts, altPts);
O2_lz = CONVERT2ILEV(O2_lz_lev, lonPts, altPts);
meanmass = CONVERT2ILEV(meanmass_lev, lonPts, altPts);  % [kg/mol]

% --------- EVERYTHING IS NOW ON ILEVS --------------

% Mass Density
mhe = He_lz.*den_lz;% Helium Mass Density [g/cm^3] 
mN2 = N2_lz.*den_lz;% N2 Mass Density
mO1 = O1_lz.*den_lz;% Oxygen Mass Density
mO2 = O2_lz.*den_lz;

%-----Finding Pressure Scale Heights for gas constituents from kT/mg-----
Hp_mean = k.*tn_lz./(meanmass/Av.*g0)/1000;
Hp_he = k.*tn_lz./(mmw_he/Av.*g0)/1000;   % [Km]
Hp_n2 = k.*tn_lz./(mmw_N2/Av.*g0)/1000;
Hp_o1 = k.*tn_lz./(mmw_O1/Av.*g0)/1000;
Hp_o2 = k.*tn_lz./(mmw_O2/Av.*g0)/1000;

% -----Calculate Scale Heights Using 3 Point Differentiation Technique-----
points = zeros(size(tn_lz));
H_He_star = points;
H_dentot = points;
H_n2 = points;
H_o1 = points;
H_o2 = points;
H_tn = points;
H_mass = points;

for l = 1:lonPts
    for z = 1:altPts
        if z == 1   % First Point gradient technique
            coeff1 = (2*zp_lz(l,1)-zp_lz(l,2)-zp_lz(l,3))/((zp_lz(l,1)-...
                zp_lz(l,2))*(zp_lz(l,1)-zp_lz(l,3)));

            coeff2 = (2*zp_lz(l,1)-zp_lz(l,1)-zp_lz(l,3))/((zp_lz(l,2)-...
                zp_lz(l,1))*(zp_lz(l,2)-zp_lz(l,3)));

            coeff3 = (2*zp_lz(l,1)-zp_lz(l,1)-zp_lz(l,2))/((zp_lz(l,3)-...
                zp_lz(l,1))*(zp_lz(l,3)-zp_lz(l,2)));

            H_He_star(l,1) = -1/mhe(l,1)*(mhe(l,1)*coeff1+mhe(l,2)*coeff2+...
                mhe(l,3)*coeff3);
            H_dentot(l,1) = -1/den_lz(l,1)*(den_lz(l,1)*coeff1+den_lz(l,2)*coeff2+...
                den_lz(l,3)*coeff3);
            H_n2(l,1) = -1/mN2(l,1)*(mN2(l,1)*coeff1+mN2(l,2)*coeff2+...
                mN2(l,3)*coeff3);
            H_o1(l,1) = -1/mO1(l,1)*(mO1(l,1)*coeff1+mO1(l,2)*coeff2+...
                mO1(l,3)*coeff3);
            H_o2(l,1) = -1/mO2(l,1)*(mO2(l,1)*coeff1+mO2(l,2)*coeff2+...
                mO2(l,3)*coeff3);
            H_tn(l,1) = 1/tn_lz(l,1)*(tn_lz(l,1)*coeff1+tn_lz(l,2)*coeff2+...
                tn_lz(l,3)*coeff3);
            H_mass(l,1) = -1/meanmass(l,1)*(meanmass(l,1)*coeff1+meanmass(l,2)*coeff2+...
                meanmass(l,3)*coeff3);


        elseif z == altPts %Last point gradient technique
            coeff1 = (2*zp_lz(l,z)-zp_lz(l,z-1)-zp_lz(l,z))/((zp_lz(l,z-2)-...
                zp_lz(l,z-1))*(zp_lz(l,z-2)-zp_lz(l,z)));

            coeff2 = (2*zp_lz(l,z)-zp_lz(l,z-2)-zp_lz(l,z))/((zp_lz(l,z-1)-...
                zp_lz(l,z-2))*(zp_lz(l,z-1)-zp_lz(l,z)));

            coeff3 = (2*zp_lz(l,z)-zp_lz(l,z-2)-zp_lz(l,z-1))/((zp_lz(l,z)-...
                zp_lz(l,z-2))*(zp_lz(l,z)-zp_lz(l,z-1)));

            H_He_star(l,z) = -1/mhe(l,z)*(mhe(l,z-2)*coeff1+mhe(l,z-1)*coeff2+...
                mhe(l,z)*coeff3);
            H_dentot(l,z) = -1/den_lz(l,z)*(den_lz(l,z-2)*coeff1+den_lz(l,z-1)*coeff2+...
                den_lz(l,z)*coeff3);
            H_n2(l,z) = -1/mN2(l,z)*(mN2(l,z-2)*coeff1+mN2(l,z-1)*coeff2+...
                mN2(l,z)*coeff3);            
            H_o1(l,z) = -1/mO1(l,z)*(mO1(l,z-2)*coeff1+mO1(l,z-1)*coeff2+...
                mO1(l,z)*coeff3); 
            H_o2(l,z) = -1/mO2(l,z)*(mO2(l,z-2)*coeff1+mO2(l,z-1)*coeff2+...
                mO2(l,z)*coeff3);
            H_tn(l,z) = 1/tn_lz(l,z)*(tn_lz(l,z-2)*coeff1+tn_lz(l,z-1)*coeff2+...
                tn_lz(l,z)*coeff3);
            H_mass(l,z) = -1/meanmass(l,z)*(meanmass(l,z-2)*coeff1+meanmass(l,z-1)*coeff2+...
                meanmass(l,z)*coeff3);


        else %Middle Points gradient technique
            coeff1 = (2*zp_lz(l,z)-zp_lz(l,z)-zp_lz(l,z+1))/((zp_lz(l,z-1)-...
                zp_lz(l,z))*(zp_lz(l,z-1)-zp_lz(l,z+1)));

            coeff2 = (2*zp_lz(l,z)-zp_lz(l,z-1)-zp_lz(l,z+1))/((zp_lz(l,z)-...
                zp_lz(l,z-1))*(zp_lz(l,z)-zp_lz(l,z+1)));

            coeff3 = (2*zp_lz(l,z)-zp_lz(l,z-1)-zp_lz(l,z))/((zp_lz(l,z+1)-...
                zp_lz(l,z-1))*(zp_lz(l,z+1)-zp_lz(l,z)));

            H_He_star(l,z) = -1/mhe(l,z)*(mhe(l,z-1)*coeff1+mhe(l,z)*coeff2+...
                mhe(l,z+1)*coeff3);
            H_dentot(l,z) = -1/den_lz(l,z)*(den_lz(l,z-1)*coeff1+den_lz(l,z)*coeff2+...
                den_lz(l,z+1)*coeff3);
            H_n2(l,z) = -1/mN2(l,z)*(mN2(l,z-1)*coeff1+mN2(l,z)*coeff2+...
                mN2(l,z+1)*coeff3);
            H_o1(l,z) = -1/mO1(l,z)*(mO1(l,z-1)*coeff1+mO1(l,z)*coeff2+...
                mO1(l,z+1)*coeff3);
            H_o2(l,z) = -1/mO2(l,z)*(mO2(l,z-1)*coeff1+mO2(l,z)*coeff2+...
                mO2(l,z+1)*coeff3);
            H_tn(l,z) = 1/tn_lz(l,z)*(tn_lz(l,z-1)*coeff1+tn_lz(l,z)*coeff2+...
                tn_lz(l,z+1)*coeff3);
            H_mass(l,z) = -1/meanmass(l,z)*(meanmass(l,z-1)*coeff1+meanmass(l,z)*coeff2+...
                meanmass(l,z+1)*coeff3);

        end
    end
end
%-----Get Scale Height from Inverse-----
H_He_star = 1./H_He_star;
H_tot_star = 1./H_dentot;
H_N2_star = 1./H_n2;
H_O1_star = 1./H_o1;
H_O2_star = 1./H_o2;
H_temp = (1./H_tn);
H_temp_he = H_temp/.62;% <---- Alpha for helium is -.38
H_mass = (1./H_mass);


%----Put Together Diffusive Profiles and Mean Mass Profile-----
H_He_diff = 1./(1./H_temp_he+1./Hp_he);% Helium diffusive profile
H_tot_diff = 1./(1./H_temp+1./Hp_mean+1./H_mass);   % Atmospheric scale height
H_N2_diff = 1./(1./H_temp+1./Hp_n2);        % N2 Diffusive profile
H_O1_diff = 1./(1./H_temp+1./Hp_o1);        % O1 Diffusive profile
H_O2_diff = 1./(1./Hp_o2+1./H_temp);        % O2 diffusie profile

end

