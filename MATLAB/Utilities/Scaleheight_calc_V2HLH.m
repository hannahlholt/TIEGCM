function [H_He_star, H_He_diff, H_tot_star, H_temp, H_N2_star, H_N2_diff,...
    H_O1_star, H_O1_diff, geop_alt, meanmass, H_He_diff_no_alpha, H_tot_diff, Hp_mean, H_mass] = Scaleheight_calc_V2HLH(lon_want, lat_want, ut_want)
%Calculates scale heights and parses .nc file data. Can be used to condense
%main program code.
% CALCULATES EVERYTHING AS AN ARRAY

filename = 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc';
linear = 1;
atom_unit=1.67e-27; % kg/unit

den = ncread(filename,'DEN');
zg = ncread(filename,'ZG'); %Geometric altitude [cm]
zp = ncread(filename, 'Z'); 
g0 = ncread(filename, 'grav')/100;          % const. gravitational acceleration [m/s]
he = ncread(filename,'HE');             % Units of mass mixing ratio
n2 = ncread(filename,'N2');
o1 = ncread(filename,'O1');
o2 = ncread(filename,'O2');
tn = ncread(filename,'TN');
mbar = 1./(he/4+n2/28+o1/16+o2/32);%Getting mean molecular mass from mmr
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
wn = ncread(filename,'WN');


% fixed UT time 
den_tp1=squeeze(den(:,:,:,ut_want+1));% Total neutral density
zg_tp1=squeeze(zg(:,:,:,ut_want+1));% Geometric height
z_tp1=squeeze(zp(:,:,:,ut_want+1));% Geopotential height
he_tp1=squeeze(he(:,:,:,ut_want+1));% Helium mass mixing ratio
N2_tp1=squeeze(n2(:,:,:,ut_want+1));
O1_tp1=squeeze(o1(:,:,:,ut_want+1));
tn_tp1=squeeze(tn(:,:,:,ut_want+1));% Neutral temperature
mbar_tp1=squeeze(mbar(:,:,:,ut_want+1));% Mean molecular mass
wn_tp1 = squeeze(wn(:,:,:,ut_want+1));


%fixed latitude and longitude  
i_index=find(lat==lat_want);

% Condense to selected latitude
den_alt_1=squeeze(den_tp1(:,i_index,:));
he_alt_1=squeeze(he_tp1(:,i_index,:));%TIEGCM output is mass mixing ratio
N2_alt_1=squeeze(N2_tp1(:,i_index,:));
O1_alt_1=squeeze(O1_tp1(:,i_index,:));
% geometric_alt_1=squeeze(zg_tp1(:,i_index,:))/1e5; %Convert from cm to km
geop_alt_1=squeeze(z_tp1(:,i_index,:))/1e5;
tn_alt_1=squeeze(tn_tp1(:,i_index,:));
mbar_alt_1=squeeze(mbar_tp1(:,i_index,:));
wn_alt_1 = squeeze(wn_tp1(:,i_index,:));

% -----Mean over the Lonitude-----
% den_alt=mean(den_alt_1);
% he_alt=mean(he_alt_1);
% N2_alt=mean(N2_alt_1);
% O1_alt=mean(O1_alt_1);
% % geop_alt=mean(geop_alt_1); 
% geom_alt=mean(geometric_alt_1);
% tn_alt=mean(tn_alt_1);
% mbar_alt=mean(mbar_alt_1);
% wn_alt = mean(wn_alt_1);

% -----Condense to selected Longitude-----
i_index = find(lon==lon_want);

den_alt = squeeze(den_alt_1(i_index,:));
he_alt = squeeze(he_alt_1(i_index,:));
N2_alt = squeeze(N2_alt_1(i_index,:));
O1_alt = squeeze(O1_alt_1(i_index,:));
geop_alt = squeeze(geop_alt_1(i_index,:));
% geom_alt = squeeze(geometric_alt_1(i_index,:));
tn_alt = squeeze(tn_alt_1(i_index,:));
mbar_alt = squeeze(mbar_alt_1(i_index,:)); %kg/kmol
wn_alt = squeeze(wn_alt_1(i_index,:));

% Number Density
nhe_bf=he_alt.*den_alt*1e3/4/atom_unit;   %helium number density
nN2=N2_alt.*den_alt*1e3/28.01/atom_unit; % N2 number density
nO1=O1_alt.*den_alt*1e3/16/atom_unit; % O1 number density

%Mass Density
mhe = he_alt.*den_alt;% Helium Mass Density g/cm^3 
mN2 = N2_alt.*den_alt;% N2 Mass Density
mO1 = O1_alt.*den_alt;% Oxygen Mass Density


% r=6372;%Earth Radius
% g_tp=9.8.*r^2./(r+geom_alt).^2; %Find g as altitude changes [m/s]
k=1.38e-23;% Boltzman's Constant
atom_unit=1.67e-27;%kg/unit
Av = 6.022141*10^23;
kmol=6.02214*10^26;
mmw_he=0.004; %Helium atomic mass
mmw_N2=0.02801; %N2 molecular mass
mmw_O1=0.016; % O1 molecular mass
points = length(tn_alt);
meanmass = mbar_alt/1000;  % [kg/mol]


%-----Finding Pressure Scale Heights for gas constituents from kT/mg-----
Hp_he=k.*tn_alt./(mmw_he/Av.*g0)/1000;   % in Km
Hp_mean=k.*tn_alt./(meanmass/Av.*g0)/1000;
Hp_n2=k.*tn_alt./(mmw_N2/Av.*g0)/1000;
Hp_o1=k.*tn_alt./(mmw_O1/Av.*g0)/1000;

% -----Calculate Scale Heights Using 3 Point Differentiation Technique-----
H_He_star = zeros(points,1);
H_dentot = zeros(points,1);
H_n2 = zeros(points,1);
H_o1 = zeros(points,1);
H_tn = zeros(points,1);
H_mass = zeros(points,1);


for i=1:points
    if i==1 %First Point gradient technique
        coeff1 = (2*geop_alt(1)-geop_alt(2)-geop_alt(3))/((geop_alt(1)-...
            geop_alt(2))*(geop_alt(1)-geop_alt(3)));
        
        coeff2 = (2*geop_alt(1)-geop_alt(1)-geop_alt(3))/((geop_alt(2)-...
            geop_alt(1))*(geop_alt(2)-geop_alt(3)));
        
        coeff3 = (2*geop_alt(1)-geop_alt(1)-geop_alt(2))/((geop_alt(3)-...
            geop_alt(1))*(geop_alt(3)-geop_alt(2)));
        
        H_He_star(1) = -1/mhe(1)*(mhe(1)*coeff1+mhe(2)*coeff2+...
            mhe(3)*coeff3);
        H_dentot(1) = -1/den_alt(1)*(den_alt(1)*coeff1+den_alt(2)*coeff2+...
            den_alt(3)*coeff3);
        H_n2(1) = -1/mN2(1)*(mN2(1)*coeff1+mN2(2)*coeff2+...
            mN2(3)*coeff3);
        H_o1(1) = -1/mO1(1)*(mO1(1)*coeff1+mO1(2)*coeff2+...
            mO1(3)*coeff3);
        H_tn(1) = 1/tn_alt(1)*(tn_alt(1)*coeff1+tn_alt(2)*coeff2+...
            tn_alt(3)*coeff3);
        H_mass(1) = -1/mbar_alt(1)*(mbar_alt(1)*coeff1+mbar_alt(2)*coeff2+...
            mbar_alt(3)*coeff3);
        
        
    elseif i==points %Last point gradient technique
        coeff1 = (2*geop_alt(i)-geop_alt(i-1)-geop_alt(i))/((geop_alt(i-2)-...
            geop_alt(i-1))*(geop_alt(i-2)-geop_alt(i)));
        
        coeff2 = (2*geop_alt(i)-geop_alt(i-2)-geop_alt(i))/((geop_alt(i-1)-...
            geop_alt(i-2))*(geop_alt(i-1)-geop_alt(i)));
        
        coeff3 = (2*geop_alt(i)-geop_alt(i-2)-geop_alt(i-1))/((geop_alt(i)-...
            geop_alt(i-2))*(geop_alt(i)-geop_alt(i-1)));
        
        H_He_star(i) = -1/mhe(i)*(mhe(i-2)*coeff1+mhe(i-1)*coeff2+...
            mhe(i)*coeff3);
        H_dentot(i) = -1/den_alt(i)*(den_alt(i-2)*coeff1+den_alt(i-1)*coeff2+...
            den_alt(i)*coeff3);
        H_n2(i) = -1/mN2(i)*(mN2(i-2)*coeff1+mN2(i-1)*coeff2+...
            mN2(i)*coeff3);
        H_o1(i) = -1/mO1(i)*(mO1(i-2)*coeff1+mO1(i-1)*coeff2+...
            mO1(i)*coeff3);
        H_tn(i) = 1/tn_alt(i)*(tn_alt(i-2)*coeff1+tn_alt(i-1)*coeff2+...
            tn_alt(i)*coeff3);
        H_mass(i) = -1/mbar_alt(i)*(mbar_alt(i-2)*coeff1+mbar_alt(i-1)*coeff2+...
            mbar_alt(i)*coeff3);

        
    else %Middle Points gradient technique
        coeff1 = (2*geop_alt(i)-geop_alt(i)-geop_alt(i+1))/((geop_alt(i-1)-...
            geop_alt(i))*(geop_alt(i-1)-geop_alt(i+1)));
        
        coeff2 = (2*geop_alt(i)-geop_alt(i-1)-geop_alt(i+1))/((geop_alt(i)-...
            geop_alt(i-1))*(geop_alt(i)-geop_alt(i+1)));
        
        coeff3 = (2*geop_alt(i)-geop_alt(i-1)-geop_alt(i))/((geop_alt(i+1)-...
            geop_alt(i-1))*(geop_alt(i+1)-geop_alt(i)));
        
        H_He_star(i) = -1/mhe(i)*(mhe(i-1)*coeff1+mhe(i)*coeff2+...
            mhe(i+1)*coeff3);
        H_dentot(i) = -1/den_alt(i)*(den_alt(i-1)*coeff1+den_alt(i)*coeff2+...
            den_alt(i+1)*coeff3);
        H_n2(i) = -1/mN2(i)*(mN2(i-1)*coeff1+mN2(i)*coeff2+...
            mN2(i+1)*coeff3);
        H_o1(i) = -1/mO1(i)*(mO1(i-1)*coeff1+mO1(i)*coeff2+...
            mO1(i+1)*coeff3);
        H_tn(i) = 1/tn_alt(i)*(tn_alt(i-1)*coeff1+tn_alt(i)*coeff2+...
            tn_alt(i+1)*coeff3);
        H_mass(i) = -1/mbar_alt(i)*(mbar_alt(i-1)*coeff1+mbar_alt(i)*coeff2+...
            mbar_alt(i+1)*coeff3);

    end
end
%-----Get Scale Height from Inverse-----
H_He_star = 1./H_He_star;
H_tot_star = 1./H_dentot;
H_N2_star = 1./H_n2;
H_O1_star = 1./H_o1(1:end);
H_temp = (1./H_tn).';
H_temp_he = H_temp/.62;% <---- Alpha for helium is -.38
H_mass = (1./H_mass).';

H_He_star = H_He_star.';
H_tot_star = H_tot_star.';
H_N2_star = H_N2_star.';
H_O1_star = H_O1_star.';


if linear==1
    %--------Linear Gradient Calculations--------
    % -----Helium Scale Height-----
    lnn_bf = reallog(nhe_bf);
    H_he_star_r_bf = zeros(points,1);
    for i=1:points
        if i==1
            H_he_star_r_bf(i)=-(lnn_bf(2)-lnn_bf(1))./...
                            (geop_alt(2)-geop_alt(1));
        elseif i==points
            H_he_star_r_bf(i)=-(lnn_bf(i)-lnn_bf(i-1))./...
                            (geop_alt(i)-geop_alt(i-1));
        else
            H_he_star_r_bf(i)=-(lnn_bf(i+1)-lnn_bf(i-1))./...
                            (geop_alt(i+1)-geop_alt(i-1));
        end
    end

    H_He_star=1./H_he_star_r_bf;
    H_He_star = H_He_star.';

    % -----Scale Height for Total Density-----
    lnnden = reallog(den_alt);
    H_dentot = zeros(points,1);
    for i=1:points
        if i==1
            H_dentot(1) = -(lnnden(2)-lnnden(1))/(geop_alt(2)-...
                geop_alt(1));
        elseif i==points
            H_dentot(i) = -(lnnden(i)-lnnden(i-1))/(geop_alt(i)-...
                geop_alt(i-1));
        else
        H_dentot(i) = -(.5*(lnnden(i)-lnnden(i-1))/(geop_alt(i)-...
            geop_alt(i-1))+.5*(lnnden(i+1)-lnnden(i))/(geop_alt(i+1)-...
            geop_alt(i)));
        end
    end
    H_den_tot = 1./H_dentot;

    % -----Scale Height for N2-----
    lnnN2 = reallog(mN2);
    H_n2 = zeros(points,1);
    for i=1:points
        if i==1
            H_n2(1) = -(lnnN2(2)-lnnN2(1))/(geop_alt(2)-...
                geop_alt(1));
        elseif i==points
            H_n2(i) = -(lnnN2(i)-lnnN2(i-1))/(geop_alt(i)-...
                geop_alt(i-1));
        else
        H_n2(i) = -(lnnN2(i+1)-lnnN2(i-1))/(geop_alt(i+1)-...
            geop_alt(i-1));
        end
    end
    H_N2_star = 1./H_n2;
    H_N2_star = H_N2_star.';

    % -----Scale Height for O1-----
    lnO1 = reallog(nO1);
    H_o1 = zeros(points,1);
    for i=1:points
        if i==1
            H_o1(1) = -(lnO1(2)-lnO1(1))/(geop_alt(2)-...
                geop_alt(1));
        elseif i==points
            H_o1(i) = -(lnO1(i)-lnO1(i-1))/(geop_alt(i)-...
                geop_alt(i-1));
        else
        H_o1(i) = -(lnO1(i+1)-lnO1(i-1))/(geop_alt(i+1)-...
            geop_alt(i-1));
        end
    end
    H_O1_star = 1./H_o1(4:end);     % changed from 4:end to 1:end ? .... HLH
end


%----Put Together Diffusive Profiles and Mean Mass Profile-----
H_He_diff = 1./(1./H_temp_he+1./Hp_he);% Helium diffusive profile
H_He_diff_no_alpha = 1./(1./H_temp+1./Hp_he);
H_tot_diff = 1./(1./H_temp+1./Hp_mean+1./H_mass);% Atmospheric scale height
H_N2_diff = 1./(1./H_temp+1./Hp_n2);% N2 Diffusive profile
H_O1_diff = 1./(1./H_temp+1./Hp_o1);% O1 Diffusive profile

end

