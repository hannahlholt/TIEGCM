function [D_i] = Di_TIEGCM(species, T, p)
% Calculates the molecular diffusion coefficient using TIEGCM values

% T = input array of tempurates [K] vs. altitude at specific lat, lon, and UT
% P = input array of pressures [Pa] vs. altitude at specific lat, lon, and UT
% D_i = output array of species diffusion coefficient in [m^2/s]

T0 = 273;       % reference Temp [K]
p0 = 1E5;       % reference pressure [Pa]

switch species
    case 'N2'
        % [O2-N2, O-N2, He-N2] 
        a = [0.18; 0.26; 0.612661];       %[cm^2/s]  
        s = [1.75; 1.75; 1.718];
    case 'O2'
        % [O2-N2, O2-O, O2-He]
        a = [0.18; 0.26; 0.648966];
        s = [1.75; 1.75; 1.71];
    case 'O1'
        % [O-N2, O-O2, O-He]
        a = [0.26; 0.26; 0.865538];
        s = [1.75; 1.75; 1.749]; 
    case 'He'
        % [He-N2, He-O2, He-O]
        a = [0.621661; 0.648966; 0.865539];
        s = [1.718; 1.710; 1.749];
end

D_i = 0;

for i=1:3
    D_i = D_i + ( a(i) .* (T./T0).^s(i) .* (p0./p) ) ./ 1E4; % [m^2/s]
end

end

