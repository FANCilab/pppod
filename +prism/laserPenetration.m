function trasmission = laserPenetration(waveL, z)

g = 0.9; %brain scattering anisotropy
a = 1.1; %scaling factor mm^-1
b= 1.37; %average brain scattering power;

scattering = (a.*(waveL/500)./(1-g)).^(-b);

trasmission = exp(-scattering*z);

for iw = 1: numel(waveL)
    trasmission(:, iw) = exp(-scattering(iw)*z);
end

% ww = [780 890 970 1020];
% 
% Lw = [0.97 1.3 1.8 2];
% 
% for iw = 1: numel(waveL)
% trasmission(:, iw) = exp(-z./(Lw(waveL(iw) == ww)*50));
% end

end