function dTsub = SahaZuber(qppw,G,D,cpL,kL)
%-------------------------------------------------------------------------%
% dTsub = SahaZuber(qppw,G,P,D,fluid)
% Subcooled temperature difference
% based on Saha and Zuber correlation
%
% Inputs
% qppw: [W/m2] Wall heat flux 
% G: [kg/m2s] Total mass flux
% D: [m] Pipe inside diameter
%-------------Properties--------------------------------------------------%
% cpL: [J/(kg K)] Heat capacity
% kL: [W/(m K)] Thermal conductivity %

% Output
% dTsub: [K] subcooled temperature difference
%-------------------------------------------------------------------------%

%--------------------------------------------------------------------------
Pe = G*D.*cpL./kL;      % [-] Peclet number
% if Pe<=70000 Contrôlé par la diffusion
dTsub(Pe<=70000,1) = 0.0022*qppw(Pe<=70000)*D./kL(Pe<=70000); % [K] Tsat-Tbulk_liquid
% Pe>70000 --> Contrôlé par l'advection   
dTsub(Pe>70000,1) = 154*qppw(Pe>70000)./(G*cpL(Pe>70000)); % [K] Tsat-Tbulk_liquid

end






