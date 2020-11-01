function [flow,propri,results] = propriIni(param,flow,data)

% Inputs:
% data: Contient les paramètres associés au système
% param: Contient les paramètres associés au tuyau ( géométrie, etc.)
% flow: Contient les paramètres associés à l'écoulement


% Outputs:
% flow: Contient les paramètres associés à l'écoulement
% propri : Contient les propriétés thermodynamiques associées à l'écoulement
% results --> Initialisation
    % epsilon_num:  [-]     Taux de vide
    % P_num:        [kPa]   Pression du réfrigérant
    % hm_num:       [J/kg]  Enthalpie du réfrigérant
    % Tf_numk:      [K]     Champs de température du réfrigérant
    % x_num:        [-]     Titre de l'écoulement du réfrigérant
    % xth_num:      [-]     Titre thermodynamique du réfrigérant
    % Vz_num:       [m/s]   Vitesse du réfrigérant

%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯%
%----------------------------- Ref Gas -----------------------------------%
%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% RefGas  flow parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                   % [J/kg] Enthalpie massique du réfrigérant à l'entrée
flow.Tf_ink = refpropm('T','P',flow.Pin,'H',flow.hf_in,data.refrigerant);   % [K] Température du réfrigérant à l'entrée
flow.hyp_dp = -param.L;                            % [kPa] hypothèse sur la perte de charge dans la conduite (basée sur la longueur de la conduite)
flow.Pout = flow.Pin + flow.hyp_dp;         % [kPa] Presion du réfrigérant à la sortie 
flow.hyp_dt = (flow.Tf_ink-param.Tinfk)/2;  % [K] hypothèse sur la différence de température entre la paroie extérieure de la conduite et le fluide utilisée dans le réservoir de stockage
flow.Ts_in_k = param.Tinfk + flow.hyp_dt;       % [K] Température de la surface intérieure de la conduite
flow.Ts_out_k = param.Tinfk + flow.hyp_dt-1;    % [K] Température de la surface extérieure de la conduite
flow.G = flow.m_in/param.As;                         % [kg/(s m2)] Flux massique du réfrigérant
flow.rho_in = refpropm('D','P',flow.Pin,'H',flow.hf_in,data.refrigerant);   % [kg/m3] Masse volumique du réfrigérant
flow.Vz_in = flow.G/flow.rho_in;                                                % [m/s] Vitesse d'entrée du réfrigérant
flow.Pcr = refpropm('P','C',0,' ',0,data.refrigerant);                       % [kPa] Pression critique du réfrigérant
flow.hcr = refpropm('H','C',0,' ',0,data.refrigerant);                       % [J/kg] Enthalpie critique du réfrigérant
if flow.Pin < flow.Pcr % L'écoulement est sous-critique
    flow.hLL_in = refpropm('H','P',flow.Pin,'Q',0,data.refrigerant);           % [J/kg] Enthalpie massique du réfrigérant HP à saturation liquide
    flow.rho_LL_in = refpropm('D','P',flow.Pin,'Q',0,data.refrigerant);        % [kg/m3] Masse volumique du réfrigérant HP à saturation liquide
    flow.hGG_in = refpropm('H','P',flow.Pin,'Q',1,data.refrigerant);           % [J/kg] Enthalpie massique du réfrigérant HP à saturation gazeuse
    flow.rho_GG_in = refpropm('D','P',flow.Pin,'Q',1,data.refrigerant);        % [kg/m3] Masse volumique du réfrigérant HP à saturation liquide
    flow.xth_in = (flow.hf_in-flow.hLL_in)/(flow.hGG_in-flow.hLL_in);       % [-] Titre thermodynamique d'entrée du réfrigérant HP
    if  flow.xth_in > 1 % L'écoulement est gazeux
        flow.x_in = 1;             % [-] Titre de l'écoulement d'entrée du réfrigérant HP
    elseif flow.xth_in < 0 % L'écoulement est liquide
        flow.x_in = 0;             % [-] Titre de l'écoulement d'entrée du réfrigérant HP
    else % L'écoulement est à l'état diphasique
        flow.x_in =flow.xth_in;  % [-] Titre de l'écoulement d'entrée du réfrigérant HP
    end
else % L'écoulement est super-critique
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% RefGas Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  flow.Pin >= flow.Pcr % L'écoulement est super-critique
    if flow.hf_in >= flow.hcr
        results.xth_num = ones(param.nb_z,2);             % [-] Initialisation du champs de titre thermodynamique du réfrigérant
        results.x_num = ones(param.nb_z,2);               % [-] Initialisation du champs de titre de l'écoulement du réfrigérant
        results.x_star = zeros(param.nb_z,1) + eps;       % [-] Initialisation de l'enthalpie de sous-refroidissement
        results.epsilon_num = ones(param.nb_z,2) - eps;   % [-] Initialisation du champs de taux de vide du réfrigérant
    else
        results.xth_num = zeros(param.nb_z,2);                   % [-] Initialisation du champs de titre thermodynamique du réfrigérant
        results.x_num = zeros(param.nb_z,2);                     % [-] Initialisation du champs de titre de l'écoulement du réfrigérant
        results.x_star = zeros(param.nb_z,1) + eps;              % [-] Initialisation de l'enthalpie de sous-refroidissement
        results.epsilon_num = zeros(param.nb_z,2) + eps;         % [-] Initialisation du champs de taux de vide du réfrigérant
    end
else % L'écoulement est sous-critique
    if flow.xth_in >= 1 % L'écoulement est gazeux
        results.xth_num = ones(param.nb_z,2) + flow.xth_in;  % [-] Initialisation du champs de titre thermodynamique du réfrigérant
        results.x_num = ones(param.nb_z,2);                    % [-] Initialisation du champs de titre de l'écoulement du réfrigérant
        results.x_star = zeros(param.nb_z,1) + eps;            % [-] Initialisation de l'enthalpie de sous-refroidissement
        results.epsilon_num = ones(param.nb_z,2) - eps;        % [-] Initialisation du champs de taux de vide du réfrigérant
    elseif flow.xth_in <= 0 % L'écoulement est liquide
        results.xth_num = zeros(param.nb_z,2) + flow.xth_in; % [-] Initialisation du champs de titre thermodynamique du réfrigérant
        results.x_num = zeros(param.nb_z,2);                   % [-] Initialisation du champs de titre de l'écoulement du réfrigérant
        results.x_star = zeros(param.nb_z,1) + eps;            % [-] Initialisation de l'enthalpie de sous-refroidissement
        results.epsilon_num = zeros(param.nb_z,2) + eps;       % [-] Initialisation du champs de taux de vide du réfrigérant
    else % L,écoulement est diphasique
        results.xth_num = zeros(param.nb_z,2) + flow.xth_in;  % [-] Initialisation du champs de titre thermodynamique du réfrigérant
        results.x_num = zeros(param.nb_z,2) + flow.x_in;      % [-] Initialisation du champs de titre de l'écoulement du réfrigérant
        results.x_star = zeros(param.nb_z,1) + eps;             % [-] Initialisation de l'enthalpie de sous-refroidissement
        results.epsilon_num = zeros(param.nb_z,2) + flow.x_in*flow.rho_LL_in/(flow.x_in*flow.rho_LL_in+(1-flow.x_in)*flow.rho_GG_in); % [-] Initialisation du champs de taux de vide du réfrigérant
    end
end
results.P_num = [linspace(flow.Pin,flow.Pout,param.nb_z)',linspace(flow.Pin,flow.Pout,param.nb_z)']; % [kPa] Initialisation de la pression du réfrigérant
results.hm_num = zeros(param.nb_z,2) + flow.hf_in;      % [J/kg] Initialisation de l'enthalpie du réfrigérant
results.Vz_num = zeros(param.nb_z,2) + flow.Vz_in;      % [m/s] Initialisation de la vitesse du réfrigérant
if flow.Pin >= flow.Pcr || flow.xth_in > 1 || flow.xth_in < 0 % L'écoulement est gazeux, liquide ou super-critique
    results.Tf_numk = [linspace(flow.Tf_ink,param.Tinfk,param.nb_z)',linspace(flow.Tf_ink,param.Tinfk,param.nb_z)'];  % [K] Initialisation de la température du réfrigérant
else % L'écoulement est diphasique
    results.Tf_numk = zeros(param.nb_z,3)+ flow.Tf_ink; % [K] Initialisation de la température du réfrigérant HP
end
results.rho_m = zeros(param.nb_z,1) + flow.rho_in;           % [kg/m3] Masse volumique du mélange
results.qp = zeros(param.nb_z,1);                                   % [W] Initialisation du flux de chaleur
results.Ts_in_numk = zeros(param.nb_z,1) + flow.Ts_in_k;     % [K] Initialisation de la température de paroie intérieure
results.Ts_out_numk = zeros(param.nb_z,1) + flow.Ts_out_k;  % [K] Initialisation de la température de paroie extérieure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Propriétés initiales du réfrigérant RefGas  %%%%%%%%%%%%%%%%%%%%%%%
% Déclaration des vecteurs associées aux propriétés du réfrigérant 
propri.rho_L = zeros(param.nb_z,1);
propri.rho_G = zeros(param.nb_z,1);
propri.mu_L = zeros(param.nb_z,1);
propri.mu_G = zeros(param.nb_z,1);
propri.h_L = zeros(param.nb_z,1);
propri.h_G = zeros(param.nb_z,1);
propri.cp_L = zeros(param.nb_z,1);
propri.cp_G = zeros(param.nb_z,1);
propri.h_LL = zeros(param.nb_z,1);  % à saturation liquide
propri.h_GG = zeros(param.nb_z,1);  % à saturation gazeuse
propri.Tsat = zeros(param.nb_z,1);
propri.Pr_L = zeros(param.nb_z,1);
propri.Pr_G = zeros(param.nb_z,1);
propri.k_L = zeros(param.nb_z,1);
propri.k_G = zeros(param.nb_z,1);
propri.sigma_L = zeros(param.nb_z,1);
propri.Psat_wall = zeros(param.nb_z,1);

for k=1:param.nb_z
    if  flow.Pin < flow.Pcr % L'écoulement est sous-critique
        [propri.h_LL(k),propri.Tsat(k)] = refpropm('HT','P',results.P_num(k,2),'Q',0,data.refrigerant);   % [J/kg] Enthalpie massique du liquide saturée; [K] Saturation temperature
        propri.h_GG(k) = refpropm('H','P',results.P_num(k,2),'Q',1,data.refrigerant);   % [J/kg] Enthalpie massique du vapeur saturée
        if  results.x_num(k,2) == 0  % L'écoulement est liquide
            [propri.rho_L(k),propri.mu_L(k),propri.h_L(k),propri.cp_L(k),propri.Pr_L(k),propri.k_L(k)] = refpropm('DVHC^L','T',results.Tf_numk(k,2),'P',results.P_num(k,2),data.refrigerant);   % liquide
            [propri.rho_G(k),propri.mu_G(k),propri.h_G(k),propri.cp_G(k),propri.Pr_G(k),propri.k_G(k)] = refpropm('DVHC^L','P',results.P_num(k,2),'Q',1,data.refrigerant); % vapeur saturée
        elseif results.x_num(k,2) == 1  % L'écoulement est gazeux
            [propri.rho_L(k),propri.mu_L(k),propri.h_L(k),propri.cp_L(k),propri.Pr_L(k),propri.k_L(k)] = refpropm('DVHC^L','P',results.P_num(k,2),'Q',0,data.refrigerant);  % liquide saturée
            [propri.rho_G(k),propri.mu_G(k),propri.h_G(k),propri.cp_G(k),propri.Pr_G(k),propri.k_G(k)] = refpropm('DVHC^L','T',results.Tf_numk(k,2),'P',results.P_num(k,2),data.refrigerant); % vapeur
        else % L'écoulement est diphasique
            [propri.rho_L(k),propri.mu_L(k),propri.h_L(k),propri.cp_L(k),propri.Pr_L(k),propri.k_L(k),propri.sigma_L(k)] = refpropm('DVHC^LI','P',results.P_num(k,2),'Q',0,data.refrigerant); % liquide saturée
            [propri.rho_G(k),propri.mu_G(k),propri.h_G(k),propri.cp_G(k),propri.Pr_G(k),propri.k_G(k)] = refpropm('DVHC^L','P',results.P_num(k,2),'Q',1,data.refrigerant); % vapeur saturée
            propri.Psat_wall(k) = refpropm('P','T',results.Ts_in_numk(k),'Q',0,data.refrigerant); % [kPa] Saturation pressure at Wall temperature
        end
    else % L'écoulement est super-critique
        [propri.rho_L(k),propri.mu_L(k),propri.h_L(k),propri.cp_L(k),propri.Pr_L(k),propri.k_L(k)] = refpropm('DVHC^L','T',results.Tf_numk(k,2),'P',results.P_num(k,2),data.refrigerant);   % liquide
        [propri.rho_G(k),propri.mu_G(k),propri.h_G(k),propri.cp_G(k),propri.Pr_G(k),propri.k_G(k)] = refpropm('DVHC^L','T',results.Tf_numk(k,2),'P',results.P_num(k,2),data.refrigerant); % vapeur
    end
    %                 fprintf('%c%c%c%c%03d%%',8,8,8,8,round(100*(k/FlowSetup.nb_z))); % Progression 0 --> 100 %
end
