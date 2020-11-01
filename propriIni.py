function [flow,propri,results] = propriIni(param,flow,data)

% Inputs:
% data: Contient les param�tres associ�s au syst�me
% param: Contient les param�tres associ�s au tuyau ( g�om�trie, etc.)
% flow: Contient les param�tres associ�s � l'�coulement


% Outputs:
% flow: Contient les param�tres associ�s � l'�coulement
% propri : Contient les propri�t�s thermodynamiques associ�es � l'�coulement
% results --> Initialisation
    % epsilon_num:  [-]     Taux de vide
    % P_num:        [kPa]   Pression du r�frig�rant
    % hm_num:       [J/kg]  Enthalpie du r�frig�rant
    % Tf_numk:      [K]     Champs de temp�rature du r�frig�rant
    % x_num:        [-]     Titre de l'�coulement du r�frig�rant
    % xth_num:      [-]     Titre thermodynamique du r�frig�rant
    % Vz_num:       [m/s]   Vitesse du r�frig�rant

%�������������������������������������������������������������������������%
%----------------------------- Ref Gas -----------------------------------%
%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% RefGas  flow parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                   % [J/kg] Enthalpie massique du r�frig�rant � l'entr�e
flow.Tf_ink = refpropm('T','P',flow.Pin,'H',flow.hf_in,data.refrigerant);   % [K] Temp�rature du r�frig�rant � l'entr�e
flow.hyp_dp = -param.L;                            % [kPa] hypoth�se sur la perte de charge dans la conduite (bas�e sur la longueur de la conduite)
flow.Pout = flow.Pin + flow.hyp_dp;         % [kPa] Presion du r�frig�rant � la sortie 
flow.hyp_dt = (flow.Tf_ink-param.Tinfk)/2;  % [K] hypoth�se sur la diff�rence de temp�rature entre la paroie ext�rieure de la conduite et le fluide utilis�e dans le r�servoir de stockage
flow.Ts_in_k = param.Tinfk + flow.hyp_dt;       % [K] Temp�rature de la surface int�rieure de la conduite
flow.Ts_out_k = param.Tinfk + flow.hyp_dt-1;    % [K] Temp�rature de la surface ext�rieure de la conduite
flow.G = flow.m_in/param.As;                         % [kg/(s m2)] Flux massique du r�frig�rant
flow.rho_in = refpropm('D','P',flow.Pin,'H',flow.hf_in,data.refrigerant);   % [kg/m3] Masse volumique du r�frig�rant
flow.Vz_in = flow.G/flow.rho_in;                                                % [m/s] Vitesse d'entr�e du r�frig�rant
flow.Pcr = refpropm('P','C',0,' ',0,data.refrigerant);                       % [kPa] Pression critique du r�frig�rant
flow.hcr = refpropm('H','C',0,' ',0,data.refrigerant);                       % [J/kg] Enthalpie critique du r�frig�rant
if flow.Pin < flow.Pcr % L'�coulement est sous-critique
    flow.hLL_in = refpropm('H','P',flow.Pin,'Q',0,data.refrigerant);           % [J/kg] Enthalpie massique du r�frig�rant HP � saturation liquide
    flow.rho_LL_in = refpropm('D','P',flow.Pin,'Q',0,data.refrigerant);        % [kg/m3] Masse volumique du r�frig�rant HP � saturation liquide
    flow.hGG_in = refpropm('H','P',flow.Pin,'Q',1,data.refrigerant);           % [J/kg] Enthalpie massique du r�frig�rant HP � saturation gazeuse
    flow.rho_GG_in = refpropm('D','P',flow.Pin,'Q',1,data.refrigerant);        % [kg/m3] Masse volumique du r�frig�rant HP � saturation liquide
    flow.xth_in = (flow.hf_in-flow.hLL_in)/(flow.hGG_in-flow.hLL_in);       % [-] Titre thermodynamique d'entr�e du r�frig�rant HP
    if  flow.xth_in > 1 % L'�coulement est gazeux
        flow.x_in = 1;             % [-] Titre de l'�coulement d'entr�e du r�frig�rant HP
    elseif flow.xth_in < 0 % L'�coulement est liquide
        flow.x_in = 0;             % [-] Titre de l'�coulement d'entr�e du r�frig�rant HP
    else % L'�coulement est � l'�tat diphasique
        flow.x_in =flow.xth_in;  % [-] Titre de l'�coulement d'entr�e du r�frig�rant HP
    end
else % L'�coulement est super-critique
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% RefGas Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  flow.Pin >= flow.Pcr % L'�coulement est super-critique
    if flow.hf_in >= flow.hcr
        results.xth_num = ones(param.nb_z,2);             % [-] Initialisation du champs de titre thermodynamique du r�frig�rant
        results.x_num = ones(param.nb_z,2);               % [-] Initialisation du champs de titre de l'�coulement du r�frig�rant
        results.x_star = zeros(param.nb_z,1) + eps;       % [-] Initialisation de l'enthalpie de sous-refroidissement
        results.epsilon_num = ones(param.nb_z,2) - eps;   % [-] Initialisation du champs de taux de vide du r�frig�rant
    else
        results.xth_num = zeros(param.nb_z,2);                   % [-] Initialisation du champs de titre thermodynamique du r�frig�rant
        results.x_num = zeros(param.nb_z,2);                     % [-] Initialisation du champs de titre de l'�coulement du r�frig�rant
        results.x_star = zeros(param.nb_z,1) + eps;              % [-] Initialisation de l'enthalpie de sous-refroidissement
        results.epsilon_num = zeros(param.nb_z,2) + eps;         % [-] Initialisation du champs de taux de vide du r�frig�rant
    end
else % L'�coulement est sous-critique
    if flow.xth_in >= 1 % L'�coulement est gazeux
        results.xth_num = ones(param.nb_z,2) + flow.xth_in;  % [-] Initialisation du champs de titre thermodynamique du r�frig�rant
        results.x_num = ones(param.nb_z,2);                    % [-] Initialisation du champs de titre de l'�coulement du r�frig�rant
        results.x_star = zeros(param.nb_z,1) + eps;            % [-] Initialisation de l'enthalpie de sous-refroidissement
        results.epsilon_num = ones(param.nb_z,2) - eps;        % [-] Initialisation du champs de taux de vide du r�frig�rant
    elseif flow.xth_in <= 0 % L'�coulement est liquide
        results.xth_num = zeros(param.nb_z,2) + flow.xth_in; % [-] Initialisation du champs de titre thermodynamique du r�frig�rant
        results.x_num = zeros(param.nb_z,2);                   % [-] Initialisation du champs de titre de l'�coulement du r�frig�rant
        results.x_star = zeros(param.nb_z,1) + eps;            % [-] Initialisation de l'enthalpie de sous-refroidissement
        results.epsilon_num = zeros(param.nb_z,2) + eps;       % [-] Initialisation du champs de taux de vide du r�frig�rant
    else % L,�coulement est diphasique
        results.xth_num = zeros(param.nb_z,2) + flow.xth_in;  % [-] Initialisation du champs de titre thermodynamique du r�frig�rant
        results.x_num = zeros(param.nb_z,2) + flow.x_in;      % [-] Initialisation du champs de titre de l'�coulement du r�frig�rant
        results.x_star = zeros(param.nb_z,1) + eps;             % [-] Initialisation de l'enthalpie de sous-refroidissement
        results.epsilon_num = zeros(param.nb_z,2) + flow.x_in*flow.rho_LL_in/(flow.x_in*flow.rho_LL_in+(1-flow.x_in)*flow.rho_GG_in); % [-] Initialisation du champs de taux de vide du r�frig�rant
    end
end
results.P_num = [linspace(flow.Pin,flow.Pout,param.nb_z)',linspace(flow.Pin,flow.Pout,param.nb_z)']; % [kPa] Initialisation de la pression du r�frig�rant
results.hm_num = zeros(param.nb_z,2) + flow.hf_in;      % [J/kg] Initialisation de l'enthalpie du r�frig�rant
results.Vz_num = zeros(param.nb_z,2) + flow.Vz_in;      % [m/s] Initialisation de la vitesse du r�frig�rant
if flow.Pin >= flow.Pcr || flow.xth_in > 1 || flow.xth_in < 0 % L'�coulement est gazeux, liquide ou super-critique
    results.Tf_numk = [linspace(flow.Tf_ink,param.Tinfk,param.nb_z)',linspace(flow.Tf_ink,param.Tinfk,param.nb_z)'];  % [K] Initialisation de la temp�rature du r�frig�rant
else % L'�coulement est diphasique
    results.Tf_numk = zeros(param.nb_z,3)+ flow.Tf_ink; % [K] Initialisation de la temp�rature du r�frig�rant HP
end
results.rho_m = zeros(param.nb_z,1) + flow.rho_in;           % [kg/m3] Masse volumique du m�lange
results.qp = zeros(param.nb_z,1);                                   % [W] Initialisation du flux de chaleur
results.Ts_in_numk = zeros(param.nb_z,1) + flow.Ts_in_k;     % [K] Initialisation de la temp�rature de paroie int�rieure
results.Ts_out_numk = zeros(param.nb_z,1) + flow.Ts_out_k;  % [K] Initialisation de la temp�rature de paroie ext�rieure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Propri�t�s initiales du r�frig�rant RefGas  %%%%%%%%%%%%%%%%%%%%%%%
% D�claration des vecteurs associ�es aux propri�t�s du r�frig�rant 
propri.rho_L = zeros(param.nb_z,1);
propri.rho_G = zeros(param.nb_z,1);
propri.mu_L = zeros(param.nb_z,1);
propri.mu_G = zeros(param.nb_z,1);
propri.h_L = zeros(param.nb_z,1);
propri.h_G = zeros(param.nb_z,1);
propri.cp_L = zeros(param.nb_z,1);
propri.cp_G = zeros(param.nb_z,1);
propri.h_LL = zeros(param.nb_z,1);  % � saturation liquide
propri.h_GG = zeros(param.nb_z,1);  % � saturation gazeuse
propri.Tsat = zeros(param.nb_z,1);
propri.Pr_L = zeros(param.nb_z,1);
propri.Pr_G = zeros(param.nb_z,1);
propri.k_L = zeros(param.nb_z,1);
propri.k_G = zeros(param.nb_z,1);
propri.sigma_L = zeros(param.nb_z,1);
propri.Psat_wall = zeros(param.nb_z,1);

for k=1:param.nb_z
    if  flow.Pin < flow.Pcr % L'�coulement est sous-critique
        [propri.h_LL(k),propri.Tsat(k)] = refpropm('HT','P',results.P_num(k,2),'Q',0,data.refrigerant);   % [J/kg] Enthalpie massique du liquide satur�e; [K] Saturation temperature
        propri.h_GG(k) = refpropm('H','P',results.P_num(k,2),'Q',1,data.refrigerant);   % [J/kg] Enthalpie massique du vapeur satur�e
        if  results.x_num(k,2) == 0  % L'�coulement est liquide
            [propri.rho_L(k),propri.mu_L(k),propri.h_L(k),propri.cp_L(k),propri.Pr_L(k),propri.k_L(k)] = refpropm('DVHC^L','T',results.Tf_numk(k,2),'P',results.P_num(k,2),data.refrigerant);   % liquide
            [propri.rho_G(k),propri.mu_G(k),propri.h_G(k),propri.cp_G(k),propri.Pr_G(k),propri.k_G(k)] = refpropm('DVHC^L','P',results.P_num(k,2),'Q',1,data.refrigerant); % vapeur satur�e
        elseif results.x_num(k,2) == 1  % L'�coulement est gazeux
            [propri.rho_L(k),propri.mu_L(k),propri.h_L(k),propri.cp_L(k),propri.Pr_L(k),propri.k_L(k)] = refpropm('DVHC^L','P',results.P_num(k,2),'Q',0,data.refrigerant);  % liquide satur�e
            [propri.rho_G(k),propri.mu_G(k),propri.h_G(k),propri.cp_G(k),propri.Pr_G(k),propri.k_G(k)] = refpropm('DVHC^L','T',results.Tf_numk(k,2),'P',results.P_num(k,2),data.refrigerant); % vapeur
        else % L'�coulement est diphasique
            [propri.rho_L(k),propri.mu_L(k),propri.h_L(k),propri.cp_L(k),propri.Pr_L(k),propri.k_L(k),propri.sigma_L(k)] = refpropm('DVHC^LI','P',results.P_num(k,2),'Q',0,data.refrigerant); % liquide satur�e
            [propri.rho_G(k),propri.mu_G(k),propri.h_G(k),propri.cp_G(k),propri.Pr_G(k),propri.k_G(k)] = refpropm('DVHC^L','P',results.P_num(k,2),'Q',1,data.refrigerant); % vapeur satur�e
            propri.Psat_wall(k) = refpropm('P','T',results.Ts_in_numk(k),'Q',0,data.refrigerant); % [kPa] Saturation pressure at Wall temperature
        end
    else % L'�coulement est super-critique
        [propri.rho_L(k),propri.mu_L(k),propri.h_L(k),propri.cp_L(k),propri.Pr_L(k),propri.k_L(k)] = refpropm('DVHC^L','T',results.Tf_numk(k,2),'P',results.P_num(k,2),data.refrigerant);   % liquide
        [propri.rho_G(k),propri.mu_G(k),propri.h_G(k),propri.cp_G(k),propri.Pr_G(k),propri.k_G(k)] = refpropm('DVHC^L','T',results.Tf_numk(k,2),'P',results.P_num(k,2),data.refrigerant); % vapeur
    end
    %                 fprintf('%c%c%c%c%03d%%',8,8,8,8,round(100*(k/FlowSetup.nb_z))); % Progression 0 --> 100 %
end
