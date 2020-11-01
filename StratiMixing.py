%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-- Heat storage tank with thermal stratification - Main routine --%%%%%
%-------------------------------------------------------------------------%
%             Par Pierre-Luc Paradis et Marie-Hélène Talbot               %
%                            September 2016                               %
%-------------------------------------------------------------------------%
% This program is the main routine that solves the 1-D temperature field
% along the height of the tank (y axis). It also solves the 1-D flow field
% of CO2 (z axis). The heat storage tank also has a heat exchanger who
% preheats domestic hot water.
%
tic
%____________________________
% General data and parameters
%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
data.refrigerant = 'CO2';            % [-] Réfrigérant
data.fluidRefGas = 'water';          % [-] Fluide utilisée dans le réservoir (Stockage de chaleur)
% data.fluidEvap = 'INCOMP::MPG-40%';  % [-] Fluide utilisée dans le réservoir (Stockage de froid)
data.theta = -pi/2;                  % [radians] Inclinaison de la conduite par rapport à la verticale (-pi/2 --> horizontal, 0 --> vertical)
data.g = 9.81;                       % [m/s2] Gravitational acceleration
data.rugosity_ratio = 1e-6;          % [-] Rugosité relative des conduites (1e-6 --> conduites lisses)
mixing = 1;                          % 1: Active l'Algorithme de mixing / 0: désactive l'algorithme de mixing

data.kss = 14.9;                     % [W/(m K)] Thermal conductivity of stainless steel at 300[K]
data.rho_ss = 7900;                  % [kg/m3] Density of stainless steel at 300[K]
data.cp_ss = 477;                    % [J/kgK] Heat capacity of stainless steel at 300[K]
data.alpha_ss = 3.95e-6;             % [m2/s] Thermal diffusivity of stainless steel at 300[K]

data.kcu = 401;                      % [W/(m K)] Thermal conductivity of copper at 300[K]
data.rho_cu = 8933;                  % [kg/m3] Density of copper at 300[K]
data.cp_cu = 385;                    % [J/kgK] Heat capacity of copper at 300[K]
data.alpha_cu = 117e-6;              % [m2/s] Thermal diffusivity of copper at 300[K]

data.kiso = 0.038;                   % [W/(m K)] Thermal conductivity of insulation material (Coated Glass Fiber Blanket; duct liner)at 300[K] 
data.rho_iso = 32;                   % [kg/m3] Density of insulation material (Coated Glass Fiber Blanket; duct liner)at 300[K] 
data.cp_iso = 835;                   % [J/kgK] Heat capacity of insulation material (Coated Glass Fiber Blanket; duct liner)at 300[K] 
data.alpha_iso = data.kiso/(data.rho_iso*data.cp_iso); % [m2/s] Thermal diffusivity of insulation material (Coated Glass Fiber Blanket; duct liner)at 300[K] 

data.k_EC = eau_prop('kf',300);          % [W/(m K)] Thermal conductivity de l'eau à 300[K]
data.rho_EC = eau_prop('vf',300)^-1;     % [kg/m3] Density de l'eau à 300[K]
data.Cp_EC = eau_prop('Cpf',300);        % [J/(kg K)] Heat capacity de l'eau à 300[K]
data.alpha_EC = data.k_EC/(data.rho_EC*data.Cp_EC); % [m^2/s] Thermal diffusivity de l'eau à 300[K]

data.Tairk = 20 + 273.15;            % [K] Température dans la salle mécanique
data.critere = 1e-6;                 % [-] Critère de convergence
% data.Pglycol = 101325;               % [Pa] Pression dans le réservoir froid

%_______________________________________
% Paramètres de la résolution temporelle
%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
nb_t = 121;                   % [-] Nombre de pas de temps
dt = 30;                     % [s] Pas de temps ( 300 [s] --> 5 [min], 900 [s] --> 15 [min])
temps = (0:dt:(nb_t-1)*dt)'; % [s] Vecteur des pas de temps
%__________________
% Réservoir #2 - EC
%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
% Res2 - Geometry
%
paramRes2.De = 24*0.0254;     % [m] Diamètre extérieure du réservoir
paramRes2.ye = 44*0.0254;     % [m] Hauteur extérieure du réservoir
paramRes2.t_side = 0.18750*0.0254;                  % [m] Épaisseur des parois de côté
paramRes2.t_top = 0.125*0.0254;                     % [m] Épaisseur de la paroi du dessus
paramRes2.t_bottom = 0.375*0.0254;                  % [m] Épaisseur de la paroi du dessous
paramRes2.Di = paramRes2.De - 2*paramRes2.t_side;   % [m] Diamètre intérieure du réservoir
paramRes2.yi = paramRes2.ye - paramRes2.t_top - paramRes2.t_bottom; % [m] Hauteur intérieure du réservoir
paramRes2.As = pi*paramRes2.Di^2/4;                                 % [m^2] Aire de la section interne du réservoir
paramRes2.As_ext = pi*paramRes2.De^2/4;                             % [m^2] Aire de la section externe du réservoir
paramRes2.Vres = paramRes2.As*paramRes2.yi;                         % [m^3] Volume interne du réservoir
paramRes2.Vres_ext = paramRes2.As_ext*paramRes2.ye;                 % [m^3] Volume externe du réservoir
paramRes2.Vss_tank = paramRes2.Vres_ext-paramRes2.Vres;             % [m^3] Volume d'acier inoxydable du réservoir
paramRes2.t_iso = 2*0.0254;                                         % [m] Épaisseur de l'isolant (2 pouces)

% Res2 - Discretization
%
paramRes2.nb_y = 16;                                             % [noeuds] Nombre de noeuds du reservoir de stockage (arbritraire)
paramRes2.dy = paramRes2.yi/paramRes2.nb_y;                      % [m] Distance entre les noeuds (selon l'axe y)
paramRes2.y_num = (paramRes2.dy/2:paramRes2.dy:paramRes2.yi)';   % [m] Vecteur de position de chacun des noeuds

% Res2 - Flow parameters
%
paramRes2.m_in = 0.13;%0.05;                         % [kg/s] Débit de la pompe CIR2 (2,1 GPM nominal)
paramRes2.intlet_EC = paramRes2.nb_y;          % [-] Numéro du noeud entrée (nb_y --> bas du réservoir)
paramRes2.outlet_EC = 1;                       % [-] Numéro du noeud sortie (1 --> haut du reservoir)
paramRes2.Ti = 60;                             % [°C] Température initiale du réservoir 2
paramRes2.Tik = paramRes2.Ti+273.15;           % [K] Conversion
paramRes2.Tin_EC = 20;                         % [°C] Température de retour du réseau de chauffage (Température d'entrée EC)
paramRes2.Tink_EC = paramRes2.Tin_EC + 273.15; % [K] Conversion

%__________________
% Échangeur EFD
%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
% HX EFD - Geometry
%
paramHXefd.Do = 0.75 * 0.0254;                      % [m] Diamètre extérieur du tuyau (3/4 pouces)
paramHXefd.t_pipe = 0.0486 * 0.0254;                % [m] épaisseur de la paroi du tuyau
paramHXefd.Di = paramHXefd.Do-2*paramHXefd.t_pipe;  % [m] diamètre intérieur du tuyau
paramHXefd.Ltot = 496*.0254;                        % [m] Longueur du tuyau
paramHXefd.As = pi*paramHXefd.Di^2/4;               % [m^2] Aire de la section interne de la conduite
paramHXefd.Aext = pi*paramHXefd.Do^2/4;             % [m^2] Aire de la section externe de la conduite
paramHXefd.D = 20*0.0254;                           % [m] Diamètre de la spirale de l'échangeur
paramHXefd.pitch = 2*0.0254;                        % [m] Pitch de la spirale de l'échangeur

% HX EFD - Flow parameters
%
paramHXefd.m_in = 0.5;%0.1;                                % [kg/s] Débit de consommation d'eau chaude domestique
paramHXefd.intlet = paramRes2.nb_y;                        % [-] Numéro du noeud entrée (nb_y --> bas du réservoir)
paramHXefd.outlet = 1;                                     % [-] Numéro du noeud sortie (1--> haut du reservoir)
paramHXefd.nb_y = paramHXefd.intlet-(paramHXefd.outlet-1); % [noeuds] Nombre de noeuds où l'échangeur EFD est présent
paramHXefd.Tin = 4;                                        % [°C] Température de l'eau potable de la ville (Température d'entrée EFD)
paramHXefd.Tink = paramHXefd.Tin + 273.15;                 % [K] Conversion

% HX EFD - Discretization
%
paramHXefd.L = paramHXefd.Ltot/paramHXefd.nb_y;              % [m] Longueur de tuyau de l'échangeur EFD dans un noeud j du réservoir 2

%__________________
% Échangeur CO2
%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
% HX CO2 - Geometry
%
paramRefGas.Do = 0.009525;                               % [m] diamètre extérieur du tuyau (3/8 pouces)
paramRefGas.t_pipe = 0.035*0.0254;                       % [m] épaisseur de la paroi du tuyau
paramRefGas.Di = paramRefGas.Do-2*paramRefGas.t_pipe;    % [m] diamètre intérieur du tuyau
paramRefGas.Ltot = 6.25*4;                               % [m] Longueur totale de l'échangeur CO2
paramRefGas.As = pi*paramRefGas.Di^2/4;                  % [m^2] Aire de la section inerne de la conduite
paramRefGas.Aext = pi*paramRefGas.Do^2/4;                % [m^2] Aire de la section externe de la conduite

% HX CO2 - Flow parameters
%
paramRefGas.m_in = 0.0146;             % [kg/s] Débit massique de réfrigérant
paramRefGas.Pin = 9000;                % [kPa] Pression de sortie du Compresseur
paramRefGas.hf_in = 470*1e3;           % [J/kg] Enthalpie massique du réfrigérant
paramRefGas.inlet = 1;                 % [-] Numéro du noeud entree
paramRefGas.outlet = paramRes2.nb_y;   % [-] Numéro du noeud sortie (1 --> haut du reservoir)
paramRefGas.nb_y = paramRefGas.outlet-paramRefGas.inlet+1; % [-] Nombre de noeud j du réservoir contenant l'échangeur CO2

% HX CO2 - Discretization
%
paramRefGas.L = paramRefGas.Ltot/paramRes2.nb_y;                         % [m] Longueur de tuyau de l'échangeur CO2 dans un noeud j du réservoir #2
paramRefGas.dz = 0.0625;                                                 % [m] Distance entre les noeuds (selon l'axe z)
paramRefGas.z_num = (paramRefGas.dz/2:paramRefGas.dz:paramRefGas.L)';    % [m] Vecteur de position de chacun des noeuds
paramRefGas.nb_z = length(paramRefGas.z_num);                            % [noeuds] Nombre de noeuds selon z
paramRefGas.Z_num = (paramRefGas.dz/2:paramRefGas.dz:paramRefGas.Ltot)'; % [m] Vecteur de position global de chacun des noeuds

%__________________
% Res2
%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
paramRes2.Veau = paramRes2.Vres -(paramRefGas.Aext*paramRefGas.Ltot+paramHXefd.Aext*paramHXefd.Ltot);     % [m^3] Volume d'eau dans le réservoir --> Calcul Volume EC = (Volume interne du réservoir - Volume des échangeurs de chaleurs)
paramRes2.dV = paramRes2.Veau/paramRes2.nb_y;                                                             % [m^3] Volume d'un noeud du réservoir de stockage
paramRes2.deltak = data.kss*(pi*paramRes2.De^2/4-paramRes2.As)/paramRes2.As;                              % Dé-stratification due au transfert de chaleur par conduction axiale dans la paroi du réservoir

%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯%
%--------- Bilan d'énergie sur volume de contrôle j du réservoir  --------%
%_________________________________________________________________________%

% Masques
%
M_topbott = zeros(paramRes2.nb_y,1);
M_topbott(1) = 1;        % Identification du noeud du haut (Condition aux frontières)
M_topbott(end) = 1;      % Identification du noeud du bas (Condition aux frontières)

M_debit_D = zeros(paramRes2.nb_y,1);
M_debit_D(paramRes2.outlet_EC:paramRes2.intlet_EC) = 1; % Identification des noeuds de débit EC pour la diagonale D (Matrice des coefficients)

M_debit_A = zeros(paramRes2.nb_y-1,1);
M_debit_A(paramRes2.outlet_EC:paramRes2.intlet_EC-1) = 1; % Identification des noeuds de débit EC pour la diagonale A (Matrice des coefficients)

M_debit_C = zeros(paramRes2.nb_y,1);
M_debit_C(paramRes2.intlet_EC) = 1; % Identification des noeuds de débit EC pour la diagonale C (Matrice des coefficients)

M_HXefd = zeros(paramRes2.nb_y,1);
M_HXefd(paramHXefd.outlet:paramHXefd.intlet) = 1; % Localisation de l'échangeur EFD

% Initialisation du vecteur de solution du RES2 - Température EC
%
Res2results.Tk = zeros(paramRes2.nb_y,nb_t) + paramRes2.Tik;   % [K] Initialisation de la matrice de température EC (ligne = #noeud et colonne = pas de temps)

% Initialisation des vecteurs de solutions de l'échangeur EFD
%
HXefdresults.Tink = zeros(paramRes2.nb_y,nb_t);                % [K] Initialisation de la température dans l'échangeur EFD à l'entrée de chaque noeud du réservoir
HXefdresults.Tink(:,1) = Res2results.Tk(:,1).*M_HXefd;         % [K] Initialisation de la température dans l'échangeur EFD (Conditions initiales --> Time step #1)
HXefdresults.Tink(paramHXefd.intlet,2:nb_t) = paramHXefd.Tink; % [K] Initialisation de la température à l'entrée de l'échangeur d'EFD (Boundary condition)

HXefdresults.Toutk = HXefdresults.Tink;                        % [K] Initialisation de la température dans l'échangeur EFD à la sortie de chaque noeud du réservoir

HXefdresults.Ts_outk = zeros(paramRes2.nb_y,nb_t);             % [K] Initialisation de la température de surface externe de l'échangeur d'EFD
HXefdresults.Ts_outk(:,1) = Res2results.Tk(:,1).*M_HXefd;      % [K] Initialisation de la température de surface externe de l'échangeur d'EFD (Conditions initiales --> Time step #1)

% Initialisation des Matrices de solutions sur un pas de temps donné pour l'échangeur CO2
%
hm_num_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y); % [J/kg] Enthalpie du réfrigérant
P_num_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y);  % [kPa] Pression du réfrigérant
qp_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y);     % [W/m] Flux de chaleur total
Tf_num_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y); % [K] Champs de température du réfrigérant
rho_m_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y);  % [kg/m3] Masse volumique du réfrigérant
Vz_num_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y); % [m/s] Vitesse du réfrigérant
cp_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y);     % [J/(kg K)] Chaleur massique du réfrigérant

% Initialisation des Matrices de solutions globales pour l'échangeur CO2
%
hm_num_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);  % [J/kg] Enthalpie du réfrigérant
P_num_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);   % [kPa] Pression du réfrigérant
qp_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);      % [W/m] Flux de chaleur total
Tf_num_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);  % [K] Champs de température du réfrigérant
rho_m_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);   % [kg/m3] Masse volumique du réfrigérant
Vz_num_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);  % [m/s] Vitesse du réfrigérant
cp_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);      % [J/(kg K)] Chaleur massique du réfrigérant

% Initialisation des vecteurs pour le calcul des pertes thermiques du réservoir
%
h_convEC = zeros(paramRes2.nb_y + 2,1);       % [W/(m^2) K] Convection heat transfer coefficient in EC
h_convair = zeros(paramRes2.nb_y + 2,1);      % [W/(m^2) K] Convection heat transfer coefficient in ambiant air
Rpp_pertestopbott = ones(paramRes2.nb_y,1);   % [K (m^2)/W] Résistance thermique résultante pour les extrémitées
qpp_pertes = zeros(paramRes2.nb_y + 2,1);     % [W/m^2] Flux de chaleur total des pertes thermiques du réservoir

% Initialisation des températures de surfaces des parois du réservoir
%
Tsk_in = zeros(paramRes2.nb_y + 2,1) + paramRes2.Tik - 4; % [K] Surface interne
Tsk_out = zeros(paramRes2.nb_y + 2,1) + data.Tairk - 4;   % [K] Surface externe

% Initialisation des vecteurs pour le calcul de l'échangeur EFD
%
h_convHXefd_in = zeros(paramRes2.nb_y,1);       % [W/(m^2) K] Convection heat transfer coefficient inside échangeur EFD (convection forcée)
h_convHXefd_out = zeros(paramRes2.nb_y,1);      % [W/(m^2) K] Convection heat transfer coefficient outside échangeur EFD (convection naturelle)
TmoyEFD = zeros(paramRes2.nb_y,1);              % [K] Température moyenne entre l'entrée et la sortie (échangeur EFD)
deltaTi = zeros(paramRes2.nb_y,1);              % [K] Pincement à l'entrée de l'échangeur EFD
deltaTo = zeros(paramRes2.nb_y,1);              % [K] Pincement à la sortie de l'échangeur EFD

% Initialisation des flux de chaleur des échangeurs
%
q_EFD = zeros(paramRes2.nb_y,1);      % [W] Flux de chaleur total échangeur EFD
q_R744 = zeros(paramRefGas.nb_y,1);   % [W] Flux de chaleur total échangeur CO2

% Calcul des Termes constant dans le temps
%
Rpp_condtanktop = paramRes2.t_top/data.kss;      % [K m^2/W] Calcul de la résistance thermique par unité de surface au travers de la paroi du dessus du réservoir
Rpp_condtankbott = paramRes2.t_bottom/data.kss;  % [K m^2/W] Calcul de la résistance thermique par unité de surface au travers de la paroi du dessous du réservoir
Rpp_condtankside = paramRes2.t_side/data.kss;    % [K m^2/W] Calcul de la résistance thermique par unité de surface au travers des parois de côté du réservoir
Rpp_condiso = paramRes2.t_iso/data.kiso;         % [K m^2/W] Calcul de la résistance thermique par unité de surface au travers de la couche d'isolant

R_condHXecd_wall = log(paramHXefd.Do/paramHXefd.Di)/(2*pi*data.kcu*paramHXefd.L);      % [K/W] Calcul de la résistance thermique totale au travers de la paroi de l'échangeur EFD

d2 = paramRes2.As*(data.kss + paramRes2.deltak)/paramRes2.dy;                          % Terme constant de la diagonale D (Matrice des coefficients)
B = ones(paramRes2.nb_y-1,1)*(-d2);                                                    % Vecteur de la diagonale inférieure B de longueur (Nb_y-1) (Matrice des coefficients)

%____________________
% RÉSOLUTION EN TEMPS
%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
for m = 2:nb_t
    iterTres2 = 1;
    limit_iter_Tres2 = 25;
    erreurTres2 = ones(limit_iter_Tres2,1);
    Res2results.Tk(:,m) = Res2results.Tk(:,m-1);  % Initilaisation de la solution du pas de temps courant à partir de la solution du pas de temps précédent
    
    % Recalcul des j températures du réservoir de stockage jusqu'à convergence pour le pas de temps courant, m
    %
    while erreurTres2(iterTres2) > data.critere
        
        Tguess = Res2results.Tk(:,m);  % Enregistrement de la dernière solution calculée (permet de calculer le critère de convergence sur la température EC)
        
        %___________________________________
        % PERTES THERMIQUES DU RÉSERVOIR
        %¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
        for i = 1:paramRes2.nb_y + 2  % La boucle se fait sur nb_y+2 afin de calculer les pertes thermiques des parois sur les nb_y noeuds du réservoir ainsi que pour le dessus et le dessous du réservoir (conditions aux frontières)
            
            % Structure "if" pour faire corresponde le bon noeud du réservoir pour le calcul des propriétés pour le dessus et le dessous du réservoir (Conditions aux frontières) à l'aide de l'indice j
            %
            if i == 1
                j = 1;
            elseif i > 1 && i < paramRes2.nb_y + 2
                j = i - 1;
            elseif i == paramRes2.nb_y + 2
                j = i - 2;
            end
            
            % Calcul des propriétés de l'eau chaude dans le réservoir
            %
            propriFluid_EC.rho_inf(i,1) = eau_prop('vf',(Res2results.Tk(j,m)+Tsk_in(i))/2)^-1;     % [kg/m3] Density
            propriFluid_EC.k_inf(i,1) = eau_prop('kf',(Res2results.Tk(j,m)+Tsk_in(i))/2);          % [W/(m K)] Thermal conductivity
            propriFluid_EC.mu_inf(i,1) = eau_prop('muf',(Res2results.Tk(j,m)+Tsk_in(i))/2);        % [Ns/m^2] Dynamic viscosity
            propriFluid_EC.Pr_inf(i,1) = eau_prop('Prf',(Res2results.Tk(j,m)+Tsk_in(i))/2);        % [-] Prandtl Number
            propriFluid_EC.beta_inf(i,1) = eau_prop('betaf',(Res2results.Tk(j,m)+Tsk_in(i))/2);    % [1/K]  Volumetric expansivity (beta)
            propriFluid_EC.Cp_inf(i,1) = eau_prop('Cpf',(Res2results.Tk(j,m)+Tsk_in(i))/2);        % [J/(kg K)] Heat capacity
            propriFluid_EC.alpha_inf(i,1) = propriFluid_EC.k_inf(i)/propriFluid_EC.rho_inf(i)/propriFluid_EC.Cp_inf(i); % [m^2/s] Thermal diffusivity
            
            % Calul des coefficients de convection entre l'eau dans le réservoir et les parois
            %
            switch i
                case 1 % Pour le haut du réservoir
                    h_convEC(i) = Churchill_Chu_plate(paramRes2.Di,propriFluid_EC.rho_inf(i),...
                        propriFluid_EC.Pr_inf(i),propriFluid_EC.k_inf(i),propriFluid_EC.beta_inf(i),...
                        propriFluid_EC.alpha_inf(i),propriFluid_EC.mu_inf(i),Tsk_in(i),Res2results.Tk(j,1),3); % [W/(m^2) K] Convection heat transfer coefficient in EC
                case (paramRes2.nb_y+2) % Pour le bas du réservoir
                    h_convEC(i) = Churchill_Chu_plate(paramRes2.Di,propriFluid_EC.rho_inf(i),...
                        propriFluid_EC.Pr_inf(i),propriFluid_EC.k_inf(i),propriFluid_EC.beta_inf(i),...
                        propriFluid_EC.alpha_inf(i),propriFluid_EC.mu_inf(i),Tsk_in(i),Res2results.Tk(j,1),2); % [W/(m^2) K] Convection heat transfer coefficient in EC
                otherwise % Pour les côtés
                    h_convEC(i) = Churchill_Chu_plate(paramRes2.yi,propriFluid_EC.rho_inf(i),...
                        propriFluid_EC.Pr_inf(i),propriFluid_EC.k_inf(i),propriFluid_EC.beta_inf(i),...
                        propriFluid_EC.alpha_inf(i),propriFluid_EC.mu_inf(i),Tsk_in(i),Res2results.Tk(j,1),1); % [W/(m^2) K] Convection heat transfer coefficient in EC
            end
            
            % Calcul des propriétés de l'air ambiant qui échange avec le réservoir
            %
            propriAir.rho(i,1) = air_prop('rho',(data.Tairk+Tsk_out(i))/2);  % [kg/m3] Density
            propriAir.k(i,1) = air_prop('k',(data.Tairk+Tsk_out(i))/2);      % [W/(m K)] Thermal conductivity
            propriAir.mu(i,1) = air_prop('mu',(data.Tairk+Tsk_out(i))/2);    % [Ns/m^2] Dynamic viscosity
            propriAir.Pr(i,1) = air_prop('Pr',(data.Tairk+Tsk_out(i))/2);    % [-] Prandtl number
            propriAir.beta(i,1) = 1/((data.Tairk+Tsk_out(i))/2);             % [1/K]  Volumetric expansivity (beta)
            propriAir.alpha(i,1) = air_prop('al',(data.Tairk+Tsk_out(i))/2); % [m^2/s] Thermal diffusivity
            
            % Calul des coefficients de convection entre les parois du réservoir et l'air ambiant
            %
            switch i
                case 1 % Pour le dessus du réservoir
                    h_convair(i) = Churchill_Chu_plate(paramRes2.De,propriAir.rho(i),...
                        propriAir.Pr(i),propriAir.k(i),propriAir.beta(i),...
                        propriAir.alpha(i),propriAir.mu(i),Tsk_out(i),data.Tairk,2); % [W/(m^2) K] Convection heat transfer coefficient
                case (paramRes2.nb_y+2) % Pour le dessous du réservoir
                    h_convair(i) = Churchill_Chu_plate(paramRes2.De,propriAir.rho(i),...
                        propriAir.Pr(i),propriAir.k(i),propriAir.beta(i),...
                        propriAir.alpha(i),propriAir.mu(i),Tsk_out(i),data.Tairk,3); % [W/(m^2) K] Convection heat transfer coefficient
                otherwise % Pour les côtés
                    h_convair(i) = Churchill_Chu_plate(paramRes2.ye,propriAir.rho(i),...
                        propriAir.Pr(i),propriAir.k(i),propriAir.beta(i),...
                        propriAir.alpha(i),propriAir.mu(i),Tsk_out(i),data.Tairk,1); % [W/(m^2) K] Convection heat transfer coefficient
            end
        end
        
        % Calcul des résistances des flux de pertes thermiques
        %
        Rpp_convEC = 1./h_convEC;      % [K (m^2)/W] Résistance thermique par convection naturelle dans l'eau chaude
        Rpp_convair = 1./h_convair;    % [K (m^2)/W] Résistance thermique par convection naturelle dans l'air ambiant
        Rpp_pertestopbott(1) = Rpp_convEC(1) + Rpp_condtanktop + Rpp_condiso + Rpp_convair(1);        % [K (m^2)/W] Résistance thermique résultante pour le dessus
        Rpp_pertestopbott(end) = Rpp_convEC(end) + Rpp_condtankbott + Rpp_condiso + Rpp_convair(end); % [K (m^2)/W] Résistance thermique résultante pour le dessous
        Rpp_pertesside = Rpp_convEC(2:end-1) + Rpp_condtankside + Rpp_condiso + Rpp_convair(2:end-1); % [K (m^2)/W] Résistance thermique résultante pour les côtés
        Rpp_pertes = [Rpp_pertestopbott(1);Rpp_pertesside;Rpp_pertestopbott(end)];                    % [K (m^2)/W] Vecteur des résistances thermiques totales pour le dessus, dessous et côtés
        
        for i = 1:paramRes2.nb_y + 2
            % Structure if pour faire corresponde le bon noeud du réservoir pour le calcul de propriétés pour le dessus et le dessous du réservoir (Conditions aux frontières) à l'aide de l'indice j
            %
            if i == 1
                j = 1;
            elseif i > 1 && i < paramRes2.nb_y + 2
                j=i-1;
            elseif i == paramRes2.nb_y + 2
                j = i - 2;
            end
            
            % Calcul des pertes thermiques du réservoir EC
            %
            qpp_pertes(i) = (Res2results.Tk(j,m) - data.Tairk)./Rpp_pertes(i);
            
            % Calcul des températures de surface à l'intérieur et à l'extérieur des parois du réservoir EC
            %
            Tsk_in(i) = Res2results.Tk(j,m) - qpp_pertes(i)*Rpp_convEC(i);
            Tsk_out(i) = data.Tairk + qpp_pertes(i)*Rpp_convair(i);
        end
        
        %______________
        % ÉCHANGEUR EFD
        %¯¯¯¯¯¯¯¯¯¯¯¯¯¯
        for i = paramHXefd.intlet:-1:paramHXefd.outlet
            
            if paramHXefd.m_in > 0
                
                if iterTres2 == 1
                    % Hypothèse sur la température de sortie du noeud (nécessaire seulement pour la première itération sur le champs de température EC de chaque pas de temps)
                    %
                    if HXefdresults.Tink(i,m) > Res2results.Tk(i,m)        % Si l'échangeur EFDS chauffe le réservoir
                        HXefdresults.Toutk(i,m) = HXefdresults.Tink(i,m)-5;
                    else % Si le réservoir chauffe l'échangeur EFD
                        HXefdresults.Toutk(i,m) = HXefdresults.Tink(i,m)+5;
                    end
                    % Hypothèse sur la température de surface externe de l'échangeur EFD
                    %
                    HXefdresults.Ts_outk(i,m) = (((HXefdresults.Tink(i,m) + HXefdresults.Toutk(i,m))/2) + Res2results.Tk(i,m))/2; % [K] On initilaise la température de surface avec la moyenne de la température moyenne entre l'entrée et la sortie et la température de l'eau dans le réservoir au noeud i
                end
                iterToutHX = 1;
                erreurTout = ones(50,1);
                % Recalcul des j températures de sortie de l'échangeur jusqu'à convergence pour le temps donné m
                %
                while erreurTout(iterToutHX) > data.critere
                    
                    Toutguess = HXefdresults.Toutk(i,m);                               % Enregistrement de la dernière solution de Température de sortie calculée (permet de valider la convergence sur ToutHX)
                    TmoyEFD(i) = (HXefdresults.Tink(i,m) + HXefdresults.Toutk(i,m))/2; % Calcul de la température moyenne entre la température d'entrée et la température de sortie de l'échangeur
                    
                    % Calcul des propriétés de l'eau dans le réservoir EC
                    %
                    propriHXefd_out.rho_inf(i,1) = eau_prop('vf',(TmoyEFD(i) + Res2results.Tk(i,m))/2)^-1;                            % [kg/m3] Density
                    propriHXefd_out.k_inf(i,1) = eau_prop('kf',(TmoyEFD(i) + Res2results.Tk(i,m))/2);                                 % [W/(m K)] Thermal conductivity
                    propriHXefd_out.mu_inf(i,1) = eau_prop('muf',(TmoyEFD(i) + Res2results.Tk(i,m))/2);                               % [Ns/m^2] Dynamic viscosity
                    propriHXefd_out.Pr_inf(i,1) = eau_prop('Prf',(TmoyEFD(i) + Res2results.Tk(i,m))/2);                               % [-] Prandtl Number
                    propriHXefd_out.beta_inf(i,1) = eau_prop('betaf',(TmoyEFD(i) + Res2results.Tk(i,m))/2);                           % [1/K]  Volumetric expansivity (beta)
                    propriHXefd_out.Cp_inf(i,1) = eau_prop('Cpf',(TmoyEFD(i) + Res2results.Tk(i,m))/2);                               % [J/(kg K)] Heat capacity
                    propriHXefd_out.alpha_inf(i,1) = propriHXefd_out.k_inf(i)/(propriHXefd_out.rho_inf(i)*propriHXefd_out.Cp_inf(i)); % [m^2/s] Thermal diffusivity
                    
                    % Calul du coefficient de convection naturelle entre l'échangeur EFD et l'eau du réservoir
                    %
                    h_convHXefd_out(i) = Churchill_Chu(paramHXefd.Do,propriHXefd_out.rho_inf(i),...
                        propriHXefd_out.Pr_inf(i),propriHXefd_out.k_inf(i),propriHXefd_out.beta_inf(i),...
                        propriHXefd_out.alpha_inf(i),propriHXefd_out.mu_inf(i),HXefdresults.Ts_outk(i,m),...
                        Res2results.Tk(i,1)); % [W/(m^2 K)] Convection heat transfer coefficient
                    
                    % Calcul des propriétés de l'eau dans l'échangeur EFD
                    %
                    propriHXefd_in.rho_inf(i,1) = eau_prop('vf',TmoyEFD(i))^-1;                                                 % [kg/m3] Density
                    propriHXefd_in.k_inf(i,1) = eau_prop('kf',TmoyEFD(i));                                                      % [W/(m K)] Thermal conductivity
                    propriHXefd_in.mu_inf(i,1) = eau_prop('muf',TmoyEFD(i));                                                    % [Ns/m^2] Dynamic viscosity
                    propriHXefd_in.Pr_inf(i,1) = eau_prop('Prf',TmoyEFD(i));                                                    % [-] Prandtl Number
                    propriHXefd_in.beta_inf(i,1) = eau_prop('betaf',TmoyEFD(i));                                                % [1/K]  Volumetric expansivity (beta)
                    propriHXefd_in.Cp_inf(i,1) = eau_prop('Cpf',TmoyEFD(i));                                                    % [J/(kg K)] Heat capacity
                    propriHXefd_in.alpha_inf(i,1) = propriHXefd_in.k_inf(i)/(propriHXefd_in.rho_inf(i)*propriHXefd_in.Cp_inf(i)); % [m^2/s] Thermal diffusivity
                    
                    % Calul du coefficient de convection forcée à l'intérieur de l'échangeur EFD
                    %
                    h_convHXefd_in(i) = Dittus_Boelter_mod(paramHXefd.Di,paramHXefd.D,paramHXefd.pitch,paramHXefd.m_in,propriHXefd_in.Cp_inf(i)...
                        ,propriHXefd_in.k_inf(i),propriHXefd_in.mu_inf(i),propriHXefd_in.Pr_inf(i));      % [W/(m^2 K)] Convection heat transfer coefficient
                    
                    % Calcul des résistances des flux de chaleur entre l'échangeur ECD et le réservoir
                    %
                    R_convHXefd_in(i,1)=(h_convHXefd_in(i)*pi*paramHXefd.Di*paramHXefd.L).^-1;    % [K/W] Résistance convection interne de l'échangeur EFD
                    R_convHXefd_out(i,1)=(h_convHXefd_out(i)*pi*paramHXefd.Do*paramHXefd.L).^-1;  % [K/W] Résistance convection externe de l'échangeur EFD
                    R_HXefd(i,1) = R_convHXefd_in(i) + R_condHXecd_wall + R_convHXefd_out(i);     % [K/W] Résistance totale de l'échangeur EFD
                    
                    % Recalcul de la température de sortie de l'échangeur pour le noeud du réservoir i
                    %
                    HXefdresults.Toutk(i,m) = Res2results.Tk(i,m) + (HXefdresults.Tink(i,m) - Res2results.Tk(i,m))*exp((-R_HXefd(i)*paramHXefd.m_in*propriHXefd_in.Cp_inf(i))^-1);
                    
                    % Calcul de DeltaTlm (Méthode LMTD)
                    %
                    deltaTo(i) = Res2results.Tk(i,m)-HXefdresults.Toutk(i,m);
                    deltaTi(i) = Res2results.Tk(i,m)-HXefdresults.Tink(i,m);
                    if deltaTi(i) == deltaTo(i)
                        deltaTlm(i) = 0;
                    else
                        deltaTlm(i) = (deltaTo(i)-deltaTi(i))/log(deltaTo(i)/deltaTi(i));
                    end
                    
                    % Calcul du flux de chaleur de l'échangeur EFD
                    %
                    q_EFD(i) = (R_HXefd(i))^-1 * deltaTlm(i);   % [W]
                    
                    % Recalcul de la température de surface externe de l'échangeur pour le noeud du réservoir i
                    %
                    HXefdresults.Ts_outk(i,m) = TmoyEFD(i) + q_EFD(i)*(R_convHXefd_in(i) + R_condHXecd_wall); % [K]
                    
                    % Calcul de l'erreur entre cette itération et l'itération précédente
                    %
                    iterToutHX = iterToutHX + 1;
                    erreurTout(iterToutHX,1) = sqrt((HXefdresults.Toutk(i,m) - Toutguess)^2);
                end
                
                % On assigne la température calculée en sortie à la température d'entrée du prochain noeud, i-1
                %
                if i == paramHXefd.outlet
                    break
                else
                    HXefdresults.Tink(i-1,m) = HXefdresults.Toutk(i,m);
                end
            end
        end
        
        %______________
        % ÉCHANGEUR CO2
        %¯¯¯¯¯¯¯¯¯¯¯¯¯¯
        if paramRefGas.m_in > 0 
            % Caractéristiques de l'écoulement de CO2 à l'entrée de l'échangeur
            %
            RefGasflow.m_in =  paramRefGas.m_in;          % [kg/s] Débit massique de réfrigérant
            RefGasflow.Pin = paramRefGas.Pin;             % [kPa] Pression de sortie du Compresseur
            RefGasflow.hf_in = paramRefGas.hf_in;         % [J/kg] Enthalpie massique du réfrigérant
            
            % Boucle sur chaque noeud du réservoir de stockage contenant l'échangeur CO2
            %
%             index = 1;
            for i = paramRefGas.inlet:paramRefGas.outlet
                % Assignation de la température Tinf pour le noeud i utilisée dans les fonctions propriIniRefGas et modelFlow1D
                % 
                paramRefGas.Tinfk = Res2results.Tk(i,m);
                
                % Initialisation des propriétés
                %
                [RefGasflow,propriRefGas,RefGasresults] = propriIni(paramRefGas,RefGasflow,data);
                
                % Calcul des résultats pour la longueur de l'échangeur dans
                % 
                [RefGasresults, propriRefGas] = modelFlow1D(data,paramRefGas,RefGasflow,RefGasresults,propriRefGas,1,'Cooling',data.fluidRefGas);

                % Calcul de la chaleur transmise à l'eau du réservoir par l'échangeur R744
                %
                q_R744(i) = sum(RefGasresults.qp)*paramRefGas.dz;
                
                % On assigne l'enthlapie et la pression calculées en sortie comme valeurs à l'entrée du prochain noeud
                %
                RefGasflow.hf_in = RefGasresults.hm_num(end,2);
                RefGasflow.Pin = RefGasresults.P_num(end,2);
                
                % Résultats pour le pas de temps donné pour tous les noeuds du réservoir
                %
                hm_num_TimeStep(:,i) = RefGasresults.hm_num(:,2);     % [J/kg] Enthalpie du réfrigérant
                P_num_TimeStep(:,i) = RefGasresults.P_num(:,2);       % [kPa] Pression du réfrigérant
                qp_TimeStep(:,i) = RefGasresults.qp(:,1);             % [W/m] Flux de chaleur total
                Tf_num_TimeStep(:,i) = RefGasresults.Tf_numk(:,2);    % [K] Champs de température du réfrigérant
                rho_m_TimeStep(:,i) = RefGasresults.rho_m(:,1);       % [kg/m3] Masse volumique du réfrigérant
                Vz_num_TimeStep(:,i) = RefGasresults.Vz_num(:,2);     % [m/s] Vitesse du réfrigérant
                cp_TimeStep(:,i) = propriRefGas.cp_L(:,1);            % [J/(kg K)] Chaleur massique du réfrigérant       
            end
            
            % Enregistrement des résultats pour l'affichage (colonne = temps & ligne = nb_z*nb_y)
            %
            hm_num_tot(:,m) = hm_num_TimeStep(:);    % [J/kg] Enthalpie du réfrigérant
            P_num_tot(:,m) = P_num_TimeStep(:);      % [kPa] Pression du réfrigérant
            qp_tot(:,m) = qp_TimeStep(:);            % [W/m] Flux de chaleur total
            Tf_num_tot(:,m) = Tf_num_TimeStep(:);    % [K] Champs de température du réfrigérant
            rho_m_tot(:,m) = rho_m_TimeStep(:);      % [kg/m3] Masse volumique du réfrigérant
            Vz_num_tot(:,m) = Vz_num_TimeStep(:);    % [m/s] Vitesse du réfrigérant
            cp_tot(:,m) = cp_TimeStep(:);            % [J/(kg K)] Chaleur massique du réfrigérant
            
        end
        
        % Calcul des coefficients de la matrice D
        %
        d1 = propriFluid_EC.rho_inf(2:end-1).*propriFluid_EC.Cp_inf(2:end-1)*paramRes2.dV/dt;
        d3 = paramRes2.m_in*propriFluid_EC.Cp_inf(2:end-1);
        d4 = pi*paramRes2.Di*paramRes2.dy./Rpp_pertesside;
        d5 = paramRes2.As./Rpp_pertestopbott ;
        
        %____________________________________________________
        % CALCUL DU CHAMPS DE TEMPÉRATURE EC, ALGORITHME TDMA
        %¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
        A = ones(paramRes2.nb_y-1,1).*(- d2 - d3(1:end-1).*M_debit_A);                                                               % Vecteur de la diagonale supérieure de longueur N-1
        D = ones(paramRes2.nb_y,1).*(d1 + d2*(2 - M_topbott) + d3.*M_debit_D + d4 + d5.*M_topbott);                                  % Vecteur de la diagonale centrale de longueur N
        C1 = ones(paramRes2.nb_y,1).*(d3.*M_debit_C*paramRes2.Tink_EC + d4.*data.Tairk + d5.*M_topbott*data.Tairk - q_EFD - q_R744); % Vecteur des constantes C (indépendante du temps m)
        C = ones(paramRes2.nb_y,1).*(d1.*Res2results.Tk(:,m-1)) + C1;                                                                % Vecteur des constantes C (dépendante du temps m)
        
        % Résolution dans l'espace des températures de l'eau du réservoir EC
        %
        Res2results.Tk(:,m) = TDMA(A,B,C,D);

        % Calcul du critère de convergence
        % 
        iterTres2 = iterTres2 + 1;
        erreurTres2(iterTres2,1) = sqrt(sum((Res2results.Tk(:,m) - Tguess).^2))/paramRes2.nb_y;
        fprintf([' \n CALCULATION COMPLETED - PAS DE TEMPS #',num2str(m),'/',num2str(nb_t),'- ITERATION #',num2str(iterTres2),'\n'])

        if iterTres2 > limit_iter_Tres2
            disp(['Temperature field not converged - ResiduTres2 = ',erreurTres2(iterTres2,1)]);
            break
        end
    end
    
    if mixing
        MixingCounter = 1;     % Initialisation du compteur associé à l'algorithme d'inversion
        Inversion = Res2results.Tk(1:end-1,m) < Res2results.Tk(2:end,m); % identifie les noeuds avec inversion
        LowerNode = find(Inversion, 1, 'last') + 1; % Noeud au bas de l'inversion
        CritereInversion = sum(Inversion);  % Initialisation du critère d'inversion
        if CritereInversion > 0
            while CritereInversion > 0
                Res2results.Tk(LowerNode-MixingCounter:LowerNode,m) = mean(Res2results.Tk(LowerNode-MixingCounter:LowerNode,m));
                Inversion = Res2results.Tk(1:end-1,m) < Res2results.Tk(2:end,m);
                CritereInversion = sum(Inversion);
                MixingCounter = MixingCounter + 1;
            end
        end        
    end
end

% Conversion des résultats de K à °C
%
HXefdresults.Tin = HXefdresults.Tink - 273.15;
HXefdresults.Tout = HXefdresults.Toutk - 273.15;
Res2results.T = Res2results.Tk - 273.15;

% Affichage des conditions d'opération
%
fprintf('La température de la pièce est de %i°C.\n',data.Tairk-273.15)
fprintf('La température initial du réservoir est %i°C.\n',paramRes2.Ti)
fprintf('La température du retour d''eau chaude est %i°C \n',paramRes2.Tin_EC)
fprintf('avec un débit de %.2f kg/s \n',paramRes2.m_in)
fprintf('La température de l''eau chaude domestique à l''entrée est %i°C.\n',paramHXefd.Tin)
fprintf('avec un débit de %.2f kg/s \n',paramHXefd.m_in)
fprintf('La pression d''entrée du CO2 est de %ikPa.\n',paramRefGas.Pin)
fprintf('L''enthalpie d''entrée du CO2 est de %ikJ/kg.\n',paramRefGas.hf_in*1e-3)
fprintf('avec un débit de %.4f kg/s \n',paramRefGas.m_in)

fich = ('ResultsWorkspaceV14mixing');
ComputationTime = toc;
save(fich);
%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯%
%------------------------ Impression des figures -------------------------%
%_________________________________________________________________________%

PrintFig;

%¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯%
%------------------- Analyse des échelles de temps -----------------------%
%_________________________________________________________________________%

disp(['Volume interne du réservoir = ',num2str(paramRes2.Vres),' [m3]'])
disp(['Volume occupé par l''échangeur EFD = ',num2str(paramHXefd.Aext*paramHXefd.Ltot),' [m3]'])
disp(['Volume occupé par l''échangeur CO2 = ',num2str(paramRefGas.Aext*paramRefGas.Ltot),' [m3]'])
Volume_eau_EC = paramRes2.Vres-paramHXefd.Aext*paramHXefd.Ltot-paramRefGas.Aext*paramRefGas.Ltot;
disp(['Volume d''eau EC = ',num2str(Volume_eau_EC),' [m3]'])
Volume_eau_EFD = paramHXefd.As*paramHXefd.Ltot;
disp(['Volume d''eau dans l''échangeur EFD = ',num2str(Volume_eau_EFD),' [m3]'])
Volume_cuivre = (paramHXefd.Aext-paramHXefd.As)*paramHXefd.Ltot;
disp(['Volume de Cuivre dans l''échangeur EFD = ',num2str(Volume_cuivre),' [m3]'])
disp(['Volume de CO2 dans l''échangeur CO2 = ',num2str(paramRefGas.As*paramRefGas.Ltot),' [m3]'])
Volume_inox = (paramRefGas.Aext-paramRefGas.As)*paramRefGas.Ltot;
disp(['Volume d''Inox dans l''échangeur CO2 = ',num2str(Volume_inox),' [m3]'])

disp(['Masse thermique d''eau EC = ',num2str(data.rho_EC*Volume_eau_EC*data.Cp_EC*1e-3),' [kJ/K]'])
disp(['Masse thermique d''eau EFD à 300[K] = ',num2str(data.rho_EC*Volume_eau_EFD*data.Cp_EC*1e-3),' [kJ/K]'])
disp(['Masse thermique du Cuivre de l''échangeur EFD = ',num2str(data.rho_cu*Volume_cuivre*data.cp_cu*1e-3),' [kJ/K]'])
Masse_thermique_co2_ini = sum(rho_m_tot(:,2).*cp_tot(:,2)*paramRefGas.As*paramRefGas.dz*1e-3); 
Masse_thermique_co2_inter = sum(rho_m_tot(:,round(nb_t/2)).*cp_tot(:,round(nb_t/2))*paramRefGas.As*paramRefGas.dz*1e-3);
Masse_thermique_co2_final = sum(rho_m_tot(:,nb_t).*cp_tot(:,nb_t)*paramRefGas.As*paramRefGas.dz*1e-3);
disp(['Masse thermique du CO2 dans l''échangeur CO2 ini/inter/final = ',num2str(Masse_thermique_co2_ini),' / ',num2str(Masse_thermique_co2_inter),' / ',num2str(Masse_thermique_co2_final), ' [kJ/K]'])
disp(['Masse thermique d''Inox de l''échangeur CO2 = ',num2str(data.rho_ss*Volume_inox*data.cp_ss*1e-3),' [kJ/K]'])
disp(['Masse thermique d''Inox du réservoir = ',num2str(data.rho_ss*paramRes2.Vss_tank*data.cp_ss*1e-3),' [kJ/K]'])

disp(['Temps de transport EFD = ',num2str(data.rho_EC*Volume_eau_EFD/paramHXefd.m_in),' [s]'])
average_velocity_CO2_ini = mean(Vz_num_tot(:,2));
average_velocity_CO2_inter = mean(Vz_num_tot(:,round(nb_t/2)));
average_velocity_CO2_final = mean(Vz_num_tot(:,nb_t));
disp(['Temps de transport CO2 ini/inter/final = ',num2str(paramRefGas.Ltot/average_velocity_CO2_ini),' / ',num2str(paramRefGas.Ltot/average_velocity_CO2_inter),' / ',num2str(paramRefGas.Ltot/average_velocity_CO2_final),' [s]'])

disp(['Temps caractéristique du réservoir = ',num2str(((paramRes2.Di/2)^2/(4*data.alpha_EC))/3600),' [h]'])



toc


