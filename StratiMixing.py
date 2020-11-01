%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-- Heat storage tank with thermal stratification - Main routine --%%%%%
%-------------------------------------------------------------------------%
%             Par Pierre-Luc Paradis et Marie-H�l�ne Talbot               %
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
%����������������������������
data.refrigerant = 'CO2';            % [-] R�frig�rant
data.fluidRefGas = 'water';          % [-] Fluide utilis�e dans le r�servoir (Stockage de chaleur)
% data.fluidEvap = 'INCOMP::MPG-40%';  % [-] Fluide utilis�e dans le r�servoir (Stockage de froid)
data.theta = -pi/2;                  % [radians] Inclinaison de la conduite par rapport � la verticale (-pi/2 --> horizontal, 0 --> vertical)
data.g = 9.81;                       % [m/s2] Gravitational acceleration
data.rugosity_ratio = 1e-6;          % [-] Rugosit� relative des conduites (1e-6 --> conduites lisses)
mixing = 1;                          % 1: Active l'Algorithme de mixing / 0: d�sactive l'algorithme de mixing

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

data.k_EC = eau_prop('kf',300);          % [W/(m K)] Thermal conductivity de l'eau � 300[K]
data.rho_EC = eau_prop('vf',300)^-1;     % [kg/m3] Density de l'eau � 300[K]
data.Cp_EC = eau_prop('Cpf',300);        % [J/(kg K)] Heat capacity de l'eau � 300[K]
data.alpha_EC = data.k_EC/(data.rho_EC*data.Cp_EC); % [m^2/s] Thermal diffusivity de l'eau � 300[K]

data.Tairk = 20 + 273.15;            % [K] Temp�rature dans la salle m�canique
data.critere = 1e-6;                 % [-] Crit�re de convergence
% data.Pglycol = 101325;               % [Pa] Pression dans le r�servoir froid

%_______________________________________
% Param�tres de la r�solution temporelle
%���������������������������������������
nb_t = 121;                   % [-] Nombre de pas de temps
dt = 30;                     % [s] Pas de temps ( 300 [s] --> 5 [min], 900 [s] --> 15 [min])
temps = (0:dt:(nb_t-1)*dt)'; % [s] Vecteur des pas de temps
%__________________
% R�servoir #2 - EC
%������������������
% Res2 - Geometry
%
paramRes2.De = 24*0.0254;     % [m] Diam�tre ext�rieure du r�servoir
paramRes2.ye = 44*0.0254;     % [m] Hauteur ext�rieure du r�servoir
paramRes2.t_side = 0.18750*0.0254;                  % [m] �paisseur des parois de c�t�
paramRes2.t_top = 0.125*0.0254;                     % [m] �paisseur de la paroi du dessus
paramRes2.t_bottom = 0.375*0.0254;                  % [m] �paisseur de la paroi du dessous
paramRes2.Di = paramRes2.De - 2*paramRes2.t_side;   % [m] Diam�tre int�rieure du r�servoir
paramRes2.yi = paramRes2.ye - paramRes2.t_top - paramRes2.t_bottom; % [m] Hauteur int�rieure du r�servoir
paramRes2.As = pi*paramRes2.Di^2/4;                                 % [m^2] Aire de la section interne du r�servoir
paramRes2.As_ext = pi*paramRes2.De^2/4;                             % [m^2] Aire de la section externe du r�servoir
paramRes2.Vres = paramRes2.As*paramRes2.yi;                         % [m^3] Volume interne du r�servoir
paramRes2.Vres_ext = paramRes2.As_ext*paramRes2.ye;                 % [m^3] Volume externe du r�servoir
paramRes2.Vss_tank = paramRes2.Vres_ext-paramRes2.Vres;             % [m^3] Volume d'acier inoxydable du r�servoir
paramRes2.t_iso = 2*0.0254;                                         % [m] �paisseur de l'isolant (2 pouces)

% Res2 - Discretization
%
paramRes2.nb_y = 16;                                             % [noeuds] Nombre de noeuds du reservoir de stockage (arbritraire)
paramRes2.dy = paramRes2.yi/paramRes2.nb_y;                      % [m] Distance entre les noeuds (selon l'axe y)
paramRes2.y_num = (paramRes2.dy/2:paramRes2.dy:paramRes2.yi)';   % [m] Vecteur de position de chacun des noeuds

% Res2 - Flow parameters
%
paramRes2.m_in = 0.13;%0.05;                         % [kg/s] D�bit de la pompe CIR2 (2,1 GPM nominal)
paramRes2.intlet_EC = paramRes2.nb_y;          % [-] Num�ro du noeud entr�e (nb_y --> bas du r�servoir)
paramRes2.outlet_EC = 1;                       % [-] Num�ro du noeud sortie (1 --> haut du reservoir)
paramRes2.Ti = 60;                             % [�C] Temp�rature initiale du r�servoir 2
paramRes2.Tik = paramRes2.Ti+273.15;           % [K] Conversion
paramRes2.Tin_EC = 20;                         % [�C] Temp�rature de retour du r�seau de chauffage (Temp�rature d'entr�e EC)
paramRes2.Tink_EC = paramRes2.Tin_EC + 273.15; % [K] Conversion

%__________________
% �changeur EFD
%������������������
% HX EFD - Geometry
%
paramHXefd.Do = 0.75 * 0.0254;                      % [m] Diam�tre ext�rieur du tuyau (3/4 pouces)
paramHXefd.t_pipe = 0.0486 * 0.0254;                % [m] �paisseur de la paroi du tuyau
paramHXefd.Di = paramHXefd.Do-2*paramHXefd.t_pipe;  % [m] diam�tre int�rieur du tuyau
paramHXefd.Ltot = 496*.0254;                        % [m] Longueur du tuyau
paramHXefd.As = pi*paramHXefd.Di^2/4;               % [m^2] Aire de la section interne de la conduite
paramHXefd.Aext = pi*paramHXefd.Do^2/4;             % [m^2] Aire de la section externe de la conduite
paramHXefd.D = 20*0.0254;                           % [m] Diam�tre de la spirale de l'�changeur
paramHXefd.pitch = 2*0.0254;                        % [m] Pitch de la spirale de l'�changeur

% HX EFD - Flow parameters
%
paramHXefd.m_in = 0.5;%0.1;                                % [kg/s] D�bit de consommation d'eau chaude domestique
paramHXefd.intlet = paramRes2.nb_y;                        % [-] Num�ro du noeud entr�e (nb_y --> bas du r�servoir)
paramHXefd.outlet = 1;                                     % [-] Num�ro du noeud sortie (1--> haut du reservoir)
paramHXefd.nb_y = paramHXefd.intlet-(paramHXefd.outlet-1); % [noeuds] Nombre de noeuds o� l'�changeur EFD est pr�sent
paramHXefd.Tin = 4;                                        % [�C] Temp�rature de l'eau potable de la ville (Temp�rature d'entr�e EFD)
paramHXefd.Tink = paramHXefd.Tin + 273.15;                 % [K] Conversion

% HX EFD - Discretization
%
paramHXefd.L = paramHXefd.Ltot/paramHXefd.nb_y;              % [m] Longueur de tuyau de l'�changeur EFD dans un noeud j du r�servoir 2

%__________________
% �changeur CO2
%������������������
% HX CO2 - Geometry
%
paramRefGas.Do = 0.009525;                               % [m] diam�tre ext�rieur du tuyau (3/8 pouces)
paramRefGas.t_pipe = 0.035*0.0254;                       % [m] �paisseur de la paroi du tuyau
paramRefGas.Di = paramRefGas.Do-2*paramRefGas.t_pipe;    % [m] diam�tre int�rieur du tuyau
paramRefGas.Ltot = 6.25*4;                               % [m] Longueur totale de l'�changeur CO2
paramRefGas.As = pi*paramRefGas.Di^2/4;                  % [m^2] Aire de la section inerne de la conduite
paramRefGas.Aext = pi*paramRefGas.Do^2/4;                % [m^2] Aire de la section externe de la conduite

% HX CO2 - Flow parameters
%
paramRefGas.m_in = 0.0146;             % [kg/s] D�bit massique de r�frig�rant
paramRefGas.Pin = 9000;                % [kPa] Pression de sortie du Compresseur
paramRefGas.hf_in = 470*1e3;           % [J/kg] Enthalpie massique du r�frig�rant
paramRefGas.inlet = 1;                 % [-] Num�ro du noeud entree
paramRefGas.outlet = paramRes2.nb_y;   % [-] Num�ro du noeud sortie (1 --> haut du reservoir)
paramRefGas.nb_y = paramRefGas.outlet-paramRefGas.inlet+1; % [-] Nombre de noeud j du r�servoir contenant l'�changeur CO2

% HX CO2 - Discretization
%
paramRefGas.L = paramRefGas.Ltot/paramRes2.nb_y;                         % [m] Longueur de tuyau de l'�changeur CO2 dans un noeud j du r�servoir #2
paramRefGas.dz = 0.0625;                                                 % [m] Distance entre les noeuds (selon l'axe z)
paramRefGas.z_num = (paramRefGas.dz/2:paramRefGas.dz:paramRefGas.L)';    % [m] Vecteur de position de chacun des noeuds
paramRefGas.nb_z = length(paramRefGas.z_num);                            % [noeuds] Nombre de noeuds selon z
paramRefGas.Z_num = (paramRefGas.dz/2:paramRefGas.dz:paramRefGas.Ltot)'; % [m] Vecteur de position global de chacun des noeuds

%__________________
% Res2
%������������������
paramRes2.Veau = paramRes2.Vres -(paramRefGas.Aext*paramRefGas.Ltot+paramHXefd.Aext*paramHXefd.Ltot);     % [m^3] Volume d'eau dans le r�servoir --> Calcul Volume EC = (Volume interne du r�servoir - Volume des �changeurs de chaleurs)
paramRes2.dV = paramRes2.Veau/paramRes2.nb_y;                                                             % [m^3] Volume d'un noeud du r�servoir de stockage
paramRes2.deltak = data.kss*(pi*paramRes2.De^2/4-paramRes2.As)/paramRes2.As;                              % D�-stratification due au transfert de chaleur par conduction axiale dans la paroi du r�servoir

%�������������������������������������������������������������������������%
%--------- Bilan d'�nergie sur volume de contr�le j du r�servoir  --------%
%_________________________________________________________________________%

% Masques
%
M_topbott = zeros(paramRes2.nb_y,1);
M_topbott(1) = 1;        % Identification du noeud du haut (Condition aux fronti�res)
M_topbott(end) = 1;      % Identification du noeud du bas (Condition aux fronti�res)

M_debit_D = zeros(paramRes2.nb_y,1);
M_debit_D(paramRes2.outlet_EC:paramRes2.intlet_EC) = 1; % Identification des noeuds de d�bit EC pour la diagonale D (Matrice des coefficients)

M_debit_A = zeros(paramRes2.nb_y-1,1);
M_debit_A(paramRes2.outlet_EC:paramRes2.intlet_EC-1) = 1; % Identification des noeuds de d�bit EC pour la diagonale A (Matrice des coefficients)

M_debit_C = zeros(paramRes2.nb_y,1);
M_debit_C(paramRes2.intlet_EC) = 1; % Identification des noeuds de d�bit EC pour la diagonale C (Matrice des coefficients)

M_HXefd = zeros(paramRes2.nb_y,1);
M_HXefd(paramHXefd.outlet:paramHXefd.intlet) = 1; % Localisation de l'�changeur EFD

% Initialisation du vecteur de solution du RES2 - Temp�rature EC
%
Res2results.Tk = zeros(paramRes2.nb_y,nb_t) + paramRes2.Tik;   % [K] Initialisation de la matrice de temp�rature EC (ligne = #noeud et colonne = pas de temps)

% Initialisation des vecteurs de solutions de l'�changeur EFD
%
HXefdresults.Tink = zeros(paramRes2.nb_y,nb_t);                % [K] Initialisation de la temp�rature dans l'�changeur EFD � l'entr�e de chaque noeud du r�servoir
HXefdresults.Tink(:,1) = Res2results.Tk(:,1).*M_HXefd;         % [K] Initialisation de la temp�rature dans l'�changeur EFD (Conditions initiales --> Time step #1)
HXefdresults.Tink(paramHXefd.intlet,2:nb_t) = paramHXefd.Tink; % [K] Initialisation de la temp�rature � l'entr�e de l'�changeur d'EFD (Boundary condition)

HXefdresults.Toutk = HXefdresults.Tink;                        % [K] Initialisation de la temp�rature dans l'�changeur EFD � la sortie de chaque noeud du r�servoir

HXefdresults.Ts_outk = zeros(paramRes2.nb_y,nb_t);             % [K] Initialisation de la temp�rature de surface externe de l'�changeur d'EFD
HXefdresults.Ts_outk(:,1) = Res2results.Tk(:,1).*M_HXefd;      % [K] Initialisation de la temp�rature de surface externe de l'�changeur d'EFD (Conditions initiales --> Time step #1)

% Initialisation des Matrices de solutions sur un pas de temps donn� pour l'�changeur CO2
%
hm_num_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y); % [J/kg] Enthalpie du r�frig�rant
P_num_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y);  % [kPa] Pression du r�frig�rant
qp_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y);     % [W/m] Flux de chaleur total
Tf_num_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y); % [K] Champs de temp�rature du r�frig�rant
rho_m_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y);  % [kg/m3] Masse volumique du r�frig�rant
Vz_num_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y); % [m/s] Vitesse du r�frig�rant
cp_TimeStep = zeros(paramRefGas.nb_z,paramRes2.nb_y);     % [J/(kg K)] Chaleur massique du r�frig�rant

% Initialisation des Matrices de solutions globales pour l'�changeur CO2
%
hm_num_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);  % [J/kg] Enthalpie du r�frig�rant
P_num_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);   % [kPa] Pression du r�frig�rant
qp_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);      % [W/m] Flux de chaleur total
Tf_num_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);  % [K] Champs de temp�rature du r�frig�rant
rho_m_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);   % [kg/m3] Masse volumique du r�frig�rant
Vz_num_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);  % [m/s] Vitesse du r�frig�rant
cp_tot = zeros(paramRefGas.nb_z*paramRes2.nb_y,nb_t);      % [J/(kg K)] Chaleur massique du r�frig�rant

% Initialisation des vecteurs pour le calcul des pertes thermiques du r�servoir
%
h_convEC = zeros(paramRes2.nb_y + 2,1);       % [W/(m^2) K] Convection heat transfer coefficient in EC
h_convair = zeros(paramRes2.nb_y + 2,1);      % [W/(m^2) K] Convection heat transfer coefficient in ambiant air
Rpp_pertestopbott = ones(paramRes2.nb_y,1);   % [K (m^2)/W] R�sistance thermique r�sultante pour les extr�mit�es
qpp_pertes = zeros(paramRes2.nb_y + 2,1);     % [W/m^2] Flux de chaleur total des pertes thermiques du r�servoir

% Initialisation des temp�ratures de surfaces des parois du r�servoir
%
Tsk_in = zeros(paramRes2.nb_y + 2,1) + paramRes2.Tik - 4; % [K] Surface interne
Tsk_out = zeros(paramRes2.nb_y + 2,1) + data.Tairk - 4;   % [K] Surface externe

% Initialisation des vecteurs pour le calcul de l'�changeur EFD
%
h_convHXefd_in = zeros(paramRes2.nb_y,1);       % [W/(m^2) K] Convection heat transfer coefficient inside �changeur EFD (convection forc�e)
h_convHXefd_out = zeros(paramRes2.nb_y,1);      % [W/(m^2) K] Convection heat transfer coefficient outside �changeur EFD (convection naturelle)
TmoyEFD = zeros(paramRes2.nb_y,1);              % [K] Temp�rature moyenne entre l'entr�e et la sortie (�changeur EFD)
deltaTi = zeros(paramRes2.nb_y,1);              % [K] Pincement � l'entr�e de l'�changeur EFD
deltaTo = zeros(paramRes2.nb_y,1);              % [K] Pincement � la sortie de l'�changeur EFD

% Initialisation des flux de chaleur des �changeurs
%
q_EFD = zeros(paramRes2.nb_y,1);      % [W] Flux de chaleur total �changeur EFD
q_R744 = zeros(paramRefGas.nb_y,1);   % [W] Flux de chaleur total �changeur CO2

% Calcul des Termes constant dans le temps
%
Rpp_condtanktop = paramRes2.t_top/data.kss;      % [K m^2/W] Calcul de la r�sistance thermique par unit� de surface au travers de la paroi du dessus du r�servoir
Rpp_condtankbott = paramRes2.t_bottom/data.kss;  % [K m^2/W] Calcul de la r�sistance thermique par unit� de surface au travers de la paroi du dessous du r�servoir
Rpp_condtankside = paramRes2.t_side/data.kss;    % [K m^2/W] Calcul de la r�sistance thermique par unit� de surface au travers des parois de c�t� du r�servoir
Rpp_condiso = paramRes2.t_iso/data.kiso;         % [K m^2/W] Calcul de la r�sistance thermique par unit� de surface au travers de la couche d'isolant

R_condHXecd_wall = log(paramHXefd.Do/paramHXefd.Di)/(2*pi*data.kcu*paramHXefd.L);      % [K/W] Calcul de la r�sistance thermique totale au travers de la paroi de l'�changeur EFD

d2 = paramRes2.As*(data.kss + paramRes2.deltak)/paramRes2.dy;                          % Terme constant de la diagonale D (Matrice des coefficients)
B = ones(paramRes2.nb_y-1,1)*(-d2);                                                    % Vecteur de la diagonale inf�rieure B de longueur (Nb_y-1) (Matrice des coefficients)

%____________________
% R�SOLUTION EN TEMPS
%��������������������
for m = 2:nb_t
    iterTres2 = 1;
    limit_iter_Tres2 = 25;
    erreurTres2 = ones(limit_iter_Tres2,1);
    Res2results.Tk(:,m) = Res2results.Tk(:,m-1);  % Initilaisation de la solution du pas de temps courant � partir de la solution du pas de temps pr�c�dent
    
    % Recalcul des j temp�ratures du r�servoir de stockage jusqu'� convergence pour le pas de temps courant, m
    %
    while erreurTres2(iterTres2) > data.critere
        
        Tguess = Res2results.Tk(:,m);  % Enregistrement de la derni�re solution calcul�e (permet de calculer le crit�re de convergence sur la temp�rature EC)
        
        %___________________________________
        % PERTES THERMIQUES DU R�SERVOIR
        %�����������������������������������
        for i = 1:paramRes2.nb_y + 2  % La boucle se fait sur nb_y+2 afin de calculer les pertes thermiques des parois sur les nb_y noeuds du r�servoir ainsi que pour le dessus et le dessous du r�servoir (conditions aux fronti�res)
            
            % Structure "if" pour faire corresponde le bon noeud du r�servoir pour le calcul des propri�t�s pour le dessus et le dessous du r�servoir (Conditions aux fronti�res) � l'aide de l'indice j
            %
            if i == 1
                j = 1;
            elseif i > 1 && i < paramRes2.nb_y + 2
                j = i - 1;
            elseif i == paramRes2.nb_y + 2
                j = i - 2;
            end
            
            % Calcul des propri�t�s de l'eau chaude dans le r�servoir
            %
            propriFluid_EC.rho_inf(i,1) = eau_prop('vf',(Res2results.Tk(j,m)+Tsk_in(i))/2)^-1;     % [kg/m3] Density
            propriFluid_EC.k_inf(i,1) = eau_prop('kf',(Res2results.Tk(j,m)+Tsk_in(i))/2);          % [W/(m K)] Thermal conductivity
            propriFluid_EC.mu_inf(i,1) = eau_prop('muf',(Res2results.Tk(j,m)+Tsk_in(i))/2);        % [Ns/m^2] Dynamic viscosity
            propriFluid_EC.Pr_inf(i,1) = eau_prop('Prf',(Res2results.Tk(j,m)+Tsk_in(i))/2);        % [-] Prandtl Number
            propriFluid_EC.beta_inf(i,1) = eau_prop('betaf',(Res2results.Tk(j,m)+Tsk_in(i))/2);    % [1/K]  Volumetric expansivity (beta)
            propriFluid_EC.Cp_inf(i,1) = eau_prop('Cpf',(Res2results.Tk(j,m)+Tsk_in(i))/2);        % [J/(kg K)] Heat capacity
            propriFluid_EC.alpha_inf(i,1) = propriFluid_EC.k_inf(i)/propriFluid_EC.rho_inf(i)/propriFluid_EC.Cp_inf(i); % [m^2/s] Thermal diffusivity
            
            % Calul des coefficients de convection entre l'eau dans le r�servoir et les parois
            %
            switch i
                case 1 % Pour le haut du r�servoir
                    h_convEC(i) = Churchill_Chu_plate(paramRes2.Di,propriFluid_EC.rho_inf(i),...
                        propriFluid_EC.Pr_inf(i),propriFluid_EC.k_inf(i),propriFluid_EC.beta_inf(i),...
                        propriFluid_EC.alpha_inf(i),propriFluid_EC.mu_inf(i),Tsk_in(i),Res2results.Tk(j,1),3); % [W/(m^2) K] Convection heat transfer coefficient in EC
                case (paramRes2.nb_y+2) % Pour le bas du r�servoir
                    h_convEC(i) = Churchill_Chu_plate(paramRes2.Di,propriFluid_EC.rho_inf(i),...
                        propriFluid_EC.Pr_inf(i),propriFluid_EC.k_inf(i),propriFluid_EC.beta_inf(i),...
                        propriFluid_EC.alpha_inf(i),propriFluid_EC.mu_inf(i),Tsk_in(i),Res2results.Tk(j,1),2); % [W/(m^2) K] Convection heat transfer coefficient in EC
                otherwise % Pour les c�t�s
                    h_convEC(i) = Churchill_Chu_plate(paramRes2.yi,propriFluid_EC.rho_inf(i),...
                        propriFluid_EC.Pr_inf(i),propriFluid_EC.k_inf(i),propriFluid_EC.beta_inf(i),...
                        propriFluid_EC.alpha_inf(i),propriFluid_EC.mu_inf(i),Tsk_in(i),Res2results.Tk(j,1),1); % [W/(m^2) K] Convection heat transfer coefficient in EC
            end
            
            % Calcul des propri�t�s de l'air ambiant qui �change avec le r�servoir
            %
            propriAir.rho(i,1) = air_prop('rho',(data.Tairk+Tsk_out(i))/2);  % [kg/m3] Density
            propriAir.k(i,1) = air_prop('k',(data.Tairk+Tsk_out(i))/2);      % [W/(m K)] Thermal conductivity
            propriAir.mu(i,1) = air_prop('mu',(data.Tairk+Tsk_out(i))/2);    % [Ns/m^2] Dynamic viscosity
            propriAir.Pr(i,1) = air_prop('Pr',(data.Tairk+Tsk_out(i))/2);    % [-] Prandtl number
            propriAir.beta(i,1) = 1/((data.Tairk+Tsk_out(i))/2);             % [1/K]  Volumetric expansivity (beta)
            propriAir.alpha(i,1) = air_prop('al',(data.Tairk+Tsk_out(i))/2); % [m^2/s] Thermal diffusivity
            
            % Calul des coefficients de convection entre les parois du r�servoir et l'air ambiant
            %
            switch i
                case 1 % Pour le dessus du r�servoir
                    h_convair(i) = Churchill_Chu_plate(paramRes2.De,propriAir.rho(i),...
                        propriAir.Pr(i),propriAir.k(i),propriAir.beta(i),...
                        propriAir.alpha(i),propriAir.mu(i),Tsk_out(i),data.Tairk,2); % [W/(m^2) K] Convection heat transfer coefficient
                case (paramRes2.nb_y+2) % Pour le dessous du r�servoir
                    h_convair(i) = Churchill_Chu_plate(paramRes2.De,propriAir.rho(i),...
                        propriAir.Pr(i),propriAir.k(i),propriAir.beta(i),...
                        propriAir.alpha(i),propriAir.mu(i),Tsk_out(i),data.Tairk,3); % [W/(m^2) K] Convection heat transfer coefficient
                otherwise % Pour les c�t�s
                    h_convair(i) = Churchill_Chu_plate(paramRes2.ye,propriAir.rho(i),...
                        propriAir.Pr(i),propriAir.k(i),propriAir.beta(i),...
                        propriAir.alpha(i),propriAir.mu(i),Tsk_out(i),data.Tairk,1); % [W/(m^2) K] Convection heat transfer coefficient
            end
        end
        
        % Calcul des r�sistances des flux de pertes thermiques
        %
        Rpp_convEC = 1./h_convEC;      % [K (m^2)/W] R�sistance thermique par convection naturelle dans l'eau chaude
        Rpp_convair = 1./h_convair;    % [K (m^2)/W] R�sistance thermique par convection naturelle dans l'air ambiant
        Rpp_pertestopbott(1) = Rpp_convEC(1) + Rpp_condtanktop + Rpp_condiso + Rpp_convair(1);        % [K (m^2)/W] R�sistance thermique r�sultante pour le dessus
        Rpp_pertestopbott(end) = Rpp_convEC(end) + Rpp_condtankbott + Rpp_condiso + Rpp_convair(end); % [K (m^2)/W] R�sistance thermique r�sultante pour le dessous
        Rpp_pertesside = Rpp_convEC(2:end-1) + Rpp_condtankside + Rpp_condiso + Rpp_convair(2:end-1); % [K (m^2)/W] R�sistance thermique r�sultante pour les c�t�s
        Rpp_pertes = [Rpp_pertestopbott(1);Rpp_pertesside;Rpp_pertestopbott(end)];                    % [K (m^2)/W] Vecteur des r�sistances thermiques totales pour le dessus, dessous et c�t�s
        
        for i = 1:paramRes2.nb_y + 2
            % Structure if pour faire corresponde le bon noeud du r�servoir pour le calcul de propri�t�s pour le dessus et le dessous du r�servoir (Conditions aux fronti�res) � l'aide de l'indice j
            %
            if i == 1
                j = 1;
            elseif i > 1 && i < paramRes2.nb_y + 2
                j=i-1;
            elseif i == paramRes2.nb_y + 2
                j = i - 2;
            end
            
            % Calcul des pertes thermiques du r�servoir EC
            %
            qpp_pertes(i) = (Res2results.Tk(j,m) - data.Tairk)./Rpp_pertes(i);
            
            % Calcul des temp�ratures de surface � l'int�rieur et � l'ext�rieur des parois du r�servoir EC
            %
            Tsk_in(i) = Res2results.Tk(j,m) - qpp_pertes(i)*Rpp_convEC(i);
            Tsk_out(i) = data.Tairk + qpp_pertes(i)*Rpp_convair(i);
        end
        
        %______________
        % �CHANGEUR EFD
        %��������������
        for i = paramHXefd.intlet:-1:paramHXefd.outlet
            
            if paramHXefd.m_in > 0
                
                if iterTres2 == 1
                    % Hypoth�se sur la temp�rature de sortie du noeud (n�cessaire seulement pour la premi�re it�ration sur le champs de temp�rature EC de chaque pas de temps)
                    %
                    if HXefdresults.Tink(i,m) > Res2results.Tk(i,m)        % Si l'�changeur EFDS chauffe le r�servoir
                        HXefdresults.Toutk(i,m) = HXefdresults.Tink(i,m)-5;
                    else % Si le r�servoir chauffe l'�changeur EFD
                        HXefdresults.Toutk(i,m) = HXefdresults.Tink(i,m)+5;
                    end
                    % Hypoth�se sur la temp�rature de surface externe de l'�changeur EFD
                    %
                    HXefdresults.Ts_outk(i,m) = (((HXefdresults.Tink(i,m) + HXefdresults.Toutk(i,m))/2) + Res2results.Tk(i,m))/2; % [K] On initilaise la temp�rature de surface avec la moyenne de la temp�rature moyenne entre l'entr�e et la sortie et la temp�rature de l'eau dans le r�servoir au noeud i
                end
                iterToutHX = 1;
                erreurTout = ones(50,1);
                % Recalcul des j temp�ratures de sortie de l'�changeur jusqu'� convergence pour le temps donn� m
                %
                while erreurTout(iterToutHX) > data.critere
                    
                    Toutguess = HXefdresults.Toutk(i,m);                               % Enregistrement de la derni�re solution de Temp�rature de sortie calcul�e (permet de valider la convergence sur ToutHX)
                    TmoyEFD(i) = (HXefdresults.Tink(i,m) + HXefdresults.Toutk(i,m))/2; % Calcul de la temp�rature moyenne entre la temp�rature d'entr�e et la temp�rature de sortie de l'�changeur
                    
                    % Calcul des propri�t�s de l'eau dans le r�servoir EC
                    %
                    propriHXefd_out.rho_inf(i,1) = eau_prop('vf',(TmoyEFD(i) + Res2results.Tk(i,m))/2)^-1;                            % [kg/m3] Density
                    propriHXefd_out.k_inf(i,1) = eau_prop('kf',(TmoyEFD(i) + Res2results.Tk(i,m))/2);                                 % [W/(m K)] Thermal conductivity
                    propriHXefd_out.mu_inf(i,1) = eau_prop('muf',(TmoyEFD(i) + Res2results.Tk(i,m))/2);                               % [Ns/m^2] Dynamic viscosity
                    propriHXefd_out.Pr_inf(i,1) = eau_prop('Prf',(TmoyEFD(i) + Res2results.Tk(i,m))/2);                               % [-] Prandtl Number
                    propriHXefd_out.beta_inf(i,1) = eau_prop('betaf',(TmoyEFD(i) + Res2results.Tk(i,m))/2);                           % [1/K]  Volumetric expansivity (beta)
                    propriHXefd_out.Cp_inf(i,1) = eau_prop('Cpf',(TmoyEFD(i) + Res2results.Tk(i,m))/2);                               % [J/(kg K)] Heat capacity
                    propriHXefd_out.alpha_inf(i,1) = propriHXefd_out.k_inf(i)/(propriHXefd_out.rho_inf(i)*propriHXefd_out.Cp_inf(i)); % [m^2/s] Thermal diffusivity
                    
                    % Calul du coefficient de convection naturelle entre l'�changeur EFD et l'eau du r�servoir
                    %
                    h_convHXefd_out(i) = Churchill_Chu(paramHXefd.Do,propriHXefd_out.rho_inf(i),...
                        propriHXefd_out.Pr_inf(i),propriHXefd_out.k_inf(i),propriHXefd_out.beta_inf(i),...
                        propriHXefd_out.alpha_inf(i),propriHXefd_out.mu_inf(i),HXefdresults.Ts_outk(i,m),...
                        Res2results.Tk(i,1)); % [W/(m^2 K)] Convection heat transfer coefficient
                    
                    % Calcul des propri�t�s de l'eau dans l'�changeur EFD
                    %
                    propriHXefd_in.rho_inf(i,1) = eau_prop('vf',TmoyEFD(i))^-1;                                                 % [kg/m3] Density
                    propriHXefd_in.k_inf(i,1) = eau_prop('kf',TmoyEFD(i));                                                      % [W/(m K)] Thermal conductivity
                    propriHXefd_in.mu_inf(i,1) = eau_prop('muf',TmoyEFD(i));                                                    % [Ns/m^2] Dynamic viscosity
                    propriHXefd_in.Pr_inf(i,1) = eau_prop('Prf',TmoyEFD(i));                                                    % [-] Prandtl Number
                    propriHXefd_in.beta_inf(i,1) = eau_prop('betaf',TmoyEFD(i));                                                % [1/K]  Volumetric expansivity (beta)
                    propriHXefd_in.Cp_inf(i,1) = eau_prop('Cpf',TmoyEFD(i));                                                    % [J/(kg K)] Heat capacity
                    propriHXefd_in.alpha_inf(i,1) = propriHXefd_in.k_inf(i)/(propriHXefd_in.rho_inf(i)*propriHXefd_in.Cp_inf(i)); % [m^2/s] Thermal diffusivity
                    
                    % Calul du coefficient de convection forc�e � l'int�rieur de l'�changeur EFD
                    %
                    h_convHXefd_in(i) = Dittus_Boelter_mod(paramHXefd.Di,paramHXefd.D,paramHXefd.pitch,paramHXefd.m_in,propriHXefd_in.Cp_inf(i)...
                        ,propriHXefd_in.k_inf(i),propriHXefd_in.mu_inf(i),propriHXefd_in.Pr_inf(i));      % [W/(m^2 K)] Convection heat transfer coefficient
                    
                    % Calcul des r�sistances des flux de chaleur entre l'�changeur ECD et le r�servoir
                    %
                    R_convHXefd_in(i,1)=(h_convHXefd_in(i)*pi*paramHXefd.Di*paramHXefd.L).^-1;    % [K/W] R�sistance convection interne de l'�changeur EFD
                    R_convHXefd_out(i,1)=(h_convHXefd_out(i)*pi*paramHXefd.Do*paramHXefd.L).^-1;  % [K/W] R�sistance convection externe de l'�changeur EFD
                    R_HXefd(i,1) = R_convHXefd_in(i) + R_condHXecd_wall + R_convHXefd_out(i);     % [K/W] R�sistance totale de l'�changeur EFD
                    
                    % Recalcul de la temp�rature de sortie de l'�changeur pour le noeud du r�servoir i
                    %
                    HXefdresults.Toutk(i,m) = Res2results.Tk(i,m) + (HXefdresults.Tink(i,m) - Res2results.Tk(i,m))*exp((-R_HXefd(i)*paramHXefd.m_in*propriHXefd_in.Cp_inf(i))^-1);
                    
                    % Calcul de DeltaTlm (M�thode LMTD)
                    %
                    deltaTo(i) = Res2results.Tk(i,m)-HXefdresults.Toutk(i,m);
                    deltaTi(i) = Res2results.Tk(i,m)-HXefdresults.Tink(i,m);
                    if deltaTi(i) == deltaTo(i)
                        deltaTlm(i) = 0;
                    else
                        deltaTlm(i) = (deltaTo(i)-deltaTi(i))/log(deltaTo(i)/deltaTi(i));
                    end
                    
                    % Calcul du flux de chaleur de l'�changeur EFD
                    %
                    q_EFD(i) = (R_HXefd(i))^-1 * deltaTlm(i);   % [W]
                    
                    % Recalcul de la temp�rature de surface externe de l'�changeur pour le noeud du r�servoir i
                    %
                    HXefdresults.Ts_outk(i,m) = TmoyEFD(i) + q_EFD(i)*(R_convHXefd_in(i) + R_condHXecd_wall); % [K]
                    
                    % Calcul de l'erreur entre cette it�ration et l'it�ration pr�c�dente
                    %
                    iterToutHX = iterToutHX + 1;
                    erreurTout(iterToutHX,1) = sqrt((HXefdresults.Toutk(i,m) - Toutguess)^2);
                end
                
                % On assigne la temp�rature calcul�e en sortie � la temp�rature d'entr�e du prochain noeud, i-1
                %
                if i == paramHXefd.outlet
                    break
                else
                    HXefdresults.Tink(i-1,m) = HXefdresults.Toutk(i,m);
                end
            end
        end
        
        %______________
        % �CHANGEUR CO2
        %��������������
        if paramRefGas.m_in > 0 
            % Caract�ristiques de l'�coulement de CO2 � l'entr�e de l'�changeur
            %
            RefGasflow.m_in =  paramRefGas.m_in;          % [kg/s] D�bit massique de r�frig�rant
            RefGasflow.Pin = paramRefGas.Pin;             % [kPa] Pression de sortie du Compresseur
            RefGasflow.hf_in = paramRefGas.hf_in;         % [J/kg] Enthalpie massique du r�frig�rant
            
            % Boucle sur chaque noeud du r�servoir de stockage contenant l'�changeur CO2
            %
%             index = 1;
            for i = paramRefGas.inlet:paramRefGas.outlet
                % Assignation de la temp�rature Tinf pour le noeud i utilis�e dans les fonctions propriIniRefGas et modelFlow1D
                % 
                paramRefGas.Tinfk = Res2results.Tk(i,m);
                
                % Initialisation des propri�t�s
                %
                [RefGasflow,propriRefGas,RefGasresults] = propriIni(paramRefGas,RefGasflow,data);
                
                % Calcul des r�sultats pour la longueur de l'�changeur dans
                % 
                [RefGasresults, propriRefGas] = modelFlow1D(data,paramRefGas,RefGasflow,RefGasresults,propriRefGas,1,'Cooling',data.fluidRefGas);

                % Calcul de la chaleur transmise � l'eau du r�servoir par l'�changeur R744
                %
                q_R744(i) = sum(RefGasresults.qp)*paramRefGas.dz;
                
                % On assigne l'enthlapie et la pression calcul�es en sortie comme valeurs � l'entr�e du prochain noeud
                %
                RefGasflow.hf_in = RefGasresults.hm_num(end,2);
                RefGasflow.Pin = RefGasresults.P_num(end,2);
                
                % R�sultats pour le pas de temps donn� pour tous les noeuds du r�servoir
                %
                hm_num_TimeStep(:,i) = RefGasresults.hm_num(:,2);     % [J/kg] Enthalpie du r�frig�rant
                P_num_TimeStep(:,i) = RefGasresults.P_num(:,2);       % [kPa] Pression du r�frig�rant
                qp_TimeStep(:,i) = RefGasresults.qp(:,1);             % [W/m] Flux de chaleur total
                Tf_num_TimeStep(:,i) = RefGasresults.Tf_numk(:,2);    % [K] Champs de temp�rature du r�frig�rant
                rho_m_TimeStep(:,i) = RefGasresults.rho_m(:,1);       % [kg/m3] Masse volumique du r�frig�rant
                Vz_num_TimeStep(:,i) = RefGasresults.Vz_num(:,2);     % [m/s] Vitesse du r�frig�rant
                cp_TimeStep(:,i) = propriRefGas.cp_L(:,1);            % [J/(kg K)] Chaleur massique du r�frig�rant       
            end
            
            % Enregistrement des r�sultats pour l'affichage (colonne = temps & ligne = nb_z*nb_y)
            %
            hm_num_tot(:,m) = hm_num_TimeStep(:);    % [J/kg] Enthalpie du r�frig�rant
            P_num_tot(:,m) = P_num_TimeStep(:);      % [kPa] Pression du r�frig�rant
            qp_tot(:,m) = qp_TimeStep(:);            % [W/m] Flux de chaleur total
            Tf_num_tot(:,m) = Tf_num_TimeStep(:);    % [K] Champs de temp�rature du r�frig�rant
            rho_m_tot(:,m) = rho_m_TimeStep(:);      % [kg/m3] Masse volumique du r�frig�rant
            Vz_num_tot(:,m) = Vz_num_TimeStep(:);    % [m/s] Vitesse du r�frig�rant
            cp_tot(:,m) = cp_TimeStep(:);            % [J/(kg K)] Chaleur massique du r�frig�rant
            
        end
        
        % Calcul des coefficients de la matrice D
        %
        d1 = propriFluid_EC.rho_inf(2:end-1).*propriFluid_EC.Cp_inf(2:end-1)*paramRes2.dV/dt;
        d3 = paramRes2.m_in*propriFluid_EC.Cp_inf(2:end-1);
        d4 = pi*paramRes2.Di*paramRes2.dy./Rpp_pertesside;
        d5 = paramRes2.As./Rpp_pertestopbott ;
        
        %____________________________________________________
        % CALCUL DU CHAMPS DE TEMP�RATURE EC, ALGORITHME TDMA
        %����������������������������������������������������
        A = ones(paramRes2.nb_y-1,1).*(- d2 - d3(1:end-1).*M_debit_A);                                                               % Vecteur de la diagonale sup�rieure de longueur N-1
        D = ones(paramRes2.nb_y,1).*(d1 + d2*(2 - M_topbott) + d3.*M_debit_D + d4 + d5.*M_topbott);                                  % Vecteur de la diagonale centrale de longueur N
        C1 = ones(paramRes2.nb_y,1).*(d3.*M_debit_C*paramRes2.Tink_EC + d4.*data.Tairk + d5.*M_topbott*data.Tairk - q_EFD - q_R744); % Vecteur des constantes C (ind�pendante du temps m)
        C = ones(paramRes2.nb_y,1).*(d1.*Res2results.Tk(:,m-1)) + C1;                                                                % Vecteur des constantes C (d�pendante du temps m)
        
        % R�solution dans l'espace des temp�ratures de l'eau du r�servoir EC
        %
        Res2results.Tk(:,m) = TDMA(A,B,C,D);

        % Calcul du crit�re de convergence
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
        MixingCounter = 1;     % Initialisation du compteur associ� � l'algorithme d'inversion
        Inversion = Res2results.Tk(1:end-1,m) < Res2results.Tk(2:end,m); % identifie les noeuds avec inversion
        LowerNode = find(Inversion, 1, 'last') + 1; % Noeud au bas de l'inversion
        CritereInversion = sum(Inversion);  % Initialisation du crit�re d'inversion
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

% Conversion des r�sultats de K � �C
%
HXefdresults.Tin = HXefdresults.Tink - 273.15;
HXefdresults.Tout = HXefdresults.Toutk - 273.15;
Res2results.T = Res2results.Tk - 273.15;

% Affichage des conditions d'op�ration
%
fprintf('La temp�rature de la pi�ce est de %i�C.\n',data.Tairk-273.15)
fprintf('La temp�rature initial du r�servoir est %i�C.\n',paramRes2.Ti)
fprintf('La temp�rature du retour d''eau chaude est %i�C \n',paramRes2.Tin_EC)
fprintf('avec un d�bit de %.2f kg/s \n',paramRes2.m_in)
fprintf('La temp�rature de l''eau chaude domestique � l''entr�e est %i�C.\n',paramHXefd.Tin)
fprintf('avec un d�bit de %.2f kg/s \n',paramHXefd.m_in)
fprintf('La pression d''entr�e du CO2 est de %ikPa.\n',paramRefGas.Pin)
fprintf('L''enthalpie d''entr�e du CO2 est de %ikJ/kg.\n',paramRefGas.hf_in*1e-3)
fprintf('avec un d�bit de %.4f kg/s \n',paramRefGas.m_in)

fich = ('ResultsWorkspaceV14mixing');
ComputationTime = toc;
save(fich);
%�������������������������������������������������������������������������%
%------------------------ Impression des figures -------------------------%
%_________________________________________________________________________%

PrintFig;

%�������������������������������������������������������������������������%
%------------------- Analyse des �chelles de temps -----------------------%
%_________________________________________________________________________%

disp(['Volume interne du r�servoir = ',num2str(paramRes2.Vres),' [m3]'])
disp(['Volume occup� par l''�changeur EFD = ',num2str(paramHXefd.Aext*paramHXefd.Ltot),' [m3]'])
disp(['Volume occup� par l''�changeur CO2 = ',num2str(paramRefGas.Aext*paramRefGas.Ltot),' [m3]'])
Volume_eau_EC = paramRes2.Vres-paramHXefd.Aext*paramHXefd.Ltot-paramRefGas.Aext*paramRefGas.Ltot;
disp(['Volume d''eau EC = ',num2str(Volume_eau_EC),' [m3]'])
Volume_eau_EFD = paramHXefd.As*paramHXefd.Ltot;
disp(['Volume d''eau dans l''�changeur EFD = ',num2str(Volume_eau_EFD),' [m3]'])
Volume_cuivre = (paramHXefd.Aext-paramHXefd.As)*paramHXefd.Ltot;
disp(['Volume de Cuivre dans l''�changeur EFD = ',num2str(Volume_cuivre),' [m3]'])
disp(['Volume de CO2 dans l''�changeur CO2 = ',num2str(paramRefGas.As*paramRefGas.Ltot),' [m3]'])
Volume_inox = (paramRefGas.Aext-paramRefGas.As)*paramRefGas.Ltot;
disp(['Volume d''Inox dans l''�changeur CO2 = ',num2str(Volume_inox),' [m3]'])

disp(['Masse thermique d''eau EC = ',num2str(data.rho_EC*Volume_eau_EC*data.Cp_EC*1e-3),' [kJ/K]'])
disp(['Masse thermique d''eau EFD � 300[K] = ',num2str(data.rho_EC*Volume_eau_EFD*data.Cp_EC*1e-3),' [kJ/K]'])
disp(['Masse thermique du Cuivre de l''�changeur EFD = ',num2str(data.rho_cu*Volume_cuivre*data.cp_cu*1e-3),' [kJ/K]'])
Masse_thermique_co2_ini = sum(rho_m_tot(:,2).*cp_tot(:,2)*paramRefGas.As*paramRefGas.dz*1e-3); 
Masse_thermique_co2_inter = sum(rho_m_tot(:,round(nb_t/2)).*cp_tot(:,round(nb_t/2))*paramRefGas.As*paramRefGas.dz*1e-3);
Masse_thermique_co2_final = sum(rho_m_tot(:,nb_t).*cp_tot(:,nb_t)*paramRefGas.As*paramRefGas.dz*1e-3);
disp(['Masse thermique du CO2 dans l''�changeur CO2 ini/inter/final = ',num2str(Masse_thermique_co2_ini),' / ',num2str(Masse_thermique_co2_inter),' / ',num2str(Masse_thermique_co2_final), ' [kJ/K]'])
disp(['Masse thermique d''Inox de l''�changeur CO2 = ',num2str(data.rho_ss*Volume_inox*data.cp_ss*1e-3),' [kJ/K]'])
disp(['Masse thermique d''Inox du r�servoir = ',num2str(data.rho_ss*paramRes2.Vss_tank*data.cp_ss*1e-3),' [kJ/K]'])

disp(['Temps de transport EFD = ',num2str(data.rho_EC*Volume_eau_EFD/paramHXefd.m_in),' [s]'])
average_velocity_CO2_ini = mean(Vz_num_tot(:,2));
average_velocity_CO2_inter = mean(Vz_num_tot(:,round(nb_t/2)));
average_velocity_CO2_final = mean(Vz_num_tot(:,nb_t));
disp(['Temps de transport CO2 ini/inter/final = ',num2str(paramRefGas.Ltot/average_velocity_CO2_ini),' / ',num2str(paramRefGas.Ltot/average_velocity_CO2_inter),' / ',num2str(paramRefGas.Ltot/average_velocity_CO2_final),' [s]'])

disp(['Temps caract�ristique du r�servoir = ',num2str(((paramRes2.Di/2)^2/(4*data.alpha_EC))/3600),' [h]'])



toc


