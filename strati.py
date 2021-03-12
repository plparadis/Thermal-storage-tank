"""
###########################################################################
####-- Heat storage tank with thermal stratification - Main routine --#####
#-------------------------------------------------------------------------#
#            Par Pierre-Luc Paradis et Marie-H�l�ne Talbot               #
#                            September 2016                               #
#-------------------------------------------------------------------------#
# This program is the main routine that solves the 1-D temperature field
# along the height of the tank (y axis).
# The heat storage tank also has a heat exchanger who preheats domestic hot water.
#
"""
import numpy as np
from MatlabStructModule import MatlabStruct
import properties_mod as properties
from TDMA_func import TDMA
from Churchill_Chu_plate_func import Churchill_Chu_plate
from Churchill_Chu_func import Churchill_Chu
from Colebrook_func import Colebrook
from Dittus_Boelter_func import Dittus_Boelter
from Dittus_Boelter_mod_func import Dittus_Boelter_mod

data = MatlabStruct()
paramRes2 = MatlabStruct()
paramHXefd = MatlabStruct()
Res2results = MatlabStruct()
HXefdresults = MatlabStruct()
propriFluid_EC = MatlabStruct()
propriAir = MatlabStruct()
propriHXefd_out = MatlabStruct()
propriHXefd_in = MatlabStruct()

# General data and parameters
data.theta = -np.pi / 2  # [radians] Inclinaison de la conduite par rapport à la verticale (-pi/2 --> horizontal, 0 --> vertical)
data.g = 9.81  # [m/s2] Gravitational acceleration
data.rugosity_ratio = 1e-6  # [-] Rugosite relative des conduites (1e-6 --> conduites lisses)

data.kss = 14.9  # [W/(m K)] Thermal conductivity of stainless steel at 300[K]
data.rho_ss = 7900  # [kg/m3] Density of stainless steel at 300[K]
data.cp_ss = 477  # [J/kgK] Heat capacity of stainless steel at 300[K]
data.alpha_ss = 3.95e-6  # [m2/s] Thermal diffusivity of stainless steel at 300[K]

data.kcu = 401  # [W/(m K)] Thermal conductivity of copper at 300[K]
data.rho_cu = 8933  # [kg/m3] Density of copper at 300[K]
data.cp_cu = 385  # [J/kgK] Heat capacity of copper at 300[K]
data.alpha_cu = 117e-6  # [m2/s] Thermal diffusivity of copper at 300[K]

data.kiso = 0.038  # [W/(m K)] Thermal conductivity of insulation material (Coated Glass Fiber Blanket; duct liner)at 300[K]
data.rho_iso = 32  # [kg/m3] Density of insulation material (Coated Glass Fiber Blanket; duct liner)at 300[K]
data.cp_iso = 835  # [J/kgK] Heat capacity of insulation material (Coated Glass Fiber Blanket; duct liner)at 300[K]
data.alpha_iso = data.kiso / (
            data.rho_iso * data.cp_iso)  # [m2/s] Thermal diffusivity of insulation material (Coated Glass Fiber Blanket; duct liner)at 300[K]

data.k_EC = properties.eau_prop('kf', 300)  # [W/(m K)] Thermal conductivity de l'eau � 300[K]
data.rho_EC = properties.eau_prop('vf', 300) ** -1  # [kg/m3] Density de l'eau � 300[K]
data.Cp_EC = properties.eau_prop('Cpf', 300)  # [J/(kg K)] Heat capacity de l'eau � 300[K]
data.alpha_EC = data.k_EC / (data.rho_EC * data.Cp_EC)  # [m^2/s] Thermal diffusivity de l'eau � 300[K]

data.Tairk = 20 + 273.15  # [K] Temp�rature dans la salle m�canique
data.critere = 1e-6  # [-] Crit�re de convergence

# Parametres de la resolution temporelle
nb_t = 10  # [-] Nombre de pas de temps
dt = 900  # [s] Pas de temps ( 300 [s] --> 5 [min], 900 [s] --> 15 [min])
temps = np.arange(start=0, stop=(nb_t) * dt, step=dt).transpose()  # [s] Vecteur des pas de temps

# Reservoir #2 - EC
# Res2 - Geometry
paramRes2.De = 24 * 0.0254  # [m] Diam�tre ext�rieure du r�servoir
paramRes2.ye = 44 * 0.0254  # [m] Hauteur ext�rieure du r�servoir
paramRes2.t_side = 0.18750 * 0.0254  # [m] �paisseur des parois de c�t�
paramRes2.t_top = 0.125 * 0.0254  # [m] �paisseur de la paroi du dessus
paramRes2.t_bottom = 0.375 * 0.0254  # [m] �paisseur de la paroi du dessous
paramRes2.Di = paramRes2.De - 2 * paramRes2.t_side  # [m] Diam�tre int�rieure du r�servoir
paramRes2.yi = paramRes2.ye - paramRes2.t_top - paramRes2.t_bottom  # [m] Hauteur int�rieure du r�servoir
paramRes2.As = np.pi * paramRes2.Di ** 2 / 4  # [m^2] Aire de la section interne du r�servoir
paramRes2.As_ext = np.pi * paramRes2.De ** 2 / 4  # [m^2] Aire de la section externe du r�servoir
paramRes2.Vres = paramRes2.As * paramRes2.yi  # [m^3] Volume interne du r�servoir
paramRes2.Vres_ext = paramRes2.As_ext * paramRes2.ye  # [m^3] Volume externe du r�servoir
paramRes2.Vss_tank = paramRes2.Vres_ext - paramRes2.Vres  # [m^3] Volume d'acier inoxydable du r�servoir
paramRes2.t_iso = 2 * 0.0254  # [m] �paisseur de l'isolant (2 pouces)

# Res2 - Discretization
paramRes2.nb_y = 5  # [noeuds] Nombre de noeuds du reservoir de stockage (arbritraire)
paramRes2.dy = paramRes2.yi / paramRes2.nb_y  # [m] Distance entre les noeuds (selon l'axe y)
paramRes2.y_num = np.arange(paramRes2.dy / 2, paramRes2.yi,
                            paramRes2.dy)  # [m] Vecteur de position de chacun des noeuds

# Res2 - Flow parameters
paramRes2.m_in = 0.13  # [kg/s] Débit de la pompe
paramRes2.intlet_EC = paramRes2.nb_y  # [-] Numéro du noeud entrée (nb_y --> bas du réservoir)
paramRes2.outlet_EC = 0  # [-] Num�ro du noeud sortie (0 --> haut du reservoir)
paramRes2.Ti = 60  # [°C] Température initiale du réservoir 2
paramRes2.Tik = paramRes2.Ti + 273.15  # [K] Conversion
paramRes2.Tin_EC = 20  # [°C] Température de retour du réseau de chauffage (Température d'entrée EC)
paramRes2.Tink_EC = paramRes2.Tin_EC + 273.15  # [K] Conversion

# Échangeur EFD
# HX EFD - Geometry
paramHXefd.Do = 0.75 * 0.0254  # [m] Diamètre extérieur du tuyau (3/4 pouces)
paramHXefd.t_pipe = 0.0486 * 0.0254  # [m] épaisseur de la paroi du tuyau
paramHXefd.Di = paramHXefd.Do - 2 * paramHXefd.t_pipe  # [m] diamètre intérieur du tuyau
paramHXefd.Ltot = 496 * .0254  # [m] Longueur du tuyau
paramHXefd.As = np.pi * paramHXefd.Di ** 2 / 4  # [m^2] Aire de la section interne de la conduite
paramHXefd.Aext = np.pi * paramHXefd.Do ** 2 / 4  # [m^2] Aire de la section externe de la conduite
paramHXefd.D = 20 * 0.0254  # [m] Diam�tre de la spirale de l'échangeur
paramHXefd.pitch = 2 * 0.0254; # [m] Pitch de la spirale de l'échangeur

# HX EFD - Flow parameters
paramHXefd.m_in = 0.5  # [kg/s] Débit de consommation d'eau chaude domestique
paramHXefd.intlet = paramRes2.nb_y  # [-] Numéro du noeud entrée (nb_y --> bas du réservoir)
paramHXefd.outlet = 0  # [-] Numéro du noeud sortie (0--> haut du reservoir)
paramHXefd.nb_y = paramHXefd.intlet - (
            paramHXefd.outlet - 1)  # [noeuds] Nombre de noeuds où l'échangeur EFD est présent
paramHXefd.Tin = 4  # [°C] Température de l'eau potable de la ville (Température d'entrée EFD)
paramHXefd.Tink = paramHXefd.Tin + 273.15  # [K] Conversion

# HX EFD - Discretization
paramHXefd.L = paramHXefd.Ltot / paramHXefd.nb_y  # [m] Longueur de tuyau de l'échangeur EFD dans un noeud j du réservoir 2

# Res2
paramRes2.Veau = paramRes2.Vres - (
            paramHXefd.Aext * paramHXefd.Ltot)  # [m^3] Volume d'eau dans le réservoir --> Calcul Volume EC = (Volume interne du r�servoir - Volume des �changeurs de chaleurs)
paramRes2.dV = paramRes2.Veau / paramRes2.nb_y  # [m^3] Volume d'un noeud du réservoir de stockage
paramRes2.deltak = data.kss * (
            np.pi * paramRes2.De ** 2 / 4 - paramRes2.As) / paramRes2.As;  # De-stratification due au transfert de chaleur par conduction axiale dans la paroi du r�servoir


# Bilan d'énergie sur volume de contrôle j du réservoir

# Masques
M_topbott = np.zeros([paramRes2.nb_y, 1])
M_topbott[0] = 1  # Identification du noeud du haut (Condition aux frontières)
M_topbott[-1] = 1  # Identification du noeud du bas (Condition aux frontières)

M_debit_D = np.zeros([paramRes2.nb_y, 1])
M_debit_D[
paramRes2.outlet_EC:paramRes2.intlet_EC] = 1  # Identification des noeuds de débit EC pour la diagonale D (Matrice des coefficients)

M_debit_A = np.zeros([paramRes2.nb_y - 1, 1])
M_debit_A[
paramRes2.outlet_EC:paramRes2.intlet_EC - 1] = 1  # Identification des noeuds de débit EC pour la diagonale A (Matrice des coefficients)

M_debit_C = np.zeros([paramRes2.nb_y, 1])
M_debit_C[
    paramRes2.intlet_EC - 1] = 1  # Identification des noeuds de débit EC pour la diagonale C (Matrice des coefficients)

M_HXefd = np.zeros([paramRes2.nb_y, 1])
M_HXefd[paramHXefd.outlet:paramHXefd.intlet] = 1  # Localisation de l'échangeur EFD

# Initialisation du vecteur de solution du RES2 - Température EC
Res2results.Tk = np.zeros([paramRes2.nb_y,
                           nb_t]) + paramRes2.Tik  # [K] Initialisation de la matrice de température EC (ligne = #noeud et colonne = pas de temps)

# Initialisation des vecteurs de solutions de l'échangeur EFD
HXefdresults.Tink = np.zeros([paramRes2.nb_y,
                              nb_t])  # [K] Initialisation de la température dans l'échangeur EFD à l'entrée de chaque noeud du réservoir
HXefdresults.Tink[:, [0]] = Res2results.Tk[:, [
                                                  0]] * M_HXefd  # [K] Initialisation de la température dans l'échangeur EFD (Conditions initiales --> Time step #1)
HXefdresults.Tink[paramHXefd.intlet - 1, np.arange(1, nb_t,
                                                   1)] = paramHXefd.Tink  # [K] Initialisation de la température à l'entrée de l'échangeur d'EFD (Boundary condition)

HXefdresults.Toutk = np.copy(
    HXefdresults.Tink)  # [K] Initialisation de la température dans l'échangeur EFD à la sortie de chaque noeud du réservoir

HXefdresults.Ts_outk = np.zeros(
    [paramRes2.nb_y, nb_t])  # [K] Initialisation de la temp�rature de surface externe de l'échangeur d'EFD
HXefdresults.Ts_outk[:, [1]] = Res2results.Tk[:, [
                                                     0]] * M_HXefd  # [K] Initialisation de la température de surface externe de l'échangeur d'EFD (Conditions initiales --> Time step #1)

# Initialisation des vecteurs pour le calcul des pertes thermiques du réservoir
h_convEC = np.zeros([paramRes2.nb_y + 2, 1]);  # [W/(m^2) K] Convection heat transfer coefficient in EC
h_convair = np.zeros([paramRes2.nb_y + 2, 1]);  # [W/(m^2) K] Convection heat transfer coefficient in ambiant air
Rpp_pertestopbott = np.ones([paramRes2.nb_y, 1]);  # [K (m^2)/W] R�sistance thermique r�sultante pour les extr�mit�es
qpp_pertes = np.zeros([paramRes2.nb_y + 2, 1]);  # [W/m^2] Flux de chaleur total des pertes thermiques du r�servoir

# Initialisation des temp�ratures de surfaces des parois du r�servoir
#
Tsk_in = np.zeros([paramRes2.nb_y + 2, 1]) + paramRes2.Tik - 4;  # [K] Surface interne
Tsk_out = np.zeros([paramRes2.nb_y + 2, 1]) + data.Tairk - 4;  # [K] Surface externe

# Initialisation des vecteurs pour le calcul de l'�changeur EFD
#
h_convHXefd_in = np.zeros(
    [paramRes2.nb_y, 1]);  # [W/(m^2) K] Convection heat transfer coefficient inside �changeur EFD (convection forc�e)
h_convHXefd_out = np.zeros([paramRes2.nb_y,
                            1]);  # [W/(m^2) K] Convection heat transfer coefficient outside �changeur EFD (convection naturelle)
TmoyEFD = np.zeros([paramRes2.nb_y, 1]);  # [K] Temp�rature moyenne entre l'entr�e et la sortie (�changeur EFD)
deltaTi = np.zeros([paramRes2.nb_y, 1]);  # [K] Pincement � l'entr�e de l'�changeur EFD
deltaTo = np.zeros([paramRes2.nb_y, 1]);  # [K] Pincement � la sortie de l'�changeur EFD
R_convHXefd_in = np.zeros([paramRes2.nb_y, 1]);  # [K/W] R�sistance convection interne de l'�changeur EFD
R_convHXefd_out = np.zeros([paramRes2.nb_y, 1]);  # [K/W] R�sistance convection externe de l'�changeur EFD
R_HXefd = np.zeros([paramRes2.nb_y, 1]);  # [K/W] R�sistance totale de l'�changeur EFD
deltaTlm = np.zeros([paramRes2.nb_y, 1]);  # Log Mean Temperature Difference

# Initialisation des flux de chaleur des �changeurs
#
q_EFD = np.zeros([paramRes2.nb_y, 1]);  # [W] Flux de chaleur total �changeur EFD

# Calcul des Termes constant dans le temps
#
Rpp_condtanktop = paramRes2.t_top / data.kss;  # [K m^2/W] Calcul de la r�sistance thermique par unit� de surface au travers de la paroi du dessus du r�servoir
Rpp_condtankbott = paramRes2.t_bottom / data.kss;  # [K m^2/W] Calcul de la r�sistance thermique par unit� de surface au travers de la paroi du dessous du r�servoir
Rpp_condtankside = paramRes2.t_side / data.kss;  # [K m^2/W] Calcul de la r�sistance thermique par unit� de surface au travers des parois de c�t� du r�servoir
Rpp_condiso = paramRes2.t_iso / data.kiso;  # [K m^2/W] Calcul de la r�sistance thermique par unit� de surface au travers de la couche d'isolant

R_condHXecd_wall = np.log(paramHXefd.Do / paramHXefd.Di) / (
            2 * np.pi * data.kcu * paramHXefd.L);  # [K/W] Calcul de la r�sistance thermique totale au travers de la paroi de l'�changeur EFD

d2 = paramRes2.As * (
            data.kss + paramRes2.deltak) / paramRes2.dy;  # Terme constant de la diagonale D (Matrice des coefficients)
B = np.ones([paramRes2.nb_y - 1, 1]) * (
    -d2);  # Vecteur de la diagonale inf�rieure B de longueur (Nb_y-1) (Matrice des coefficients)

# %%
# _________________________
# Initalisation des propri
# ��������������������
propriFluid_EC.rho_inf = np.zeros([paramRes2.nb_y + 2, 1])
propriFluid_EC.k_inf = np.zeros([paramRes2.nb_y + 2, 1])
propriFluid_EC.mu_inf = np.zeros([paramRes2.nb_y + 2, 1])
propriFluid_EC.Pr_inf = np.zeros([paramRes2.nb_y + 2, 1])
propriFluid_EC.beta_inf = np.zeros([paramRes2.nb_y + 2, 1])
propriFluid_EC.Cp_inf = np.zeros([paramRes2.nb_y + 2, 1])
propriFluid_EC.alpha_inf = np.zeros([paramRes2.nb_y + 2, 1])

propriAir.rho = np.zeros([paramRes2.nb_y + 2, 1])
propriAir.k = np.zeros([paramRes2.nb_y + 2, 1])
propriAir.mu = np.zeros([paramRes2.nb_y + 2, 1])
propriAir.Pr = np.zeros([paramRes2.nb_y + 2, 1])
propriAir.beta = np.zeros([paramRes2.nb_y + 2, 1])
propriAir.alpha = np.zeros([paramRes2.nb_y + 2, 1])

propriHXefd_out.rho_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_out.k_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_out.mu_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_out.Pr_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_out.beta_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_out.Cp_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_out.alpha_inf = np.zeros([paramHXefd.nb_y, 1])

propriHXefd_in.rho_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_in.k_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_in.mu_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_in.Pr_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_in.beta_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_in.Cp_inf = np.zeros([paramHXefd.nb_y, 1])
propriHXefd_in.alpha_inf = np.zeros([paramHXefd.nb_y, 1])

# %%
# ____________________
# R�SOLUTION EN TEMPS
# ��������������������

for m in range(1, nb_t, 1):
    iterTres2 = 0;
    limit_iter_Tres2 = 25;
    erreurTres2 = np.ones([limit_iter_Tres2, 1]);
    Res2results.Tk[:, [m]] = Res2results.Tk[:, [
                                                   m - 1]];  # Initilaisation de la solution du pas de temps courant � partir de la solution du pas de temps pr�c�dent

    # Recalcul des j temp�ratures du r�servoir de stockage jusqu'� convergence pour le pas de temps courant, m
    #
    while erreurTres2[iterTres2] > data.critere:

        Tguess = Res2results.Tk[:, [
                                       m]];  # Enregistrement de la derni�re solution calcul�e (permet de calculer le crit�re de convergence sur la temp�rature EC)

        # ___________________________________
        # PERTES THERMIQUES DU R�SERVOIR
        # �����������������������������������
        for i in range(
                paramRes2.nb_y + 2):  # La boucle se fait sur nb_y+2 afin de calculer les pertes thermiques des parois sur les nb_y noeuds du r�servoir ainsi que pour le dessus et le dessous du r�servoir (conditions aux fronti�res)
            # Structure "if" pour faire corresponde le bon noeud du r�servoir pour le calcul des propri�t�s pour le dessus et le dessous du r�servoir (Conditions aux fronti�res) � l'aide de l'indice j
            #
            if i == 0:
                j = 0;
            elif i > 1 and i < paramRes2.nb_y + 1:
                j = i - 1;
            elif i == paramRes2.nb_y + 1:
                j = i - 2;

            # Calcul des propri�t�s de l'eau chaude dans le r�servoir
            #
            propriFluid_EC.rho_inf[i, 0] = properties.eau_prop('vf', (
                        Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2) ** -1;  # [kg/m3] Density
            propriFluid_EC.k_inf[i, 0] = properties.eau_prop('kf', (
                        Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2);  # [W/(m K)] Thermal conductivity
            propriFluid_EC.mu_inf[i, 0] = properties.eau_prop('muf', (
                        Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2);  # [Ns/m^2] Dynamic viscosity
            propriFluid_EC.Pr_inf[i, 0] = properties.eau_prop('Prf', (
                        Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2);  # [-] Prandtl Number
            propriFluid_EC.beta_inf[i, 0] = properties.eau_prop('betaf', (
                        Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2);  # [1/K]  Volumetric expansivity (beta)
            propriFluid_EC.Cp_inf[i, 0] = properties.eau_prop('Cpf', (
                        Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2);  # [J/(kg K)] Heat capacity
            propriFluid_EC.alpha_inf[i, 0] = propriFluid_EC.k_inf[i, 0] / propriFluid_EC.rho_inf[i, 0] / \
                                             propriFluid_EC.Cp_inf[i, 0];  # [m^2/s] Thermal diffusivity

            # Calul des coefficients de convection entre l'eau dans le r�servoir et les parois
            #
            if i == 0:  # Pour le haut du r�servoir
                h_convEC[i, 0] = Churchill_Chu_plate(paramRes2.Di, propriFluid_EC.rho_inf[i, 0],
                                                     propriFluid_EC.Pr_inf[i, 0], propriFluid_EC.k_inf[i, 0],
                                                     propriFluid_EC.beta_inf[i, 0],
                                                     propriFluid_EC.alpha_inf[i, 0], propriFluid_EC.mu_inf[i, 0],
                                                     Tsk_in[i, 0], Res2results.Tk[j, 1],
                                                     3);  # [W/(m^2) K] Convection heat transfer coefficient in EC
            elif i == (paramRes2.nb_y + 1):  # Pour le bas du r�servoir
                h_convEC[i, 0] = Churchill_Chu_plate(paramRes2.Di, propriFluid_EC.rho_inf[i, 0],
                                                     propriFluid_EC.Pr_inf[i, 0], propriFluid_EC.k_inf[i, 0],
                                                     propriFluid_EC.beta_inf[i, 0],
                                                     propriFluid_EC.alpha_inf[i, 0], propriFluid_EC.mu_inf[i, 0],
                                                     Tsk_in[i, 0], Res2results.Tk[j, 1],
                                                     2);  # [W/(m^2) K] Convection heat transfer coefficient in EC
            else:  # Pour les c�t�s
                h_convEC[i, 0] = Churchill_Chu_plate(paramRes2.yi, propriFluid_EC.rho_inf[i, 0],
                                                     propriFluid_EC.Pr_inf[i, 0], propriFluid_EC.k_inf[i, 0],
                                                     propriFluid_EC.beta_inf[i, 0],
                                                     propriFluid_EC.alpha_inf[i, 0], propriFluid_EC.mu_inf[i, 0],
                                                     Tsk_in[i, 0], Res2results.Tk[j, 1],
                                                     1);  # [W/(m^2) K] Convection heat transfer coefficient in EC

            # Calcul des propri�t�s de l'air ambiant qui �change avec le r�servoir
            #
            propriAir.rho[i, 0] = properties.air_prop('rho', (data.Tairk + Tsk_out[i, 0]) / 2);  # [kg/m3] Density
            propriAir.k[i, 0] = properties.air_prop('k',
                                                    (data.Tairk + Tsk_out[i, 0]) / 2);  # [W/(m K)] Thermal conductivity
            propriAir.mu[i, 0] = properties.air_prop('mu',
                                                     (data.Tairk + Tsk_out[i, 0]) / 2);  # [Ns/m^2] Dynamic viscosity
            propriAir.Pr[i, 0] = properties.air_prop('Pr', (data.Tairk + Tsk_out[i, 0]) / 2);  # [-] Prandtl number
            propriAir.beta[i, 0] = 1 / ((data.Tairk + Tsk_out[i, 0]) / 2);  # [1/K]  Volumetric expansivity (beta)
            propriAir.alpha[i, 0] = properties.air_prop('al', (
                        data.Tairk + Tsk_out[i, 0]) / 2);  # [m^2/s] Thermal diffusivity

            # Calul des coefficients de convection entre les parois du r�servoir et l'air ambiant
            #
            if i == 0:  # Pour le dessus du r�servoir
                h_convair[i, 0] = Churchill_Chu_plate(paramRes2.De, propriAir.rho[i, 0],
                                                      propriAir.Pr[i, 0], propriAir.k[i, 0], propriAir.beta[i, 0],
                                                      propriAir.alpha[i, 0], propriAir.mu[i, 0], Tsk_out[i, 0],
                                                      data.Tairk,
                                                      2);  # [W/(m^2) K] Convection heat transfer coefficient
            elif i == (paramRes2.nb_y + 1):  # Pour le dessous du r�servoir
                h_convair[i, 0] = Churchill_Chu_plate(paramRes2.De, propriAir.rho[i, 0],
                                                      propriAir.Pr[i, 0], propriAir.k[i, 0], propriAir.beta[i, 0],
                                                      propriAir.alpha[i, 0], propriAir.mu[i, 0], Tsk_out[i, 0],
                                                      data.Tairk,
                                                      3);  # [W/(m^2) K] Convection heat transfer coefficient
            else:  # Pour les c�t�s
                h_convair[i, 0] = Churchill_Chu_plate(paramRes2.ye, propriAir.rho[i, 0],
                                                      propriAir.Pr[i, 0], propriAir.k[i, 0], propriAir.beta[i, 0],
                                                      propriAir.alpha[i, 0], propriAir.mu[i, 0], Tsk_out[i, 0],
                                                      data.Tairk,
                                                      1);  # [W/(m^2) K] Convection heat transfer coefficient

        # Calcul des r�sistances des flux de pertes thermiques
        #
        Rpp_convEC = 1 / h_convEC;  # [K (m^2)/W] R�sistance thermique par convection naturelle dans l'eau chaude
        Rpp_convair = 1 / h_convair;  # [K (m^2)/W] R�sistance thermique par convection naturelle dans l'air ambiant
        Rpp_pertestopbott[0] = Rpp_convEC[0] + Rpp_condtanktop + Rpp_condiso + Rpp_convair[
            0];  # [K (m^2)/W] R�sistance thermique r�sultante pour le dessus
        Rpp_pertestopbott[-1] = Rpp_convEC[-1] + Rpp_condtankbott + Rpp_condiso + Rpp_convair[
            -1];  # [K (m^2)/W] R�sistance thermique r�sultante pour le dessous
        Rpp_pertesside = Rpp_convEC[1:-1] + Rpp_condtankside + Rpp_condiso + Rpp_convair[
                                                                             1:-1];  # [K (m^2)/W] R�sistance thermique r�sultante pour les c�t�s
        Rpp_pertes = np.vstack((Rpp_pertestopbott[0], Rpp_pertesside, Rpp_pertestopbott[
            -1]));  # [K (m^2)/W] Vecteur des r�sistances thermiques totales pour le dessus, dessous et c�t�s

        for i in range(paramRes2.nb_y + 2):
            # Structure if pour faire corresponde le bon noeud du r�servoir pour le calcul de propri�t�s pour le dessus et le dessous du r�servoir (Conditions aux fronti�res) � l'aide de l'indice j
            #
            if i == 0:
                j = 0;
            elif i > 1 and i < paramRes2.nb_y + 1:
                j = i - 1;
            elif i == paramRes2.nb_y + 1:
                j = i - 2;

            # Calcul des pertes thermiques du r�servoir EC
            #
            qpp_pertes[i, 0] = (Res2results.Tk[j, m] - data.Tairk) / Rpp_pertes[i, 0];

            # Calcul des temp�ratures de surface � l'int�rieur et � l'ext�rieur des parois du r�servoir EC
            #
            Tsk_in[i, 0] = Res2results.Tk[j, m] - qpp_pertes[i, 0] * Rpp_convEC[i, 0];
            Tsk_out[i, 0] = data.Tairk + qpp_pertes[i, 0] * Rpp_convair[i, 0];

        # ______________
        # �CHANGEUR EFD
        # ��������������
        for ii in range(paramHXefd.intlet, paramHXefd.outlet, -1):
            i = ii - 1
            if paramHXefd.m_in > 0:

                if iterTres2 == 0:
                    # Hypoth�se sur la temp�rature de sortie du noeud (n�cessaire seulement pour la premi�re it�ration sur le champs de temp�rature EC de chaque pas de temps)
                    #
                    if HXefdresults.Tink[i, m] > Res2results.Tk[i, m]:  # Si l'�changeur EFDS chauffe le r�servoir
                        HXefdresults.Toutk[i, m] = HXefdresults.Tink[i, m] - 5;
                    else:  # Si le r�servoir chauffe l'�changeur EFD
                        HXefdresults.Toutk[i, m] = HXefdresults.Tink[i, m] + 5;
                    # Hypoth�se sur la temp�rature de surface externe de l'�changeur EFD
                    #
                    HXefdresults.Ts_outk[i, m] = (((HXefdresults.Tink[i, m] + HXefdresults.Toutk[i, m]) / 2) +
                                                  Res2results.Tk[
                                                      i, m]) / 2;  # [K] On initilaise la temp�rature de surface avec la moyenne de la temp�rature moyenne entre l'entr�e et la sortie et la temp�rature de l'eau dans le r�servoir au noeud i
                iterToutHX = 0;
                erreurTout = np.ones([50, 1]);
                # Recalcul des j temp�ratures de sortie de l'�changeur jusqu'� convergence pour le temps donn� m
                #
                while erreurTout[iterToutHX] > data.critere:

                    Toutguess = HXefdresults.Toutk[
                        i, m];  # Enregistrement de la derni�re solution de Temp�rature de sortie calcul�e (permet de valider la convergence sur ToutHX)
                    TmoyEFD[i, 0] = (HXefdresults.Tink[i, m] + HXefdresults.Toutk[
                        i, m]) / 2;  # Calcul de la temp�rature moyenne entre la temp�rature d'entr�e et la temp�rature de sortie de l'�changeur

                    # Calcul des propri�t�s de l'eau dans le r�servoir EC
                    #
                    propriHXefd_out.rho_inf[i, 0] = properties.eau_prop('vf', (
                                TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2) ** -1;  # [kg/m3] Density
                    propriHXefd_out.k_inf[i, 0] = properties.eau_prop('kf', (
                                TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2);  # [W/(m K)] Thermal conductivity
                    propriHXefd_out.mu_inf[i, 0] = properties.eau_prop('muf', (
                                TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2);  # [Ns/m^2] Dynamic viscosity
                    propriHXefd_out.Pr_inf[i, 0] = properties.eau_prop('Prf', (
                                TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2);  # [-] Prandtl Number
                    propriHXefd_out.beta_inf[i, 0] = properties.eau_prop('betaf', (
                                TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2);  # [1/K]  Volumetric expansivity (beta)
                    propriHXefd_out.Cp_inf[i, 0] = properties.eau_prop('Cpf', (
                                TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2);  # [J/(kg K)] Heat capacity
                    propriHXefd_out.alpha_inf[i, 0] = propriHXefd_out.k_inf[i, 0] / (
                                propriHXefd_out.rho_inf[i, 0] * propriHXefd_out.Cp_inf[
                            i, 0]);  # [m^2/s] Thermal diffusivity

                    # Calul du coefficient de convection naturelle entre l'�changeur EFD et l'eau du r�servoir
                    #
                    h_convHXefd_out[i, 0] = Churchill_Chu(paramHXefd.Do, propriHXefd_out.rho_inf[i, 0],
                                                          propriHXefd_out.Pr_inf[i, 0], propriHXefd_out.k_inf[i, 0],
                                                          propriHXefd_out.beta_inf[i, 0],
                                                          propriHXefd_out.alpha_inf[i, 0], propriHXefd_out.mu_inf[i, 0],
                                                          HXefdresults.Ts_outk[i, m], Res2results.Tk[
                                                              i, 0]);  # [W/(m^2 K)] Convection heat transfer coefficient

                    # Calcul des propri�t�s de l'eau dans l'�changeur EFD
                    #
                    propriHXefd_in.rho_inf[i, 0] = properties.eau_prop('vf', TmoyEFD[i, 0]) ** -1;  # [kg/m3] Density
                    propriHXefd_in.k_inf[i, 0] = properties.eau_prop('kf',
                                                                     TmoyEFD[i, 0]);  # [W/(m K)] Thermal conductivity
                    propriHXefd_in.mu_inf[i, 0] = properties.eau_prop('muf',
                                                                      TmoyEFD[i, 0]);  # [Ns/m^2] Dynamic viscosity
                    propriHXefd_in.Pr_inf[i, 0] = properties.eau_prop('Prf', TmoyEFD[i, 0]);  # [-] Prandtl Number
                    propriHXefd_in.beta_inf[i, 0] = properties.eau_prop('betaf', TmoyEFD[
                        i, 0]);  # [1/K]  Volumetric expansivity (beta)
                    propriHXefd_in.Cp_inf[i, 0] = properties.eau_prop('Cpf', TmoyEFD[i, 0]);  # [J/(kg K)] Heat capacity
                    propriHXefd_in.alpha_inf[i, 0] = propriHXefd_in.k_inf[i, 0] / (
                                propriHXefd_in.rho_inf[i, 0] * propriHXefd_in.Cp_inf[
                            i, 0]);  # [m^2/s] Thermal diffusivity

                    # Calul du coefficient de convection forc�e � l'int�rieur de l'�changeur EFD
                    #
                    h_convHXefd_in[i, 0] = Dittus_Boelter_mod(paramHXefd.Di, paramHXefd.D, paramHXefd.pitch,
                                                              paramHXefd.m_in, propriHXefd_in.Cp_inf[i, 0],
                                                              propriHXefd_in.k_inf[i, 0], propriHXefd_in.mu_inf[i, 0],
                                                              propriHXefd_in.Pr_inf[
                                                                  i, 0]);  # [W/(m^2 K)] Convection heat transfer coefficient

                    # Calcul des r�sistances des flux de chaleur entre l'�changeur ECD et le r�servoir
                    #
                    R_convHXefd_in[i, 0] = (h_convHXefd_in[
                                                i, 0] * np.pi * paramHXefd.Di * paramHXefd.L) ** -1;  # [K/W] R�sistance convection interne de l'�changeur EFD
                    R_convHXefd_out[i, 0] = (h_convHXefd_out[
                                                 i, 0] * np.pi * paramHXefd.Do * paramHXefd.L) ** -1;  # [K/W] R�sistance convection externe de l'�changeur EFD
                    R_HXefd[i, 0] = R_convHXefd_in[i, 0] + R_condHXecd_wall + R_convHXefd_out[
                        i, 0];  # [K/W] R�sistance totale de l'�changeur EFD

                    # Recalcul de la temp�rature de sortie de l'�changeur pour le noeud du r�servoir i
                    #
                    HXefdresults.Toutk[i, m] = Res2results.Tk[i, m] + (
                                HXefdresults.Tink[i, m] - Res2results.Tk[i, m]) * np.exp(
                        (-R_HXefd[i, 0] * paramHXefd.m_in * propriHXefd_in.Cp_inf[i, 0]) ** -1)

                    # Calcul de DeltaTlm (Méthode LMTD)
                    deltaTo[i, 0] = Res2results.Tk[i, m] - HXefdresults.Toutk[i, m]
                    deltaTi[i, 0] = Res2results.Tk[i, m] - HXefdresults.Tink[i, m]
                    if deltaTi[i, 0] == deltaTo[i, 0]:
                        deltaTlm[i, 0] = 0
                    else:
                        deltaTlm[i, 0] = (deltaTo[i, 0] - deltaTi[i, 0]) / np.log(deltaTo[i, 0] / deltaTi[i, 0])

                    # Calcul du flux de chaleur de l'échangeur EFD
                    q_EFD[i, 0] = R_HXefd[i, 0] ** -1 * deltaTlm[i, 0]  # [W]

                    # Recalcul de la temp�rature de surface externe de l'�changeur pour le noeud du r�servoir i
                    HXefdresults.Ts_outk[i, m] = TmoyEFD[i, 0] + q_EFD[i, 0] * (
                                R_convHXefd_in[i, 0] + R_condHXecd_wall);  # [K]

                    # Calcul de l'erreur entre cette it�ration et l'it�ration pr�c�dente
                    #
                    iterToutHX = iterToutHX + 1;
                    erreurTout[iterToutHX, 0] = np.sqrt((HXefdresults.Toutk[i, m] - Toutguess) ** 2);

                # On assigne la temp�rature calcul�e en sortie � la temp�rature d'entr�e du prochain noeud, i-1
                #
                if i == paramHXefd.outlet:
                    break
                else:
                    HXefdresults.Tink[i - 1, m] = HXefdresults.Toutk[i, m];

        # Calcul des coefficients de la matrice D
        #
        d1 = propriFluid_EC.rho_inf[1:-1] * propriFluid_EC.Cp_inf[1:-1] * paramRes2.dV / dt;
        d3 = paramRes2.m_in * propriFluid_EC.Cp_inf[1:-1];
        d4 = np.pi * paramRes2.Di * paramRes2.dy / Rpp_pertesside;
        d5 = paramRes2.As / Rpp_pertestopbott;

        # ____________________________________________________
        # CALCUL DU CHAMPS DE TEMP�RATURE EC, ALGORITHME TDMA
        # ����������������������������������������������������
        A = np.ones([paramRes2.nb_y - 1, 1]) * (
                    -d2 - d3[0:-1] * M_debit_A);  # Vecteur de la diagonale sup�rieure de longueur N-1
        D = np.ones([paramRes2.nb_y, 1]) * (d1 + d2 * (
                    2 - M_topbott) + d3 * M_debit_D + d4 + d5 * M_topbott);  # Vecteur de la diagonale centrale de longueur N
        C1 = np.ones([paramRes2.nb_y, 1]) * (
                    d3 * M_debit_C * paramRes2.Tink_EC + d4 * data.Tairk + d5 * M_topbott * data.Tairk - q_EFD);  # Vecteur des constantes C (ind�pendante du temps m)
        C = np.ones([paramRes2.nb_y, 1]) * (
                    d1 * Res2results.Tk[:, [m - 1]]) + C1;  # Vecteur des constantes C (d�pendante du temps m)

        # R�solution dans l'espace des temp�ratures de l'eau du r�servoir EC
        #
        Res2results.Tk[:, [m]] = TDMA(A, B, C, D);

        # Calcul du crit�re de convergence
        # 
        iterTres2 = iterTres2 + 1;
        erreurTres2[iterTres2, 0] = np.sqrt(sum((Res2results.Tk[:, [m]] - Tguess) ** 2)) / paramRes2.nb_y;
        print('\n CALCULATION COMPLETED - PAS DE TEMPS # {}/{} - ITERATION # {} \n'.format(m, nb_t, iterTres2))

        if iterTres2 >= limit_iter_Tres2:
            print('Temperature field not converged - ResiduTres2 = {}'.format(erreurTres2[iterTres2, 0]));
            break
# %%
# Conversion des r�sultats de K � �C
#
HXefdresults.Tin = HXefdresults.Tink - 273.15;
HXefdresults.Tout = HXefdresults.Toutk - 273.15;
Res2results.T = Res2results.Tk - 273.15;

# Affichage des conditions d'op�ration
#
print('La temp�rature de la pi�ce est de {0}�C.\n'.format(data.Tairk - 273.15))
print('La temp�rature initial du r�servoir est {0}�C.\n'.format(paramRes2.Ti))
print('La temp�rature du retour d\'eau chaude est {0}�C \n'.format(paramRes2.Tin_EC))
print('avec un d�bit de {0:.2f} kg/s \n'.format(paramRes2.m_in))
print('La temp�rature de l\'eau chaude domestique � l\'entr�e est {0}�C.\n'.format(paramHXefd.Tin))
print('avec un d�bit de {0:.2f} kg/s \n'.format(paramHXefd.m_in))

# %%
# �������������������������������������������������������������������������%
# ------------------------ Impression des figures -------------------------%
# _________________________________________________________________________%

from matplotlib import pyplot as plt

fig1, (ax1, ax11) = plt.subplots(2, 1, num='Montly Loads')
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,nb_t-1]), label='$\itt$={} h'.format(temps[nb_t-1]/3600))
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,round(nb_t/2)]), label='$\itt$={} h'.format(temps[round(nb_t/2)]/3600))
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,round(nb_t/3.25)]), label='$\itt$={} h'.format(temps[round(nb_t/3.25)]/3600))
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,round(nb_t/5)]), label='$\itt$={} h'.format(temps[round(nb_t/5)]/3600))
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,round(nb_t/10)]), label='$\itt$={} h'.format(temps[round(nb_t/10)]/3600))
ax1.plot(paramRes2.y_num,np.flipud(Res2results.T[:,0]), label='$\itt$={} h'.format(temps[0]/3600))

ax1.set_xlabel('$\ity$ [m]', color='k', fontname='Lucida Bright',fontsize=9);
ax1.set_ylabel('$\itT_{SHW}$ [°C]', color='k', fontname='Lucida Bright',fontsize=9);

ax11.plot(temps/3600, Res2results.T[0,:], label='$\ity$={0:.2f} m'.format(paramRes2.y_num[paramRes2.nb_y-1]))        # Haut du réservoir
ax11.plot(temps/3600, Res2results.T[round(paramRes2.nb_y/2),:], label='$\ity$={0:.2f} m'.format(paramRes2.y_num[round(paramRes2.nb_y/2)])) # Mi hauteur
ax11.plot(temps/3600, Res2results.T[paramRes2.nb_y-1,:], label='$\ity$={0:.2f} m'.format(paramRes2.y_num[0]))    #Bas du réservoir
ax11.plot(temps/3600, HXefdresults.Tout[0,:], label='$\itDCW_{out}$')         # sortie DHW HX

ax11.set_xlabel('Time [h]', color='k', fontname='Lucida Bright',fontsize=9);
ax11.set_ylabel('$\itT_{SHW}$ or $\itT_{DCW_{out}}$ [°C]', color='k', fontname='Lucida Bright',fontsize=9);

ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=3, prop = {'family':'Lucida Bright', 'size':9})
ax11.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20), ncol=2, prop = {'family':'Lucida Bright', 'size':9})

fig1.set_figheight(fig1.get_figheight()*1.5)
fig1.set_figwidth(fig1.get_figwidth())

fig1.tight_layout()  # otherwise the right y-label is slightly clipped


fig1.savefig('Tank Simulation results.jpg', format='JPG', dpi=300)

# %%
# �������������������������������������������������������������������������#
# ------------------- Analyse des �chelles de temps -----------------------#
# _________________________________________________________________________#

print('Volume interne du r�servoir = {0:.2f} [m3]'.format(paramRes2.Vres))
print('Volume occup� par l''�changeur EFD = {0:.2f} [m3]'.format(paramHXefd.Aext * paramHXefd.Ltot))
Volume_eau_EC = paramRes2.Vres - paramHXefd.Aext * paramHXefd.Ltot;
print('Volume d''eau EC = {0:.2f} [m3]'.format(Volume_eau_EC))
Volume_eau_EFD = paramHXefd.As * paramHXefd.Ltot;
print('Volume d''eau dans l''�changeur EFD = {0:.2f} [m3]'.format(Volume_eau_EFD))
Volume_cuivre = (paramHXefd.Aext - paramHXefd.As) * paramHXefd.Ltot;
print('Volume de Cuivre dans l''�changeur EFD = {0:.2f} [m3]'.format(Volume_cuivre))

print('Masse thermique d''eau EC = {0:.2f} [kJ/K]'.format(data.rho_EC * Volume_eau_EC * data.Cp_EC * 1e-3))
print('Masse thermique d''eau EFD � 300[K] = {0:.2f} [kJ/K]'.format(data.rho_EC * Volume_eau_EFD * data.Cp_EC * 1e-3))
print('Masse thermique du Cuivre de l''�changeur EFD = {0:.2f} [kJ/K]'.format(
    data.rho_cu * Volume_cuivre * data.cp_cu * 1e-3))

print('Temps de transport SHW = {0:.2f} [min]'.format(data.rho_EC * Volume_eau_EC / paramRes2.m_in / 60))
print('Temps de transport EFD = {0:.2f} [s]'.format(data.rho_EC * Volume_eau_EFD / paramHXefd.m_in))

print('Temps caract�ristique du r�servoir = {0:.2f} [h]'.format((paramRes2.Di / 2) ** 2 / (4 * data.alpha_EC) / 3600))
