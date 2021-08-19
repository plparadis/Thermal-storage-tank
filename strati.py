"""
###########################################################################
####-- Heat storage tank with thermal stratification - Main routine --#####
#-------------------------------------------------------------------------#
#            Par Pierre-Luc Paradis et Marie-Hélène Talbot               #
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
from Dittus_Boelter_mod_func import Dittus_Boelter_mod
import streamlit as st
from matplotlib import pyplot as plt
import UnitConversionModule as conversion

TextColor = "#404041"
plt.rcParams['text.color'] = TextColor
plt.rcParams['axes.labelcolor'] = TextColor
plt.rcParams['xtick.color'] = TextColor
plt.rcParams['ytick.color'] = TextColor
fontname = 'Helvetica LT Std'

st.set_page_config(page_title="Thermal Storage Tank",
                   page_icon="static/favicon.ico",
                   layout="centered",
                   initial_sidebar_state="auto")

col1, col2 = st.columns((3, 1))
col1.title("Thermal Storage Tank")
with col2:
    col2.image("static/logo.png", use_column_width=True, output_format='PNG')

with st.expander("Tool description", expanded=True):
    st.markdown("A simple simulation tool to simulate stratified thermal storage tank as seen below.")
    st.image("static/Schematic.JPG", use_column_width=True, output_format='JPG')

    st.markdown("A heat balance on a slice of the tank gives")
    st.write(r'''
        $$
        \underbrace{(rho*c_p)_{SHW} * 
        \frac{\partial T_{SHW}}{\partial t}}_{\text{Transient term}}
        =
        \underbrace{\frac{\partial }{\partial y} \lbrack (k_{SHW}+\Delta k) \frac{\partial T_{SHW}}{\partial y} \rbrack}_{\text{1-D conduction term}} -
        \underbrace{\frac{(\dot{m}c_p)_{SHW}}{A_{SHW}}\frac{\partial T_{SHW}}{\partial y}}_{\text{Advection term}} +
        $$
        ''')
    st.write(r'''
        $$
        \underbrace{\frac{dq_{DCW} - dq_{loss}}{d\forall_{SHW}}}_{\text{Source terms}}     
        $$
        ''')

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
mixing = st.sidebar.checkbox("Handle Temperature inversion",value=False)  # Active/désactive l'Algorithme de mixing

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

data.k_EC = properties.eau_prop('kf', 300)  # [W/(m K)] Thermal conductivity de l'eau à 300[K]
data.rho_EC = properties.eau_prop('vf', 300) ** -1  # [kg/m3] Density de l'eau à 300[K]
data.Cp_EC = properties.eau_prop('Cpf', 300)  # [J/(kg K)] Heat capacity de l'eau à 300[K]
data.alpha_EC = data.k_EC / (data.rho_EC * data.Cp_EC)  # [m^2/s] Thermal diffusivity de l'eau à 300[K]

data.Tairk = conversion.FtoC(st.sidebar.number_input("Mechanical room Temp. [°F]", 40, 100,
                                     70, 2)) + 273.15  # [K] Température dans la salle mécanique
data.critere = 1e-6  # [-] Critère de convergence

# Parametres de la resolution temporelle
st.sidebar.markdown('** &#10112 Transient solution parameter **', unsafe_allow_html=False)
nb_t = st.sidebar.number_input("Number of time steps", 1, 3600, 10)  # [-] Nombre de pas de temps
dt = st.sidebar.number_input("Time step, [s]", 15, 3600, 900,
                             300)  # [s] Pas de temps ( 300 [s] --> 5 [min], 900 [s] --> 15 [min])
temps = np.arange(start=0, stop=nb_t * dt, step=dt).transpose()  # [s] Vecteur des pas de temps

# Reservoir #2 - EC
# Res2 - Geometry
st.sidebar.markdown('** &#10113 Tank **', unsafe_allow_html=False)
st.sidebar.markdown('Geometry')
paramRes2.De = st.sidebar.number_input("Tank outer diameter, [in]", 12, 60, 24,
                                       6) * 0.0254  # [m] Diamètre extérieure du réservoir
paramRes2.ye = st.sidebar.number_input("Tank height, [in]", 12, 96, 44,
                                       6) * 0.0254  # [m] Hauteur extérieure du réservoir
paramRes2.t_side = st.sidebar.number_input("Tank wall tickness, [in]", 0.0625, 0.5, 0.18750, 0.0625,
                                           format="%.4f") * 0.0254  # [m] épaisseur des parois de côté
paramRes2.t_top = st.sidebar.number_input("Tank top wall tickness, [in]", 0.0625, 0.5, 0.125, 0.0625,
                                          format="%.4f") * 0.0254  # [m] épaisseur de la paroi du dessus
paramRes2.t_bottom = st.sidebar.number_input("Tank bottom wall tickness, [in]", 0.0625, 0.5, 0.375, 0.0625,
                                             format="%.4f") * 0.0254  # [m] épaisseur de la paroi du dessous
paramRes2.Di = paramRes2.De - 2 * paramRes2.t_side  # [m] Diamètre intérieure du réservoir
paramRes2.yi = paramRes2.ye - paramRes2.t_top - paramRes2.t_bottom  # [m] Hauteur intérieure du réservoir
paramRes2.As = np.pi * paramRes2.Di ** 2 / 4  # [m^2] Aire de la section interne du réservoir
paramRes2.As_ext = np.pi * paramRes2.De ** 2 / 4  # [m^2] Aire de la section externe du réservoir
paramRes2.Vres = paramRes2.As * paramRes2.yi  # [m^3] Volume interne du réservoir
paramRes2.Vres_ext = paramRes2.As_ext * paramRes2.ye  # [m^3] Volume externe du réservoir
paramRes2.Vss_tank = paramRes2.Vres_ext - paramRes2.Vres  # [m^3] Volume d'acier inoxydable du réservoir
paramRes2.t_iso = st.sidebar.number_input("Tank insulation thickness, [in]", 0.25, 4., 2., 0.25,
                                          format="%.2f") * 0.0254  # [m] épaisseur de l'isolant (2 pouces)

# Res2 - Discretization
st.sidebar.markdown('Discretization')
paramRes2.nb_y = st.sidebar.number_input("Number of vertical Tank nodes", 1, 15, 5, 1)  # [noeuds] Nombre de noeuds du reservoir de stockage (arbritraire)
paramRes2.dy = paramRes2.yi / paramRes2.nb_y  # [m] Distance entre les noeuds (selon l'axe y)
paramRes2.y_num = np.arange(paramRes2.dy / 2, paramRes2.yi,
                            paramRes2.dy)  # [m] Vecteur de position de chacun des noeuds

# Res2 - Flow parameters
st.sidebar.markdown('Flow parameters')
paramRes2.m_in = st.sidebar.number_input("HW flow, [gpm]", 0., 50., 2.5, 0.5)/60*3.78541/1000*data.rho_EC  # [kg/s] Débit de la pompe -->60sec/hr | 3.78541liters/gal | 1000liters/m3
paramRes2.intlet_EC = paramRes2.nb_y  # [-] Numéro du noeud entrée (nb_y --> bas du réservoir)
paramRes2.outlet_EC = 0  # [-] Numéro du noeud sortie (0 --> haut du reservoir)
paramRes2.Ti = conversion.FtoC(st.sidebar.number_input("initial tank temperature, [°F]", 40, 160, 140, 10))  # [°C] Température initiale du réservoir 2
paramRes2.Tik = paramRes2.Ti + 273.15  # [K] Conversion
paramRes2.Tin_EC = conversion.FtoC(st.sidebar.number_input("HWR temperature, [°F]",60, 160, 120, 10))  # [°C] Température de retour du réseau de chauffage (Température d'entrée EC)
paramRes2.Tink_EC = paramRes2.Tin_EC + 273.15  # [K] Conversion

# Échangeur EFD
st.sidebar.markdown('** &#10114 Heat Exchanger **', unsafe_allow_html=False)
# HX EFD - Geometry
st.sidebar.markdown('Geometry')
paramHXefd.Do = st.sidebar.number_input("Coil pipe outside diameter, [in]", 0., 3., 0.75, 0.25) * 0.0254  # [m] Diamètre extérieur du tuyau (3/4 pouces)
paramHXefd.t_pipe = st.sidebar.number_input("Coil pipe wall thickness, [in]", 0., 0.125, 0.0486, 0.0001) * 0.0254  # [m] épaisseur de la paroi du tuyau
paramHXefd.Di = paramHXefd.Do - 2 * paramHXefd.t_pipe  # [m] diamètre intérieur du tuyau
paramHXefd.Ltot = st.sidebar.number_input("Total length of the coiled pipe, [in]", 0, 2500, 496, 6) * .0254  # [m] Longueur du tuyau
paramHXefd.As = np.pi * paramHXefd.Di ** 2 / 4  # [m^2] Aire de la section interne de la conduite
paramHXefd.Aext = np.pi * paramHXefd.Do ** 2 / 4  # [m^2] Aire de la section externe de la conduite
paramHXefd.D = st.sidebar.number_input("Coil internal diameter, [in]", 0, 60, 20, 2) * 0.0254  # [m] Diamètre de la spirale de l'échangeur
paramHXefd.pitch = st.sidebar.number_input("Coil pitch, [in]", 0., 6., 2., 0.5) * 0.0254  # [m] Pitch de la spirale de l'échangeur

# HX EFD - Flow parameters
st.sidebar.markdown('Flow parameters')
paramHXefd.m_in = st.sidebar.number_input("HX flow, [gpm]", 0., 30., 0.5, 0.5)/60*3.78541/1000*data.rho_EC  # [kg/s] Débit de consommation d'eau chaude domestique
paramHXefd.intlet = paramRes2.nb_y  # [-] Numéro du noeud entrée (nb_y --> bas du réservoir)
paramHXefd.outlet = 0  # [-] Numéro du noeud sortie (0--> haut du reservoir)
paramHXefd.nb_y = paramHXefd.intlet - (
        paramHXefd.outlet - 1)  # [noeuds] Nombre de noeuds où l'échangeur EFD est présent
paramHXefd.Tin = conversion.FtoC(st.sidebar.number_input("HX inlet temperature, [°F]", 40, 80, 40, 5))  # [°C] Température de l'eau potable de la ville (Température d'entrée EFD)
paramHXefd.Tink = paramHXefd.Tin + 273.15  # [K] Conversion

# HX EFD - Discretization
paramHXefd.L = paramHXefd.Ltot / paramHXefd.nb_y  # [m] Longueur de tuyau de l'échangeur EFD dans un noeud j du réservoir 2

# Res2
paramRes2.Veau = paramRes2.Vres - (
        paramHXefd.Aext * paramHXefd.Ltot)  # [m^3] Volume d'eau dans le réservoir --> Calcul Volume EC = (Volume interne du r�servoir - Volume des �changeurs de chaleurs)
paramRes2.dV = paramRes2.Veau / paramRes2.nb_y  # [m^3] Volume d'un noeud du réservoir de stockage
paramRes2.deltak = data.kss * (
        np.pi * paramRes2.De ** 2 / 4 - paramRes2.As) / paramRes2.As  # De-stratification due au transfert de chaleur par conduction axiale dans la paroi du r�servoir

# Bilan d'énergie sur volume de contrôle j du réservoir

# Masques
M_topbott = np.zeros([paramRes2.nb_y, 1])
M_topbott[0] = 1  # Identification du noeud du haut (Condition aux frontières)
M_topbott[-1] = 1  # Identification du noeud du bas (Condition aux frontières)

M_debit_D = np.zeros([paramRes2.nb_y, 1])
M_debit_D[paramRes2.outlet_EC:paramRes2.intlet_EC] = 1  # Identification des noeuds de débit EC pour la diagonale D (Matrice des coefficients)

M_debit_A = np.zeros([paramRes2.nb_y - 1, 1])
M_debit_A[paramRes2.outlet_EC:paramRes2.intlet_EC - 1] = 1  # Identification des noeuds de débit EC pour la diagonale A (Matrice des coefficients)

M_debit_C = np.zeros([paramRes2.nb_y, 1])
M_debit_C[paramRes2.intlet_EC - 1] = 1  # Identification des noeuds de débit EC pour la diagonale C (Matrice des coefficients)

M_HXefd = np.zeros([paramRes2.nb_y, 1])
M_HXefd[paramHXefd.outlet:paramHXefd.intlet] = 1  # Localisation de l'échangeur EFD

# Initialisation du vecteur de solution du RES2 - Température EC
Res2results.Tk = np.zeros([paramRes2.nb_y, nb_t]) + paramRes2.Tik  # [K] Initialisation de la matrice de température EC (ligne = #noeud et colonne = pas de temps)

# Initialisation des vecteurs de solutions de l'échangeur EFD
HXefdresults.Tink = np.zeros([paramRes2.nb_y, nb_t])  # [K] Initialisation de la température dans l'échangeur EFD à l'entrée de chaque noeud du réservoir
HXefdresults.Tink[:, [0]] = Res2results.Tk[:, [0]] * M_HXefd  # [K] Initialisation de la température dans l'échangeur EFD (Conditions initiales --> Time step #1)
HXefdresults.Tink[paramHXefd.intlet - 1, np.arange(1, nb_t, 1)] = paramHXefd.Tink  # [K] Initialisation de la température à l'entrée de l'échangeur d'EFD (Boundary condition)

HXefdresults.Toutk = np.copy(HXefdresults.Tink)  # [K] Initialisation de la température dans l'échangeur EFD à la sortie de chaque noeud du réservoir

HXefdresults.Ts_outk = np.zeros([paramRes2.nb_y, nb_t])  # [K] Initialisation de la température de surface externe de l'échangeur d'EFD
HXefdresults.Ts_outk[:, [1]] = Res2results.Tk[:, [0]] * M_HXefd  # [K] Initialisation de la température de surface externe de l'échangeur d'EFD (Conditions initiales --> Time step #1)

# Initialisation des vecteurs pour le calcul des pertes thermiques du réservoir
h_convEC = np.zeros([paramRes2.nb_y + 2, 1])  # [W/(m^2) K] Convection heat transfer coefficient in EC
h_convair = np.zeros([paramRes2.nb_y + 2, 1])  # [W/(m^2) K] Convection heat transfer coefficient in ambiant air
Rpp_pertestopbott = np.ones([paramRes2.nb_y, 1])  # [K (m^2)/W] Résistance thermique résultante pour les extrémitées
qpp_pertes = np.zeros([paramRes2.nb_y + 2, 1])  # [W/m^2] Flux de chaleur total des pertes thermiques du réservoir

# Initialisation des temp�ratures de surfaces des parois du réservoir
Tsk_in = np.zeros([paramRes2.nb_y + 2, 1]) + paramRes2.Tik - 4  # [K] Surface interne
Tsk_out = np.zeros([paramRes2.nb_y + 2, 1]) + data.Tairk - 4  # [K] Surface externe

# Initialisation des vecteurs pour le calcul de l'échangeur EFD
h_convHXefd_in = np.zeros(
    [paramRes2.nb_y, 1])  # [W/(m^2) K] Convection heat transfer coefficient inside échangeur EFD (convection forcée)
h_convHXefd_out = np.zeros([paramRes2.nb_y,
                            1])  # [W/(m^2) K] Convection heat transfer coefficient outside échangeur EFD (convection naturelle)
TmoyEFD = np.zeros([paramRes2.nb_y, 1])  # [K] Température moyenne entre l'entrée et la sortie (échangeur EFD)
deltaTi = np.zeros([paramRes2.nb_y, 1])  # [K] Pincement à l'entrée de l'échangeur EFD
deltaTo = np.zeros([paramRes2.nb_y, 1])  # [K] Pincement à la sortie de l'échangeur EFD
R_convHXefd_in = np.zeros([paramRes2.nb_y, 1])  # [K/W] Résistance convection interne de l'échangeur EFD
R_convHXefd_out = np.zeros([paramRes2.nb_y, 1])  # [K/W] Résistance convection externe de l'échangeur EFD
R_HXefd = np.zeros([paramRes2.nb_y, 1])  # [K/W] Résistance totale de l'échangeur EFD
deltaTlm = np.zeros([paramRes2.nb_y, 1])  # Log Mean Temperature Difference

# Initialisation des flux de chaleur des échangeurs
q_EFD = np.zeros([paramRes2.nb_y, 1])  # [W] Flux de chaleur total échangeur EFD

# Calcul des Termes constant dans le temps
Rpp_condtanktop = paramRes2.t_top / data.kss  # [K m^2/W] Calcul de la résistance thermique par unité de surface au travers de la paroi du dessus du réservoir
Rpp_condtankbott = paramRes2.t_bottom / data.kss  # [K m^2/W] Calcul de la résistance thermique par unité de surface au travers de la paroi du dessous du réservoir
Rpp_condtankside = paramRes2.t_side / data.kss  # [K m^2/W] Calcul de la résistance thermique par unité de surface au travers des parois de côté du réservoir
Rpp_condiso = paramRes2.t_iso / data.kiso  # [K m^2/W] Calcul de la résistance thermique par unité de surface au travers de la couche d'isolant

R_condHXecd_wall = np.log(paramHXefd.Do / paramHXefd.Di) / (
        2 * np.pi * data.kcu * paramHXefd.L)  # [K/W] Calcul de la résistance thermique totale au travers de la paroi de l'échangeur EFD

d2 = paramRes2.As * (
        data.kss + paramRes2.deltak) / paramRes2.dy  # Terme constant de la diagonale D (Matrice des coefficients)
B = np.ones([paramRes2.nb_y - 1, 1]) * (
    -d2)  # Vecteur de la diagonale inférieure B de longueur (Nb_y-1) (Matrice des coefficients)

# Initalisation des propri
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

# RÉSOLUTION EN TEMPS
for m in range(1, nb_t, 1):
    iterTres2 = 0
    limit_iter_Tres2 = 25
    erreurTres2 = np.ones([limit_iter_Tres2, 1])
    Res2results.Tk[:, [m]] = Res2results.Tk[:, [
                                                   m - 1]]  # Initilaisation de la solution du pas de temps courant à partir de la solution du pas de temps précédent

    # Recalcul des j températures du réservoir de stockage jusqu'à convergence pour le pas de temps courant, m
    while erreurTres2[iterTres2] > data.critere:

        Tguess = Res2results.Tk[:, [
                                       m]]  # Enregistrement de la dernière solution calculée (permet de calculer le critère de convergence sur la température EC)

        # PERTES THERMIQUES DU RÉSERVOIR
        for i in range(
                paramRes2.nb_y + 2):  # La boucle se fait sur nb_y+2 afin de calculer les pertes thermiques des parois sur les nb_y noeuds du réservoir ainsi que pour le dessus et le dessous du réservoir (conditions aux frontières)
            # Structure "if" pour faire corresponde le bon noeud du réservoir pour le calcul des propriétés pour le dessus et le dessous du réservoir (Conditions aux frontières) à l'aide de l'indice j
            if i == 0:
                j = 0
            elif i > 1 and i < paramRes2.nb_y + 1:
                j = i - 1
            elif i == paramRes2.nb_y + 1:
                j = i - 2

            # Calcul des propriétés de l'eau chaude dans le réservoir
            propriFluid_EC.rho_inf[i, 0] = properties.eau_prop('vf', (
                    Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2) ** -1  # [kg/m3] Density
            propriFluid_EC.k_inf[i, 0] = properties.eau_prop('kf', (
                    Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2)  # [W/(m K)] Thermal conductivity
            propriFluid_EC.mu_inf[i, 0] = properties.eau_prop('muf', (
                    Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2)  # [Ns/m^2] Dynamic viscosity
            propriFluid_EC.Pr_inf[i, 0] = properties.eau_prop('Prf', (
                    Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2)  # [-] Prandtl Number
            propriFluid_EC.beta_inf[i, 0] = properties.eau_prop('betaf', (
                    Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2)  # [1/K]  Volumetric expansivity (beta)
            propriFluid_EC.Cp_inf[i, 0] = properties.eau_prop('Cpf', (
                    Res2results.Tk[j, m] + Tsk_in[i, 0]) / 2)  # [J/(kg K)] Heat capacity
            propriFluid_EC.alpha_inf[i, 0] = propriFluid_EC.k_inf[i, 0] / propriFluid_EC.rho_inf[i, 0] / \
                                             propriFluid_EC.Cp_inf[i, 0]  # [m^2/s] Thermal diffusivity

            # Calul des coefficients de convection entre l'eau dans le réservoir et les parois
            if i == 0:  # Pour le haut du réservoir
                h_convEC[i, 0] = Churchill_Chu_plate(paramRes2.Di, propriFluid_EC.rho_inf[i, 0],
                                                     propriFluid_EC.Pr_inf[i, 0], propriFluid_EC.k_inf[i, 0],
                                                     propriFluid_EC.beta_inf[i, 0],
                                                     propriFluid_EC.alpha_inf[i, 0], propriFluid_EC.mu_inf[i, 0],
                                                     Tsk_in[i, 0], Res2results.Tk[j, 1],
                                                     3)  # [W/(m^2) K] Convection heat transfer coefficient in EC
            elif i == (paramRes2.nb_y + 1):  # Pour le bas du réservoir
                h_convEC[i, 0] = Churchill_Chu_plate(paramRes2.Di, propriFluid_EC.rho_inf[i, 0],
                                                     propriFluid_EC.Pr_inf[i, 0], propriFluid_EC.k_inf[i, 0],
                                                     propriFluid_EC.beta_inf[i, 0],
                                                     propriFluid_EC.alpha_inf[i, 0], propriFluid_EC.mu_inf[i, 0],
                                                     Tsk_in[i, 0], Res2results.Tk[j, 1],
                                                     2)  # [W/(m^2) K] Convection heat transfer coefficient in EC
            else:  # Pour les côtés
                h_convEC[i, 0] = Churchill_Chu_plate(paramRes2.yi, propriFluid_EC.rho_inf[i, 0],
                                                     propriFluid_EC.Pr_inf[i, 0], propriFluid_EC.k_inf[i, 0],
                                                     propriFluid_EC.beta_inf[i, 0],
                                                     propriFluid_EC.alpha_inf[i, 0], propriFluid_EC.mu_inf[i, 0],
                                                     Tsk_in[i, 0], Res2results.Tk[j, 1],
                                                     1)  # [W/(m^2) K] Convection heat transfer coefficient in EC

            # Calcul des propriétés de l'air ambiant qui échange avec le réservoir
            propriAir.rho[i, 0] = properties.air_prop('rho', (data.Tairk + Tsk_out[i, 0]) / 2)  # [kg/m3] Density
            propriAir.k[i, 0] = properties.air_prop('k',
                                                    (data.Tairk + Tsk_out[i, 0]) / 2)  # [W/(m K)] Thermal conductivity
            propriAir.mu[i, 0] = properties.air_prop('mu',
                                                     (data.Tairk + Tsk_out[i, 0]) / 2)  # [Ns/m^2] Dynamic viscosity
            propriAir.Pr[i, 0] = properties.air_prop('Pr', (data.Tairk + Tsk_out[i, 0]) / 2)  # [-] Prandtl number
            propriAir.beta[i, 0] = 1 / ((data.Tairk + Tsk_out[i, 0]) / 2)  # [1/K]  Volumetric expansivity (beta)
            propriAir.alpha[i, 0] = properties.air_prop('al', (
                    data.Tairk + Tsk_out[i, 0]) / 2)  # [m^2/s] Thermal diffusivity

            # Calul des coefficients de convection entre les parois du réservoir et l'air ambiant
            if i == 0:  # Pour le dessus du réservoir
                h_convair[i, 0] = Churchill_Chu_plate(paramRes2.De, propriAir.rho[i, 0],
                                                      propriAir.Pr[i, 0], propriAir.k[i, 0], propriAir.beta[i, 0],
                                                      propriAir.alpha[i, 0], propriAir.mu[i, 0], Tsk_out[i, 0],
                                                      data.Tairk,
                                                      2)  # [W/(m^2) K] Convection heat transfer coefficient
            elif i == (paramRes2.nb_y + 1):  # Pour le dessous du réservoir
                h_convair[i, 0] = Churchill_Chu_plate(paramRes2.De, propriAir.rho[i, 0],
                                                      propriAir.Pr[i, 0], propriAir.k[i, 0], propriAir.beta[i, 0],
                                                      propriAir.alpha[i, 0], propriAir.mu[i, 0], Tsk_out[i, 0],
                                                      data.Tairk,
                                                      3)  # [W/(m^2) K] Convection heat transfer coefficient
            else:  # Pour les côtés
                h_convair[i, 0] = Churchill_Chu_plate(paramRes2.ye, propriAir.rho[i, 0],
                                                      propriAir.Pr[i, 0], propriAir.k[i, 0], propriAir.beta[i, 0],
                                                      propriAir.alpha[i, 0], propriAir.mu[i, 0], Tsk_out[i, 0],
                                                      data.Tairk,
                                                      1)  # [W/(m^2) K] Convection heat transfer coefficient

        # Calcul des résistances des flux de pertes thermiques
        Rpp_convEC = 1 / h_convEC  # [K (m^2)/W] Résistance thermique par convection naturelle dans l'eau chaude
        Rpp_convair = 1 / h_convair  # [K (m^2)/W] Résistance thermique par convection naturelle dans l'air ambiant
        Rpp_pertestopbott[0] = Rpp_convEC[0] + Rpp_condtanktop + Rpp_condiso + Rpp_convair[
            0]  # [K (m^2)/W] Résistance thermique résultante pour le dessus
        Rpp_pertestopbott[-1] = Rpp_convEC[-1] + Rpp_condtankbott + Rpp_condiso + Rpp_convair[
            -1]  # [K (m^2)/W] Résistance thermique résultante pour le dessous
        Rpp_pertesside = Rpp_convEC[1:-1] + Rpp_condtankside + Rpp_condiso + Rpp_convair[
                                                                             1:-1]  # [K (m^2)/W] Résistance thermique résultante pour les côtés
        Rpp_pertes = np.vstack((Rpp_pertestopbott[0], Rpp_pertesside, Rpp_pertestopbott[
            -1]))  # [K (m^2)/W] Vecteur des résistances thermiques totales pour le dessus, dessous et côtés

        for i in range(paramRes2.nb_y + 2):
            # Structure if pour faire corresponde le bon noeud du réservoir pour le calcul de propriétés pour le dessus et le dessous du réservoir (Conditions aux frontières) à l'aide de l'indice j
            if i == 0:
                j = 0
            elif i > 1 and i < paramRes2.nb_y + 1:
                j = i - 1
            elif i == paramRes2.nb_y + 1:
                j = i - 2

            # Calcul des pertes thermiques du réservoir EC
            qpp_pertes[i, 0] = (Res2results.Tk[j, m] - data.Tairk) / Rpp_pertes[i, 0]

            # Calcul des températures de surface à l'intérieur et à l'extérieur des parois du réservoir EC
            Tsk_in[i, 0] = Res2results.Tk[j, m] - qpp_pertes[i, 0] * Rpp_convEC[i, 0]
            Tsk_out[i, 0] = data.Tairk + qpp_pertes[i, 0] * Rpp_convair[i, 0]

        # ÉCHANGEUR EFD
        for ii in range(paramHXefd.intlet, paramHXefd.outlet, -1):
            i = ii - 1
            if paramHXefd.m_in > 0:
                if iterTres2 == 0:
                    # Hypothèse sur la température de sortie du noeud (nécessaire seulement pour la première itération sur le champs de température EC de chaque pas de temps)
                    if HXefdresults.Tink[i, m] > Res2results.Tk[i, m]:  # Si l'échangeur EFDS chauffe le réservoir
                        HXefdresults.Toutk[i, m] = HXefdresults.Tink[i, m] - 5
                    else:  # Si le réservoir chauffe l'échangeur EFD
                        HXefdresults.Toutk[i, m] = HXefdresults.Tink[i, m] + 5
                    # Hypothèse sur la température de surface externe de l'échangeur EFD
                    HXefdresults.Ts_outk[i, m] = (((HXefdresults.Tink[i, m] + HXefdresults.Toutk[i, m]) / 2) +
                                                  Res2results.Tk[
                                                      i, m]) / 2  # [K] On initilaise la température de surface avec la moyenne de la température moyenne entre l'entrée et la sortie et la température de l'eau dans le réservoir au noeud i
                iterToutHX = 0
                limit_iter_ToutHX = 50
                erreurTout = np.ones([limit_iter_ToutHX, 1])
                # Recalcul des j températures de sortie de l'échangeur jusqu'à convergence pour le temps donné m
                while erreurTout[iterToutHX] > data.critere:

                    Toutguess = HXefdresults.Toutk[
                        i, m]  # Enregistrement de la derniére solution de Temp�rature de sortie calculée (permet de valider la convergence sur ToutHX)
                    TmoyEFD[i, 0] = (HXefdresults.Tink[i, m] + HXefdresults.Toutk[
                        i, m]) / 2  # Calcul de la température moyenne entre la température d'entrée et la température de sortie de l'échangeur

                    # Calcul des propriétés de l'eau dans le réservoir EC
                    propriHXefd_out.rho_inf[i, 0] = properties.eau_prop('vf', (
                            TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2) ** -1  # [kg/m3] Density
                    propriHXefd_out.k_inf[i, 0] = properties.eau_prop('kf', (
                            TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2)  # [W/(m K)] Thermal conductivity
                    propriHXefd_out.mu_inf[i, 0] = properties.eau_prop('muf', (
                            TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2)  # [Ns/m^2] Dynamic viscosity
                    propriHXefd_out.Pr_inf[i, 0] = properties.eau_prop('Prf', (
                            TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2)  # [-] Prandtl Number
                    propriHXefd_out.beta_inf[i, 0] = properties.eau_prop('betaf', (
                            TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2)  # [1/K]  Volumetric expansivity (beta)
                    propriHXefd_out.Cp_inf[i, 0] = properties.eau_prop('Cpf', (
                            TmoyEFD[i, 0] + Res2results.Tk[i, m]) / 2)  # [J/(kg K)] Heat capacity
                    propriHXefd_out.alpha_inf[i, 0] = propriHXefd_out.k_inf[i, 0] / (
                            propriHXefd_out.rho_inf[i, 0] * propriHXefd_out.Cp_inf[
                        i, 0])  # [m^2/s] Thermal diffusivity

                    # Calul du coefficient de convection naturelle entre l'échangeur EFD et l'eau du réservoir
                    #
                    h_convHXefd_out[i, 0] = Churchill_Chu(paramHXefd.Do, propriHXefd_out.rho_inf[i, 0],
                                                          propriHXefd_out.Pr_inf[i, 0], propriHXefd_out.k_inf[i, 0],
                                                          propriHXefd_out.beta_inf[i, 0],
                                                          propriHXefd_out.alpha_inf[i, 0], propriHXefd_out.mu_inf[i, 0],
                                                          HXefdresults.Ts_outk[i, m], Res2results.Tk[
                                                              i, 0])  # [W/(m^2 K)] Convection heat transfer coefficient

                    # Calcul des propriétés de l'eau dans l'échangeur EFD
                    propriHXefd_in.rho_inf[i, 0] = properties.eau_prop('vf', TmoyEFD[i, 0]) ** -1  # [kg/m3] Density
                    propriHXefd_in.k_inf[i, 0] = properties.eau_prop('kf',
                                                                     TmoyEFD[i, 0])  # [W/(m K)] Thermal conductivity
                    propriHXefd_in.mu_inf[i, 0] = properties.eau_prop('muf',
                                                                      TmoyEFD[i, 0])  # [Ns/m^2] Dynamic viscosity
                    propriHXefd_in.Pr_inf[i, 0] = properties.eau_prop('Prf', TmoyEFD[i, 0])  # [-] Prandtl Number
                    propriHXefd_in.beta_inf[i, 0] = properties.eau_prop('betaf', TmoyEFD[
                        i, 0])  # [1/K]  Volumetric expansivity (beta)
                    propriHXefd_in.Cp_inf[i, 0] = properties.eau_prop('Cpf', TmoyEFD[i, 0])  # [J/(kg K)] Heat capacity
                    propriHXefd_in.alpha_inf[i, 0] = propriHXefd_in.k_inf[i, 0] / (
                            propriHXefd_in.rho_inf[i, 0] * propriHXefd_in.Cp_inf[
                        i, 0])  # [m^2/s] Thermal diffusivity

                    # Calul du coefficient de convection forcée à l'intérieur de l'échangeur EFD
                    h_convHXefd_in[i, 0] = Dittus_Boelter_mod(paramHXefd.Di, paramHXefd.D, paramHXefd.pitch,
                                                              paramHXefd.m_in, propriHXefd_in.Cp_inf[i, 0],
                                                              propriHXefd_in.k_inf[i, 0], propriHXefd_in.mu_inf[i, 0],
                                                              propriHXefd_in.Pr_inf[
                                                                  i, 0])  # [W/(m^2 K)] Convection heat transfer coefficient

                    # Calcul des résistances des flux de chaleur entre l'échangeur ECD et le réservoir
                    R_convHXefd_in[i, 0] = (h_convHXefd_in[
                                                i, 0] * np.pi * paramHXefd.Di * paramHXefd.L) ** -1  # [K/W] Résistance convection interne de l'échangeur EFD
                    R_convHXefd_out[i, 0] = (h_convHXefd_out[
                                                 i, 0] * np.pi * paramHXefd.Do * paramHXefd.L) ** -1  # [K/W] Résistance convection externe de l'échangeur EFD
                    R_HXefd[i, 0] = R_convHXefd_in[i, 0] + R_condHXecd_wall + R_convHXefd_out[
                        i, 0]  # [K/W] Résistance totale de l'échangeur EFD

                    # Recalcul de la température de sortie de l'échangeur pour le noeud du réservoir i
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

                    # Recalcul de la température de surface externe de l'échangeur pour le noeud du réservoir i
                    HXefdresults.Ts_outk[i, m] = TmoyEFD[i, 0] + q_EFD[i, 0] * (
                            R_convHXefd_in[i, 0] + R_condHXecd_wall)  # [K]

                    # Calcul de l'erreur entre cette itération et l'itération précédente
                    iterToutHX = iterToutHX + 1
                    erreurTout[iterToutHX, 0] = np.sqrt((HXefdresults.Toutk[i, m] - Toutguess) ** 2)
                    if iterToutHX >= limit_iter_ToutHX - 1:
                        print('Heat Exchanger Outlet Temperature not converged - ResiduToutHX = {}'.format(erreurTout[iterToutHX, 0]))
                    break
                # On assigne la température calculée en sortie à la température d'entrée du prochain noeud, i-1
                if i == paramHXefd.outlet:
                    break
                else:
                    HXefdresults.Tink[i - 1, m] = HXefdresults.Toutk[i, m]

        # Calcul des coefficients de la matrice D
        d1 = propriFluid_EC.rho_inf[1:-1] * propriFluid_EC.Cp_inf[1:-1] * paramRes2.dV / dt
        d3 = paramRes2.m_in * propriFluid_EC.Cp_inf[1:-1]
        d4 = np.pi * paramRes2.Di * paramRes2.dy / Rpp_pertesside
        d5 = paramRes2.As / Rpp_pertestopbott

        # CALCUL DU CHAMPS DE TEMPÉRATURE EC, ALGORITHME TDMA
        A = np.ones([paramRes2.nb_y - 1, 1]) * (
                -d2 - d3[0:-1] * M_debit_A)  # Vecteur de la diagonale supérieure de longueur N-1
        D = np.ones([paramRes2.nb_y, 1]) * (d1 + d2 * (
                2 - M_topbott) + d3 * M_debit_D + d4 + d5 * M_topbott)  # Vecteur de la diagonale centrale de longueur N
        C1 = np.ones([paramRes2.nb_y, 1]) * (
                d3 * M_debit_C * paramRes2.Tink_EC + d4 * data.Tairk + d5 * M_topbott * data.Tairk - q_EFD)  # Vecteur des constantes C (indépendante du temps m)
        C = np.ones([paramRes2.nb_y, 1]) * (
                d1 * Res2results.Tk[:, [m - 1]]) + C1  # Vecteur des constantes C (dépendante du temps m)

        # Résolution dans l'espace des températures de l'eau du réservoir EC
        Res2results.Tk[:, [m]] = TDMA(A, B, C, D)

        # Calcul du critère de convergence
        iterTres2 = iterTres2 + 1
        erreurTres2[iterTres2, 0] = np.sqrt(sum((Res2results.Tk[:, [m]] - Tguess) ** 2)) / paramRes2.nb_y
        print('\n CALCULATION COMPLETED - PAS DE TEMPS # {}/{} - ITERATION # {} \n'.format(m, nb_t, iterTres2))

        if iterTres2 >= limit_iter_Tres2-1:
            print('Temperature field not converged - ResiduTres2 = {}'.format(erreurTres2[iterTres2, 0]))
            break
    if mixing:
        MixingCounter = 1  # Initialisation du compteur associé à l'algorithme d'inversion
        Inversion = Res2results.Tk[0:-1, [m]] < Res2results.Tk[1::, [m]]  # identifie les noeuds avec inversion
        CritereInversion = np.sum(Inversion)  # Initialisation du critère d'inversion
        if CritereInversion > 0:
            LowerNode = np.where(Inversion)[0][-1] + 1  # index du Noeud au bas de l'inversion
            while CritereInversion > 0:
                Res2results.Tk[(LowerNode - MixingCounter):LowerNode+1, [m]] = np.mean(Res2results.Tk[(LowerNode - MixingCounter):LowerNode+1, [m]])
                Inversion = Res2results.Tk[0:-1, [m]] < Res2results.Tk[1::, [m]]
                CritereInversion = np.sum(Inversion)
                MixingCounter = MixingCounter + 1

# Conversion des résultats de K à °C
HXefdresults.Tin = HXefdresults.Tink - 273.15
HXefdresults.Tout = HXefdresults.Toutk - 273.15
Res2results.T = Res2results.Tk - 273.15


# Affichage des conditions d'opération
with st.expander("Operating Conditions", expanded=True):
    st.markdown('Temperature in the mechanical room: **{:,.0f}°F**.'.format(conversion.CtoF(data.Tairk - 273.15)))
    st.markdown('Initial tank temperature: **{:,.0f}°F**.'.format(conversion.CtoF(paramRes2.Ti)))
    st.markdown('Temperature of Hot Water Return (HWR): **{:,.0f}°F**.'.format(conversion.CtoF(paramRes2.Tin_EC)))
    st.markdown('Hot water flow: **{0:.2f} gpm**.'.format(paramRes2.m_in*60/3.78541*1000/data.rho_EC))
    st.markdown('Coil Heat Exchanger (HX) inlet temperature: **{:,.0f}°F.**'.format(conversion.CtoF(paramHXefd.Tin)))
    st.markdown('Coil Heat exchanger flow: **{0:.2f} gpm**.'.format(paramHXefd.m_in*60/3.78541*1000/data.rho_EC))

#  Impression des figures
with st.expander("Transient Results", expanded=True):
    fig1, (ax1, ax11) = plt.subplots(2, 1, num='Montly Loads')
    ax1.plot(conversion.mtoft(paramRes2.y_num), conversion.CtoF(np.flipud(Res2results.T[:, nb_t - 1])), label='$\itt$={0:.2f} h'.format(temps[nb_t - 1] / 3600))
    ax1.plot(conversion.mtoft(paramRes2.y_num), conversion.CtoF(np.flipud(Res2results.T[:, round(nb_t / 2)])),
             label='$\itt$={0:.2f} h'.format(temps[round(nb_t / 2)] / 3600))
    ax1.plot(conversion.mtoft(paramRes2.y_num), conversion.CtoF(np.flipud(Res2results.T[:, round(nb_t / 3.25)])),
             label='$\itt$={0:.2f} h'.format(temps[round(nb_t / 3.25)] / 3600))
    ax1.plot(conversion.mtoft(paramRes2.y_num), conversion.CtoF(np.flipud(Res2results.T[:, round(nb_t / 5)])),
             label='$\itt$={0:.2f} h'.format(temps[round(nb_t / 5)] / 3600))
    ax1.plot(conversion.mtoft(paramRes2.y_num), conversion.CtoF(np.flipud(Res2results.T[:, round(nb_t / 10)])),
             label='$\itt$={0:.2f} h'.format(temps[round(nb_t / 10)] / 3600))
    ax1.plot(conversion.mtoft(paramRes2.y_num), conversion.CtoF(np.flipud(Res2results.T[:, 0])), label='$\itt$={0:.2f} h'.format(temps[0] / 3600))
    ax1.set_xlabel('$\ity$ [ft]', color='k', fontname='Lucida Bright', fontsize=9)
    ax1.set_ylabel('$\itT_{HW}$ [°F]', color='k', fontname='Lucida Bright', fontsize=9)
    ax11.plot(temps / 3600, conversion.CtoF(Res2results.T[0, :]),
              label='$\ity$={0:.2f} ft'.format(conversion.mtoft(paramRes2.y_num[paramRes2.nb_y - 1])))  # Haut du réservoir
    ax11.plot(temps / 3600, conversion.CtoF(Res2results.T[round(paramRes2.nb_y / 2), :]),
              label='$\ity$={0:.2f} ft'.format(conversion.mtoft(paramRes2.y_num[round(paramRes2.nb_y / 2)])))  # Mi hauteur
    ax11.plot(temps / 3600, conversion.CtoF(Res2results.T[paramRes2.nb_y - 1, :]),
              label='$\ity$={0:.2f} ft'.format(conversion.mtoft(paramRes2.y_num[0])))  # Bas du réservoir
    if paramHXefd.m_in>0:
        ax11.plot(temps / 3600, conversion.CtoF(HXefdresults.Tout[0, :]), label='$\itHX_{out}$')  # sortie DHW HX
    ax11.set_xlabel('Time [h]', color='k', fontname='Lucida Bright', fontsize=9)
    ax11.set_ylabel('$\itT_{HW}$ or $\itT_{HX_{out}}$ [°F]', color='k', fontname='Lucida Bright', fontsize=9)
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=3, prop={'family': 'Lucida Bright', 'size': 9})
    ax11.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20), ncol=2, prop={'family': 'Lucida Bright', 'size': 9})
    fig1.set_figheight(fig1.get_figheight() * 1.5)
    fig1.set_figwidth(fig1.get_figwidth())
    fig1.tight_layout()  # otherwise the right y-label is slightly clipped
    #fig1.savefig('Tank Simulation results.jpg', format='JPG', dpi=300)
    st.pyplot(fig1)

# Analyse des Échelles de temps
with st.expander("Time Scale Analysis", expanded=True):
    st.markdown('### Volume')
    st.markdown('Tank internal volume: **{0:.0f} [Gal]**'.format(conversion.m3toGal(paramRes2.Vres)))
    st.markdown("Heat exchange area of Coil heat Exchanger **{0:.2f} [ft^2]**".format(conversion.m2toft2(np.pi * paramHXefd.Do * paramHXefd.Ltot)))
    st.markdown("Volume of the Coil Heat Exchanger in the tank: **{0:.0f} [Gal]**".format(conversion.m3toGal(paramHXefd.Aext * paramHXefd.Ltot)))
    Volume_eau_EC = paramRes2.Vres - paramHXefd.Aext * paramHXefd.Ltot
    st.markdown("Volume of Hot Water in the tank (HW): **{0:.0f} [Gal]**".format(conversion.m3toGal(Volume_eau_EC)))
    Volume_eau_EFD = paramHXefd.As * paramHXefd.Ltot
    st.markdown("Volume of Water in the Coil Heat Exchanger: **{0:.0f} [Gal]**".format(conversion.m3toGal(Volume_eau_EFD)))
    Volume_cuivre = (paramHXefd.Aext - paramHXefd.As) * paramHXefd.Ltot
    st.markdown("Volume of Copper in the Coil heat Exchanger **{0:.0f} [Gal]**".format(conversion.m3toGal(Volume_cuivre)))
    st.markdown('### Thermal capacity')
    st.markdown("Total Thermal capacity of the water in the tank: **{0:.2f} [kJ/K]**".format(data.rho_EC * Volume_eau_EC * data.Cp_EC * 1e-3))
    st.markdown("Total Thermal capacity of the water in the coil heat exchanger @300K: **{0:.2f} [kJ/K]**".format(data.rho_EC * Volume_eau_EFD * data.Cp_EC * 1e-3))
    st.markdown("Total thermal capacity of Copper: **{0:.2f} [kJ/K]**".format(
        data.rho_cu * Volume_cuivre * data.cp_cu * 1e-3))
    st.markdown('### Transport time')
    st.markdown('Transport time of Hot Water (HW): **{0:.2f} [min]**'.format(data.rho_EC * Volume_eau_EC / paramRes2.m_in / 60))
    st.markdown('Transport time of water in the coil heat exchanger: **{0:.2f} [s]**'.format(data.rho_EC * Volume_eau_EFD / paramHXefd.m_in))

    st.markdown('Tank characteristic time: **{0:.2f} [h]**'.format((paramRes2.Di / 2) ** 2 / (4 * data.alpha_EC) / 3600))

print("**** Calculation Successfull ****")