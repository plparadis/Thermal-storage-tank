# -*- coding: utf-8 -*-
'''
# Unit Conversion Module 
# Developed by PLParadis

'''

def FtoC(TF):
    # Convertie les Fahrenheit en Celsius
    # TC: Temperature [°C]
    # TF: Temperature [°F]
    # note: 0°C = 32°F et 100°C = 212°F
    # Conversion: °F = (9/5)*°C + 32
    TC = (TF - 32)*(5/9);
    return TC
def CtoF(TC):
    # Convertie les Fahrenheit en Celsius
    # TC: Temperature [°C]
    # TF: Temperature [°F]
    # note: 0°C = 32°F et 100°C = 212°F
    # Conversion: °F = (9/5)*°C + 32
    TF = TC*(9/5)+32;
    return TF
def RtoK(DTr):
    # Convertie une différence de température en Rankine vers Kelvin
    # DTk: Différence de Temperature [K]
    # DTr: Différence de Temperature [R]
    # note: 0°C = 32°F et 100°C = 212°F
    # Conversion: °F = (9/5)*°C + 32
    # T1[°F] - T2[°F] = (9/5)*T1[°C] + 32 - ((9/5)*T2[°C] + 32)
    # (T1 - T2)[R] = (9/5)*(T1 - T2)[K]
    DTk = (5/9)*DTr;
    return DTk
def KtoR(DTk):
    # Convertie une différence de température en Rankine vers Kelvin
    # DTk: Différence de Temperature [K]
    # DTr: Différence de Temperature [R]
    # note: 0°C = 32°F et 100°C = 212°F
    # Conversion: °F = (9/5)*°C + 32
    # T1[°F] - T2[°F] = (9/5)*T1[°C] + 32 - ((9/5)*T2[°C] + 32)
    # (T1 - T2)[R] = (9/5)*(T1 - T2)[K]
    DTr = (9/5)*DTk;
    return DTr

def PSItoKPA(Ppsig):
    # Convertie une Pression en psig vers kPa
    # Pkpa: Pression [kPa]
    # Ppsig: Pression [psig]
    # note: 0.14503773800722 psig = 1 kPa
    Pkpa = Ppsig/0.14503773800722;
    return Pkpa

def KPAtoPSI(Pkpa):
    # Convertie une Pression en kpa vers psi
    # Pkpa: Pression [kPa]
    # Ppsig: Pression [psig]
    # note: 0.14503773800722 psig = 1 kPa
    Ppsig = Pkpa*0.14503773800722;
    return Ppsig

def LBHtoKGS(mlbh):
    # Convertie un dÃ©bit massique en lb/h vers kg/s
    # mlbh: Débit massique [lb/h]
    # mkgs: Débit massique [kg/s]
    # note: 1 kg/s = 7936.641438 lb/h
    mkgs = mlbh/7936.641438;
    return mkgs

def LBtoKG(mlb):
    # Convertie une masse en lb vers kg
    # mlb: masse [lb]
    # mkg: masse [kg]
    # note: 1 kg = 2.20462 lb
    mkg = mlb/2.20462;
    return mkg

def KGtoLB(mkg):
    # Convertie une masse en kg vers lb
    # mlb: masse [lb]
    # mkg: masse [kg]
    # note: 1 kg = 2.20462 lb
    mlb = mkg*2.20462;
    return mlb

def KJKGtoBTULB(hkjkg):
    # Convertie l'enthalpie en kJ/kg vers btu/lb
    # hkjkg: Enthalpie [kJ/kg]
    # hbtulb: Enthalpie [btu/lb]
    # Conversion: 1 kJ/kg = 0.429923 Btu/lb
    hbtulb = hkjkg*0.429923;
    return hbtulb

def BTULBtoKJKG(hbtulb):
    # Convertie l'enthalpie en btu/lb vers kJ/kg
    # hkjkg: Enthalpie [kJ/kg]
    # hbtulb: Enthalpie [btu/lb]
    # Conversion: 1 kJ/kg = 0.429923 Btu/lb
    hkjkg = hbtulb/0.429923;
    return hkjkg

def WtoBTUhr(qW):
    # Convertie la puissance en W vers btu/hr
    # qW: Puissance [W]
    # qBTUhr: Puissance [btu/hr]
    # Conversion: 1 W = 3.412142 BTU/hr
    qBTUhr = qW*3.412142;
    return qBTUhr

def BTUhrtoW(qBTUhr):
    # Convertie la puissance en W vers btu/hr
    # qW: Puissance [W]
    # qBTUhr: Puissance [btu/hr]
    # Conversion: 1 W = 3.412142 BTU/hr
    qW = qBTUhr/3.412142;
    return qW

def BTUtoJ(eBTU):
    # Convertie l'energie en Btu vers Joules
    # eJ: Energie [J]
    # eBTU: Energie [btu]
    # Conversion: 1 Btu = 1055 Joules
    eJ = eBTU*1055
    return eJ

def JtoBTU(eJ):
    # Convertie l'energie en Joules vers Btu
    # eJ: Energie [J]
    # eBTU: Energie [btu]
    # Conversion: 1 Btu = 1055 Joules
    eBTU = eJ/1055
    return eBTU

def CFMtom3sec(VCFM):
    # Convertie le debit volumique en CFM vers m3/sec
    # VCFM: Debit volumique [CFM]
    # Vm3sec: Debit volumique [m3/sec]
    # Conversion: 2118.8799727597 CFM = 1 m3/sec
    Vm3sec = VCFM/2118.8799727597
    return Vm3sec


def m3sectoCFM(Vm3sec):
    # Convertie le debit volumique en m3/sec vers CFM
    # Vm3sec: Debit volumique [m3/sec]
    # VCFM: Debit volumique [CFM]
    # Conversion: 1 m3/sec = 2 118.88 [CFM]
    VCFM = Vm3sec*2118.8799727597
    return VCFM

def GPMtoLPS(Vgpm):
    # Convertie le debit volumique en gpm vers l/sec
    # Vgpm: Debit volumique [gpm]
    # Vlps: Debit volumique [l/sec]
    # Conversion: 3.7854118 l = 1 gallon
    Vlps = Vgpm*3.7854118/60
    return Vlps

def LPStoGPM(Vlps):
    # Convertie le debit volumique en gpm vers l/sec
    # Vgpm: Debit volumique [gpm]
    # Vlps: Debit volumique [l/sec]
    # Conversion: 3.7854118 l = 1 gallon
    Vgpm = Vlps/3.7854118*60
    return Vgpm

def ft2tom2(ft2):
    # Convertie les square feets en square meters
    # ft2: area [ft2]
    # m2: area [m2]
    # note: 12 [in] = 1 [ft] and 1 [in] = 25.4 [mm] and 1000 [mm] = 1 [m]
    # Conversion: 
    m2 = ft2*(12*25.4/1000)**2
    return m2

def m2toft2(m2):
    # Convertie les square meters en square feets
    # ft2: area [ft2]
    # m2: area [m2]
    # note: 12 [in] = 1 [ft] and 1 [in] = 25.4 [mm] and 1000 [mm] = 1 [m]
    # Conversion: 
    ft2 = m2*(12*25.4/1000)**-2
    return ft2

def mtoft(m):
    # Convertie les meters en feets
    # ft: length [ft]
    # m: length [m]
    # note: 12 [in] = 1 [ft] and 1 [in] = 25.4 [mm] and 1000 [mm] = 1 [m]
    # Conversion: 
    ft = m*(12*25.4/1000)**-1
    return ft

def fttom(ft):
    # Convertie les feets en meters
    # ft: length [ft]
    # m: length [m]
    # note: 12 [in] = 1 [ft] and 1 [in] = 25.4 [mm] and 1000 [mm] = 1 [m]
    # Conversion: 
    m = ft*(12*25.4/1000)
    return m

def m3kgtoft3lb(m3kg):
    # Convertie les m3/kg en ft3/lb
    # ft3lb: density [ft]
    # m3kg: density [m]
    # note: 1m3/kg = 16.0185ft3/lb
    # Conversion: 
    ft3lb = m3kg*(16.0185)
    return ft3lb


def m3toGal(m3):
    # Convertie les m3 en Gal
    # Gal: Volume [Gal]
    # m3: Volume [m3]
    # note: 1m3 = 1000liters and 3.78541liters = 1Gal
    # Conversion:
    Gal = m3*1000/3.78541
    return Gal









