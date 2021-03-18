# -*- coding: utf-8 -*-
"""
Unit Conversion Module
--> Developed by PLParadis
"""


def FtoC(TF):
    """
    Convertie les Fahrenheit en Celsius
        note: 0°C = 32°F et 100°C = 212°F
        Conversion: °F = (9/5)*°C + 32
    :param TF: Temperature [°F]
    :return TC: Temperature [°C]
    """
    TC = (TF - 32) * (5 / 9)
    return TC


def CtoF(TC):
    """
    Convertie les Fahrenheit en Celsius
        note: 0°C = 32°F et 100°C = 212°F
        Conversion: °F = (9/5)*°C + 32
    :param TC: Temperature [°C]
    :return TF: Temperature [°F]
    """
    TF = TC * (9 / 5) + 32
    return TF


def RtoK(DTr):
    """
    Convertie une différence de température en Rankine vers Kelvin
        note: 0°C = 32°F et 100°C = 212°F
        Conversion: °F = (9/5)*°C + 32
        T1[°F] - T2[°F] = (9/5)*T1[°C] + 32 - ((9/5)*T2[°C] + 32)
        (T1 - T2)[R] = (9/5)*(T1 - T2)[K]
    :param DTr: Différence de Temperature [R]
    :return DTk: Différence de Temperature [K]
    """
    DTk = (5 / 9) * DTr
    return DTk


def KtoR(DTk):
    """
    Convertie une différence de température en Rankine vers Kelvin
        note: 0°C = 32°F et 100°C = 212°F
        Conversion: °F = (9/5)*°C + 32
        T1[°F] - T2[°F] = (9/5)*T1[°C] + 32 - ((9/5)*T2[°C] + 32)
        (T1 - T2)[R] = (9/5)*(T1 - T2)[K]
    :param DTk: Différence de Temperature [K]
    :return DTr: Différence de Temperature [R]
    """
    DTr = (9 / 5) * DTk
    return DTr


def PSItoKPA(Ppsig):
    """
    Convertie une Pression en psig vers kPa
        note: 0.14503773800722 psig = 1 kPa
    :param Ppsig: Pression [psig]
    :return Pkpa: Pression [kPa]
    """
    Pkpa = Ppsig / 0.14503773800722
    return Pkpa


def KPAtoPSI(Pkpa):
    """
    Convertie une Pression en kpa vers psi
        note: 0.14503773800722 psig = 1 kPa
    :param Pkpa: Pression [kPa]
    :return Ppsig: Pression [psig]
    """
    Ppsig = Pkpa * 0.14503773800722
    return Ppsig


def LBHtoKGS(mlbh):
    """
    Convertie un débit massique en lb/h vers kg/s
        note: 1 kg/s = 7936.641438 lb/h
    :param mlbh: Débit massique [lb/h]
    :return mkgs: Débit massique [kg/s]
    """
    mkgs = mlbh / 7936.641438
    return mkgs


def LBtoKG(mlb):
    """
    Convertie une masse en lb vers kg
        note: 1 kg = 2.20462 lb
    :param mlb: masse [lb]
    :return mkg: masse [kg]
    """
    mkg = mlb / 2.20462
    return mkg


def KGtoLB(mkg):
    """
    Convertie une masse en kg vers lb
        note: 1 kg = 2.20462 lb
    :param mkg: masse [kg]
    :return mlb: masse [lb]
    """
    mlb = mkg * 2.20462
    return mlb


def KJKGtoBTULB(hkjkg):
    """
    Convertie l'enthalpie en kJ/kg vers btu/lb
        Conversion: 1 kJ/kg = 0.429923 Btu/lb
    :param hkjkg: Enthalpie [kJ/kg]
    :return hbtulb: Enthalpie [btu/lb]
    """
    hbtulb = hkjkg * 0.429923
    return hbtulb


def BTULBtoKJKG(hbtulb):
    """
    Convertie l'enthalpie en btu/lb vers kJ/kg
        Conversion: 1 kJ/kg = 0.429923 Btu/lb
    :param hbtulb: Enthalpie [btu/lb]
    :return hkjkg: Enthalpie [kJ/kg]
    """
    hkjkg = hbtulb / 0.429923
    return hkjkg


def WtoBTUhr(qW):
    """
    Convertie la puissance en W vers btu/hr
        Conversion: 1 W = 3.412142 BTU/hr
    :param qW: Puissance [W]
    :return qBTUhr: Puissance [btu/hr]
    """
    qBTUhr = qW * 3.412142
    return qBTUhr


def BTUhrtoW(qBTUhr):
    """
    Convertie la puissance en W vers btu/hr
        Conversion: 1 W = 3.412142 BTU/hr
    :param qBTUhr: Puissance [btu/hr]
    :return qW: Puissance [W]
    """
    qW = qBTUhr / 3.412142
    return qW


def BTUtoJ(eBTU):
    """
    Convertie l'energie en Btu vers Joules
        Conversion: 1 Btu = 1055 Joules
    :param eBTU: Energie [btu]
    :return eJ: Energie [J]
    """
    eJ = eBTU * 1055
    return eJ


def JtoBTU(eJ):
    """
    Convertie l'energie en Joules vers Btu
        Conversion: 1 Btu = 1055 Joules
    :param eJ: Energie [J]
    :return eBTU: Energie [btu]
    """
    eBTU = eJ / 1055
    return eBTU


def CFMtom3sec(VCFM):
    """
    Convertie le debit volumique en CFM vers m3/sec
        Conversion: 2118.8799727597 CFM = 1 m3/sec
    :param VCFM: Debit volumique [CFM]
    :return Vm3sec: Debit volumique [m3/sec]
    """
    Vm3sec = VCFM / 2118.8799727597
    return Vm3sec


def m3sectoCFM(Vm3sec):
    """
    Convertie le debit volumique en m3/sec vers CFM
        Conversion: 1 m3/sec = 2 118.88 [CFM]
    :param Vm3sec: Debit volumique [m3/sec]
    :return VCFM: Debit volumique [CFM]
    """
    VCFM = Vm3sec * 2118.8799727597
    return VCFM


def GPMtoLPS(Vgpm):
    """
    Convertie le debit volumique en gpm vers l/sec
        Conversion: 3.7854118 l = 1 gallon
    :param Vgpm: Debit volumique [gpm]
    :return Vlps: Debit volumique [l/sec]
    """
    Vlps = Vgpm * 3.7854118 / 60
    return Vlps


def LPStoGPM(Vlps):
    """
    Convertie le debit volumique en gpm vers l/sec
        Conversion: 3.7854118 l = 1 gallon
    :param Vlps: Debit volumique [l/sec]
    :return Vgpm: Debit volumique [gpm]
    """
    Vgpm = Vlps / 3.7854118 * 60
    return Vgpm


def ft2tom2(ft2):
    """
    Convertie les square feets en square meters
    note: 12 [in] = 1 [ft] and 1 [in] = 25.4 [mm] and 1000 [mm] = 1 [m]
    :param ft2: area [ft2]
    :return m2: area [m2]
    """
    m2 = ft2 * (12 * 25.4 / 1000) ** 2
    return m2


def m2toft2(m2):
    """
    Convertie les square meters en square feets
    note: 12 [in] = 1 [ft] and 1 [in] = 25.4 [mm] and 1000 [mm] = 1 [m]
    :param m2: area [m2]
    :return ft2: area [ft2]
    """
    ft2 = m2 * (12 * 25.4 / 1000) ** -2
    return ft2


def mtoft(m):
    """
    Convertie les meters en feets
    note: 12 [in] = 1 [ft] and 1 [in] = 25.4 [mm] and 1000 [mm] = 1 [m]
    :param m: length [m]
    :return ft: length [ft]
    """
    ft = m * (12 * 25.4 / 1000) ** -1
    return ft


def fttom(ft):
    """
    Convertie les feets en meters
    note: 12 [in] = 1 [ft] and 1 [in] = 25.4 [mm] and 1000 [mm] = 1 [m]
    :param ft: length [ft]
    :return m: length [m]
    """
    m = ft * (12 * 25.4 / 1000)
    return m


def m3kgtoft3lb(m3kg):
    """
    Convertie les m3/kg en ft3/lb
    note: 1m3/kg = 16.0185ft3/lb
    :param m3kg: density [m]
    :return ft3lb: density [ft]
    """
    ft3lb = m3kg * 16.0185
    return ft3lb


def m3toGal(m3):
    """
    Convertie les m3 en Gal
    note: 1m3 = 1000liters and 3.78541liters = 1Gal
    :param m3: Volume [m3]
    :return Gal: Volume [Gal]
    """
    Gal = m3 * 1000 / 3.78541
    return Gal
