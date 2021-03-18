import numpy as np
import sys
from scipy.special import erf


def ierf(x):
    y = x * erf(x) - (1 - np.exp(-x * x)) / np.sqrt(np.pi)
    return y


def f_Eim(xx):
    x = -xx
    if x < -5.0:
        y = Continued_Fraction_Ei(x)
    elif x == 0.0:
        y = realmax
    elif x < 6.8:
        y = Power_Series_Ei(x)
    elif x < 50.0:
        y = Argument_Addition_Series_Ei(x)
    else:
        y = Continued_Fraction_Ei(x)
    return -y


# ////////////////////////////////////////////////////////////////////////////////
# //  Continued_Fraction_Ei( double x )                  //
# //                                                                            //
# //  Description:                                                              //
# //     For x < -5 or x > 50, the continued fraction representation of Ei      //
# //     converges fairly rapidly.                                              //
# //                                                                            //
# //     The continued fraction expansion of Ei(x) is:                          //
# //        Ei(x) = -exp(x) { 1/(-x+1-) 1/(-x+3-) 4/(-x+5-) 9/(-x+7-) \ }.    //
# //                                                                            //
# //                                                                            //
# //  Arguments:                                                                //
# //     double  x                                                         //
# //                The argument of the exponential integral Ei().              //
# //                                                                            //
# //  Return Value:                                                             //
# //     The value of the exponential integral Ei evaluated at x.               //
# ////////////////////////////////////////////////////////////////////////////////

realmax = sys.float_info.max
eps = np.spacing(1)


def Continued_Fraction_Ei(x):
    Am1 = 1.0
    A0 = 0.0
    Bm1 = 0.0
    B0 = 1.0
    a = np.exp(x)
    b = -x + 1.0
    Ap1 = b * A0 + a * Am1
    Bp1 = b * B0 + a * Bm1
    a = 1.0
    j = 1
    while (abs(Ap1 * B0 - A0 * Bp1) > eps * abs(A0 * Bp1)):
        if abs(Bp1) > 1.0:
            Am1 = A0 / Bp1
            A0 = Ap1 / Bp1
            Bm1 = B0 / Bp1
            B0 = 1.0
        else:
            Am1 = A0
            A0 = Ap1
            Bm1 = B0
            B0 = Bp1
        a = -j * j
        b = b + 2.0
        Ap1 = b * A0 + a * Am1
        Bp1 = b * B0 + a * Bm1
        j = j + 1
    y = (-Ap1 / Bp1)
    return y


# ////////////////////////////////////////////////////////////////////////////////
# // static double Power_Series_Ei( double x )                        //
# //                                                                            //
# //  Description:                                                              //
# //     For -5 < x < 6.8, the power series representation for                  //
# //     (Ei(x) - gamma - ln|x|)/exp(x) is used, where gamma is Euler's gamma   //
# //     constant.                                                              //
# //     Note that for x = 0.0, Ei is -inf.  In which case -DBL_MAX is          //
# //     returned.                                                              //
# //                                                                            //
# //     The power series expansion of (Ei(x) - gamma - ln|x|) / exp(x) is      //
# //        - Sum(1 + 1/2 + \ + 1/j) (-x)^j / j!, where the Sum extends       //
# //        from j = 1 to inf.                                                  //
# //                                                                            //
# //  Arguments:                                                                //
# //     double  x                                                         //
# //                The argument of the exponential integral Ei().              //
# //                                                                            //
# //  Return Value:                                                             //
# //     The value of the exponential integral Ei evaluated at x.               //
# ////////////////////////////////////////////////////////////////////////////////
#
def Power_Series_Ei(x):
    xn = -x
    Sn = -x
    Sm1 = 0.0
    hsum = 1.0
    g = 0.5772156649015328606065121
    y = 1.0
    factorial = 1.0
    if x == 0.0:
        y = realmax
    else:
        while (abs(Sn - Sm1) > eps * abs(Sm1)):
            Sm1 = Sn
            y = y + 1.0
            xn = xn * (-x)
            factorial = factorial * y
            hsum = hsum + (1.0 / y)
            Sn = Sn + hsum * xn / factorial
        y = (g + np.log(abs(x)) - np.exp(x) * Sn)
    return y


# ////////////////////////////////////////////////////////////////////////////////
# // static double Argument_Addition_Series_Ei(double x)              //
# //                                                                            //
# //  Description:                                                              //
# //     For 6.8 < x < 50.0, the argument addition series is used to calculate  //
# //     Ei.                                                                    //
# //                                                                            //
# //     The argument addition series for Ei(x) is:                             //
# //     Ei(x+dx) = Ei(x) + exp(x) Sum j! [exp(j) expj(-dx) - 1] / x^(j+1),     //
# //     where the Sum extends from j = 0 to inf, |x| > |dx| and expj(y) is     //
# //     the exponential polynomial expj(y) = Sum y^k / k!, the Sum extending   //
# //     from k = 0 to k = j.                                                   //
# //                                                                            //
# //  Arguments:                                                                //
# //     double  x                                                         //
# //                The argument of the exponential integral Ei().              //
# //                                                                            //
# //  Return Value:                                                             //
# //     The value of the exponential integral Ei evaluated at x.               //
# ////////////////////////////////////////////////////////////////////////////////
def Argument_Addition_Series_Ei(x):
    ei = np.array([
        1.915047433355013959531e2, 4.403798995348382689974e2,
        1.037878290717089587658e3, 2.492228976241877759138e3,
        6.071406374098611507965e3, 1.495953266639752885229e4,
        3.719768849068903560439e4, 9.319251363396537129882e4,
        2.349558524907683035782e5, 5.955609986708370018502e5,
        1.516637894042516884433e6, 3.877904330597443502996e6,
        9.950907251046844760026e6, 2.561565266405658882048e7,
        6.612718635548492136250e7, 1.711446713003636684975e8,
        4.439663698302712208698e8, 1.154115391849182948287e9,
        3.005950906525548689841e9, 7.842940991898186370453e9,
        2.049649711988081236484e10, 5.364511859231469415605e10,
        1.405991957584069047340e11, 3.689732094072741970640e11,
        9.694555759683939661662e11, 2.550043566357786926147e12,
        6.714640184076497558707e12, 1.769803724411626854310e13,
        4.669055014466159544500e13, 1.232852079912097685431e14,
        3.257988998672263996790e14, 8.616388199965786544948e14,
        2.280446200301902595341e15, 6.039718263611241578359e15,
        1.600664914324504111070e16, 4.244796092136850759368e16,
        1.126348290166966760275e17, 2.990444718632336675058e17,
        7.943916035704453771510e17, 2.111342388647824195000e18,
        5.614329680810343111535e18, 1.493630213112993142255e19,
        3.975442747903744836007e19, 1.058563689713169096306e20])
    k = round(x + 0.5)
    j = 0
    xx = k
    dx = x - xx
    xxj = xx
    edx = np.exp(dx)
    Sm = 1.0
    Sn = (edx - 1.0) / xxj
    term = realmax
    factorial = 1.0
    dxj = 1.0
    while (abs(term) > eps * np.fabs(Sn)):
        j = j + 1
        factorial = factorial * j
        xxj = xxj * xx
        dxj = dxj * (-dx)
        Sm = Sm + (dxj / factorial)
        term = (factorial * (edx * Sm - 1.0)) / xxj
        Sn = Sn + term
    y = ei[k - 7] + Sn * np.exp(xx)
    return y


def cherche_index(xi, x):
    """ cherche l'index où x(i) <= xi < x(i+1)"""
    err = 0
    ok = 1
    i = 0
    if (xi > 0 and xi < x[0]):
        i = 0
    elif xi == x[len(x) - 1]:
        i = len(x) - 1
    else:
        i = i + 1
        while ok:
            if (xi >= x[i - 1]) & (xi < x[i]):
                ok = 0
            else:
                i = i + 1
                if i >= len(x):
                    ok = 0
                    err = -1
                    i = np.nan
    return i


def cherche_index2(xi, x):
    """ cherche l'index où x(i) <= xi < x(i+1)"""
    err = 0
    ok = 1
    i = 0
    if (xi > 0 and xi <= x[0]):
        i = 0
    elif xi == x[len(x) - 1]:
        i = len(x) - 1
    else:
        while ok:
            if (xi >= x[i]) & (xi < x[i + 1]):
                ok = 0
            else:
                i = i + 1
                if i >= len(x):
                    ok = 0
                    err = -1
                    i = np.nan
    return i


def mon_interp2(x, y, z, x1, y1):
    i = cherche_index2(x1, x)
    j = cherche_index2(y1, y)
    A = np.array([[x[i], y[j], x[i] * y[j], 1],
                  [x[i + 1], y[j], x[i + 1] * y[j], 1],
                  [x[i], y[j + 1], x[i] * y[j + 1], 1],
                  [x[i + 1], y[j + 1], x[i + 1] * y[j + 1], 1]])
    B = np.array([z[i, j], z[i + 1, j], z[i, j + 1], z[i + 1, j + 1]])
    d = np.linalg.solve(A, B)
    z1 = d[0] * x1 + d[1] * y1 + d[2] * x1 * y1 + d[3]
    return z1


def combinaisons(n, k):
    if k == 0:
        C = 1
    elif k == 1:
        C = n
    else:
        p = n - k + 1
        for i in range(2, k + 1):
            p = p * (n - k + i) / (i * 1.0)
        C = p
    return C


def tic():
    # Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()


def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print('Elapsed time is' + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set")
