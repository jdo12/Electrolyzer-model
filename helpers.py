# -*- coding: utf-8 -*-
# !/usr/bin/env python3

from math import log, log10, exp
from numpy import polyfit, RankWarning
from pandas import read_csv
from warnings import simplefilter
simplefilter('ignore', RankWarning)


def fit_pdat(filename, time_interval=None):
    """
    Read demand profile and check if timestamp interval is at least in seconds.
    Execute only if freq does not match: milli, micro, or nano
    :param filename: Time-stamped power curve input data as CSV
    :param time_interval: Specify the data time interval in seconds
    :return p: Curve-fit results using Least square curve-fitting method
    """
    if time_interval is None:
        try:
            df = read_csv(filename, header=None, parse_dates=[0], index_col=0, nrows=10)
            interval = df.index.inferred_freq
        except (AttributeError, TypeError):
            raise Exception("\tUnable to infer time interval from input data\n"
                            "\tPlease specify the time_interval parameter in seconds!!\n")
        if any(i in interval for i in ('L', 'ms', 'U', 'us', 'N')) is False:
            print("Demand data time-interval is at least in seconds:\t\tOKAY")
            df = read_csv(filename, header=None, parse_dates=[0], names=['datetime', 'Pval'])
            df['sec'] = df['datetime'].sub(df['datetime'].iloc[0]).dt.total_seconds()
        else:
            raise Exception("\tTime-series data should be at least in seconds\n"
                            "\te.g.:\n"
                            "\t\t2019-01-21 00:00:00  646.969697\n"
                            "\t\t2019-01-21 00:05:00  642.818182\n"
                            "\t\t...................  ..........\n"
                            "\t\tYYYY-MM-DD HH:MM:SS  Value\n")
    else:
        df = read_csv(filename, header=None, names=['Pval'])
        df['sec'] = [i for i in range(0, len(df)*time_interval, time_interval)]

    # Least-square curve-fitting
    return polyfit(df['sec'].values, df['Pval'].values, 21), int(df['sec'].iloc[-1])


def vi_calc(x, Pri, Nc, V_rev, T, r1, r2, s1, s2, s3, t1, t2, t3, A):
    fa = list([x[0] - V_rev - (r1+r2*T)*x[1]/A - (s1+s2*T+s3*T**2)*log10((t1+t2/T+t3/T**2)*x[1]/A+1)])
    fa.append((Nc*x[0])*x[1]-Pri)
    return fa


def vifc_calc(x, Pri, V_rev, T, R, alpha, F, i_n, i_0, R_tot, a, b, pH2, pO2):
    #print(R, T, alpha, F, x[1], i_n, i_0)
    Vact = (R*T/alpha/F)*log((x[1]+i_n)/i_0)
    #print(Vact)
    Vohm = (x[1] + i_n)*R_tot
    Vcon = a*exp(b*x[1])
    E_cell = V_rev+R*T/2/F*log(pH2*pO2**0.5)
    fb = list([x[0]-E_cell+Vact+Vohm-Vcon])
    fb.append(x[0]*x[1]-Pri)
    return fb
