import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
import time
from typing import Any
import random

@dataclass
class Enzyme:
    name: str
    c: float
    kinetic_factors: dict = field(default_factory=dict)

    def __post_init__(self):
        self.kinetic_factors = pd.Series(self.kinetic_factors)


class Simulation:
    def __init__(self):
        #Enzymes
        self.ras = Enzyme('RAS', 0, {'v_1': 2.5,
                                     'k_mm_1': 10.})
        self.mkkk = Enzyme('MKKK', 100.)
        self.mkp_mkkk = Enzyme('MKP_MKKK', 0, {'v_2': 0.25,
                                               'k_mm_2': 8.})
        self.mkkk_p = Enzyme('MKKK_P', 0., {'k_cat_3': 0.025,
                                            'k_mm_3': 15.,
                                            'k_cat_4': 0.025,
                                            'k_mm_4': 15.})
        self.mkp_mkk = Enzyme('MKP_MKK', 0, {'v_5': 0.75,
                                             'k_mm_5': 15.,
                                             'v_6': 0.75,
                                             'k_mm_6': 15.})
        self.mkk = Enzyme('MKK', 300.)
        self.mkk_p = Enzyme('MKK_P', 0.)
        self.mkk_pp = Enzyme('MKK_PP', 0., {'k_cat7': 0.025,
                                            'k_mm7': 15.,
                                            'k_cat8': 0.025,
                                            'k_mm8': 15.})
        self.mkp_mapk = Enzyme('MKP_MAPK', 0, {'v_9': 0.5,
                                               'k_mm_9': 15.,
                                               'v_10': 0.5,
                                               'k_mm_10': 15.})
        self.mapk = Enzyme('MAPK', 300.)
        self.mapk_p = Enzyme('MAPK_P', 0)
        self.mapk_pp = Enzyme('MAPK_PP', 0, {'n': 0.,
                                             'ki': 9.})
        self.mapk_list = [self.mkkk,self.mkkk_p, self.mkk, self.mkk_p,
                          self.mkk_pp, self.mapk, self.mapk_p, self.mapk_pp]
        # constants timeframe
        self.t0 = 0
        self.t_end_minutes = 140
        self.t_end = self.t_end_minutes * 60
        self.stepps = 400
        self.t_list = np.linspace(self.t0, self.t_end, self.stepps)
        # Todo: FFunktion timeframe
        self.solution = pd.DataFrame()


    def simulate_ode_ruck(self, t, dc_dt_list):
        self.mkkk.c, self.mkkk_p.c, self.mkk.c, self.mkk_p.c, self.mkk_pp.c,\
            self.mapk.c, self.mapk_p.c, self.mapk_pp.c = dc_dt_list
        ki = self.mapk_pp.kinetic_factors['ki']
        n = self.mapk_pp.kinetic_factors['n']
        v1 = self.ras.kinetic_factors['v_1']
        k_mm1 = self.ras.kinetic_factors['k_mm_1']
        v2, k_mm2 = self.mkp_mkkk.kinetic_factors
        k_cat3, k_mm3, k_cat4, k_mm4 = self.mkkk_p.kinetic_factors
        v5, k_mm5, v6, k_mm6 = self.mkp_mkk.kinetic_factors
        k_cat7, k_mm7, k_cat8, k_mm8 = self.mkk_pp.kinetic_factors
        v9, k_mm9, v10, k_mm10 = self.mkp_mapk.kinetic_factors

        # rates
        rate1 = v1 * self.mkkk.c / ((1 + (self.mapk_pp.c / ki) ** n) * (k_mm1 + self.mkkk.c))
        rate2 = v2 * self.mkkk_p.c / (k_mm2 + self.mkkk_p.c)
        rate3 = k_cat3 * self.mkkk_p.c * self.mkk.c / (k_mm3 + self.mkk.c)
        rate4 = k_cat4 * self.mkkk_p.c * self.mkk_p.c / (k_mm4 + self.mkk_p.c)
        rate5 = v5 * self.mkk_pp.c / (k_mm5 + self.mkk_pp.c)
        rate6 = v6 * self.mkk_p.c / (k_mm6 + self.mkk_p.c)
        rate7 = k_cat7 * self.mkk_pp.c * self.mapk.c / (k_mm7 + self.mapk.c)
        rate8 = k_cat8 * self.mkk_pp.c * self.mapk_p.c / (k_mm8 + self.mapk_p.c)
        rate9 = v9 * self.mapk_pp.c / (k_mm9 + self.mapk_pp.c)
        rate10 = v10 * self.mapk_p.c / (k_mm10 + self.mapk_p.c)

        # d_kin/d_t
        d_mkkk = rate2 - rate1
        d_mkkk_p = rate1 - rate2
        d_mkk = rate6 - rate3
        d_mkk_p = rate3 + rate5 - rate4 - rate6
        d_mkk_pp = rate4 - rate5
        d_mapk = rate10 - rate7
        d_mapk_p = rate7 + rate9 - rate8 - rate10
        d_mapk_pp = rate8 - rate9

        return [d_mkkk, d_mkkk_p, d_mkk, d_mkk_p, d_mkk_pp, d_mapk, d_mapk_p, d_mapk_pp]

    def scipy_rk45(self, order=True):
        c_t0 = [self.mkkk.c, self.mkkk_p.c, self.mkk.c, self.mkk_p.c,
                self.mkk_pp.c, self.mapk.c, self.mapk_p.c, self.mapk_pp.c]
        scipy_rk45 = solve_ivp(self.simulate_ode_ruck,
                               [self.t0, self.t_end],
                               c_t0,
                               dense_output=True,
                               rtol=1e-3,
                               atol=1e-6,
                               method='RK45')
        solution_array = scipy_rk45.sol(self.t_list)
        df = pd.DataFrame()
        for i, value in enumerate(solution_array):
            df[self.mapk_list[i].name] = value
        if order:
            self.solution = df
            # return df
        else:
            df = df.sample(frac=1, axis=1)
            self.solution = df
            # return df

    def plot(self, list):
        t_list_min = self.t_list / 60
        if len(list) == 0:
            sns.lineplot(data=self.solution, legend='full')
        else:
            labels = []
            for i in list:
                sns.lineplot(x=t_list_min, y=i, data=self.solution, label=i)
                labels.append(i)
            plt.xlabel('t [min]')
            plt.ylabel('c [nM]')
            # plt.legend(bbox_to_anchor=(1, 1))
            plt.legend()
        plt.show()




sim1 = Simulation()
sim1.mapk_pp.kinetic_factors['n'] = 1.1
sim1.t_end_minutes = 100
print(sim1.t_end_minutes)
sim1.scipy_rk45(order=True)
sim1.plot(['MAPK', 'MAPK_PP'])
#sim1.plot(['MKKK', 'MKKK_P', 'MKK', 'MKK_P', 'MKK_PP', 'MAPK', 'MAPK_P', 'MAPK_PP'])