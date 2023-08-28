import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import random
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
import time
from typing import Any
import random
import os

def random_expdist_boundary(lower_bound, upper_bound, scale_parameter = 1.):
    random_number = np.random.exponential(scale=scale_parameter)
    while random_number < lower_bound or random_number > upper_bound:
        random_number = np.random.exponential(scale=scale_parameter)
    return random_number


@dataclass
class Enzyme:
    name: str
    c: float
    kinetic_factors: dict = field(default_factory=dict)

    def __post_init__(self):
        self.kinetic_factors = pd.Series(self.kinetic_factors, dtype='float64')


class Simulation:
    def __init__(self):
        # Enzymes
        self.signal_lyase = Enzyme('Lyase', 100., {'k_cat': 0.025,
                                                   'k_mm': 10.,
                                                   'i_c': 1.,
                                                   'i_nc': 1.,
                                                   'i_uc': 1.})
        self.signal = Enzyme('Signal', 0, {'k_cat': 0.1,
                                           'k_mm': 10.,
                                           'i_c': 1.,
                                           'i_nc': 1.,
                                           'i_uc': 1.})
        self.signal_active = Enzyme('Signal_A', 200, {'k_cat': 0.1,
                                                      'k_mm': 10.,
                                                      'i_c': 1.,
                                                      'i_nc': 1.,
                                                      'i_uc': 1.})
        self.ras = Enzyme('RAS', 100.)
        self.ras_active = Enzyme('RAS_A', 0, {'k_cat': 0.025,
                                              'k_mm': 10.,
                                              'i_c': 1.,
                                              'i_nc': 1.,
                                              'i_uc': 1.})
        self.mkkk = Enzyme('MKKK', 100.)
        self.mkp_mkkk = Enzyme('MKP_MKKK', 10, {'k_cat': 0.025,
                                                'k_mm': 8.,
                                                'i_c': 1.,
                                                'i_nc': 1.,
                                                'i_uc': 1.})
        self.mkkk_p = Enzyme('MKKK_P', 0., {'k_cat': 0.025,
                                            'k_mm': 15.,
                                            'i_c': 1.,
                                            'i_nc': 1.,
                                            'i_uc': 1.})
        self.mkp_mkk = Enzyme('MKP_MKK', 30, {'k_cat': 0.025,
                                              'k_mm': 15.,
                                              'i_c': 1.,
                                              'i_nc': 1.,
                                              'i_uc': 1.})
        self.mkk = Enzyme('MKK', 300.)
        self.mkk_p = Enzyme('MKK_P', 0.)
        self.mkk_pp = Enzyme('MKK_PP', 0., {'k_cat1': 0.025,
                                            'k_mm1': 15.,
                                            'i_c': 1.,
                                            'i_nc': 1.,
                                            'i_uc': 1.})
        self.mkp_mapk = Enzyme('MKP_MAPK', 20, {'v_9': 0.025,
                                                'k_mm_9': 15.,
                                                'i_c': 1.,
                                                'i_nc': 1.,
                                                'i_uc': 1.})
        self.mapk = Enzyme('MAPK', 300.)
        self.mapk_p = Enzyme('MAPK_P', 0)
        self.mapk_pp = Enzyme('MAPK_PP', 0, {'n': 0.,
                                             'ki': 9.})
        self.mapk_list = [self.mkkk, self.mkkk_p, self.mkk, self.mkk_p,
                          self.mkk_pp, self.mapk, self.mapk_p, self.mapk_pp]
        self.mapk_list_normal = [self.mkkk, self.mkkk_p, self.mkk, self.mkk_p,
                                 self.mkk_pp, self.mapk, self.mapk_p, self.mapk_pp,
                                 self.signal, self.signal_active, self.ras, self.ras_active]
        # constants timeframe
        self.t0 = 0
        self.t_end_minutes = 30
        self.t_end = self.t_end_minutes * 60
        self.stepps = 361
        self.t_list = np.linspace(self.t0, self.t_end, self.stepps)

        self.solution = pd.DataFrame()

    def michaelis_menten_kin_inh(self, enzyme=Enzyme, substrate=Enzyme):
        k_cat, k_mm, i_c, i_nc, i_uc = enzyme.kinetic_factors
        c_enzyme = enzyme.c
        c_substrate = substrate.c
        rate = k_cat * c_enzyme * c_substrate / (c_substrate + (k_mm / i_uc) * i_c) * (1 / i_nc) * (1 / i_uc)
        return rate

    def simulate_ode_normal(self, t, dc_dt_list):
        self.mkkk.c, self.mkkk_p.c, self.mkk.c, self.mkk_p.c, self.mkk_pp.c, \
            self.mapk.c, self.mapk_p.c, self.mapk_pp.c, self.signal.c, self.signal_active.c, \
            self.ras.c, self.ras_active.c = dc_dt_list

        rate_minus2 = self.michaelis_menten_kin_inh(self.signal, self.ras_active)
        rate_minus1 = self.michaelis_menten_kin_inh(self.signal_lyase, self.signal_active)
        rate0 = self.michaelis_menten_kin_inh(self.signal_active, self.ras)
        rate1 = self.michaelis_menten_kin_inh(self.ras_active, self.mkkk)
        rate2 = self.michaelis_menten_kin_inh(self.mkp_mkkk, self.mkkk_p)
        rate3 = self.michaelis_menten_kin_inh(self.mkkk_p, self.mkk)
        rate4 = self.michaelis_menten_kin_inh(self.mkkk_p, self.mkk_p)
        rate5 = self.michaelis_menten_kin_inh(self.mkp_mkk, self.mkk_pp)
        rate6 = self.michaelis_menten_kin_inh(self.mkp_mkk, self.mkk_p)
        rate7 = self.michaelis_menten_kin_inh(self.mkk_pp, self.mapk)
        rate8 = self.michaelis_menten_kin_inh(self.mkk_pp, self.mapk_p)
        rate9 = self.michaelis_menten_kin_inh(self.mkp_mapk, self.mapk_pp)
        rate10 = self.michaelis_menten_kin_inh(self.mkp_mapk, self.mapk_p)

        # d_kin/d_t
        d_signal_activ = - rate_minus1
        d_signal = + rate_minus1
        d_ras = -rate0 + rate_minus2
        d_ras_active = rate0 - rate_minus2
        d_mkkk = rate2 - rate1
        d_mkkk_p = rate1 - rate2
        d_mkk = rate6 - rate3
        d_mkk_p = rate3 + rate5 - rate4 - rate6
        d_mkk_pp = rate4 - rate5
        d_mapk = rate10 - rate7
        d_mapk_p = rate7 + rate9 - rate8 - rate10
        d_mapk_pp = rate8 - rate9

        return [d_mkkk, d_mkkk_p, d_mkk, d_mkk_p, d_mkk_pp, d_mapk, d_mapk_p, d_mapk_pp, d_signal, d_signal_activ,
                d_ras, d_ras_active]

    def scipy_rk45(self, order=True):
        c_t0 = [i.c for i in self.mapk_list_normal]
        scipy_rk45 = solve_ivp(self.simulate_ode_normal,
                               [self.t0, self.t_end],
                               c_t0,
                               dense_output=True,
                               rtol=1e-3,
                               atol=1e-6,
                               method='RK45')
        solution_array = scipy_rk45.sol(self.t_list)
        df = pd.DataFrame()
        for i, value in enumerate(solution_array):
            df[self.mapk_list_normal[i].name] = value
        df['t_sec'] = self.t_list
        if order:
            self.solution = df
        else:
            df = df.sample(frac=1, axis=1)
            self.solution = df

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
            plt.legend(bbox_to_anchor=(1, 0.86), loc='upper right')
        plt.show()

    def noisy(self, data_list, sigma=6):
        noisy_list = []
        for i, value in enumerate(data_list):
            noisy_list.append(random.normalvariate(value, sigma))
        return noisy_list


def gen_data(inhibition_type: str, enzymes: str, plot=False, noisy=False, healthy_variance=0.025,
             n_samples=25, min_inhib = 0, set='nonoise_h25_test'):
    '''
    :param inhibition_type: types = lowc, competitive, noncompetitive, uncompetitive,
                             else generates healthy dataset.
    :param enzymes: signal_lyase, signal, signal_active, ras, ras_active, mkkk, mkp_mkkk, mkkk_p, mkp_mkk, mkk, mkk_p,
                    mkk_pp, mkp_mapk, mapk, mapk_p, mapk_pp, mapk_list.
    :param plot:
    :param healthy_variance:
    :param n_samples:
    :return:
    '''
    sim1 = Simulation()
    if not os.path.exists(f'dataset_mapk{sim1.stepps}_n{n_samples}_inh{min_inhib}_{set}.csv'):
        sim1.mapk.c = 300
        sim1.scipy_rk45(order=True)
        for i in sim1.solution:
            if i == 't_sec':
                continue
            sim1.solution[i] = sim1.noisy(sim1.solution[i])
        # init dataframe
        sim1.solution['id'] = 0
        sim1.solution['status'] = 'g'
        sim1.solution['desc'] = '-'
        sim1.solution['inhibition'] = '-'
        sim1.solution['inh_strength'] = '-'
        df = sim1.solution
        df = df.set_index(['id', 'status', 'desc', 'inhibition', 'inh_strength', 't_sec'])
        df.to_csv(f'dataset_mapk{sim1.stepps}_n{n_samples}_inh{min_inhib}_{set}.csv', mode='a', header=True)

    enzymes_list = enzymes.split(',')
    global id_nr
    # solution_list = []
    for i in range(n_samples):
        if i%10 == 0:
            print(i)
        # id = df['id'].max() + 1
        id_nr += 1
        id = id_nr
        sim = Simulation()
        sim_unchanged = Simulation()
        sim.ras.c = np.random.randint(sim.ras.c * (1 - healthy_variance), sim.ras.c * (1 + healthy_variance))
        sim.mapk.c = np.random.randint(sim.mapk.c * (1 - healthy_variance), sim.mapk.c * (1 + healthy_variance))
        sim.mkkk.c = np.random.randint(sim.mkkk.c * (1 - healthy_variance), sim.mkkk.c * (1 + healthy_variance))
        sim.mkk.c = np.random.randint(sim.mkk.c * (1 - healthy_variance), sim.mkk.c * (1 + healthy_variance))
        sim.mkp_mapk.c = np.random.randint(sim.mkp_mapk.c * (1 - healthy_variance),
                                           sim.mkp_mapk.c * (1 + healthy_variance))
        sim.mkp_mkk.c = np.random.randint(sim.mkp_mkk.c * (1 - healthy_variance),
                                          sim.mkp_mkk.c * (1 + healthy_variance))
        sim.mkp_mkkk.c = np.random.randint(sim.mkp_mkkk.c * (1 - healthy_variance),
                                           sim.mkp_mkkk.c * (1 + healthy_variance))

        healthy_inhib = 1
        match inhibition_type:
            case 'lowc':
                inhibition_strength = np.random.rand()
                if inhibition_strength <= min_inhib:
                    inhibition_strength = np.random.rand()
                parameter_string_list = [f'sim.{i}.c' for i in enzymes_list]
                parameter_string_list_unchanged = [f'sim_unchanged.{enzyme}.c' for enzyme in enzymes_list]
                healthy_values_c = []
                for parameter in parameter_string_list_unchanged:
                    healthy = eval(parameter)
                    healthy_values_c.append(healthy)
                for index, c in enumerate(parameter_string_list):
                    exec(f'{c} = inhibition_strength * {c}')
                if eval(parameter_string_list[0]) < healthy_values_c[0]:
                    status = 'k'
                    desc = enzymes
                    inhibition = inhibition_type
                    strength = inhibition_strength
                else:
                    status = 'g'
                    desc = '-'
                    inhibition = '-'
                    strength = inhibition_strength
            case 'competitive':
                parameter_string_list = [f'sim.{enzyme}.kinetic_factors["i_c"]' for enzyme in enzymes_list]
                parameter_string_list_unchanged = [f'sim_unchanged.{enzyme}.kinetic_factors["i_c"]' for enzyme in enzymes_list]
                healthy_values_i_c = [healthy_inhib]
                # inhibition_strength = np.random.uniform(1, 10)
                inhibition_strength = random_expdist_boundary((1+min_inhib*6), 10, 3)
                for index, i_c in enumerate(parameter_string_list):
                    exec(f'{i_c} = inhibition_strength * {i_c}')
                if eval(parameter_string_list[0]) > healthy_values_i_c[0]:
                    status = 'k'
                    desc = enzymes
                    inhibition = inhibition_type
                    strength = inhibition_strength
                else:
                    status = 'g'
                    desc = '-'
                    inhibition = '-'
                    strength = inhibition_strength
            case 'noncompetitive':
                parameter_string_list = [f'sim.{enzyme}.kinetic_factors["i_nc"]' for enzyme in enzymes_list]
                parameter_string_list_unchanged = [f'sim_unchanged.{enzyme}.kinetic_factors["i_nc"]' for enzyme in enzymes_list]
                healthy_values_i_nc = [healthy_inhib]
                # inhibition_strength = np.random.uniform(1, 5)
                inhibition_strength = random_expdist_boundary(1+min_inhib, 5, 1.5)
                for index, i_nc in enumerate(parameter_string_list):
                    exec(f'{i_nc} = inhibition_strength * {i_nc}')
                if eval(parameter_string_list[0]) > healthy_values_i_nc[0]:
                    status = 'k'
                    desc = enzymes
                    inhibition = inhibition_type
                    strength = inhibition_strength
                else:
                    status = 'g'
                    desc = '-'
                    inhibition = '-'
                    strength = inhibition_strength
            case 'uncompetitive':
                parameter_string_list = [f'sim.{enzyme}.kinetic_factors["i_uc"]' for enzyme in enzymes_list]
                parameter_string_list_unchanged = [f'sim_unchanged.{enzyme}.kinetic_factors["i_uc"]' for enzyme in enzymes_list]
                healthy_values_i_uc = [healthy_inhib]
                # inhibition_strength = np.random.uniform(1, 5)
                inhibition_strength = random_expdist_boundary(1+min_inhib, 5, 1.5)
                for index, i_uc in enumerate(parameter_string_list):
                    exec(f'{i_uc} = inhibition_strength * {i_uc}')
                if eval(parameter_string_list[0]) > healthy_values_i_uc[0]:
                    status = 'k'
                    desc = enzymes
                    inhibition = inhibition_type
                    strength = inhibition_strength
                else:
                    status = 'g'
                    desc = '-'
                    inhibition = '-'
                    strength = inhibition_strength
            case _:
                status = 'g'
                desc = '-'
                inhibition = '-'
                strength = 0

        sim.scipy_rk45(order=True)

        if noisy:
            for col in sim.solution:
                if col == 't_sec':
                    continue
                sim.solution[col] = sim.noisy(sim.solution[col])

        if plot:
            # sim.plot(['MKKK', 'MKKK_P', 'MKK', 'MKK_P', 'MKK_PP', 'MAPK', 'MAPK_P', 'MAPK_PP'])
            # sim.plot(['MAPK', 'MAPK_PP'])
            # sim.plot([])
            # sim.plot(['Signal', 'Signal_A', 'RAS', 'RAS_A'])
            sim.plot(
                ['MKKK', 'MKKK_P', 'MKK', 'MKK_P', 'MKK_PP', 'MAPK', 'MAPK_P', 'MAPK_PP', 'Signal', 'RAS', 'RAS_A'])

        sim.solution['id'] = id
        sim.solution['status'] = status
        sim.solution['desc'] = desc
        sim.solution['inhibition'] = inhibition
        sim.solution['inh_strength'] = strength
        sim.solution = sim.solution.set_index(['id', 'status', 'desc', 'inhibition', 'inh_strength', 't_sec'])

        # df = pd.concat([df, sim.solution])
        # solution_list.append(sim.solution)
        sim.solution.to_csv(f'dataset_mapk{sim1.stepps}_n{n_samples}_inh{min_inhib}_{set}.csv', mode='a', header=False)
        # global gen_data_list
        # gen_data_list.append(sim.solution)
    print(f'fin {inhibition_type}, {enzymes}')
    # global gen_data_list
    # gen_data_list.append(solution_list)
    # return solution_list


if __name__ == '__main__':

    sim1 = Simulation()
    sim1.mkkk_p.kinetic_factors['i_nc'] = 1.5
    sim1.scipy_rk45(order=True)
    # for i in sim1.solution:
    #     if i == 't_sec':
    #         continue
    #     sim1.solution[i] = sim1.noisy(sim1.solution[i])
    # sim1.solution['status'] = 'g'
    # sim1.plot(['MAPK', 'MAPK_PP'])
    # sim1.plot([])
    # sim1.plot(['Signal', 'Signal_A', 'RAS', 'RAS_A'])
    sim1.plot(['MKKK', 'MKKK_P', 'MKK', 'MKK_P', 'MKK_PP', 'MAPK', 'MAPK_P', 'MAPK_PP', 'Signal', 'RAS', 'RAS_A'])

    # init dataframe
    sim1.solution['id'] = 0
    sim1.solution['status'] = 'g'
    sim1.solution['desc'] = '-'
    sim1.solution['inhibition'] = '-'
    sim1.solution['inh_strength'] = '-'
    df = sim1.solution
    df = df.set_index(['id', 'status', 'desc', 'inhibition', 'inh_strength', 't_sec'])
    print(df)
    # df.to_csv(f'dataset_mapk{sim1.stepps}_n_test.csv')

    # gen_data_test
    # df = gen_data(df, '--', 'mkk,mkkk', plot=True)
    # df = gen_data(df, 'lowc', 'mkk,mkkk', plot=True)
    # df = gen_data(df, 'noncompetitive', 'ras_active', plot=True)

    # Generate the Dataset 29
    id_nr = 0
    gen_data_list = [df]
    gen_data('lowc', 'mapk, mkk, mkkk')
    gen_data('lowc', 'mkp_mapk, mkp_mkk, mkp_mkkk')
    gen_data('uncompetitive', 'mkp_mapk')
    gen_data('noncompetitive', 'mkp_mapk')
    gen_data('competitive', 'mkp_mapk')
    gen_data('lowc', 'mkp_mapk')
    gen_data('uncompetitive', 'mkp_mkk')
    gen_data('noncompetitive', 'mkp_mkk')
    gen_data('competitive', 'mkp_mkk')
    gen_data('lowc', 'mkp_mkk')
    gen_data('uncompetitive', 'mkp_mkkk')
    gen_data('noncompetitive', 'mkp_mkkk')
    gen_data('competitive', 'mkp_mkkk')
    gen_data('lowc', 'mkp_mkkk')
    gen_data('lowc', 'mapk')
    gen_data('uncompetitive', 'mkk_pp')
    gen_data('noncompetitive', 'mkk_pp')
    gen_data('competitive', 'mkk_pp')
    gen_data('lowc', 'mkk')
    gen_data('uncompetitive', 'mkkk_p')
    gen_data('noncompetitive', 'mkkk_p')
    gen_data('competitive', 'mkkk_p')
    gen_data('lowc', 'mkkk')
    gen_data('uncompetitive', 'ras_active')
    gen_data('noncompetitive', 'ras_active')
    gen_data('lowc', 'ras')
    gen_data('competitive', 'ras_active')
    gen_data('lowc', 'signal_active')
    gen_data('-', '-')
    # for i in range(20):
    #     print(i)
    #     gen_data('-', '-')
    print('*****************')

    # df = pd.concat(gen_data_list, axis=0)
    # create multiindex after full df
    # df = df.set_index(['id', 'status', 'desc', 'inhibition', 'inh_strength', 't_sec'])
    # print(df)

    # df.to_csv(f'dataset_mapk{sim1.stepps}_t.csv')
