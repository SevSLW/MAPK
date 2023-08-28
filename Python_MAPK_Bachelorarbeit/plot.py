import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.colors as mcolors
import os

import mapk_sim_neu

#%%
simm = mapk_sim_neu.Simulation()
simm.mkp_mkk.kinetic_factors['i_nc'] = 1.5
simm.scipy_rk45(order=True)
print(simm.mkk.kinetic_factors)
simm.plot(['MKKK', 'MKKK_P', 'MKK', 'MKK_P', 'MKK_PP', 'MAPK', 'MAPK_P', 'MAPK_PP', 'Signal', 'RAS', 'RAS_A'])
plt.show()

#%%
data = 'dataset_mapk361_n25_inh0_nonoise_h25_train.csv'
df = pd.read_csv(data, sep=',', index_col=['id', 'desc'])
df = df.drop('status', axis=1)
df = df.drop('inhibition', axis=1)
df = df.drop('inh_strength', axis=1)
df = df.drop(0, level='id')
# df = df.drop('desc', axis=1)
print(df)

for i, row in df.groupby(level='desc'):
    print('i')
    print(i)
    for j, innerrow in row.groupby(level='id'):
        print(j)
        list = ['MKKK', 'MKKK_P', 'MKK', 'MKK_P', 'MKK_PP', 'MAPK', 'MAPK_P', 'MAPK_PP', 'Signal', 'RAS', 'RAS_A']
        col_list = [f for f in mcolors.TABLEAU_COLORS]
        col_list.append('lime')
        for n, val in enumerate(list):
            sns.lineplot(x=innerrow['t_sec'], y=innerrow[val], data=df, color=col_list[n])
        plt.xlabel('t [min]')
        plt.ylabel('c [nM]')
        plt.title(i)
    output_folder = os.path.join('plots', data.strip(".csv"))
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    path = os.path.join(output_folder, f'{i}.png')
    plt.savefig(path)
    # plt.show()

print('**** FIN ****')
