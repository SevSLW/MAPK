import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.colors as mcolors
from sklearn.metrics import confusion_matrix
import os
import ast

def gen_labellist(columns):
    global df
    label_list =[]
    for value, group_df in df.groupby(level='id'):
        inner = []
        for col in columns.split(','):
            id_series = group_df[col]
            first_element = id_series.iloc[0]
            inner.append(first_element)
        label_list.append(inner)
    label_list = np.array(label_list)
    for col in columns.split(','):
        df = df.drop(col, axis=1)
    return label_list

data = 'mapk361_n25_inh0_nonoise.csv'
df_pred = pd.read_csv('mapk361_m25_inh0_nonoise.csv')
df = pd.read_csv('../../dataset_mapk361_n25_inh0_nonoise_h25_test.csv', sep=',', index_col=['id', 't_sec'])
y_predict = df_pred['y']
print(y_predict)
y_pred = y_predict.iloc[0]
y_pred = y_pred.strip('[]')
y_pred = y_pred.replace('\n', '')

n = 0
y_pred_list = []
desc = ''
for i, char in enumerate(y_pred):
    if char == "'":
        n += 1
        if n % 2 == 0:
            y_pred_list.append(desc[1:])
            desc = ''
        continue
    desc += char
y_pred_list = np.array(y_pred_list)

status_array = gen_labellist('desc')
y_test = np.array([status_array[i][0] for i, val in enumerate(status_array)])

n = 0
id_list = []
for i, value in enumerate(y_pred_list):
    if y_test[i] != y_pred_list[i]:
        if y_pred_list[i] == '-':
            n += 1
            id_list.append(i)
            print(n)
            print(i)
            print(y_test[i])
            print(y_pred[i])
            print()
print(id_list)

# fault_id_list = [
#     233, 248, 326, 327, 334, 469, 573, 636, 642, 646
# ]

fault_id_list = id_list


df = pd.read_csv('../../dataset_mapk361_n25_inh0_nonoise_h25_test.csv', sep=',', index_col=['id', 't_sec'])
fault_tseries_list = []

for i in fault_id_list:
    tseries = df.xs(i, level='id', drop_level=False)
    fault_tseries_list.append(tseries)

df_fault_tseries = pd.concat(fault_tseries_list)
print('unique desc')
print(df_fault_tseries['inhibition'].unique())
df_fault_tseries = df_fault_tseries.reset_index()
df_fault_tseries = df_fault_tseries.set_index(['id', 'desc'])
# df_fault_tseries = df_fault_tseries.drop(0, level='id')

# data = 'dataset_mapk361_n25_inh0_nonoise_h25_train.csv'
# df = pd.read_csv(data, sep=',', index_col=['id', 'desc'])
# df = df.drop('status', axis=1)
# df = df.drop('inhibition', axis=1)
# df = df.drop('inh_strength', axis=1)
# df = df.drop(0, level='id')
# # df = df.drop('desc', axis=1)
# print(df)

for i, row in df_fault_tseries.groupby(level='desc'):
    print('i')
    print(i)
    for j, innerrow in row.groupby(level='id'):
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
    plt.show()
