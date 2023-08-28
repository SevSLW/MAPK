import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.colors as mcolors
from sklearn.metrics import confusion_matrix
import os
import ast

df_list_inh = []
df_list_n = []
for datei_name in os.listdir():
    print(datei_name)
    if datei_name.endswith('.csv'):
        if 'n25' in datei_name:
            df = pd.read_csv(datei_name, index_col='name')
            df_list_inh.append(df)
        elif 'm25' in datei_name:
            df_nonoise = pd.read_csv(datei_name)
        else:
            df = pd.read_csv(datei_name, index_col='name')
            df_list_n.append(df)
df = pd.read_csv('mapk361_n25_inh0.csv', index_col='name')
df_list_n.append(df)

df_inh = pd.concat(df_list_inh)
df_n = pd.concat(df_list_n)

df = pd.read_csv('mapk361_n100_inh0.csv', sep=',')
for i, col in enumerate(df_inh.columns):
    print(i, '  ', col)


df_prec = df_inh.iloc[:, 1:15]
df_prec = df_prec.reset_index()
df_prec = df_prec.drop('name', axis=1)
print(df_prec.at[0, '-'])
inh_strengts = np.arange(0, 1.2, 0.2)
inh_strengts += 1
print(inh_strengts)
inh_strengts = inh_strengts.round(1)
label = []
for col in df_prec:
    if df_prec.at[0, col] < 0.8:
        print(col)
        label.append(col)
        sns.lineplot(x=inh_strengts, y=df_prec[col], data=df_inh, label=col)
plt.legend(bbox_to_anchor=(1, 0.86), loc='upper right')
plt.xlabel('Inhibitorfaktor')
plt.ylabel('Precision score')
plt.show()

print('INH')
prec_latex = ''
for i in range(len(df_prec)):
    for col in df_prec:
        prec = df_prec.at[i, col]
        prec = np.round(prec, 2)
        prec_latex += f'{prec}&'
    prec_latex += r'\\'
    prec_latex += '\n'
print('prec_latex')
print(prec_latex)
for i in df_prec.columns:
    print(i, end='&')



print()
print('********* df_n')
print(df_n)

df_prec = df_n.iloc[:, 1:15]
df_prec = df_prec.reset_index()
df_prec = df_prec.drop('name', axis=1)
print(df_prec.at[0, '-'])
n = [1, 2, 3, 4, 5, 6]
inh_strengts = inh_strengts.round(1)
label = []
for col in df_prec:
    if col in ['-', 'mkk', 'mkp_mapk', 'ras']:
        print(col)
        label.append(col)
        ax = sns.lineplot(x=n, y=df_prec[col], data=df_inh, label=col)
plt.legend(bbox_to_anchor=(1, 0.86), loc='upper right')
plt.xlabel('n')
plt.ylabel('Precision score')
plt.xticks(n, labels=['2', '5', '10', '25', '50', '100'])
plt.show()

print()
print('PREC N')
prec_latex = ''
for i in range(len(df_prec)):
    for col in df_prec:
        prec = df_prec.at[i, col]
        prec = np.round(prec, 2)
        prec_latex += f'{prec}&'
    prec_latex += r'\\'
    prec_latex += '\n'
print('prec_latex')
print(prec_latex)
for i in df_prec.columns:
    print(i, end='&')

#%%
dataframe = pd.read_csv('mapk361_n25_inh1.csv')
print('nonoise')
print(dataframe['cm'])
print(dataframe.at[0, 'cm'])
cm = dataframe.at[0, 'cm']
matrix_string = ''

n = 1
matrix = []
row = [n]
number = ''
for char in cm:
    if char.isnumeric():
        number += char
    elif char is ' ':
        if len(number) > 0:
            row.append(int(number))
            number = ''
    elif char is ']':
        if len(number) > 0:
            row.append(int(number))
            number = ''
        matrix.append(row)
        n += 1
        row = [n]
print(matrix)

matrix_string = r''
for i in range(len(matrix)):
    matrix_string += f'{i}&'
matrix_string = matrix_string.rstrip('&')
matrix_string += r'\\'
matrix_string += '\n'
matrix_string += r'\hline'
matrix_string += '\n'
row_str = r''
for row in matrix:
    for number in row:
        row_str += str(number)
        row_str += '&'
    row_str = row_str.rstrip('&')
    row_str += r'\\'
    row_str += '\n'
    matrix_string += row_str
    row_str = ''
print(matrix_string)

prec_latex = ''
prec_noise = df_nonoise.iloc[:, 2:16]
prec_noise = prec_noise.reset_index()
# prec_noise = prec_noise.drop('name', axis=1)
for col in prec_noise:
    prec = prec_noise.at[0, col]
    prec = np.round(prec, 2)
    prec_latex += f'{prec}&'
prec_latex += r'\\'
prec_latex += '\n'

print(prec_latex)

print(prec_latex)
print(df_nonoise.columns)
print(df_nonoise.iloc[:, 2:15])





