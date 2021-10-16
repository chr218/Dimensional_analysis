# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 00:48:20 2021

@author: cvr52
"""
import numpy as np
import pandas as pd
import sympy
import sys

#%%
'''
Functions and Subroutines
'''

def make_var_str(vars_):
    #Returns comma delimited joined string.
    return ",".join(vars_)

def find_lone_var(sol_):
    chars_list =[]
    for line in sol_.args[0]:
        if str(line).isalpha():#returns True if all the characters are alphabet letters (a-z).
            chars_list.append(str(line))
    return chars_list

def pretty_output(sol_,vars_):
    for set_ in sol_:
        str_list = []
        for i,arg_ in enumerate(set_):
            if not arg_ == 0:
                SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
                str_list.append((vars_[i]+str(int(arg_))).translate(SUP))
        print(''.join(str_list))

    

#%%
'''
Read file
'''


df = pd.read_excel('input_test.xlsx')
cols = df.columns
df_arr = df.to_numpy(dtype=str, copy=False)

'''
Initialize sub-arrays
'''

var_names = df_arr[:,0]
var_sym = df_arr[:,1]
var_dim = df_arr[:,2:]

'''
Perform matrix operations
'''
#Convert numpy array from str to float
var_dim = var_dim.astype(np.float)

#Define coefficient numpy array.
A_np = var_dim.T#transpose

#Convert to sympy matrix
A_sympy = sympy.Matrix(A_np)
b_sympy = sympy.Matrix([0,0,0])

#initialize symbolic variables
#a,b,c,d,e,f,g,h = sympy.symbols("a,b,c,d,e,f,g,h")#find better way of doing this!
syms = sympy.symbols(make_var_str(var_sym))#See above comment for how this works.
syms = list(syms)

#Solve system
sol_sympy = sympy.linsolve((A_sympy,b_sympy),syms)

lone_chars = find_lone_var(sol_sympy)#Find variables which must be substituted with values. "lone_chars": V,Vo,l etc.

'''
Report Solution
'''

#Substitute exhuastive list of values for variables.
print("Variables:")
print(var_sym)
print()
print("Dimensionless Quantities:\n")
for i in range(len(lone_chars)):
    vals = np.zeros(len(lone_chars))
    vals[i] = 1

    subbed_sol_sympy = sol_sympy.subs((list(zip(lone_chars,vals))))
    pretty_output(subbed_sol_sympy,var_sym)


