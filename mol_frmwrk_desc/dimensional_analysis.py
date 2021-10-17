# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 00:48:20 2021

@author: cvr52
"""
import numpy as np
import pandas as pd
import sympy
import sys

from itertools import combinations, chain
from scipy.special import comb

#%%
'''
Functions and Subroutines
'''

def equal_ignore_order(a, b):
    unmatched = b
    
    print(a,b)

    for element in a:
        try:
            unmatched.remove(element)
        except ValueError:
            return False
    return not unmatched

def roll_arr(var_names_, var_sym_, var_dim_, roll):
    
    var_dim_col_list = []
    for col in var_dim_.T:
        var_dim_col_list.append(np.roll(col,roll))
        
    var_names_rolled = np.roll(var_names_,roll)
    var_sym_rolled = np.roll(var_sym_,roll)
    var_dim_rolled = np.stack(var_dim_col_list,axis=1)
    
    return var_names_rolled, var_sym_rolled, var_dim_rolled
        

def comb_index(n, k):
    count = comb(n, k, exact=True)
    index = np.fromiter(chain.from_iterable(combinations(range(n), k)), 
                        int, count=count*k)
    return index.reshape(-1, k)

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

    tot_var_list_ = []
    for set_ in sol_:
        str_list = []
        var_list_ = []
        for i,arg_ in enumerate(set_):
            if not arg_ == 0:
                
                SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
                str_list.append((vars_[i]+str(round(arg_))).translate(SUP))
                var_list_.append(vars_[i])
                
        tot_var_list_.append(var_list_)
        output_str_ = (''.join(str_list))
        
    return tot_var_list_, output_str_
            
        
    
def find_dimensionless(var_dim_,var_sym_):
    #Convert numpy array from str to float
    var_dim_ = var_dim_.astype(np.float)
    
    #Define coefficient numpy array.
    A_np = var_dim_.T#transpose
    
    #Convert to sympy matrix
    A_sympy = sympy.Matrix(A_np)
    b_sympy = sympy.Matrix([0,0,0])
    
    #initialize symbolic variables
    syms = sympy.symbols(make_var_str(var_sym_))#See above comment for how this works.
    syms = list(syms)
    
    #Solve system
    sol_sympy_ = sympy.linsolve((A_sympy,b_sympy),syms)
    
    lone_chars_ = find_lone_var(sol_sympy_)#Find variables which must be substituted with values. "lone_chars": V,Vo,l etc.
    
    return sol_sympy_, lone_chars_, var_sym_


if __name__ == "__main__":
    
    '''
    Read file
    '''
    
    
    df = pd.read_excel('input_Greenberg.xlsx')
    cols = df.columns
    df_arr = df.to_numpy(dtype=str, copy=False)
    
    '''
    Initialize sub-arrays
    '''
    
    var_names = df_arr[:,0]
    var_sym = df_arr[:,1]
    var_dim = df_arr[:,2:]    

    #Roll each array. (i.e shift position of array elements by n places,
    # such that the final element goes in front.) This is necessary to find
    # exhuastive list of dimensionless quantities.
    
    unique_vars_list = []
    for iter_ in range(len(var_sym)):
        
        var_names_tmp, var_sym_tmp, var_dim_tmp = roll_arr(var_names,var_sym,var_dim,iter_)
    
    
        '''
        Perform matrix operations
        '''
    
        sol_sympy, lone_chars, var_sym_tmp = find_dimensionless(var_dim_tmp,var_sym_tmp)
        
        
        '''
        Report Solution
        '''
        
        #Obtain exhuastive list of indices for "vals" variable.
        #i.e within Greenberg pg 503, we set alpha_1 =1, ..., alpha_N = 0
        #But we can also set alpha_1 =1 , alpha_2 = 1, ..., alpha_N = 0.
        #Therefore the "comb_index" function returns all unique combinations
        # of indices for setting alphas = 1.
        
        idx = comb_index(len(lone_chars), 1)
        
        #Substitute list of values for variables.
        print("Variables:")
        print(var_sym_tmp)
        print()
        print("Dimensionless Quantities:\n")
        
        #for i in range(len(lone_chars)):
            #vals = np.zeros(len(lone_chars))
            #vals[i] = 1
    
    
        for pair in idx:
            vals = np.zeros(len(lone_chars))
            
            for element in pair:
                vals[element] = 1
            subbed_sol_sympy = sol_sympy.subs((list(zip(lone_chars,vals))))
            tot_vars_list, output_str = pretty_output(subbed_sol_sympy,var_sym_tmp)
            

            
            
            #print(equal_ignore_order(tot_vars_list,unique_vars_list))
            if not tot_vars_list in unique_vars_list:
                unique_vars_list.append(tot_vars_list)
                print(output_str)
             
            
    #print(unique_vars_list)
            
            


