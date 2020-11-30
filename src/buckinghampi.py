"""buckinghampi.py: a symbolic module that generates the pi terms based on some variables by applying the pi-theorem."""

__author__ = "Mokbel Karam"
__copyright__ = "Copyright (c) 2020, Mokbel Karam"

__credits__ = ["University of Utah Department of Chemical Engineering"]
__license__ = "Apache 2.0"
__version__ = "0.1.2"
__maintainer__ = "Mokbel Karam"
__email__ = "mokbel.karam@chemeng.utah.edu"
__status__ = "Production"

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from sympy.core.mul import Mul, Pow
from sympy.core.expr import Expr
import numpy as np
from itertools import combinations

class BuckinghamPi:
    def __init__(self, physical_dimensions:str,sep=' '):
        '''
        Construct an instance of the BuckinghamPi theorem
        :param physical_dimensions: (string) of the physical dimensions used for this instance. ex: physical_dimensions='m l t'
        :param sep: (string) used to separate the physical_dimensions, by default it is a white space ' '
        '''
        self.__all_physical_dimensions = ('a','k','t','l','m','cd','mol')
        physical_dimensions_list = [x.lower() for x in physical_dimensions.split(sep=sep)]
        if not(all(x in self.__all_physical_dimensions for x in physical_dimensions_list)):
            raise Exception('physical_dimensions has to be a subset of the all physical dimensions {}'.format(self.__all_physical_dimensions))

        self.__physical_dimensions = {v:sp.symbols(v) for v  in physical_dimensions_list}
        self.__physical_dimensions_list = [self.__physical_dimensions[key] for key in self.__physical_dimensions.keys()]
        self.num_physical_dimensions = len(self.__physical_dimensions_list)

        self.__var_from_idx={}
        self.__idx_from_var = {}
        self.__variables={}
        self.__sym_variables={}
        self.__flagged_var = {'var_name':None, 'var_index':None,'selected':False}

        self.__null_spaces = []

    @property
    def physical_dimensions(self):
        '''
        :return: the physical dimensions being used
        '''
        return self.__physical_dimensions

    @property
    def all_physical_dimensions(self):
        '''
        :return: all possible physical dimensions. they are 7 in total:  ampere (a), kelvin (k), second (t), metre (l), kilogram (m), candela (cd) and mole (mol)
        '''
        return self.__all_physical_dimensions

    @property
    def variables(self):
        '''
        :return: a dict of the variables added by the user.
        '''
        return self.__variables


    def __parse_expression(self,string:str):
        expr = parse_expr(string.lower(),local_dict=self.physical_dimensions,global_dict={"Integer":sp.Integer})
        if not (isinstance(expr,Mul) or isinstance(expr,Pow) or isinstance(expr,sp.Symbol)):
            raise Exception('expression of type {} is not of the accepted types ({}, {}, {})'.format(type(expr), Mul, Pow, sp.Symbol))
        if expr.as_coeff_Mul()[0] != 1:
            raise Exception('cannot have coefficients, {}, that multiply the expression'.format(expr.as_coeff_Mul()[0]))

        return expr

    def __extract_exponents(self,expr:Expr):
        vect = np.zeros(self.num_physical_dimensions)
        args = list(expr.args) if list(expr.args) else [expr]
        # print(args)
        if isinstance(expr, Pow):
            vect[self.__physical_dimensions_list.index(args[0])] = int(args[1])
        else:
            for e in args:
                if isinstance(expr, sp.Symbol):
                    vect[self.__physical_dimensions_list.index(e)]= int(1)
                    # print('({}, {})'.format(e, 1))
                else:
                    var, exponent= e.as_base_exp()
                    vect[self.__physical_dimensions_list.index(var)] = int(exponent)
                    # print('({}, {})'.format(var, exponent))

        return vect

    def add_variable(self, name:str, expression:str, select=False):
        '''
        Add variables to use for the pi-theorem
        :param name: (string) name of the variable to be added
        :param expression: (string) expression of the independent physical variable expressed in terms of the k independent physical units.
        :return: (Boolean) True if done perfectly
        '''
        expr =  self.__parse_expression(expression)
        self.__variables.update({name:expr})
        var_idx = len(list(self.__variables.keys()))-1
        self.__var_from_idx[var_idx]= name
        self.__idx_from_var[name] = var_idx
        if select and (self.__flagged_var['selected']==False):
            self.__flagged_var['var_name'] = name
            self.__flagged_var['var_index'] = var_idx
            self.__flagged_var['selected'] = True
        elif select and (self.__flagged_var['selected']==True):
            raise Exception("you cannot select more than one variable at a time.")

        return True

    def __create_M(self):
        self.num_variable = len(list(self.__variables.keys()))
        if self.num_variable < self.num_physical_dimensions:
            raise Exception('The number of variables has to be greater than the number of physical dimensions.')

        self.M = np.zeros(shape=(self.num_variable, self.num_physical_dimensions))
        # fill M
        for var_name in self.__variables.keys():
            expr = self.__variables[var_name]
            vect = self.__extract_exponents(expr)
            row = self.__idx_from_var[var_name]
            self.M[row, :] = vect

        self.M = self.M.transpose()

    def __create_symbolic_variables(self):
        for var_name in self.__variables.keys():
            self.__sym_variables[var_name] = sp.symbols(var_name)

    def __solve_null_spaces(self):
        if self.__flagged_var['selected']==True:
            self.__solve_null_spaces_for_flagged_variables()

        else:
            for idx in self.__var_from_idx.keys():
                self.__flagged_var['var_name'] = self.__var_from_idx[idx]
                self.__flagged_var['var_index'] = idx
                self.__flagged_var['selected'] = True

                self.__solve_null_spaces_for_flagged_variables()

    def __solve_null_spaces_for_flagged_variables(self):

        assert self.__flagged_var['selected']==True, " you need to select a variable"

        n = self.num_variable
        m = self.num_physical_dimensions

        original_indicies = list(range(0, n))
        all_idx = original_indicies.copy()
        if self.__flagged_var['selected']:
            del all_idx[self.__flagged_var['var_index']]

        # print(all_idx)
        all_combs = list(combinations(all_idx,m))
        # print(all_combs)

        num_det_0 = 0
        for comb in all_combs:
            temp_comb = list(comb).copy()
            extra_vars = [i for i in original_indicies if i not in temp_comb ]
            b_ns = []
            for extra_var in extra_vars:
                new_order = {}
                temp_comb.append(extra_var)
                A = self.M[:,temp_comb].copy()
                for num,var_idx in enumerate(temp_comb):
                    new_order[num] =  self.__var_from_idx[var_idx]
                B = sp.Matrix(A)
                test_mat = B[:,:m]
                if sp.det(test_mat) !=0:
                    ns = B.nullspace()[0]
                    b_ns.append({'order': new_order, 'power': ns.tolist()})

                else:
                    num_det_0+=1
                temp_comb = list(comb).copy()
            self.__null_spaces.append(b_ns)
        # print("num of det 0 : ",num_det_0)

    def __construct_symbolic_pi_terms(self):
        self.__allpiterms = []
        for space in self.__null_spaces:
            spacepiterms = []
            for term in space:
                expr = 1
                idx = 0
                for order,power in zip(term['order'].keys(),term['power']):
                    expr *= self.__sym_variables[term['order'][order]] ** sp.nsimplify(sp.Rational(power[0]))
                    idx += 1
                spacepiterms.append(expr)
            # check for already existing pi terms in previous null-spaces
            already_exists = False
            for previouspiterms in self.__allpiterms:
                if all(x in previouspiterms for x in spacepiterms):
                    already_exists = True
                    break
            if not already_exists:
                self.__allpiterms.append(spacepiterms)

    def generate_pi_terms(self):
        '''
        Generates all the possible pi terms
        :return: Boolean true if done perfectly
        '''
        self.__create_M()

        self.__create_symbolic_variables()

        self.__solve_null_spaces()

        self.__construct_symbolic_pi_terms()

        return True

    @property
    def pi_terms(self):
        '''
        :return: a list with all the symbolic dimensionless terms for all permutation of the dimensional Matrix M
        '''
        return self.__allpiterms
