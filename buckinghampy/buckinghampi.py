"""buckinghampi.py: a symbolic module that generates the pi terms based on some variables by applying the pi-theorem."""

__author__ = "Mokbel Karam"
__copyright__ = "Copyright (c) 2021, Mokbel Karam"

__credits__ = ["University of Utah Department of Chemical Engineering"]
__license__ = "MIT"
__version__ = "1.0.3"
__maintainer__ = "Mokbel Karam"
__email__ = "karammokbel@gmail.com"
__status__ = "Production"

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from sympy.core.mul import Mul, Pow
from sympy.core.expr import Expr
import numpy as np
from itertools import combinations,permutations
from tabulate import tabulate

try:
    from IPython.display import display, clear_output, Math, Markdown
except:
    pass

class BuckinghamPi:
    def __init__(self):
        '''
        Construct an instance of the BuckinghamPi theorem
        '''
        self.__var_from_idx={}
        self.__idx_from_var = {}
        self.__variables={}
        self.__sym_variables={}
        self.__flagged_var = {'var_name':None, 'var_index':None,'selected':False}

        self.__null_spaces = []

        self.__fundamental_vars_used = [] # list of fundamental variables being used

        self.__prefixed_dimensionless_terms = []

    @property
    def fundamental_variables(self):
        '''
        :return: a list of the fundamental variables being used
        '''
        return self.__fundamental_vars_used

    @property
    def variables(self):
        '''
        :return: a dict of the variables added by the user.
        '''
        return self.__variables


    def __parse_expression(self,string:str):
        if '^' in string:
            # convert the xor operator to power operator
            string = string.replace('^','**')

        expr = parse_expr(string.lower())

        if not (isinstance(expr,Mul) or isinstance(expr,Pow) or isinstance(expr,sp.Symbol)):
            raise Exception('expression of type {} is not of the accepted types ({}, {}, {})'.format(type(expr), Mul, Pow, sp.Symbol))
        if expr.as_coeff_Mul()[0] != 1:
            raise Exception('cannot have coefficients, {}, that multiply the expression {}'.format(expr.as_coeff_Mul()[0],expr.as_coeff_Mul()[1]))

        #extract the physical dimensions from the dimensions expressions
        used_symbols = list(expr.free_symbols)
        for sym in used_symbols:
            if not sym in self.__fundamental_vars_used:
                self.__fundamental_vars_used.append(sym)

        return expr

    def __extract_exponents(self,expr:Expr):
        num_physical_dimensions = len(self.__fundamental_vars_used)
        vect = np.zeros(num_physical_dimensions)
        args = list(expr.args) if list(expr.args) else [expr]
        # print(args)
        if isinstance(expr, Pow):
            vect[self.__fundamental_vars_used.index(args[0])] = int(args[1])
        else:
            for e in args:
                if isinstance(expr, sp.Symbol):
                    vect[self.__fundamental_vars_used.index(e)]= int(1)
                    # print('({}, {})'.format(e, 1))
                else:
                    var, exponent= e.as_base_exp()
                    vect[self.__fundamental_vars_used.index(var)] = int(exponent)
                    # print('({}, {})'.format(var, exponent))

        return vect

    def add_variable(self, name: str, dimensions: str, non_repeating=False):
        '''
        Add variables to use for the pi-theorem
        :param name: (string) name of the variable to be added
        :param dimensions: (string) expression of the independent physical variable expressed in terms of the k independent fundamental dimensions.
        :param non_repeating: (boolean) select a variable to belong to the non-repeating variables matrix. This will ensure that the selected variable
                                        only shows up in one dimensionless group.
        '''
        if dimensions!="1":
            expr =  self.__parse_expression(dimensions)
            self.__variables.update({name:expr})
            var_idx = len(list(self.__variables.keys()))-1
            self.__var_from_idx[var_idx]= name
            self.__idx_from_var[name] = var_idx
            if non_repeating and (self.__flagged_var['selected'] == False):
                self.__flagged_var['var_name'] = name
                self.__flagged_var['var_index'] = var_idx
                self.__flagged_var['selected'] = True
            elif non_repeating and (self.__flagged_var['selected'] == True):
                raise Exception("you cannot select more than one variable at a time to be a non_repeating.")
        else:
            self.__prefixed_dimensionless_terms.append(sp.symbols(name))

    def __create_M(self):
        self.num_variable = len(list(self.__variables.keys()))
        num_physical_dimensions = len(self.__fundamental_vars_used)
        if self.num_variable <= num_physical_dimensions:
            raise Exception('The number of variables has to be greater than the number of physical dimensions.')

        self.M = np.zeros(shape=(self.num_variable, num_physical_dimensions))
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

        assert self.__flagged_var['selected']==True, " you need to select a variable to be explicit"

        n = self.num_variable
        m = len(self.__fundamental_vars_used)

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
            if b_ns: # if b_ns is not empty add it to the nullspaces list
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

    def __rm_duplicated_powers(self):
        # this algorithm rely on the fact that the nullspace function
        # in sympy set one free variable to 1 and the all other to zero
        # then solve the system by back substitution.
        duplicate = []
        dummy_other_terms = self.__allpiterms.copy()
        for num_set, pi_set in enumerate(self.__allpiterms):
            dummy_other_terms.remove(pi_set)
            for num_other, other in enumerate(dummy_other_terms):
                permutations_sets = permutations(pi_set)
                for p_set in permutations_sets:
                    # create a permutation vector from the permutation set
                    p_V = sp.Matrix(list(p_set))
                    # create a vector from the other set of dimensionless groups that we are comparing to.
                    o_V = sp.Matrix(other)
                    # create an element wise inverse of the vector of dimensionless groups
                    o_V_inv = o_V.applyfunc(lambda x:x**(-1))

                    result = sp.matrix_multiply_elementwise(p_V, o_V)
                    # obtain the index of numerical value in the result vector.
                    # numerical values indicates that one dimensionless group is the inverse of the other group
                    # in this algorithm the numerical value will be equal to 1 (this is a result of the nullspace function in sympy)
                    idx_num_result = [x for x in range(len(p_set)) if isinstance(result[x,0],sp.Number)]
                    # also repeat the multiplication with the inverse vector
                    result_inv = sp.matrix_multiply_elementwise(p_V, o_V_inv)
                    # check for the index of the numerical values in the result vector
                    idx_num_result_inv = [x for x in range(len(p_set)) if isinstance(result_inv[x,0],sp.Number)]
                    # concatinate the indices into one list
                    all_indices = idx_num_result + idx_num_result_inv
                    # compare if the two vector are duplicates
                    if set(all_indices) == set(list(range(len(p_set)))):
                        duplicate.append(pi_set)

        # remove duplicates from the main dict of all pi terms
        for dup in duplicate:
            if dup in self.__allpiterms:
                self.__allpiterms.remove(dup)
        return duplicate

    def __populate_prefixed_dimensionless_groups(self):
        for num_set, pi_set in enumerate(self.__allpiterms):
            for pre_fixed_dimensionless_group in self.__prefixed_dimensionless_terms:
                self.__allpiterms[num_set].append(pre_fixed_dimensionless_group)

    def generate_pi_terms(self):
        '''
        Generates all the possible pi terms
        '''
        self.__create_M()

        self.__create_symbolic_variables()

        self.__solve_null_spaces()

        self.__construct_symbolic_pi_terms()

        self.__rm_duplicated_powers()

        self.__populate_prefixed_dimensionless_groups()

    @property
    def pi_terms(self):
        '''
        :return: a list with all the symbolic dimensionless terms for all permutation of the dimensional Matrix M
        '''
        return self.__allpiterms


    def __Jupyter_print(self):
        ''' print the rendered Latex format in Jupyter cell'''
        for set_num, space in enumerate(self.__allpiterms):
            latex_str= '\\text{Set }'
            latex_str+='{}: \\quad'.format(set_num+1)
            for num, term in enumerate(space):
                latex_str += '\\pi_{} = '.format(num+1)+sp.latex(term)
                latex_str += '\\quad'
            display(Math(latex_str))
            display(Markdown('---'))

    def __tabulate_print(self,latex_string=False):
        ''' print the dimensionless sets in a tabulated format'''

        latex_form = []
        for pi_set in self.__allpiterms:
            latex_set = []
            for pi in pi_set:
                if latex_string:
                    if latex_string:
                        latex_set.append(sp.latex(pi))
                    else:
                        latex_set.append(pi)
                else:
                    latex_set.append(pi)
            latex_form.append(latex_set)

        num_of_pi_terms = len(latex_form[0])

        headers = ['sets']
        for num in range(num_of_pi_terms):
            headers.append('Pi {}'.format(num + 1))

        for num, set in enumerate(latex_form):
            set.insert(0, num + 1)

        print(tabulate(latex_form, headers=headers))

    def print_all(self, latex_string=False):
        '''
        print all the sets of dimensionless groups in latex or symbolic form.
        :latex_string: optional boolean. If set to True the function will print the latex string of the
                        dimensionless groups. if set to False the function will print the symbolic form of the
                        dimensionless groups.
        '''
        try:
            ''' Try to render the latex in Jupyter cell'''
            self.__Jupyter_print()
        except:
            ''' print the dimensionless sets in a tabulated format when in terminal session'''
            self.__tabulate_print(latex_string)