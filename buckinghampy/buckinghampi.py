"""buckinghampi.py: a symbolic module that generates the pi terms based on some variables by applying the pi-theorem."""

__author__ = "Mokbel Karam"
__copyright__ = "Copyright (c) 2021, Mokbel Karam"

__credits__ = ["University of Utah Department of Chemical Engineering"]
__license__ = "MIT"
__version__ = "1.0.4"
__maintainer__ = "Mokbel Karam"
__email__ = "karammokbel@gmail.com"
__status__ = "Production"

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from sympy.core.mul import Mul, Pow
from sympy.core.expr import Expr
import numpy as np
from fractions import Fraction
from itertools import combinations
from tabulate import tabulate

# scipy.linalg.null_space uses a robust SVD-based algorithm and is orders of
# magnitude faster than sympy's exact rational nullspace for our purposes.
from scipy.linalg import null_space as scipy_null_space

try:
    from IPython.display import display, clear_output, Math, Markdown
except:
    pass


class BuckinghamPi:
    def __init__(self):
        '''
        Construct an instance of the BuckinghamPi theorem
        '''
        self.__var_from_idx = {}
        self.__idx_from_var = {}
        self.__variables = {}
        self.__sym_variables = {}
        self.__flagged_var = {'var_name': None, 'var_index': None, 'selected': False}

        self.__null_spaces = []

        self.__fundamental_vars_used = []  # list of fundamental variables being used

        self.__prefixed_dimensionless_terms = []

        self.__flagged_var_max_sets = 20

    @property
    def fundamental_variables(self):
        return self.__fundamental_vars_used

    @property
    def variables(self):
        return self.__variables

    def __parse_expression(self, string: str):
        if '^' in string:
            string = string.replace('^', '**')

        expr = parse_expr(string.lower())

        if not (isinstance(expr, Mul) or isinstance(expr, Pow) or isinstance(expr, sp.Symbol)):
            raise Exception('expression of type {} is not of the accepted types ({}, {}, {})'.format(
                type(expr), Mul, Pow, sp.Symbol))
        if expr.as_coeff_Mul()[0] != 1:
            raise Exception('cannot have coefficients, {}, that multiply the expression {}'.format(
                expr.as_coeff_Mul()[0], expr.as_coeff_Mul()[1]))

        used_symbols = list(expr.free_symbols)
        for sym in used_symbols:
            if not sym in self.__fundamental_vars_used:
                self.__fundamental_vars_used.append(sym)

        return expr

    def __extract_exponents(self, expr: Expr):
        num_physical_dimensions = len(self.__fundamental_vars_used)
        vect = np.zeros(num_physical_dimensions)
        args = list(expr.args) if list(expr.args) else [expr]
        if isinstance(expr, Pow):
            vect[self.__fundamental_vars_used.index(args[0])] = int(args[1])
        else:
            for e in args:
                if isinstance(expr, sp.Symbol):
                    vect[self.__fundamental_vars_used.index(e)] = int(1)
                else:
                    var, exponent = e.as_base_exp()
                    vect[self.__fundamental_vars_used.index(var)] = int(exponent)
        return vect

    def add_variable(self, name: str, dimensions: str, non_repeating=False):
        '''
        Add variables to use for the pi-theorem
        :param name: (string) name of the variable to be added
        :param dimensions: (string) expression of the independent physical variable expressed in terms of
                           the k independent fundamental dimensions.
        :param non_repeating: (boolean) select a variable to belong to the non-repeating variables matrix.
                              This will ensure that the selected variable only shows up in one dimensionless group.
        '''
        if dimensions != "1":
            expr = self.__parse_expression(dimensions)
            self.__variables.update({name: expr})
            var_idx = len(list(self.__variables.keys())) - 1
            self.__var_from_idx[var_idx] = name
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
        if self.__flagged_var['selected'] == True:
            self.__solve_null_spaces_for_flagged_variables()
        else:
            for idx in self.__var_from_idx.keys():
                self.__flagged_var['var_name'] = self.__var_from_idx[idx]
                self.__flagged_var['var_index'] = idx
                self.__flagged_var['selected'] = True
                self.__solve_null_spaces_for_flagged_variables()

    def __rationalize_vector(self, vec, max_denominator=1000):
        """Convert a float numpy vector to a list of sp.Rational using stdlib fractions."""
        result = []
        for v in vec:
            frac = Fraction(float(v)).limit_denominator(max_denominator)
            result.append(sp.Rational(frac.numerator, frac.denominator))
        return result

    def __solve_null_spaces_for_flagged_variables(self):
        assert self.__flagged_var['selected'] == True, "you need to select a variable to be explicit"

        n = self.num_variable
        m = len(self.__fundamental_vars_used)

        original_indices = list(range(0, n))
        all_idx = original_indices.copy()
        if self.__flagged_var['selected']:
            del all_idx[self.__flagged_var['var_index']]

        all_combs = list(combinations(all_idx, m))[:self.__flagged_var_max_sets]

        num_all_pi_terms = (n - m) * len(all_combs)
        num_singular = 0

        for comb in all_combs:
            temp_comb = list(comb).copy()
            extra_vars = [i for i in original_indices if i not in temp_comb]
            b_ns = []

            for extra_var in extra_vars:
                new_order = {}
                temp_comb.append(extra_var)
                A = self.M[:, temp_comb].copy()  # shape: (m, n)

                for num, var_idx in enumerate(temp_comb):
                    new_order[num] = self.__var_from_idx[var_idx]

                # FAST singularity check via numpy determinant (float arithmetic)
                test_mat = A[:, :m]
                if abs(np.linalg.det(test_mat)) > 1e-10:
                    # FAST null space via scipy SVD — much faster than sympy rational arithmetic
                    ns_float = scipy_null_space(A)  # shape: (n_cols, nullity)

                    for col in range(ns_float.shape[1]):
                        col_vec = ns_float[:, col]
                        # Normalize so last-column component is positive (matches SymPy convention)
                        pivot = col_vec[-1]
                        if abs(pivot) > 1e-12:
                            col_vec = col_vec / pivot
                        rational_vec = self.__rationalize_vector(col_vec)
                        b_ns.append({'order': new_order, 'power': [[v] for v in rational_vec]})
                else:
                    num_singular += 1

                temp_comb = list(comb).copy()

            if b_ns:
                self.__null_spaces.append(b_ns)

        if num_singular == num_all_pi_terms:
            raise Exception(
                "All the P matrices in the possible sets of dimensionless groups were singular, resulting in no pi terms.")

    def __construct_symbolic_pi_terms(self):
        self.__allpiterms = []
        for space in self.__null_spaces:
            spacepiterms = []
            for term in space:
                expr = 1
                for order, power in zip(term['order'].keys(), term['power']):
                    # power[0] is already sp.Rational — no nsimplify needed
                    expr *= self.__sym_variables[term['order'][order]] ** power[0]
                spacepiterms.append(expr)
            self.__allpiterms.append(spacepiterms)

    def __rm_duplicated_powers(self):
        """
        Hash-based deduplication: O(n_sets) instead of O(n_sets^2 * k!).

        Each pi-term set is represented as a frozenset of canonical keys where
        each key is the lexicographically smaller of latex(term) and latex(1/term).
        This makes the check invariant to both ordering within the set and
        inversion of individual terms — the same equivalences the original
        permutation-based algorithm was checking.
        """
        seen = set()
        unique = []
        for pi_set in self.__allpiterms:
            key_parts = set()
            for term in pi_set:
                s = sp.latex(term)
                s_inv = sp.latex(term ** -1)
                key_parts.add(min(s, s_inv))
            key = frozenset(key_parts)
            if key not in seen:
                seen.add(key)
                unique.append(pi_set)
        self.__allpiterms = unique

    def __populate_prefixed_dimensionless_groups(self):
        for num_set, pi_set in enumerate(self.__allpiterms):
            for pre_fixed_dimensionless_group in self.__prefixed_dimensionless_terms:
                self.__allpiterms[num_set].append(pre_fixed_dimensionless_group)

    def generate_pi_terms(self):
        '''
        Generates all the possible pi terms.
        Note: this function can throw exceptions.
        '''
        self.__create_M()
        self.__create_symbolic_variables()
        self.__solve_null_spaces()
        self.__construct_symbolic_pi_terms()
        self.__rm_duplicated_powers()
        self.__populate_prefixed_dimensionless_groups()

    @property
    def pi_terms(self):
        return self.__allpiterms

    def __Jupyter_print(self):
        for set_num, space in enumerate(self.__allpiterms):
            latex_str = '\\text{Set }'
            latex_str += '{}: \\quad'.format(set_num + 1)
            for num, term in enumerate(space):
                latex_str += '\\pi_{} = '.format(num + 1) + sp.latex(term)
                latex_str += '\\quad'
            display(Math(latex_str))
            display(Markdown('---'))

    def __get_latex_form(self, latex_string=False):
        latex_form = []
        for pi_set in self.__allpiterms:
            latex_set = []
            for pi in pi_set:
                if latex_string:
                    latex_set.append(sp.latex(pi))
                else:
                    latex_set.append(pi)
            latex_form.append(latex_set)

        for num, set in enumerate(latex_form):
            set.insert(0, num + 1)

        return latex_form

    def __tabulate_print(self, latex_string=False):
        latex_sets = self.__get_latex_form(latex_string)
        n = self.num_variable
        m = len(self.__fundamental_vars_used)
        num_of_pi_terms = n - m
        headers = ['sets']
        for num in range(num_of_pi_terms):
            headers.append('Pi {}'.format(num + 1))
        print(tabulate(latex_sets, headers=headers))

    def print_all(self, latex_string=False):
        '''
        print all the sets of dimensionless groups in latex or symbolic form.
        '''
        try:
            self.__Jupyter_print()
        except:
            self.__tabulate_print(latex_string)

    def return_all(self, latex_string=False):
        return self.__get_latex_form(latex_string)
