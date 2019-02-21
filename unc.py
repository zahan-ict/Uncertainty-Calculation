#!/usr/bin/env python 3.6.4
# -*- coding: utf-8 -*-
"""
@developer: Md Sarwar Zahan
Matrix: 01461419
University of Klagenfurt
"""
# !/usr/bin/python 3.6.4
import math
import numpy as np
from info import info_text
from numpy.linalg import inv
from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import inspect


class Unc:
    value = []
    dict_logic = {}
    dict_value = {}
    std_unc = 0
    cat = ''
    name = ''
    dep = ''
    eq = ''

    info_text()
    input("Enter any key to continue ...........")

    # ################## Constructor definition ##############################################
    def __init__(self, value, std_unc=None, cat=None, name=None, dep=None, eq=None):
        self.value = value
        self.std_unc = std_unc
        self.cat = cat
        self.name = name
        self.dep = dep
        self.eq = eq

        # print((inspect.stack()[1].code_context[-1].split()[0]))

        # ************************ Update equation **************************************************
        if self.eq is not None:
            Unc.dict_logic[str(id(self))] = self.eq
            # print(Unc.dict_logic)

        # ************************ ********************************************************************

        if isinstance(self.value, np.matrix) and isinstance(self.std_unc, np.ndarray):
            # convert to array and do list

            if self.value.shape == self.std_unc.shape:
                a = self.value
                e = self.std_unc
                value_b = np.array(a).flatten().tolist()
                std_unc_f = np.array(e).flatten().tolist()

                data_r = getStdUnc(value_b, std_unc_f, self.cat)
                data = np.array(data_r)
                shape = a.shape
                Unc.nameANDdep(self.name, self.dep)
                print(data.reshape(shape))

            else:
                print('Error: Matrix dimension must agree')

    def __repr__(self):
        return "Unc(%s,%s)" % (self.value, self.std_unc)

    # ####################################### Output result format ##################################################
    def __str__(self):

        if self.cat == "L":

            # --------------structured 3 digits standard uncertainty---------------------
            if (self.std_unc == 0) or (np.log10(abs(self.value) / abs(self.std_unc)) > 12):
                digit = -12

            else:
                digit = np.floor(np.log10(abs(self.std_unc)))
            std_unc_digits = abs(round(self.std_unc / math.pow(10, digit - 2)))
            # --------------- Show Only the last 3 values--------------------------------------
            value_digits = round(self.value / math.pow(10, digit - 2)) * math.pow(10, digit - 2)
            if std_unc_digits == 0:
                Unc.nameANDdep(self.name, self.dep)
                return "%f%s" % (value_digits, '(0)')
            else:
                if digit < 0:
                    t = int(abs(digit) + 2)  # number of decimal point on value, must be int value
                    Unc.nameANDdep(self.name, self.dep)
                    return "%.*f(%d)" % (t, value_digits, std_unc_digits)

                else:
                    if digit < 2:

                        std_digits = std_unc_digits * math.pow(10, digit - 2)
                        tc = int(2 - digit)  # number of decimal point on value, must be int value
                        Unc.nameANDdep(self.name, self.dep)
                        return "%.*f(%.*f)" % (tc, value_digits, tc, std_digits)
                    else:

                        std_digits = std_unc_digits * math.pow(10, digit - 2)
                        Unc.nameANDdep(self.name, self.dep)
                        return "%s(%d)" % (value_digits, std_digits)

        elif self.cat == "R":

            # ************************************ format nominal value ******************************************************
            if self.value == 0.0:
                check_value = self.value
            elif self.value < 0:
                check_value = self.value
            elif self.value > 0:
                check_value = math.log10(self.value)  # check if value start with 0.0001 or 1.0001 or -1.0001
            else:
                check_value = self.value

            if check_value > 0:
                value = float("{0:.1f}".format(self.value))  # Print values till 1 decimal point.ex: 12 = 12.0
            else:
                value = float("{0:.4f}".format(self.value))  # Print values till 4 decimal point.ex: 12 = 12.0001
            # ***************************************************************************************************************

            us = format(self.std_unc, '.4f')  # Print  uncertainties till 3 decimal point.ex: 12 = 12.000
            uf = float(us)  # Convert string into float for uncertainty
            std_unc = abs(uf)  # Print absolute uncertainties. remove negative uncertainty number.ex: -12.333 = 12.333
            unc_relative = round(std_unc * 100 / value, 2)
            Unc.nameANDdep(self.name, self.dep)

            # send output
            return '%s(%s%s)' % (value, unc_relative, '%')

        else:

            # ************************************ format nominal value ******************************************************
            if self.value == 0.0:
                check_value = self.value
            elif self.value < 0:
                check_value = self.value
            elif self.value > 0:
                check_value = math.log10(self.value)  # check if value start with 0.0001 or 1.0001 or -1.0001
            else:
                check_value = self.value

            if check_value > 0:
                value = float("{0:.1f}".format(self.value))  # Print values till 1 decimal point.ex: 12 = 12.0
            else:
                value = float("{0:.4f}".format(self.value))  # Print values till 4 decimal point.ex: 12 = 12.0001
            # ***************************************************************************************************************

            # ************************************ format absolute uncertainty ******************************************************
            if self.std_unc == 0.0:
                check_std = self.std_unc
            elif self.std_unc < 0:
                check_std = self.std_unc
            elif self.std_unc > 0:
                check_std = math.log10(self.std_unc)  # check if number start with 0.0001 or 1.0001 or -1.0001
            else:
                check_std = self.std_unc

            if check_std > 0:
                std = float("{0:.2f}".format(abs(self.std_unc)))  # Print uncertainty till 1 decimal point.ex: 12 = 12.0
            else:
                std = float(
                    "{0:.4f}".format(abs(self.std_unc)))  # Print uncertainty till 4 decimal point.ex: 12 = 12.0001
            # ***************************************************************************************************************

            Unc.nameANDdep(self.name, self.dep)

            # send output
            return "%s(%s)" % (value, std)

    # ######################################################################################################################
    # ############################ ARITHMETIC CALCULATION  ##########################

    # Overloading "+" Addition Operator
    def __add__(self, other):

        #  Check the value if it is  matrix, vector or scalar
        if isinstance(self.value, np.matrix) and isinstance(self.std_unc, np.ndarray):
            if self.value.shape == self.std_unc.shape:
                # ========== Get the matrix values and sum it up ==================
                M1 = self.value
                M2 = other.value
                M_sum = M1 + M2
                M_sum = np.array(M_sum).flatten().tolist()

                # ========== Get the uncertainties and sum it up ====================
                U1 = self.std_unc
                U2 = other.std_unc
                U1_1D = np.array(U1).flatten().tolist()  # Convert uncertainties into 1 dimension
                U2_1D = np.array(U2).flatten().tolist()  # Convert uncertainties into 1 dimension

                # Perform Matrix Addition
                U_DATA = []
                for M1_U, M2_U in zip(U1_1D, U2_1D):
                    U_DATA.append(math.sqrt(math.pow(M1_U, 2) + math.pow(M2_U, 2)))

                data_r = getStdUnc(M_sum, U_DATA, self.cat)
                data = np.array(data_r)

                shape = M1.shape
                out = data.reshape(shape)
                return out

            else:
                print("Error: Matrix dimension must agree")

        elif isinstance(self.value, np.ndarray) and isinstance(self.std_unc, np.ndarray):

            # Multi array uncertainty addition calculation
            cal = 0
            for U in self.std_unc:
                a = math.pow(U, 2)  # perform sum calculation
                cal += a
            std_unc = math.sqrt(cal)

            cal2 = 0
            for U2 in other.std_unc:
                a2 = math.pow(U2, 2)  # perform sum calculation
                cal2 += a2
            other_unc = math.sqrt(cal2)

            new_unc = math.sqrt(math.pow(std_unc, 2) + math.pow(other_unc, 2))
            v_sum = np.sum(self.value) + np.sum(other.value)
            return Unc(v_sum, new_unc, self.cat, self.name, self.dep)

        else:

            # Check if there is a constant value
            # print('--add-- can work like z = Unc(2,3)+2')

            self.std_unc = chkConst(self.std_unc)
            other.std_unc = chkConst(other.std_unc)

            # Perform an addition calculation
            # ############################################### start basic loop control for correlation ######################
            opera = '+'
            nominal_value = self.value + other.value

            ans = correlation_check(inspect.currentframe(), self, other, opera)

            return Unc(nominal_value, ans[0], self.cat, self.name, self.dep, ans[1])

            # ############################################### end  basic loop control for correlation ######################


    # def __radd__(self, other):
    # print('Assignment  possible by class object. exe: x=Unc(2,0)')

    # Overloading "-" Subtraction Operator
    def __sub__(self, other):

        #  Check the value if it is scalar or matrix
        if isinstance(self.value, np.matrix) and isinstance(self.std_unc, np.ndarray):
            if self.value.shape == self.std_unc.shape:
                # ========== Get the matrix values and sum it up ==================
                M1 = self.value
                M2 = other.value
                M_sum = M1 - M2
                M_sum = OneDList(M_sum)
                # ========== Get the uncertainties and perform subtraction rule  ====================
                U1 = self.std_unc
                U2 = other.std_unc
                U1_1D = OneDList(U1)  # Convert uncertainties into 1 dimension
                U2_1D = OneDList(U2)  # Convert uncertainties into 1 dimension

                # Perform Matrix subtraction
                U_DATA = []
                for M1_U, M2_U in zip(U1_1D, U2_1D):
                    U_DATA.append(math.sqrt(math.pow(M1_U, 2) + math.pow(M2_U, 2)))

                data_r = getStdUnc(M_sum, U_DATA, self.cat)
                data = np.array(data_r)
                shape = M1.shape
                return data.reshape(shape)
            else:
                print("Error: Matrix dimension must agree")

        elif isinstance(self.value, np.ndarray) and isinstance(self.std_unc, np.ndarray):

            cal = 0
            for U in self.std_unc:
                a = math.pow(U, 2)  # perform sum calculation
                cal += a
            std_unc = math.sqrt(cal)

            cal2 = 0
            for U2 in other.std_unc:
                a2 = math.pow(U2, 2)  # perform sum calculation
                cal2 += a2
            other_unc = math.sqrt(cal2)

            new_unc = math.sqrt(math.pow(std_unc, 2) + math.pow(other_unc, 2))
            v_sub = np.sum(self.value) - np.sum(other.value)
            return Unc(v_sub, new_unc, self.cat, self.name, self.dep)

        else:
            # Check if there is a Constant value
            self.std_unc = chkConst(self.std_unc)
            other.std_unc = chkConst(other.std_unc)

            # ############################################### start basic loop control for correlation ######################
            opera = "-"
            nominal_value = self.value - other.value

            ans = correlation_check(inspect.currentframe(), self, other, opera)

            return Unc(nominal_value, ans[0], self.cat, self.name, self.dep, ans[1])

            # ############################################### end  basic loop control for correlation ######################

    # Overloading "*"  Multiplication Operator
    def __mul__(self, other):
        #  Check the value if it is scalar or matrix
        if isinstance(self.value, np.matrix) and isinstance(self.std_unc, np.ndarray):
            if self.value.shape == self.std_unc.shape:
                # ========== Get the matrix values and start multiplication process ==================
                M1 = self.value
                M2 = other.value
                U1 = self.std_unc
                U2 = other.std_unc

                # Perform value multiplication
                M_sum = M1 * M2
                M_sum = OneDList(M_sum)

                # Perform uncertainties multiplication in matrix form
                MatShap = self.std_unc.shape[1]

                MAT_UNC = []
                for RC in range(MatShap):
                    M1_1D_R = OneDList(M1[RC,])  # Convert  value matrix 1 into 1 dimension
                    U1_1D_R = OneDList(U1[RC,])  # Convert uncertainties of first matrix into 1 dimension

                    for RC2 in range(MatShap):
                        M2_1D_C = OneDList(M2[:, RC2])  # Convert value matrix 2 into 1 dimension
                        U2_1D_C = OneDList(
                            U2[:, RC2])  # Convert uncertainties uncertainties of 2nd matrix into 1 dimension

                        SUM = 0
                        for M1_U, M2_U, M1_V1, M2_V2 in zip(U1_1D_R, U2_1D_C, M1_1D_R, M2_1D_C):
                            SUM_V = M1_V1 * M2_V2
                            gen = SUM_V * (math.sqrt(math.pow(M1_U / M1_V1, 2) + math.pow(M2_U / M2_V2,
                                                                                          2)))  # General Formula of Uncertainty for Multiplication
                            a1 = math.pow(gen, 2)
                            SUM += a1
                        MAT_UNC.append(math.sqrt(SUM))

                # print(MAT_UNC)
                data_r = getStdUnc(M_sum, MAT_UNC, self.cat)
                data = np.array(data_r)
                # print(data_r)
                shape = M1.shape
                return data.reshape(shape)

            else:
                print("Error: Matrix dimension must agree")

        elif isinstance(self.value, np.ndarray) and isinstance(self.std_unc, np.ndarray):
            # Multi array uncertainty addition calculation

            add = 0
            for U in self.std_unc:
                a = math.pow(U, 2)  # perform sum calculation
                add += a
            std_unc = math.sqrt(add)
            v_sub = np.sum(self.value) * other.value

            if self.cat == 'R':
                new_unc = math.sqrt(math.pow(std_unc, 2))
            else:
                new_unc = math.sqrt(math.pow(std_unc, 2)) * other.value

            return Unc(v_sub, new_unc, self.cat, self.name, self.dep)

        else:
            # Check if there is a constant value
            self.std_unc = chkConst(self.std_unc)
            other.std_unc = chkConst(other.std_unc)

            # ############################################### start basic loop control for correlation ######################
            opera = "*"
            nominal_value = self.value * other.value

            ans = correlation_check(inspect.currentframe(), self, other, opera)

            return Unc(nominal_value, ans[0], self.cat, self.name, self.dep, ans[1])

            # ############################################### end  basic loop control for correlation ######################

    # Overloading "/" Division Operator
    def __truediv__(self, other):
        # Check  vector or scalar values and uncertainties
        if isinstance(self.value, np.ndarray) and isinstance(self.std_unc, np.ndarray):

            con = 0
            for U in self.std_unc:
                a = math.pow(U, 2)  # perform sum calculation
                con += a
            std_unc = math.sqrt(con)
            v_sub = np.sum(self.value) / other.value

            if self.cat == 'R':
                new_unc = math.sqrt(math.pow(std_unc, 2))
            else:
                new_unc = math.sqrt(math.pow(std_unc, 2)) / other.value

            return Unc(v_sub, new_unc, self.cat, self.name, self.dep)

        else:
            # Check if there is a constant value
            self.std_unc = chkConst(self.std_unc)
            other.std_unc = chkConst(other.std_unc)

            # ############################################### start basic loop control for correlation ######################
            opera = "/"
            nominal_value = self.value / other.value

            ans = correlation_check(inspect.currentframe(), self, other, opera)

            return Unc(nominal_value, ans[0], self.cat, self.name, self.dep, ans[1])

            # ############################################### end  basic loop control for correlation ######################

    # Overloading "power" as a Polynomial functions. Formula if R = X^n    delta(R) = (|n|.delta(X).|R|)/X
    def __pow__(self, other):

        frame = inspect.currentframe()
        raw_d = frame.f_back.f_locals
        # key_list = list(raw_d.keys())
        # print(key_list)

        # key_value = list(raw_d.values())
        # print(key_value)

        tmp_equ = []

        # ***************************** if self is other ****************************
        # search dict_logic dictionary if any correlation exits as self id
        if str(id(self)) in Unc.dict_logic:
            tmp_equ.append(Unc.dict_logic[str(id(self))])

        for q_src in raw_d:
            if raw_d[q_src] == self:
                tmp_equ.append(q_src)

                # build global dict_value
                Unc.dict_value[q_src] = self.value
                Unc.dict_value[q_src + '_u'] = self.std_unc
                Unc.dict_value[q_src + '_id'] = id(self)

        # print(tmp_equ)
        sym = symbols(tmp_equ[0])
        sym_eq = '(%s)' % (sym)

        # build equation
        equation = parse_expr(sym_eq)**other
        eq = str(equation)
        # ***************** power calculation ******************************

        nominal_value = self.value ** other
        std_unc = (other * self.std_unc * self.value ** other) / self.value

        return Unc(nominal_value, std_unc, self.cat, self.name, self.dep, eq)

    # Calculation of  Square root on given function
    def sqrt(self):
        value = self.value ** 0.5
        std_unc = (0.5 * self.std_unc * value) / self.value
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # Calculation of natural logarithm - ln(x) on given function
    def ln(self):
        value = math.log(self.value)  # ln(Value)
        std_unc = self.std_unc / self.value
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # Calculation of logarithm - log10 on given function
    def ulog(self):
        value = math.log10(self.value)
        std_unc = 0.434 * (self.std_unc / self.value)
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # Calculation of Antilog 10^x
    def tenPower(self):
        value = 10 ** self.value
        ln10 = 2.3026
        std_unc = value * ln10 * self.std_unc  # ln10=2.3026
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # Exponential function (e^x) calculation
    def uexp(self):
        value = math.exp(self.value)  # e^1=2.718
        std_unc = value * self.std_unc
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # ############################  STATISTICAL CALCULATION  ###########################
    # ======== !!Calculations are performed in radian !! =================
    def sum(self):
        V_sum = np.sum(self.value)
        cal = 0
        for U in self.std_unc:
            a = math.pow(U, 2)  # perform sum calculation
            cal += a
        std_unc = math.sqrt(cal)
        return Unc(V_sum, std_unc, self.cat, self.name, self.dep)

    def mean(vs):
        mean = np.mean(vs)
        return mean

    def stdev(vs):
        std = np.std(vs)
        return std

    def statUnc(self):
        value = self.value
        count = np.count_nonzero(value)
        std = Unc.stdev(value)
        mean = Unc.mean(value)
        meanAvg = std / math.sqrt(count)
        return Unc(mean, meanAvg, self.cat, self.name, self.dep)

    def nameANDdep(name, dep):

        if isinstance(name, str) and isinstance(dep, str):
            print("Category is:", name, "  & Department is: ", dep)
        elif isinstance(name, str):
            print("Category is:", name)
        elif isinstance(dep, str):
            print("Department is:", dep)
        else:
            return None

    # ############################ TRIGONOMETRIC FUNCTIONS CALCULATION  ###########################
    # ======== !!Calculations are performed in radian !! =================
    ####################################################################################

    # sinus function calculation.
    def sin(self):
        value = math.sin(self.value)
        std_unc = self.std_unc * math.cos(self.value)  # if y = sin(x) than U(y) = U(x)cos(x)
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # cosine function calculation
    def cos(self):
        value = math.cos(self.value)
        std_unc = - self.std_unc * math.sin(self.value)  # if y = cos(x) than U(y) = -U(x)sin(x)
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # tan function calculation
    def tan(self):
        value = math.tan(self.value)
        secSquared = (2 / (math.cos(2 * self.value)) + 1)
        std_unc = self.std_unc * secSquared  # if y = tan^2(x) than U(y) = -U(x)sec^2(x)
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # cot function calculation
    def cot(self):
        value = 1 / math.tan(self.value)
        csecSquared = -(2 / (1 - math.cos(2 * self.value)))
        std_unc = self.std_unc * csecSquared  # if y = cot^2(x) than U(y) = -U(x) csc^2(x)
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # ########## Inverse Trigonometric Calculation ###########################################

    # arcsin function calculation
    def arcsin(self):
        value = math.asin(self.value)
        dx = (1 / math.sqrt(1 - self.value ** 2))
        std_unc = self.std_unc * dx  # if y = sin^-1(x) than U(y) = -U(x) 1/sqrt(1-x^2)
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # arcos function calculation
    def arccos(self):
        value = math.acos(self.value)
        dx = math.sqrt(1 - self.value ** 2)
        dxr = -1 / dx
        std_unc = self.std_unc * dxr  # if y = cos^-1(x) than U(y) = -U(x) -1/sqrt(1-x^2)
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # arctan function calculation
    def arctan(self):
        value = math.atan(self.value)
        dx = 1 + self.value ** 2
        dxr = 1 / dx
        std_unc = self.std_unc * dxr  # if y = tan^-1(x) than U(y) = -U(x) 1/1+x^2
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # ########## Hyperbolic Trigonometric Calculation ###########################################

    # sinhx function calculation
    def sinh(self):
        value = math.sinh(self.value)
        dxr = math.cosh(self.value)
        std_unc = self.std_unc * dxr  # if y = sinhx than U(y) = U(x)coshx
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # coshx function calculation
    def cosh(self):
        value = math.cosh(self.value)
        dxr = math.sinh(self.value)
        std_unc = self.std_unc * dxr  # if y = coshx than U(y) = U(x)sinhx
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # tanhx function calculation
    def tanh(self):
        value = math.tanh(self.value)
        dx1 = 1 - math.cosh(2 * self.value)
        dx2 = 1 + math.cosh(2 * self.value)
        dx3 = dx1 * dx2
        dxr = dx3 / 4
        dxrf = (1 - dxr)
        std_unc = self.std_unc * dxrf  # if y = tanhx than U(y) = U(x)(1-tanh^2x)
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # ########## Inverse Hyperbolic Trigonometric Calculation ###########################################

    # asinhx function calculation

    def arcsinh(self):
        value = math.asinh(self.value)
        dx1 = math.sqrt(self.value ** 2) + math.sqrt(1)
        dxr = 1 / dx1
        std_unc = self.std_unc * dxr  # if y = asinh(x) than U(y) = U(x) 1/sqrt(x^2+1)
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # acoshx function calculation
    def arccosh(self):
        value = math.acosh(self.value)
        dx1 = math.sqrt(self.value ** 2) - math.sqrt(1)
        dxr = 1 / dx1
        std_unc = self.std_unc * dxr  # if y = acosh(x) than U(y) = U(x) 1/sqrt(x^2-1)
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    # atanhx function calculation
    def arctanh(self):
        value = math.atanh(self.value)
        dx1 = 1 - self.value ** 2
        dxr = 1 / dx1
        std_unc = self.std_unc * dxr  # if y = atanh(x) than U(y) = U(x) 1/1-x^2
        return Unc(value, std_unc, self.cat, self.name, self.dep)

    ####################################################################################
    # ======== !!Calculations are performed in degree !! =================
    ####################################################################################

    def degree(self, str):

        if str == "sin":
            # ========= Example of Sin values in degree ==============================
            # np.sin(np.array((0., 30., 45., 60., 90.)) * np.pi / 180. )
            # array([ 0.        ,  0.5       ,  0.70710678,  0.8660254 ,  1. ])

            # ================ sin(x) is calculated in degree  ==================
            value = np.sin(self.value * np.pi / 180)
            std_unc = (self.std_unc * np.pi / 180) * np.cos(self.value * np.pi / 180)
            return "%s(%s)" % (value, std_unc)

        elif str == "cos":
            # ================ cos(x) is calculated in degree ==================
            value = np.cos(self.value * np.pi / 180)
            std_unc = - (self.std_unc * np.pi / 180) * np.sin(self.value * np.pi / 180)
            return "%s(%s)" % (value, std_unc)

        elif str == "tan":
            # ================ tan(x) is calculated in degree ==================
            valueN = self.value * np.pi / 180
            std_uncN = self.std_unc * np.pi / 180
            value = np.tan(self.value * np.pi / 180)
            secSquared = (2 / (math.cos(2 * valueN)) + 1)
            std_unc = std_uncN * secSquared  # if y = tan^2(x) than U(y) = -U(x)sec^2(x)
            return "%s(%s)" % (value, std_unc)

        elif str == "cot":
            # ================ cot(x) is calculated in degree ==================
            valueN = self.value * np.pi / 180
            std_uncN = self.std_unc * np.pi / 180
            value = 1 / np.tan(self.value * np.pi / 180)
            csecSquared = -(2 / (1 - math.cos(2 * valueN)))
            std_unc = std_uncN * csecSquared  # if y = cot^2(x) than U(y) = -U(x) csc^2(x)
            return "%s(%s)" % (value, std_unc)
        else:
            print(" Format is not correct")

    def inv_matrix(self, value):
        print('\n Inverse Matrix is:')
        invsMatrix = inv(self.value)
        print(invsMatrix)


####################################################################################
# ======== Correlation handler =================
####################################################################################


def correlation_check(frame, self_data, other_data, opera):

    # print("-----------------correlation_covariance function start---------------------------")

    # print('self id:', id(self_data), 'other id:', id(other_data))
    # print('self data:', self_data, 'other data:', other_data)

    raw_d = frame.f_back.f_locals
    # key_list = list(raw_d.keys())
    # print(key_list)

    # key_value = list(raw_d.values())
    # print(key_value)

    tmp_equ = []

    # ***************** start inspecting frame and dict_logic dictionary. *****************************************

    # search dict_logic dictionary if any correlation exits as self id
    if str(id(self_data)) in Unc.dict_logic:
        tmp_equ.append(Unc.dict_logic[str(id(self_data))])

    else:
        # two same loop keep variable sequence correct
        for q_src in raw_d:
            if raw_d[q_src] == self_data:
                tmp_equ.append(q_src)

                # build global dict_value
                Unc.dict_value[q_src] = self_data.value
                Unc.dict_value[q_src + '_u'] = self_data.std_unc
                Unc.dict_value[q_src + '_id'] = id(self_data)
                break

    # search dict_logic dictionary if any correlation exits as other id
    if str(id(other_data)) in Unc.dict_logic:
        tmp_equ.append(Unc.dict_logic[str(id(other_data))])
    else:
        for q_src in raw_d:

            if raw_d[q_src] == other_data:
                tmp_equ.append(q_src)
                # build global dict_value
                Unc.dict_value[q_src] = other_data.value
                Unc.dict_value[q_src + '_u'] = other_data.std_unc
                Unc.dict_value[q_src + '_id'] = id(other_data)
                break

    # *************************************************************************************
    # ************* correct equation building *********************************************

    tmp = []
    for lm in tmp_equ:
        element_in_list = '(%s)' % (lm)
        tmp.append(element_in_list)
    equation = opera.join(tmp)  # build equation as string

    # build equation
    equ = parse_expr(equation)

    # *************************************************************************************
    # ***************** start partial derivative calculation ******************************

    # step1: symbol extraction form equation
    sym_list = str(srepr(equ)).split(",")
    sym = []

    for sign in sym_list:

        if 'Symbol' in sign:
            sz = sign[sign.find("'") + 1:]
            sym.append(sz[:sz.find(")") - 1])

    # step:2 remove duplicate from symbol list
    sym = list(set(sym))

    # step: 3 perform partial derivative
    cal_r = 0
    for symbol in sym:
        derivative = diff(equ, symbol)
        sym_list2 = str(srepr(derivative)).split(",")

        # step: 4 start numerical calculation

        sym2 = []
        for sign2 in sym_list2:

            if 'Symbol' in sign2:
                sz2 = sign2[sign2.find("'") + 1:]
                v_value = Unc.dict_value[sz2[:sz2.find(")") - 1]]
                list_r = sz2[:sz2.find(")") - 1], v_value
                sym2.append(list_r)

        # step : 5 apply uncertainty correlation formula

        pde = derivative.subs(sym2)
        std = Unc.dict_value[symbol + '_u']
        std_sqr = np.power(std, 2)
        pde_sqr = np.power(pde, 2)
        sol = pde_sqr * std_sqr
        cal_r += sol
    final_std_unc = math.sqrt(cal_r)

    return final_std_unc, equation

# ********************************************************************************************

# ########################## Convert  matrix or array into single dimension ##################

def OneDList(params):
    return np.array(params).flatten().tolist()


# ###################################################################################
# ======== !!Check the constant value  !! =================
# ####################################################################################
def chkConst(std_check):
    if std_check is None:
        return 0
    else:
        return std_check


# ########################### remove cot########################################################
def rmParen(comp):
    comp = str(comp)
    comp = comp.strip(')')
    comp = comp.strip('(')
    return comp


# ###################################################################################
# ======== !!Matrix or array Calculation  !! =================
# ####################################################################################

def getStdUnc(value_b, std_unc_f, cat):
    if cat == 'L':

        data = []
        for value_item, std_item in zip(value_b, std_unc_f):
            # --------------structured 3 digits standard uncertainty---------------------
            if (std_item == 0) or (np.log10(abs(value_item) / abs(std_item)) > 12):
                digit = -12
            else:
                digit = np.floor(np.log10(abs(std_item)))
            std_unc_digits = round(std_item / math.pow(10, digit - 2))
            # --------------- Show Only the last 3 values--------------------------------------
            value_digits = round(value_item / math.pow(10, digit - 2)) * math.pow(10, digit - 2)

            if std_unc_digits == 0:
                return "%f%s" % (value_digits, '(0)')
            else:
                if digit < 0:
                    t = int(abs(digit) + 2)  # number of decimal point on value, must be int value

                    data.append("%.*f(%d)" % (t, value_digits, std_unc_digits))
                else:
                    if digit < 2:

                        std_digits = std_unc_digits * math.pow(10, digit - 2)
                        tc = int(2 - digit)  # number of decimal point on value, must be int value
                        data.append('%.*f(%.*f)' % (tc, value_digits, tc, std_digits))

                    else:
                        std_digits = std_unc_digits * math.pow(10, digit - 2)

                        data.append('%s(%d)' % (value_digits, std_digits))
        return data

    elif cat == 'R':

        data = []
        for value_item, std_item in zip(value_b, std_unc_f):
            vs = format(value_item, '.4f')  # Print values till 3 decimal point.ex: 12 = 12.000
            vf = float(vs)  # Convert string into float for values
            value = abs(vf)  # Print absolute values. remove negative numbers.ex:-12.333 = 12.333

            us = format(std_item, '.3f')  # Print  uncertainties till 3 decimal point.ex: 12 = 12.000
            uf = float(us)  # Convert string into float for uncertainty
            std_unc = abs(uf)  # Print absolute uncertainties. remove negative uncertainty number.ex: -12.333 = 12.333
            unc_relative = round(std_unc * 100 / value, 2)
            data.append("%s(%s%s)" % (value, unc_relative, '%'))

        return data

    else:
        data = []
        for value_item, std_item in zip(value_b, std_unc_f):
            vs = format(value_item, '.4f')  # Print values till 3 decimal point.ex: 12 = 12.000
            vf = float(vs)  # Convert string into float for values
            value = abs(vf)  # Print absolute values. remove negative numbers.ex:-12.333 = 12.333

            us = format(std_item, '.3f')  # Print  uncertainties till 3 decimal point.ex: 12 = 12.000
            uf = float(us)  # Convert string into float for uncertainty
            std_unc = abs(uf)  # Print absolute uncertainties. remove negative uncertainty number.ex: -12.333 = 12.333

            data.append("%s(%s)" % (value, std_unc))

        return data
