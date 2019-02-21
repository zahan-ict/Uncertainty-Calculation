#!/usr/bin/env python 3.6.4
# -*- coding: utf-8 -*-
"""
Complex Uncertainty Calculation
@developer: Md Sarwar Zahan
Matrix: 01461419
University of Klagenfurt
"""
import inspect
from unc import Unc
import numpy as np
import math
import cmath
import sympy


class Ucom(Unc):
    comp_value = 0
    id_check = [0]
    comp_unc = 0
    name = ''
    dep = ''

    # Constructor definition ##############################################
    def __init__(self, comp_value, comp_unc=None, name=None, dep=None):
        self.comp_value = comp_value
        self.comp_unc = comp_unc
        self.name = name
        self.dep = dep

    # Default string function for formatting output ###############################
    def __str__(self):
        comp_value = rmParen(self.comp_value)
        comp_unc = rmParen(self.comp_unc)
        Unc.nameANDdep(self.name, self.dep)

        return "%s[%s]" % (comp_value, comp_unc)

    #  Overloading "+" Addition Operator For Complex Calculation
    def __add__(self, other):

        #  Check the value if it is scalar or vector
        if isinstance(self.comp_value, np.ndarray) and isinstance(self.comp_unc, np.ndarray):
            if self.comp_value.shape[0] == self.comp_unc.shape[0]:

                C_V1 = self.comp_value
                C_U1 = self.comp_unc
                C_V2 = other.comp_value
                C_U2 = other.comp_unc

                M1 = 0
                for U in C_U1:
                    a = np.power(U, 2)
                    M1 += a
                comp_unc = np.sqrt(M1)

                M2 = 0
                for U2 in C_U2:
                    a2 = np.power(U2, 2)
                    M2 += a2
                other_comp_unc = np.sqrt(M2)

                new_comp_unc = np.sqrt(np.power(comp_unc, 2) + np.power(other_comp_unc, 2))
                v_sum = np.sum(C_V1) + np.sum(C_V2) # Complex Uncertainty calculation  for addition

                return Ucom(v_sum, new_comp_unc, self.name, self.dep)

        else:
            # Perform complex addition calculation
            comp_value = np.add(self.comp_value,other.comp_value)

            # check variable correlation on complex addition

            if self is other:
                comp_unc = np.subtract(self.comp_unc, other.comp_unc)
                Ucom.id_check = id(self), id(other)
            elif Ucom.id_check[0] == id(other):
                comp_unc = np.subtract(self.comp_unc, other.comp_unc)
            else:
                comp_unc = np.sqrt(np.power(self.comp_unc, 2) + np.power(other.comp_unc, 2))  # General Formula of Uncertainty for Subtraction

            return Ucom(comp_value, comp_unc, self.name, self.dep)

    # Overloading "-" Subtraction Operator For Complex Calculation
    def __sub__(self, other):

        #  Check the value if it is scalar or vector
        if isinstance(self.comp_value, np.ndarray) and isinstance(self.comp_unc, np.ndarray):
            if self.comp_value.shape[0] == self.comp_unc.shape[0]:

                C_V1 = self.comp_value
                C_U1 = self.comp_unc
                C_V2 = other.comp_value
                C_U2 = other.comp_unc

                M1 = 0
                for U in C_U1:
                    a = np.power(U, 2)
                    M1 += a
                comp_unc = np.sqrt(M1)

                v1 = 0
                for U2 in C_U2:
                    a2 = np.power(U2, 2)
                    v1 += a2
                other_comp_unc = np.sqrt(v1)

                new_comp_unc = np.sqrt(np.power(comp_unc, 2) + np.power(other_comp_unc, 2))
                v_sub = np.sum(C_V1) - np.sum(C_V2) # Complex Uncertainty calculation  for subtraction

                return Ucom(v_sub, new_comp_unc, self.name, self.dep)

        else:
            # Perform complex subtraction calculation
            comp_value = np.subtract(self.comp_value, other.comp_value)

            # check variable correlation on complex subtraction
            if self is other:
                comp_unc = np.subtract(self.comp_unc, other.comp_unc)
                Ucom.id_check= id(self), id(other)
            elif Ucom.id_check[0] == id(other):
                comp_unc = np.subtract(self.comp_unc, other.comp_unc)
            else:
                comp_unc = np.sqrt(np.power(self.comp_unc, 2) + np.power(other.comp_unc,2))  # General Formula of Uncertainty for Subtraction

            return Ucom(comp_value, comp_unc, self.name, self.dep)

    # Overloading "*" multiplication Operator For Complex Calculation
    def __mul__(self, other):

        #  Check the value if it is scalar or vector
        if isinstance(self.comp_value, np.ndarray) and isinstance(self.comp_unc, np.ndarray):
            if self.comp_value.shape[0] == self.comp_unc.shape[0]:
                C_V1 = self.comp_value
                C_U1 = self.comp_unc
                C_V2 = other.comp_value
                C_U2 = other.comp_unc
                # Perform value multiplication
                mul = 0
                for U in C_U1:
                    a = np.power(U, 2)  # do sum
                    mul += a
                comp_unc = np.sqrt(mul)
                comp_value = np.sum(C_V1) * C_V2

                return Ucom(comp_value, comp_unc, self.name, self.dep)
            else:
                print("Error: Number of element of values and uncertainties must be same")

        else:

            # Perform complex multiplication calculation
            comp_value = np.multiply(self.comp_value, other.comp_value)

            # check variable correlation on complex division
            if self is other:
                comp_unc = np.divide(self.comp_unc, other.comp_unc)

                Ucom.id_check = id(self), id(other)

            elif Ucom.id_check[0] == id(other):

                comp_unc = np.divide(self.comp_unc, other.comp_unc)
            else:
                comp_unc = comp_value * (np.sqrt(np.power(self.comp_unc / self.comp_value, 2) + np.power(other.comp_unc / other.comp_value, 2)))  # General Formula of Uncertainty for multiplication

            return Ucom(comp_value, comp_unc, self.name, self.dep)

    # Overloading "/" Division Operator For Complex Calculation
    def __truediv__(self, other):

        # Check the value if it is scalar or vector
        if isinstance(self.comp_value, np.ndarray) and isinstance(self.comp_unc, np.ndarray):
            if self.comp_value.shape[0] == self.comp_unc.shape[0]:
                C_V1 = self.comp_value
                C_U1 = self.comp_unc
                C_V2 = other.comp_value
                C_U2 = other.comp_unc

                # Perform complex division calculation
                div = 0
                for U in C_U1:
                    a = np.power(U, 2)  # do sum
                    div += a
                comp_unc = np.sqrt(div)
                comp_value = np.sum(C_V1) / C_V2

                return Ucom(comp_value, comp_unc, self.name, self.dep)

            else:
                print("Error: Number of element of values and uncertainties must be same")

        else:
            # Perform complex division calculation

            comp_value = np.divide(self.comp_value, other.comp_value)

            # check variable correlation on complex division
            if self is other:

                comp_unc = np.divide(self.comp_unc, other.comp_unc)
                Ucom.id_check = id(self), id(other)
            elif Ucom.id_check[0] == id(other):

                comp_unc = np.divide(self.comp_unc, other.comp_unc)
            else:
                comp_unc = comp_value * (np.sqrt(np.power(self.comp_unc / self.comp_value, 2) + np.power(other.comp_unc / other.comp_value, 2)))  # General Formula of Uncertainty for Division

            return Ucom(comp_value, comp_unc, self.name, self.dep)

    # Overloading "power" as a Polynomial functions for complex. Formula if R = X^n    delta(R) = (|n|.delta(X).|R|)/X
    def __pow__(self, other):
        # Perform an power calculation
        comp_value = self.comp_value ** other
        comp_unc = (other * self.comp_unc * self.comp_value ** other) / self.comp_value
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # Calculation of  Square root on given function
    def sqrt(self):
        comp_value = self.comp_value ** 0.5
        comp_unc = (0.5 * self.comp_unc * comp_value) / self.comp_value
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # Calculation of natural logarithm - ln(x) on given function
    def ln(self):
        comp_value = math.log(self.comp_value)  # ln(Value)
        comp_unc = self.comp_unc / self.comp_value
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # Calculation of logarithm - log10 on given function
    def ulog(self):
        comp_value = math.log10(self.comp_value)
        comp_unc = 0.434 * (self.comp_unc / self.comp_value)
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # Calculation of Antilog 10^x
    def tenPower(self):
        comp_value = 10 ** self.comp_value
        ln10 = 2.3026
        comp_unc = comp_value * ln10 * self.comp_unc  # ln10=2.3026
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # Exponential function (e^x) calculation
    def uexp(self):
        comp_value = cmath.exp(self.comp_value)  # e^1=2.718
        comp_unc = comp_value * self.comp_unc
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # ############################ COMPLEX TRIGONOMETRIC FUNCTIONS CALCULATION  ###########################
    # ======== !!Calculations are performed in radian !! =================
    ####################################################################################
    # sinus function calculation.
    def sin(self):
        comp_value = cmath.sin(self.comp_value)
        comp_unc = self.comp_unc * cmath.cos(self.comp_value)  # if y = sin(x) than U(y) = U(x)cos(x)

        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # cosine function calculation
    def cos(self):
        comp_value = cmath.cos(self.comp_value)
        comp_unc = self.comp_unc * cmath.sin(self.comp_value) # if y = sin(x) than U(y) = U(x)cos(x)

        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # tan function calculation
    def tan(self):
        comp_value = cmath.tan(self.comp_value)
        secSquared = (2 / (cmath.cos(2 * self.comp_value)) + 1)
        comp_unc = self.comp_unc * secSquared  # if y = tan^2(x) than U(y) = -U(x)sec^2(x)
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # cot function calculation
    def cot(self):
        comp_value = 1 / cmath.tan(self.comp_value)
        csecSquared = -(2 / (1 - cmath.cos(2 * self.comp_value)))
        comp_unc = self.comp_unc * csecSquared  # if y = cot^2(x) than U(y) = -U(x) csc^2(x)
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # ########## Inverse Trigonometric Calculation (Complex) ###########################################
    # arcsin function calculation
    def arcsin(self):
        comp_value = cmath.asin(self.comp_value)
        dx = (1 / cmath.sqrt(1 - self.comp_value ** 2))
        comp_unc = self.comp_unc * dx  # if y = sin^-1(x) than U(y) = -U(x) 1/sqrt(1-x^2)
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # arcos function calculation
    def arccos(self):
        comp_value = cmath.acos(self.comp_value)
        dx = cmath.sqrt(1 - self.comp_value ** 2)
        dxr = -1 / dx
        comp_unc = self.comp_unc * dxr  # if y = cos^-1(x) than U(y) = -U(x) -1/sqrt(1-x^2)
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # arctan function calculation
    def arctan(self):
        comp_value = cmath.atan(self.comp_value)
        dx = 1 + self.comp_value ** 2
        dxr = 1 / dx
        comp_unc = self.comp_unc * dxr  # if y = tan^-1(x) than U(y) = -U(x) 1/1+x^2
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # ########## Hyperbolic Trigonometric Calculation ###########################################
    # sinhx function calculation
    def sinh(self):
        comp_value = cmath.sinh(self.comp_value)
        dxr = cmath.cosh(self.comp_value)
        comp_unc = self.comp_unc * dxr  # if y = sinhx than U(y) = U(x)coshx
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # coshx function calculation
    def cosh(self):
        comp_value = cmath.cosh(self.comp_value)
        dxr = cmath.sinh(self.comp_value)
        comp_unc = self.comp_unc * dxr  # if y = coshx than U(y) = U(x)sinhx
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # tanhx function calculation
    def tanh(self):
        comp_value = cmath.tanh(self.comp_value)
        dx1 = 1 - cmath.cosh(2 * self.comp_value)
        dx2 = 1 + cmath.cosh(2 * self.comp_value)
        dx3 = dx1 * dx2
        dxr = dx3 / 4
        dxrf = (1 - dxr)
        comp_unc = self.comp_unc * dxrf  # if y = tanhx than U(y) = U(x)(1-tanh^2x)
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # ########## Complex Inverse Hyperbolic Trigonometric Calculation ###########################################
    # asinhx function calculation
    def arcsinh(self):
        comp_value = cmath.asinh(self.comp_value)
        dx1 = cmath.sqrt(self.comp_value ** 2) + cmath.sqrt(1)
        dxr = 1 / dx1
        comp_unc = self.comp_unc * dxr  # if y = asinh(x) than U(y) = U(x) 1/sqrt(x^2+1)
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # acoshx function calculation
    def arccosh(self):
        comp_value = cmath.acosh(self.comp_value)
        dx1 = cmath.sqrt(self.comp_value ** 2) - cmath.sqrt(1)
        dxr = 1 / dx1
        comp_unc = self.comp_unc * dxr  # if y = acosh(x) than U(y) = U(x) 1/sqrt(x^2-1)
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # atanhx function calculation
    def arctanh(self):
        comp_value = cmath.atanh(self.comp_value)
        dx1 = 1 - self.comp_value ** 2
        dxr = 1 / dx1
        comp_unc = self.comp_unc * dxr  # if y = atanh(x) than U(y) = U(x) 1/1-x^2
        return Ucom(comp_value, comp_unc, self.name, self.dep)

    # Complex sum calculation
    def Csum(self):
        CV_sum = np.sum(self.comp_value)
        cal = 0
        for U in self.comp_unc:
            a = np.power(U, 2)  # perform sum calculation
            cal += a
        comp_unc = np.sqrt(cal)
        return Ucom(CV_sum, comp_unc, self.name, self.dep)


    def mean(vs):
        mean = np.mean(vs)
        return mean

    def stdev(vs):
        std = np.std(vs)
        return std

    def comStatisticUnc(self):
        value = self.comp_value
        count = np.count_nonzero(value)
        std = Unc.stdev(value)
        mean = Unc.mean(value)
        meanAvg = std/math.sqrt(count)
        return Ucom(mean, meanAvg, self.name, self.dep)


# Complex numbers without enclosing	parentheses. (6+4j) = 6+4j ==========
def rmParen(comp):
    comp = str(comp)
    comp = comp.strip(')')
    comp = comp.strip('(')
    return comp


# ########## format complex number before print ###########################################

def comformat(com_num):

    if com_num.real > 0 or com_num.real < 0:
        if com_num.real > 0:

            if checktype(com_num.real) == 'int':

                ncom_num = "%.1f%s%.1f%s" % (com_num.real, '+', com_num.imag, 'j')
                ncom_num= ncom_num.replace(" ", "")
                ncom_num = complex(ncom_num)
            else:
                ncom_num = "%.3f%s%.3f%s" % (com_num.real, '+', com_num.imag, 'j')
                ncom_num= ncom_num.replace(" ", "")
                ncom_num = complex(ncom_num)

        elif com_num.real < 0:

            if checktype(com_num.real) == 'int':

                ncom_num = "%.1f%.1f%s" % (com_num.real, com_num.imag, 'j')
                ncom_num= ncom_num.replace(" ", "")
                ncom_num = complex(ncom_num)
            else:
                ncom_num = "%.3f%s%.3f%s" % (com_num.real, '+', com_num.imag, 'j')
                ncom_num= ncom_num.replace(" ", "")
                ncom_num = complex(ncom_num)

        return ncom_num
    else:
        vs = "%.3f%s" % (com_num.imag, 'j')
        vs = vs.replace(" ", "")
        vs = complex(vs)
        return vs


def checktype(num):
    num = str(num-int(num))[1:]
    count = len(num)
    if count == 2:
        nv = num[1]
        if nv == '0':
            num = "int"

    return num





