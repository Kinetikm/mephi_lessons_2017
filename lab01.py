# coding=utf-8
import math
import unittest


def bisection(a, b, delta, eps, f):
    while (b-a)/2. >= eps:
        x_one = (a + b - delta)/2.
        x_two = (a + b + delta)/2.
        if f(x_one) > f(x_two):
            a = x_one
        if f(x_two) > f(x_one):
            b = x_two
    return (a + b)/2.


def golden_ratio(a, b, eps, f):
    phi = (1+math.sqrt(5))/2.
    while (b-a)/2. >= eps:
        x_one = b - (b - a) / phi
        x_two = a + (b - a) / phi
        if f(x_one) >= f(x_two):
            a = x_one
        if f(x_two) > f(x_one):
            b = x_two
    return (a + b)/2.


class Tester(unittest.TestCase):

    def test_bisection(self):
        self.assertTrue(abs(bisection(-1, 2, 0.05, 0.05, lambda x: x**2)-0) < 0.05)
        self.assertTrue(abs(bisection(-7, 1, 0.05, 0.05, lambda x: abs(x))-0) < 0.05)
        self.assertTrue(abs(bisection(0.5, 1.5, 0.05, 0.05, lambda x: ((x-1)**2)*((x-3)**2)) - 1) < 0.05)
        self.assertTrue(abs(bisection(2.4, 5, 0.05, 0.05, lambda x: 2*x**3 - 15*x**2+ 36*x - 14) - 3) < 0.05)
        self.assertTrue(abs(bisection(0.5, 6, 0.05, 0.05, lambda x: 16/x + x + 3) - 4) < 0.05)
        self.assertTrue(abs(bisection(1.1, 5, 0.05, 0.05, lambda x: x**3 / (x**2 - 1)) - math.sqrt(3)) < 0.05)
        self.assertTrue(abs(bisection(0, 15, 0.05, 0.05, lambda x: x**2 -34*x + 140*math.log(x)-10) - 10) < 0.05)
        self.assertTrue(abs(bisection(-20, 17, 0.05, 0.05, lambda x: (x**4)/4 - 2 * x**2 + 1) + 2) < 0.05)
        self.assertTrue(abs(bisection(0, 5, 0.05, 0.05, lambda x: math.cos(x)) - math.pi) < 0.05)
        self.assertTrue(abs(bisection(-7, 5, 0.05, 0.05, lambda x: (x-1)*math.exp(3*x)) - 2/3.) < 0.05)

    def test_golden(self):
        self.assertTrue(abs(golden_ratio(-1, 2, 0.05, lambda x: x ** 2) - 0) < 0.15)
        self.assertTrue(abs(golden_ratio(-7, 1, 0.05, lambda x: abs(x)) - 0) < 0.15)
        self.assertTrue(abs(golden_ratio(0.5, 1.5, 0.05, lambda x: ((x - 1) ** 2) * ((x - 3) ** 2)) - 1) < 0.15)
        self.assertTrue(abs(golden_ratio(2.4, 5, 0.05, lambda x: 2 * x ** 3 - 15 * x ** 2 + 36 * x - 14) - 3) < 0.15)
        self.assertTrue(abs(golden_ratio(0.5, 6, 0.05, lambda x: 16 / x + x + 3) - 4) < 0.15)
        self.assertTrue(abs(golden_ratio(1.1, 5, 0.05, lambda x: x ** 3 / (x ** 2 - 1)) - math.sqrt(3)) < 0.15)
        self.assertTrue(abs(golden_ratio(8, 15, 0.05, lambda x: x ** 2 - 34 * x + 140 * math.log(x) - 10) - 10) < 0.15)
        self.assertTrue(abs(golden_ratio(-20, 17, 0.05, lambda x: (x ** 4) / 4 - 2 * x ** 2 + 1) + 2) < 0.15)
        self.assertTrue(abs(golden_ratio(0, 5, 0.05, lambda x: math.cos(x)) - math.pi) < 0.15)
        self.assertTrue(abs(golden_ratio(-7, 5, 0.05, lambda x: (x - 1) * math.exp(3 * x)) - 2 / 3.) < 0.15)

    def test_compare(self):
        self.assertTrue(abs(bisection(-1, 2, 0.05, 0.05, lambda x: x ** 2) -
                            golden_ratio(-1, 2, 0.05, lambda x: x ** 2)) < 0.3)
        self.assertTrue(abs(bisection(-7, 1, 0.05, 0.05, lambda x: abs(x)) -
                            golden_ratio(-7, 1, 0.05, lambda x: abs(x))) < 0.3)
        self.assertTrue(abs(bisection(0.5, 1.5, 0.05, 0.05, lambda x: ((x - 1) ** 2) * ((x - 3) ** 2)) -
                            golden_ratio(0.5, 1.5, 0.05, lambda x: ((x - 1) ** 2) * ((x - 3) ** 2))) < 0.3)
        self.assertTrue(abs(bisection(2.4, 5, 0.05, 0.05, lambda x: 2 * x ** 3 - 15 * x ** 2 + 36 * x - 14) -
                            golden_ratio(2.4, 5, 0.05, lambda x: 2 * x ** 3 - 15 * x ** 2 + 36 * x - 14)) < 0.3)
        self.assertTrue(abs(bisection(0.5, 6, 0.05, 0.05, lambda x: 16 / x + x + 3) -
                            golden_ratio(0.5, 6, 0.05, lambda x: 16 / x + x + 3)) < 0.3)
        self.assertTrue(abs(bisection(1.1, 5, 0.05, 0.05, lambda x: x ** 3 / (x ** 2 - 1)) -
                            golden_ratio(1.1, 5, 0.05, lambda x: x ** 3 / (x ** 2 - 1))) < 0.3)
        self.assertTrue(abs(bisection(0, 15, 0.05, 0.05, lambda x: x ** 2 - 34 * x + 140 * math.log(x) - 10) -
                            golden_ratio(8, 15, 0.05, lambda x: x ** 2 - 34 * x + 140 * math.log(x) - 10)) < 0.3)
        self.assertTrue(abs(bisection(-20, 17, 0.05, 0.05, lambda x: (x ** 4) / 4 - 2 * x ** 2 + 1) -
                            golden_ratio(-20, 17, 0.05, lambda x: (x ** 4) / 4 - 2 * x ** 2 + 1)) < 0.3)
        self.assertTrue(abs(bisection(0, 5, 0.05, 0.05, lambda x: math.cos(x)) -
                            golden_ratio(0, 5, 0.05, lambda x: math.cos(x))) < 0.3)
        self.assertTrue(abs(bisection(-7, 5, 0.05, 0.05, lambda x: (x - 1) * math.exp(3 * x)) -
                            golden_ratio(-7, 5, 0.05, lambda x: (x - 1) * math.exp(3 * x))) < 0.3)


if __name__ == "__main__":
    unittest.main()

