# coding=utf-8

import unittest
import numpy as np
import warnings
warnings.filterwarnings('ignore')

class Func_integrator:
    
    import numpy as np
    
    class Optimizer:
    
        def __init__(self, deriv_precise=1e-6, eps=1e-5, points=1e+5):
            self.__deriv_precise = deriv_precise
            self.__eps = eps
            self.__linespace_points = points

        def derivative(self, func):
            def derivative_func(x, precise=self.__deriv_precise):
                return (func(x + precise/2) - func(x - precise/2)) / precise
            return derivative_func

        def second_derivative(self, func):
            return derivative(derivative(func))
        
        def third_derivative(self, func):
            return derivative(derivative(derivative(func)))
        
        def near_maximizer(self, func, a, b):
    
            def new_point_maximizer(func, x, delta, eps, max_old):
                points = np.linspace(x-delta, x+delta, self.__linespace_points)[1:-2]
                deltas = points[1]-points[0]

                ys = []
                for i in points:
                    ys.append(func(i))
                ys.sort()
                themax=[]
                if abs(ys[-1] - max_old) < eps or ys[-1] < max_old:
                    return max_old
                for i in ys[-5:]:
                    themax.append(new_point_maximizer(func, i, deltas, self.__eps, ys[-1]))
                return max(themax)

            first_points_init = np.linspace(a, b, self.__linespace_points*10)[1:-2]
            deltas = first_points_init[1]-first_points_init[0]

            ys = []
            for i in first_points_init:
                ys.append(func(i))
            ys.sort()

            themax = []
            for i in ys[:5]:
                themax.append(new_point_maximizer(func, i, deltas, self.__eps, ys[-1]))
            return max(themax), deltas


        def near_minimizer(self, func, a, b):
    
            def new_point_minimizer(func, x, delta, eps, min_old):
                points = np.linspace(x-delta, x+delta, self.__linespace_points)[1:-2]
                deltas = points[1]-points[0]

                ys = []
                for i in points:
                    ys.append(func(i))
                ys.sort()
                if abs(ys[0] - min_old) < eps or ys[0] > min_old:
                    return min_old
                for i in ys[:5]:
                    themins.append(new_point_minimizer(func, i, deltas, self.__eps, ys[0]))
                return min(themins)

            first_points_init = np.linspace(a, b, self.__linespace_points*10)[1:-2]
            deltas = first_points_init[1]-first_points_init[0]

            ys = []
            for i in first_points_init:
                ys.append(func(i))
            ys.sort()

            themins = []
            for i in ys[:5]:
                themins.append(new_point_minimizer(func, i, deltas, self.__eps, ys[0]))

            return min(themins), deltas

    
    def __init__(self, step=1e+5, deriv_precise=1e-6, eps=1e-5, points=1e+5):
        self.__step = step
        self.__deriv_precise = deriv_precise
        self.__eps = eps
        self.__points = points
        self.optimizer = self.Optimizer(self.__deriv_precise, eps = self.__eps, points = self.__points)

    def rect_integrator(self, func, a, b):
        h = (b - a) / float(self.__step)
        result = 0
        for i in range(int(self.__step)):
            result += func(a + h * (i + 0.5))
        result *= h
        return result
    
    def trap_integrator(self, func, a, b):
        h = (b - a) / float(self.__step)
        result = 0
        for i in range(int(self.__step)):
            result += (func(a + h * (i+1)) + func(a + h * i)) * h / 2.
        return result
    
    def sympson_integrator(self, func, a, b):
        h = (b - a) / float(self.__step)
        result = func(a) + func(b)
        for i in range(1, int(self.__step)-1):
            if i%2!=0:
                result += func(a + h * i) * 4
            else:
                result += func(a + h * i) * 2
        result *= h / 3.
        return result
    
    def trap_integrator_derriv(self, func, a, b):
        h = (b - a) / float(self.__step)
        result = 0
        for i in range(int(self.__step)):
            result += (func(a + h * (i+1)) + func(a + h * i)) * h / 2.
        max_val = self.maximize_deriv_two_func(func, a, b)
        int_val = (b - a) / 24. * h**2 * (abs(max_val[0]) + max_val[1] + self.__eps)
        return {"result": result,
                "delta": int_val,
                "conf_int": (result - int_val, result + int_val)}
    
    def sympson_integrator_quadr(self, func, a, b):
        h = (b - a) / float(self.__step)
        result = func(a) + func(b)
        for i in range(1, int(self.__step)-1):
            if i%2!=0:
                result += func(a + h * i) * 4
            else:
                result += func(a + h * i) * 2
        result *= h / 3.
        max_val = self.maximize_deriv_three_func(func, a, b)
        int_val = (b - a) / 288. * h**3 * (abs(max_val[0]) + max_val[1] + self.__eps)
        return {"result": result,
                "delta": int_val,
                "conf_int": (result - int_val, result + int_val)}
    
    def rect_integrator_quadr(self, func, a, b):
        h = (b - a) / float(self.__step)
        result = 0
        for i in range(int(self.__step)):
            result += func(a + h * (i + 0.5))
        result *= h
        int_val = self.count_diff_rect(func, a, b, h)
        return {"result": result,
                "delta": int_val,
                "conf_int": (result - int_val, result + int_val)}
    
    def trap_integrator_quadr(self, func, a, b):
        h = (b - a) / float(self.__step)
        result = 0
        for i in range(int(self.__step)):
            result += (func(a + h * (i+1)) + func(a + h * i)) * h / 2.
        int_val = self.count_diff_trap(func, a, b, h)
        return {"result": result,
                "delta": int_val,
                "conf_int": (result - int_val, result + int_val)}  
    
    
    def count_diff_trap(self, func, a, b, h):
        I_h_2 = 0
        I_h = 0
        for i in range(int(self.__step)):
            I_h += (func(a + h * (i+1)) - func(a + h * i)) * h / 2.
            I_h_2 += (func(a + 0.5 * h * (i+1)) - func(a + h * i)) * h / 4.
            I_h_2 += (func(a + h * (i+1)) - func(a + 0.5 * h * (i+1))) * h / 4.
        return abs(I_h - I_h_2) / float(2**self.__step - 1)
    
    def count_diff_rect(self, func, a, b, h):
        I_h_2 = 0
        I_h = 0
        for i in range(int(self.__step)):
            I_h += func(a + h * (i + 0.5)) * h
            I_h_2 += func(a + h * (i + 0.125)) * h / 2
            I_h_2 += func(a + h * (i + 0.375)) * h / 2
        return abs(I_h - I_h_2) / float(2**self.__step - 1)
    
    def sympson_integrator_deriv(self, func, a, b):
        h = (b - a) / float(self.__step)
        result = func(a) + func(b)
        for i in range(1, int(self.__step)-1):
            if i%2!=0:
                result += func(a + h * i) * 4
            else:
                result += func(a + h * i) * 2
        result *= h / 3.
        max_val = self.maximize_deriv_three_func(func, a, b)
        int_val = (b - a) / 288. * h**3 * (abs(max_val[0]) + max_val[1] + self.__eps)
        return {"result": result,
                "delta": int_val,
                "conf_int": (result - int_val, result + int_val)}
    
    def maximize_deriv_two_func(self, func, a, b):
        max_val = self.optimizer.second_derivative(func)
        return self.optimizer.near_maximizer(max_val, a, b)
    
    def maximize_deriv_three_func(self, func, a, b):
        max_val = self.optimizer.third_derivative(func)
        return self.optimizer.near_maximizer(max_val, a, b)
    
class Tester(unittest.TestCase):

    def test_rect(self):
        itegrator = Func_integrator()
        self.assertTrue(abs(itegrator.rect_integrator(lambda x: x, 0, 2) - 2) < 0.05)
        self.assertTrue(abs(itegrator.rect_integrator(np.sin, -1, 1) - 0) < 0.05)
        self.assertTrue(abs(itegrator.rect_integrator(np.cos, -1, 1) - 1.68) < 0.05)
        self.assertTrue(abs(itegrator.rect_integrator(lambda x: x**2, 0, 2) - 2.67) < 0.05)
        self.assertTrue(abs(itegrator.rect_integrator(lambda x: x**3, -2, 2) - 0) < 0.05)
        self.assertTrue(abs(itegrator.rect_integrator(np.log, 1, 3) - 1.296) < 0.05)
        self.assertTrue(abs(itegrator.rect_integrator(lambda x: np.sin(x)+np.cos(x)+np.sin(2*x)+np.cos(2*x), -1, 3) - 2.13) < 0.05)
        self.assertTrue(abs(itegrator.rect_integrator(lambda x: 1/x, 2, 6) - 1.099) < 0.05)
        self.assertTrue(abs(itegrator.rect_integrator(lambda x: np.exp(x), 0, 2) - 6.389) < 0.05)
        self.assertTrue(abs(itegrator.rect_integrator(lambda x: np.sin(x)/x, 4, 8) + 0.184) < 0.05)

    def test_trap(self):
        itegrator = Func_integrator()
        self.assertTrue(abs(itegrator.trap_integrator(lambda x: x, 0, 2) - 2) < 0.05)
        self.assertTrue(abs(itegrator.trap_integrator(np.sin, -1, 1) - 0) < 0.05)
        self.assertTrue(abs(itegrator.trap_integrator(np.cos, -1, 1) - 1.68) < 0.05)
        self.assertTrue(abs(itegrator.trap_integrator(lambda x: x**2, 0, 2) - 2.67) < 0.05)
        self.assertTrue(abs(itegrator.trap_integrator(lambda x: x**3, -2, 2) - 0) < 0.05)
        self.assertTrue(abs(itegrator.trap_integrator(np.log, 1, 3) - 1.296) < 0.05)
        self.assertTrue(abs(itegrator.trap_integrator(lambda x: np.sin(x)+np.cos(x)+np.sin(2*x)+np.cos(2*x), -1, 3) - 2.13) < 0.05)
        self.assertTrue(abs(itegrator.trap_integrator(lambda x: 1/x, 2, 6) - 1.099) < 0.05)
        self.assertTrue(abs(itegrator.trap_integrator(lambda x: np.exp(x), 0, 2) - 6.389) < 0.05)
        self.assertTrue(abs(itegrator.trap_integrator(lambda x: np.sin(x)/x, 4, 8) + 0.184) < 0.05)

    def test_parabol(self):
        itegrator = Func_integrator()
        self.assertTrue(abs(itegrator.sympson_integrator(lambda x: x, 0, 2) - 2) < 0.05)
        self.assertTrue(abs(itegrator.sympson_integrator(np.sin, -1, 1) - 0) < 0.05)
        self.assertTrue(abs(itegrator.sympson_integrator(np.cos, -1, 1) - 1.68) < 0.05)
        self.assertTrue(abs(itegrator.sympson_integrator(lambda x: x**2, 0, 2) - 2.67) < 0.05)
        self.assertTrue(abs(itegrator.sympson_integrator(lambda x: x**3, -2, 2) - 0) < 0.05)
        self.assertTrue(abs(itegrator.sympson_integrator(np.log, 1, 3) - 1.296) < 0.05)
        self.assertTrue(abs(itegrator.sympson_integrator(lambda x: np.sin(x)+np.cos(x)+np.sin(2*x)+np.cos(2*x), -1, 3) - 2.13) < 0.05)
        self.assertTrue(abs(itegrator.sympson_integrator(lambda x: 1/x, 2, 6) - 1.099) < 0.05)
        self.assertTrue(abs(itegrator.sympson_integrator(lambda x: np.exp(x), 0, 2) - 6.389) < 0.05)
        self.assertTrue(abs(itegrator.sympson_integrator(lambda x: np.sin(x)/x, 4, 8) + 0.184) < 0.05)


if __name__ == "__main__":
    unittest.main()