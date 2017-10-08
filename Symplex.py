# coding=utf-8

import unittest
import numpy as np
import warnings
warnings.filterwarnings('ignore')

class Symplex:
    
    import numpy as np
    
    def __init__(self):
        self.__matrix = np.array([])
    
    def read_from_file(self, filename, act = "max"):
        tmp_mx = np.array([])
        with open(filename, "r") as f:
            flag = True
            for line in f.readlines():
                append = np.array([float(i) for i in line.split()])
                if flag:
                    tar_func = append
                    flag = False
                    continue
                if len(tmp_mx) == 0:
                    tmp_mx = append
                    continue
                if len(append)!=0:
                    tmp_mx = np.vstack((tmp_mx, append))
            if act == "max":
                for i in range(len(tar_func)):
                    tar_func[i]*=-1
        add_mx = np.diag(np.ones((tmp_mx.shape[0]+1)))
        tmp = np.vstack((tmp_mx, tar_func))
        result = np.hstack((np.hstack((tmp[:, :-1], add_mx)), tmp[:, -1:]))
        self.__matrix = result
        self.__transf = False
        
    def transform(self):
        
        self.__row_names = {}
        while min(self.__matrix[-1, :-1]) < 0:
            neg_col = np.argmin(self.__matrix[-1, :-1])
            temper = []
            for i in range(self.__matrix.shape[0]-1):
                temper.append(self.__matrix[i, -1] / self.__matrix[i, neg_col])
            neg_row = np.argmin(temper)
            self.__row_names[neg_row] = neg_col
            divider = self.__matrix[neg_row, neg_col]
            for i in range(self.__matrix.shape[1]):
                self.__matrix[neg_row, i] /= divider
            for i in range(self.__matrix.shape[0]):
                if i == neg_row:
                    continue
                coef = self.__matrix[i, neg_col]
                for j in range(self.__matrix.shape[1]):
                    self.__matrix[i, j] -= self.__matrix[neg_row, j]*coef
        self.__transf = True
        return self.__matrix       
        
    def get_profit(self):
        if self.__transf:
            return self.__matrix[-1, -1]
        else:
            print("Not fitted")
            return np.nan
    
    def get_report(self):
        if self.__transf:
            res = {}
            res['profit'] = self.__matrix[-1, -1]
            for i in self.__row_names.keys():
                res['x_'+str(self.__row_names[i])] = self.__matrix[i, -1]
            return res
        else:
            print("Not fitted")
            return np.nan
        
    def get_matrix(self):
        return self.__matrix
    
class Tester(unittest.TestCase):

    def test_sympl(self):
        smpl = Symplex()
        smpl.read_from_file("table.txt")
        smpl.transform()
        self.assertTrue(abs(smpl.get_profit() - 4400) < 10)
        smpl.read_from_file("table_one.txt")
        smpl.transform()
        self.assertTrue(abs(smpl.get_profit() - 13) < 10)

if __name__ == "__main__":
    unittest.main()