# import numpy as np
import math
import sympy

class Robotic:

    def __init__(self, DH):
        self.DH = DH
        self.Ts = list()
        self.Rs = list()
        self.Ws = list()
        self.Vs = list()
        self.THp = list()
        
    def makesingleTandR(self,DH):
        T = sympy.Array([[sympy.cos(DH[2]),                -sympy.sin(DH[2]),  0,  DH[1]],
                     [sympy.sin(DH[2])*sympy.cos(DH[0]), sympy.cos(DH[2])*sympy.cos(DH[0]), -sympy.sin(DH[0]), -DH[3]*sympy.sin(DH[0])],
                     [sympy.sin(DH[2])*sympy.sin(DH[0]), sympy.cos(DH[2])*sympy.sin(DH[0]),  sympy.cos(DH[0]),  DH[3]*sympy.cos(DH[0])],
                     [0,                                0,                                0,                1]]).tomatrix()
        self.Ts.append(T)
        self.Rs.append(T[0:3, 0:3].T)


    def makeTandR(self):
        for i in range(self.DH.shape[0]):
            self.makesingleTandR(self.DH[i])


    def PInZero(self):
        T_all = sympy.eye(4)
        for T in self.Ts:
            T_all = T_all @ (T)
        self.pinzero = T_all @ sympy.Matrix([0,0,0,1])
        return self.pinzero

    def computeWandV(self):
        for i in range(self.DH.shape[0]):
            if self.DH[i,2] != 0:
                self.THp.append(sympy.Symbol(f"th{i+1}_p"))
            if self.DH[i,2] == 0:
                self.THp.append(0)
        
        self.Ws.append(sympy.Matrix([0,0,0]))
        self.Vs.append(sympy.Matrix([0,0,0]))
        for i in range(self.DH.shape[0]):
            if self.DH[i,4] == 0:
                W_new = (self.Rs[i] @ self.Ws[i]) + sympy.Matrix([0,0,self.THp[i]])
                print(W_new)
                self.Ws.append(W_new)
                V_new = self.Rs[i] @ (self.Vs[i] + sympy.matrix_multiply_elementwise(self.Ws[i], sympy.Matrix([self.DH[i,1], 0, 0])))
                print(V_new)
                self.Vs.append(V_new)
            if self.DH[i,4] == 1:
                pass



def run():
    L2 = sympy.Symbol("L2")
    L3 = sympy.Symbol("L3")
    Th1 = sympy.Symbol("Th1")
    Th2 = sympy.Symbol("Th2")
    Th3 = sympy.Symbol("Th3")

    DH = sympy.Array([[0, 0, Th1, 0, 0], [math.radians(90), 0, Th2, 0, 0,],
                [0, L2, Th3, 0, 0], [0, L3, 0, 0, 0]])
    robot = Robotic(DH)
    # la = robot.fkin([30,30,30])
    robot.makeTandR()
    robot.PInZero()
    robot.computeWandV()
    


if __name__ == "__main__":
    run()