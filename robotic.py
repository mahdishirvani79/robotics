# import numpy as np
import math
import sympy

class Robotic:

    def __init__(self, DH):
        self.DH = DH
        self.Ts = list()
        self.Rs = list()
        
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
        print(self.pinzero)

    # def computeWandV(self):
    #     th1_p = sympy("th1_p")
    #     th1_p = sympy("th1_p")
    #     th1_p = sympy("th1_p")


    # def fkinCalc(self):
    #     TT = np.eye(4)
    #     for l in range(np.size(self.links,0)):
    #         T = self.makeT(self.links[l])
    #         TT = TT.dot(T)
            
    #     return TT
    
    # def fkin(self,joints):
    #     # joints = list(map(np.deg2rad, joints))
    #     for i in range(np.size(joints)):
    #         if self.links[i][4] == 0:
    #             self.links[i][3] = joints[i]
    #         else:
    #             self.links[i][2] = joints[i]
        
    #     T = self.fkinCalc()
    #     return T



def run():
    L2 = sympy.Symbol("L1")
    L3 = sympy.Symbol("L2")
    Th1 = sympy.Symbol("Th1")
    Th2 = sympy.Symbol("Th2")
    Th3 = sympy.Symbol("Th3")

    DH = sympy.Array([[0, 0, Th1, 0, 0], [math.radians(90), 0, Th2, 0, 0,],
                [0, L2, Th3, 0, 0], [0, L3, 0, 0, 0]])
    robot = Robotic(DH)
    # la = robot.fkin([30,30,30])
    robot.makeTandR()
    robot.PInZero()
    


if __name__ == "__main__":
    run()