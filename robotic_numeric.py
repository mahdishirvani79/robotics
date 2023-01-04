import numpy as np
import math
import sympy

class Robotic:

    def __init__(self, DH, F, M):
        self.DH = DH
        self.F = F
        self.M = M
        self.Ts = list()
        self.Rs = list()
        self.Ws = list()
        self.Vs = list()
        self.THp = list()
        

    def makesingleTandR(self,DH):
        T = np.array([[math.cos(DH[2]),                -math.sin(DH[2]),  0,  DH[1]],
                     [math.sin(DH[2])*math.cos(DH[0]), math.cos(DH[2])*math.cos(DH[0]), -math.sin(DH[0]), -DH[3]*math.sin(DH[0])],
                     [math.sin(DH[2])*math.sin(DH[0]), math.cos(DH[2])*math.sin(DH[0]),  math.cos(DH[0]),  DH[3]*math.cos(DH[0])],
                     [0,                                0,                                0,                1]]).astype(np.float16)
        self.Ts.append(T)
        self.Rs.append(T[0:3, 0:3].T)


    def makeTandR(self):
        for i in range(self.DH.shape[0]):
            self.makesingleTandR(self.DH[i])


    def PInZero(self):
        T_all = np.eye(4)
        for T in self.Ts:
            T_all = T_all @ (T)
        self.pinzero = T_all @ np.array([0,0,0,1]).reshape(4,1)
        return self.pinzero


    def computeWandV(self):
        for i in range(self.DH.shape[0]):
            if self.DH[i,2] != 0:
                self.THp.append(sympy.Symbol(f"th{i+1}_p"))
            elif self.DH[i,2] == 0:
                self.THp.append(0)
        
        self.Ws.append(sympy.Matrix([0,0,0]))
        self.Vs.append(sympy.Matrix([0,0,0]))
        for i in range(self.DH.shape[0]):
            if self.DH[i,4] == 0:
                W_new = (self.Rs[i] @ self.Ws[i]) + sympy.Matrix([0,0,self.THp[i]])
                # print(W_new)
                self.Ws.append(W_new)
                V_new = self.Rs[i] @ (self.Vs[i] + self.Ws[i].cross(sympy.Matrix([self.DH[i,1], 0, 0])))
                # print(V_new)
                self.Vs.append(V_new)
            elif self.DH[i,4] == 1:
                W_new = (self.Rs[i] @ self.Ws[i]) 
                # print(W_new)
                self.Ws.append(W_new)
                V_new = self.Rs[i] @ (self.Vs[i] + self.Ws[i].cross(sympy.Matrix([self.DH[i,1], 0, 0])))
                # print(V_new)
                self.Vs.append(V_new)


    def computeWandVinZero(self):
        WInZero = sympy.eye(3)
        VInZero = sympy.eye(3)
        for i in self.Rs:
            i = i.T
            WInZero = WInZero @ i
            VInZero = VInZero @ i
        self.WInZero = WInZero @ self.Ws[-1]
        self.VInZero = VInZero @ self.Vs[-1]


    def computeJacobian(self):
        Jcalc = sympy.zeros(6,1)
        Jcalc[0:3, :] = self.VInZero
        Jcalc[3:6, :] = self.WInZero
        self.Jacobian = sympy.zeros(6,3)
        for i in range(self.Jacobian.shape[0]):
            for j in range(self.Jacobian.shape[1]):
                self.Jacobian[i,j] = Jcalc[i].diff(self.THp[j])
        self.Jacobian = np.array(self.Jacobian.tolist())


    def computeTaw(self):
        FInZero = np.eye(3)
        MInZero = np.eye(3)
        for i in self.Rs:
            i = i.T
            FInZero = FInZero @ i
            MInZero = MInZero @ i
        FInZero = FInZero @ self.F
        MInZero = MInZero @ self.M
        tempMat = np.zeros((6,1))
        tempMat[0:3, :] = FInZero
        tempMat[3:6, :] = MInZero
        self.Taw = self.Jacobian.T @ tempMat
    

def run():
    L2 = 10
    L3 = 10
    Th1 = math.radians(30)
    Th2 = math.radians(30)
    Th3 = math.radians(30)

    DH = np.array([[0, 0, Th1, 0, 0], [math.radians(90), 0, Th2, 0, 0,],
                [0, L2, Th3, 0, 0], [0, L3, 0, 0, 0]])

    F = np.array([-500, 0, 0]).reshape(3,1)
    M = np.array([200, 0, 0]).reshape(3,1)

    robot = Robotic(DH, F, M)
    robot.makeTandR()
    robot.PInZero()
    robot.computeWandV()
    robot.computeWandVinZero()
    robot.computeJacobian()
    robot.computeTaw()


if __name__ == "__main__":
    run()