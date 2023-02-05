import numpy as np
import math
import sympy

class Robotic:

    def __init__(self, DH, F, M, path):
        self.DH = DH
        self.F = F
        self.M = M
        self.Ts = list()
        self.Rs = list()
        self.Ws = list()
        self.Vs = list()
        self.THp = list()
        self.f = open(path, 'w')

        
    def makesingleTandR(self,DH):
        T = np.array([[math.cos(DH[2]),                -math.sin(DH[2]),  0,  DH[1]],
                     [math.sin(DH[2])*math.cos(DH[0]), math.cos(DH[2])*math.cos(DH[0]), -math.sin(DH[0]), -DH[3]*math.sin(DH[0])],
                     [math.sin(DH[2])*math.sin(DH[0]), math.cos(DH[2])*math.sin(DH[0]),  math.cos(DH[0]),  DH[3]*math.cos(DH[0])],
                     [0,                                0,                                0,                1]]).astype(np.float16)
        self.Ts.append(T)
        print(str(len(self.Ts)) + "th T is" + str(np.round_(T, decimals=2)) + "\n")
        self.f.write(str(len(self.Ts)) + "th T is" + str((np.round_(T, decimals=2))) + "\n")
        R = T[0:3, 0:3].T
        self.Rs.append(R)
        print(str(len(self.Rs)) + "th R is" + str((np.round_(R, decimals=2))) + "\n")
        self.f.write(str(len(self.Rs)) + "th R is" + str((np.round_(R, decimals=2))) + "\n")


    def makeTandR(self):
        for i in range(self.DH.shape[0]):
            self.makesingleTandR(self.DH[i])


    def PInZero(self):
        T_all = np.eye(4)
        for T in self.Ts:
            T_all = T_all @ (T)
        self.pinzero = T_all @ np.array([0,0,0,1]).reshape(4,1)
        self.pinzero = self.pinzero[:-1]
        print("P in zero coordinates is: " + str((np.round_(self.pinzero, decimals=2))) + "\n")
        self.f.write("P in zero coordinates is: " + str((np.round_(self.pinzero, decimals=2))) + "\n")
        return self.pinzero


    def computeWandV(self):
        for i in range(self.DH.shape[0]):
            if self.DH[i,2] != 0:
                self.THp.append(sympy.Symbol(f"th{i+1}_p"))
            elif self.DH[i,2] == 0:
                self.THp.append(sympy.Symbol(f"d{i+1}_p"))
        
        self.Ws.append(sympy.Matrix([0,0,0]))
        self.Vs.append(sympy.Matrix([0,0,0]))
        for i in range(self.DH.shape[0]):
            if self.DH[i,4] == 0:
                W_new = (self.Rs[i] @ self.Ws[i]) + sympy.Matrix([0,0,self.THp[i]])
                self.Ws.append(W_new)
                V_new = self.Rs[i] @ (self.Vs[i] + self.Ws[i].cross(sympy.Matrix([self.DH[i,1], 0, 0])))
                self.Vs.append(V_new)
            elif self.DH[i,4] == 1:
                W_new = (self.Rs[i] @ self.Ws[i]) 
                self.Ws.append(W_new)
                V_new = self.Rs[i] @ (self.Vs[i] + self.Ws[i].cross(sympy.Matrix([self.DH[i,1], 0, 0])))  + sympy.Matrix([0,0,self.THp[i]])
                self.Vs.append(V_new)
        print("last W is: " + str(self.Ws[-1]) + "\n")
        print("last V is: " + str(self.Vs[-1]) + "\n")
        self.f.write("last W is: " + str(self.Ws[-1]) + "\n")
        self.f.write("last V is: " + str(self.Vs[-1]) + "\n")


    def computeWandVinZero(self):
        WInZero = sympy.eye(3)
        VInZero = sympy.eye(3)
        for i in self.Rs:
            i = i.T
            WInZero = WInZero @ i
            VInZero = VInZero @ i
        self.WInZero = WInZero @ self.Ws[-1]
        self.VInZero = VInZero @ self.Vs[-1]
        print("W in zero is:", str(self.WInZero) + "\n")
        print("V in zero is:", str(self.VInZero) + "\n")


    def computeJacobian(self):
        Jcalc = sympy.zeros(6,1)
        Jcalc[0:3, :] = self.VInZero
        Jcalc[3:6, :] = self.WInZero
        self.Jacobian = sympy.zeros(6,3)
        for i in range(self.Jacobian.shape[0]):
            for j in range(self.Jacobian.shape[1]):
                self.Jacobian[i,j] = Jcalc[i].diff(self.THp[j])
        self.Jacobian = np.array(self.Jacobian.tolist(), dtype=np.float16)
        print("Jacobian matrix is: " + str((np.round_(self.Jacobian, decimals=2))) + "\n")
        self.f.write("Jacobian matrix is: " + str((np.round_(self.Jacobian, decimals=2))) + "\n")


    def computeTaw(self):
        FInZero = np.eye(3)
        MInZero = np.eye(3)
        for i in self.Rs:
            i = i.T
            FInZero = FInZero @ i
            MInZero = MInZero @ i
        FInZero = FInZero @ self.F
        MInZero = MInZero @ self.M
        print("F in zero is:", str(np.round_(FInZero, decimals=2)), "\n")
        self.f.write("F in zero is:" + str(np.round_(FInZero, decimals=2)) + "\n")
        print("M in zero is:", str(np.round_(MInZero, decimals=2)), "\n")
        self.f.write("M in zero is:" + str(np.round_(MInZero, decimals=2)) + "\n")
        tempMat = np.zeros((6,1))
        tempMat[0:3, :] = FInZero
        tempMat[3:6, :] = MInZero
        self.Taw = self.Jacobian.T @ tempMat
        print("Taw is: " + str(np.round_(self.Taw, decimals=2)) + "\n")
        self.f.write("Taw is: " + str(np.round_(self.Taw, decimals=2)) + "\n")
    

def run():
    print("number of joints: ")
    n = int(input())
    print("inter denavit hartenberg:")
    DH = list()
    for _ in range(n):
        lst = input().split()
        DH.append(lst)
    for i in range(n):
        for j in range(5):
            if DH[i][j].isnumeric():
                if j == 0 or j == 2:
                    DH[i][j] = math.radians(int(DH[i][j]))
                else:
                    DH[i][j] = float(DH[i][j])
    file_name = "./out_numeric.txt"
    DH = np.array(DH)
    print("inter F parameters:")
    F = input().split()
    F = [float(f) for f in F]
    F = np.array(F).reshape(3,1)
    print("inter M parameters:")
    M = input().split()
    M = [float(m) for m in M]
    M = np.array(M).reshape(3,1)

    print("\n")
    robot = Robotic(DH, F, M, file_name)
    robot.makeTandR()
    robot.PInZero()
    robot.computeWandV()
    robot.computeWandVinZero()
    robot.computeJacobian()
    robot.computeTaw()


if __name__ == "__main__":
    run()