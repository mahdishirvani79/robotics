import numpy as np
import math


class SerialLink:

    def __init__(self, name, links):
        self.name  = name
        self.links = links
        
    def makeT(self,DH):
        T = np.array([[math.cos(DH[3]),                -math.sin(DH[3]),  0,  DH[1]],
                     [math.sin(DH[3])*math.cos(DH[0]), math.cos(DH[3])*math.cos(DH[0]), -math.sin(DH[0]), -DH[2]*math.sin(DH[0])],
                     [math.sin(DH[3])*math.sin(DH[0]), math.cos(DH[3])*math.sin(DH[0]),  math.cos(DH[0]),  DH[2]*math.cos(DH[0])],
                     [0,                                0,                                0,                1]])
        return T

    def fkinCalc(self):
        TT = np.eye(4)
        for l in range(np.size(self.links,0)):
            T = self.makeT(self.links[l])
            TT = TT.dot(T)
            
        return TT
    
    def fkin(self,joints):
        joints = list(map(np.deg2rad, joints))
        for i in range(np.size(joints)):
            if self.links[i][4] == 0:
                self.links[i][3] = joints[i]
            else:
                self.links[i][2] = joints[i]
        
        T = self.fkinCalc()
        return T



def run():
    robot_links = np.array([[0.0, 0.0, 17.0, 0.0, 0], [math.radians(90), 0.0, 0.0, 0.0, 0],
    [0, 12, 0.0, 0.0, 0]])
    robot = SerialLink("robot", robot_links)
    la = robot.fkin([30,30,30])


if __name__ == "__main__":
    run()