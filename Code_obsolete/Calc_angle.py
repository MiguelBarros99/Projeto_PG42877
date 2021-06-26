import math

class Calc_angles:

    def __init__(self, Coords):
        '''
        :param Coords: recebe uma lista com coordenadas contendo as coordenas em X;Y;Z para cada atomo
        p.exp,:[[1,2,3],[4,5,6],[7,8,9]]
        '''
        self.coords = Coords
        Co1, Co2, Co3, Co4 = Coords
        self.Co1 = self.get_xyz(Co1)
        self.Co2 = self.get_xyz(Co2)
        self.Co3 = self.get_xyz(Co3)
        self.Co4 = self.get_xyz(Co4)

    def get_xyz(self, Co):  # separar todos os pontos das coords
        p_x,p_y,p_z  = Co[0],Co[1],Co[2]
        return p_x,p_y,p_z

    def angulo(self, DD1,DD3):
        lenA = self.Distance3d([0,0,0], DD1)
        lenB = self.Distance3d([0,0,0], DD3)
        DProd = self.DotProduct3d(DD1, DD3)

        CosA = self.CosAfromDotProduct(DProd, lenA, lenB)
        return math.acos(CosA) * (180 / math.pi)

    def Distance3d(self, P1, P2):
        tmp = (P2[0] - P1[0])**2 + (P2[1] - P1[1])**2 + (P2[2] - P1[2])**2
        return math.sqrt(tmp)

    def CosAfromDotProduct(self, aDP, aLine1Len, aLine2Len):
        try:
            return aDP / (aLine1Len * aLine2Len)
        except ZeroDivisionError as errormsg:
            print('An error occured: ', errormsg)
            return 0

    def VectorSubtract(self, P1,P2):
        return [P1[0] - P2[0], P1[1] - P2[1], P1[2] - P2[2]]

    def CrossProduct(self, V1,V2):
        return [V1[1] * V2[2] - V1[2] * V2[1], V1[2] * V2[0] - V1[0] * V2[2], V1[0] * V2[1] - V1[1] * V2[0]]

    def DotProduct3d(self, P1,P2):# Return Dot product of two vectors
        return P1[0] * P2[0] + P1[1] * P2[1] + P1[2] * P2[2]

    def DihedralAngle(self):
        V_dist32 = self.VectorSubtract(self.Co3, self.Co2)
        V_dist12 = self.VectorSubtract(self.Co1, self.Co2)
        V_dist43 = self.VectorSubtract(self.Co4, self.Co3)

        arr_dd1 = self.CrossProduct(V_dist32 , V_dist12)
        arr_dd3 = self.CrossProduct(V_dist32 , V_dist43)

        r = self.angulo(arr_dd1, arr_dd3)
        arr_pos_d = self.CrossProduct(V_dist32, arr_dd1)

        if (self.DotProduct3d(arr_dd3, arr_pos_d) < 0):
            r = -r
        return r
