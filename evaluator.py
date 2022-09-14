import numpy as np

SPEED_OF_LIGHT = 299792458

class Mirror:
    def __init__(self, coord, angle, radius):
        self.coord = coord
        self.angle = angle
        self.radius = radius
        self.terminating = (self.angle < 0.000000001)

    def __lt__(self, x):
        return self.coord < x.coord

    def get_coord(self):
        return self.coord

    def get_angle(self):
        return self.angle

    def get_radius(self):
        return self.radius

    def get_matrix_sagittal(self):
        return np.matrix([[1, 0],[-1 * np.cos(np.pi / 180 * self.angle) / self.radius, 1]])

    def get_matrix_tangential(self):
        return np.matrix([[1, 0],[-1 / np.cos(np.pi / 180 * self.angle) / self.radius, 1]])

    def is_terminating(self):
        return self.terminating


class System:
    def __init__ (self, num_of_mirrors, *args):
        self.num_of_mirrors = num_of_mirrors
        self.elems = []
        for arg in args:
            self.elems.append(arg)
        self.elems.sort()

    def is_consistent(self):
        mrkr = 0
        for i in range(self.num_of_mirrors):
            mrkr += self.elems[i].is_terminating()
        return (mrkr == 2 and self.elems[0].is_terminating() and self.elems[-1].is_terminating())

    def st_matrix_sagittal(self):
        res = self.elems[0].get_matrix_sagittal()
        for i in range(1, self.num_of_mirrors - 1):
            res = np.matmul(np.matrix([[1, (self.elems[i].get_coord() - self.elems[i - 1].get_coord())], [0, 1]]), res)
            res = np.matmul(self.elems[i].get_matrix_sagittal(), res)
            res = np.matmul(self.elems[i].get_matrix_sagittal(), res)
        res = np.matmul(np.matrix([[1, (self.elems[-1].get_coord() - self.elems[-2].get_coord())], [0, 1]]), res)
        res = np.matmul(self.elems[-1].get_matrix_sagittal(), res)
        return res

    def st_matrix_tangential(self):
        res = self.elems[0].get_matrix_tangential()
        for i in range(1, self.num_of_mirrors - 1):
            res = np.matmul(np.matrix([[1, (self.elems[i].get_coord() - self.elems[i - 1].get_coord())], [0, 1]]), res)
            res = np.matmul(self.elems[i].get_matrix_tangential(), res)
            res = np.matmul(self.elems[i].get_matrix_tangential(), res)
        res = np.matmul(np.matrix([[1, (self.elems[-1].get_coord() - self.elems[-2].get_coord())], [0, 1]]), res)
        res = np.matmul(self.elems[-1].get_matrix_tangential(), res)
        return res

    def is_g_stable_sagittal(self):
        mx = self.st_matrix_sagittal()
        return mx[0,0] * mx[1, 0] * mx[0, 1] * mx[1, 1] < 0

    def is_g_stable_tangential(self):
        mx = self.st_matrix_tangential()
        return mx[0,0] * mx[1, 0] * mx[0, 1] * mx[1, 1] < 0

    def get_length(self):
        return self.elems[-1].get_coord() - self.elems[0].get_coord()

    def first_lambda_approx(self, zero_lambda_approx):
        resonator_length = self.get_length()
        n_approx = np.around(resonator_length * 2 / zero_lambda_approx)
        return 2 * resonator_length / n_approx

    def transverse_split (self):
        mx_sagittal = self.st_matrix_sagittal()
        mx_tangential = self.st_matrix_tangential()
        delta_phy = (np.arctan2(2 * np.sqrt(-1 * mx_sagittal[1, 0] * mx_sagittal[1, 1] * mx_sagittal[0, 0] * mx_sagittal[0, 1]), 
                                (mx_sagittal[0, 0] * mx_sagittal[1, 1] + mx_sagittal[0, 1] * mx_sagittal[1, 0])) +
                    np.arctan2(2 * np.sqrt(-1 * mx_tangential[1, 0] * mx_tangential[1, 1] * mx_tangential[0, 0] * mx_tangential[0, 1]), 
                                (mx_tangential[0, 0] * mx_tangential[1, 1] + mx_tangential[0, 1] * mx_tangential[1, 0]))) / 2
        return SPEED_OF_LIGHT / 2 / self.get_length() * delta_phy

    def longitude_split(self):
        return np.pi * SPEED_OF_LIGHT / self.get_length()

m1 = Mirror(0, 0, 1)
m2 = Mirror(1.99, 0, 1)


Sys = System(2, m1, m2)

print(Sys.is_consistent())
print(Sys.st_matrix_sagittal())
print(Sys.is_g_stable_sagittal())
print(Sys.st_matrix_tangential())
print(Sys.is_g_stable_tangential())
print(Sys.get_length())
print(Sys.longitude_split())
print(Sys.transverse_split())
