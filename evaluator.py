from math import nan
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

SPEED_OF_LIGHT = 299792458
MIRROR_AMGLE_SIZE = 10

def transrorm_waist(waist2, radius, lmbd, transform_mx):
    res_waist2 = 0
    res_radius = 0
    if transform_mx[0, 1] < 0.000000001:
        res_waist2 = transform_mx[0, 0] * waist2
        res_radius = 1 / (transform_mx[1, 0] / transform_mx[0, 0] + transform_mx[1, 1] / transform_mx[0, 0] / radius)
    else:
        res_waist2 = waist2 * ((transform_mx[0, 0] + transform_mx[0, 1] / radius) ** 2 + (lmbd * transform_mx[0, 1] / np.pi / waist2) ** 2)
        res_radius = transform_mx[0, 1] / (transform_mx[1, 1] - waist2 / res_waist2 * (transform_mx[0, 0] + transform_mx[0, 1] / radius))
    return res_waist2, res_radius

class Mirror:
    def __init__(self, coord, angle, radius):
        self.coord = coord
        self.angle = angle
        self.radius = radius
        self.z_central_coord = nan
        self.x_central_coord = nan
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

    def set_central_coord(self, z_central_coord, x_central_coord):
        self.x_central_coord = x_central_coord
        self.z_central_coord = z_central_coord

    def get_central_coord(self):
        return np.array([self.z_central_coord, self.x_central_coord])

class System:
    def __init__ (self, num_of_mirrors, *args):
        self.num_of_mirrors = num_of_mirrors
        self.elems = []
        for arg in args:
            self.elems.append(arg)
        self.elems.sort()
        tmp_coord = np.array([0, 0])
        tmp_angle_rad = self.elems[0].get_angle() * np.pi / 180
        self.elems[0].set_central_coord(np.cos(tmp_angle_rad) * self.elems[0].get_radius(), np.sin(tmp_angle_rad) * self.elems[0].get_radius())
        tmp_diraction = np.array([1, 0])
        for i in range(1, self.num_of_mirrors):
            tmp_angle_rad = self.elems[i].get_angle() * np.pi / 180
            tmp_coord = tmp_coord + (self.elems[i].get_coord() - self.elems[i - 1].get_coord()) * tmp_diraction
            tmp_diraction = np.array([-1 * tmp_diraction[0] * np.cos(tmp_angle_rad) + tmp_diraction[1] * np.sin(tmp_angle_rad), -1 * (tmp_diraction[0] * np.sin(tmp_angle_rad) + tmp_diraction[1] * np.cos(tmp_angle_rad))])
            tmp_center = tmp_coord + tmp_diraction * self.elems[i].get_radius()
            tmp_diraction = np.array([tmp_diraction[0] * np.cos(tmp_angle_rad) - tmp_diraction[1] * np.sin(tmp_angle_rad), tmp_diraction[0] * np.sin(tmp_angle_rad) + tmp_diraction[1] * np.cos(tmp_angle_rad)])
            self.elems[i].set_central_coord(tmp_center[0], tmp_center[1])

    def is_consistent(self):
        mrkr = 0
        for i in range(self.num_of_mirrors):
            mrkr += self.elems[i].is_terminating()
        return (mrkr == 2 and self.elems[0].is_terminating() and self.elems[-1].is_terminating() and self.elems[0].get_coord() < 0.000000001)

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

    def transverse_split(self):
        mx_sagittal = self.st_matrix_sagittal()
        mx_tangential = self.st_matrix_tangential()
        delta_phy = (np.arctan2(2 * np.sqrt(-1 * mx_sagittal[1, 0] * mx_sagittal[1, 1] * mx_sagittal[0, 0] * mx_sagittal[0, 1]), 
                                (mx_sagittal[0, 0] * mx_sagittal[1, 1] + mx_sagittal[0, 1] * mx_sagittal[1, 0])) +
                    np.arctan2(2 * np.sqrt(-1 * mx_tangential[1, 0] * mx_tangential[1, 1] * mx_tangential[0, 0] * mx_tangential[0, 1]), 
                                (mx_tangential[0, 0] * mx_tangential[1, 1] + mx_tangential[0, 1] * mx_tangential[1, 0]))) / 2
        return SPEED_OF_LIGHT / 2 / self.get_length() * delta_phy

    def longitude_split(self):
        return np.pi * SPEED_OF_LIGHT / self.get_length()

    def system_scheme(self):
        fig, ax = plt.subplots()

        clctd_angle = 0
        center = np.array([0, 0])
        position = np.array([0, 0])

        for i in range(self.num_of_mirrors - 1):
            clctd_angle += self.elems[i].get_angle()
            center = position + self.elems[i].get_radius() * np.array([np.cos(clctd_angle * np.pi / 180), np.sin(clctd_angle * np.pi / 180)])
            clctd_angle += 180
            tmp_arc = patches.Arc(center, 2 * self.elems[i].get_radius(), 2 * self.elems[i].get_radius(), 0, clctd_angle - MIRROR_AMGLE_SIZE, clctd_angle + MIRROR_AMGLE_SIZE, color= 'blue')
            ax.add_patch(tmp_arc)
            if i == 0:
                d_position = (self.elems[i + 1].get_coord() - self.elems[i].get_coord()) * np.array([1, 0])
                clctd_angle -= self.elems[i].get_angle()
            else:
                d_position = -1 * (self.elems[i + 1].get_coord() - self.elems[i].get_coord()) * np.array([np.cos(np.pi * clctd_angle / 180), np.sin(np.pi * clctd_angle / 180)])
                clctd_angle += self.elems[i].get_angle()
            tmp_line = mlines.Line2D([position[0], position[0] + d_position[0]], [position[1], position[1] + d_position[1]])
            ax.add_line(tmp_line)
            position = position + d_position
        
        clctd_angle += self.elems[-1].get_angle()
        center = position + self.elems[-1].get_radius() * np.array([np.cos(clctd_angle * np.pi / 180), np.sin(clctd_angle * np.pi / 180)])
        clctd_angle += 180
        tmp_arc = patches.Arc(center, 2 * self.elems[-1].get_radius(), 2 * self.elems[-1].get_radius(), 0, clctd_angle - MIRROR_AMGLE_SIZE, clctd_angle + MIRROR_AMGLE_SIZE, color= 'blue')
        ax.add_patch(tmp_arc)

        plt.axis('equal')
        plt.axis('off')
        plt.show()

    def fund_lambda_choice(self, zero_lambda_approx):
        zero_omega_approx = 2 * np.pi * SPEED_OF_LIGHT / zero_lambda_approx
        n_closest = np.around((zero_omega_approx - self.transverse_split()) / self.longitude_split())
        return 2 * np.pi * SPEED_OF_LIGHT / (self.longitude_split() * n_closest + self.transverse_split())

    def waist_search(self, zero_lambda_approx):
        using_lambda = self.fund_lambda_choice(zero_lambda_approx)
        mx_tangential = self.st_matrix_tangential()
        mx_sagittal = self.st_matrix_sagittal()
        init_waists2 = np.array([using_lambda / np.pi * np.sqrt(-1 * mx_tangential[1, 1] * mx_tangential[0, 1] / mx_tangential[1, 0] / mx_tangential[0, 0]),
                                using_lambda / np.pi * np.sqrt(-1 * mx_sagittal[1, 1] * mx_sagittal[0, 1] / mx_sagittal[1, 0] / mx_sagittal[0, 0])])
        init_curv_radius = np.array([-1, -1]) * self.elems[0].get_radius()
        waist = np.array([[[0, 0], [0, 0]]])
        for i in range(self.num_of_mirrors - 1):
            beta_tangential = (using_lambda * init_curv_radius[0] / np.pi / init_waists2[0]) ** 2
            beta_sagittal = (using_lambda * init_curv_radius[1] / np.pi / init_waists2[1]) ** 2
            waist = np.concatenate((waist, np.array([[[- init_curv_radius[0] / (1 + beta_tangential), beta_tangential * np.sqrt(init_waists2[0]) / (1 + beta_tangential)],
                                                [- init_curv_radius[1] / (1 + beta_sagittal), beta_sagittal * np.sqrt(init_waists2[1]) / (1 + beta_sagittal)]]])))
            transform_mx_tangential = np.matrix([[1, (self.elems[i + 1].get_coord() - self.elems[i].get_coord())], [0, 1]])
            transform_mx_tangential = np.matmul(self.elems[i + 1].get_matrix_tangential(), transform_mx_tangential)
            transform_mx_tangential = np.matmul(self.elems[i + 1].get_matrix_tangential(), transform_mx_tangential)
            res_waist2, res_radius = transrorm_waist(init_waists2[0], init_curv_radius[0], using_lambda, transform_mx_tangential)
            init_waists2[0] = res_waist2
            init_curv_radius[0] = res_radius
            transform_mx_sagittal = np.matrix([[1, (self.elems[i + 1].get_coord() - self.elems[i].get_coord())], [0, 1]])
            transform_mx_sagittal = np.matmul(self.elems[i + 1].get_matrix_sagittal(), transform_mx_sagittal)
            transform_mx_sagittal = np.matmul(self.elems[i + 1].get_matrix_sagittal(), transform_mx_sagittal)
            res_waist2, res_radius = transrorm_waist(init_waists2[1], init_curv_radius[1], using_lambda, transform_mx_sagittal)
            init_waists2[1] = res_waist2
            init_curv_radius[1] = res_radius
        return waist[1:]
            
            

m1 = Mirror(0, 10, 0.005)
m2 = Mirror(0.009, 0, 0.005)


Sys = System(2, m1, m2)

print(Sys.is_consistent())
Sys.system_scheme()

