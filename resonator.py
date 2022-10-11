from math import nan
from math import isnan
from turtle import color
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

SPEED_OF_LIGHT = 299792458
MIRROR_ANGLE_SIZE = 10
ANGLE_ACCURACY = 10 ** -13
ZERO_ACCURACY = 10 ** -13

def transrorm_waist(waist2, radius, lmbd, transform_mx):
    res_waist2 = 0
    res_radius = 0
    if transform_mx[0, 1] < ZERO_ACCURACY:
        res_waist2 = transform_mx[0, 0] * waist2
        res_radius = 1 / (transform_mx[1, 0] / transform_mx[0, 0] + transform_mx[1, 1] / transform_mx[0, 0] / radius)
    else:
        res_waist2 = waist2 * ((transform_mx[0, 0] + transform_mx[0, 1] / radius) ** 2 + (lmbd * transform_mx[0, 1] / np.pi / waist2) ** 2)
        res_radius = transform_mx[0, 1] / (transform_mx[1, 1] - waist2 / res_waist2 * (transform_mx[0, 0] + transform_mx[0, 1] / radius))
    return res_waist2, res_radius

def quadratic_solver(A, B, C):
    if np.abs(A) < ZERO_ACCURACY:
        return np.array([-1 * C / B, nan])
    D = B * B - 4 * A * C
    if D >= 0:
        x1 = (-1 * B + np.sqrt(D)) / 2 / A
        x2 = (-1 * B - np.sqrt(D)) / 2 / A
        return np.array([x1, x2])
    else:
        return np.array([nan, nan])

class Mirror:
    def __init__(self, coord, angle, radius):
        self.coord = coord
        self.angle = angle
        self.radius = radius
        self.z_central_coord = nan
        self.x_central_coord = nan
        self.terminating = (np.abs(self.angle) < ANGLE_ACCURACY)

    def __lt__(self, x):
        return self.coord < x.coord

    def get_matrix_sagittal(self):
        return np.matrix([[1, 0],[-1 * np.cos(self.angle) / self.radius, 1]])

    def get_matrix_tangential(self):
        return np.matrix([[1, 0],[-1 / np.cos(self.angle) / self.radius, 1]])

    def set_central_coord(self, z_central_coord, x_central_coord):
        self.x_central_coord = x_central_coord
        self.z_central_coord = z_central_coord

    def get_central_coord(self):
        return np.array([self.z_central_coord, self.x_central_coord])

class Resonator:
    def __init__ (self, num_of_mirrors, *args):
        self.num_of_mirrors = num_of_mirrors
        self.elems = []
        for arg in args:
            self.elems.append(arg)
        self.elems.sort()
        self.refresh()
        

    def refresh(self):
        tmp_coord = np.array([0, 0])
        tmp_angle_rad = self.elems[0].angle
        self.elems[0].set_central_coord(np.cos(tmp_angle_rad) * self.elems[0].radius, np.sin(tmp_angle_rad) * self.elems[0].radius)
        self.elems[0].terminating = (np.abs(self.elems[0].angle) < ANGLE_ACCURACY)
        tmp_diraction = np.array([1, 0])
        for i in range(1, self.num_of_mirrors):
            tmp_angle_rad = self.elems[i].angle
            tmp_coord = tmp_coord + (self.elems[i].coord - self.elems[i - 1].coord) * tmp_diraction
            tmp_diraction = np.array([-1 * tmp_diraction[0] * np.cos(tmp_angle_rad) + tmp_diraction[1] * np.sin(tmp_angle_rad), -1 * (tmp_diraction[0] * np.sin(tmp_angle_rad) + tmp_diraction[1] * np.cos(tmp_angle_rad))])
            tmp_center = tmp_coord + tmp_diraction * self.elems[i].radius
            tmp_diraction = np.array([tmp_diraction[0] * np.cos(tmp_angle_rad) - tmp_diraction[1] * np.sin(tmp_angle_rad), tmp_diraction[0] * np.sin(tmp_angle_rad) + tmp_diraction[1] * np.cos(tmp_angle_rad)])
            self.elems[i].set_central_coord(tmp_center[0], tmp_center[1])
            self.elems[i].terminating = (np.abs(self.elems[i].angle) < ANGLE_ACCURACY)

    def is_consistent(self):
        mrkr = 0
        for i in range(self.num_of_mirrors):
            mrkr += self.elems[i].terminating
        return (mrkr == 2 and self.elems[0].terminating and self.elems[-1].terminating and self.elems[0].coord < 0.000000001)

    def st_matrix_sagittal(self):
        res = self.elems[0].get_matrix_sagittal()
        for i in range(1, self.num_of_mirrors - 1):
            res = np.matmul(np.matrix([[1, (self.elems[i].coord - self.elems[i - 1].coord)], [0, 1]]), res)
            res = np.matmul(self.elems[i].get_matrix_sagittal(), res)
            res = np.matmul(self.elems[i].get_matrix_sagittal(), res)
        res = np.matmul(np.matrix([[1, (self.elems[-1].coord - self.elems[-2].coord)], [0, 1]]), res)
        res = np.matmul(self.elems[-1].get_matrix_sagittal(), res)
        return res

    def st_matrix_tangential(self):
        res = self.elems[0].get_matrix_tangential()
        for i in range(1, self.num_of_mirrors - 1):
            res = np.matmul(np.matrix([[1, (self.elems[i].coord - self.elems[i - 1].coord)], [0, 1]]), res)
            res = np.matmul(self.elems[i].get_matrix_tangential(), res)
            res = np.matmul(self.elems[i].get_matrix_tangential(), res)
        res = np.matmul(np.matrix([[1, (self.elems[-1].coord - self.elems[-2].coord)], [0, 1]]), res)
        res = np.matmul(self.elems[-1].get_matrix_tangential(), res)
        return res

    def is_g_stable_sagittal(self):
        mx = self.st_matrix_sagittal()
        return mx[0,0] * mx[1, 0] * mx[0, 1] * mx[1, 1] < ZERO_ACCURACY

    def is_g_stable_tangential(self):
        mx = self.st_matrix_tangential()
        return mx[0,0] * mx[1, 0] * mx[0, 1] * mx[1, 1] < ZERO_ACCURACY

    def get_length(self):
        return self.elems[-1].coord - self.elems[0].coord

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
            clctd_angle += self.elems[i].angle
            center = position + self.elems[i].radius * np.array([np.cos(clctd_angle), np.sin(clctd_angle)])
            clctd_angle += 180
            tmp_arc = patches.Arc(center, 2 * self.elems[i].radius, 2 * self.elems[i].radius, 0, clctd_angle - MIRROR_ANGLE_SIZE, clctd_angle + MIRROR_ANGLE_SIZE, color= 'blue')
            ax.add_patch(tmp_arc)
            if i == 0:
                clctd_angle -= self.elems[i].angle
                d_position = (self.elems[i + 1].coord - self.elems[i].coord) * np.array([1, 0])
            else:
                clctd_angle += self.elems[i].angle
                d_position = -1 * (self.elems[i + 1].coord - self.elems[i].coord) * np.array([np.cos(np.pi * clctd_angle / 180), np.sin(np.pi * clctd_angle / 180)])
            tmp_line = mlines.Line2D([position[0], position[0] + d_position[0]], [position[1], position[1] + d_position[1]])
            ax.add_line(tmp_line)
            position = position + d_position
        
        clctd_angle += self.elems[-1].angle
        center = position + self.elems[-1].radius * np.array([np.cos(clctd_angle), np.sin(clctd_angle)])
        clctd_angle += 180
        tmp_arc = patches.Arc(center, 2 * self.elems[-1].radius, 2 * self.elems[-1].radius, 0, clctd_angle - MIRROR_ANGLE_SIZE, clctd_angle + MIRROR_ANGLE_SIZE, color= 'blue')
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
        init_curv_radius = np.array([-1, -1]) * self.elems[0].radius
        waist = np.array([[[0, 0, 0], [0, 0, 0]]])
        for i in range(self.num_of_mirrors - 1):
            beta_tangential = (using_lambda * init_curv_radius[0] / np.pi / init_waists2[0]) ** 2
            beta_sagittal = (using_lambda * init_curv_radius[1] / np.pi / init_waists2[1]) ** 2
            waist = np.concatenate((waist, np.array([[[- init_curv_radius[0] / (1 + beta_tangential), np.sqrt(beta_tangential) * np.sqrt(init_waists2[0]) / np.sqrt(1 + beta_tangential), 0],
                                                [- init_curv_radius[1] / (1 + beta_sagittal), np.sqrt(beta_sagittal) * np.sqrt(init_waists2[1]) / np.sqrt(1 + beta_sagittal), 0]]])))
            waist[i + 1, 0, 2] = using_lambda / np.pi / waist[i + 1, 0, 1]
            waist[i + 1, 1, 2] = using_lambda / np.pi / waist[i + 1, 1, 1]
            transform_mx_tangential = np.matrix([[1, (self.elems[i + 1].coord - self.elems[i].coord)], [0, 1]])
            transform_mx_tangential = np.matmul(self.elems[i + 1].get_matrix_tangential(), transform_mx_tangential)
            transform_mx_tangential = np.matmul(self.elems[i + 1].get_matrix_tangential(), transform_mx_tangential)
            res_waist2, res_radius = transrorm_waist(init_waists2[0], init_curv_radius[0], using_lambda, transform_mx_tangential)
            init_waists2[0] = res_waist2
            init_curv_radius[0] = res_radius
            transform_mx_sagittal = np.matrix([[1, (self.elems[i + 1].coord - self.elems[i].coord)], [0, 1]])
            transform_mx_sagittal = np.matmul(self.elems[i + 1].get_matrix_sagittal(), transform_mx_sagittal)
            transform_mx_sagittal = np.matmul(self.elems[i + 1].get_matrix_sagittal(), transform_mx_sagittal)
            res_waist2, res_radius = transrorm_waist(init_waists2[1], init_curv_radius[1], using_lambda, transform_mx_sagittal)
            init_waists2[1] = res_waist2
            init_curv_radius[1] = res_radius
        return waist[1:]    #[waist coord along axis, waist, Numerical Aperture]

    def realign(self):
        flag = True
        rotation = 0
        while not self.is_consistent():

            r_transform_mx = np.matrix([[1, 0], [0, 1]])
            for i in range(1, self.num_of_mirrors - 1):
                r_transform_mx = np.matmul(np.matrix([[1, self.elems[i].coord - self.elems[i - 1].coord], [0, 1]]), r_transform_mx)
                r_transform_mx = np.matmul(self.elems[i].get_matrix_tangential(), r_transform_mx)
                r_transform_mx = np.matmul(self.elems[i].get_matrix_tangential(), r_transform_mx)
            r_transform_mx = np.matmul(np.matrix([[1, (self.elems[-1].coord - self.elems[-2].coord)], [0, 1]]), r_transform_mx)

            z_central_start = self.elems[0].radius * np.cos(self.elems[0].angle)
            x_central_start = self.elems[0].radius * np.sin(self.elems[0].angle)

            self.elems[0].angle = 0
            self.elems[0].coord = 0

            term_radius = -1 * self.elems[-1].radius

            x_central_term = (-1) ** self.num_of_mirrors * self.elems[-1].radius * np.sin(self.elems[-1].angle) / (r_transform_mx[1, 0] * term_radius + r_transform_mx[0, 0])
            new_radius = (r_transform_mx[1, 1] * term_radius + r_transform_mx[0, 1]) / (r_transform_mx[1, 0] * term_radius + r_transform_mx[0, 0])
            z_central_term = np.sqrt(new_radius ** 2 - x_central_term ** 2)

            direct = np.array([np.abs(z_central_term - z_central_start), x_central_start - x_central_term])
            direct = 1 / np.sqrt(direct[0] ** 2 + direct[1] ** 2) * direct
            opt_path = self.elems[0].radius
            start_coord = np.array([z_central_start, x_central_start])

            if flag:
                rotation = np.arctan2(direct[1], direct[0])
                flag = False

            for i in range(1, self.num_of_mirrors):
                central_coord = self.elems[i].get_central_coord()

                if np.abs(direct[1]) < 10 ** -9:
                    x_solutions = np.array([start_coord[1], start_coord[1]])
                    if np.abs(x_solutions[0] - central_coord[1]) > self.elems[i].radius:
                        raise Exception("Missing mirror")
                    z_diff = np.sqrt(self.elems[i].radius ** 2 - (x_solutions[0] - central_coord[1]) ** 2)
                    z_solutions = np.array([central_coord[0] + z_diff, central_coord[0] - z_diff])
                else:
                    A_coef = 1
                    B_coef = 2 * (direct[0] * direct[1] * (start_coord[0] - central_coord[0]) - direct[0] ** 2 * start_coord[1] - direct[1] ** 2 * central_coord[1])
                    C_coef = (direct[1] ** 2 * (start_coord[0] - central_coord[0]) ** 2 - 2 * direct[0] * direct[1] * (start_coord[0] - central_coord[0]) * start_coord[1]
                             + direct[0] ** 2 * start_coord[1] ** 2 + direct[1] ** 2 * central_coord[1] ** 2 - direct[1] ** 2 * self.elems[i].radius ** 2)
                    x_solutions = quadratic_solver(A_coef, B_coef, C_coef)
                    if np.isnan(x_solutions[0]):
                        raise Exception("Missing mirror")
                    z_solutions = start_coord[0] * np.array([1, 1]) + direct[0] / direct[1] * (x_solutions - start_coord[1] * np.array([1, 1]))
            
                cntrl_direct = np.array([[central_coord[0] - z_solutions[0], central_coord[1] - x_solutions[0]], [central_coord[0] - z_solutions[1], central_coord[1] - x_solutions[1]]])
                cntrl_direct[0] = cntrl_direct[0] / np.sqrt(cntrl_direct[0, 0] ** 2 + cntrl_direct[0, 1] ** 2)
                cntrl_direct[1] = cntrl_direct[1] / np.sqrt(cntrl_direct[1, 0] ** 2 + cntrl_direct[1, 1] ** 2)

                angle = 0

                if direct[0] * cntrl_direct[0, 0] + direct[1] * cntrl_direct[0, 1] < 0:
                    new_start_coord = np.array([z_solutions[0], x_solutions[0]])
                    angle = np.arctan2(direct[1] * cntrl_direct[0, 0] - direct[0] * cntrl_direct[0, 1], -1 * direct[0] * cntrl_direct[0, 0] - direct[1] * cntrl_direct[0, 1])
                    direct = direct - 2 * (direct[0] * cntrl_direct[0, 0] + direct[1] * cntrl_direct[0, 1]) * cntrl_direct[0]
                else:
                    new_start_coord = np.array([z_solutions[1], x_solutions[1]])
                    angle = np.arctan2(direct[1] * cntrl_direct[1, 0] - direct[0] * cntrl_direct[1, 1], -1 * direct[0] * cntrl_direct[1, 0] - direct[1] * cntrl_direct[1, 1])
                    direct = direct - 2 * (direct[0] * cntrl_direct[1, 0] + direct[1] * cntrl_direct[1, 1]) * cntrl_direct[1]

                d_start_coord = new_start_coord - start_coord
                opt_path += np.sqrt(d_start_coord[0] ** 2 + d_start_coord[1] ** 2)
                start_coord = new_start_coord
                self.elems[i].coord = opt_path
                self.elems[i].angle = angle
            self.refresh()
        return rotation
    
    def waist_scheme(self, zero_lambda_approx):
        using_lambda = self.fund_lambda_choice(zero_lambda_approx)
        waist = self.waist_search(zero_lambda_approx)

        number_of_steps = 1000
        
        fig, ax = plt.subplots()

        tmp_line = mlines.Line2D([0, 0], [1 / 3, -1 / 3])
        ax.add_line(tmp_line)
        tmp_line = mlines.Line2D([0, 0], [-2 / 3, -4 / 3])
        ax.add_line(tmp_line)
        length = self.get_length()
        tmp_line = mlines.Line2D([length, length], [1 / 3, -1 / 3])
        ax.add_line(tmp_line)
        tmp_line = mlines.Line2D([length, length], [-2 / 3, -4 / 3])
        ax.add_line(tmp_line)

        start = 0
        stop = self.elems[1].coord
        step = (stop - start) / (number_of_steps - 1)
        x_coords = np.arange(start, stop + step, step)
        
        x_0 = waist[0, 0, 0]
        z_R_0 = waist[0, 0, 1] ** 2 * np.pi / using_lambda
        upper_curv = waist[0, 0, 1] * np.sqrt(1 + ((x_coords - x_0) / z_R_0) ** 2)
        lower_curv = -upper_curv
        plt.plot(x_coords, upper_curv, color= 'r')
        plt.plot(x_coords, lower_curv, color= 'r')

        y_offset_1 = -1
        x_1 = waist[0, 1, 0]
        z_R_1 = waist[0, 1, 1] ** 2 * np.pi / using_lambda
        upper_curv = waist[0, 1, 1] * np.sqrt(1 + ((x_coords - x_1) / z_R_1) ** 2)
        lower_curv = -upper_curv
        upper_curv = y_offset_1 + upper_curv
        lower_curv = y_offset_1 + lower_curv
        plt.plot(x_coords, upper_curv, color= 'r')
        plt.plot(x_coords, lower_curv, color= 'r')

        for i in range(1, self.num_of_mirrors - 1):
            x_init = self.elems[i].coord
            x_term = self.elems[i + 1].coord
            tmp_line = mlines.Line2D([x_init, x_init], [1 / 3, -1 / 3])
            ax.add_line(tmp_line)
            tmp_line = mlines.Line2D([x_init, x_init], [-2 / 3, -4 / 3])
            ax.add_line(tmp_line)

            start = x_init
            stop = x_term
            step = (stop - start) / (number_of_steps - 1)
            x_coords = np.arange(start, stop + step, step)
        
            x_0 = waist[i, 0, 0]
            z_R_0 = waist[i, 0, 1] ** 2 * np.pi / using_lambda
            upper_curv = waist[i, 0, 1] * np.sqrt(1 + ((x_coords - x_0 - x_init) / z_R_0) ** 2)
            lower_curv = -upper_curv
            plt.plot(x_coords, upper_curv, color= 'r')
            plt.plot(x_coords, lower_curv, color= 'r')

            x_1 = waist[i, 1, 0]
            z_R_1 = waist[i, 1, 1] ** 2 * np.pi / using_lambda
            upper_curv = waist[i, 1, 1] * np.sqrt(1 + ((x_coords - x_1 - x_init) / z_R_1) ** 2)
            lower_curv = -upper_curv
            upper_curv = y_offset_1 + upper_curv
            lower_curv = y_offset_1 + lower_curv
            plt.plot(x_coords, upper_curv, color= 'r')
            plt.plot(x_coords, lower_curv, color= 'r')

        plt.axis('off')
        plt.show()  