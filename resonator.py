__version__ = '0.1'
__author__ = 'Nikolay Davydov'

from math import isnan, nan
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

SPEED_OF_LIGHT = 299792458
MIRROR_ANGLE_SIZE = 10
ANGLE_ACCURACY = 10 ** -13
ZERO_ACCURACY = 10 ** -13
number_of_plots_G = 0


def transrorm_waist(waist2, radius, lmbd, transform_mx):
    res_waist2 = 0
    res_radius = 0
    if transform_mx[0, 1] < ZERO_ACCURACY:
        res_waist2 = transform_mx[0, 0] * waist2
        res_radius = 1 / (transform_mx[1, 0]/transform_mx[0, 0]
                          + transform_mx[1, 1]/transform_mx[0, 0]/radius)
    else:
        res_waist2 = waist2 * ((transform_mx[0, 0] + transform_mx[0, 1]
                                /radius)**2 + (lmbd * transform_mx[0, 1]
                                / np.pi / waist2)**2)
        res_radius = transform_mx[0, 1] / (transform_mx[1, 1] - waist2
                                           /res_waist2*(transform_mx[0, 0]
                                           + transform_mx[0, 1]/radius))
    return res_waist2, res_radius


def quadratic_solver(A, B, C):
    if np.abs(A) < ZERO_ACCURACY:
        return np.array([-1 * C / B, nan])
    D = B*B - 4*A*C
    if D >= 0:
        x1 = (-1*B + np.sqrt(D)) / 2 / A
        x2 = (-1*B - np.sqrt(D)) / 2 / A
        return np.array([x1, x2])
    else:
        return np.array([nan, nan])


class Mirror:
    def __init__(self, coord, in_plane_angle, radius,
                 in_plane_angle_deviation=0,
                 out_of_plane_angle_deviation=0,
                 in_plane_coord_deviation=0,
                 along_axis_coord_deviation=0,
                 out_of_plane_coord_deviation=0):
        self.coord = coord
        self.in_plane_angle = in_plane_angle
        self.in_plane_angle_deviation = in_plane_angle_deviation
        self.out_of_plane_angle_deviation = out_of_plane_angle_deviation
        self.in_plane_coord_deviation = in_plane_coord_deviation
        self.along_axis_coord_deviation = along_axis_coord_deviation
        self.out_of_plane_coord_deviation = out_of_plane_coord_deviation
        self.radius = radius
        self.z_central_coord = nan
        self.x_central_coord = nan
        self.y_central_coord = nan

    def __lt__(self, x):
        return self.coord < x.coord

    def is_terminating(self):
        res = self.in_plane_angle < ANGLE_ACCURACY
        res = res and self.in_plane_angle_deviation < ANGLE_ACCURACY
        res = res and self.out_of_plane_angle_deviation < ANGLE_ACCURACY
        res = res and self.in_plane_coord_deviation < ZERO_ACCURACY
        res = res and self.along_axis_coord_deviation < ZERO_ACCURACY
        res = res and self.out_of_plane_coord_deviation < ZERO_ACCURACY
        return res

    def get_matrix_sagittal(self):
        return np.matrix([[1, 0],
                          [-1 * np.cos(self.in_plane_angle) / self.radius, 1]])

    def get_matrix_tangential(self):
        return np.matrix([[1, 0],
                          [-1 / np.cos(self.in_plane_angle) / self.radius, 1]])

    def set_central_coord(self, z_central_coord,
                          x_central_coord, y_central_coord):
        self.x_central_coord = x_central_coord
        self.z_central_coord = z_central_coord
        self.y_central_coord = y_central_coord

    def get_central_coord(self):
        return np.array([self.z_central_coord,
                         self.x_central_coord,
                         self.y_central_coord])


class Resonator:
    def __init__ (self, num_of_mirrors, *args):
        self.num_of_mirrors = num_of_mirrors
        self.elems = []
        for arg in args:
            self.elems.append(arg)
        self.elems.sort()
        self.refresh()
        
    def refresh(self):
        tmp_coord = np.array([0, 0, 0])
        tmp_angle_rad = (self.elems[0].in_plane_angle
                         + self.elems[0].in_plane_angle_deviation)
        self.elems[0].set_central_coord(np.cos(tmp_angle_rad)
                                        *self.elems[0].radius
                                        + self.elems[0].along_axis_coord_deviation,
                                        np.sin(tmp_angle_rad)
                                        *self.elems[0].radius
                                        + self.elems[0].in_plane_coord_deviation,
                                        np.sin(self.elems[0].out_of_plane_angle_deviation)
                                        *self.elems[0].radius
                                        + self.elems[0].out_of_plane_coord_deviation)
        tmp_diraction = np.array([1, 0, 0])
        for i in range(1, self.num_of_mirrors):
            tmp_angle_rad = self.elems[i].in_plane_angle
            tmp_coord = (tmp_coord
                         + (self.elems[i].coord - self.elems[i - 1].coord)
                         *tmp_diraction)
            tmp_diraction_cr = np.array([-1*tmp_diraction[0]
                                         *np.cos(tmp_angle_rad
                                         + self.elems[i].in_plane_angle_deviation)
                                         + tmp_diraction[1]
                                         *np.sin(tmp_angle_rad
                                         + self.elems[i].in_plane_angle_deviation),
                                         -1 * (tmp_diraction[0]*np.sin(tmp_angle_rad
                                         + self.elems[i].in_plane_angle_deviation)
                                         + tmp_diraction[1]*np.cos(tmp_angle_rad
                                         + self.elems[i].in_plane_angle_deviation)),
                                         np.sin(self.elems[i].out_of_plane_coord_deviation)])
            tmp_diraction_cr = (tmp_diraction_cr
                                / np.sqrt(tmp_diraction_cr[0]**2
                                + tmp_diraction_cr[1]**2
                                + tmp_diraction_cr[2]**2))
            tmp_center = (tmp_coord + tmp_diraction_cr*self.elems[i].radius
                          + tmp_diraction
                          *self.elems[0].along_axis_coord_deviation)
            tmp_center[0] += (-1 * tmp_diraction[1]
                              * self.elems[0].in_plane_coord_deviation)
            tmp_center[1] += (tmp_diraction[0]
                              * self.elems[0].in_plane_coord_deviation)
            tmp_center[2] +=  self.elems[0].out_of_plane_angle_deviation
            tmp_diraction = (-1 * np.array([tmp_diraction[0]
                                            *np.cos(2 * tmp_angle_rad)
                                            - tmp_diraction[1]
                                            *np.sin(2 * tmp_angle_rad),
                                            tmp_diraction[0]
                                            *np.sin(2 * tmp_angle_rad)
                                            + tmp_diraction[1]
                                            *np.cos(2 * tmp_angle_rad), 0]))
            self.elems[i].set_central_coord(tmp_center[0], tmp_center[1]
                                            + self.elems[i].in_plane_coord_deviation,
                                            tmp_center[2]
                                            + self.elems[i].out_of_plane_coord_deviation)
            self.elems[i].terminating = (np.abs(self.elems[i].in_plane_angle)
                                         < ANGLE_ACCURACY
                                         and np.abs(self.elems[i].in_plane_angle_deviation)
                                         < ANGLE_ACCURACY
                                         and np.abs(self.elems[i].out_of_plane_angle_deviation)
                                         < ANGLE_ACCURACY)

    def is_consistent(self):
        res = True
        for i in range(self.num_of_mirrors):
            res = (res and self.elems[i].in_plane_angle_deviation
                   < ANGLE_ACCURACY)
            res = (res and self.elems[i].out_of_plane_angle_deviation
                   < ANGLE_ACCURACY)
            res = (res and self.elems[i].in_plane_coord_deviation
                   < ZERO_ACCURACY)
            res = (res and self.elems[i].along_axis_coord_deviation
                   < ZERO_ACCURACY)
            res = (res and self.elems[i].out_of_plane_coord_deviation
                   < ZERO_ACCURACY)
        return (self.elems[0].is_terminating()
                and self.elems[-1].is_terminating() and self.elems[0].coord
                < ZERO_ACCURACY and res)

    def st_matrix_sagittal(self):
        res = self.elems[0].get_matrix_sagittal()
        for i in range(1, self.num_of_mirrors - 1):
            res = np.matmul(np.matrix([[1, (self.elems[i].coord
                            - self.elems[i - 1].coord)], [0, 1]]), res)
            res = np.matmul(self.elems[i].get_matrix_sagittal(), res)
            res = np.matmul(self.elems[i].get_matrix_sagittal(), res)
        res = np.matmul(np.matrix([[1, (self.elems[-1].coord
                        - self.elems[-2].coord)], [0, 1]]), res)
        res = np.matmul(self.elems[-1].get_matrix_sagittal(), res)
        return res

    def st_matrix_tangential(self):
        res = self.elems[0].get_matrix_tangential()
        for i in range(1, self.num_of_mirrors - 1):
            res = np.matmul(np.matrix([[1, (self.elems[i].coord
                            - self.elems[i - 1].coord)], [0, 1]]), res)
            res = np.matmul(self.elems[i].get_matrix_tangential(), res)
            res = np.matmul(self.elems[i].get_matrix_tangential(), res)
        res = np.matmul(np.matrix([[1, (self.elems[-1].coord
                        - self.elems[-2].coord)], [0, 1]]), res)
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
        delta_phy = (np.arctan2(2 * np.sqrt(-1 * mx_sagittal[1, 0]
                     * mx_sagittal[1, 1] * mx_sagittal[0, 0]
                     * mx_sagittal[0, 1]),
                     (mx_sagittal[0, 0]*mx_sagittal[1, 1] + mx_sagittal[0, 1]
                     *mx_sagittal[1, 0]))
                     + np.arctan2(2 * np.sqrt(-1 * mx_tangential[1, 0]
                     * mx_tangential[1, 1] * mx_tangential[0, 0]
                     * mx_tangential[0, 1]), 
                     (mx_tangential[0, 0]*mx_tangential[1, 1]
                     + mx_tangential[0, 1]*mx_tangential[1, 0]))) / 2
        return SPEED_OF_LIGHT / 2 / self.get_length() * delta_phy

    def longitude_split(self):
        return np.pi * SPEED_OF_LIGHT / self.get_length()

    def system_scheme(self):
        fig, ax = plt.subplots()

        clctd_angle = 0
        center = np.array([0, 0])
        position = np.array([0, 0])

        for i in range(self.num_of_mirrors - 1):
            clctd_angle += self.elems[i].in_plane_angle
            center = (position + self.elems[i].radius
                      *np.array([np.cos(clctd_angle), np.sin(clctd_angle)]))
            clctd_angle += 180
            tmp_arc = patches.Arc(center, 2 * self.elems[i].radius,
                                  2 * self.elems[i].radius, 0,
                                  clctd_angle - MIRROR_ANGLE_SIZE,
                                  clctd_angle + MIRROR_ANGLE_SIZE,
                                  color='blue')
            ax.add_patch(tmp_arc)
            if i == 0:
                clctd_angle -= self.elems[i].in_plane_angle
                d_position = ((self.elems[i + 1].coord - self.elems[i].coord)
                              * np.array([1, 0]))
            else:
                clctd_angle += self.elems[i].in_plane_angle
                d_position = (-1 * (self.elems[i + 1].coord - self.elems[i].coord)
                              * np.array([np.cos(np.pi * clctd_angle / 180),
                                          np.sin(np.pi * clctd_angle / 180)]))
            tmp_line = mlines.Line2D([position[0], position[0] + d_position[0]],
                                     [position[1], position[1] + d_position[1]])
            ax.add_line(tmp_line)
            position = position + d_position
        
        clctd_angle += self.elems[-1].in_plane_angle
        center = (position + self.elems[-1].radius
                  *np.array([np.cos(clctd_angle), np.sin(clctd_angle)]))
        clctd_angle += 180
        tmp_arc = patches.Arc(center, 2 * self.elems[-1].radius,
                              2 * self.elems[-1].radius, 0,
                              clctd_angle - MIRROR_ANGLE_SIZE,
                              clctd_angle + MIRROR_ANGLE_SIZE,
                              color='blue')
        ax.add_patch(tmp_arc)
        plt.axis('equal')
        plt.axis('off')
        plt.show()

    def fund_lambda_choice(self, zero_lambda_approx):
        zero_omega_approx = 2 * np.pi * SPEED_OF_LIGHT / zero_lambda_approx
        n_closest = np.round((zero_omega_approx - self.transverse_split())
                             / self.longitude_split())
        return (2 * np.pi * SPEED_OF_LIGHT
                / (self.longitude_split() * n_closest + self.transverse_split()))

    def waist_search(self, zero_lambda_approx):
        using_lambda = self.fund_lambda_choice(zero_lambda_approx)
        mx_tangential = self.st_matrix_tangential()
        mx_sagittal = self.st_matrix_sagittal()
        init_waists2 = np.array([using_lambda / np.pi
                                * np.sqrt(-1 * mx_tangential[1, 1]
                                * mx_tangential[0, 1] / mx_tangential[1, 0]
                                / mx_tangential[0, 0]),
                                using_lambda / np.pi
                                * np.sqrt(-1 * mx_sagittal[1, 1]
                                * mx_sagittal[0, 1] / mx_sagittal[1, 0]
                                / mx_sagittal[0, 0])])
        init_curv_radius = np.array([-1, -1]) * self.elems[0].radius
        waist = np.array([[[0, 0, 0], [0, 0, 0]]])
        for i in range(self.num_of_mirrors - 1):
            beta_tangential = ((using_lambda * init_curv_radius[0] / np.pi
                               / init_waists2[0]) ** 2)
            beta_sagittal = ((using_lambda * init_curv_radius[1] / np.pi
                             / init_waists2[1]) ** 2)
            waist = np.concatenate((waist, np.array([[[- init_curv_radius[0]
                                    / (1 + beta_tangential),
                                    np.sqrt(beta_tangential)
                                    * np.sqrt(init_waists2[0])
                                    / np.sqrt(1 + beta_tangential), 0],
                                    [- init_curv_radius[1] / (1 + beta_sagittal),
                                    np.sqrt(beta_sagittal)
                                    * np.sqrt(init_waists2[1])
                                    / np.sqrt(1 + beta_sagittal), 0]]])))
            waist[i + 1, 0, 2] = using_lambda / np.pi / waist[i + 1, 0, 1]
            waist[i + 1, 1, 2] = using_lambda / np.pi / waist[i + 1, 1, 1]
            transform_mx_tangential = np.matrix([[1, (self.elems[i + 1].coord
                                                - self.elems[i].coord)],
                                                [0, 1]])
            transform_mx_tangential = np.matmul(self.elems[i + 1].get_matrix_tangential(),
                                                transform_mx_tangential)
            transform_mx_tangential = np.matmul(self.elems[i + 1].get_matrix_tangential(),
                                                transform_mx_tangential)
            res_waist2, res_radius = transrorm_waist(init_waists2[0],
                                                     init_curv_radius[0],
                                                     using_lambda,
                                                     transform_mx_tangential)
            init_waists2[0] = res_waist2
            init_curv_radius[0] = res_radius
            transform_mx_sagittal = np.matrix([[1, (self.elems[i + 1].coord
                                              - self.elems[i].coord)],
                                              [0, 1]])
            transform_mx_sagittal = np.matmul(self.elems[i + 1].get_matrix_sagittal(),
                                              transform_mx_sagittal)
            transform_mx_sagittal = np.matmul(self.elems[i + 1].get_matrix_sagittal(),
                                              transform_mx_sagittal)
            res_waist2, res_radius = transrorm_waist(init_waists2[1],
                                                     init_curv_radius[1],
                                                     using_lambda,
                                                     transform_mx_sagittal)
            init_waists2[1] = res_waist2
            init_curv_radius[1] = res_radius
        return waist[1:]    # [waist coord along axis, waist, Numerical Aperture]

    def realign(self):
        flag = True
        rotation = 0
        while not self.is_consistent():

            r_transform_mx_tangential = np.matrix([[1, 0], [0, 1]])
            for i in range(1, self.num_of_mirrors - 1):
                r_transform_mx_tangential = np.matmul(np.matrix([[1,
                                                      self.elems[i].coord
                                                      - self.elems[i - 1].coord],
                                                      [0, 1]]),
                                                      r_transform_mx_tangential)
                r_transform_mx_tangential = np.matmul(self.elems[i].get_matrix_tangential(),
                                                      r_transform_mx_tangential)
                r_transform_mx_tangential = np.matmul(self.elems[i].get_matrix_tangential(),
                                                      r_transform_mx_tangential)
            r_transform_mx_tangential = np.matmul(np.matrix([[1,
                                                  (self.elems[-1].coord
                                                  - self.elems[-2].coord)],
                                                  [0, 1]]),
                                                  r_transform_mx_tangential)

            r_transform_mx_sagittal = np.matrix([[1, 0], [0, 1]])
            for i in range(1, self.num_of_mirrors - 1):
                r_transform_mx_sagittal = np.matmul(np.matrix([[1,
                                                    self.elems[i].coord
                                                    - self.elems[i - 1].coord],
                                                    [0, 1]]),
                                                    r_transform_mx_sagittal)
                r_transform_mx_sagittal = np.matmul(self.elems[i].get_matrix_sagittal(),
                                                    r_transform_mx_sagittal)
                r_transform_mx_sagittal = np.matmul(self.elems[i].get_matrix_sagittal(),
                                                    r_transform_mx_sagittal)
            r_transform_mx_sagittal = np.matmul(np.matrix([[1,
                                                (self.elems[-1].coord
                                                - self.elems[-2].coord)],
                                                [0, 1]]),
                                                r_transform_mx_sagittal)

            central_start = self.elems[0].get_central_coord()
            z_central_start = central_start[0]
            x_central_start = central_start[1]
            y_central_start = central_start[2]

            self.elems[0].in_plane_angle = 0
            self.elems[0].coord = 0
            self.elems[0].in_plane_angle_deviation = 0
            self.elems[0].out_of_plane_angle_deviation = 0
            self.elems[0].in_plane_coord_deviation = 0
            self.elems[0].out_of_plane_coord_deviation = 0

            term_radius = -1 * self.elems[-1].radius

            x_central_term_tangential = (-1 * (-1)**self.num_of_mirrors
                                         * (self.elems[-1].radius
                                         * np.sin(self.elems[-1].in_plane_angle
                                         + self.elems[-1].in_plane_angle_deviation)
                                         + self.elems[-1].in_plane_coord_deviation)
                                         / (r_transform_mx_tangential[1, 0]
                                         *term_radius
                                         + r_transform_mx_tangential[0, 0]))
            new_radius = ((r_transform_mx_tangential[1, 1]*term_radius
                          + r_transform_mx_tangential[0, 1])
                          / (r_transform_mx_tangential[1, 0]*term_radius
                          + r_transform_mx_tangential[0, 0]))
            z_central_term_tangential = new_radius

            y_central_term_sagittal = (-1 * (-1)**self.num_of_mirrors
                                       * (self.elems[-1].radius
                                       *np.sin(self.elems[-1].out_of_plane_angle_deviation)
                                       + self.elems[-1].out_of_plane_coord_deviation)
                                       / (r_transform_mx_sagittal[1, 0]
                                       *term_radius
                                       + r_transform_mx_sagittal[0, 0]))
            new_radius = ((r_transform_mx_sagittal[1, 1]*term_radius
                          + r_transform_mx_sagittal[0, 1])
                          / (r_transform_mx_sagittal[1, 0]*term_radius
                          + r_transform_mx_sagittal[0, 0]))
            z_central_term_sagittal = new_radius

            direct = np.array([1, (x_central_start - x_central_term_tangential)
                              / np.abs(z_central_term_tangential - z_central_start),
                              (y_central_start - y_central_term_sagittal)
                              / np.abs(z_central_term_sagittal - z_central_start)])
            
            if flag:
                rotation = float(np.sqrt(direct[1]**2 + direct[2]**2))
                flag = False            
            
            direct = direct / np.sqrt(direct[0]**2
                                      + direct[1]**2
                                      + direct[2]**2)
            opt_path = self.elems[0].radius
            start_coord = np.array([z_central_start, x_central_start,
                                    y_central_start])

            z_solutions = np.array([0, 0])
            x_solutions = np.array([0, 0])
            y_solutions = np.array([0, 0])

            for i in range(1, self.num_of_mirrors):
                central_coord = self.elems[i].get_central_coord()

                if np.abs(direct[1]) < ZERO_ACCURACY:
                    if np.abs(direct[2]) < ZERO_ACCURACY:
                        x_solutions = np.array([start_coord[1], start_coord[1]])
                        y_solutions = np.array([start_coord[2], start_coord[2]])
                        if ((x_solutions[0] - central_coord[1])**2
                            + (y_solutions[0] - central_coord[2])**2
                            > self.elems[i].radius**2):
                            raise Exception("Misnp.sing mirror")
                        z_diff = np.sqrt(self.elems[i].radius**2
                                         - (x_solutions[0] - central_coord[1])**2
                                         - (y_solutions[0] - central_coord[2])**2)
                        z_solutions = np.array([central_coord[0] + z_diff,
                                                central_coord[0] - z_diff])
                    else:
                        x_solutions = np.array([start_coord[1], start_coord[1]])
                        A_coef = 1
                        B_coef = 2 * ((direct[2]*start_coord[0]
                                      - direct[0]*start_coord[2]
                                      - direct[2]*central_coord[0])*direct[0]
                                      - direct[2]**2 * central_coord[2])
                        C_coef = (direct[2]**2 * central_coord[2]**2
                                  + (direct[2]*start_coord[0]
                                  - direct[0]*start_coord[2]
                                  - direct[2]*central_coord[0])**2
                                  - direct[2]**2 * self.elems[i].radius**2
                                  + direct[2]**2
                                  * (start_coord[1] - central_coord[1])**2)
                        y_solutions = quadratic_solver(A_coef, B_coef, C_coef)
                        if isnan(y_solutions[0]):
                            raise Exception("Misnp.sing mirror")
                        z_solutions = (start_coord[0]*np.array([1, 1])
                                       + direct[0]/direct[2]
                                       *(y_solutions - start_coord[2]
                                       *np.array([1, 1])))
                else:
                    A_coef = 1
                    B_coef = 2 * ((direct[1]*start_coord[2]
                                  - direct[2]*start_coord[1]
                                  - direct[1]*central_coord[2])*direct[2]
                                  + (direct[1]*start_coord[0]
                                  - direct[0]*start_coord[1]
                                  - direct[1]*central_coord[0])*direct[0]
                                  - direct[1]**2 * central_coord[1])
                    C_coef = (direct[1]**2 * central_coord[1]**2
                              + (direct[1]*start_coord[2]
                              - direct[2]*start_coord[1]
                              - direct[1]*central_coord[2])**2
                              + (direct[1]*start_coord[0]
                              - direct[0]*start_coord[1]
                              - direct[1]*central_coord[0])**2
                              - direct[1]**2 * self.elems[i].radius**2)
                    x_solutions = quadratic_solver(A_coef, B_coef, C_coef)
                    if isnan(x_solutions[0]):
                        raise Exception("Misnp.sing mirror")
                    z_solutions = (start_coord[0]*np.array([1, 1])
                                   + direct[0]/direct[1]*(x_solutions
                                   - start_coord[1]*np.array([1, 1])))
                    y_solutions = (start_coord[2]*np.array([1, 1])
                                   + direct[2]/direct[1]*(x_solutions
                                   - start_coord[1]*np.array([1, 1])))
            
                cntrl_direct = np.array([[central_coord[0] - z_solutions[0],
                                         central_coord[1] - x_solutions[0],
                                         central_coord[2] - y_solutions[0]],
                                         [central_coord[0] - z_solutions[1],
                                         central_coord[1] - x_solutions[1],
                                         central_coord[2] - y_solutions[1]]])
                cntrl_direct[0] = cntrl_direct[0] / np.sqrt(cntrl_direct[0, 0]**2
                                                            + cntrl_direct[0, 1]**2
                                                            + cntrl_direct[0, 2]**2)
                cntrl_direct[1] = cntrl_direct[1] / np.sqrt(cntrl_direct[1, 0]**2
                                                            + cntrl_direct[1, 1]**2
                                                            + cntrl_direct[1, 2]**2)

                angle = 0

                if (direct[0]*cntrl_direct[0, 0] + direct[1]*cntrl_direct[0, 1]
                    + direct[2]*cntrl_direct[0, 2] < 0):
                    new_start_coord = np.array([z_solutions[0],
                                                x_solutions[0],
                                                y_solutions[0]])
                    cosx = (direct[0]*cntrl_direct[0, 0]
                            + direct[1]*cntrl_direct[0, 1]
                            + direct[2]*cntrl_direct[0, 2])
                    sinx = direct - cosx*cntrl_direct[0]
                    angle = np.arctan2(np.sqrt(sinx[0]**2 + sinx[1]**2 + sinx[2]**2),
                                       -1*cosx)
                    direct = direct - 2*cosx*cntrl_direct[0]
                else:
                    new_start_coord = np.array([z_solutions[1],
                                                x_solutions[1],
                                                y_solutions[1]])
                    cosx = (direct[0]*cntrl_direct[1, 0]
                            + direct[1]*cntrl_direct[1, 1]
                            + direct[2]*cntrl_direct[1, 2])
                    sinx = direct - cosx*cntrl_direct[1]
                    angle = np.arctan2(np.sqrt(sinx[0]**2
                                       + sinx[1]**2
                                       + sinx[2]**2),
                                       -1*cosx)
                    direct = direct - 2*cosx*cntrl_direct[1]

                d_start_coord = new_start_coord - start_coord
                opt_path += np.sqrt(d_start_coord[0]**2 + d_start_coord[1]**2
                                    + d_start_coord[2]**2)
                start_coord = new_start_coord
                self.elems[i].coord = opt_path
                self.elems[i].in_plane_angle = angle
                self.elems[i].in_plane_angle_deviation = 0
                self.elems[i].out_of_plane_angle_deviation = 0
                self.elems[i].in_plane_coord_deviation = 0
                self.elems[i].out_of_plane_coord_deviation = 0
            self.refresh()
        return rotation
    
    def waist_scheme(self, zero_lambda_approx, equal_axis=False):
        using_lambda = self.fund_lambda_choice(zero_lambda_approx)
        waist = self.waist_search(zero_lambda_approx)

        number_of_steps = 10000

        global number_of_plots_G
        number_of_plots_G += 1

        plt.figure(number_of_plots_G)
        ax = plt.gca()

        tmp_line = mlines.Line2D([0, 0], [1, -1])
        ax.add_line(tmp_line)
        length = self.get_length()
        tmp_line = mlines.Line2D([length, length], [1, -1])
        ax.add_line(tmp_line)
        tmp_line = mlines.Line2D([0, length], [0, 0], color='green')
        ax.add_line(tmp_line)
        tmp_line = mlines.Line2D([self.elems[0].radius,
                                 self.elems[0].radius - 0.05],
                                 [-1 / 6, -1 / 6], color='thistle')
        ax.add_line(tmp_line)
        tmp_line = mlines.Line2D([length - self.elems[-1].radius,
                                 length - self.elems[-1].radius + 0.05],
                                 [-1 / 6, -1 / 6], color='thistle')
        ax.add_line(tmp_line)
        tmp_line = mlines.Line2D([self.elems[0].radius - 0.05,
                                 self.elems[0].radius, self.elems[0].radius],
                                 [1 / 6, 1 / 6, -1 / 6], color='indigo')
        ax.add_line(tmp_line)
        tmp_line = mlines.Line2D([length - self.elems[-1].radius + 0.05,
                                 length - self.elems[-1].radius,
                                 length - self.elems[-1].radius],
                                 [1 / 6, 1 / 6, -1 / 6], color='indigo')
        ax.add_line(tmp_line)

        start = 0
        stop = self.elems[1].coord
        step = (stop - start) / (number_of_steps - 1)
        x_coords = np.arange(start, stop + step, step)

        x_1 = waist[0, 1, 0]
        tmp_line = mlines.Line2D([x_1, x_1], [1 / 12, -1 / 12],
                                 color='lightcoral')
        ax.add_line(tmp_line)
        z_R_1 = waist[0, 1, 1]**2 * np.pi / using_lambda
        upper_curv = np.empty(number_of_steps)
        lower_curv = np.empty(number_of_steps)
        for i in range(number_of_steps):
            upper_curv[i] = (waist[0, 1, 1]
                             * np.sqrt(1 + ((x_coords[i] - x_1) / z_R_1)**2))
        lower_curv = -upper_curv
        plt.plot(x_coords, upper_curv, color='lightcoral')
        plt.plot(x_coords, lower_curv, color='lightcoral')
        
        x_0 = waist[0, 0, 0]
        tmp_line = mlines.Line2D([x_0, x_0], [1 / 12 , -1 / 12], color='red')
        ax.add_line(tmp_line)
        z_R_0 = waist[0, 0, 1] ** 2 * np.pi / using_lambda
        for i in range(x_coords.size):
            upper_curv[i] = (waist[0, 0, 1]
                             * np.sqrt(1 + ((x_coords[i] - x_0) / z_R_0)**2))
        lower_curv = -upper_curv
        plt.plot(x_coords, upper_curv, color='r')
        plt.plot(x_coords, lower_curv, color='r')

        for i in range(1, self.num_of_mirrors - 1):
            x_init = self.elems[i].coord
            x_term = self.elems[i + 1].coord
            tmp_line = mlines.Line2D([x_init, x_init], [1 / 3, -1 / 3])
            ax.add_line(tmp_line)
            focal_length_tg = (self.elems[i].radius / 2
                               * np.cos(self.elems[i].in_plane_angle))
            focal_length_sag = (self.elems[i].radius / 2
                                / np.cos(self.elems[i].in_plane_angle))
            tmp_line = mlines.Line2D([x_init - focal_length_sag,
                                     x_init - focal_length_sag,
                                     x_init - focal_length_sag + 0.05],
                                     [1 / 6, -1 / 6, -1 / 6], color='thistle')
            ax.add_line(tmp_line)
            tmp_line = mlines.Line2D([x_init + focal_length_sag,
                                     x_init + focal_length_sag,
                                     x_init + focal_length_sag - 0.05],
                                     [1 / 6, -1 / 6, -1 / 6], color='thistle')
            ax.add_line(tmp_line)
            tmp_line = mlines.Line2D([x_init - focal_length_tg + 0.05,
                                     x_init - focal_length_tg,
                                     x_init - focal_length_tg],
                                     [1 / 6, 1 / 6, -1 / 6], color='indigo')
            ax.add_line(tmp_line)
            tmp_line = mlines.Line2D([x_init + focal_length_tg - 0.05,
                                     x_init + focal_length_tg,
                                     x_init + focal_length_tg],
                                     [1 / 6, 1 / 6, -1 / 6], color='indigo')
            ax.add_line(tmp_line)
            
            start = x_init
            stop = x_term
            step = (stop - start) / (number_of_steps - 1)
            x_coords = np.arange(start, stop + step, step)

            x_1 = waist[i, 1, 0]
            tmp_line = mlines.Line2D([x_1 + x_init, x_1 + x_init],
                                     [1 / 12, -1 / 12], color='lightcoral')
            ax.add_line(tmp_line)
            z_R_1 = waist[i, 1, 1]**2 * np.pi / using_lambda
            for j in range(x_coords.size):
                upper_curv[j] = (waist[i, 1, 1] * np.sqrt(1
                                 + ((x_coords[j] - x_1 - x_init) / z_R_1)**2))
            lower_curv = -upper_curv
            plt.plot(x_coords, upper_curv, color='lightcoral')
            plt.plot(x_coords, lower_curv, color='lightcoral')
        
            x_0 = waist[i, 0, 0]
            tmp_line = mlines.Line2D([x_0 + x_init, x_0 + x_init],
                                     [1 / 12, -1 / 12], color='red')
            ax.add_line(tmp_line)
            z_R_0 = waist[i, 0, 1] ** 2 * np.pi / using_lambda
            for j in range(x_coords.size):
                upper_curv[j] = (waist[i, 0, 1] * np.sqrt(1
                                 + ((x_coords[j] - x_0 - x_init) / z_R_0)**2))
            lower_curv = -upper_curv
            plt.plot(x_coords, upper_curv, color='r')
            plt.plot(x_coords, lower_curv, color='r')

                    
        plt.get_current_fig_manager().window.state('zoomed')
        if equal_axis:
            plt.axis('equal')
        plt.axis('off')
