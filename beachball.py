import numpy as np
from matplotlib import pyplot as plt

from focal_mechanism import plane2mt, plane2nd, plane2tbp, mt2plane, fibonacci_sphere


class Beachball(object):
    def __init__(
        self, strike=None, dip=None, rake=None, mt=None, pure_dp=False, N=2000
    ):
        self.strike = strike
        self.dip = dip
        self.rake = rake
        self.mt = mt
        self.pure_dp = pure_dp
        self.N = N

        self.p = None
        self.b = None
        self.t = None
        self.z = None
        self.y = None
        self.x = None
        self.pl1cross_ball = None
        self.pl2cross_ball = None
        self.tension_sphere = None

    def cal_plane_sphere(self):
        if self.mt is None:
            self.mt = plane2mt(1, self.strike, self.dip, self.rake)
            n1, n2 = plane2nd(self.strike, self.dip, self.rake)
            self.t, self.b, self.p = plane2tbp(self.strike, self.dip, self.rake)
        else:
            [_, _, n1, _, n2, _, self.t, self.b, self.p, _] = mt2plane(self.mt)

        sphere = fibonacci_sphere(r0=1, N=self.N)

        index_pl1 = np.abs(np.dot(sphere, np.array([n1]).T)) <= 10 ** (-2)
        self.pl1cross_ball = sphere[index_pl1[:, 0]]

        index_pl2 = np.abs(np.dot(sphere, np.array([n2]).T)) <= 10 ** (-2)
        self.pl2cross_ball = sphere[index_pl2[:, 0]]

        if self.pure_dp:
            index_tension1 = (np.dot(sphere, np.array([n1]).T) >= 0) & (
                np.dot(sphere, np.array([n2]).T) >= 0
            )
            index_tension2 = (np.dot(sphere, np.array([n1]).T) <= 0) & (
                np.dot(sphere, np.array([n2]).T) <= 0
            )
            self.tension_sphere = sphere[index_tension1[:, 0] | index_tension2[:, 0]]
        else:
            M = np.array(
                [
                    [self.mt[0], self.mt[1], self.mt[2]],
                    [self.mt[1], self.mt[3], self.mt[4]],
                    [self.mt[2], self.mt[4], self.mt[5]],
                ]
            )

            index_tension = []
            for i in range(len(sphere)):
                if np.dot(sphere[i], np.dot(M, sphere[i].T)) > 0:
                    index_tension.append(True)
                else:
                    index_tension.append(False)
            # print(index_tension)
            self.tension_sphere = sphere[np.array(index_tension)]

    @staticmethod
    def wulff_projection(points):
        points = points[points[:, 2] >= 0]
        x_ = points[:, 0] * 2 / (points[:, 2] - (-1)) / 2
        y_ = points[:, 1] * 2 / (points[:, 2] - (-1)) / 2
        return [x_, y_]

    @staticmethod
    def wulff_projection_point(point):
        if np.abs(point[2] - (-1)) <= 10 ** (-4):
            return [0, 0]
        if point[2] < 0:
            point = -point
        x_ = point[0] * 2 / (point[2] - (-1)) / 2
        y_ = point[1] * 2 / (point[2] - (-1)) / 2
        return [x_, y_]

    def plot_beachball_3d(self):
        if self.tension_sphere is None:
            self.cal_plane_sphere()
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        # 画图时转化至ENU坐标系
        ax.set_xlabel("x (E)")
        ax.set_ylabel("y (N)")
        ax.set_zlabel("z (U)")

        theta_s = np.linspace(0, np.pi, 100, endpoint=True)
        phi_s = np.linspace(0, 2 * np.pi, 100, endpoint=True)
        theta_s, phi_s = np.meshgrid(theta_s, phi_s)
        x_s = np.sin(theta_s) * np.cos(phi_s)
        y_s = np.sin(theta_s) * np.sin(phi_s)
        z_s = np.cos(theta_s)

        ax.plot_surface(
            y_s, x_s, -z_s, antialiased=True, alpha=0.5, color="gray", shade=False
        )
        ax.scatter3D(
            self.pl1cross_ball[:, 1],
            self.pl1cross_ball[:, 0],
            -self.pl1cross_ball[:, 2],
            color="black",
        )
        ax.scatter3D(
            self.pl2cross_ball[:, 1],
            self.pl2cross_ball[:, 0],
            -self.pl2cross_ball[:, 2],
            color="black",
        )
        print(self.pl1cross_ball.shape)
        print(self.pl2cross_ball.shape)
        ax.scatter3D(
            self.tension_sphere[:, 1],
            self.tension_sphere[:, 0],
            -self.tension_sphere[:, 2],
            alpha=0.25,
            color="blue",
        )

        ax.scatter3D(0, 0, 0, color="black")
        ax.text3D(0, 0, 0, "O", color="black", fontsize=20)

        ax.quiver(0, 0, 0, self.p[1], self.p[0], -self.p[2], color="red")
        ax.text3D(
            self.p[1] * 1.2,
            self.p[0] * 1.2,
            -self.p[2] * 1.2,
            "P",
            color="red",
            fontsize=20,
        )

        ax.quiver(0, 0, 0, self.b[1], self.b[0], -self.b[2], color="black")
        ax.text3D(
            self.b[1] * 1.2,
            self.b[0] * 1.2,
            -self.b[2] * 1.2,
            "B",
            color="black",
            fontsize=20,
        )

        ax.quiver(0, 0, 0, self.t[1], self.t[0], -self.t[2], color="green")
        ax.text3D(
            self.t[1] * 1.2,
            self.t[0] * 1.2,
            -self.t[2] * 1.2,
            "T",
            color="green",
            fontsize=20,
        )
        ax.view_init(elev=30, azim=45)
        plt.show()

    def cal_wulff_points(self):
        if self.tension_sphere is None:
            self.cal_plane_sphere()
        theta_p = np.linspace(0, 2 * np.pi, self.N)
        big_circle = [np.cos(theta_p), np.sin(theta_p)]
        pl1_wulff_projection = self.wulff_projection(self.pl1cross_ball)
        pl2_wulff_projection = self.wulff_projection(self.pl2cross_ball)
        tension_sphere_wulff_projection = self.wulff_projection(self.tension_sphere)

        p_p = self.wulff_projection_point(self.p)
        b_p = self.wulff_projection_point(self.b)
        t_p = self.wulff_projection_point(self.t)
        return [
            pl1_wulff_projection,
            pl2_wulff_projection,
            tension_sphere_wulff_projection,
            p_p,
            b_p,
            t_p,
            big_circle,
        ]

    def plot_beachball_2d(self):
        [
            pl1_wulff_projection,
            pl2_wulff_projection,
            tension_sphere_wulff_projection,
            p_p,
            b_p,
            t_p,
            big_circle,
        ] = self.cal_wulff_points()

        fig = plt.figure()
        ax = fig.add_subplot()
        # 画图时转化至东北EN坐标系
        ax.scatter(big_circle[1], big_circle[0], color="black", s=0.5)
        ax.scatter(
            pl1_wulff_projection[1], pl1_wulff_projection[0], color="black", s=0.5
        )
        ax.scatter(
            pl2_wulff_projection[1], pl2_wulff_projection[0], color="black", s=0.5
        )
        ax.scatter(
            tension_sphere_wulff_projection[1],
            tension_sphere_wulff_projection[0],
            marker="+",
            color="gray",
        )
        ax.text(p_p[1], p_p[0], "P", color="red", fontsize=20)
        ax.text(b_p[1], b_p[0], "B", color="black", fontsize=20)
        ax.text(t_p[1], t_p[0], "T", color="green", fontsize=20)
        ax.axis("equal")

        for i in range(36):
            s = i * 10
            ax.text(
                np.sin(i / 36 * 2 * np.pi) * 1.2,
                np.cos(i / 36 * 2 * np.pi) * 1.2,
                str(s),
            )
        ax.axis("off")
        plt.show()
    
    @staticmethod
    def map_to_canvas(x, y, min_x, max_x, min_y, max_y, width, height):
        col = round((x - min_x) / (max_x - min_x) * (width - 1))
        row = round((max_y - y) / (max_y - min_y) * (height - 1))
        return row, col

    def plot_beachball_ascii(self, width=50, height=16):
        [
            pl1_wulff_projection,
            pl2_wulff_projection,
            tension_sphere_wulff_projection,
            p_p,
            b_p,
            t_p,
            big_circle,
        ] = self.cal_wulff_points()

        x_circle = np.array(big_circle[1])
        y_circle = np.array(big_circle[0])
        x_tension = np.array(tension_sphere_wulff_projection[1])
        y_tension = np.array(tension_sphere_wulff_projection[0])
        
        
        margin = 0
        x_all = np.concatenate([x_circle, x_tension])
        y_all = np.concatenate([y_circle, y_tension])
        min_x, max_x = x_all.min() - margin, x_all.max() + margin
        min_y, max_y = y_all.min() - margin, y_all.max() + margin

        canvas = [[" " for _ in range(width)] for _ in range(height)]

        for x, y in zip(x_circle, y_circle):
            row, col = self.map_to_canvas(x, y, min_x, max_x, min_y, max_y, width, height)
            if 0 <= row < height and 0 <= col < width:
                canvas[row][col] = "."

        for x, y in zip(x_tension, y_tension):
            row, col = self.map_to_canvas(x, y, min_x, max_x, min_y, max_y, width, height)
            if 0 <= row < height and 0 <= col < width:
                canvas[row][col] = "+"

        for row in canvas:
            print("".join(row))

if __name__ == "__main__":
    import sys
    paras = [float(item) for item in sys.argv[1:]]
    bb=Beachball(*paras, N=10000)
    bb.plot_beachball_ascii()
