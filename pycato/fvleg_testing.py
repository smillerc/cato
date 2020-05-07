# -*- coding: utf-8 -*-
import numpy as np

π = np.pi


def E0(rho, u, v, p, theta_ib, theta_ie, rho_tilde, a_tilde, rho_p_prime, p_p_prime):

    π = np.pi

    u_P = (1 / π) * np.sum(
        (-p / (rho_tilde * a_tilde)) * (np.sin(theta_ie) - np.sin(theta_ib))
        + u
        * (
            ((theta_ie - theta_ib) / 2)
            + ((np.sin(2 * theta_ie) - np.sin(2 * theta_ib)) / 4)
        )
        - v * ((np.cos(2 * theta_ie) - np.cos(2 * theta_ib)) / 4)
    )

    v_P = (1 / π) * np.sum(
        (p / (rho_tilde * a_tilde)) * (np.cos(theta_ie) - np.cos(theta_ib))
        - u * ((np.cos(2 * theta_ie) - np.cos(2 * theta_ib)) / 4)
        + v
        * (
            ((theta_ie - theta_ib) / 2)
            + ((np.sin(2 * theta_ie) - np.sin(2 * theta_ib)) / 4)
        )
    )

    p_P = (1 / (2 * π)) * np.sum(
        p * (theta_ie - theta_ib)
        - rho_tilde * a_tilde * u * (np.sin(theta_ie) - np.sin(theta_ib))
        + rho_tilde * a_tilde * v * (np.cos(theta_ie) - np.cos(theta_ib))
    )

    rho_P = rho_p_prime + (p_P / a_tilde ** 2) - (p_p_prime / a_tilde ** 2)

    return rho_P, u_P, v_P, p_P


def mach_cone(rho_i, u_i, v_i, p_i, plot=False):
    gamma = 5 / 3
    rho_i = np.array(rho_i)
    u_i = np.array(u_i)
    v_i = np.array(v_i)
    p_i = np.array(p_i)
    rho_tilde = np.sum(rho_i) / len(rho_i)
    u_tilde = np.sum(u_i) / len(u_i)
    v_tilde = np.sum(v_i) / len(v_i)
    p_tilde = np.sum(p_i) / len(p_i)
    a_tilde = np.sqrt(gamma * p_tilde / rho_tilde)
    mach_tilde = np.sqrt(u_tilde ** 2 + v_tilde ** 2) / a_tilde

    y1 = 1.0
    y2 = -1.0
    x_s = mach
    x1 = x_s
    x2 = x_s

    def intersect(x1, y1, x2, y2, r):
        ax = 0
        bx = 0
        ay = 0
        by = 0
        n_intersections = 0
        A = y2 - y1
        B = x1 - x2
        C = x1 * y2 - x2 * y1
        EPS = 1e-20
        if C ** 2 > (r ** 2 * (A ** 2 + B ** 2)):
            #             print('no intersections')
            n_intersections = 0

        elif abs(C ** 2 - r ** 2 * (A ** 2 + B ** 2)) < EPS:
            n_intersections = 1
            x0 = (A * C) / (A ** 2 + B ** 2)
            y0 = (B * C) / (A ** 2 + B ** 2)
            ax = x0
            ay = y0
        else:
            n_intersections = 2
            x0 = (A * C) / (A ** 2 + B ** 2)
            y0 = (B * C) / (A ** 2 + B ** 2)
            d = np.sqrt(r ** 2 - ((C ** 2) / (A ** 2 + B ** 2)))
            m = np.sqrt(d ** 2 / (A ** 2 + B ** 2))
            ax = x0 + B * m
            bx = x0 - B * m
            ay = y0 - A * m
            by = y0 + A * m

        return n_intersections, (ax, ay), (bx, by)

    r = 1
    theta = np.linspace(0, 2 * np.pi, 100)
    circle_x = r * np.cos(theta)
    circle_y = r * np.sin(theta)

    n_intersections, (ax, ay), (bx, by) = intersect(x1, y1, x2, y2, r)

    theta_ib = np.array([np.arctan2(ay, ax), 2 * np.pi + np.arctan2(by, bx)])
    theta_ie = np.array([np.arctan2(by, bx), np.arctan2(ay, ax)])
    dtheta = 2 * np.pi - np.abs(theta_ie - theta_ib)
    #     print(f'dtheta: {np.rad2deg(dtheta)}')

    if plot:
        plt.figure(figsize=(3, 3))

        plt.plot(ax, ay, "o", label="a")
        plt.plot(bx, by, "s", label="b")
        plt.plot([x1, x2], [y1, y2])
        plt.plot(0, 0, "s")
        plt.plot(circle_x, circle_y)
        plt.gca().set_aspect("equal")
        # plt.legend()
        plt.show()

    rho, u, v, p = E0(
        rho_i,
        u_i,
        v_i,
        p_i,
        theta_ib,
        theta_ie,
        rho_tilde,
        a_tilde,
        rho_p_prime,
        p_p_prime,
    )
    return rho, u, v, p, theta_ie, theta_ib


mach_cone(rho_i, u_i, v_i, p_i, plot=True)


plt.plot(mach, (theta_ie_c2 - theta_ib_c2) / 2, label="(theta_ie_c2 - theta_ib_c2)/2")
plt.plot(
    mach,
    np.sin(theta_ie_c2) - np.sin(theta_ib_c2),
    label="sin(theta_ie) - sin(theta_ib)",
)
plt.plot(
    mach,
    (np.sin(2 * theta_ie_c2) - np.sin(2 * theta_ib_c2)) / 4,
    label="(sin(2*theta_ie) - sin(2*theta_ib))/4",
)
plt.plot(
    mach,
    (np.cos(2 * theta_ie_c2) - np.cos(2 * theta_ib_c2)) / 4,
    label="(cos(2*theta_ie) - cos(2*theta_ib))/4",
)

plt.plot(mach, np.sin(theta_ie_c2), label="sin(theta_ie)")
plt.plot(mach, np.sin(theta_ib_c2), label="sin(theta_ib)")
plt.plot(mach, np.cos(theta_ie_c2), label="cos(theta_ie)")
plt.plot(mach, np.cos(theta_ib_c2), label="cos(theta_ib)")
plt.plot(mach, np.cos(2 * theta_ie_c2), label="cos(2theta_ie)")
plt.plot(mach, np.cos(2 * theta_ib_c2), label="cos(2theta_ib)")
plt.plot(mach, np.sin(2 * theta_ie_c2), label="sin(2theta_ie)")
plt.plot(mach, np.sin(2 * theta_ib_c2), label="sin(2theta_ib)")
