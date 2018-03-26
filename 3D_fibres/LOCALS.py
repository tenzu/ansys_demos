import math, random
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

r_disk = 37.5  # disk radius (mm)
td_ratio = 0.4
min_gap = 0.1  # minimum gap between blocks (mm)
f_num = 4  # number of fibres

# original segment center for cylindars (fibres)
r = 1.5  # spiral fibre rotation radius (mm)
s = 8  # spiral fibre pitch (mm)
n = 16  # division in single pitch (n equals to 2**?)
p = 1.75  # spiral fibre total pitch number (should equal to an integer plus integral multiple of 1/n)
theta = p * 360  # spiral fibre total rotation angle (DEG)
margin = 0.2  # margin from fibre to outer cylindar (mm)
r_cylindar = math.sqrt((r + margin) ** 2 + (
        s * p / 2) ** 2) * 1.2  # radius of sphere which covers whole cylindar (mm), 1.2 for additional margin
r_sph = 2 * r ** 2 * (1 - math.cos(math.pi / n))  # radius of spheres in cylindar (mm)
n_sph = int(n * p + 1)  # number of spheres in one cylindar
cld_ctr = np.zeros((f_num, 3))  # center of cylindars
cld_ctr_1 = np.zeros((f_num, 3))  # center of cylindars
sph_ctr = np.zeros((f_num * n_sph, 3))  # centers of spheres
rd = np.linspace(0, 2 * math.pi, n + 1)  # angles of sphere centers
local_cs = []  # local cordinate systems for fibres
b1_num = 5  # number of aggregate type 1
b1_ctr = np.zeros((b1_num, 4))  # aggregate 1 info (x, y, z, radius)
b2_num = 12  # number of aggregate type 2
b2_ctr = np.zeros((b2_num, 4))  # aggregate 2 info (x, y, z, radius)
b3_num = 32  # number of aggregate type 3
b3_ctr = np.zeros((b3_num, 4))  # aggregate 3 info (x, y, z, radius)
b4_num = 64  # number of aggregate type 4
b4_ctr = np.zeros((b4_num, 4))  # aggregate 4 info (x, y, z, radius)


# generate a random local cordinate system
def rdm_cs():
    rdm_cs = [0 for i in range(0, 6, 1)]  # random local cordinate systems (x,y,z,yz,zx,xy)
    for i in range(0, 6, 1):
        if i == 0:
            x_ = r_disk - 2 * r_cylindar - min_gap
            rdm_cs[i] = random.uniform(-x_, x_)
        elif i == 1:
            rdm_cs[i] = random.uniform(-math.sqrt(x_ ** 2 - rdm_cs[0] ** 2),
                                       math.sqrt(x_ ** 2 - rdm_cs[0] ** 2))
        elif i == 2:
            z_ = 2 * r_disk * td_ratio - 4 * r_cylindar - 2 * min_gap
            rdm_cs[i] = random.uniform(-z_ / 2, z_ / 2)
        else:
            rdm_cs[i] = random.uniform(0, 360)
    return (rdm_cs)


# generate original centers of spheres
def s_ctrs():
    s_ctr = np.zeros((n_sph, 3))  # sphere centers in cylindar (original zeros)
    for i in range(0, n_sph, 1):
        s_ctr[i][0] = r * math.cos(rd[i % n])
        s_ctr[i][1] = r * math.sin(rd[i % n])
        s_ctr[i][2] = s / n * i
    return (s_ctr)


# locations for fibres
def f_ctrs():
    f_cont = 0  # fibre counter
    while f_cont < f_num:
        f_avlb = 1
        tmp_cs = rdm_cs()  # generate a temporary local system
        theta1 = tmp_cs[3]
        theta2 = tmp_cs[4]
        theta3 = tmp_cs[5]
        T1_matrix = np.array(
            [[1, 0, 0], [0, math.cos(theta1), -math.sin(theta1)], [0, math.sin(theta1), math.cos(theta1)]])
        T2_matrix = np.array(
            [[math.cos(theta2), 0, math.sin(theta2)], [0, 1, 0], [-math.sin(theta2), 0, math.cos(theta2)]])
        T3_matrix = np.array(
            [[math.cos(theta3), -math.sin(theta3), 0], [math.sin(theta3), math.cos(theta3), 0], [0, 0, 1]])
        # before_rotation = np.array([0, 0, s * p / 2])
        before_rotation = np.array([0, 0, 0])
        before_rotation_1 = np.array([0, 0, s * p])
        after_rotation = np.dot(T1_matrix, np.dot(T2_matrix, np.dot(T3_matrix, before_rotation)))
        after_rotation_1 = np.dot(T1_matrix, np.dot(T2_matrix, np.dot(T3_matrix, before_rotation_1)))
        cld_ctr[f_cont, 0] = after_rotation[0] + tmp_cs[0]
        cld_ctr[f_cont, 1] = after_rotation[1] + tmp_cs[1]
        cld_ctr[f_cont, 2] = after_rotation[2] + tmp_cs[2]
        cld_ctr_1[f_cont, 2] = after_rotation_1[2] + tmp_cs[2]
        for i in range(0, f_cont, 1):
            if (cld_ctr[i, 0] - cld_ctr[f_cont, 0]) ** 2 + (cld_ctr[i, 1] - cld_ctr[f_cont, 1]) ** 2 + (
                    cld_ctr[i, 2] - cld_ctr[f_cont, 2]) ** 2 > (r_cylindar + r_cylindar) ** 2 + min_gap and abs(
                cld_ctr_1[f_cont, 2]) < abs(r_disk * td_ratio - r - min_gap):
                f_avlb = f_avlb * 1
            else:
                f_avlb = f_avlb * 0
        if f_avlb == 1:
            local_cs.append(tmp_cs)
            for j in range(0, n_sph, 1):
                before_rotation_2 = np.array([org_ctr[j, 0], org_ctr[j, 1], org_ctr[j, 2]])
                after_rotation_2 = np.dot(T1_matrix, np.dot(T2_matrix, np.dot(T3_matrix, before_rotation_2)))
                sph_ctr[j + f_cont * n_sph, 0] = after_rotation_2[0] + tmp_cs[0]
                sph_ctr[j + f_cont * n_sph, 1] = after_rotation_2[1] + tmp_cs[1]
                sph_ctr[j + f_cont * n_sph, 2] = after_rotation_2[2] + tmp_cs[2]
            f_cont += 1
        else:
            f_cont = f_cont


# Location for aggregate type 1
def b1_ctrs():
    b1_cont = 0  # aggregate counter
    while b1_cont < b1_num:
        b1_avlb = 1
        r_b1 = random.uniform(0.004e3, 0.005e3)
        x_ = r_disk - r_b1 - min_gap
        x_b1 = random.uniform(-x_, x_)
        y_b1 = random.uniform(-math.sqrt(x_ ** 2 - x_b1 ** 2), math.sqrt(x_ ** 2 - x_b1 ** 2))
        z_ = 2 * r_disk * td_ratio - 2 * r_b1 - 2 * min_gap
        z_b1 = random.uniform(-z_ / 2, z_ / 2)
        for i in range(0, int(f_num * (n * p + 1)), 1):
            if (sph_ctr[i, 0] - x_b1) ** 2 + (sph_ctr[i, 1] - y_b1) ** 2 + (sph_ctr[i, 2] - z_b1) ** 2 > (
                    r_sph + r_b1) ** 2 + min_gap:
                b1_avlb = b1_avlb * 1
            else:
                b1_avlb = b1_avlb * 0
        for i in range(0, b1_cont, 1):
            if (b1_ctr[i, 0] - x_b1) ** 2 + (b1_ctr[i, 1] - y_b1) ** 2 + (b1_ctr[i, 2] - z_b1) ** 2 > (
                    b1_ctr[i, 3] + r_b1) ** 2 + min_gap:
                b1_avlb = b1_avlb * 1
            else:
                b1_avlb = b1_avlb * 0
        if b1_avlb == 1:
            b1_ctr[b1_cont, 0] = x_b1
            b1_ctr[b1_cont, 1] = y_b1
            b1_ctr[b1_cont, 2] = z_b1
            b1_ctr[b1_cont, 3] = r_b1
            b1_cont += 1
        else:
            b1_cont = b1_cont


# Location for aggregate type 2
def b2_ctrs():
    b2_cont = 0  # aggregate counter
    while b2_cont < b2_num:
        b2_avlb = 1
        r_b2 = random.uniform(0.003e3, 0.004e3)
        x_ = r_disk - r_b2 - min_gap
        x_b2 = random.uniform(-x_, x_)
        y_b2 = random.uniform(-math.sqrt(x_ ** 2 - x_b2 ** 2), math.sqrt(x_ ** 2 - x_b2 ** 2))
        z_ = 2 * r_disk * td_ratio - 2 * r_b2 - 2 * min_gap
        z_b2 = random.uniform(-z_ / 2, z_ / 2)
        for i in range(0, int(f_num * (n * p + 1)), 1):
            if (sph_ctr[i, 0] - x_b2) ** 2 + (sph_ctr[i, 1] - y_b2) ** 2 + (sph_ctr[i, 2] - z_b2) ** 2 > (
                    r_sph + r_b2) ** 2 + min_gap:
                b2_avlb = b2_avlb * 1
            else:
                b2_avlb = b2_avlb * 0
        for i in range(0, b1_num, 1):
            if (b1_ctr[i, 0] - x_b2) ** 2 + (b1_ctr[i, 1] - y_b2) ** 2 + (b1_ctr[i, 2] - z_b2) ** 2 > (
                    b1_ctr[i, 3] + r_b2) ** 2 + min_gap:
                b2_avlb = b2_avlb * 1
            else:
                b2_avlb = b2_avlb * 0
        for i in range(0, b2_cont, 1):
            if (b2_ctr[i, 0] - x_b2) ** 2 + (b2_ctr[i, 1] - y_b2) ** 2 + (b2_ctr[i, 2] - z_b2) ** 2 > (
                    b2_ctr[i, 3] + r_b2) ** 2 + min_gap:
                b2_avlb = b2_avlb * 1
            else:
                b2_avlb = b2_avlb * 0
        if b2_avlb == 1:
            b2_ctr[b2_cont, 0] = x_b2
            b2_ctr[b2_cont, 1] = y_b2
            b2_ctr[b2_cont, 2] = z_b2
            b2_ctr[b2_cont, 3] = r_b2
            b2_cont += 1
        else:
            b2_cont = b2_cont


# Location for aggregate type 3
def b3_ctrs():
    b3_cont = 0  # aggregate counter
    while b3_cont < b3_num:
        b3_avlb = 1
        r_b3 = random.uniform(0.002e3, 0.003e3)
        x_ = r_disk - r_b3 - min_gap
        x_b3 = random.uniform(-x_, x_)
        y_b3 = random.uniform(-math.sqrt(x_ ** 2 - x_b3 ** 2), math.sqrt(x_ ** 2 - x_b3 ** 2))
        z_ = 2 * r_disk * td_ratio - 2 * r_b3 - 2 * min_gap
        z_b3 = random.uniform(-z_ / 2, z_ / 2)
        for i in range(0, int(f_num * (n * p + 1)), 1):
            if (sph_ctr[i, 0] - x_b3) ** 2 + (sph_ctr[i, 1] - y_b3) ** 2 + (sph_ctr[i, 2] - z_b3) ** 2 > (
                    r_sph + r_b3) ** 2 + min_gap:
                b3_avlb = b3_avlb * 1
            else:
                b3_avlb = b3_avlb * 0
        for i in range(0, b1_num, 1):
            if (b1_ctr[i, 0] - x_b3) ** 2 + (b1_ctr[i, 1] - y_b3) ** 2 + (b1_ctr[i, 2] - z_b3) ** 2 > (
                    b1_ctr[i, 3] + r_b3) ** 2 + min_gap:
                b3_avlb = b3_avlb * 1
            else:
                b3_avlb = b3_avlb * 0
        for i in range(0, b2_num, 1):
            if (b2_ctr[i, 0] - x_b3) ** 2 + (b2_ctr[i, 1] - y_b3) ** 2 + (b2_ctr[i, 2] - z_b3) ** 2 > (
                    b2_ctr[i, 3] + r_b3) ** 2 + min_gap:
                b3_avlb = b3_avlb * 1
            else:
                b3_avlb = b3_avlb * 0
        for i in range(0, b3_cont, 1):
            if (b3_ctr[i, 0] - x_b3) ** 2 + (b3_ctr[i, 1] - y_b3) ** 2 + (b3_ctr[i, 2] - z_b3) ** 2 > (
                    b3_ctr[i, 3] + r_b3) ** 2 + min_gap:
                b3_avlb = b3_avlb * 1
            else:
                b3_avlb = b3_avlb * 0
        if b3_avlb == 1:
            b3_ctr[b3_cont, 0] = x_b3
            b3_ctr[b3_cont, 1] = y_b3
            b3_ctr[b3_cont, 2] = z_b3
            b3_ctr[b3_cont, 3] = r_b3
            b3_cont += 1
        else:
            b3_cont = b3_cont


# Location for aggregate type 4
def b4_ctrs():
    b4_cont = 0  # aggregate counter
    while b4_cont < b4_num:
        b4_avlb = 1
        r_b4 = random.uniform(0.001e3, 0.002e3)
        x_ = r_disk - r_b4 - min_gap
        x_b4 = random.uniform(-x_, x_)
        y_b4 = random.uniform(-math.sqrt(x_ ** 2 - x_b4 ** 2), math.sqrt(x_ ** 2 - x_b4 ** 2))
        z_ = 2 * r_disk * td_ratio - 2 * r_b4 - 2 * min_gap
        z_b4 = random.uniform(-z_ / 2, z_ / 2)
        for i in range(0, int(f_num * (n * p + 1)), 1):
            if (sph_ctr[i, 0] - x_b4) ** 2 + (sph_ctr[i, 1] - y_b4) ** 2 + (sph_ctr[i, 2] - z_b4) ** 2 > (
                    r_sph + r_b4) ** 2 + min_gap:
                b4_avlb = b4_avlb * 1
            else:
                b4_avlb = b4_avlb * 0
        for i in range(0, b1_num):
            if (b1_ctr[i, 0] - x_b4) ** 2 + (b1_ctr[i, 1] - y_b4) ** 2 + (b1_ctr[i, 2] - z_b4) ** 2 > (
                    b1_ctr[i, 3] + r_b4) ** 2 + min_gap:
                b4_avlb = b4_avlb * 1
            else:
                b4_avlb = b4_avlb * 0
        for i in range(0, b2_num, 1):
            if (b2_ctr[i, 0] - x_b4) ** 2 + (b2_ctr[i, 1] - y_b4) ** 2 + (b2_ctr[i, 2] - z_b4) ** 2 > (
                    b2_ctr[i, 3] + r_b4) ** 2 + min_gap:
                b4_avlb = b4_avlb * 1
            else:
                b4_avlb = b4_avlb * 0
        for i in range(0, b3_num, 1):
            if (b3_ctr[i, 0] - x_b4) ** 2 + (b3_ctr[i, 1] - y_b4) ** 2 + (b3_ctr[i, 2] - z_b4) ** 2 > (
                    b3_ctr[i, 3] + r_b4) ** 2 + min_gap:
                b4_avlb = b4_avlb * 1
            else:
                b4_avlb = b4_avlb * 0
        for i in range(0, b4_cont, 1):
            if (b4_ctr[i, 0] - x_b4) ** 2 + (b4_ctr[i, 1] - y_b4) ** 2 + (b4_ctr[i, 2] - z_b4) ** 2 > (
                    b4_ctr[i, 3] + r_b4) ** 2 + min_gap:
                b4_avlb = b4_avlb * 1
            else:
                b4_avlb = b4_avlb * 0
        if b4_avlb == 1:
            b4_ctr[b4_cont, 0] = x_b4
            b4_ctr[b4_cont, 1] = y_b4
            b4_ctr[b4_cont, 2] = z_b4
            b4_ctr[b4_cont, 3] = r_b4
            b4_cont += 1
        else:
            b4_cont = b4_cont


org_ctr = s_ctrs()
f_ctrs()
b1_ctrs()
b2_ctrs()
b3_ctrs()
b4_ctrs()

# Record fibre locations
f1 = open('sph_ctrs.txt', 'w')
for i in range(len(sph_ctr)):
    f1.write('%11.5f' % sph_ctr[i][0] + ',' + '%11.5f' % sph_ctr[i][1] + ',' + '%11.5f' % sph_ctr[i][2] + '\n')
f1.close()
f2 = open('CS_trans.txt', 'w')
for i in range(len(local_cs)):
    f2.write('%11.5f' % local_cs[i][0] + ',' + '%11.5f' % local_cs[i][1] + ',' + '%11.5f' % local_cs[i][2] + '\n')
f2.close()
f3 = open('CS_rotat.txt', 'w')
for i in range(len(local_cs)):
    f3.write('%11.5f' % local_cs[i][3] + ',' + '%11.5f' % local_cs[i][4] + ',' + '%11.5f' % local_cs[i][5] + '\n')
f3.close()
# Record b1 location
f1 = open('b1_ctrs.txt', 'w')
for i in range(0, len(b1_ctr), 1):
    f1.write(
        '%11.5f' % b1_ctr[i, 0] + ' ' + '%11.5f' % b1_ctr[i, 1] + ' ' + '%11.5f' % b1_ctr[i, 2] + '%11.5f' % b1_ctr[
            i, 3] + '\n')
f1.close()
# Record b2 location
f1 = open('b2_ctrs.txt', 'w')
for i in range(0, len(b2_ctr), 1):
    f1.write(
        '%11.5f' % b2_ctr[i, 0] + ' ' + '%11.5f' % b2_ctr[i, 1] + ' ' + '%11.5f' % b2_ctr[i, 2] + '%11.5f' % b2_ctr[
            i, 3] + '\n')
f1.close()
# Record b3 location
f1 = open('b3_ctrs.txt', 'w')
for i in range(0, len(b3_ctr), 1):
    f1.write(
        '%11.5f' % b3_ctr[i, 0] + ' ' + '%11.5f' % b3_ctr[i, 1] + ' ' + '%11.5f' % b3_ctr[i, 2] + '%11.5f' % b3_ctr[
            i, 3] + '\n')
f1.close()
# Record b4 location
f1 = open('b4_ctrs.txt', 'w')
for i in range(0, len(b4_ctr), 1):
    f1.write(
        '%11.5f' % b4_ctr[i, 0] + ' ' + '%11.5f' % b4_ctr[i, 1] + ' ' + '%11.5f' % b4_ctr[i, 2] + '%11.5f' % b4_ctr[
            i, 3] + '\n')
f1.close()


# plot spheres in matplotlib
def plt_all():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim([-r_disk, r_disk])
    ax.set_ylim([-r_disk, r_disk])
    ax.set_zlim([-(2 * r_disk * td_ratio) / 2, (2 * r_disk * td_ratio) / 2])  # for plotting unstreched model
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    sphere_color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k']
    for i in range(0, f_num, 1):
        for j in range(i * len(org_ctr), (i + 1) * len(org_ctr), 1):
            ax.scatter(sph_ctr[j][0], sph_ctr[j][1], sph_ctr[j][2], c=sphere_color[i])  # draw spheres (scatters)
    for i in range(0, len(b1_ctr), 1):
        ax.scatter(b1_ctr[i][0], b1_ctr[i][1], b1_ctr[i][2], c='b')  # plot aggregate 1 (scatters)
    for i in range(0, len(b2_ctr), 1):
        ax.scatter(b2_ctr[i][0], b2_ctr[i][1], b2_ctr[i][2], c='g')  # plot aggregate 2 (scatters)
    for i in range(0, len(b3_ctr), 1):
        ax.scatter(b3_ctr[i][0], b3_ctr[i][1], b3_ctr[i][2], c='r')  # plot aggregate 3 (scatters)
    for i in range(0, len(b4_ctr), 1):
        ax.scatter(b4_ctr[i][0], b4_ctr[i][1], b4_ctr[i][2], c='c')  # plot aggregate 3 (scatters)
    plt.show()


plt_all()
