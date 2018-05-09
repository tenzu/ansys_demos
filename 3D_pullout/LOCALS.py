import math, random
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

r_disk = 37.5  # disk radius (mm)
td_ratio = 0.4
min_gap = 0.1  # minimum gap between blocks (mm)
f_num = 1  # number of fibres

# original segment center for cylindars (fibres)
r = 1.5  # spiral fibre rotation radius (mm)
s = 5  # spiral fibre pitch (mm)
n = 16  # division in single pitch (n equals to 2**?)
p = 2.75  # spiral fibre total pitch number (should equal to an integer plus integral multiple of 1/n)
theta = p * 360  # spiral fibre total rotation angle (DEG)
margin = 0.2  # margin from fibre to outer cylindar (mm)
r_cylindar = math.sqrt((r + margin) ** 2 + (
        s * p / 2) ** 2)  # radius of sphere which covers whole cylindar (mm)
r_sph = 2 * r ** 2 * (1 - math.cos(math.pi / n))  # radius of spheres in cylindar (mm)
n_sph = int(n * p + 1)  # number of spheres in one cylindar
cld_ctr = np.zeros((f_num, 3))  # center of cylindars
cld_ctr_1 = np.zeros((f_num, 3))  # center of cylindars
cld_ctr_2 = np.zeros((f_num, 3))  # center of cylindars
sph_ctr = np.zeros((f_num * n_sph, 3))  # centers of spheres
rd = np.linspace(0, 2 * math.pi, n + 1)  # angles of sphere centers
local_cs = []  # local cordinate systems for fibres



# generate a random local cordinate system
def rdm_cs():
    rdm_cs = [0 for i in range(0, 6, 1)]  # random local cordinate systems (x,y,z,yz,zx,xy)
    for i in range(0, 6, 1):
        if i == 0:
            x_ = 0
            rdm_cs[i] = random.uniform(-x_, x_)
        elif i == 1:
            rdm_cs[i] = random.uniform(0, 0)
        elif i == 2:
            z_ = 0
            rdm_cs[i] = random.uniform(-z_ / 2, z_ / 2)
        else:
            rdm_cs[i] = random.uniform(0, 0)
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
        # row 1 for user CS
        # row 2 for overlap judgement of cylindar in cylindar axis direction
        # row 3 for overlap judgement of cylindar in disk thickness
        before_rotation = np.array([[0, 0, 0], [0, 0, s * p / 2], [0, 0, s * p]])
        after_rotation = np.dot(T1_matrix, np.dot(T2_matrix, np.dot(T3_matrix, before_rotation)))
        cld_ctr[f_cont, 0] = after_rotation[0][0] + tmp_cs[0]
        cld_ctr[f_cont, 1] = after_rotation[0][1] + tmp_cs[1]
        cld_ctr[f_cont, 2] = after_rotation[0][2] + tmp_cs[2]
        cld_ctr_1[f_cont, 2] = after_rotation[2][2] + tmp_cs[2]
        cld_ctr_2[f_cont, 0] = after_rotation[1][0] + tmp_cs[0]
        cld_ctr_2[f_cont, 1] = after_rotation[1][1] + tmp_cs[1]
        cld_ctr_2[f_cont, 2] = after_rotation[1][2] + tmp_cs[2]
        for i in range(0, f_cont, 1):
            if (cld_ctr[i, 0] - cld_ctr[f_cont, 0]) ** 2 + (cld_ctr[i, 1] - cld_ctr[f_cont, 1]) ** 2 + (
                    cld_ctr[i, 2] - cld_ctr[f_cont, 2]) ** 2 > (r_cylindar + r_cylindar) ** 2 + min_gap and abs(
                cld_ctr_1[f_cont, 2]) < abs(r_disk * td_ratio - r - min_gap) and (
                    cld_ctr_2[i, 0] - cld_ctr_2[f_cont, 0]) ** 2 + (cld_ctr_2[i, 1] - cld_ctr_2[f_cont, 1]) ** 2 + (
                    cld_ctr_2[i, 2] - cld_ctr_2[f_cont, 2]) ** 2 > (r_cylindar + r_cylindar) ** 2 + min_gap:
                f_avlb = f_avlb * 1
            else:
                f_avlb = f_avlb * 0
        if f_avlb == 1:
            local_cs.append(tmp_cs)
            for j in range(0, n_sph, 1):
                sph_ctr_before_rotation = np.array([org_ctr[j, 0], org_ctr[j, 1], org_ctr[j, 2]])
                sph_ctr_after_rotation = np.dot(T1_matrix,
                                                np.dot(T2_matrix, np.dot(T3_matrix, sph_ctr_before_rotation)))
                sph_ctr[j + f_cont * n_sph, 0] = sph_ctr_after_rotation[0] + tmp_cs[0]
                sph_ctr[j + f_cont * n_sph, 1] = sph_ctr_after_rotation[1] + tmp_cs[1]
                sph_ctr[j + f_cont * n_sph, 2] = sph_ctr_after_rotation[2] + tmp_cs[2]
            f_cont += 1
        else:
            f_cont = f_cont


org_ctr = s_ctrs()
f_ctrs()

# Record fibre locations
f1 = open('sph_ctrs.txt', 'w')
for i in range(len(sph_ctr)):
    f1.write('%11.5f' % sph_ctr[i][0] + ',' + '%11.5f' % sph_ctr[i][1] + ',' + '%11.5f' % sph_ctr[i][2] + '\n')
f1.close()
f1 = open('CS_trans.txt', 'w')
for i in range(len(local_cs)):
    f1.write('%11.5f' % local_cs[i][0] + ',' + '%11.5f' % local_cs[i][1] + ',' + '%11.5f' % local_cs[i][2] + '\n')
f1.close()
f1 = open('CS_rotat.txt', 'w')
for i in range(len(local_cs)):
    f1.write('%11.5f' % local_cs[i][3] + ',' + '%11.5f' % local_cs[i][4] + ',' + '%11.5f' % local_cs[i][5] + '\n')
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
    plt.show()


plt_all()
