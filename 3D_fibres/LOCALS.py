import math
import random

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


# generate an odd number (to make sure L_sphere is an odd number)
def odd_int(x):
    if x % 2 == 0:
        x += 1
    return (x)


r_disk = 37.5  # disk radius (mm)
td_ratio = 0.4
thickness = 2 * r_disk * td_ratio  # disk thickness (mm)
min_gap = 0.1  # minimum gap between blocks (mm)
f_num = 2  # number of fibres

# original segment center for cylindars (fibres)
r = 2.5  # spiral fibre rotation radius (mm)
s = 5  # spiral fibre pitch (mm)
n = 8  # division in single pitch (n equals to 2**?)
p = 1.5  # spiral fibre total pitch number (should equal to an integer plus integral multiple of 1/n)
theta = p * 360  # spiral fibre total rotation angle (DEG)
bargin = r * 0.1  # bargin from fibre to outer cylindar (mm)
r_cylindar = math.sqrt((r / 2) ** 2 + (s * p / 2) ** 2)  # radius of sphere which covers whole cylindar (mm)
r_sphere = 2 * r ** 2 * (1 - math.cos(math.pi / n))  # radius of spheres in cylindar (mm)
L_sphere = math.ceil((s * p + 2 * r_sphere) / r_sphere)  # number of spheres in cylindar longitudinal direction
n_sphs = n * odd_int(L_sphere)  # number of spheres in one cylindar
sph_ctrs = np.zeros((f_num * n_sphs, 3))  # centers of spheres
rd = np.linspace(0, 2 * math.pi, n + 1)  # angles of sphere centers
local_CSs = []  # local cordinate systems for fibres


# generate a random local cordinate system
def rdm_CS():
    rdm_CS = [0 for i in range(0, 6, 1)]  # random local cordinate systems (x,y,z,yz,zx,xy)
    for i in range(0, 6, 1):
        if i == 0:
            # r~=r_disk-r_cylindar-min_gap, x between (-r~,r~)
            r_ = r_disk - r_cylindar - min_gap
            rdm_CS[i] = random.uniform(-r_, r_)
        elif i == 1:
            rdm_CS[i] = random.uniform(-math.sqrt(r_ ** 2 - rdm_CS[0] ** 2),
                                       math.sqrt(r_ ** 2 - rdm_CS[0] ** 2))
        elif i == 2:
            # t~=2*r_disk*td_ratio, z between (-(t~/2-s*p/2),t~/2-s*p/2)
            t_ = 2 * r_disk * td_ratio - min_gap
            rdm_CS[i] = random.uniform(-(t_ / 3 - s * p / 2), t_ / 2 - s * p / 2)
        else:
            rdm_CS[i] = random.uniform(0, 360)
    return (rdm_CS)


# generate original centers of spheres
def s_ctrs():
    s_ctrs = np.zeros((n_sphs, 3))  # sphere centers in cylindar (original zeros)
    # generate list of [0,1,-1,2,-2,3,-3,...]
    tmp_list = [i for i in range(-(n_sphs // n // 2), (n_sphs // n // 2) + 1, 1)]
    for i in range(0, n_sphs, 1):
        s_ctrs[i][0] = (r + r_sphere) * math.cos(rd[i % n])
        s_ctrs[i][1] = (r + r_sphere) * math.sin(rd[i % n])
        s_ctrs[i][2] = s * p / L_sphere * tmp_list[i // n]
    return (s_ctrs)


# locations for fibres
org_ctrs = s_ctrs()
f_cont = 0  # fibre counter
while f_cont < f_num:
    f_avlb = 1
    tmp_CS = rdm_CS()  # generate a temporary local system
    for i in range(0, n_sphs, 1):
        theta1 = tmp_CS[3]
        theta2 = tmp_CS[4]
        theta3 = tmp_CS[5]
        T1_matrix = np.array(
            [[1, 0, 0], [0, math.cos(theta1), -math.sin(theta1)], [0, math.sin(theta1), math.cos(theta1)]])
        T2_matrix = np.array(
            [[math.cos(theta2), 0, math.sin(theta2)], [0, 1, 0], [-math.sin(theta2), 0, math.cos(theta2)]])
        T3_matrix = np.array(
            [[math.cos(theta3), -math.sin(theta3), 0], [math.sin(theta3), math.cos(theta3), 0], [0, 0, 1]])
        before_rotation = np.array([org_ctrs[i, 0], org_ctrs[i, 1], org_ctrs[i, 2]])
        after_rotation = np.dot(T1_matrix, np.dot(T2_matrix, np.dot(T3_matrix, before_rotation)))
        sph_ctrs[i + (f_cont) * (n_sphs), 0] = after_rotation[0] + tmp_CS[0]
        sph_ctrs[i + (f_cont) * (n_sphs), 1] = after_rotation[1] + tmp_CS[1]
        sph_ctrs[i + (f_cont) * (n_sphs), 2] = after_rotation[2] + tmp_CS[2]
    if f_cont < f_num:
        for j in range((f_cont - 1) * (n_sphs), f_cont * (n_sphs)):
            for k in range(1, (f_cont - 1) * (n_sphs)):
                if (sph_ctrs[j, 0] - sph_ctrs[k, 0]) ** 2 + (sph_ctrs[j, 1] - sph_ctrs[k, 1]) ** 2 + (
                        sph_ctrs[j, 2] - sph_ctrs[k, 2]) ** 2 > (r_sphere + r_sphere) ** 2 + min_gap:
                    f_avlb = f_avlb * 1
                else:
                    f_avlb = f_avlb * 0
    else:
        break
    if f_avlb == 1:
        local_CSs.append(tmp_CS)
        f_cont += 1
    else:
        f_cont = f_cont

# Location for fibres
f1 = open('sph_ctrs.txt', 'w')
for i in range(len(sph_ctrs)):
    f1.write('%11.5f' % sph_ctrs[i][0] + '\t' + '%11.5f' % sph_ctrs[i][1] + '\t' + '%11.5f' % sph_ctrs[i][2] + '\n')
f1.close()
f2 = open('local_CSs.txt', 'w')
for i in range(len(local_CSs)):
    f2.write('%11.5f' % local_CSs[i][0] + '\t' + '%11.5f' % local_CSs[i][1] + '\t' + '%11.5f' % local_CSs[i][
        2] + '\t' + '%11.5f' % local_CSs[i][3] + '\t' + '%11.5f' % local_CSs[i][4] + '\t' + '%11.5f' % local_CSs[i][
                 5] + '\n')
f2.close()
# plot fibres in matplob
RD1 = np.linspace(0, 2 * np.pi, 8)
RD2 = np.linspace(0, np.pi, 8)
# prepare for matlibplot plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([-r_disk, r_disk])
ax.set_ylim([-r_disk, r_disk])
ax.set_zlim([-thickness / 2, thickness / 2])  # for plotting unstreched model
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
for i in range(0, len(sph_ctrs) // 2):
    ax.scatter(sph_ctrs[i][0], sph_ctrs[i][1], sph_ctrs[i][2], c='r')  # 绘制数据点
for i in range(len(sph_ctrs) // 2, 2 * len(sph_ctrs) // 2):
    ax.scatter(sph_ctrs[i][0], sph_ctrs[i][1], sph_ctrs[i][2], c='b')  # 绘制数据点
# for i in range(len(sph_ctrs)//2, 2*len(sph_ctrs)//2):
#    cx = sph_ctrs[i][0] + r_sphere * np.outer(np.cos(RD1), np.sin(RD2))
#    cy = sph_ctrs[i][1] + r_sphere * np.outer(np.sin(RD1), np.sin(RD2))
#    cz = sph_ctrs[i][2] + r_sphere * np.outer(np.ones(np.size(RD1)), np.cos(RD2))
#    ax.plot_surface(cx, cy, cz, color='b')

plt.show()
