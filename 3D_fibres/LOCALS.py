import math, random
import numpy as np
import matplotlib.pyplot as plt

r_disk = 37.5  # disk radius (mm)
td_ratio = 0.4
thickness = 2 * r_disk * td_ratio  # disk thickness (mm)

# Original segment center for cylindars (fibres)
r = 2.5  # spiral fibre rotation radius (mm)
s = 10  # spiral fibre pitch (mm)
n = 8  # division in single pitch (n equals to 2**?)
p = 1.125  # spiral fibre total pitch number (should equal to an integer plus integral multiple of 1/n)
theta = p * 360  # spiral fibre total rotation angle (DEG)
bargin = r * 0.1  # bargin from fibre to outer cylindar (mm)
r_sphere = math.pi * 2 * (r + r + bargin) / 2 / n / 2  # radius of sphere in cylindar (mm)
L_sphere = math.ceil((s * p + 2 * r_sphere) / r_sphere)  # number of spheres in cylindar longitudinal direction


# generate an odd number (to make sure L_sphere is an odd number)
def odd_int(x):
    if x % 2 == 0:
        x += 1
    return (x)


n_sphs = n * odd_int(L_sphere)  # number of spheres in one cylindar
org_ctrs = np.zeros((n_sphs, 3))  # sphere centers in cylindar (original zeros)
min_gap = 1  # minimum gap between blocks


# generate local cordinate systems
def local_CS():
    tmp_CS = [0 for i in range(0, 6, 1)]  # local cordinate systems (x,y,z,yz,zx,xy)
    for i in range(0, 6, 1):
        if i == 0:
            # r~=r_disk-s*p-min_gap, x between (-r~,r~)
            r_ = r_disk - s * p - min_gap
            tmp_CS[i] = random.uniform(-r_, r_)
        elif i == 1:
            tmp_CS[i] = random.uniform(-math.sqrt(r_ ** 2 - tmp_CS[0] ** 2),
                                       math.sqrt(r_ ** 2 - tmp_CS[0] ** 2))
        elif i == 2:
            # t~=2*r_disk*td_ratio, z between (-(t~/2-s*p/2),t~/2-s*p/2)
            t_ = 2 * r_disk * td_ratio - 2 * min_gap
            tmp_CS[i] = random.uniform(-(t_ / 2 - s * p / 2), t_ / 2 - s * p / 2)
        else:
            tmp_CS[i] = random.uniform(0, 360)
    return (tmp_CS)


# define original centers of spheres
sph_ctrs = np.zeros((f_num * n_sphs, 3))  # centers of spheres
rd = np.linspace(0, 2 * math.pi, n + 1)  # angles of sphere centers
# generate list of [0,1,-1,2,-2,3,-3,...]
tmp_list = np.linspace(0, n_sphs / n, n_sphs / n + 1)
tmp_list = [((i + 1) // 2) * (-1) ** (i % 2) for i in tmp_list]
for i in range(0, n_sphs, 1):
    org_ctrs[i][0] = (r + r_sphere) * math.cos(rd[i % n])
    org_ctrs[i][1] = (r + r_sphere) * math.sin(rd[i % n])
    org_ctrs[i][2] = r_sphere * tmp_list[i // n]

# Locations for fibres
f_num = 4  # number of fibres
cld_cont = 1  # cylindar counter
while cld_cont <= f_num:
    cld_avlb = 1
    for i in range(0, n_sphs + 1, 1):
        sph_ctrs[i + (cld_cont - 1) * (n_sphs + 1), 0] = org_ctrs[i, 0] * math.cos(
            local_CS()[4] * math.pi / 180) * math.cos(local_CS()[5] * math.pi / 180) + local_CS()[0]
        sph_ctrs[i + (cld_cont - 1) * (n_sphs + 1), 1] = org_ctrs[i, 0] * math.cos(
            local_CS()[3] * math.pi / 180) * math.cos(local_CS()[5] * math.pi / 180) + local_CS()[1]
        sph_ctrs[i + (cld_cont - 1) * (n_sphs + 1), 2] = org_ctrs[i, 0] * math.cos(
            local_CS()[3] * math.pi / 180) * math.cos(local_CS()[4] * math.pi / 180) + local_CS()[2]
        if cld_cont <= f_num:
            for j in range():

        ##Location for fibres
        # NUM_F = 12   #fibre number
        # FC_CTR = numpy.zeros((NUM_F*(SPL_SEG+1),2)) #circle center on translated fibre (x, y)
        # F_LOCATIONS = numpy.zeros((NUM_F,3))    #fibre translation and rotation (x, y, angle)
        # F_CONT = 1  #fibre counter
        # while F_CONT <= NUM_F:
        #    F_AVLB = 1
        #    theta_F = random.randint(0.,360.)
        #    x_F = random.uniform(-(DISK_R-F_LENG/2.-0.00024133e3),(DISK_R-F_LENG/2.-0.00024133e3))
        #    y_F = random.uniform(-math.sqrt(0.02975867e3**2-x_F**2),math.sqrt(0.02975867e3**2-x_F**2))
        #    for i in range(0,SPL_SEG+1):
        #        FC_CTR[i+(F_CONT-1)*(SPL_SEG+1),0] = ORG_CTR[i,0]*math.cos(theta_F*pi/180.)-ORG_CTR[i,1]*math.sin(theta_F*pi/180.)+x_F
        #        FC_CTR[i+(F_CONT-1)*(SPL_SEG+1),1] = ORG_CTR[i,0]*math.sin(theta_F*pi/180.)+ORG_CTR[i,1]*math.cos(theta_F*pi/180.)+y_F
        #    if F_CONT <= NUM_F:
        #        for j in range((F_CONT-1)*(SPL_SEG+1),F_CONT*(SPL_SEG+1)):
        #            for k in range(1,(F_CONT-1)*(SPL_SEG+1)):
        #                if (FC_CTR[j,0]-FC_CTR[k,0])**2+(FC_CTR[j,1]-FC_CTR[k,1])**2 > (SPL_CR+SPL_CR)**2+MIN_GAP:
        #                    F_AVLB = F_AVLB*1
        #                else:
        #                    F_AVLB = F_AVLB*0
        #    else:
        #       break
        #    if F_AVLB == 1:
        #        F_LOCATIONS[F_CONT-1,0] = x_F
        #        F_LOCATIONS[F_CONT-1,1] = y_F
        #        F_LOCATIONS[F_CONT-1,2] = theta_F
        #        F_CONT += 1
        #    else:
        #        F_CONT = F_CONT
        #
        ##Location for fibres
        # f1 = open('ORG_CTR.txt','w')
        # for i in range(len(ORG_CTR)):
        #    f1.write('%11.5f' % ORG_CTR[i,0] +' '+'%11.5f' % ORG_CTR[i,1] +'\n')
        # f1.close()
        # f2 = open('F_LOCATIONS.txt','w')
        # for i in range(len(F_LOCATIONS)):
        #    f2.write('%11.5f' % F_LOCATIONS[i,0] +' '+'%11.5f' % F_LOCATIONS[i,1]+' '+'%11.5f' % F_LOCATIONS[i,2] +'\n')
        # f2.close()
        # f3 = open('FC_CTR.txt','w')
        # for i in range(len(FC_CTR)):
        #    f3.write('%11.5f' % FC_CTR[i,0] +' '+'%11.5f' % FC_CTR[i,1]+'\n')
        # f3.close()
        ##plot fibres in matplob
        ##for i in range(0,len(ORG_CTR)):
        ##    cx = ORG_CTR[i,0]+SPL_CR*numpy.cos(RD)
        ##    cy = ORG_CTR[i,1]+SPL_CR*numpy.sin(RD)
        ##    plt.plot(cx,cy,'-b')
        # for i in range(0,len(FC_CTR)):
        #    cx = FC_CTR[i,0]+SPL_CR*numpy.cos(RD)
        #    cy = FC_CTR[i,1]+SPL_CR*numpy.sin(RD)
        #    plt.plot(cx,cy,'-r')