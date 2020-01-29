*CONTROL_SHELL
$#  wrpang     esort     irnxx    istupd    theory       bwc     miter      proj
 20.000000         1        -1         0         2         2         1         0
$# rotascl    intgrd    lamsht    cstyp6    tshell
  1.000000         0         0         1         0
$# psstupd   sidt4tu     cntco    itsflg    irquad
         0         0         0         0         2
$#  nfail1    nfail4   psnfail    keepcs     delfr   drcpsid    drcprm
         0         0         0         0         0         0  1.000000
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas
  10000.00         0     0.000     0.000     0.000
*CONTROL_TIMESTEP
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
     0.000  0.900000         0     0.000     0.000         0         0         0
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl
     0.000         0         0                         0.000
$$$*DATABASE_ELOUT
$$$$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
$$$ 1.0000E-3         1         0         1         0         0         0         0
$$$*DATABASE_NODOUT
$$$$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
$$$ 1.0000E-3         1         0         1         0         0         0         0
*DATABASE_BINARY_D3PLOT
100.00    
*DATABASE_BINARY_D3THDT
25.000
*DATABASE_BINARY_D3DUMP
5000000
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
         6         6         3         1         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         4         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0  1.000000         0         0         0STRESS    STRESS
$#    dtdt    resplt
         0         0
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$                              PARTS DEFINITIONS                               $
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*SECTION_SHELL
$#   secid    elform      shrf       nip     propt   qr/irid     icomp     setyp
         1        13  1.000000         1         1         0         0         1
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
  1.000000  1.000000  1.000000  1.000000     0.000     0.000     0.000         
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$                             MATERIAL DEFINITIONS                             $
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*MAT_PIECEWISE_LINEAR_PLASTICITY_TITLE
fibre
$#     mid        ro         e        pr      sigy      etan      fail      tdel
         1 7.8000E-3 2.0000E+5  0.280000 1300.0000     0.000  0.200000     0.000
$#       c         p      lcss      lcsr        vp
     0.000     0.000         0         0     0.000
$#    eps1      eps2      eps3      eps4      eps5      eps6      eps7      eps8
     0.000 1.4000E-2 1.7000E-2 2.7000E-2 4.3000E-2 6.1000E-2 8.5000E-2  0.105000
$#     es1       es2       es3       es4       es5       es6       es7       es8
 535.00000 540.00000 560.00000 595.00000 625.00000 640.00000 655.00000 660.00000
*MAT_CONCRETE_DAMAGE_REL3
$#     mid        ro        pr
         2 2.2000E-3  0.180000
$#      ft        a0        a1        a2        b1     omega       a1f
  3.000000-35.000000     0.000     0.000     0.000     0.000     0.000
$# slambda      nout     edrop     rsize       ucf    lcrate  locwidth      npts
     0.000     0.000     0.000  0.393700 145.00000      2222     0.000     0.000
$# lambda1   lambda2   lambda3   lambda4   lambda5   lambda6   lambda7   lambda8
     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000
$#lambda09  lambda10  lambda11  lambda12  lambda13        b3       a0y       a1y
     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000
$#    eta1      eta2      eta3      eta4      eta5      eta6      eta7      eta8
     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000
$#   eta09     eta10     eta11     eta12     eta13        b2       a2f       a2y
     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000
*MAT_ADD_EROSION
$#     mid      excl    mxpres     mneps    effeps    voleps    numfip       ncs
         1     0.000     0.000     0.000     0.000     0.000  1.000000  1.000000
$#  mnpres     sigp1     sigvm     mxeps     epssh     sigth   impulse    failtm
     0.000     0.000     210.0  0.250000     0.000     0.000     0.000     0.000
$#    idam    dmgtyp     lcsdg     ecrit    dmgexp     dcrit    fadexp    lcregd
         0     0.000         0     0.000  1.000000     0.000  1.000000     0.000
$#   lcfld             epsthin    engcrt    radcrt
         0         0     0.000     0.000     0.000
*MAT_ADD_EROSION
$#     mid      excl    mxpres     mneps    effeps    voleps    numfip       ncs
         2     0.000     0.000     0.000   0.00050    -0.001  1.000000  1.000000
$#  mnpres     sigp1     sigvm     mxeps     epssh     sigth   impulse    failtm
     0.000     0.000     1.710    0.0010     0.000     0.000     0.000     0.000
$#    idam    dmgtyp     lcsdg     ecrit    dmgexp     dcrit    fadexp    lcregd
         0     0.000         0     0.000  1.000000     0.000  1.000000     0.000
$#   lcfld             epsthin    engcrt    radcrt
         0         0     0.000     0.000     0.000
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$                               LOAD DEFINITIONS                               $
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$ Strain rate effect for aggregates (Hao and Hao 2013 in Rock Mechanics and Rock Engineering)
$ Strain rate effect for mortar matrix
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp
      1111         0  1.000000  1.000000     0.000     0.000         0
$#                a1                  o1
          -30.000000            6.100000
           -0.300000            6.100000
           -0.100000            6.100000
      -5.000000e-002            6.100000
      -3.500000e-002            5.600000
      -2.000000e-002            4.880000
      -1.000000e-002            4.070000
      -1.000000e-004            1.300000
      -1.000000e-006            1.180000
      -1.000000e-009            1.000000
               0.000            1.000000
       3.000000e-008            1.000000
       1.000000e-007            1.220000
       1.000000e-003            1.290000
            0.220000            1.336000
            0.300000            1.475000
            1.000000            2.655500
           30.000000            2.655500
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp
      2222         0  1.000000  1.000000     0.000     0.000         0
$#                a1                  o1
          -30.000000           10.780000
           -1.000000           10.780000
           -0.300000           10.780000
           -0.100000            7.470000
      -5.000000e-002            5.930000
      -3.000000e-002            5.000000
      -1.000000e-002            3.470000
      -1.000000e-003            1.600000
      -1.000000e-004            1.500000
      -1.000000e-005            1.380000
      -1.000000e-008            1.080000
               0.000            1.000000
       3.000000e-008            1.000000
       1.000000e-007            1.050000
       1.000000e-003            1.220000
       5.000000e-002            1.285000
            0.100000            1.435000
            0.200000            1.748000
            0.500000            2.412000
            0.700000            2.727000
            1.000000            3.100000
            3.000000            3.100000
           30.000000            3.100000
