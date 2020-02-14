$# LS-DYNA Keyword file created by LS-PrePost(R) V4.6.7 - 09May2019
$# Created on Feb-03-2020 (15:58:59)
*KEYWORD
*TITLE
$#                                                                         title
LS-DYNA keyword deck by LS-PrePost
*CONTROL_SHELL
$#  wrpang     esort     irnxx    istupd    theory       bwc     miter      proj
      20.0         1        -1         0         2         2         1         0
$# rotascl    intgrd    lamsht    cstyp6    thshel      
       1.0         0         0         1         0
$# psstupd   sidt4tu     cntco    itsflg    irquad    w-mode   stretch      icrq
         0         0         0         0         2       0.0       0.0         0
$#  nfail1    nfail4   psnfail    keepcs     delfr   drcpsid    drcprm   intperr
         0         0         0         0         0         0       1.0         0
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol     
    1000.0         0       0.0       0.0       0.0         0
*CONTROL_TIMESTEP
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
       0.0       0.9         0       0.0       0.0         0         0         0
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl     
       0.0         0         0                           0.0
*DATABASE_BINARY_D3DUMP
$$$5000000
$#    cycl      lcdt      beam     npltc    psetid      
   1000000         0         0         0         0
$$$*DATABASE_ELOUT
$$$$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
$$$ 1.0000E-3         1         0         1         0         0         0         0
$$$*DATABASE_NODOUT
$$$$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
$$$ 1.0000E-3         1         0         1         0         0         0         0
*DATABASE_BINARY_D3PLOT
$$$100.00
$#      dt      lcdt      beam     npltc    psetid      
      10.0         0         0         0         0
*DATABASE_BINARY_D3THDT
$$$25.000
$#      dt      lcdt      beam     npltc    psetid      
       5.0         0         0         0         0
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
         6         6         3         1         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         4         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0       1.0         0         0         0STRESS    STRESS
$#    dtdt    resplt     neipb     quadr     cubic     
         0         0         0         0         0
$$$*BOUNDARY_PRESCRIBED_MOTION_SET
$$$$#    nsid       dof       vad      lcid        sf       vid     death     birth
$$$         2         1         2         1       1.0         0       0.0       0.0
$$$*BOUNDARY_SPC_SET
$$$$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
$$$         1         0         1         1         1         1         1         1
*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE_ID
$#     cid                                                                 title
         1                                                                     
$#    ssid      msid     sstyp     mstyp    sboxid    mboxid       spr       mpr
         1         2         3         3         0         0         0         0
$#      fs        fd        dc        vc       vdc    penchk        bt        dt
      0.15       0.1       0.0       0.0       0.2         0       0.09.000000E7
$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf
       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0
*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE_ID
$#     cid                                                                 title
         2                                                                     
$#    ssid      msid     sstyp     mstyp    sboxid    mboxid       spr       mpr
         2         3         3         3         0         0         0         0
$#      fs        fd        dc        vc       vdc    penchk        bt        dt
      0.15       0.1       0.0       0.0       0.2         0       0.09.000000E7
$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf
       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0
*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE_ID
$#     cid                                                                 title
         3                                                                     
$#    ssid      msid     sstyp     mstyp    sboxid    mboxid       spr       mpr
         1         3         3         3         0         0         0         0
$#      fs        fd        dc        vc       vdc    penchk        bt        dt
      0.15       0.1       0.0       0.0       0.2         0       0.09.000000E7
$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf
       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0
$$$*CONTACT_TIEBREAK_SURFACE_TO_SURFACE_ID
$$$$#     cid                                                                 title
$$$         1
$$$$#    ssid      msid     sstyp     mstyp    sboxid    mboxid       spr       mpr
$$$         1         2         3         3         0         0         0         0
$$$$#      fs        fd        dc        vc       vdc    penchk        bt        dt
$$$      0.15       0.1       0.0       0.0       0.2         0       0.09.000000E7
$$$$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf
$$$       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0
$$$$#    nfls      sfls    tblcid    thkoff
$$$     100.0       3.5         0         0
$$$*CONTACT_TIEBREAK_SURFACE_TO_SURFACE_ID
$$$$#     cid                                                                 title
$$$         2
$$$$#    ssid      msid     sstyp     mstyp    sboxid    mboxid       spr       mpr
$$$         2         3         3         3         0         0         0         0
$$$$#      fs        fd        dc        vc       vdc    penchk        bt        dt
$$$      0.15       0.1       0.0       0.0       0.2         0       0.09.000000E7
$$$$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf
$$$       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0
$$$$#    nfls      sfls    tblcid    thkoff
$$$     100.0       3.5         0         0
$$$*CONTACT_TIEBREAK_SURFACE_TO_SURFACE_ID
$$$$#     cid                                                                 title
$$$         3
$$$$#    ssid      msid     sstyp     mstyp    sboxid    mboxid       spr       mpr
$$$         1         3         3         3         0         0         0         0
$$$$#      fs        fd        dc        vc       vdc    penchk        bt        dt
$$$      0.15       0.1       0.0       0.0       0.2         0       0.09.000000E7
$$$$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf
$$$       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0
$$$$#    nfls      sfls    tblcid    thkoff
$$$     100.0       3.5         0         0
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$                              PARTS DEFINITIONS                               $
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*SECTION_SHELL
$#   secid    elform      shrf       nip     propt   qr/irid     icomp     setyp
         1        13       1.0         1       1.0         0         0         1
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
       1.0       1.0       1.0       1.0       0.0       0.0       0.0         0
*MAT_ADD_EROSION
$#     mid      excl    mxpres     mneps    effeps    voleps    numfip       ncs
         1       0.0       0.0       0.0       0.0       0.0       1.0       1.0
$#  mnpres     sigp1     sigvm     mxeps     epssh     sigth   impulse    failtm
       0.0       0.0     210.0      0.25       0.0       0.0       0.0       0.0
$#    idam         -         -         -         -         -         -    lcregd
         0         0         0         0         1         0         1         0
$#   lcfld      nsff   epsthin    engcrt    radcrt   lceps12   lceps13   lcepsmx
         0         0       0.0       0.0       0.0         0         0         0
$#  dteflt    unused     mxtmp     dtmin     
       0.0                 0.0       0.0
$$$*MAT_ADD_EROSION
$$$$#     mid      excl    mxpres     mneps    effeps    voleps    numfip       ncs
$$$         2       0.0       0.0       0.05.00000E-4    -0.001       1.0       1.0
$$$$#  mnpres     sigp1     sigvm     mxeps     epssh     sigth   impulse    failtm
$$$       0.0       0.0      40.0       0.2       0.0       0.0       0.0       0.0
$$$$#    idam         -         -         -         -         -         -    lcregd
$$$         0         0         0         0         1         0         1         0
$$$$#   lcfld      nsff   epsthin    engcrt    radcrt   lceps12   lceps13   lcepsmx
$$$         0         0       0.0       0.0       0.0         0         0         0
$$$$#  dteflt    unused     mxtmp     dtmin     
$$$       0.0                 0.0       0.0
*MAT_ADD_EROSION
$#     mid      excl    mxpres     mneps    effeps    voleps    numfip       ncs
         3       0.0       0.0       0.05.00000E-4    -0.001       1.0       1.0
$#  mnpres     sigp1     sigvm     mxeps     epssh     sigth   impulse    failtm
       0.0       0.0      1.71     0.001       0.0       0.0       0.0       0.0
$#    idam         -         -         -         -         -         -    lcregd
         0         0         0         0         1         0         1         0
$#   lcfld      nsff   epsthin    engcrt    radcrt   lceps12   lceps13   lcepsmx
         0         0       0.0       0.0       0.0         0         0         0
$#  dteflt    unused     mxtmp     dtmin     
       0.0                 0.0       0.0
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$                             MATERIAL DEFINITIONS                             $
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*MAT_PIECEWISE_LINEAR_PLASTICITY_TITLE
fibre
$#     mid        ro         e        pr      sigy      etan      fail      tdel
         1    0.0078  200000.0      0.28    1300.0       0.0       0.2       0.0
$#       c         p      lcss      lcsr        vp  
       0.0       0.0         0         0       0.0
$#    eps1      eps2      eps3      eps4      eps5      eps6      eps7      eps8
       0.0     0.014     0.017     0.027     0.043     0.061     0.085     0.105
$#     es1       es2       es3       es4       es5       es6       es7       es8
     535.0     540.0     560.0     595.0     625.0     640.0     655.0     660.0
*MAT_CONCRETE_DAMAGE_REL3
$#     mid        ro        pr  
         2    0.0022      0.18
$#      ft        a0        a1        a2        b1     omega       a1f   
      35.0     -35.0       0.0       0.0       0.0       0.0       0.0
$# slambda      nout     edrop     rsize       ucf    lcrate  locwidth      npts
       0.0       0.0       0.0    0.3937     145.0      2222       0.0       0.0
$# lambda1   lambda2   lambda3   lambda4   lambda5   lambda6   lambda7   lambda8
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#lambda09  lambda10  lambda11  lambda12  lambda13        b3       a0y       a1y
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#    eta1      eta2      eta3      eta4      eta5      eta6      eta7      eta8
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#   eta09     eta10     eta11     eta12     eta13        b2       a2f       a2y
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
*MAT_CONCRETE_DAMAGE_REL3
$#     mid        ro        pr  
         3    0.0022      0.18
$#      ft        a0        a1        a2        b1     omega       a1f   
       3.0     -35.0       0.0       0.0       0.0       0.0       0.0
$# slambda      nout     edrop     rsize       ucf    lcrate  locwidth      npts
       0.0       0.0       0.0    0.3937     145.0      2222       0.0       0.0
$# lambda1   lambda2   lambda3   lambda4   lambda5   lambda6   lambda7   lambda8
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#lambda09  lambda10  lambda11  lambda12  lambda13        b3       a0y       a1y
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#    eta1      eta2      eta3      eta4      eta5      eta6      eta7      eta8
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#   eta09     eta10     eta11     eta12     eta13        b2       a2f       a2y
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$                               LOAD DEFINITIONS                               $
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$ Strain rate effect for aggregates (Hao and Hao 2013 in Rock Mechanics and Rock Engineering)
$ Strain rate effect for mortar matrix
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0                 0.0
             60000.0                10.0
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
      1111         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
               -30.0                 6.1
                -0.3                 6.1
                -0.1                 6.1
               -0.05                 6.1
              -0.035                 5.6
               -0.02                4.88
               -0.01                4.07
   -9.9999997474e-05                 1.3
   -9.9999999748e-07                1.18
   -9.9999997172e-10                 1.0
                 0.0                 1.0
    2.9999998929e-08                 1.0
    1.0000000117e-07                1.22
               0.001                1.29
                0.22               1.336
                 0.3               1.475
                 1.0              2.6555
                30.0              2.6555
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
      2222         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
               -30.0               10.78
                -1.0               10.78
                -0.3               10.78
                -0.1                7.47
               -0.05                5.93
               -0.03                 5.0
               -0.01                3.47
              -0.001                 1.6
   -9.9999997474e-05                 1.5
   -9.9999997474e-06                1.38
   -9.9999999392e-09                1.08
                 0.0                 1.0
    2.9999998929e-08                 1.0
    1.0000000117e-07                1.05
               0.001                1.22
                0.05               1.285
                 0.1               1.435
                 0.2               1.748
                 0.5               2.412
                 0.7               2.727
                 1.0                 3.1
                 3.0                 3.1
                30.0                 3.1
*END
