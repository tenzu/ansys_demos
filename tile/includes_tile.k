$# LS-DYNA Keyword file created by LS-PrePost(R) V4.3 (Beta) - 02Jan2016(09:00)
$# Created on Sep-05-2017 (09:31:26)
*KEYWORD
*PARAMETER
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R err     0.28                                                                  
R sigc    40.0                                                                  
R gb      30.0                                                                  
R smax    0.07                                                                  
R exp     1.6                                                                   
*TITLE
$#                                                                         title
LS-DYNA keyword deck by LS-PrePost
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
         6         6         3         1         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         4         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0       1.0         0         0         0STRESS    STRESS
$#    dtdt    resplt     neipb     
         0         0         0
*SECTION_BEAM
$#   secid    elform      shrf   qr/irid       cst     scoor       nsm   
         2         1       0.0         1         1       0.0       0.0
$#     ts1       ts2       tt1       tt2     nsloc     ntloc     
      0.56      0.56       0.0       0.0       0.0       0.0
*MAT_ADD_EROSION
$#     mid      excl    mxpres     mneps    effeps    voleps    numfip       ncs
         1       0.0       0.0       0.0       0.0       0.0       1.0       1.0
$#  mnpres     sigp1     sigvm     mxeps     epssh     sigth   impulse    failtm
       0.0       0.0       0.0       0.2       0.0       0.0       0.0       0.0
$#    idam    dmgtyp     lcsdg     ecrit    dmgexp     dcrit    fadexp    lcregd
         0       0.0         0       0.0       1.0       0.0       1.0         0
$#   lcfld             epsthin    engcrt    radcrt      
         0         0       0.0       0.0       0.0
*MAT_CONCRETE_DAMAGE_REL3
$#     mid        ro        pr  
         1   0.00275       0.2
$#      ft        a0        a1        a2        b1     omega       a1f   
       5.0    -140.0       0.0       0.0       0.0       0.0       0.0
$# slambda      nout     edrop     rsize       ucf    lcrate  locwidth      npts
       0.0       0.0       0.0    0.3937     145.0      1111       0.0       0.0
$# lambda1   lambda2   lambda3   lambda4   lambda5   lambda6   lambda7   lambda8
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#lambda09  lambda10  lambda11  lambda12  lambda13        b3       a0y       a1y
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#    eta1      eta2      eta3      eta4      eta5      eta6      eta7      eta8
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#   eta09     eta10     eta11     eta12     eta13        b2       a2f       a2y
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
*MAT_PLASTIC_KINEMATIC
$#     mid        ro         e        pr      sigy      etan      beta    
         2    0.0078  200000.0      0.28    1280.0      81.0       1.0
$#     src       srp        fs        vp  
       0.0       0.0       0.0       0.0
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
      1111         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
               -30.0           6.0999999
         -0.30000001           6.0999999
                -0.1           6.0999999
               -0.05           6.0999999
              -0.035           5.5999999
               -0.02          4.88000011
               -0.01          4.07000017
      -1.000000e-004          1.29999995
      -1.000000e-006          1.17999995
      -1.000000e-009                 1.0
                 0.0                 1.0
       3.000000e-008                 1.0
       1.000000e-007          1.22000003
               0.001          1.28999996
                0.22          1.33599997
          0.30000001          1.47500002
                 1.0          2.65549994
                30.0          2.65549994
$*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE_ID
$$#     cid                                                                 title
$         1                                                                      
$$#    ssid      msid     sstyp     mstyp    sboxid    mboxid       spr       mpr
$         2         1         3         3         0         0         0         0
$$#      fs        fd        dc        vc       vdc    penchk        bt        dt
$       0.0       0.0       0.0       0.0       0.0         0       0.01.00000E20
$$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf
$       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0
*CONTROL_BULK_VISCOSITY
$#      q1        q2      type     btype     
       1.5      0.06         1         0
*CONTROL_COUPLING
$#  unleng    untime    unforc    timidl     flipx     flipy     flipz    subcyl
       1.0       1.0       1.0       0.0         0         0         0         1
*CONTROL_CPU
$#  cputim      
       0.0
*CONTROL_DYNAMIC_RELAXATION
$#  nrcyck     drtol    drfctr    drterm    tssfdr    irelal     edttl    idrflg
       250     0.001     0.9951.00000E300.89999998         0      0.04         0
*CONTROL_ENERGY
$#    hgen      rwen    slnten     rylen     
         1         2         1         1
$$#     ihq        qh
$         1  0.100000
*CONTROL_SHELL
$#  wrpang     esort     irnxx    istupd    theory       bwc     miter      proj
      20.0         2        -1         0         2         2         1         0
$# rotascl    intgrd    lamsht    cstyp6    tshell      
       1.0         0         0         1         0
$# psstupd   sidt4tu     cntco    itsflg    irquad      
         0         0         0         0         2
$#  nfail1    nfail4   psnfail    keepcs     delfr   drcpsid    drcprm      
         0         0         0         0         0         0       1.0
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas      
       4.0         0       0.0       0.0       0.0
*CONTROL_TIMESTEP
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
       0.00.89999998         0       0.0       0.0         0         0         0
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl     
       0.0         0         0                           0.0
*BOUNDARY_SPC_SET
$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
         1         0         1         1         1         0         0         0
*END
