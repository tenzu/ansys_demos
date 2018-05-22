*KEYWORD
*PARAMETER
R err     0.28      
R sigc    40.0      
R gb      30.0      
R smax    0.07      
R exp     1.6       
*CONTROL_BULK_VISCOSITY
$#      q1        q2      type     btype
  1.500000 6.0000E-2         1         0
*CONTROL_CONTACT
$#  slsfac    rwpnal    islchk    shlthk    penopt    thkchg     orien    enmass
  0.100000     0.000         2         0         1         1         1         0
$#  usrstr    usrfrc     nsbcs    interm     xpene     ssthk      ecdt   tiedprj
         0         0        10         0  4.000000         0         0         0
$#   sfric     dfric       edc       vfc        th     th_sf    pen_sf
     0.000     0.000     0.000     0.000     0.000     0.000     0.000
$#  ignore    frceng   skiprwg    outseg   spotstp   spotdel   spothin
         0         0         0         0         0         0     0.000
$#    isym    nserod    rwgaps    rwgdth     rwksf      icov    swradf    ithoff
         0         0         1     0.000  1.000000         0     0.000         0
$#  shledg    pstiff    ithcnt    tdcnof     ftall    unused    shltrw
         0         0         0         0         0               0.000
*CONTROL_COUPLING
$#  unleng    untime    unforc    timidl     flipx     flipy     flipz    subcyl
  1.000000  1.000000  1.000000     0.000         0         0         0         1
*CONTROL_CPU
$#  cputim
     0.000
*CONTROL_DYNAMIC_RELAXATION
$#  nrcyck     drtol    drfctr    drterm    tssfdr    irelal     edttl    idrflg
       250 1.0000E-3  0.9950001.0000E+30  0.900000         0 4.0000E-2         0
*CONTROL_ENERGY
$#    hgen      rwen    slnten     rylen
         1         2         1         1
$*CONTROL_HOURGLASS
$$#     ihq        qh
$         1  0.100000
*CONTROL_OUTPUT
$#   npopt    neecho    nrefup    iaccop     opifs    ipnint    ikedit    iflush
         0         3         0         0     0.000         0         0         0
$#   iprtf    ierode     tet10    msgmax    ipcurv      gmdt   ip1dblt      eocs
         0         0         2        50         0     0.000         0         0
*CONTROL_SHELL
$#  wrpang     esort     irnxx    istupd    theory       bwc     miter      proj
 20.000000         2        -1         0         2         2         1         0
$# rotascl    intgrd    lamsht    cstyp6    tshell
  1.000000         0         0         1         0         0         0         0
$# psstupd   sidt4tu     cntco    itsflg    irquad
         0         0         0         0         2
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas
  3.000000         0     0.000     0.000     0.000
*CONTROL_TIMESTEP
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
     0.000  0.900000         0     0.000     0.000         0         0         0
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl
     0.000         0         0                         0.000
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
         6         6         3         1         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         4         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0  1.000000         0         0         0STRESS    STRESS
$#    dtdt    resplt
         0         0
*CONTACT_1D
$#   nsids     nsidm       err      sigc        gb      smax       exp
         1         2      &err     &sigc       &gb     &smax      &exp 
*MAT_CONCRETE_DAMAGE_REL3
$#     mid        ro        pr
         1 2.7500E-3  0.200000
$#      ft        a0        a1        a2        b1     omega       a1f
  5.000000-140.00000     0.000     0.000     0.000     0.000     0.000
$# slambda      nout     edrop     rsize       ucf    lcrate  locwidth      npts
     0.000     0.000     0.000  0.393700 145.00000      1111     0.000     0.000
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
     0.000     0.000     0.000  0.200000     0.000     0.000     0.000     0.000
$#    idam    dmgtyp     lcsdg     ecrit    dmgexp     dcrit    fadexp    lcregd
         0     0.000         0     0.000  1.000000     0.000  1.000000     0.000
$#   lcfld             epsthin    engcrt    radcrt
         0         0     0.000     0.000     0.000
*SECTION_BEAM
$#   secid    elform      shrf   qr/irid       cst     scoor       nsm
         2         1     0.000         0         1     0.000     0.000
$#     ts1       ts2       tt1       tt2     nsloc     ntloc
  0.560000  0.560000     0.000     0.000     0.000     0.000
*MAT_PLASTIC_KINEMATIC
$#     mid        ro         e        pr      sigy      etan      beta
         2 7.8000E-3 2.0000E+5  0.280000 1280.0000 81.000000  1.000000
$#     src       srp        fs        vp
     0.000     0.000     0.000     0.000
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
