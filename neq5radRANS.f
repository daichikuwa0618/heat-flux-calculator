*************************************************************
c     *   Thermochemical nonequilibrium axisymmetroc or 2D NS solver    *
c     *                      for Earh entry on 2D structured mesh       *
c     *                                                                 *
c     *   Thermochemical model                                          *
c     *     - Air neutral 5-species (O,N,NO,O2,N2)                      *
c     *       subscript: 1:O, 2:N, 3:NO, 4:O2, 5:N2                     *
c     *     - Park's two-temparature model                              *
c     *     - Chemical reaction                                         *
c     *       Forward rate coefficients: Park, JTHT, Vol.15, No1, 2001  *
c     *       O2 + M = O + O + M                                        *
c     *       N2 + M = N + N + M                                        *
c     *       NO + M = N + N + M                                        *
c     *       N2 + O = NO + N  Bose and Candler                         *
c     *       NO + O = O2 + N  Bose and Candler                         *
c     *       (M: third body)                                           *
c     *       Coefficients for equibrium constant: Park's Book          *
c     *       "nonequilibrium hypersonic aerothermodynamics",           *
c     *        at number density of 1.0e17 cm-3, for reaction 1,3 and 5,*
c     *        a2 is increased by 13.8155 for mole/m3 unit              *
c     *                                                                 *
c     *     - Vibrational energy source term                            *
c     *       Translational-vibrational Energy relaxation:              *
c     *       Landau-Teller + Park's two modification                   *
c     *       (1) Park's average collision time is added to             *
c     *           Millikan-White time                                   *
c     *           The coefficients for relaxation time:                 *
c     *           Park, JTHT Vol.7, No.3, 1993                          *
c     *       (2) Bridging model is added to Landau-Teller equation     *
c     *           this model is more valid when Tshock > 10000K         *
c     *           ON(nrlx=1),OFF(nrlx=0)                                *
c     *           In one dimensional flow, the determination of the post*
c     *           shock temperatures is trivial, However, in the        *
c     *           multidimensional flow, a treatment of the model calls *
c     *           for tracing the stream line back to the point of      *
c     *           origin at the shock which becomes higher cost         *
c     *           So in multidimensional flow over complex geometry,    *
c     *           Bridgine model is set to be OFF(nflx=0) in NASA LAURA *
c     *                                                                 *
c     *     - Transport properties                                      *
c     *       1. species viscocities: Blottner's curv fitting           *
c     *          species conductivities of translational-rotational     *
c     *          temperature and vibrational- electronic temperature    *
c     *          : Eucken's relation                                    *
c     *          total viscosity and conductivity                       *
c     *          : Wilke's semi-empirical mixing rule                   *
c     *          diffusion coefficient: constant Schmidt number (nvis=0)*
c     *       2. total viscosity and conductivity:Collision Integrals   *
c     *          Gupta curve fitting for collision integrals            *
c     *          diffusion coefficients: binary diffusion               *
c     *          (nvis=1)                                               *
c     *                                                                 *
c     *   Numerical methods                                             *
c     *     - Mass conservation equations for O, N, NO are solved       *
c     *       Mass for O2, N2 is determind by mass of O,N,NO and inflow *
c     *       fraction                                                  *
c     *     - Higher order spatial accuracy is archieved by MUSCL for   *
c     *       primitive value (nacu=2)                                  *
c     *       artificial compression parameter == 1                     *
c     *     - Time integration                                          *
c     *       Lower-upper symmetric-Gauss-Seidel (LU-SGS)               *
c     *       The diagonal point implicit method is employed            *
c     *       for chemical source term                                  *
c     *     - local time step(ntim=1) or global time step(ntim=0)       *
c     *     - numerical flux                                            *
c     *       convective :                                              *
c     *       1. AUSM-DV (shock-fix)                                    *
c     *       2. SLAU                                                   *
c     *       3. SLAU2                                                  *
c     *       viscous    : 2nd order central difference                 *
c     *                                                                 *
c     *     - RANS model : Baldwin-Lomax (nprc=2)                       *
c     *       fully-turbulence or with transition model                 *
c     *                                                                 *
c     *     - Radiation loose-coupled calculation for 5 species         *
c     *       1-dimensional radiative transfer calculation              *
c     *       with tangent-slab approximation for 5 species             *
c     *       spectrally detailed absorption using multi-band model     *
c     *       (sub. rad1,loose,radheat)                                 *
c     *     - output file for linebyline  (./output/lbl_rt.dat)         *
c     *             (when irad=2)                                       *
c     *             radiaion for 5 species                              *
c     *             (sub. lbldata)                                      *
c     *                                                                 *
c     *   output                                                        *
c     *     - res.his    : history of residual for density and          *
c     *                    convective heat flux                         *
c     *     - acf.his    : history of aerocoefficients                  *
c     *     - field data     (field.dat)                                *
c     *     - tecplot data   (tec.dat)                                  *
c     *     - heat flux data (qwxy.dat)                                 *
c     *                                                                 *
c     *******************************************************************
c     *   Tricks for more robust                                        *
c     *     -sub. tsol                                                  *
c     *      Restrict to tv or tt which is determined by newton iter.   *
c     *     -sub. sorc                                                  *
c     *      Restrict to effective temperature, Ta, Te                  *
c     *      change the telim in input file                             *
c     *     -sub. intg                                                  *
c     *      Restrict to ev                                             *
c     *      qv(i,j,1)  =dmax1(evinf,qv(i,j,1)+dvq(i,j,1))              *
c     *                                                                 *
c     *******************************************************************
c     *     for application                                             *
c     *     -Geometry                                                   *
c     *      Modify boundary conditions                                 *
c     *      (subroutine metr, subroutine rbdy and subroutine sbdy)     *
c     *     -Default boundary conditions                                *
c     *      i=1:axis  , i=imax: supersonic outflow                     *
c     *      j=1:wall  , j=jmax: inflow                                 *
c     *     -for 2D calculation                                         *
c     *      1. comment out 'call saxi'                                 *
c     *      2. comment out ' qcr(:,:,1-4) in Sub. intg'                *
c     *                                                                 *
c     *                                   update 2014.10.14 T.Ishihara  *
c     *******************************************************************
c
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /aray3/ qcl(nx,ny,4),qcr(nx,ny,4),
     &               qsl(nx,ny,3),qsr(nx,ny,3),
     &               qvl(nx,ny,1),qvr(nx,ny,1)
      common /aray4/ dcq(nx,ny,4),dsq(nx,ny,3),dvq(0:nx,0:ny,1)
      common /aray5/ ish(0:nx,0:ny)
      common /aray7/ scs(nx,ny,3),scv(nx,ny,1)
      common /aray8/ ajs(nx,ny,3),ajv(nx,ny,1)
      common /aray9/ tvev(nx,0:ny)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /coor3/ xg(0:nx,0:ny),yg(0:nx,0:ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /dtst1/ dt(nx,ny),dt0,dti
      common /cntl1/ istart
      common /cntl2/ cflmx,grwmx
      common /cntl3/ tlmt
      common /cntl4/ lodt
      common /cntl5/ nprc,iprc
      common /cntl6/ nacu
      common /cntl8/ nchm
      common /cntl9/ qq,fmul3,fmul4
      common /cntl10/nvis,nrlx
      common /parm1/ nstp,nout
      common /parm3/ mdd
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem2/ tss(nx)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /chem6/ telim
      common /wbcd1/ twall,cefn,cefo
      common /fscd1/ fsv,alpha
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /fscd3/ vmueinf,vkapinf,vkavinf,
     &               ds1inf,ds2inf,ds3inf,ds4inf,ds5inf,
     &               evs1inf,evs2inf,evs3inf,evs4inf,evs5inf,
     &               hs1inf,hs2inf,hs3inf,hs4inf,hs5inf
      common /stat1/ loop,locl
      common /stat2/ time
      common /cntrd/ irad  !
      common /file1/ file_mesh,file_field
      character*50 file_mesh,file_field
      character field*50

c read run_time parameters
      call parm
c read initial mesh data / metrices
      open(12,file=file_mesh,form='formatted')
      call grid
      close(12)
      call metr
c set distance from wall for Baldwin-Lomax
c this sub. is used in aero_coeff
      call wdst
c set initial condition
      call init
c read disk data (if requested)
      if(istart.ne.0) then
         open(15,file=file_field,form='formatted')
         call dred
         close(15)
      endif
c* history of density residual
      open(333,file='./res.his',form='formatted',status='unknown')
      rewind 333
      if (istart.ne.0) then
 800     read(333,'(i10,1p4e15.7,i4)',end=801) idummy
         go to 800
 801     backspace 333
      endif
c* history of aero_coeff
c if you want to calculate aerodynamic coefficient,
c please commen in following 5 lines
      open(77,file='acf.his',form='formatted')
      write(77,767)
      write(77,777) '#','loop','Cd','Pressure drag','Friction drag'
 767  format('# aero coefficient')
 777  format(a,a9,1x,3(a14))
c
c* initial setup sequence
c* entry/main sequence
      locl   =0
 1000 locl =locl+1
      loop =loop+1
c translational-rotational and vibrational-electronic temperatures
      call tsol
c highest temperature behind shock wave
      call htbs
c transport coefficients
      if(nprc.ne.0) call trpc
c RANS (Baldwin-Lomax)
      if(nprc.ge.2) then
         call rbdy(1)
         call rbdy(2)
         call baldwin_lomax
      end if
c calc. aerodynamic coefficient
c if you want to calculate aerodynamic coefficient,
c please comment in Sub. aero_coeff
c you can use this subroutine, only when axisymmetric calc.
      mw =mod(loop,nout)
      if(mw.eq.0) then
         call aero_coeff
         call cpcf_coeff
         call knudsen_gll
      end if
c shock detection using ausm-dv,rad1
      call shock
c i-integration
      call rbdy(1)
      call musl(1)
c select flux scheme
      call flux(1)
      call vflx(1)
c j-integraton
      call rbdy(2)
      call musl(2)
c select flux scheme
      call flux(2)
      call vflx(2)
c axisymmetric term
      call saxi
c chemical source term
      call sorc
c set time step
      call dtsz
c for rmshflx
      call heatrms
c integration
      call intg
c
c* post processing sequence
c save data
      m =mod(locl,nout)
      if(m.eq.0) then
         open(15,file=file_field,form='formatted')
         call wrtd
         close(15)
         call heat
         call tecp
      endif
c
      nouten = nout*10
      kl =mod(locl,nout)
      if(kl.eq.0) then
      k = loop/nout
      write (field,'(a,i6.6,a)') "./output/field",k,".dat"
      open(16,file=field,form='formatted')
         call wrtd2
      close(16)
      endif
c judge (continue / stop)
      if(locl.lt.nstp) then
         if(time.lt.tlmt) then
            goto 1000
         else
            if(m.ne.0) then
               open(15,file=file_field,form='formatted')
               call wrtd
               close(15)
               call heat
               call tecp
            endif
         endif
      endif
c
      WRITE(6,*) ' '
      WRITE(6,*) '<< END OF COMPUTATION >>'
c
c* termination
      stop
      end
c
      subroutine grid
c     **************************************************************
c     *     read grid system                                       *
c     **************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
c* read grid data
      write(6,*) ''
      write(6,*) '*** grid ***'
      read(12,*) imax,jmax
      write(6,*) 'imax=',imax,' jmax=',jmax
      write(6,*) ''
!      read(12,*) ((x(i,j),i=1,imax),j=1,jmax)
!      read(12,*) ((y(i,j),i=1,imax),j=1,jmax)
      read(12,*) ((x(i,j),y(i,j),i=1,imax),j=1,jmax)
c
c* termination
      return
      end
c
      subroutine dred
c     ************************************************************
c     *     read disk data  (istart.ne.0)                        *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /coor1/ imax,jmax
      common /stat1/ loop,locl
      common /stat2/ time
      common /hfrad/ qsh(nx,ny),qwa(nx,ny) !
      dimension icom(20),pcom(20)
c* main
      read(15,*) (icom(i),i=1,20)
      read(15,*) (pcom(i),i=1,20)
      read(15,*) ((qs(i,j,1),qs(i,j,2),qs(i,j,3),qv(i,j,1),
     &             qc(i,j,1),qc(i,j,2),qc(i,j,3),qc(i,j,4),
     &             i=1,imax),j=1,jmax)
      read(15,*) ((tt(i,j),tv(i,j),i=1,imax),j=1,jmax)
      read(15,*) ((vmue(i,j),vkap(i,j),vkav(i,j),i=1,imax),j=1,jmax)
      read(15,*) ((qsh(i,j),qwa(i,j),i=1,imax-1),j=1,jmax)  !
c* disp
      write(6,600) (icom(i),i=1,20)
  600 format(/1x,'dred:disp) icom(1 - 20) ...'/
     &       (1x,4i13))
      write(6,610) (pcom(i),i=1,20)
  610 format(/1x,'dred:disp) pcom(1 - 20) ...'/
     &       (1x,4e13.6))
c* set parameters
      loop =icom(5)
      time =pcom(1)
      write(6,620) loop,time
  620 format(/1x,'dred:disp) loop/time are reset to followings.'/
     &        1x,'                     loop = ',i13/
     &        1x,'                     time = ',e13.6)
c* termination
      return
      end
c
      subroutine wrtd
c     ************************************************************
c     *     write data                                           *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /dtst1/ dt(nx,ny),dt0,dti
      common /cntl1/ istart
      common /cntl2/ cflmx,grwmx
      common /cntl3/ tlmt
      common /cntl4/ lodt
      common /cntl5/ nprc,iprc
      common /cntl6/ nacu
      common /parm1/ nstp,nout
      common /parm3/ mdd
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem2/ tss(nx)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /wbcd1/ twall,cefn,cefo
      common /fscd1/ fsv,alpha
      common /stat1/ loop,locl
      common /stat2/ time
      common /hfrad/ qsh(nx,ny),qwa(nx,ny) !
      dimension icom(20),pcom(20)
      data icom/20*0/pcom/20*0.0/
c* set headder
      icom(1)  =imax
      icom(2)  =jmax
      icom(5)  =loop
      icom(6)  =nstp
      icom(7)  =nout
      icom(8)  =lodt
      icom(9)  =nprc
      icom(10) =iprc
      icom(11) =nacu
      icom(12) =mdd
      icom(13) =istart
      pcom(1)  =time
      pcom(2)  =cflmx
      pcom(3)  =grwmx
      pcom(4)  =dt0
      pcom(5)  =tlmt
      pcom(6)  =fsv
      pcom(7)  =alpha
      pcom(8)  =twall
      pcom(9)  =cefn
      pcom(10) =cefo
      pcom(11) =twall
      pcom(12) =bet1
      pcom(13) =bet2
      pcom(14) =bet3
      pcom(15) =bet4
      pcom(16) =bet5
      pcom(17) =bet6
c* write data
      write(15,*) (icom(i),i=1,20)
      write(15,*) (pcom(i),i=1,20)
      write(15,*) ((qs(i,j,1),qs(i,j,2),qs(i,j,3),qv(i,j,1),
     &              qc(i,j,1),qc(i,j,2),qc(i,j,3),qc(i,j,4),
     &              i=1,imax),j=1,jmax)
      write(15,*) ((tt(i,j),tv(i,j),i=1,imax),j=1,jmax)
      write(15,*) ((vmue(i,j),vkap(i,j),vkav(i,j),i=1,imax),j=1,jmax)
      write(15,*) ((qsh(i,j),qwa(i,j),i=1,imax-1),j=1,jmax)   !
c* termination
      return
      end
c
      subroutine wrtd2
c     ************************************************************
c     *     write data                                           *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /dtst1/ dt(nx,ny),dt0,dti
      common /cntl1/ istart
      common /cntl2/ cflmx,grwmx
      common /cntl3/ tlmt
      common /cntl4/ lodt
      common /cntl5/ nprc,iprc
      common /cntl6/ nacu
      common /parm1/ nstp,nout
      common /parm3/ mdd
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem2/ tss(nx)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /wbcd1/ twall,cefn,cefo
      common /fscd1/ fsv,alpha
      common /stat1/ loop,locl
      common /stat2/ time
      common /hfrad/ qsh(nx,ny),qwa(nx,ny) !
      dimension icom(20),pcom(20)
      data icom/20*0/pcom/20*0.0/
c* set headder
      icom(1)  =imax
      icom(2)  =jmax
      icom(5)  =loop
      icom(6)  =nstp
      icom(7)  =nout
      icom(8)  =lodt
      icom(9)  =nprc
      icom(10) =iprc
      icom(11) =nacu
      icom(12) =mdd
      icom(13) =istart
      pcom(1)  =time
      pcom(2)  =cflmx
      pcom(3)  =grwmx
      pcom(4)  =dt0
      pcom(5)  =tlmt
      pcom(6)  =fsv
      pcom(7)  =alpha
      pcom(8)  =twall
      pcom(9)  =cefn
      pcom(10) =cefo
      pcom(11) =twall
      pcom(12) =bet1
      pcom(13) =bet2
      pcom(14) =bet3
      pcom(15) =bet4
      pcom(16) =bet5
      pcom(17) =bet6
c* write data
      write(16,*) (icom(i),i=1,20)
      write(16,*) (pcom(i),i=1,20)
      write(16,*) ((qs(i,j,1),qs(i,j,2),qs(i,j,3),qv(i,j,1),
     &              qc(i,j,1),qc(i,j,2),qc(i,j,3),qc(i,j,4),
     &              i=1,imax),j=1,jmax)
      write(16,*) ((tt(i,j),tv(i,j),i=1,imax),j=1,jmax)
      write(16,*) ((vmue(i,j),vkap(i,j),vkav(i,j),i=1,imax),j=1,jmax)
      write(16,*) ((qsh(i,j),qwa(i,j),i=1,imax-1),j=1,jmax)   !
c* termination
      return
      end
c
      subroutine parm
c     **************************************************************
c     *     run_time parameters                                    *
c     **************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      parameter(narray=10,nsp=5) ! for radiation
      common /cntl1/ istart
      common /cntl2/ cflmx,grwmx
      common /cntl3/ tlmt
      common /cntl4/ lodt
      common /cntl5/ nprc,iprc
      common /cntl6/ nacu
      common /cntl8/ nchm
      common /cntl9/ qq,fmul3,fmul4
      common /cntl10/nvis,nrlx
      common /coor1/ imax,jmax                        !
      common /parm1/ nstp,nout
      common /parm3/ mdd
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /chem6/ telim
      common /wbcd1/ twall,cefn,cefo
      common /fscd1/ fsv,alpha
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /fscd3/ vmueinf,vkapinf,vkavinf,
     &               ds1inf,ds2inf,ds3inf,ds4inf,ds5inf,
     &               evs1inf,evs2inf,evs3inf,evs4inf,evs5inf,
     &               hs1inf,hs2inf,hs3inf,hs4inf,hs5inf
      common /stat1/ loop,locl
      common /stat2/ time
!!--- for radiation----------------------------------------------------
      common /cntrd/ irad
      common /absb1/ wavel(narray),dlam(0:narray),absbx(narray,5,nsp),
     &               absbk(nx,ny,narray)
      common /absb2/ nspec,nwavel
      common /efnct/ e1tbl(1501),e2tbl(1501),e3tbl(1501),atbl(1501),
     &               ccofe(1501,6)
      common /hfrad/ qsh(nx,ny),qwa(nx,ny)
!!---------------------------------------------------------------------
      common /file1/ file_mesh,file_field
      common /file2/ file_absb,file_eftb  ! for radiation
      character*50 file_mesh,file_field
      character*80 file_absb,file_eftb    ! for radiation
c Blottner curve fitting parameters
      data a1/0.0203144/,b1/ 0.4294404/,c1/-11.6031403/,
     &     a2/0.0115572/,b2/ 0.6031679/,c2/-12.4327495/,
     &     a3/0.0436378/,b3/-0.0335511/,c3/ -9.5767430/,
     &     a4/0.0449290/,b4/-0.0826158/,c4/ -9.2019475/,
     &     a5/0.0268142/,b5/ 0.3177838/,c5/-11.3155513/
c array for collisional integrals
      dimension xr(5),xcinf(5)
      dimension dlsr1(5,5),dlsr2(5,5),dlsrv(5,5)
      dimension alsr(5,5),ddsr(5,5)
      dimension xxdl1(5),xxdl2(5),xxdlv(5),axdl2(5)
      dimension dsb(5)
c coefficients for collisional integrals by Gupta curve fitting
      dimension aog(5,5,2),bog(5,5,2),cog(5,5,2),dog(5,5,2)
c Omega_sr(1,1) and Omega_sr(2,2)
      data aog/50*0.d0/
      data bog/-0.0034d0, 0.0048d0,-0.0179d0,-0.0226d0,-0.0139d0, ! b(1-5,1,1) ! for Omega_sr(1,1)
     &          0.0048d0,-0.0033d0,-0.0185d0,-0.0179d0,-0.0194d0, ! b(1-5,2,1)
     &         -0.0179d0,-0.0185d0,-0.0364d0,-0.0438d0,-0.0291d0, ! b(1-5,3,1)
     &         -0.0226d0,-0.0179d0,-0.0438d0,-0.0410d0,-0.0465d0, ! b(1-5,4,1)
     &         -0.0139d0,-0.0194d0,-0.0291d0,-0.0465d0,-0.0112d0, ! b(1-5,5,1)
     &         -0.0207d0, 0.0065d0,-0.0203d0,-0.0247d0,-0.0169d0, ! b(1-5,1,2) ! for Omega_sr(2,2)
     &          0.0065d0,-0.0118d0,-0.0196d0,-0.0203d0,-0.0190d0, ! b(1-5,2,2)
     &         -0.0203d0,-0.0196d0,-0.0453d0,-0.0522d0,-0.0385d0, ! b(1-5,3,2)
     &         -0.0247d0,-0.0203d0,-0.0522d0,-0.0485d0,-0.0558d0, ! b(1-5,4,2)
     &         -0.0169d0,-0.0190d0,-0.0385d0,-0.0558d0,-0.0203d0/ ! b(1-5,5,2)
      data cog/-0.0572d0,-0.4195d0, 0.0152d0, 0.1300d0,-0.0825d0, ! c(1-5,1,1) ! for Omega_sr(1,1)
     &         -0.4195d0,-0.0572d0, 0.0118d0, 0.0152d0, 0.0119d0, ! c(1-5,2,1)
     &          0.0152d0, 0.0118d0, 0.3825d0, 0.5352d0, 0.2324d0, ! c(1-5,3,1)
     &          0.1300d0, 0.0152d0, 0.5352d0, 0.4977d0, 0.5729d0, ! c(1-5,4,1)
     &         -0.0825d0, 0.0119d0, 0.2324d0, 0.5729d0,-0.1182d0, ! c(1-5,5,1)
     &          0.0780d0,-0.4467d0, 0.0730d0, 0.1783d0, 0.0143d0, ! c(1-5,1,2) ! for Omega_sr(2,2)
     &         -0.4467d0,-0.0960d0, 0.0478d0, 0.0730d0, 0.0239d0, ! c(1-5,2,2)
     &          0.0730d0, 0.0478d0, 0.5624d0, 0.7045d0, 0.4226d0, ! c(1-5,3,2)
     &          0.1783d0, 0.0730d0, 0.7045d0, 0.6475d0, 0.7590d0, ! c(1-5,4,2)
     &          0.0143d0, 0.0239d0, 0.4226d0, 0.7590d0, 0.0683d0/ ! c(1-5,5,2)
      data dog/ 4.9901d0, 5.7774d0, 3.9996d0, 3.3363d0, 4.5785d0, ! d(1-5,1,1) ! for Omega_sr(1,1)
     &          5.7774d0, 5.0452d0, 4.0590d0, 3.9996d0, 4.1055d0, ! d(1-5,2,1)
     &          3.9996d0, 4.0590d0, 2.4718d0, 1.7252d0, 3.2082d0, ! d(1-5,3,1)
     &          3.3363d0, 3.9996d0, 1.7252d0, 1.8302d0, 1.6185d0, ! d(1-5,4,1)
     &          4.5785d0, 4.1055d0, 3.2082d0, 1.6185d0, 4.8464d0, ! d(1-5,5,1)
     &          3.5658d0, 6.0426d0, 3.8818d0, 3.2517d0, 4.4195d0, ! d(1-5,1,2) ! for Omega_sr(2,2)
     &          6.0426d0, 4.3252d0, 4.0321d0, 3.8818d0, 4.1782d0, ! d(1-5,2,2)
     &          3.8818d0, 4.0321d0, 1.7669d0, 1.0738d0, 2.4507d0, ! d(1-5,3,2)
     &          3.2517d0, 3.8818d0, 1.0738d0, 1.2607d0, 0.8955d0, ! d(1-5,4,2)
     &          4.4195d0, 4.1782d0, 2.4507d0, 0.8955d0, 4.0900d0/ ! d(1-5,5,2)
c
c* i/o files
      write(6,*) '?input  : file name of mesh data           ...'
      read(5,500) file_mesh
  500 format(a50)
      write(6,*) '?input  : file name of field data          ...'
      read(5,500) file_field
c* entry control
      write(6,*) '?select : start from initial(0) / disk(1)  ...'
      read(5,*) istart
c* procedure control
      write(6,*) '?select : euler(0) NS_lam(1) NS_RANS_BL    ...'
      read(5,*) nprc
      nprc =max(nprc,0)
!      nprc =min(nprc,1)
      nprc =min(nprc,2)
      write(6,*) '?select : implicit(0)  explicit(1) scheme  ...'
      read(5,*) iprc
      write(6,*) '?select : frozen(0)  non-equilibrium(1)    ...'
      read(5,*) nchm
c* step control
      write(6,*) '?input  : number of steps (nstp)           ...'
      read(5,*) nstp
      write(6,*) '?input  : interval to store data (nout)    ...'
      read(5,*) nout
      write(6,*) '?input  : mdd (interval for residual est.) ...'
      read(5,*) mdd
      write(6,*) '?input  : maximum time to continue (tlmt)  ...'
      read(5,*) tlmt
c* step_size control
      write(6,*) '?select : global(0) or local(1) time step  ...'
      read(5,*) lodt
      lodt =max(lodt,0)
      lodt =min(lodt,1)
      write(6,*) '?input  : maximum allowable cfl limit      ...'
      read(5,*) cflmx
      write(6,*) '?input  : maximum growth rate              ...'
      read(5,*) grwmx
c* accuracy control
      write(6,*) '?input  : 1st order(1)  2nd order(2)       ...'
      read(5,*) nacu
c* run_time parameters
      write(6,*) '?input  : free stream density (kg/m^3)     ...'
      read(5,*) rinf
      write(6,*) '?input  : fraction of species mass in freestream'
      write(6,*) '     O, N, NO, O2, N2'
      read(5,*) x1inf,x2inf,x3inf,x4inf,x5inf
c
      x5inf=dmax1(0.d0,1.d0-(x1inf+x2inf+x3inf+x4inf))
c
      write(6,*) '?input  : free stream translational temperature(K)'
      read(5,*) ttinf
      write(6,*) '?input  : free stream vibrational temperature(K)'
      read(5,*) tvinf
      write(6,*) '?input  : free stream velocity (m/s)       ...'
      read(5,*) fsv
      write(6,*) '?input  : alpha (deg.)                     ...'
      read(5,*) alpha
      write(6,*) '?input  : wall temperature (0=adiabatic)   ...'
      read(5,*) twall
      if(nprc.eq.0) twall =-1.0
      write(6,*) '?input  : wall catalytic efficiency (N)    ...'
      read(5,*) cefn
      write(6,*) '?input  : wall catalytic efficiency (O)    ...'
      read(5,*) cefo
      write(6,*) '?input  : contribution of vibrational mode, q '
      write(6,*) '          ta = (T)**q * (Tv)**(1-q)'
      read(5,*) qq
      write(6,*) '?input  : limit temperature for Keq curve fitting '
      read(5,*) telim
      write(6,*) '?input  : acceleration parameter for reactions '
      read(5,*) fmul3
      write(6,*) '?input  : acceleration parameter for relaxations '
      read(5,*) fmul4
      write(6,*) '?input  : tranport calculation type'
      write(6,*) '        : Wilke,Blottner,Eucken =0,Gupta=1 '
      read(5,*) nvis
      write(6,*) '?input  : Park bridging model ON(1)/OFF(0) '
      read(5,*) nrlx
c radiation input
      write(6,*) '?input  : radiation coupling off(0) / on(1)...'
      write(6,*) '?input  : output data for linebyline(2)    ...'
      read(5,*) irad
c
c reset radiative heat flux
      do j=1,jmax
         do i=1,imax-1
            qsh(i,j)=0.0d0
            qwa(i,j)=0.0d0
         end do
      end do
c read data for radiation calculation
      if(irad.eq.1) then
         read(5,510) file_absb
         read(5,510) file_eftb
c absorption coefficients
         open(16,file=file_absb,form='formatted')
         read(16,300) nwavel
         read(16,310) (wavel(m),m=1,nwavel)
         read(16,510)
         read(16,320) nspec
         do j=1,nspec
            read(16,510)
            read(16,510)
            read(16,310) ((absbx(l,i,j),l=1,nwavel),i=1,5)
         end do
         close(16)
c if trapezoidal integration for radiative heat flux
         dlam(0) =0.0d0
         do i=1,nwavel-1
            dlam(i)=wavel(i+1)-wavel(i)
         end do
         dlam(nwavel) =0.0d0
c exponential function
         itble=1501
         open(17,file=file_eftb,form='formatted')
         read(17,400) (atbl(i),e1tbl(i),e2tbl(i),e3tbl(i),i=1,itble)
         close(17)
         do i=1,itble-1
            ccofe(i,1)=e1tbl(i)-(e1tbl(i+1)-e1tbl(i))
     &                 /(atbl(i+1)-atbl(i))*atbl(i)
            ccofe(i,2)=(e1tbl(i+1)-e1tbl(i))/(atbl(i+1)-atbl(i))
            ccofe(i,3)=e2tbl(i)-(e2tbl(i+1)-e2tbl(i))
     &                 /(atbl(i+1)-atbl(i))*atbl(i)
            ccofe(i,4)=(e2tbl(i+1)-e2tbl(i))/(atbl(i+1)-atbl(i))
            ccofe(i,5)=e3tbl(i)-(e3tbl(i+1)-e3tbl(i))
     &                 /(atbl(i+1)-atbl(i))*atbl(i)
            ccofe(i,6)=(e3tbl(i+1)-e3tbl(i))/(atbl(i+1)-atbl(i))
         end do
      end if
c
 300  format(i10)
 310  format(1p,6e13.6)
 320  format(i5)
 400  format(4e12.4)
 510  format(a80)
c* set constatns
      alp =alpha*acos(-1.0)/180.0
c* set constants related to chemistry
c amass(is) is molecular weight in kg/mole
      amass(1) =0.016
      amass(2) =0.0140067
      amass(3) =0.0300067
      amass(4) =0.032
      amass(5) =0.0280134
      amas1(1) =1.0/amass(1)
      amas1(2) =1.0/amass(2)
      amas1(3) =1.0/amass(3)
      amas1(4) =1.0/amass(4)
      amas1(5) =1.0/amass(5)
c cvi(is) is specific heat at constant volume in j/mole
      cvi(1)   =12.472
      cvi(2)   =12.472
      cvi(3)   =20.786
      cvi(4)   =20.786
      cvi(5)   =20.786
c h0i(is) is formation energy in j/mole
      h0i(1)   =246814.0
      h0i(2)   =470700.0
      h0i(3)   = 89789.0
      gasc     =     8.314
      r1inf    =rinf*x1inf
      r2inf    =rinf*x2inf
      r3inf    =rinf*x3inf
      r4inf    =rinf*x4inf
      r5inf    =rinf*x5inf
      ano      =r1inf*amas1(1)
      ann      =r2inf*amas1(2)
      anno     =r3inf*amas1(3)
      ano2     =r4inf*amas1(4)
      ann2     =r5inf*amas1(5)
      rhor     =(ano+ann+anno+ano2+ann2)*gasc
      gascs    =rhor/rinf
      pinf     =rhor*ttinf
      cpbyr    =2.5*(ano+ann)+3.5*(anno+ano2+ann2)
      cvbyr    =1.5*(ano+ann)+2.5*(anno+ano2+ann2)
      gaminf   =cpbyr/cvbyr
      cinf     =sqrt(gaminf*gascs*ttinf)
      xminf    =fsv/cinf
      uinf     =fsv*cos(alp)
      vinf     =fsv*sin(alp)
      ruinf    =rinf*uinf
      rvinf    =rinf*vinf
      tv1      =1.0/tvinf
      xa       =exp(  226.0*tv1)
      xb       =exp(-3785.5*tv1)
      xa2      =xa*xa
      xa4      =xa2*xa2
      xa8      =xa4*xa4
      xa10     =xa8*xa2
      xa12     =xa10*xa2
      xa15     =xa12*xa2*xa
      xb2      =xb*xb
      xb3      =xb2*xb
      xb4      =xb2*xb2
      xb5      =xb4*xb
      xb6      =xb5*xb
      xb7      =xb6*xb
      ev1      =2260.0/(xa10-1.0)
      ev2      =3390.0/(xa15-1.0)
      ev3      =2712.0/(xa12-1.0)
      ev4      =(22713.0*xb3+18927.5*xb5)/(3.0+2.0*xb3+xb5)
      ev5      =113565.0*xb6/(9.0+ 5.0*xb6)
      ev6      =264985.0*xb7/(4.0+10.0*xb7)
      evinf    =8.314*(ano2*ev1+ann2*ev2+anno*ev3
     &                +ano2*ev4+ ano*ev5+ ann*ev6)
      ehinf    =h0i(1)*ano+h0i(2)*ann+h0i(3)*anno
      ekinf    =0.5*rinf*(uinf**2+vinf**2)
      cvav     =ano*cvi(1)+ann*cvi(2)+anno*cvi(3)+ano2*cvi(4)
     &                                           +ann2*cvi(5)
      einf     =cvav*ttinf+evinf+ehinf+ekinf
c transport coefficients
      sc       =0.5
c
      rs1      =gasc*amas1(1)
      cvts1    =1.5d0*rs1
      cvrs1    =0.0d0
      cvvs1    =0.0d0
      rs2      =gasc*amas1(2)
      cvts2    =1.5d0*rs2
      cvrs2    =0.0d0
      cvvs2    =0.0d0
      rs3      =gasc*amas1(3)
      cvts3    =1.5d0*rs3
      cvrs3    =rs3
      cvvs3    =rs3
      rs4      =gasc*amas1(4)
      cvts4    =1.5d0*rs4
      cvrs4    =rs4
      cvvs4    =rs4
      rs5      =gasc*amas1(5)
      cvts5    =1.5d0*rs5
      cvrs5    =rs5
      cvvs5    =rs5
c Wilke/Blottner/Eucken
      if(nvis.eq.0) then
         am       =1.0/(x1inf*amas1(1)+x2inf*amas1(2)+x3inf*amas1(3)
     &                 +x4inf*amas1(4)+x5inf*amas1(5))
         alnt     =log(ttinf)
c
         tv1      =1./tvinf
         xa       =exp(226.0*tv1)
         xa2      =xa*xa
         xa4      =xa2*xa2
         xa8      =xa4*xa4
         xa10     =xa8*xa2
         xa12     =xa10*xa2
         xa15     =xa12*xa2*xa
         ev1     =2260.0/(xa10-1.0)
         ev2     =3390.0/(xa15-1.0)
         ev3     =2712.0/(xa12-1.0)
c s1(O)
         x1       =x1inf*am*amas1(1)
         vmus1    =0.1*exp((a1*alnt+b1)*alnt+c1)
         vkas1    =vmus1*(2.5*cvts1+cvrs1)
         vkvs1    =vmus1*cvvs1
c s2(N)
         x2       =x2inf*am*amas1(2)
         vmus2    =0.1*exp((a2*alnt+b2)*alnt+c2)
         vkas2    =vmus2*(2.5*cvts2+cvrs2)
         vkvs2    =vmus2*cvvs2
c s3(NO)
         x3       =x3inf*am*amas1(3)
         vmus3    =0.1*exp((a3*alnt+b3)*alnt+c3)
         vkas3    =vmus3*(2.5*cvts3+cvrs3)
         vkvs3    =vmus3*cvvs3
c s4(O2)
         x4       =x4inf*am*amas1(4)
         vmus4    =0.1*exp((a4*alnt+b4)*alnt+c4)
         vkas4    =vmus4*(2.5*cvts4+cvrs4)
         vkvs4    =vmus4*cvvs4
c s5(N2)
         x5       =x5inf*am*amas1(5)
         vmus5    =0.1*exp((a5*alnt+b5)*alnt+c5)
         vkas5    =vmus5*(2.5*cvts5+cvrs5)
         vkvs5    =vmus5*cvvs5
c mixing (1-*)
         phs11t   =1.0+sqrt(vmus1/vmus1)*(amass(1)*amas1(1))**0.25
         phs11    =phs11t*phs11t/sqrt(8.0*(1.0+amass(1)*amas1(1)))
         phs12t   =1.0+sqrt(vmus1/vmus2)*(amass(2)*amas1(1))**0.25
         phs12    =phs12t*phs12t/sqrt(8.0*(1.0+amass(1)*amas1(2)))
         phs13t   =1.0+sqrt(vmus1/vmus3)*(amass(3)*amas1(1))**0.25
         phs13    =phs13t*phs13t/sqrt(8.0*(1.0+amass(1)*amas1(3)))
         phs14t   =1.0+sqrt(vmus1/vmus4)*(amass(4)*amas1(1))**0.25
         phs14    =phs14t*phs14t/sqrt(8.0*(1.0+amass(1)*amas1(4)))
         phs15t   =1.0+sqrt(vmus1/vmus5)*(amass(5)*amas1(1))**0.25
         phs15    =phs15t*phs15t/sqrt(8.0*(1.0+amass(1)*amas1(5)))
         phis1    =1.0/(x1*phs11+x2*phs12+x3*phs13+x4*phs14+x5*phs15)
c mixing (2-*)
         phs21t   =1.0+sqrt(vmus2/vmus1)*(amass(1)*amas1(2))**0.25
         phs21    =phs21t*phs21t/sqrt(8.0*(1.0+amass(2)*amas1(1)))
         phs22t   =1.0+sqrt(vmus2/vmus2)*(amass(2)*amas1(2))**0.25
         phs22    =phs22t*phs22t/sqrt(8.0*(1.0+amass(2)*amas1(2)))
         phs23t   =1.0+sqrt(vmus2/vmus3)*(amass(3)*amas1(2))**0.25
         phs23    =phs23t*phs23t/sqrt(8.0*(1.0+amass(2)*amas1(3)))
         phs24t   =1.0+sqrt(vmus2/vmus4)*(amass(4)*amas1(2))**0.25
         phs24    =phs24t*phs24t/sqrt(8.0*(1.0+amass(2)*amas1(4)))
         phs25t   =1.0+sqrt(vmus2/vmus5)*(amass(5)*amas1(2))**0.25
         phs25    =phs25t*phs25t/sqrt(8.0*(1.0+amass(2)*amas1(5)))
         phis2    =1.0/(x1*phs21+x2*phs22+x3*phs23+x4*phs24+x5*phs25)
c mixing (3-*)
         phs31t   =1.0+sqrt(vmus3/vmus1)*(amass(1)*amas1(3))**0.25
         phs31    =phs31t*phs31t/sqrt(8.0*(1.0+amass(3)*amas1(1)))
         phs32t   =1.0+sqrt(vmus3/vmus2)*(amass(2)*amas1(3))**0.25
         phs32    =phs32t*phs32t/sqrt(8.0*(1.0+amass(3)*amas1(2)))
         phs33t   =1.0+sqrt(vmus3/vmus3)*(amass(3)*amas1(3))**0.25
         phs33    =phs33t*phs33t/sqrt(8.0*(1.0+amass(3)*amas1(3)))
         phs34t   =1.0+sqrt(vmus3/vmus4)*(amass(4)*amas1(3))**0.25
         phs34    =phs34t*phs34t/sqrt(8.0*(1.0+amass(3)*amas1(4)))
         phs35t   =1.0+sqrt(vmus3/vmus5)*(amass(5)*amas1(3))**0.25
         phs35    =phs35t*phs35t/sqrt(8.0*(1.0+amass(3)*amas1(5)))
         phis3    =1.0/(x1*phs31+x2*phs32+x3*phs33+x4*phs34+x5*phs35)
c mixing (4-*)
         phs41t   =1.0+sqrt(vmus4/vmus1)*(amass(1)*amas1(4))**0.25
         phs41    =phs41t*phs41t/sqrt(8.0*(1.0+amass(4)*amas1(1)))
         phs42t   =1.0+sqrt(vmus4/vmus2)*(amass(2)*amas1(4))**0.25
         phs42    =phs42t*phs42t/sqrt(8.0*(1.0+amass(4)*amas1(2)))
         phs43t   =1.0+sqrt(vmus4/vmus3)*(amass(3)*amas1(4))**0.25
         phs43    =phs43t*phs43t/sqrt(8.0*(1.0+amass(4)*amas1(3)))
         phs44t   =1.0+sqrt(vmus4/vmus4)*(amass(4)*amas1(4))**0.25
         phs44    =phs44t*phs44t/sqrt(8.0*(1.0+amass(4)*amas1(4)))
         phs45t   =1.0+sqrt(vmus4/vmus5)*(amass(5)*amas1(4))**0.25
         phs45    =phs45t*phs45t/sqrt(8.0*(1.0+amass(4)*amas1(5)))
         phis4    =1.0/(x1*phs41+x2*phs42+x3*phs43+x4*phs44+x5*phs45)
c mixing (5-*)
         phs51t   =1.0+sqrt(vmus5/vmus1)*(amass(1)*amas1(5))**0.25
         phs51    =phs51t*phs51t/sqrt(8.0*(1.0+amass(5)*amas1(1)))
         phs52t   =1.0+sqrt(vmus5/vmus2)*(amass(2)*amas1(5))**0.25
         phs52    =phs52t*phs52t/sqrt(8.0*(1.0+amass(5)*amas1(2)))
         phs53t   =1.0+sqrt(vmus5/vmus3)*(amass(3)*amas1(5))**0.25
         phs53    =phs53t*phs53t/sqrt(8.0*(1.0+amass(5)*amas1(3)))
         phs54t   =1.0+sqrt(vmus5/vmus4)*(amass(4)*amas1(5))**0.25
         phs54    =phs54t*phs54t/sqrt(8.0*(1.0+amass(5)*amas1(4)))
         phs55t   =1.0+sqrt(vmus5/vmus5)*(amass(5)*amas1(5))**0.25
         phs55    =phs55t*phs55t/sqrt(8.0*(1.0+amass(5)*amas1(5)))
         phis5    =1.0/(x1*phs51+x2*phs52+x3*phs53+x4*phs54+x5*phs55)
c combined values
         vmueinf  =x1*vmus1*phis1+x2*vmus2*phis2+x3*vmus3*phis3
     &            +x4*vmus4*phis4+x5*vmus5*phis5
         vkapinf  =x1*vkas1*phis1+x2*vkas2*phis2+x3*vkas3*phis3
     &            +x4*vkas4*phis4+x5*vkas5*phis5
         vkavinf  =x1*vkvs1*phis1+x2*vkvs2*phis2+x3*vkvs3*phis3
     &            +x4*vkvs4*phis4+x5*vkvs5*phis5
c diffusion coefficient
         ds1inf   =vmueinf/(rinf*sc)
         ds2inf   =ds1inf
         ds3inf   =ds1inf
         ds4inf   =ds1inf
         ds5inf   =ds1inf
      end if
c Collisonal integrals by Gupta curv fitting
c Diffusion coefficent is determined by binary diffusion coefficient
c 1:O, 2:N, 3:NO, 4:O2, 5:N2
      if(nvis.eq.1) then
         pi    = dacos(-1.0d0)
         avoga1= 1.d0/6.02214d23 ! avogadro constant  [1/mol]
         boltz = 1.3806503d-23  ! Boltzmann constant [m2kg/s2K]
         dct1  = 8.d0/3.d0      ! constant for delta_sr(1)
         dct2  = 16.d0/5.d0     ! constant for delta_sr(2)
         te    =ttinf
         ttv   =tvinf
         pp    =pinf
         alte  =dlog(te)
         alt2  =alte*alte
         prt1  =1.d0/pi/te/gasc
         altv  =dlog(ttv)
         altv2 =altv*altv
         prtv  =1.d0/pi/ttv/gasc
         xcinf(1)=x1inf
         xcinf(2)=x2inf
         xcinf(3)=x3inf
         xcinf(4)=x4inf
         xcinf(5)=x5inf
         do n=1,5
            xr(n)=xcinf(n)*amas1(n)
         end do
         xt=0.d0
         do n=1,5
            xt=xt+xr(n)
         end do
c pi*Omega_sr(1,1)(=pomega1) and Omega_sr(2,2)(=pomega2)
c unit A^2(=10^{-20}[m])
         do m=1,5               ! s
            do l=1,5            ! r
               pomega1   =dexp(dog(l,m,1))*te**(aog(l,m,1)*alt2
     &                        +bog(l,m,1)*alte+cog(l,m,1))
               pomega2   =dexp(dog(l,m,2))*te**(aog(l,m,2)*alt2
     &                        +bog(l,m,2)*alte+cog(l,m,2))
               pomegav   =dexp(dog(l,m,1))*ttv**(aog(l,m,1)
     &                  *altv2+bog(l,m,1)*altv+cog(l,m,1))
c delta_sr^1(=dlsr1) and delta_sr^2(=dlsr2)
               dlsr1(l,m)=dct1*dsqrt(2.d0*amass(m)*amass(l)
     &                   /(amass(m)+amass(l))*prt1)*pomega1*1.d-20
               dlsr2(l,m)=dct2*dsqrt(2.d0*amass(m)*amass(l)
     &                   /(amass(m)+amass(l))*prt1)*pomega2*1.d-20
               dlsrv(l,m)=dct1*dsqrt(2.d0*amass(m)*amass(l)
     &                   /(amass(m)+amass(l))*prtv)*pomegav*1.d-20
c alpha_sr^1(=alsr)
               ramass    =amass(m)/amass(l)
               alsr(l,m) = 1.d0 + ((1.d0-ramass)*(0.45d0-2.54d0*ramass)
     &                  /((1.d0+ramass)**2))
               ddsr(l,m)  = (boltz*te)/(pp*dlsr1(l,m))
            end do
         end do
c
         xxdl1(:)=0.d0
         xxdl2(:)=0.d0
         xxdlv(:)=0.d0
         axdl2(:)=0.d0
         do m = 1,5             ! s
            do l = 1,5          ! r
               xxdl1(m) = xxdl1(m)+          xr(l)*dlsr1(l,m)
               xxdl2(m) = xxdl2(m)+          xr(l)*dlsr2(l,m)
               xxdlv(m) = xxdlv(m)+          xr(l)*dlsrv(l,m)
               axdl2(m) = axdl2(m)+alsr(l,m)*xr(l)*dlsr2(l,m)
            end do
         end do
         vics=0.d0
         thct=0.d0
         thcr=0.d0
         thcv=0.d0
         do m=1,5               !O,N,NO,O2,N2
            vics=vics+xcinf(m)/xxdl2(m)
            thct=thct+xr(m)/axdl2(m)
         end do
         do m=3,5               ! NO,O2,N2
            thcr=thcr+xr(m)/xxdl1(m)
            thcv=thcv+xr(m)/xxdlv(m)
         end do
         vmueinf=vics*avoga1
         vkapinf=boltz*(3.75d0*thct+thcr)
         vkavinf=boltz*thcv
c binary diffusion coefficient
         dsb(:)=0.d0
         do m=1,5
            do l=1,5
               if(m.ne.l) then
                  dsb(m) = dsb(m) + xr(l)/ddsr(l,m)
               end if
            end do
         end do
         ds1inf = ((xt**2)*amass(1)*(1.d0-xcinf(1)))/dsb(1)
         ds2inf = ((xt**2)*amass(2)*(1.d0-xcinf(2)))/dsb(2)
         ds3inf = ((xt**2)*amass(3)*(1.d0-xcinf(3)))/dsb(3)
         ds4inf = ((xt**2)*amass(4)*(1.d0-xcinf(4)))/dsb(4)
         ds5inf = ((xt**2)*amass(5)*(1.d0-xcinf(5)))/dsb(5)
      end if
c
      write(6,*) 'vmueinf=',vmueinf
      write(6,*) 'vkapinf=',vkapinf
      write(6,*) 'vkavinf=',vkavinf
c transport coefficient for vibration-electronic energy
      ev1g      =8.314*ev1
      ev2g      =8.314*ev2
      ev3g      =8.314*ev3
      ev4g      =8.314*ev4
      ev5g      =8.314*ev5
      ev6g      =8.314*ev6
      evs1inf   =      ev5g *amas1(1)
      evs2inf   =      ev6g *amas1(2)
      evs3inf   =      ev3g *amas1(3)
      evs4inf   =(ev1g+ev4g)*amas1(4)
      evs5inf   =      ev2g *amas1(5)
      hs1inf    =(cvts1+cvrs1+rs1)*ttinf+(h0i(1)+ev5g)*amas1(1)
      hs2inf    =(cvts2+cvrs2+rs2)*ttinf+(h0i(2)+ev6g)*amas1(2)
      hs3inf    =(cvts3+cvrs3+rs3)*ttinf+(h0i(3)+ev3g)*amas1(3)
      hs4inf    =(cvts4+cvrs4+rs4)*ttinf+(  ev1g+ev4g)*amas1(4)
      hs5inf    =(cvts5+cvrs5+rs5)*ttinf+        ev2g *amas1(5)
c* coefficients
      fr        =(ano+anno+2.0*ano2)/(ann+anno+2.0*ann2)
      alp1      =1.0/(amass(4)+amass(5)/fr)
      alp2      =-alp1*(1.0+amass(5)/(2.0*fr*amass(1)))
      alp3      =-alp1*(1.0+0.5*(1.0/fr-1.0)*amass(5)/amass(3))
      alp4      =1.0/(amass(5)+fr*amass(4))
      alp5      =-alp4*(1.0+0.5*fr*amass(4)/amass(2))
      alp6      =-alp4*(1.0-0.5*(1.0-fr)*amass(4)/amass(3))
      bet1      =amass(4)*alp1
      bet2      =amass(4)*alp2
      bet3      =amass(4)*alp3
      bet4      =amass(5)*alp4
      bet5      =amass(5)*alp5
      bet6      =amass(5)*alp6
c* display setup information
      write(6,*) ' '
      write(6,*) '*** procedures ***'
      if(istart.eq.0) then
         write(6,*) 'parm) calculation starts from initial state.'
      else
         write(6,*) 'parm) calculation starts from the disk file.'
         write(6,*) '           file name of field data=',file_field
      endif
      write(6,*) ' '
      if(nprc.eq.0) then
         write(6,*) 'parm) euler simulation is chosen'
      else if (nprc.eq.1) then
         write(6,*) 'parm) NS(laminar) simulation is chosen'
      else
         write(6,*) 'parm) NS(RANS(Baldwin-Lomax)) simulation is chosen'
      end if
      if(nchm.eq.0) then
         write(6,*) 'parm) frozen flow calculation initiated'
      else
         write(6,*) 'parm) nonequilibrium reactive flow calculation'
      endif
      if(iprc.eq.0) then
         write(6,*) 'parm) implicit time integration is chosen'
      else
         write(6,*) 'parm) explicit time integration is chosen'
      endif
      if(nacu.eq.1) then
         write(6,*) 'parm) first order accuracy in space'
      else
         write(6,*) 'parm) higher order accuracy in space'
      endif
      write(6,*) ' '
      write(6,*) '*** run_time parameters ***'
      write(6,*) '           nstp  = ',nstp
      write(6,*) '           nout  = ',nout
      write(6,*) '           mdd   = ',mdd
      if(lodt.eq.0) then
         write(6,*) 'parm) global time stepping is employed'
         write(6,*) '           cflmx  = ',cflmx
         write(6,*) '           grwmx  = ',grwmx
      else
         write(6,*) 'parm) local time stepping is employed'
         write(6,*) '           cflmx  = ',cflmx
         write(6,*) '           grwmx  = ',grwmx
      endif
      write(6,*) ' '
      write(6,*) '*** free stream conditions ***'
      write(6,*) '   fsv=',     fsv,'[m/s]'
      write(6,*) ' xminf=',   xminf
      write(6,*) '  rinf=',    rinf,'[Kg/m3]'
      write(6,*) '  uinf=',    uinf,'[m/s]'
      write(6,*) '  vinf=',    vinf,'[m/s]'
      write(6,*) '  pinf=',    pinf,'[N/m2]'
      write(6,*) '  einf=',    einf,'[J/m3]'
      write(6,*) ' evinf=',   evinf,'[J/m3]'
      write(6,*) ' ttinf=',   ttinf,'[K]'
      write(6,*) ' tvinf=',   tvinf,'[K]'
      write(6,*) '  cinf=',    cinf,'[m/s]'
      write(6,*) 'gaminf=',  gaminf
      write(6,*) ' x1inf=',   x1inf
      write(6,*) ' x2inf=',   x2inf
      write(6,*) ' x3inf=',   x3inf
      write(6,*) ' x4inf=',   x4inf
      write(6,*) ' x5inf=',   x5inf
      write(6,*) ' r1inf=',   r1inf,'[Kg/m3]'
      write(6,*) ' r2inf=',   r2inf,'[Kg/m3]'
      write(6,*) ' r3inf=',   r3inf,'[Kg/m3]'
      write(6,*) ' r4inf=',   r4inf,'[Kg/m3]'
      write(6,*) ' r5inf=',   r5inf,'[Kg/m3]'
      write(6,*) ' '
      write(6,*) '*** wall boundary condition ***'
      if(twall.le.0.0) then
         write(6,*) 'parm) adiabatic wall condition is chosen'
      else
         write(6,*) 'parm) isothermal wall condition is chosen'
         write(6,*) '             wall temperature =',twall,'[K]'
      endif
      write(6,*) ' surface catalytic efficicency for N =',cefn
      write(6,*) ' surface catalytic efficicency for O =',cefo
      write(6,*) ''
      write(6,*) '*** radiation coupling      ***'
      if(irad.eq.0) then
         write(6,*) 'parm) non coupling between radiation and flow '
      else if(irad.eq.1) then
         write(6,*) 'parm) coupling between radiation and flow'
      else
         write(6,*) 'parm) output data for lignbylign (non coupling)'
      end if
      if(nvis.eq.0) then
         write(6,*) 'Wilke,Blottner,Eucken is chosen for tranport calc.'
      else
         write(6,*) 'Gupta curve fitting is chosen for tranport calc.'
      end if
      if(nrlx.eq.0) then
         write(6,*) 'Park bridging model is OFF'
      else
         write(6,*) 'Park bridging model is ON'
      end if
c* termination
      return
      end
c
      subroutine metr
c     **************************************************************
c     *     metrices of cell system                                *
c     **************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /coor3/ xg(0:nx,0:ny),yg(0:nx,0:ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
c* area vector (i=const)
      do j=1,jmax-1
         do i=1,imax
            x1        =x(i,j)
            y1        =y(i,j)
            x2        =x(i,j+1)
            y2        =y(i,j+1)
            si(i,j,1) =  y2-y1
            si(i,j,2) =-(x2-x1)
            si(i,j,3) =sqrt(si(i,j,1)**2+si(i,j,2)**2)
         end do
      end do
c* area vector (j=const)
      do j=1,jmax
         do i=1,imax-1
            x1        =x(i  ,j)
            y1        =y(i  ,j)
            x2        =x(i+1,j)
            y2        =y(i+1,j)
            sj(i,j,1) =-(y2-y1)
            sj(i,j,2) =  x2-x1
            sj(i,j,3) =sqrt(sj(i,j,1)**2+sj(i,j,2)**2)
         end do
      end do
c* cell volume / cell centroid
      do j=1,jmax-1
         do i=1,imax-1
            xg1      =(x(i,j)+x(i+1,j)+x(i+1,j+1))/3.0
            yg1      =(y(i,j)+y(i+1,j)+y(i+1,j+1))/3.0
            sg1      =0.5*((x(i+1,j)-x(i,j))*(y(i+1,j+1)-y(i,j))
     &                    -(y(i+1,j)-y(i,j))*(x(i+1,j+1)-x(i,j)))
            xg2      =(x(i,j)+x(i,j+1)+x(i+1,j+1))/3.0
            yg2      =(y(i,j)+y(i,j+1)+y(i+1,j+1))/3.0
            sg2      =0.5*((x(i+1,j+1)-x(i,j))*(y(i,j+1)-y(i,j))
     &                    -(y(i+1,j+1)-y(i,j))*(x(i,j+1)-x(i,j)))
            vol(i,j) =sg1+sg2
            xg(i,j)  =(sg1*xg1+sg2*xg2)/(sg1+sg2)
            yg(i,j)  =(sg1*yg1+sg2*yg2)/(sg1+sg2)
         end do
      end do
      do j=1,jmax-1
         vol(0   ,j) =vol(1     ,j)
         vol(imax,j) =vol(imax-1,j)
      end do
      do i=1,imax-1
         vol(i,0   ) =vol(i,1     )
         vol(i,jmax) =vol(i,jmax-1)
      end do
c* termination
      return
      end
c
      subroutine wdst
c     ************************************************************
c     *     set distance from wall                               *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /coor3/ xg(0:nx,0:ny),yg(0:nx,0:ny)
      common /coor4/ sl(0:nx,0:ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
c* reset sl array
      do j=1,jmax
         do i=1,imax
            sl(i,j)=1.0d5
         end do
      end do
c*sawada ver.
c brute force attack
!      do m=1,imax-1
!         xwl =        x(m,1)
!         ywl =        y(m,1)
!         xwm = 0.5d0*(x(m,1)+x(m+1,1))
!         ywm = 0.5d0*(y(m,1)+y(m+1,1))
!         xwr =               x(m+1,1)
!         ywr =               y(m+1,1)
!         do j=1,jmax-1
!            do i=1,imax-1
!               ds2l    = (xwl-xg(i,j))**2+(ywl-yg(i,j))**2
!               ds2     = (xwm-xg(i,j))**2+(ywm-yg(i,j))**2
!               ds2r    = (xwr-xg(i,j))**2+(ywr-yg(i,j))**2
!               sl(i,j) = dmin1(sl(i,j),ds2,ds2l,ds2r)
!            end do
!         end do
!      end do
!      do j=1,jmax-1
!         do i=1,imax-1
!            sl(i,j)=dsqrt(sl(i,j))
!         end do
!      end do
c*furudate ver.
      do l=1,imax-1
         do j=1,jmax-1
            do i=1,imax-1
               aex      =sj(l,1,1)/sj(l,1,3)
               aey      =sj(l,1,2)/sj(l,1,3)
               aes      =sj(l,1,3)
               x1       =xg(i,j)-x(l,1)
               y1       =yg(i,j)-y(l,1)
               dd1      =sqrt(x1**2+y1**2)
               x2       =xg(i,j)-x(l+1,1)
               y2       =yg(i,j)-y(l+1,1)
               dd2      =sqrt(x2**2+y2**2)
               dd       =abs(x1*aex+y1*aey)
               d1d      =abs(x1*aey-y1*aex)
               d2d      =abs(x2*aey-y2*aex)
               if (d1d.gt.aes) dd =1.0d5
               if (d2d.gt.aes) dd =1.0d5
               sl(i,j) =min(dd1,dd2,dd,sl(i,j))
            end do
         end do
      end do
c* termination
      return
      end
c
      subroutine init
c     ************************************************************
c     *     set initial condition                                *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /cntl1/ istart
      common /cntl5/ nprc,iprc
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /fscd3/ vmueinf,vkapinf,vkavinf,
     &               ds1inf,ds2inf,ds3inf,ds4inf,ds5inf,
     &               evs1inf,evs2inf,evs3inf,evs4inf,evs5inf,
     &               hs1inf,hs2inf,hs3inf,hs4inf,hs5inf
      common /stat1/ loop,locl
      common /stat2/ time
c* set initial condition (if istart.eq.0)
      do j=0,jmax
         do i=0,imax
            qs(i,j,1)  =r1inf
            qs(i,j,2)  =r2inf
            qs(i,j,3)  =r3inf
            qv(i,j,1)  =evinf
            qc(i,j,1)  =rinf
            qc(i,j,2)  =ruinf
            qc(i,j,3)  =rvinf
            qc(i,j,4)  =einf
            pc(i,j)    =pinf
            cc(i,j)    =cinf
            tt(i,j)    =ttinf
            tv(i,j)    =tvinf
            uc(i,j)    =uinf
            vc(i,j)    =vinf
            xc(i,j, 1) =x1inf
            xc(i,j, 2) =x2inf
            xc(i,j, 3) =x3inf
            xc(i,j, 4) =x4inf
            xc(i,j, 5) =x5inf
         end do
      end do
      if(istart.eq.0) then
         loop  =0
         time  =0.0
      endif
c* reset transport coefficients
      if(nprc.eq.0) then
         do j=0,jmax
            do i=0,imax
               vmue(i,j)  =0.0
               vkap(i,j)  =0.0
               vkav(i,j)  =0.0
               ds(i,j,1)  =0.0
               ds(i,j,2)  =0.0
               ds(i,j,3)  =0.0
               ds(i,j,4)  =0.0
               ds(i,j,5)  =0.0
               evs(i,j,1) =0.0
               evs(i,j,2) =0.0
               evs(i,j,3) =0.0
               evs(i,j,4) =0.0
               evs(i,j,5) =0.0
               hs(i,j,1)  =0.0
               hs(i,j,2)  =0.0
               hs(i,j,3)  =0.0
               hs(i,j,4)  =0.0
               hs(i,j,5)  =0.0
            end do
         end do
      else
         do j=0,jmax
            do i=0,imax
               vmue(i,j)  =vmueinf
               vkap(i,j)  =vkapinf
               vkav(i,j)  =vkavinf
               ds(i,j,1)  =ds1inf
               ds(i,j,2)  =ds2inf
               ds(i,j,3)  =ds3inf
               ds(i,j,4)  =ds4inf
               ds(i,j,5)  =ds5inf
               evs(i,j,1) =evs1inf
               evs(i,j,2) =evs2inf
               evs(i,j,3) =evs3inf
               evs(i,j,4) =evs4inf
               evs(i,j,5) =evs5inf
               hs(i,j,1)  =hs1inf
               hs(i,j,2)  =hs2inf
               hs(i,j,3)  =hs3inf
               hs(i,j,4)  =hs4inf
               hs(i,j,5)  =hs5inf
            end do
         end do
      endif
c* termination
      return
      end
c
      subroutine dtsz
c     ************************************************************
c     *     step_size estimation                                 *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /aray4/ dcq(nx,ny,4),dsq(nx,ny,3),dvq(0:nx,0:ny,1)
      common /coor1/ imax,jmax
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /dtst1/ dt(nx,ny),dt0,dti
      common /cntl2/ cflmx,grwmx
      common /cntl4/ lodt
      common /cntl5/ nprc,iprc
      common /parm3/ mdd
      common /stat1/ loop,locl
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
c* print control
      mw  =mod(loop,mdd)
      if(loop.le.100) mw =0
c* step size estimation
      dt0 =1.0e30
      dti =1.0e30
      do j=1,jmax-1
         do i=1,imax-1
            volv    =1.0/vol(i,j)
            coef    =vmue(i,j)/(qc(i,j,1)*vol(i,j))
            qx      =si(i,j,1)*uc(i,j)+si(i,j,2)*vc(i,j)
            dtx     =vol(i,j)/(abs(qx)
     &                        +(cc(i,j)+coef*si(i,j,3))*si(i,j,3))
            qy      =sj(i,j,1)*uc(i,j)+sj(i,j,2)*vc(i,j)
            dty     =vol(i,j)/(abs(qy)
     &                        +(cc(i,j)+coef*sj(i,j,3))*sj(i,j,3))
            dt(i,j) =cflmx*dtx*dty/(dtx+dty)
            if(dt(i,j).lt.dt0) then
               dt0 =dt(i,j)
               i0  =i
               j0  =j
            endif
            dti    =min(dti,dtx)
         end do
      end do
c* growth control
      dtmn =1.0e20
      cof1 =amas1(1)*h0i(1)
      cof2 =amas1(2)*h0i(2)
      cof3 =amas1(3)*h0i(3)
      do j=1,jmax-1
         do i=1,imax-1
            qs4       =xc(i,j,4)*qc(i,j,1)
            qs5       =xc(i,j,5)*qc(i,j,1)
            ano       =qs(i,j,1)*amas1(1)
            ann       =qs(i,j,2)*amas1(2)
            anno      =qs(i,j,3)*amas1(3)
            ano2      =qs4      *amas1(4)
            ann2      =qs5      *amas1(5)
            ansum     =ano+ann+anno+ano2+ann2
            cvav      =ano*cvi(1)+ann*cvi(2)+anno*cvi(3)+ano2*cvi(4)
     &                                                  +ann2*cvi(5)
            rbar      =8.314*ansum/qc(i,j,1)
            cvf       =cvi(1)*xc(i,j,1)+cvi(2)*xc(i,j,2)
     &                +cvi(3)*xc(i,j,3)+cvi(4)*xc(i,j,4)
     &                +cvi(5)*xc(i,j,5)
            bet       =rbar/cvf
            thc       =0.5*(uc(i,j)*uc(i,j)+vc(i,j)*vc(i,j))
            gam1      =8.314*amas1(1)*tt(i,j)+bet*(thc-cvi(1)*tt(i,j)
     &                                            -h0i(1))
            gam2      =8.314*amas1(2)*tt(i,j)+bet*(thc-cvi(2)*tt(i,j)
     &                                            -h0i(2))
            gam3      =8.314*amas1(3)*tt(i,j)+bet*(thc-cvi(3)*tt(i,j)
     &                                            -h0i(3))
            gam4      =8.314*amas1(4)*tt(i,j)+bet*(thc-cvi(4)*tt(i,j))
            gam5      =8.314*amas1(5)*tt(i,j)+bet*(thc-cvi(5)*tt(i,j))
            dsq4      =bet1*dcq(i,j,1)+bet2*dsq(i,j,1)+bet3*dsq(i,j,3)
            dsq5      =bet4*dcq(i,j,1)+bet5*dsq(i,j,2)+bet6*dsq(i,j,3)
            delp      =(gam1*dsq(i,j,1)+gam2*dsq(i,j,2)
     &                 +gam3*dsq(i,j,3)+gam4*dsq4+gam5*dsq5
     &                 -bet*(uc(i,j)*dcq(i,j,2)+vc(i,j)*dcq(i,j,3)
     &                 +dvq(i,j,1)-dcq(i,j,4)))/vol(i,j)
            dt(i,j)   =min(dt(i,j),grwmx*pc(i,j)/(abs(delp)+1.0e-15),
     &                     dti)
            if(dt(i,j).lt.dtmn) then
               imn  =i
               jmn  =j
               dtmn =dt(i,j)
            endif
         end do
      end do
      cflt =cflmx*dtmn/dt0
      dt0  =dtmn
      if(lodt.eq.0) then
         do j=1,jmax-1
            do i=1,imax-1
               dt(i,j) =dtmn
            end do
         end do
      endif
c* termination
      return
      end
c
      subroutine intg
c     ************************************************************
c     *     integration control                                  *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /aray3/ qcl(nx,ny,4),qcr(nx,ny,4),
     &               qsl(nx,ny,3),qsr(nx,ny,3),
     &               qvl(nx,ny,1),qvr(nx,ny,1)
      common /aray4/ dcq(nx,ny,4),dsq(nx,ny,3),dvq(0:nx,0:ny,1)
      common /aray7/ scs(nx,ny,3),scv(nx,ny,1)
      common /aray9/ tvev(nx,0:ny)
      common /coor1/ imax,jmax
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /dtst1/ dt(nx,ny),dt0,dti
      common /parm3/ mdd
      common /cntl2/ cflmx,grwmx
      common /cntl4/ lodt
      common /cntl5/ nprc,iprc
      common /cntl8/ nchm
      common /cntl9/ qq,fmul3,fmul4
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /stat1/ loop,locl
      common /stat2/ time
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /rms/    rmshflx   !
      dimension aj(ny),bj(ny),cj(ny),dj(ny),ej(ny)

c* disp control
      mw =mod(loop,mdd)
      if(loop.le.100) mw=0
c* increment
      rms  =0.0
      rmx  =0.0
c      do j=1,jmax-1
c         do i=1,imax-1
c            rmt =abs(dcq(i,j,1))/vol(i,j)
c            rms =rms+rmt*rmt
c            rmx =max(rmx,rmt)
c         end do
c      end do
      do j=1,jmax-1
         do i=1,imax-1
            dtt        =dt(i,j)/vol(i,j)
            scs(i,j,1) =dt(i,j)*scs(i,j,1)
            scs(i,j,2) =dt(i,j)*scs(i,j,2)
            scs(i,j,3) =dt(i,j)*scs(i,j,3)
            scv(i,j,1) =dt(i,j)*scv(i,j,1)
            dsq(i,j,1) =dtt*dsq(i,j,1)+dt(i,j)*qsr(i,j,1)
            dsq(i,j,2) =dtt*dsq(i,j,2)+dt(i,j)*qsr(i,j,2)
            dsq(i,j,3) =dtt*dsq(i,j,3)+dt(i,j)*qsr(i,j,3)
            dvq(i,j,1) =dtt*dvq(i,j,1)+dt(i,j)*qvr(i,j,1)
            dcq(i,j,1) =dtt*dcq(i,j,1)+dt(i,j)*qcr(i,j,1)
            dcq(i,j,2) =dtt*dcq(i,j,2)+dt(i,j)*qcr(i,j,2)
            dcq(i,j,3) =dtt*dcq(i,j,3)+dt(i,j)*qcr(i,j,3)
            dcq(i,j,4) =dtt*dcq(i,j,4)+dt(i,j)*qcr(i,j,4)
            dcq(i,j,4) =dcq(i,j,4)-dvq(i,j,1)
         end do
      end do
      do j=1,jmax-1
         do i=1,imax-1
            rmt =abs(dcq(i,j,1))/vol(i,j)
            rms =rms+rmt*rmt
            rmx =max(rmx,rmt)
         end do
      end do
c* implicit treatment for heat conduction in vibrational energy eq.
c in (j)
      do i=1,imax-1
         do j=1,jmax-1
            ds    =2.0*vol(i,j)/(sj(i,j,3)+sj(i,j+1,3))
            vkp   =0.5*(vkav(i,j+1)+vkav(i,j  ))
            vkm   =0.5*(vkav(i,j  )+vkav(i,j-1))
            cof   =dt(i,j)/(ds*ds)
            aj(j) =-cof*vkm*tvev(i,j-1)
            cj(j) =-cof*vkp*tvev(i,j+1)
            bj(j) =1.0+cof*(vkm+vkp)*tvev(i,j)
         end do
         aj(1     ) =0.0
         cj(jmax-1) =0.0
         dj(1     ) =bj(1)
         ej(1     ) =dvq(i,1,1)/dj(1)
         do j=2,jmax-1
            dj(j) =bj(j)-aj(j)*cj(j-1)/dj(j-1)
            ej(j) =(-aj(j)*ej(j-1)+dvq(i,j,1))/dj(j)
         end do
         dvq(i,jmax-1,1) =ej(jmax-1)
         do j=jmax-2,1,-1
            dvq(i,j,1) =ej(j)-dvq(i,j+1,1)*cj(j)/dj(j)
         end do
      end do
cc in (i)
c      do j=1,jmax-1
c         do i=1,imax-1
c            ds    =2.0*vol(i,j)/(si(i,j,3)+si(i+1,j,3))
c            vkp   =0.5*(vkav(i+1,j)+vkav(i  ,j))
c            vkm   =0.5*(vkav(i  ,j)+vkav(i-1,j))
c            cof   =dt(i,j)/(ds*ds)
c            aj(i) =-cof*vkm*tvev(i-1,j)
c            cj(i) =-cof*vkp*tvev(i+1,j)
c            bj(i) =1.0+cof*(vkm+vkp)*tvev(i,j)
c         end do
c         aj(1     ) =0.0
c         cj(imax-1) =0.0
c         dj(1     ) =bj(1)
c         ej(1     ) =dvq(1,j,1)/dj(1)
c         do i=2,imax-1
c            dj(i) =bj(i)-aj(i)*cj(i-1)/dj(i-1)
c            ej(i) =(-aj(i)*ej(i-1)+dvq(i,j,1))/dj(i)
c         end do
c         dvq(imax-1,j,1) =ej(imax-1)
c         do i=imax-2,1,-1
c            dvq(i,j,1) =ej(i)-dvq(i+1,j,1)*cj(i)/dj(i)
c         end do
c      end do
      do j=1,jmax-1
         do i=1,imax-1
            dcq(i,j,4) =dcq(i,j,4)+dvq(i,j,1)
         end do
      end do
c* display information
      rms =sqrt(rms/float((imax-1)*(jmax-1)))
      rmshflx  =dsqrt(rmshflx/dfloat(imax-1))
      time =time+dt0
      if(mw.eq.0)  then
         write(6,600) loop,time,dt0,rmx,rms,rmshflx
         write(333,601) loop,time,dt0,rmx,log10(rms),log10(rmshflx)
      end if
 600  format(' intg: ',i10,' ',f13.7,4(2x,e13.6))
 601  format(i10,' ',f13.7,4(2x,e13.6))
c
c* add source terms
      do j=1,jmax-1
         do i=1,imax-1
            dsq(i,j,1) =dsq(i,j,1)+scs(i,j,1)
            dsq(i,j,2) =dsq(i,j,2)+scs(i,j,2)
            dsq(i,j,3) =dsq(i,j,3)+scs(i,j,3)
            dvq(i,j,1) =dvq(i,j,1)+scv(i,j,1)
         end do
      end do
c* implicit scheme
      if(iprc.eq.0) call impc_pnt
c*
      if (nchm.eq.0)then
      do j=1,jmax-1
         do i=1,imax-1
           dsq(i,j,1)=0.0
           dsq(i,j,2)=0.0
           dsq(i,j,3)=0.0
           dvq(i,j,1)=0.0
         enddo
      enddo
      endif

c* set new values
      qs1mn = 1.0e20
      qs1mx =-0.1e20
      qs2mn = 1.0e20
      qs2mx =-0.1e20
      qs3mn = 1.0e20
      qs3mx =-0.1e20
      qs4mn = 1.0e20
      qs4mx =-0.1e20
      qs5mn = 1.0e20
      qs5mx =-0.1e20
      qv1mn = 1.0e20
      qv1mx =-0.1e20
      do j=1,jmax-1
         do i=1,imax-1
            qs(i,j,1) =qs(i,j,1)+dsq(i,j,1)
            qs(i,j,2) =qs(i,j,2)+dsq(i,j,2)
            qs(i,j,3) =qs(i,j,3)+dsq(i,j,3)
!            qv(i,j,1) =qv(i,j,1)+dvq(i,j,1)
!            qv(i,j,1)  =dmax1(evinf,qv(i,j,1)+dvq(i,j,1))
            qv(i,j,1) =dmax1(0.0d-30,qv(i,j,1)+dvq(i,j,1))
            qc(i,j,1) =qc(i,j,1)+dcq(i,j,1)
            qc(i,j,2) =qc(i,j,2)+dcq(i,j,2)
            qc(i,j,3) =qc(i,j,3)+dcq(i,j,3)
            qc(i,j,4) =qc(i,j,4)+dcq(i,j,4)
            qs1       =max(0.0d0,qs(i,j,1))
            qs2       =max(0.0d0,qs(i,j,2))
            qs3       =max(0.0d0,qs(i,j,3))
            qs4       =bet1*qc(i,j,1)+bet2*qs1+bet3*qs3
            qs5       =bet4*qc(i,j,1)+bet5*qs2+bet6*qs3
            qs(i,j,1) =qs1
            qs(i,j,2) =qs2
            qs(i,j,3) =qs3
            qs1mn =min(qs1mn,qs(i,j,1))
            qs1mx =max(qs1mx,qs(i,j,1))
            qs2mn =min(qs2mn,qs(i,j,2))
            qs2mx =max(qs2mx,qs(i,j,2))
            qs3mn =min(qs3mn,qs(i,j,3))
            qs3mx =max(qs3mx,qs(i,j,3))
            qs4mn =min(qs4mn,qs4      )
            qs4mx =max(qs4mx,qs4      )
            qs5mn =min(qs5mn,qs5      )
            qs5mx =max(qs5mx,qs5      )
            qv1mn =min(qv1mn,qv(i,j,1))
            qv1mx =max(qv1mx,qv(i,j,1))
         end do
      end do
c
      if (nchm.eq.0)then
      do j=1,jmax-1
         do i=1,imax-1
           qs(i,j,1) =xc(i,j,1)*qc(i,j,1)
           qs(i,j,2) =xc(i,j,2)*qc(i,j,1)
           qs(i,j,3) =xc(i,j,3)*qc(i,j,1)
c
           rv        =1.0/qc(i,j,1)
           qs4       =bet1*qc(i,j,1)+bet2*qs(i,j,1)+bet3*qs(i,j,3)
           qs5       =bet4*qc(i,j,1)+bet5*qs(i,j,2)+bet6*qs(i,j,3)
           ano       =qs(i,j,1)*amas1(1)
           ann       =qs(i,j,2)*amas1(2)
           anno      =qs(i,j,3)*amas1(3)
           ano2      =qs4      *amas1(4)
           ann2      =qs5      *amas1(5)
           tv1       =1.0/tv(i,j)
           xa       =exp(226.0*tv1)
           xb       =exp(-3785.5*tv1)
           xa2      =xa  *xa
           xa4      =xa2 *xa2
           xa8      =xa4 *xa4
           xa9      =xa8 *xa
           xa10     =xa8 *xa2
           xa11     =xa10*xa
           xa12     =xa10*xa2
           xa14     =xa12*xa2
           xa15     =xa12*xa2*xa
           xb2      =xb  *xb
           xb3      =xb2 *xb
           xb4      =xb2 *xb2
           xb5      =xb4 *xb
           xb6      =xb5 *xb
           xb7      =xb6 *xb
           ev1       =2260.0/(xa10-1.0)
           ev2       =3390.0/(xa15-1.0)
           ev3       =2712.0/(xa12-1.0)
           ev4       =(22713.0*xb3+18927.5*xb5)/(3.0+2.0*xb3+xb5)
           ev5       =113565.0*xb6/(9.0+5.0*xb6)
           ev6       =264985.0*xb7/(4.0+10.0*xb7)
           qv(i,j,1) =8.314*(ano2*ev1+ann2*ev2+anno*ev3
     &                      +ano2*ev4+ano *ev5+ann *ev6)
         enddo
      enddo
      endif
c
c* termination
      return
      end
      subroutine impc_pnt
c     ************************************************************
c     *     implicit for chemical source term                    *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray4/ dcq(nx,ny,4),dsq(nx,ny,3),dvq(0:nx,0:ny,1)
      common /aray8/ ajs(nx,ny,3),ajv(nx,ny,1)
      common /coor1/ imax,jmax
      common /dtst1/ dt(nx,ny),dt0,dti
c* forward sequence
      do j=1,jmax-1
         do i=1,imax-1
            dsq(i,j,1) =dsq(i,j,1)/(1.0+dt(i,j)*ajs(i,j,1))
            dsq(i,j,2) =dsq(i,j,2)/(1.0+dt(i,j)*ajs(i,j,2))
            dsq(i,j,3) =dsq(i,j,3)/(1.0+dt(i,j)*ajs(i,j,3))
            dvq(i,j,1) =dvq(i,j,1)/(1.0+dt(i,j)*ajv(i,j,1))
         end do
      end do
c* termination
      return
      end

      subroutine trpc
c     **************************************************************
c     *     transport coefficient                                  *
c     **************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /cntl10/nvis,nrlx
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
c Blottner's curve fitting parameters
      data a1/0.0203144/,b1/ 0.4294404/,c1/-11.6031403/,
     &     a2/0.0115572/,b2/ 0.6031679/,c2/-12.4327495/,
     &     a3/0.0436378/,b3/-0.0335511/,c3/ -9.5767430/,
     &     a4/0.0449290/,b4/-0.0826158/,c4/ -9.2019475/,
     &     a5/0.0268142/,b5/ 0.3177838/,c5/-11.3155513/
c array for collisional integrals
      dimension xr(5),xcinf(5)
      dimension dlsr1(5,5),dlsr2(5,5),dlsrv(5,5)
      dimension alsr(5,5),ddsr(5,5)
      dimension xxdl1(5),xxdl2(5),xxdlv(5),axdl2(5)
      dimension dsb(5)
c coefficients for collisional integrals by Gupta curve fitting
      dimension aog(5,5,2),bog(5,5,2),cog(5,5,2),dog(5,5,2)
c Omega_sr(1,1) and Omega_sr(2,2)
      data aog/50*0.d0/
      data bog/-0.0034d0, 0.0048d0,-0.0179d0,-0.0226d0,-0.0139d0, ! b(1-5,1,1) ! for Omega_sr(1,1)
     &          0.0048d0,-0.0033d0,-0.0185d0,-0.0179d0,-0.0194d0, ! b(1-5,2,1)
     &         -0.0179d0,-0.0185d0,-0.0364d0,-0.0438d0,-0.0291d0, ! b(1-5,3,1)
     &         -0.0226d0,-0.0179d0,-0.0438d0,-0.0410d0,-0.0465d0, ! b(1-5,4,1)
     &         -0.0139d0,-0.0194d0,-0.0291d0,-0.0465d0,-0.0112d0, ! b(1-5,5,1)
     &         -0.0207d0, 0.0065d0,-0.0203d0,-0.0247d0,-0.0169d0, ! b(1-5,1,2) ! for Omega_sr(2,2)
     &          0.0065d0,-0.0118d0,-0.0196d0,-0.0203d0,-0.0190d0, ! b(1-5,2,2)
     &         -0.0203d0,-0.0196d0,-0.0453d0,-0.0522d0,-0.0385d0, ! b(1-5,3,2)
     &         -0.0247d0,-0.0203d0,-0.0522d0,-0.0485d0,-0.0558d0, ! b(1-5,4,2)
     &         -0.0169d0,-0.0190d0,-0.0385d0,-0.0558d0,-0.0203d0/ ! b(1-5,5,2)
      data cog/-0.0572d0,-0.4195d0, 0.0152d0, 0.1300d0,-0.0825d0, ! c(1-5,1,1) ! for Omega_sr(1,1)
     &         -0.4195d0,-0.0572d0, 0.0118d0, 0.0152d0, 0.0119d0, ! c(1-5,2,1)
     &          0.0152d0, 0.0118d0, 0.3825d0, 0.5352d0, 0.2324d0, ! c(1-5,3,1)
     &          0.1300d0, 0.0152d0, 0.5352d0, 0.4977d0, 0.5729d0, ! c(1-5,4,1)
     &         -0.0825d0, 0.0119d0, 0.2324d0, 0.5729d0,-0.1182d0, ! c(1-5,5,1)
     &          0.0780d0,-0.4467d0, 0.0730d0, 0.1783d0, 0.0143d0, ! c(1-5,1,2) ! for Omega_sr(2,2)
     &         -0.4467d0,-0.0960d0, 0.0478d0, 0.0730d0, 0.0239d0, ! c(1-5,2,2)
     &          0.0730d0, 0.0478d0, 0.5624d0, 0.7045d0, 0.4226d0, ! c(1-5,3,2)
     &          0.1783d0, 0.0730d0, 0.7045d0, 0.6475d0, 0.7590d0, ! c(1-5,4,2)
     &          0.0143d0, 0.0239d0, 0.4226d0, 0.7590d0, 0.0683d0/ ! c(1-5,5,2)
      data dog/ 4.9901d0, 5.7774d0, 3.9996d0, 3.3363d0, 4.5785d0, ! d(1-5,1,1) ! for Omega_sr(1,1)
     &          5.7774d0, 5.0452d0, 4.0590d0, 3.9996d0, 4.1055d0, ! d(1-5,2,1)
     &          3.9996d0, 4.0590d0, 2.4718d0, 1.7252d0, 3.2082d0, ! d(1-5,3,1)
     &          3.3363d0, 3.9996d0, 1.7252d0, 1.8302d0, 1.6185d0, ! d(1-5,4,1)
     &          4.5785d0, 4.1055d0, 3.2082d0, 1.6185d0, 4.8464d0, ! d(1-5,5,1)
     &          3.5658d0, 6.0426d0, 3.8818d0, 3.2517d0, 4.4195d0, ! d(1-5,1,2) ! for Omega_sr(2,2)
     &          6.0426d0, 4.3252d0, 4.0321d0, 3.8818d0, 4.1782d0, ! d(1-5,2,2)
     &          3.8818d0, 4.0321d0, 1.7669d0, 1.0738d0, 2.4507d0, ! d(1-5,3,2)
     &          3.2517d0, 3.8818d0, 1.0738d0, 1.2607d0, 0.8955d0, ! d(1-5,4,2)
     &          4.4195d0, 4.1782d0, 2.4507d0, 0.8955d0, 4.0900d0/ ! d(1-5,5,2)
c
c* main
      sc =0.5
      rs1        =8.314*amas1(1)
      cvts1      =1.5*rs1
      cvrs1      =0.0
      cvvs1      =0.0
      rs2        =8.314*amas1(2)
      cvts2      =1.5*rs2
      cvrs2      =0.0
      cvvs2      =0.0
      rs3        =8.314*amas1(3)
      cvts3      =1.5*rs3
      cvrs3      =rs3
      cvvs3      =rs3
      rs4        =8.314*amas1(4)
      cvts4      =1.5*rs4
      cvrs4      =rs4
      cvvs4      =rs4
      rs5        =8.314*amas1(5)
      cvts5      =1.5*rs5
      cvrs5      =rs5
      cvvs5      =rs5
c Viscosity and conductivities by Wilke/Blottner/Eucken
c Diffusion coefficients by constant Schmidt number
      if(nvis.eq.0) then

         do j=1,jmax-1
            do i=1,imax-1
c* viscosity
               am         =1.0/(xc(i,j,1)*amas1(1)+xc(i,j,2)*amas1(2)
     &                    +xc(i,j,3)*amas1(3)+xc(i,j,4)*amas1(4)
     &                    +xc(i,j,5)*amas1(5))
               alnt       =log(tt(i,j))
c s1(O)
               x1         =xc(i,j,1)*am*amas1(1)
               vmus1      =0.1*exp((a1*alnt+b1)*alnt+c1)
               vkas1      =vmus1*(2.5*cvts1+cvrs1)
               vkvs1      =vmus1*cvvs1
c s2(N)
               x2         =xc(i,j,2)*am*amas1(2)
               vmus2      =0.1*exp((a2*alnt+b2)*alnt+c2)
               vkas2      =vmus2*(2.5*cvts2+cvrs2)
               vkvs2      =vmus2*cvvs2
c s3(NO)
               x3         =xc(i,j,3)*am*amas1(3)
               vmus3      =0.1*exp((a3*alnt+b3)*alnt+c3)
               vkas3      =vmus3*(2.5*cvts3+cvrs3)
               vkvs3      =vmus3*cvvs3
c s4(O2)
               x4         =xc(i,j,4)*am*amas1(4)
               vmus4      =0.1*exp((a4*alnt+b4)*alnt+c4)
               vkas4      =vmus4*(2.5*cvts4+cvrs4)
               vkvs4      =vmus4*cvvs4
c s5(N2)
               x5         =xc(i,j,5)*am*amas1(5)
               vmus5      =0.1*exp((a5*alnt+b5)*alnt+c5)
               vkas5      =vmus5*(2.5*cvts5+cvrs5)
               vkvs5      =vmus5*cvvs5
c mixing (1-*)
               phs11t  =1.0+sqrt(vmus1/vmus1)*(amass(1)*amas1(1))**0.25
               phs11   =phs11t*phs11t/sqrt(8.0*(1.0+amass(1)*amas1(1)))
               phs12t  =1.0+sqrt(vmus1/vmus2)*(amass(2)*amas1(1))**0.25
               phs12   =phs12t*phs12t/sqrt(8.0*(1.0+amass(1)*amas1(2)))
               phs13t  =1.0+sqrt(vmus1/vmus3)*(amass(3)*amas1(1))**0.25
               phs13   =phs13t*phs13t/sqrt(8.0*(1.0+amass(1)*amas1(3)))
               phs14t  =1.0+sqrt(vmus1/vmus4)*(amass(4)*amas1(1))**0.25
               phs14   =phs14t*phs14t/sqrt(8.0*(1.0+amass(1)*amas1(4)))
               phs15t  =1.0+sqrt(vmus1/vmus5)*(amass(5)*amas1(1))**0.25
               phs15   =phs15t*phs15t/sqrt(8.0*(1.0+amass(1)*amas1(5)))
               phis1   =1.0/(x1*phs11+x2*phs12+x3*phs13+x4*phs14
     &                                                 +x5*phs15)
c mixing (2-*)
               phs21t  =1.0+sqrt(vmus2/vmus1)*(amass(1)*amas1(2))**0.25
               phs21   =phs21t*phs21t/sqrt(8.0*(1.0+amass(2)*amas1(1)))
               phs22t  =1.0+sqrt(vmus2/vmus2)*(amass(2)*amas1(2))**0.25
               phs22   =phs22t*phs22t/sqrt(8.0*(1.0+amass(2)*amas1(2)))
               phs23t  =1.0+sqrt(vmus2/vmus3)*(amass(3)*amas1(2))**0.25
               phs23   =phs23t*phs23t/sqrt(8.0*(1.0+amass(2)*amas1(3)))
               phs24t  =1.0+sqrt(vmus2/vmus4)*(amass(4)*amas1(2))**0.25
               phs24   =phs24t*phs24t/sqrt(8.0*(1.0+amass(2)*amas1(4)))
               phs25t  =1.0+sqrt(vmus2/vmus5)*(amass(5)*amas1(2))**0.25
               phs25   =phs25t*phs25t/sqrt(8.0*(1.0+amass(2)*amas1(5)))
               phis2   =1.0/(x1*phs21+x2*phs22+x3*phs23+x4*phs24
     &                                                 +x5*phs25)
c mixing (3-*)
               phs31t  =1.0+sqrt(vmus3/vmus1)*(amass(1)*amas1(3))**0.25
               phs31   =phs31t*phs31t/sqrt(8.0*(1.0+amass(3)*amas1(1)))
               phs32t  =1.0+sqrt(vmus3/vmus2)*(amass(2)*amas1(3))**0.25
               phs32   =phs32t*phs32t/sqrt(8.0*(1.0+amass(3)*amas1(2)))
               phs33t  =1.0+sqrt(vmus3/vmus3)*(amass(3)*amas1(3))**0.25
               phs33   =phs33t*phs33t/sqrt(8.0*(1.0+amass(3)*amas1(3)))
               phs34t  =1.0+sqrt(vmus3/vmus4)*(amass(4)*amas1(3))**0.25
               phs34   =phs34t*phs34t/sqrt(8.0*(1.0+amass(3)*amas1(4)))
               phs35t  =1.0+sqrt(vmus3/vmus5)*(amass(5)*amas1(3))**0.25
               phs35   =phs35t*phs35t/sqrt(8.0*(1.0+amass(3)*amas1(5)))
               phis3   =1.0/(x1*phs31+x2*phs32+x3*phs33+x4*phs34
     &                                                 +x5*phs35)
c mixing (4-*)
               phs41t  =1.0+sqrt(vmus4/vmus1)*(amass(1)*amas1(4))**0.25
               phs41   =phs41t*phs41t/sqrt(8.0*(1.0+amass(4)*amas1(1)))
               phs42t  =1.0+sqrt(vmus4/vmus2)*(amass(2)*amas1(4))**0.25
               phs42   =phs42t*phs42t/sqrt(8.0*(1.0+amass(4)*amas1(2)))
               phs43t  =1.0+sqrt(vmus4/vmus3)*(amass(3)*amas1(4))**0.25
               phs43   =phs43t*phs43t/sqrt(8.0*(1.0+amass(4)*amas1(3)))
               phs44t  =1.0+sqrt(vmus4/vmus4)*(amass(4)*amas1(4))**0.25
               phs44   =phs44t*phs44t/sqrt(8.0*(1.0+amass(4)*amas1(4)))
               phs45t  =1.0+sqrt(vmus4/vmus5)*(amass(5)*amas1(4))**0.25
               phs45   =phs45t*phs45t/sqrt(8.0*(1.0+amass(4)*amas1(5)))
               phis4   =1.0/(x1*phs41+x2*phs42+x3*phs43+x4*phs44
     &                                                 +x5*phs45)
c mixing (5-*)
               phs51t  =1.0+sqrt(vmus5/vmus1)*(amass(1)*amas1(5))**0.25
               phs51   =phs51t*phs51t/sqrt(8.0*(1.0+amass(5)*amas1(1)))
               phs52t  =1.0+sqrt(vmus5/vmus2)*(amass(2)*amas1(5))**0.25
               phs52   =phs52t*phs52t/sqrt(8.0*(1.0+amass(5)*amas1(2)))
               phs53t  =1.0+sqrt(vmus5/vmus3)*(amass(3)*amas1(5))**0.25
               phs53   =phs53t*phs53t/sqrt(8.0*(1.0+amass(5)*amas1(3)))
               phs54t  =1.0+sqrt(vmus5/vmus4)*(amass(4)*amas1(5))**0.25
               phs54   =phs54t*phs54t/sqrt(8.0*(1.0+amass(5)*amas1(4)))
               phs55t  =1.0+sqrt(vmus5/vmus5)*(amass(5)*amas1(5))**0.25
               phs55   =phs55t*phs55t/sqrt(8.0*(1.0+amass(5)*amas1(5)))
               phis5   =1.0/(x1*phs51+x2*phs52+x3*phs53+x4*phs54
     &                                                 +x5*phs55)
c combined values
               vmue(i,j)  =x1*vmus1*phis1+x2*vmus2*phis2+x3*vmus3*phis3
     &                    +x4*vmus4*phis4+x5*vmus5*phis5
               vkap(i,j)  =x1*vkas1*phis1+x2*vkas2*phis2+x3*vkas3*phis3
     &                    +x4*vkas4*phis4+x5*vkas5*phis5
               vkav(i,j)  =x1*vkvs1*phis1+x2*vkvs2*phis2+x3*vkvs3*phis3
     &                    +x4*vkvs4*phis4+x5*vkvs5*phis5
c* diffusion coefficient
               d          =vmue(i,j)/(qc(i,j,1)*sc)
               ds(i,j,1)  =d
               ds(i,j,2)  =d
               ds(i,j,3)  =d
               ds(i,j,4)  =d
               ds(i,j,5)  =d
            end do
         end do
      end if
      if(nvis.eq.1) then
c Viscosity and conductivities by collisonal integrals
c using Gupta curve fitting
c Diffusion coefficients by binary diffusion coefficient
c 1:O, 2:N, 3:NO, 4:O2, 5:N2
         pi    = dacos(-1.0d0)
         avoga1= 1.d0/6.02214d23
         boltz = 1.3806503d-23  ! [m2kg/s2K]
         dct1  = 8.d0/3.d0      ! constant for delta_sr(1)
         dct2  = 16.d0/5.d0     ! constant for delta_sr(2)
         gasc  = 8.314d0
         do j=1,jmax-1
            do i=1,imax-1
               te  =tt(i,j)
               ttv =tv(i,j)
               pp  =pc(i,j)
               alte=dlog(te)
               alt2=alte*alte
               prt1=1.d0/pi/te/gasc
               altv=dlog(ttv)
               altv2=altv*altv
               prtv=1.d0/pi/ttv/gasc
               do n=1,5
                  xr(n)=xc(i,j,n)*amas1(n)
               end do
               xt=0.d0
               do n=1,5
                  xt=xt+xr(n)
               end do
c pi*Omega_sr(1,1)(=pomega1) and Omega_sr(2,2)(=pomega2)
c unit A^2(=10^{-20}[m])
                     do m=1,5   ! s
                        do l=1,5 ! r
                           pomega1   =dexp(dog(l,m,1))*te**(aog(l,m,1)
     &                               *alt2+bog(l,m,1)*alte+cog(l,m,1))
                           pomega2   =dexp(dog(l,m,2))*te**(aog(l,m,2)
     &                               *alt2+bog(l,m,2)*alte+cog(l,m,2))
                           pomegav   =dexp(dog(l,m,1))*ttv**(aog(l,m,1)
     &                              *altv2+bog(l,m,1)*altv+cog(l,m,1))
c delta_sr^1(=dlsr1) and delta_sr^2(=dlsr2)
                           dlsr1(l,m)=dct1*dsqrt(2.d0*amass(m)*amass(l)
     &                         /(amass(m)+amass(l))*prt1)*pomega1*1.d-20
                           dlsr2(l,m)=dct2*dsqrt(2.d0*amass(m)*amass(l)
     &                         /(amass(m)+amass(l))*prt1)*pomega2*1.d-20
                           dlsrv(l,m)=dct1*dsqrt(2.d0*amass(m)*amass(l)
     &                         /(amass(m)+amass(l))*prtv)*pomegav*1.d-20
c alpha_sr^1(=alsr)
                           ramass    =amass(m)/amass(l)
                           alsr(l,m) = 1.d0 + ((1.d0-ramass)
     &                               *(0.45d0-2.54d0*ramass)
     &                               /((1.d0+ramass)**2))
                           ddsr(l,m)  = (boltz*te)/(pp*dlsr1(l,m))
                        end do
                     end do
c
                     xxdl1(:)=0.d0
                     xxdl2(:)=0.d0
                     xxdlv(:)=0.d0
                     axdl2(:)=0.d0
                     do m = 1,5 ! s
                        do l = 1,5 ! r
                           xxdl1(m) = xxdl1(m)+xr(l)*dlsr1(l,m)
                           xxdl2(m) = xxdl2(m)+xr(l)*dlsr2(l,m)
                           xxdlv(m) = xxdlv(m)+xr(l)*dlsrv(l,m)
                           axdl2(m) = axdl2(m)+xr(l)*dlsr2(l,m)
     &                               *alsr(l,m)
                        end do
                     end do
                     vics=0.d0
                     thct=0.d0
                     thcr=0.d0
                     thcv=0.d0
                     do m=1,5   !O,N,NO,O2,N2
                        vics=vics+xc(i,j,m)/xxdl2(m)
                        thct=thct+xr(m)/axdl2(m)
                     end do
                     do m=3,5   ! NO,O2,N2
                        thcr=thcr+xr(m)/xxdl1(m)
                        thcv=thcv+xr(m)/xxdlv(m)
                     end do
                     vmue(i,j)=vics*avoga1
                     vkap(i,j)=boltz*(3.75d0*thct+thcr)
                     vkav(i,j)=boltz*thcv
c binary diffusion coefficient
                     dsb(:)=0.d0
                     do m=1,5
                        do l=1,5
                           if(m.ne.l) then
                              dsb(m) = dsb(m) + xr(l)/ddsr(l,m)
                           end if
                        end do
                     end do
                     do m=1,5
                        ds(i,j,m) = ((xt**2)*amass(m)
     &                           *(1.d0-xc(i,j,m)))/dsb(m)
                     end do
                  end do
               end do
            end if
c* transport coefficient for vibration-electronic energy
            do j=1,jmax-1
               do i=1,imax-1
                  tv1        =1.0/tv(i,j)
                  xa         =exp(226.0*tv1)
                  xb         =exp(-3785.5*tv1)
            xa2        =xa*xa
            xa4        =xa2*xa2
            xa8        =xa4*xa4
            xa10       =xa8*xa2
            xa12       =xa10*xa2
            xa15       =xa12*xa2*xa
            xb2        =xb*xb
            xb3        =xb2*xb
            xb4        =xb2*xb2
            xb5        =xb4*xb
            xb6        =xb5*xb
            xb7        =xb6*xb
            ev1        =2260.0/(xa10-1.0)
            ev2        =3390.0/(xa15-1.0)
            ev3        =2712.0/(xa12-1.0)
            ev4        =(22713.0*xb3+18927.5*xb5)/(3.0+2.0*xb3+xb5)
            ev5        =113565.0*xb6/(9.0+5.0*xb6)
            ev6        =264985.0*xb7/(4.0+10.0*xb7)
            ev1g       =8.314*ev1
            ev2g       =8.314*ev2
            ev3g       =8.314*ev3
            ev4g       =8.314*ev4
            ev5g       =8.314*ev5
            ev6g       =8.314*ev6
            evs(i,j,1) =      ev5g *amas1(1)
            evs(i,j,2) =      ev6g *amas1(2)
            evs(i,j,3) =      ev3g *amas1(3)
            evs(i,j,4) =(ev1g+ev4g)*amas1(4)
            evs(i,j,5) =      ev2g *amas1(5)
            hs(i,j,1)  =(cvts1+cvrs1+rs1)*tt(i,j)
     &                 +(h0i(1)+ev5g)*amas1(1)
            hs(i,j,2)  =(cvts2+cvrs2+rs2)*tt(i,j)
     &                 +(h0i(2)+ev6g)*amas1(2)
            hs(i,j,3)  =(cvts3+cvrs3+rs3)*tt(i,j)
     &                 +(h0i(3)+ev3g)*amas1(3)
            hs(i,j,4)  =(cvts4+cvrs4+rs4)*tt(i,j)
     &                 +(  ev1g+ev4g)*amas1(4)
            hs(i,j,5)  =(cvts5+cvrs5+rs5)*tt(i,j)
     &                 +        ev2g *amas1(5)
         end do
      end do
c* termination
      return
      end
c
      subroutine baldwin_lomax
c     ****************************************************************
c     *     RANS 0-eq turbulance model -Baldwin-Lomax model-
c     *       ref. Baldwin, B. S., and Lomax, H.,
c     *            "Thin Layer Approximation and Algebraic Model
c     *             for Separated Turbulent Flows", AIAA Paper 78-257
c     ****************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /coor4/ sl(0:nx,0:ny) !
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      dimension dudx(0:nx,0:ny),dvdx(0:nx,0:ny)
     &         ,dudy(0:nx,0:ny),dvdy(0:nx,0:ny)
      dimension vmul(0:nx,0:ny),vmut(0:nx,0:ny),yti(0:nx,0:ny)
     &         ,yto(0:nx,0:ny),  yt(0:nx,0:ny), ff(0:nx,0:ny)

c*preset
      do j=1,jmax-1
         do i=1,imax-1
            vmul(i,j)=vmue(i,j)
         end do
      end do
c*constant
      aps  = 26.d0      ! A
      ckar =  0.4d0     ! Karman constant
      cle  =  0.3d0
      cwk  =  0.25d0
      alk  =  0.0168d0
      ccp  =  1.6d0

      prnl =  0.72d0
      prnt =  0.90d0

      sc   =  0.5d0

c*set du/dy,dv/dx and vorticity
      do i=1, imax-1
         do j=1, jmax-1
            duip      =uc(i+1,j  )-uc(i  ,j  )
            duim      =uc(i  ,j  )-uc(i-1,j  )
            dujp      =uc(i  ,j+1)-uc(i  ,j  )
            dujm      =uc(i  ,j  )-uc(i  ,j-1)
            dvip      =vc(i+1,j  )-vc(i  ,j  )
            dvim      =vc(i  ,j  )-vc(i-1,j  )
            dvjp      =vc(i  ,j+1)-vc(i  ,j  )
            dvjm      =vc(i  ,j  )-vc(i  ,j-1)
!
!            duci      =0.5*(duip+duim)
!            dvci      =0.5*(dvip+dvim)
!            ducj      =0.5*(dujp+dujm)
!            dvcj      =0.5*(dvjp+dvjm)
!
            duci    =(sign(0.5d0,duip)+sign(0.5d0,duim))
     &                *dmin1(dabs(duip),     dabs(duim))
            dvci    =(sign(0.5d0,dvip)+sign(0.5d0,dvim))
     &                *dmin1(dabs(dvip),     dabs(dvim))
            ducj    =(sign(0.5d0,dujp)+sign(0.5d0,dujm))
     &                *dmin1(dabs(dujp),     dabs(dujm))
            dvcj    =(sign(0.5d0,dvjp)+sign(0.5d0,dvjm))
     &                *dmin1(dabs(dvjp),     dabs(dvjm))
            volv      =1.0/vol(i,j)
            sxx       =0.5d0*(si(i+1,j  ,1)+si(i,j,1))
            sxy       =0.5d0*(si(i+1,j  ,2)+si(i,j,2))
            sex       =0.5d0*(sj(i  ,j+1,1)+sj(i,j,1))
            sey       =0.5d0*(sj(i  ,j+1,2)+sj(i,j,2))
            sxxv      =sxx*volv
            sxyv      =sxy*volv
            sexv      =sex*volv
            seyv      =sey*volv
            dudx(i,j) =sxxv*duci+sexv*ducj
            dudy(i,j) =sxyv*duci+seyv*ducj
            dvdx(i,j) =sxxv*dvci+sexv*dvcj
            dvdy(i,j) =sxyv*dvci+seyv*dvcj
!            write(*,*) "i,j,duci,ducj",duci,ducj
!            write(*,*)"sxyv,seyv",seyv,sexv
!            write(*,*) "i,j,dvdx,dudy",dvdx(i,j),dudy(i,j)
         end do
      end do
c
c* yt (inner)
c* calc rho_w,vmulw, u_w, tau_w, u_tau
      do i=1,imax-1
         rhow  = 0.5d0*(  qc(i,0,1)+  qc(i,1,1))
         vmuw  = 0.5d0*(vmul(i,0  )+vmul(i,1  ))
         aex   = si(i,1,1)/si(i,1,3)
         aey   = si(i,1,2)/si(i,1,3)
         uw    = uc(i,1)*aey-vc(i,1)*aex    !uw is velocity along the wall
!
!         tauw  = vmuw*dabs(uw/sl(i,1))
         tauw  = vmul(i,1)*dabs(uw/sl(i,1))  !!!
!
         utau  = dsqrt(tauw/rhow)
         do j=1,jmax-1
            rho1     = qc(i,j,1)
!
!            yplus    = rhow*utau*sl(i,j)/vmuw           ! original in ref
            yplus    = rhow*utau*sl(i,j)/vmul(i,1)       !!!
!
!            dd       = 1.d0-dexp(-yplus/aps)            !van driest
            dd       =1.0-dexp(-yplus*dsqrt(rho1/rhow)
     &                          *vmuw/vmul(i,j)/aps)     !van driest for hypersonic

            dl       = ckar*sl(i,j)*dd                   !mixing length
            voab     = dabs(dvdx(i,j)-dudy(i,j))
            yti(i,j) = rho1*dl*dl*voab
!
!            ff(i,j)  = sl(i,j)*voab*dd                  ! original for yt(outer)
            ff(i,j)  = sl(i,j)*voab*dd+1.d-20            !!!
!
         end do
      end do
c* yt (outer)
c* search fmax,ymax,umax,umin
      do i=1,imax-1
         fmax = -1.0d10
         ymax =  0.d0
         umax = -1.0d10
         umin =  1.0d10
         do j=1,jmax-1
            if(ff(i,j).gt.fmax) then
               ymax = sl(i,j)
               fmax = ff(i,j)
            end if
            umax = dmax1(umax, dsqrt(uc(i,j)**2+vc(i,j)**2))
            umin = dmin1(umin, dsqrt(uc(i,j)**2+vc(i,j)**2))
         end do

         udif = umax-umin
         if(udif.eq.0.d0) udif=umax

         do j=1,jmax-1
!
!            fkleb    = 1.d0/(1.d0+5.5d0*(cle*sl(i,j)/ymax)**6) ! original
            fkleb    = 1.d0/(1.d0+5.5d0*(cle*sl(i,j)/(ymax+1.d-30))**6)  !!
!
            fwake    = ymax*dmin1(fmax,cwk*udif**2/fmax)
            rho2     = qc(i,j,1)
            yto(i,j) = alk*ccp*rho2*fwake*fkleb
         end do
      end do

c* choose yt(inner) or yt(outer)
      do i=1,imax-1
         isw = 0
         do j=1,jmax-1
            if (isw.eq.0) then
               if(yti(i,j).lt.yto(i,j)) then
                  yt(i,j) = yti(i,j)
               else
                  yt(i,j) = yto(i,j)
                  isw = 1
               end if
            else
            yt(i,j) = yto(i,j)
            end if
         end do
      end do
c* update vmue,vkap,vkav,ds
      do i=1,imax-1
         do j=1, jmax-1
c*viscosity
            vmut(i,j) = yt(i,j)
            vmue(i,j) = vmul(i,j)+vmut(i,j)
c*conductivity
            vkapl = vkap(i,j)
            vkavl = vkav(i,j)
            vkapt = vkapl*prnl*vmut(i,j)/prnt/vmul(i,j)
            vkavt = vkavl*prnl*vmut(i,j)/prnt/vmul(i,j)
            vkap(i,j) = vkapl + vkapt
            vkav(i,j) = vkavl + vkavt
c*diffusion coefficient
            rho     = qc(i,j,1)
            d       = vmue(i,j)/(rho*sc)
            ds(i,j,1) = d
            ds(i,j,2) = d
            ds(i,j,3) = d
            ds(i,j,4) = d
            ds(i,j,5) = d
         end do
      end do
c
      return
      end
c
      subroutine aero_coeff
c     **************************************************************
c     *     calculate aerodynamics coefficient                     *
c     *     accounting for only pressure drag                      *
c     *     you can use only when axisymmetrics calculation,       *
c     *     that is, AOA=0 deg                                     *
c     **************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /coor4/ sl(0:nx,0:ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /parm3/ mdd
      common /stat1/ loop,locl
      common /ref/ yRef
c* set constant
      pi = dacos(-1.0d0)
      pd = 0.d0
      fd = 0.d0
      dd = 0.d0
      cd = 0.d0
c  reference area
      ! 代表長さを出す
      yRef = 0.0d0
      do i=1,imax
         if (y(i,1).gt.yRef) then
            yRef = y(i,1)
         endif
      end do
      ra = pi*yRef*yRef
      qqinf = 0.5d0*rinf*uinf*uinf
c* calc. drag
      do i=1,imax-1
c angle between normal direction to wall and x-direction (freestream direction)
c calculate area of a truncated cone
         csth = -sj(i,1,1)/sj(i,1,3)
         sth  =  sj(i,1,2)/sj(i,1,3)
c differential area of rotaling body, ds
         dxm = x(i+1,1)-x(i,1)
         dym = y(i+1,1)-y(i,1)
         dyp = y(i+1,1)+y(i,1)
         ds = pi*dyp*dxm*dsqrt(1.0d0+(dym/dxm)*(dym/dxm))
c pressure drag
         pd = pd+pc(i,1)*csth*ds
c if you want to be friction drag= 0, please comment out following 8 lines
c friction drag
         rhow  = 0.5d0*(  qc(i,0,1)+  qc(i,1,1))
         vmuw  = 0.5d0*(vmue(i,0  )+vmue(i,1  ))
         aex   = si(i,1,1)/si(i,1,3)
         aey   = si(i,1,2)/si(i,1,3)
         uw    = uc(i,1)*aey-vc(i,1)*aex !uw is velocity along the wall
         tauw  = vmuw*dabs(uw/sl(i,1))
c         tauw  = vmul(i,1)*dabs(uw/sl(i,1))
         fd = fd+tauw*sth*ds
      end do
c* calc. drag coefficient
      pd =pd/(ra*qqinf)
      fd =fd/(ra*qqinf)
      cd = pd +fd
c
      write(6,*) ''
      write(6,*) 'yRef = ',yRef
      write(6,*) 'ra*qqinf = ',ra*qqinf
      write(6,600) ra,cd,pd,fd
      write(77,601) loop,cd,pd,fd
c
 600  format(' aero:',' area',e13.6,' Cd',f10.6,' pd',f10.6,' fd',f10.6)
 601  format(i10,' ',4(1x,e13.6))
c* termination
      return
      end
c
      subroutine cpcf_coeff
c     *******************************************************************
c     *     calculate distribution of pressure and friction coefficient  *
c     *******************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /coor4/ sl(0:nx,0:ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /parm3/ mdd
      common /stat1/ loop,locl
      dimension sw(0:nx,2),s(0:nx)
c open file
      open(79,file='./output/cpcf_coeff.dat',form='formatted')
      write(79,778) '#','x','s','cp','cf'
 778  format(a,a12,3(a13))
c* set constant
      rwall   = dsqrt(x(1,1)**2+y(1,1)**2)
      sw(0,1) = -rwall
      sw(0,2) = 0.d0
      s(0) = 0.d0
      do i=1,imax-1
         sw(i,1) = (x(i,1)+x(i+1,1))*0.5d0
         sw(i,2) = (y(i,1)+y(i+1,1))*0.5d0
         dl      = dsqrt((sw(i,1)-sw(i-1,1))**2+(sw(i,2)-sw(i-1,2))**2)
         s(i)    = s(i-1)+dl
      end do
      pi  = dacos(-1.0d0)
      pw  = 0.d0
      fw  = 0.d0
      ps  = pinf
      cp  = 0.d0   ! pressure coefficient
      cf  = 0.d0   ! skin friction coefficient
      pqq = 0.5d0*rinf*uinf*uinf
c calc pressure and fiction drag
      do i=1,imax-1
c pressure at wall
         pw    = 0.5d0*(pc(i,1)+pc(i,0))
c skin friction at wall
         rhow  = 0.5d0*(  qc(i,0,1)+  qc(i,1,1))
         vmuw  = 0.5d0*(vmue(i,0  )+vmue(i,1  ))
         aex   = si(i,1,1)/si(i,1,3)
         aey   = si(i,1,2)/si(i,1,3)
         uw    = uc(i,1)*aey-vc(i,1)*aex !uw is velocity along the wall
         tauw  = vmuw*dabs(uw/sl(i,1))
         fw    = tauw
c pressure and skin friction at wall
         cp = (pw-ps)/pqq
         cf = fw/pqq
c
         pstag = 0.5d0*(pc(1,1)+pc(1,0))
c
         write(79,602) (x(i,1)-x(1,1))*1.d3,s(i),cp,cf,
     &                 (x(i,1)-x(1,1))/(-x(1,1)),pw/pstag
       end do
       close(79)
 602   format(6(e13.6))
c
c termination
       return
       end
c
      subroutine calc_knudsen
!     ************************************************************
!     *     calculate local Knudsen number                       *
!     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /cntl5/ nprc,iprc
      common /cntl10/nvis,nrlx
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /ref/   yRef
      common /knud1/ aLocalKnud(0:nx,0:ny),frPth(0:nx,0:ny)
!
      dimension thermVelo(0:nx,0:ny)
      double precision totalWeight, meanMass, knudFlag

      ! constants
      pi = dacos(-1.0d0)
      boltz = 1.3806503d-23 ! boltzmann const. [m2kg/s2K]
      avoga = 6.02214d23 ! avogadro const. [1/mole]

      !initialize valuables
      totalWeight = 0.0d0 ! total mean mokecular weight [kg/mole]
      meanWeight = 0.0d0 ! mean molecular mass [kg]
      knudFlag = 0
      do m=1,5
         totalWeight = totalWeight + amass(m)
      enddo
      meanMass = 0.2d0*totalWeight/avoga
      do i=1,imax-1
         do j=1,jmax-1
            ! calc. Thermal velocity
            thermVelo(i,j) = dsqrt(8.0d0*boltz*tt(i,j)
     &                       /(pi*meanMass))
            ! calc. Mean free path
            frPth(i,j) = 3.0d0*vmue(i,j)/(qc(i,j,1)*thermVelo(i,j))
            ! calc. local Knudsen number
            aLocalKnud(i,j) = frPth(i,j)/(2.0d0*yRef)
         enddo
      enddo
      write(*,*)pi,boltz,avoga,amass(1),yRef,vmue(1,1),qc(1,1,1),tt(1,1)
      write(*,*)totalWeight,meanMass,thermVelo(1,1),frPth(1,1)
!
! termination
      return
      end
!
      subroutine knudsen_gll
!     ************************************************************
!     *     calculate gradient-length local Knudsen number       *
!     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /cntl5/ nprc,iprc
      common /cntl10/nvis,nrlx
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /ref/   yRef
      common /knud1/ aLocalKnud(0:nx,0:ny),frPth(0:nx,0:ny),
     &               thermVelo(0:nx,0:ny)
!
      double precision,dimension(3) :: scaleTt
      double precision,dimension(3) :: scaleTv
      double precision,dimension(3) :: scaleRoe
      double precision,dimension(3) :: scaleVelo
      double precision meanWeight,meanMass,knudFlag,minScale

      ! constants
      pi = dacos(-1.0d0)
      boltz = 1.3806503d-23 ! boltzmann const. [m2kg/s2K]
      avoga = 6.02214d23 ! avogadro const. [1/mole]
      gasc  = 8.314d0 ! gas constant

      !initialize valuables
      meanMass = 0.0d0 ! mean molecular mass [kg]
      knudFlag = 0
      open(123,file='./knud.dat',status='replace')
      do i=1,imax-1
         do j=1,jmax-1
            ! mean mass
            meanMass = 0.0d0 !initialize
            do m=1,5
               meanMass = meanMass + xc(i,j,m)*amass(m)
            enddo
            ! calc. Thermal velocity
            thermVelo(i,j) = dsqrt(8.0d0*gasc*tt(i,j)
     &                       /(pi*meanMass))
            ! calc. Mean free path
            frPth(i,j) = 3.0d0*vmue(i,j)/(qc(i,j,1)*thermVelo(i,j))
            ! gradients
            scaleTt(1)  = 2*sj(i,j,3)/(tt(i+1,j)-tt(i-1,j))*tt(i,j)
            scaleTt(2)  = 2*si(i,j,3)/(tt(i,j+1)-tt(i,j-1))*tt(i,j)
            scaleTv(1)  = 2*sj(i,j,3)/(tv(i+1,j)-tv(i-1,j))*tv(i,j)
            scaleTv(2)  = 2*si(i,j,3)/(tv(i,j+1)-tv(i,j-1))*tv(i,j)
            scaleRoe(1) = 2*sj(i,j,3)
     &                    /(qc(i+1,j,1)-qc(i-1,j,1))*qc(i,j,1)
            scaleRoe(2) = 2*si(i,j,3)
     &                    /(qc(i,j+1,1)-qc(i,j-1,1))*qc(i,j,1)
            ! calc. gradient of Velo
            scaleVelo(1)= 2*sj(i,j,3)
     &                    /(dsqrt(uc(i+1,j)**2+vc(i+1,j)**2)
     &                    -dsqrt(uc(i-1,j)**2+vc(i-1,j)**2))
     &                    *dsqrt(uc(i,j)**2+vc(i,j))
            scaleVelo(2)= 2*si(i,j,3)
     &                    /(dsqrt(uc(i,j+1)**2+vc(i,j+1)**2)
     &                    -dsqrt(uc(i,j-1)**2+vc(i,j-1)**2))
     &                    *dsqrt(uc(i,j)**2+vc(i,j))
            ! index 3 is RMS
            scaleTt(3)  = dsqrt(scaleTt(1)**2 + scaleTt(2)**2)
            scaleTv(3)  = dsqrt(scaleTv(1)**2 + scaleTv(2)**2)
            scaleRoe(3) = dsqrt(scaleRoe(1)**2 + scaleRoe(2)**2)
            scaleVelo(3)= dsqrt(scaleVelo(1)**2 + scaleVelo(2)**2)
            ! minimam scale
            minScale =
     &      min(scaleTt(3),scaleTv(3),scaleRoe(3),scaleVelo(3))
            ! calc. local Knudsen number
            aLocalKnud(i,j) = frPth(i,j)/minScale
            write(123,*) i, j, aLocalKnud(i,j)
         enddo
      enddo
      close(123)
!
! termination
      return
      end
!
      subroutine rbdy(mode)
c     ************************************************************
c     *     set boundary values                                  *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /cntl5/ nprc,iprc
      common /cntl10/nvis,nrlx
      common /stat1/ loop,locl
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /fscd3/ vmueinf,vkapinf,vkavinf,                       !
     &               ds1inf,ds2inf,ds3inf,ds4inf,ds5inf,            !
     &               evs1inf,evs2inf,evs3inf,evs4inf,evs5inf,       !
     &               hs1inf,hs2inf,hs3inf,hs4inf,hs5inf             !
      common /wbcd1/ twall,cefn,cefo
      data a1/0.0203144/,b1/ 0.4294404/,c1/-11.6031403/,
     &     a2/0.0115572/,b2/ 0.6031679/,c2/-12.4327495/,
     &     a3/0.0436378/,b3/-0.0335511/,c3/ -9.5767430/,
     &     a4/0.0449290/,b4/-0.0826158/,c4/ -9.2019475/,
     &     a5/0.0268142/,b5/ 0.3177838/,c5/-11.3155513/
c array for collisional integrals
      dimension xr(5),xcinf(5)
      dimension dlsr1(5,5),dlsr2(5,5),dlsrv(5,5)
      dimension alsr(5,5),ddsr(5,5)
      dimension xxdl1(5),xxdl2(5),xxdlv(5),axdl2(5)
      dimension dsb(5)
c coefficients for collisional integrals by Gupta curve fitting
      dimension aog(5,5,2),bog(5,5,2),cog(5,5,2),dog(5,5,2)
c Omega_sr(1,1) and Omega_sr(2,2)
      data aog/50*0.d0/
      data bog/-0.0034d0, 0.0048d0,-0.0179d0,-0.0226d0,-0.0139d0, ! b(1-5,1,1) ! for Omega_sr(1,1)
     &          0.0048d0,-0.0033d0,-0.0185d0,-0.0179d0,-0.0194d0, ! b(1-5,2,1)
     &         -0.0179d0,-0.0185d0,-0.0364d0,-0.0438d0,-0.0291d0, ! b(1-5,3,1)
     &         -0.0226d0,-0.0179d0,-0.0438d0,-0.0410d0,-0.0465d0, ! b(1-5,4,1)
     &         -0.0139d0,-0.0194d0,-0.0291d0,-0.0465d0,-0.0112d0, ! b(1-5,5,1)
     &         -0.0207d0, 0.0065d0,-0.0203d0,-0.0247d0,-0.0169d0, ! b(1-5,1,2) ! for Omega_sr(2,2)
     &          0.0065d0,-0.0118d0,-0.0196d0,-0.0203d0,-0.0190d0, ! b(1-5,2,2)
     &         -0.0203d0,-0.0196d0,-0.0453d0,-0.0522d0,-0.0385d0, ! b(1-5,3,2)
     &         -0.0247d0,-0.0203d0,-0.0522d0,-0.0485d0,-0.0558d0, ! b(1-5,4,2)
     &         -0.0169d0,-0.0190d0,-0.0385d0,-0.0558d0,-0.0203d0/ ! b(1-5,5,2)
      data cog/-0.0572d0,-0.4195d0, 0.0152d0, 0.1300d0,-0.0825d0, ! c(1-5,1,1) ! for Omega_sr(1,1)
     &         -0.4195d0,-0.0572d0, 0.0118d0, 0.0152d0, 0.0119d0, ! c(1-5,2,1)
     &          0.0152d0, 0.0118d0, 0.3825d0, 0.5352d0, 0.2324d0, ! c(1-5,3,1)
     &          0.1300d0, 0.0152d0, 0.5352d0, 0.4977d0, 0.5729d0, ! c(1-5,4,1)
     &         -0.0825d0, 0.0119d0, 0.2324d0, 0.5729d0,-0.1182d0, ! c(1-5,5,1)
     &          0.0780d0,-0.4467d0, 0.0730d0, 0.1783d0, 0.0143d0, ! c(1-5,1,2) ! for Omega_sr(2,2)
     &         -0.4467d0,-0.0960d0, 0.0478d0, 0.0730d0, 0.0239d0, ! c(1-5,2,2)
     &          0.0730d0, 0.0478d0, 0.5624d0, 0.7045d0, 0.4226d0, ! c(1-5,3,2)
     &          0.1783d0, 0.0730d0, 0.7045d0, 0.6475d0, 0.7590d0, ! c(1-5,4,2)
     &          0.0143d0, 0.0239d0, 0.4226d0, 0.7590d0, 0.0683d0/ ! c(1-5,5,2)
      data dog/ 4.9901d0, 5.7774d0, 3.9996d0, 3.3363d0, 4.5785d0, ! d(1-5,1,1) ! for Omega_sr(1,1)
     &          5.7774d0, 5.0452d0, 4.0590d0, 3.9996d0, 4.1055d0, ! d(1-5,2,1)
     &          3.9996d0, 4.0590d0, 2.4718d0, 1.7252d0, 3.2082d0, ! d(1-5,3,1)
     &          3.3363d0, 3.9996d0, 1.7252d0, 1.8302d0, 1.6185d0, ! d(1-5,4,1)
     &          4.5785d0, 4.1055d0, 3.2082d0, 1.6185d0, 4.8464d0, ! d(1-5,5,1)
     &          3.5658d0, 6.0426d0, 3.8818d0, 3.2517d0, 4.4195d0, ! d(1-5,1,2) ! for Omega_sr(2,2)
     &          6.0426d0, 4.3252d0, 4.0321d0, 3.8818d0, 4.1782d0, ! d(1-5,2,2)
     &          3.8818d0, 4.0321d0, 1.7669d0, 1.0738d0, 2.4507d0, ! d(1-5,3,2)
     &          3.2517d0, 3.8818d0, 1.0738d0, 1.2607d0, 0.8955d0, ! d(1-5,4,2)
     &          4.4195d0, 4.1782d0, 2.4507d0, 0.8955d0, 4.0900d0/ ! d(1-5,5,2)
c* (i)
      if(mode.eq.1) then
c*axis boundary condition (velocity mirror)
         do j=1,jmax-1
            qs(0,j,1)  = qs(1,j,1)
            qs(0,j,2)  = qs(1,j,2)
            qs(0,j,3)  = qs(1,j,3)
            qv(0,j,1)  = qv(1,j,1)
            qc(0,j,2)  = qc(1,j,2)
            qc(0,j,3)  =-qc(1,j,3)
            qc(0,j,4)  = qc(1,j,4)
            pc(0,j)    = pc(1,j)
            cc(0,j)    = cc(1,j)
            tt(0,j)    = tt(1,j)
            tv(0,j)    = tv(1,j)
            qc(0,j,1)  = qc(1,j,1)
            uc(0,j)    = uc(1,j)
            vc(0,j)    =-vc(1,j)
            xc(0,j,1)  =xc(1,j,1)
            xc(0,j,2)  =xc(1,j,2)
            xc(0,j,3)  =xc(1,j,3)
            xc(0,j,4)  =xc(1,j,4)
            xc(0,j,5)  =xc(1,j,5)
            vmue(0,j)  =vmue(1,j)
            vkap(0,j)  =vkap(1,j)
            vkav(0,j)  =vkav(1,j)
            ds(0,j,1)  =ds(1,j,1)
            ds(0,j,2)  =ds(1,j,2)
            ds(0,j,3)  =ds(1,j,3)
            ds(0,j,4)  =ds(1,j,4)
            ds(0,j,5)  =ds(1,j,5)
            hs(0,j,1)  =hs(1,j,1)
            hs(0,j,2)  =hs(1,j,2)
            hs(0,j,3)  =hs(1,j,3)
            hs(0,j,4)  =hs(1,j,4)
            hs(0,j,5)  =hs(1,j,5)
            evs(0,j,1) =evs(1,j,1)
            evs(0,j,2) =evs(1,j,2)
            evs(0,j,3) =evs(1,j,3)
            evs(0,j,4) =evs(1,j,4)
            evs(0,j,5) =evs(1,j,5)
         end do
c*axis boundary condition (velocity mirror)
         do j=1,jmax-1
            qs(imax,j,1)  = qs(imax-1,j,1)
            qs(imax,j,2)  = qs(imax-1,j,2)
            qs(imax,j,3)  = qs(imax-1,j,3)
            qv(imax,j,1)  = qv(imax-1,j,1)
            qc(imax,j,2)  = qc(imax-1,j,2)
            qc(imax,j,3)  = -qc(imax-1,j,3)
            qc(imax,j,4)  = qc(imax-1,j,4)
            pc(imax,j)    = pc(imax-1,j)
            cc(imax,j)    = cc(imax-1,j)
            tt(imax,j)    = tt(imax-1,j)
            tv(imax,j)    = tv(imax-1,j)
            qc(imax,j,1)  = qc(imax-1,j,1)
            uc(imax,j)    = uc(imax-1,j)
            vc(imax,j)    = -vc(imax-1,j)
            xc(imax,j,1)  = xc(imax-1,j,1)
            xc(imax,j,2)  = xc(imax-1,j,2)
            xc(imax,j,3)  = xc(imax-1,j,3)
            xc(imax,j,4)  = xc(imax-1,j,4)
            xc(imax,j,5)  = xc(imax-1,j,5)
            vmue(imax,j)  = vmue(imax-1,j)
            vkap(imax,j)  = vkap(imax-1,j)
            vkav(imax,j)  = vkav(imax-1,j)
            ds(imax,j,1)  = ds(imax-1,j,1)
            ds(imax,j,2)  = ds(imax-1,j,2)
            ds(imax,j,3)  = ds(imax-1,j,3)
            ds(imax,j,4)  = ds(imax-1,j,4)
            ds(imax,j,5)  = ds(imax-1,j,5)
            hs(imax,j,1)  = hs(imax-1,j,1)
            hs(imax,j,2)  = hs(imax-1,j,2)
            hs(imax,j,3)  = hs(imax-1,j,3)
            hs(imax,j,4)  = hs(imax-1,j,4)
            hs(imax,j,5)  = hs(imax-1,j,5)
            evs(imax,j,1) = evs(imax-1,j,1)
            evs(imax,j,2) = evs(imax-1,j,2)
            evs(imax,j,3) = evs(imax-1,j,3)
            evs(imax,j,4) = evs(imax-1,j,4)
            evs(imax,j,5) = evs(imax-1,j,5)
         end do
      end if
c* (j)
      if(mode.eq.2) then
c Euler
         if(nprc.eq.0) then
            do i=1,imax-1
               sjv       =1.0/sj(i,1,3)
               sjx       =sjv*sj(i,1,1)
               sjy       =sjv*sj(i,1,2)
               rqn       =2.0*(qc(i,1,2)*sjx+qc(i,1,3)*sjy)
               qs(i,0,1) =qs(i,1,1)
               qs(i,0,2) =qs(i,1,2)
               qs(i,0,3) =qs(i,1,3)
               qv(i,0,1) =qv(i,1,1)
               qc(i,0,1) =qc(i,1,1)
               qc(i,0,2) =qc(i,1,2)-rqn*sjx
               qc(i,0,3) =qc(i,1,3)-rqn*sjy
               qc(i,0,4) =qc(i,1,4)
               pc(i,0)   =pc(i,1)
               cc(i,0)   =cc(i,1)
               tt(i,0)   =tt(i,1)
               tv(i,0)   =tv(i,1)
               rv        =1.0/qc(i,0,1)
               uc(i,0)   =rv*qc(i,0,2)
               vc(i,0)   =rv*qc(i,0,3)
               xc(i,0,1) =xc(i,1,1)
               xc(i,0,2) =xc(i,1,2)
               xc(i,0,3) =xc(i,1,3)
               xc(i,0,4) =xc(i,1,4)
               xc(i,0,5) =xc(i,1,5)
            end do
c NS(laminar)
         else
c adiabatic wall
            if(twall.le.0.0) then
               do i=1,imax-1
                  qs(i,0,1)  = qs(i,1,1)
                  qs(i,0,2)  = qs(i,1,2)
                  qs(i,0,3)  = qs(i,1,3)
                  qv(i,0,1)  = qv(i,1,1)
                  qc(i,0,1)  = qc(i,1,1)
                  qc(i,0,2)  =-qc(i,1,2)
                  qc(i,0,3)  =-qc(i,1,3)
                  qc(i,0,4)  = qc(i,1,4)
                  pc(i,0)    =pc(i,1)
                  cc(i,0)    =cc(i,1)
                  tt(i,0)    =tt(i,1)
                  tv(i,0)    =tv(i,1)
                  uc(i,0)    =-uc(i,1)
                  vc(i,0)    =-vc(i,1)
                  xc(i,0,1)  =xc(i,1,1)
                  xc(i,0,2)  =xc(i,1,2)
                  xc(i,0,3)  =xc(i,1,3)
                  xc(i,0,4)  =xc(i,1,4)
                  xc(i,0,5)  =xc(i,1,5)
                  vmue(i,0)  =vmue(i,1)
                  vkap(i,0)  =vkap(i,1)
                  vkav(i,0)  =vkav(i,1)
                  ds(i,0,1)  =ds(i,1,1)
                  ds(i,0,2)  =ds(i,1,2)
                  ds(i,0,3)  =ds(i,1,3)
                  ds(i,0,4)  =ds(i,1,4)
                  ds(i,0,5)  =ds(i,1,5)
                  hs(i,0,1)  =hs(i,1,1)
                  hs(i,0,2)  =hs(i,1,2)
                  hs(i,0,3)  =hs(i,1,3)
                  hs(i,0,4)  =hs(i,1,4)
                  hs(i,0,5)  =hs(i,1,5)
                  evs(i,0,1) =evs(i,1,1)
                  evs(i,0,2) =evs(i,1,2)
                  evs(i,0,3) =evs(i,1,3)
                  evs(i,0,4) =evs(i,1,4)
                  evs(i,0,5) =evs(i,1,5)
               end do
c isothermal wall
            else
               pai2v =0.5/acos(-1.0)
               cof =float(max(0,min(1,100-loop)))
               do i=1,imax-1
c wall values
                  ub        =0.0
                  vb        =0.0
                  pb        =pc(i,1)
                  ttb       =twall
                  tvb       =twall
                  akwo      =cefo*sqrt(8.314*ttb*amas1(1)*pai2v)
                  ako       =akwo/ds(i,1,1)
                  akwn      =cefn*sqrt(8.314*ttb*amas1(2)*pai2v)
                  akn       =akwn/ds(i,1,2)
                  dni       =0.5*vol(i,1)/sj(i,1,3)
                  cow       =xc(i,1,1)/(1.0+ako*dni)
                  cnw       =xc(i,1,2)/(1.0+akn*dni)
                  cnow      =xc(i,1,3)
                  co2w      =xc(i,1,4)+cow*ako*dni
                  cn2w      =1.0-cow-cnw-cnow-co2w
                  cimi      = cow*amas1(1)+ cnw*amas1(2)+cnow*amas1(3)
     &                      +co2w*amas1(4)+cn2w*amas1(5)
                  rob       =pb/(8.314*cimi*ttb)
                  rsb1      =cow *rob
                  rsb2      =cnw *rob
                  rsb3      =cnow*rob
                  rsb4      =co2w*rob
                  rsb5      =cn2w*rob
                  robv      =1.0/rob
                  xcb1      =robv*rsb1
                  xcb2      =robv*rsb2
                  xcb3      =robv*rsb3
                  xcb4      =robv*rsb4
                  xcb5      =robv*rsb5
                  bno       =rsb1*amas1(1)
                  bnn       =rsb2*amas1(2)
                  bnno      =rsb3*amas1(3)
                  bno2      =rsb4*amas1(4)
                  bnn2      =rsb5*amas1(5)
                  bnsum     =bno+bnn+bnno+bno2+bnn2
                  cvav      = bno*cvi(1)+ bnn*cvi(2)+bnno*cvi(3)
     &                      +bno2*cvi(4)+bnn2*cvi(5)
                  gm        =1.0+8.314*bnsum/cvav
                  xa        =exp(  226.0/tvb)
                  xb        =exp(-3785.5/tvb)
                  xa2       =xa*xa
                  xa4       =xa2*xa2
                  xa8       =xa4*xa4
                  xa10      =xa8*xa2
                  xa12      =xa10*xa2
                  xa15      =xa12*xa2*xa
                  xb2       =xb*xb
                  xb3       =xb2*xb
                  xb4       =xb2*xb2
                  xb5       =xb4*xb
                  xb6       =xb5*xb
                  xb7       =xb6*xb
                  ev1       =2260.0/(xa10-1.0)
                  ev2       =3390.0/(xa15-1.0)
                  ev3       =2712.0/(xa12-1.0)
                  ev4       =(22713.0*xb3+18927.5*xb5)/(3.0+2.0*xb3
     &                                                         +xb5)
                  ev5       =113565.0*xb6/(9.0+5.0*xb6)
                  ev6       =264985.0*xb7/(4.0+10.0*xb7)
                  evw       =8.314*(bno2*ev1+bnn2*ev2+bnno*ev3
     &                             +bno2*ev4+bno *ev5+bnn *ev6)
                  eh        =h0i(1)*bno+h0i(2)*bnn+h0i(3)*bnno
                  eb        =cvav*ttb+evw+eh
                  cb        =sqrt(gm*pb/rob)
                  sc        =0.5
c
                  uc(i,0)   =2.0*ub-uc(i,1)
                  vc(i,0)   =2.0*vb-vc(i,1)
                  tt(i,0)   =2.0*ttb-tt(i,1)
                  tv(i,0)   =2.0*tvb-tv(i,1)
                  pc(i,0)   =2.0*pb-pc(i,1)
                  qc(i,0,1) =2.0*rob-qc(i,1,1)
                  qs(i,0,1) =2.0*rsb1-qs(i,1,1)
                  qs(i,0,2) =2.0*rsb2-qs(i,1,2)
                  qs(i,0,3) =2.0*rsb3-qs(i,1,3)
                  qc(i,0,2) =qc(i,0,1)*uc(i,0)
                  qc(i,0,3) =qc(i,0,1)*vc(i,0)
                  xc(i,0,1) =2.0*xcb1-xc(i,1,1)
                  xc(i,0,2) =2.0*xcb2-xc(i,1,2)
                  xc(i,0,3) =2.0*xcb3-xc(i,1,3)
                  xc(i,0,4) =2.0*xcb4-xc(i,1,4)
                  xc(i,0,5) =2.0*xcb5-xc(i,1,5)
                  qc(i,0,4) =2.0*eb-qc(i,1,4)
                  qv(i,0,1) =2.0*evw-qv(i,1,1)
                  cc(i,0)   =2.0*cb-cc(i,1)
c
                  rs1       =8.314*amas1(1)
                  cvts1     =1.5*rs1
                  cvrs1     =0.0
                  cvvs1     =0.0
                  rs2       =8.314*amas1(2)
                  cvts2     =1.5*rs2
                  cvrs2     =0.0
                  cvvs2     =0.0
                  rs3       =8.314*amas1(3)
                  cvts3     =1.5*rs3
                  cvrs3     =rs3
                  cvvs3     =rs3
                  rs4       =8.314*amas1(4)
                  cvts4     =1.5*rs4
                  cvrs4     =rs4
                  cvvs4     =rs4
                  rs5       =8.314*amas1(5)
                  cvts5     =1.5*rs5
                  cvrs5     =rs5
                  cvvs5     =rs5
                  if(nvis.eq.0) then
                     am        =1.0/(xcb1*amas1(1)+xcb2*amas1(2)
     &                              +xcb3*amas1(3)+xcb4*amas1(4)
     &                              +xcb5*amas1(5))
                     alnt      =log(ttb)
                     x1        =xcb1*am*amas1(1)

                     vmus1     =0.1*exp((a1*alnt+b1)*alnt+c1)
                     vkas1     =vmus1*(2.5*cvts1+cvrs1)
                     vkvs1     =vmus1*cvvs1
                     x2        =xcb2*am*amas1(2)
                     vmus2     =0.1*exp((a2*alnt+b2)*alnt+c2)
                     vkas2     =vmus2*(2.5*cvts2+cvrs2)
                     vkvs2     =vmus2*cvvs2
                     x3        =xcb3*am*amas1(3)
                     vmus3     =0.1*exp((a3*alnt+b3)*alnt+c3)
                     vkas3     =vmus3*(2.5*cvts3+cvrs3)
                     vkvs3     =vmus3*cvvs3
                     x4        =xcb4*am*amas1(4)
                     vmus4     =0.1*exp((a4*alnt+b4)*alnt+c4)
                     vkas4     =vmus4*(2.5*cvts4+cvrs4)
                     vkvs4     =vmus4*cvvs4
                     x5        =xcb5*am*amas1(5)
                     vmus5     =0.1*exp((a5*alnt+b5)*alnt+c5)
                     vkas5     =vmus5*(2.5*cvts5+cvrs5)
                     vkvs5     =vmus5*cvvs5
                     phs11t    =1.0+sqrt(vmus1/vmus1)
     &                          *(amass(1)*amas1(1))**0.25
                     phs11     =phs11t*phs11t
     &                          /sqrt(8.0*(1.0+amass(1)*amas1(1)))
                     phs12t    =1.0+sqrt(vmus1/vmus2)
     &                          *(amass(2)*amas1(1))**0.25
                     phs12     =phs12t*phs12t
     &                          /sqrt(8.0*(1.0+amass(1)*amas1(2)))
                     phs13t    =1.0+sqrt(vmus1/vmus3)
     &                          *(amass(3)*amas1(1))**0.25
                     phs13     =phs13t*phs13t
     &                          /sqrt(8.0*(1.0+amass(1)*amas1(3)))
                     phs14t    =1.0+sqrt(vmus1/vmus4)
     &                          *(amass(4)*amas1(1))**0.25
                     phs14     =phs14t*phs14t
     &                          /sqrt(8.0*(1.0+amass(1)*amas1(4)))
                     phs15t    =1.0+sqrt(vmus1/vmus5)
     &                          *(amass(5)*amas1(1))**0.25
                     phs15     =phs15t*phs15t
     &                          /sqrt(8.0*(1.0+amass(1)*amas1(5)))
                     phis1     =1.0/(x1*phs11+x2*phs12+x3*phs13
     &                           +x4*phs14+x5*phs15)
                     phs21t    =1.0+sqrt(vmus2/vmus1)
     &                          *(amass(1)*amas1(2))**0.25
                     phs21     =phs21t*phs21t
     &                          /sqrt(8.0*(1.0+amass(2)*amas1(1)))
                     phs22t    =1.0+sqrt(vmus2/vmus2)
     &                          *(amass(2)*amas1(2))**0.25
                     phs22     =phs22t*phs22t
     &                          /sqrt(8.0*(1.0+amass(2)*amas1(2)))
                     phs23t    =1.0+sqrt(vmus2/vmus3)
     &                          *(amass(3)*amas1(2))**0.25
                     phs23     =phs23t*phs23t
     &                          /sqrt(8.0*(1.0+amass(2)*amas1(3)))
                     phs24t    =1.0+sqrt(vmus2/vmus4)
     &                          *(amass(4)*amas1(2))**0.25
                     phs24     =phs24t*phs24t
     &                          /sqrt(8.0*(1.0+amass(2)*amas1(4)))
                     phs25t    =1.0+sqrt(vmus2/vmus5)
     &                          *(amass(5)*amas1(2))**0.25
                     phs25     =phs25t*phs25t
     &                          /sqrt(8.0*(1.0+amass(2)*amas1(5)))
                     phis2     =1.0/(x1*phs21+x2*phs22+x3*phs23
     &                           +x4*phs24+x5*phs25)
                     phs31t    =1.0+sqrt(vmus3/vmus1)
     &                          *(amass(1)*amas1(3))**0.25
                     phs31     =phs31t*phs31t
     &                          /sqrt(8.0*(1.0+amass(3)*amas1(1)))
                     phs32t    =1.0+sqrt(vmus3/vmus2)
     &                          *(amass(2)*amas1(3))**0.25
                     phs32     =phs32t*phs32t
     &                          /sqrt(8.0*(1.0+amass(3)*amas1(2)))
                     phs33t    =1.0+sqrt(vmus3/vmus3)
     &                          *(amass(3)*amas1(3))**0.25
                     phs33     =phs33t*phs33t
     &                          /sqrt(8.0*(1.0+amass(3)*amas1(3)))
                     phs34t    =1.0+sqrt(vmus3/vmus4)
     &                          *(amass(4)*amas1(3))**0.25
                     phs34     =phs34t*phs34t
     &                          /sqrt(8.0*(1.0+amass(3)*amas1(4)))
                     phs35t    =1.0+sqrt(vmus3/vmus5)
     &                          *(amass(5)*amas1(3))**0.25
                     phs35     =phs35t*phs35t
     &                          /sqrt(8.0*(1.0+amass(3)*amas1(5)))
                     phis3     =1.0/(x1*phs31+x2*phs32+x3*phs33
     &                           +x4*phs34+x5*phs35)
                     phs41t    =1.0+sqrt(vmus4/vmus1)
     &                          *(amass(1)*amas1(4))**0.25
                     phs41     =phs41t*phs41t
     &                          /sqrt(8.0*(1.0+amass(4)*amas1(1)))
                     phs42t    =1.0+sqrt(vmus4/vmus2)
     &                          *(amass(2)*amas1(4))**0.25
                     phs42     =phs42t*phs42t
     &                          /sqrt(8.0*(1.0+amass(4)*amas1(2)))
                     phs43t    =1.0+sqrt(vmus4/vmus3)
     &                          *(amass(3)*amas1(4))**0.25
                     phs43     =phs43t*phs43t
     &                          /sqrt(8.0*(1.0+amass(4)*amas1(3)))
                     phs44t    =1.0+sqrt(vmus4/vmus4)
     &                          *(amass(4)*amas1(4))**0.25
                     phs44     =phs44t*phs44t
     &                          /sqrt(8.0*(1.0+amass(4)*amas1(4)))
                     phs45t    =1.0+sqrt(vmus4/vmus5)
     &                          *(amass(5)*amas1(4))**0.25
                     phs45     =phs45t*phs45t
     &                          /sqrt(8.0*(1.0+amass(4)*amas1(5)))
                     phis4     =1.0/(x1*phs41+x2*phs42+x3*phs43
     &                           +x4*phs44+x5*phs45)
                     phs51t    =1.0+sqrt(vmus5/vmus1)
     &                          *(amass(1)*amas1(5))**0.25
                     phs51     =phs51t*phs51t
     &                          /sqrt(8.0*(1.0+amass(5)*amas1(1)))
                     phs52t    =1.0+sqrt(vmus5/vmus2)
     &                          *(amass(2)*amas1(5))**0.25
                     phs52     =phs52t*phs52t
     &                          /sqrt(8.0*(1.0+amass(5)*amas1(2)))
                     phs53t    =1.0+sqrt(vmus5/vmus3)
     &                          *(amass(3)*amas1(5))**0.25
                     phs53     =phs53t*phs53t
     &                          /sqrt(8.0*(1.0+amass(5)*amas1(3)))
                     phs54t    =1.0+sqrt(vmus5/vmus4)
     &                          *(amass(4)*amas1(5))**0.25
                     phs54     =phs54t*phs54t
     &                          /sqrt(8.0*(1.0+amass(5)*amas1(4)))
                     phs55t    =1.0+sqrt(vmus5/vmus5)
     &                          *(amass(5)*amas1(5))**0.25
                     phs55     =phs55t*phs55t
     &                          /sqrt(8.0*(1.0+amass(5)*amas1(5)))
                     phis5     =1.0/(x1*phs51+x2*phs52+x3*phs53
     &                           +x4*phs54+x5*phs55)
                     vmueb     =x1*vmus1*phis1+x2*vmus2*phis2
     &                         +x3*vmus3*phis3+x4*vmus4*phis4
     &                         +x5*vmus5*phis5
                     vkapb     =x1*vkas1*phis1+x2*vkas2*phis2
     &                         +x3*vkas3*phis3+x4*vkas4*phis4
     &                         +x5*vkas5*phis5
                     vkavb     =x1*vkvs1*phis1+x2*vkvs2*phis2
     &                         +x3*vkvs3*phis3+x4*vkvs4*phis4
     &                         +x5*vkvs5*phis5
                     ds1b      =vmueb/(rob*sc)
                     ds2b      =ds1b
                     ds3b      =ds1b
                     ds4b      =ds1b
                     ds5b      =ds1b
                     vmue(i,0) =2.0*vmueb-vmue(i,1)
                     vkap(i,0) =2.0*vkapb-vkap(i,1)
                     vkav(i,0) =2.0*vkavb-vkav(i,1)
                     ds(i,0,1) =2.0*ds1b-ds(i,1,1)
                     ds(i,0,2) =2.0*ds2b-ds(i,1,2)
                     ds(i,0,3) =2.0*ds3b-ds(i,1,3)
                     ds(i,0,4) =2.0*ds4b-ds(i,1,4)
                     ds(i,0,5) =2.0*ds5b-ds(i,1,5)
                  end if
                  if(nvis.eq.1) then
c Collisonal integrals by Gupta curv fitting
c Diffusion coefficent is determined by binary diffusion coefficient
c 1:O, 2:N, 3:NO, 4:O2, 5:N2
                     pi    = dacos(-1.0d0)
                     avoga1= 1.d0/6.02214d23 ! Avogadro constant [1/mol]
                     boltz = 1.3806503d-23   ! Boltzman constant [m2kg/s2K]
                     gasc  = 8.314d0
                     dct1  = 8.d0/3.d0       ! constant for delta_sr(1)
                     dct2  = 16.d0/5.d0      ! constant for delta_sr(2)
                     te    =tt(i,0)
                     ttv   =tv(i,0)
                     pp    =pc(i,0)
                     alte  =dlog(te)
                     alt2  =alte*alte
                     prt1  =1.d0/pi/te/gasc
                     altv  =dlog(ttv)
                     altv2 =altv*altv
                     prtv  =1.d0/pi/ttv/gasc
                     do n=1,5
                        xr(n)=xc(i,0,n)*amas1(n)
                     end do
                     xt=0.d0
                     do n=1,5
                        xt=xt+xr(n)
                     end do
c pi*Omega_sr(1,1)(=pomega1) and Omega_sr(2,2)(=pomega2)
                     do m=1,5   ! s
                        do l=1,5 ! r
c unit A^2(=10^{-20}[m])
                           pomega1   =dexp(dog(l,m,1))*te**
     &                                    (aog(l,m,1)*alt2
     &                                    +bog(l,m,1)*alte+cog(l,m,1))
                           pomega2   =dexp(dog(l,m,2))*te**
     &                                    (aog(l,m,2)*alt2
     &                                    +bog(l,m,2)*alte+cog(l,m,2))
                           pomegav   =dexp(dog(l,m,1))*ttv**
     &                                    (aog(l,m,1) *altv2
     &                                    +bog(l,m,1)*altv+cog(l,m,1))
c delta_sr^1(=dlsr1) and delta_sr^2(=dlsr2)
                           dlsr1(l,m)=dct1*dsqrt(2.d0*amass(m)*amass(l)
     &                        /(amass(m)+amass(l))*prt1)*pomega1*1.d-20
                           dlsr2(l,m)=dct2*dsqrt(2.d0*amass(m)*amass(l)
     &                        /(amass(m)+amass(l))*prt1)*pomega2*1.d-20
                           dlsrv(l,m)=dct1*dsqrt(2.d0*amass(m)*amass(l)
     &                        /(amass(m)+amass(l))*prtv)*pomegav*1.d-20
c alpha_sr^1(=alsr)
                           ramass    =amass(m)/amass(l)
                           alsr(l,m) = 1.d0 + ((1.d0-ramass)*(0.45d0
     &                              -2.54d0*ramass)/((1.d0+ramass)**2))
                           ddsr(l,m)  = (boltz*te)/(pp*dlsr1(l,m))
                        end do
                     end do
c
                     xxdl1(:)=0.d0
                     xxdl2(:)=0.d0
                     xxdlv(:)=0.d0
                     axdl2(:)=0.d0
                     do m = 1,5 ! s
                        do l = 1,5 ! r
                           xxdl1(m) = xxdl1(m)+ xr(l)*dlsr1(l,m)
                           xxdl2(m) = xxdl2(m)+ xr(l)*dlsr2(l,m)
                           xxdlv(m) = xxdlv(m)+ xr(l)*dlsrv(l,m)
                           axdl2(m) = axdl2(m)+ xr(l)*dlsr2(l,m)
     &                               *alsr(l,m)
                        end do
                     end do
                     vics=0.d0
                     thct=0.d0
                     thcr=0.d0
                     thcv=0.d0
                     do m=1,5   !O,N,NO,O2,N2
                        vics=vics+xc(i,0,m)/xxdl2(m)
                        thct=thct+xr(m)/axdl2(m)
                     end do
                     do m=3,5   ! NO,O2,N2
                        thcr=thcr+xr(m)/xxdl1(m)
                        thcv=thcv+xr(m)/xxdlv(m)
                     end do
                     vmue(i,0)=vics*avoga1
                     vkap(i,0)=boltz*(3.75d0*thct+thcr)
                     vkav(i,0)=boltz*thcv
c binary diffusion coefficient
                     dsb(:)=0.d0
                     do m=1,5
                        do l=1,5
                           if(m.ne.l) then
                              dsb(m) = dsb(m) + xr(l)/ddsr(l,m)
                           end if
                        end do
                     end do
                     do m=1,5
                        ds(i,0,m) = ((xt**2)*amass(m)
     &                            *(1.d0-xc(i,0,m)))/dsb(m)
                     end do
                  end if
c set
                  tv1       =1.0/tvb
                  xa        =exp(226.0*tv1)
                  xb        =exp(-3785.5*tv1)
                  xa2       =xa*xa
                  xa4       =xa2*xa2
                  xa8       =xa4*xa4
                  xa10      =xa8*xa2
                  xa12      =xa10*xa2
                  xa15      =xa12*xa2*xa
                  xb2       =xb*xb
                  xb3       =xb2*xb
                  xb4       =xb2*xb2
                  xb5       =xb4*xb
                  xb6       =xb5*xb
                  xb7       =xb6*xb
                  ev1       =2260.0/(xa10-1.0)
                  ev2       =3390.0/(xa15-1.0)
                  ev3       =2712.0/(xa12-1.0)
                  ev4       =(22713.0*xb3+18927.5*xb5)
     &                      /(3.0+2.0*xb3+xb5)
                  ev5       =113565.0*xb6/(9.0+5.0*xb6)
                  ev6       =264985.0*xb7/(4.0+10.0*xb7)
                  ev1g      =8.314*ev1
                  ev2g      =8.314*ev2
                  ev3g      =8.314*ev3
                  ev4g      =8.314*ev4
                  ev5g      =8.314*ev5
                  ev6g      =8.314*ev6
                  evs1b     =      ev5g *amas1(1)
                  evs2b     =      ev6g *amas1(2)
                  evs3b     =      ev3g *amas1(3)
                  evs4b     =(ev1g+ev4g)*amas1(4)
                  evs5b     =      ev2g *amas1(5)
                  hs1b      =(cvts1+cvrs1+rs1)*ttb
     &                      +(h0i(1)+ev5g)*amas1(1)
                  hs2b      =(cvts2+cvrs2+rs2)*ttb
     &                      +(h0i(2)+ev6g)*amas1(2)
                  hs3b      =(cvts3+cvrs3+rs3)*ttb
     &                      +(h0i(3)+ev3g)*amas1(3)
                  hs4b      =(cvts4+cvrs4+rs4)*ttb
     &                      +(  ev1g+ev4g)*amas1(4)
                  hs5b      =(cvts5+cvrs5+rs5)*ttb
     &                      +        ev2g *amas1(5)
                  evs(i,0,1)=2.0*evs1b-evs(i,1,1)
                  evs(i,0,2)=2.0*evs2b-evs(i,1,2)
                  evs(i,0,3)=2.0*evs3b-evs(i,1,3)
                  evs(i,0,4)=2.0*evs4b-evs(i,1,4)
                  evs(i,0,5)=2.0*evs5b-evs(i,1,5)
                  hs(i,0,1) =2.0*hs1b-hs(i,1,1)
                  hs(i,0,2) =2.0*hs2b-hs(i,1,2)
                  hs(i,0,3) =2.0*hs3b-hs(i,1,3)
                  hs(i,0,4) =2.0*hs4b-hs(i,1,4)
                  hs(i,0,5) =2.0*hs5b-hs(i,1,5)
               end do
            endif
         endif
      endif
c* termination
      return
      end
c
      subroutine musl(mode)
c     ************************************************************
c     *     preprocessing by MUSCL                               *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /aray3/ qcl(nx,ny,4),qcr(nx,ny,4),
     &               qsl(nx,ny,3),qsr(nx,ny,3),
     &               qvl(nx,ny,1),qvr(nx,ny,1)
      common /aray5/ ish(0:nx,0:ny)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /coor3/ xg(0:nx,0:ny),yg(0:nx,0:ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /cntl5/ nprc,iprc
      common /cntl6/ nacu
      common /wbcd1/ twall,cefn,cefo
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
c* first order scheme
      if(mode.eq.1) then
         do j=1,jmax-1
            do i=1,imax-1
               qsl(i+1,j,1) =qs(i,j,1)
               qsl(i+1,j,2) =qs(i,j,2)
               qsl(i+1,j,3) =qs(i,j,3)
               qvl(i+1,j,1) =qv(i,j,1)
               qcl(i+1,j,1) =qc(i,j,1)
               qcl(i+1,j,2) =uc(i,j  )
               qcl(i+1,j,3) =vc(i,j  )
               qcl(i+1,j,4) =pc(i,j  )
               qsr(i  ,j,1) =qs(i,j,1)
               qsr(i  ,j,2) =qs(i,j,2)
               qsr(i  ,j,3) =qs(i,j,3)
               qvr(i  ,j,1) =qv(i,j,1)
               qcr(i  ,j,1) =qc(i,j,1)
               qcr(i  ,j,2) =uc(i,j  )
               qcr(i  ,j,3) =vc(i,j  )
               qcr(i  ,j,4) =pc(i,j  )
            end do
         end do
c*boundary
         do j=1,jmax-1
c*axis boundary condition (velocty mirror)
            qsl(1   ,j,1) = qsr(1,j,1)
            qsl(1   ,j,2) = qsr(1,j,2)
            qsl(1   ,j,3) = qsr(1,j,3)
            qvl(1   ,j,1) = qvr(1,j,1)
            qcl(1   ,j,1) = qcr(1,j,1)
            qcl(1   ,j,2) = qcr(1,j,2)
            qcl(1   ,j,3) =-qcr(1,j,3)
            qcl(1   ,j,4) = qcr(1,j,4)
c*outflow boundary condition (supersonic outflow)
            qsr(imax,j,1) = qsl(imax,j,1)
            qsr(imax,j,2) = qsl(imax,j,2)
            qsr(imax,j,3) = qsl(imax,j,3)
            qvr(imax,j,1) = qvl(imax,j,1)
            qcr(imax,j,1) = qcl(imax,j,1)
            qcr(imax,j,2) = qcl(imax,j,2)
            qcr(imax,j,3) = -qcl(imax,j,3)
            qcr(imax,j,4) = qcl(imax,j,4)
         end do
c
      elseif(mode.eq.2) then
         do i=1,imax-1
            do j=1,jmax-1
               qsl(i,j+1,1) =qs(i,j,1)
               qsl(i,j+1,2) =qs(i,j,2)
               qsl(i,j+1,3) =qs(i,j,3)
               qvl(i,j+1,1) =qv(i,j,1)
               qcl(i,j+1,1) =qc(i,j,1)
               qcl(i,j+1,2) =uc(i,j  )
               qcl(i,j+1,3) =vc(i,j  )
               qcl(i,j+1,4) =pc(i,j  )
               qsr(i,j  ,1) =qs(i,j,1)
               qsr(i,j  ,2) =qs(i,j,2)
               qsr(i,j  ,3) =qs(i,j,3)
               qvr(i,j  ,1) =qv(i,j,1)
               qcr(i,j  ,1) =qc(i,j,1)
               qcr(i,j  ,2) =uc(i,j  )
               qcr(i,j  ,3) =vc(i,j  )
               qcr(i,j  ,4) =pc(i,j  )
            end do
         end do
c
c*boundary
         do i=1,imax-1
c*wall boundary condition (for Euler, slip wall)
            if (nprc.eq.0) then
               sjv           =1.0/sj(i,1,3)
               sjx           =sjv*sj(i,1,1)
               sjy           =sjv*sj(i,1,2)
               qn            =2.0*(qcr(i,1,2)*sjx+qcr(i,1,3)*sjy)
               qsl(i,1   ,1) =qsr(i,1,1)
               qsl(i,1   ,2) =qsr(i,1,2)
               qsl(i,1   ,3) =qsr(i,1,3)
               qvl(i,1   ,1) =qvr(i,1,1)
               qcl(i,1   ,1) =qcr(i,1,1)
               qcl(i,1   ,2) =qcr(i,1,2)-qn*sjx
               qcl(i,1   ,3) =qcr(i,1,3)-qn*sjy
               qcl(i,1   ,4) =qcr(i,1,4)
c*wall boundary condition (for NS, non-slip wall)
            else
               qsl(i,1   ,1) =qsr(i,1,1)
               qsl(i,1   ,2) =qsr(i,1,2)
               qsl(i,1   ,3) =qsr(i,1,3)
               qvl(i,1   ,1) =qvr(i,1,1)
               qcl(i,1   ,1) =qcr(i,1,1)
               qcl(i,1   ,2) =-qcr(i,1,2)
               qcl(i,1   ,3) =-qcr(i,1,3)
               qcl(i,1   ,4) =qcr(i,1,4)
            end if
c*inflow boundary confition
c 流入・流出を x 座標によって変える (x<0:流入)
            if(xg(i, jmax).le.0.d0)then
               qsr(i,jmax,1) =r1inf
               qsr(i,jmax,2) =r2inf
               qsr(i,jmax,3) =r3inf
               qvr(i,jmax,1) =evinf
               qcr(i,jmax,1) =rinf
               qcr(i,jmax,2) =uinf
               qcr(i,jmax,3) =vinf
               qcr(i,jmax,4) =pinf
            else
               qsr(i,jmax,1) = qsr(i,jmax-1,1)
               qsr(i,jmax,2) = qsr(i,jmax-1,2)
               qsr(i,jmax,3) = qsr(i,jmax-1,3)
               qvr(i,jmax,1) = qvr(i,jmax-1,1)
               qcr(i,jmax,1) = qcr(i,jmax-1,1)
               qcr(i,jmax,2) = qcr(i,jmax-1,2)
               qcr(i,jmax,3) = qcr(i,jmax-1,3)
               qcr(i,jmax,4) = qcr(i,jmax-1,4)
            end if
         end do
      endif
      if(nacu.eq.1) return
c
c* higher order scheme
      if(mode.eq.1) then
         do j=1,jmax-1
            do i=1,imax-1
               ds1p        =qs(i+1,j,1)-qs(i  ,j,1)
               ds2p        =qs(i+1,j,2)-qs(i  ,j,2)
               ds3p        =qs(i+1,j,3)-qs(i  ,j,3)
               dv1p        =qv(i+1,j,1)-qv(i  ,j,1)
               dc1p        =qc(i+1,j,1)-qc(i  ,j,1)
               dc2p        =uc(i+1,j  )-uc(i  ,j  )
               dc3p        =vc(i+1,j  )-vc(i  ,j  )
               dc4p        =pc(i+1,j  )-pc(i  ,j  )
               ds1m        =qs(i  ,j,1)-qs(i-1,j,1)
               ds2m        =qs(i  ,j,2)-qs(i-1,j,2)
               ds3m        =qs(i  ,j,3)-qs(i-1,j,3)
               dv1m        =qv(i  ,j,1)-qv(i-1,j,1)
               dc1m        =qc(i  ,j,1)-qc(i-1,j,1)
               dc2m        =uc(i  ,j  )-uc(i-1,j  )
               dc3m        =vc(i  ,j  )-vc(i-1,j  )
               dc4m        =pc(i  ,j  )-pc(i-1,j  )
               sgns1       =sign(0.5d0,ds1m)+sign(0.5d0,ds1p)
               sgns2       =sign(0.5d0,ds2m)+sign(0.5d0,ds2p)
               sgns3       =sign(0.5d0,ds3m)+sign(0.5d0,ds3p)
               sgnv1       =sign(0.5d0,dv1m)+sign(0.5d0,dv1p)
               sgnc1       =sign(0.5d0,dc1m)+sign(0.5d0,dc1p)
               sgnc2       =sign(0.5d0,dc2m)+sign(0.5d0,dc2p)
               sgnc3       =sign(0.5d0,dc3m)+sign(0.5d0,dc3p)
               sgnc4       =sign(0.5d0,dc4m)+sign(0.5d0,dc4p)
               es1p        =sgns1*min(abs(ds1m),0.5*abs(ds1p))
               es2p        =sgns2*min(abs(ds2m),0.5*abs(ds2p))
               es3p        =sgns3*min(abs(ds3m),0.5*abs(ds3p))
               ev1p        =sgnv1*min(abs(dv1m),0.5*abs(dv1p))
               ec1p        =sgnc1*min(abs(dc1m),0.5*abs(dc1p))
               ec2p        =sgnc2*min(abs(dc2m),0.5*abs(dc2p))
               ec3p        =sgnc3*min(abs(dc3m),0.5*abs(dc3p))
               ec4p        =sgnc4*min(abs(dc4m),0.5*abs(dc4p))
               es1m        =sgns1*min(abs(ds1p),0.5*abs(ds1m))
               es2m        =sgns2*min(abs(ds2p),0.5*abs(ds2m))
               es3m        =sgns3*min(abs(ds3p),0.5*abs(ds3m))
               ev1m        =sgnv1*min(abs(dv1p),0.5*abs(dv1m))
               ec1m        =sgnc1*min(abs(dc1p),0.5*abs(dc1m))
               ec2m        =sgnc2*min(abs(dc2p),0.5*abs(dc2m))
               ec3m        =sgnc3*min(abs(dc3p),0.5*abs(dc3m))
               ec4m        =sgnc4*min(abs(dc4p),0.5*abs(dc4m))
               qsl(i+1,j,1) =qs(i,j,1)+es1p
               qsl(i+1,j,2) =qs(i,j,2)+es2p
               qsl(i+1,j,3) =qs(i,j,3)+es3p
               qvl(i+1,j,1) =qv(i,j,1)+ev1p
               qcl(i+1,j,1) =qc(i,j,1)+ec1p
               qcl(i+1,j,2) =uc(i,j  )+ec2p
               qcl(i+1,j,3) =vc(i,j  )+ec3p
               qcl(i+1,j,4) =pc(i,j  )+ec4p
               qsr(i  ,j,1) =qs(i,j,1)-es1m
               qsr(i  ,j,2) =qs(i,j,2)-es2m
               qsr(i  ,j,3) =qs(i,j,3)-es3m
               qvr(i  ,j,1) =qv(i,j,1)-ev1m
               qcr(i  ,j,1) =qc(i,j,1)-ec1m
               qcr(i  ,j,2) =uc(i,j  )-ec2m
               qcr(i  ,j,3) =vc(i,j  )-ec3m
               qcr(i  ,j,4) =pc(i,j  )-ec4m
            end do
         end do
c*boundary
         do j=1,jmax-1
c*axis boundary condition (velocity mirror)
            qsl(1   ,j,1) = qsr(1,j,1)
            qsl(1   ,j,2) = qsr(1,j,2)
            qsl(1   ,j,3) = qsr(1,j,3)
            qvl(1   ,j,1) = qvr(1,j,1)
            qcl(1   ,j,1) = qcr(1,j,1)
            qcl(1   ,j,2) = qcr(1,j,2)
            qcl(1   ,j,3) =-qcr(1,j,3)
            qcl(1   ,j,4) = qcr(1,j,4)
c*outflow boundary condition (supersonic outflow)
            qsr(imax,j,1) = qsl(imax,j,1)
            qsr(imax,j,2) = qsl(imax,j,2)
            qsr(imax,j,3) = qsl(imax,j,3)
            qvr(imax,j,1) = qvl(imax,j,1)
            qcr(imax,j,1) = qcl(imax,j,1)
            qcr(imax,j,2) = qcl(imax,j,2)
            qcr(imax,j,3) = -qcl(imax,j,3)
            qcr(imax,j,4) = qcl(imax,j,4)
         end do
c
      elseif(mode.eq.2) then
         do j=1,jmax-1
            do i=1,imax-1
               ds1p        =qs(i,j+1,1)-qs(i,j  ,1)
               ds2p        =qs(i,j+1,2)-qs(i,j  ,2)
               ds3p        =qs(i,j+1,3)-qs(i,j  ,3)
               dv1p        =qv(i,j+1,1)-qv(i,j  ,1)
               dc1p        =qc(i,j+1,1)-qc(i,j  ,1)
               dc2p        =uc(i,j+1  )-uc(i,j    )
               dc3p        =vc(i,j+1  )-vc(i,j    )
               dc4p        =pc(i,j+1  )-pc(i,j    )
               ds1m        =qs(i,j  ,1)-qs(i,j-1,1)
               ds2m        =qs(i,j  ,2)-qs(i,j-1,2)
               ds3m        =qs(i,j  ,3)-qs(i,j-1,3)
               dv1m        =qv(i,j  ,1)-qv(i,j-1,1)
               dc1m        =qc(i,j  ,1)-qc(i,j-1,1)
               dc2m        =uc(i,j    )-uc(i,j-1  )
               dc3m        =vc(i,j    )-vc(i,j-1  )
               dc4m        =pc(i,j    )-pc(i,j-1  )
               sgns1       =sign(0.5d0,ds1m)+sign(0.5d0,ds1p)
               sgns2       =sign(0.5d0,ds2m)+sign(0.5d0,ds2p)
               sgns3       =sign(0.5d0,ds3m)+sign(0.5d0,ds3p)
               sgnv1       =sign(0.5d0,dv1m)+sign(0.5d0,dv1p)
               sgnc1       =sign(0.5d0,dc1m)+sign(0.5d0,dc1p)
               sgnc2       =sign(0.5d0,dc2m)+sign(0.5d0,dc2p)
               sgnc3       =sign(0.5d0,dc3m)+sign(0.5d0,dc3p)
               sgnc4       =sign(0.5d0,dc4m)+sign(0.5d0,dc4p)
               es1p        =sgns1*min(abs(ds1m),0.5*abs(ds1p))
               es2p        =sgns2*min(abs(ds2m),0.5*abs(ds2p))
               es3p        =sgns3*min(abs(ds3m),0.5*abs(ds3p))
               ev1p        =sgnv1*min(abs(dv1m),0.5*abs(dv1p))
               ec1p        =sgnc1*min(abs(dc1m),0.5*abs(dc1p))
               ec2p        =sgnc2*min(abs(dc2m),0.5*abs(dc2p))
               ec3p        =sgnc3*min(abs(dc3m),0.5*abs(dc3p))
               ec4p        =sgnc4*min(abs(dc4m),0.5*abs(dc4p))
               es1m        =sgns1*min(abs(ds1p),0.5*abs(ds1m))
               es2m        =sgns2*min(abs(ds2p),0.5*abs(ds2m))
               es3m        =sgns3*min(abs(ds3p),0.5*abs(ds3m))
               ev1m        =sgnv1*min(abs(dv1p),0.5*abs(dv1m))
               ec1m        =sgnc1*min(abs(dc1p),0.5*abs(dc1m))
               ec2m        =sgnc2*min(abs(dc2p),0.5*abs(dc2m))
               ec3m        =sgnc3*min(abs(dc3p),0.5*abs(dc3m))
               ec4m        =sgnc4*min(abs(dc4p),0.5*abs(dc4m))
               qsl(i,j+1,1) =qs(i,j,1)+es1p
               qsl(i,j+1,2) =qs(i,j,2)+es2p
               qsl(i,j+1,3) =qs(i,j,3)+es3p
               qvl(i,j+1,1) =qv(i,j,1)+ev1p
               qcl(i,j+1,1) =qc(i,j,1)+ec1p
               qcl(i,j+1,2) =uc(i,j  )+ec2p
               qcl(i,j+1,3) =vc(i,j  )+ec3p
               qcl(i,j+1,4) =pc(i,j  )+ec4p
               qsr(i,j  ,1) =qs(i,j,1)-es1m
               qsr(i,j  ,2) =qs(i,j,2)-es2m
               qsr(i,j  ,3) =qs(i,j,3)-es3m
               qvr(i,j  ,1) =qv(i,j,1)-ev1m
               qcr(i,j  ,1) =qc(i,j,1)-ec1m
               qcr(i,j  ,2) =uc(i,j  )-ec2m
               qcr(i,j  ,3) =vc(i,j  )-ec3m
               qcr(i,j  ,4) =pc(i,j  )-ec4m
            end do
         end do
c*boundary
         do i=1,imax-1
c*wall boundary condition (for Euler, slip wall)
            if (nprc.eq.0) then
              sjv           =1.0/sj(i,1,3)
                sjx           =sjv*sj(i,1,1)
               sjy           =sjv*sj(i,1,2)
               qn            =2.0*(qcr(i,1,2)*sjx+qcr(i,1,3)*sjy)
               qsl(i,1   ,1) =qsr(i,1,1)
               qsl(i,1   ,2) =qsr(i,1,2)
               qsl(i,1   ,3) =qsr(i,1,3)
               qvl(i,1   ,1) =qvr(i,1,1)
               qcl(i,1   ,1) =qcr(i,1,1)
               qcl(i,1   ,2) =qcr(i,1,2)-qn*sjx
               qcl(i,1   ,3) =qcr(i,1,3)-qn*sjy
               qcl(i,1   ,4) =qcr(i,1,4)
c*wall boundary condition (for NS, non-slip wall)
            else
               qsl(i,1   ,1) =qsr(i,1,1)
               qsl(i,1   ,2) =qsr(i,1,2)
               qsl(i,1   ,3) =qsr(i,1,3)
               qvl(i,1   ,1) =qvr(i,1,1)
               qcl(i,1   ,1) =qcr(i,1,1)
               qcl(i,1   ,2) =-qcr(i,1,2)
               qcl(i,1   ,3) =-qcr(i,1,3)
               qcl(i,1   ,4) =qcr(i,1,4)
            end if
c*inflow boudary condition
! 流入・流出を x 座標によって変える
            if (xg(i,jmax).le.0.0d0) then
               qsr(i,jmax,1) =r1inf
               qsr(i,jmax,2) =r2inf
               qsr(i,jmax,3) =r3inf
               qvr(i,jmax,1) =evinf
               qcr(i,jmax,1) =rinf
               qcr(i,jmax,2) =uinf
               qcr(i,jmax,3) =vinf
               qcr(i,jmax,4) =pinf
            else
               qsr(i,jmax,1) = qsr(i,jmax-1,1)
               qsr(i,jmax,2) = qsr(i,jmax-1,2)
               qsr(i,jmax,3) = qsr(i,jmax-1,3)
               qvr(i,jmax,1) = qvr(i,jmax-1,1)
               qcr(i,jmax,1) = qcr(i,jmax-1,1)
               qcr(i,jmax,2) = qcr(i,jmax-1,2)
               qcr(i,jmax,3) = qcr(i,jmax-1,3)
               qcr(i,jmax,4) = qcr(i,jmax-1,4)
            end if
         end do
      endif
c* termination
      return
      end
c
      subroutine shock
c     *********************************************************
c     *     identification of shock-fix position              *
c     *********************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /aray5/ ish(0:nx,0:ny)
      common /coor1/ imax,jmax
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
c* reset array
      do j=0,jmax
         do i=0,imax
            ish(i,j) =0
         end do
      end do
c* in-normal sweep
      do j=2,jmax-1
         do i=1,imax-1
            sjv        =1.0/sj(i,j,3)
            unl        =sjv*(sj(i,j,1)*uc(i,j-1)+sj(i,j,2)*vc(i,j-1))
            unr        =sjv*(sj(i,j,1)*uc(i,j  )+sj(i,j,2)*vc(i,j  ))
            dd1l       =max(0.0d0, sign(1.0d0,unl-cc(i,j-1)))
            dd1r       =max(0.0d0,-sign(1.0d0,unr-cc(i,j)))
            nd1        =nint(dd1l*dd1r)
            dd2l       =max(0.0d0, sign(1.0d0,unl+cc(i,j-1)))
            dd2r       =max(0.0d0,-sign(1.0d0,unr+cc(i,j)))
            nd2        =nint(dd2l*dd2r)
            nss        =nd1+nd2
            ish(i,j-2) =ish(i,j-2)+nss
            ish(i,j-1) =ish(i,j-1)+nss
            ish(i,j  ) =ish(i,j  )+nss
            ish(i,j+1) =ish(i,j+1)+nss
         end do
      end do
c* in-tangential direction
c     do j=1,jmax-1
c        do i=2,imax-2
c           siv        =1.0/si(i,j,3)
c           unl        =siv*(si(i,j,1)*uc(i-1,j)+si(i,j,2)*vc(i-1,j))
c           unr        =siv*(si(i,j,1)*uc(i  ,j)+si(i,j,2)*vc(i  ,j))
c           dd1l       =max(0.0, sign(1.0,unl-cc(i-1,j)))
c           dd1r       =max(0.0,-sign(1.0,unr-cc(i  ,j)))
c           nd1        =nint(dd1l*dd1r)
c           dd2l       =max(0.0, sign(1.0,unl+cc(i-1,j)))
c           dd2r       =max(0.0,-sign(1.0,unr+cc(i  ,j)))
c           nd2        =nint(dd2l*dd2r)
c           nss        =nd1+nd2
c           ish(i-1,j) =ish(i-1,j)+nss
c           ish(i  ,j) =ish(i  ,j)+nss
c        end do
c     end do
c* termination
      return
      end
c
      subroutine flux(mode)
c     ************************************************************
c     *     AUSM-DV + shock fix                                  *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray3/ qcl(nx,ny,4),qcr(nx,ny,4),
     &               qsl(nx,ny,3),qsr(nx,ny,3),
     &               qvl(nx,ny,1),qvr(nx,ny,1)
      common /aray4/ dcq(nx,ny,4),dsq(nx,ny,3),dvq(0:nx,0:ny,1)
      common /aray5/ ish(0:nx,0:ny)
      common /coor1/ imax,jmax
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
c* (i)
      if(mode.eq.1) then
         do j=1,jmax-1
            do i=1,imax
               as        =si(i,j,3)
               asv       =1.0/as
               ax        =asv*si(i,j,1)
               ay        =asv*si(i,j,2)
c left state
               rl        =qcl(i,j,1)
               rlv       =1.0/rl
               ul        =qcl(i,j,2)
               vl        =qcl(i,j,3)
               pl        =qcl(i,j,4)
               unl       =ax*ul+ay*vl
               x1l       =rlv*qsl(i,j,1)
               x2l       =rlv*qsl(i,j,2)
               x3l       =rlv*qsl(i,j,3)
               qsl4      =bet1*qcl(i,j,1)+bet2*qsl(i,j,1)
     &                                   +bet3*qsl(i,j,3)
               qsl5      =bet4*qcl(i,j,1)+bet5*qsl(i,j,2)
     &                                   +bet6*qsl(i,j,3)
               ano       =amas1(1)*qsl(i,j,1)
               ann       =amas1(2)*qsl(i,j,2)
               anno      =amas1(3)*qsl(i,j,3)
               ano2      =amas1(4)*qsl4
               ann2      =amas1(5)*qsl5
               ansum     =ano+ann+anno+ano2+ann2
               tl        =pl/(8.314*ansum)
               cvav      =ano*cvi(1)+ann*cvi(2)+anno*cvi(3)+ano2*cvi(4)
     &                                                     +ann2*cvi(5)
               gml       =1.0+8.314*ansum/cvav
               cl        =sqrt(gml*pl*rlv)
               evw       =qvl(i,j,1)
               x6l       =rlv*evw
               eh        =h0i(1)*ano+h0i(2)*ann+h0i(3)*anno
               ek        =0.5*rl*(ul*ul+vl*vl)
               el        =cvav*tl+evw+eh+ek
               hl        =(el+pl)*rlv
c right state
               rr        =qcr(i,j,1)
               rrv       =1.0/rr
               ur        =qcr(i,j,2)
               vr        =qcr(i,j,3)
               pr        =qcr(i,j,4)
               unr       =ax*ur+ay*vr
               x1r       =rrv*qsr(i,j,1)
               x2r       =rrv*qsr(i,j,2)
               x3r       =rrv*qsr(i,j,3)
               qsr4      =bet1*qcr(i,j,1)+bet2*qsr(i,j,1)
     &                                   +bet3*qsr(i,j,3)
               qsr5      =bet4*qcr(i,j,1)+bet5*qsr(i,j,2)
     &                                   +bet6*qsr(i,j,3)
               ano       =amas1(1)*qsr(i,j,1)
               ann       =amas1(2)*qsr(i,j,2)
               anno      =amas1(3)*qsr(i,j,3)
               ano2      =amas1(4)*qsr4
               ann2      =amas1(5)*qsr5
               ansum     =ano+ann+anno+ano2+ann2
               tr        =pr/(8.314*ansum)
               cvav      =ano*cvi(1)+ann*cvi(2)+anno*cvi(3)+ano2*cvi(4)
     &                                                     +ann2*cvi(5)
               gmr       =1.0+8.314*ansum/cvav
!               cr2       =gmr*pr*rrv
               cr        =sqrt(gmr*pr*rrv)
               evw       =qvr(i,j,1)
               x6r       =rrv*evw
               eh        =h0i(1)*ano+h0i(2)*ann+h0i(3)*anno
               ek        =0.5*rr*(ur*ur+vr*vr)
               er        =cvav*tr+evw+eh+ek
               hr        =(er+pr)*rrv
c u_half(+-) and p_half (+-) using a common speed of sound on L/R state
            swvv      =(1.d0+min(1.d0,10.d0*abs(pr-pl)/min(pr,pl)))*0.5
c            swdps     =min(float(ish(i,j)+ish(i-1,j)),1.d0)
            swdps     =min(dble(ish(i,j)+ish(i-1,j)),1.d0)
               if(swdps.ne.0.0) swvv=1.0
c  sc(ro,p) = 1      iff strong shock -> Hanel 87
c           = p/ro   otherwise
               prav      =min(pl*rlv,pr*rrv)
               scl       =(1.0-swdps)*pl*rlv+swdps*prav
               scr       =(1.0-swdps)*pr*rrv+swdps*prav
               cmax      =max(cl,cr)
               cml       =(1.0-swdps)*cmax+swdps*cl
               cmr       =(1.0-swdps)*cmax+swdps*cr
               rml       =unl/cml
               rmr       =unr/cmr
               al        =2.0*scl/(scl+scr)
               ar        =2.0*scr/(scl+scr)
               if(abs(rml).le.1.0)then
                  plp    =pl*(rml+1.0)**2*(2.0-rml)*0.25
                  ulp0   =(unl+abs(unl))*0.5
                  ulp1   =(unl+cml)**2/(4.0*cml)
                  ulp    =al*(ulp1-ulp0)+ulp0
               elseif(rml.gt.1.0)then
                  plp    =pl
                  ulp    =unl
               else
                  plp    =0.0
                  ulp    =0.0
               endif
               if(abs(rmr).le.1.0)then
                  prm    = pr*(rmr-1.0)**2*(2.0+rmr)*0.25
                  urm0   = (unr-abs(unr))*0.5
                  urm1   =-(unr-cmr)**2/(4.0*cmr)
                  urm    = ar*(urm1-urm0)+urm0
               elseif(rmr.lt.-1.0)then
                  prm    =pr
                  urm    =unr
               else
                  prm    =0.0
                  urm    =0.0
               endif
c numerical flux
               rul       =ulp*rl
               rur       =urm*rr
               fm        =rul+rur
               f1d       =0.5*(fm*(x1r+x1l)-abs(fm)*(x1r-x1l))
               f1v       =rul*x1l+rur*x1r
               f1        =swdps*f1v+(1.0-swdps)*f1d
               f2d       =0.5*(fm*(x2r+x2l)-abs(fm)*(x2r-x2l))
               f2v       =rul*x2l+rur*x2r
               f2        =swdps*f2v+(1.0-swdps)*f2d
               f3d       =0.5*(fm*(x3r+x3l)-abs(fm)*(x3r-x3l))
               f3v       =rul*x3l+rur*x3r
               f3        =swdps*f3v+(1.0-swdps)*f3d
               f6d       =0.5*(fm*(x6r+x6l)-abs(fm)*(x6r-x6l))
               f6v       =rul*x6l+rur*x6r
               f6        =swdps*f6v+(1.0-swdps)*f6d
               f7d       =0.5*(fm*(unr+unl)-abs(fm)*(unr-unl))
               f7v       =rul*unl+rur*unr
               f7        =swvv*f7v+(1.0-swvv)*f7d+(plp+prm)
               f9d       =0.5*(fm*(hr+hl)-abs(fm)*(hr-hl))
               f9v       =rul*hl+rur*hr
               f9        =swdps*f9v+(1.0-swdps)*f9d
c tangential velocity components
               ult       =ul-ax*unl
               vlt       =vl-ay*unl
               urt       =ur-ax*unr
               vrt       =vr-ay*unr
c numerical fluxes
               qsl(i,j,1) =as*f1
               qsl(i,j,2) =as*f2
               qsl(i,j,3) =as*f3
               qvl(i,j,1) =as*f6
               qcl(i,j,1) =as*fm
               qcl(i,j,2) =as*(ax*f7+swdps*(rul*ult+rur*urt)
     &                            +(1.0-swdps)*0.5*(fm*(ult+urt)
     &                                        -abs(fm)*(urt-ult)))
               qcl(i,j,3) =as*(ay*f7+swdps*(rul*vlt+rur*vrt)
     &                            +(1.0-swdps)*0.5*(fm*(vlt+vrt)
     &                                        -abs(fm)*(vrt-vlt)))
               qcl(i,j,4) =as*f9
            end do
         end do
c increment
         do i=1,imax-1
            do j=1,jmax-1
               dsq(i,j,1) =-qsl(i+1,j,1)+qsl(i,j,1)
               dsq(i,j,2) =-qsl(i+1,j,2)+qsl(i,j,2)
               dsq(i,j,3) =-qsl(i+1,j,3)+qsl(i,j,3)
               dvq(i,j,1) =-qvl(i+1,j,1)+qvl(i,j,1)
               dcq(i,j,1) =-qcl(i+1,j,1)+qcl(i,j,1)
               dcq(i,j,2) =-qcl(i+1,j,2)+qcl(i,j,2)
               dcq(i,j,3) =-qcl(i+1,j,3)+qcl(i,j,3)
               dcq(i,j,4) =-qcl(i+1,j,4)+qcl(i,j,4)
            end do
         end do
c* interface (in-normal)
      elseif(mode.eq.2) then
         do j=1,jmax
            do i=1,imax-1
               as        =sj(i,j,3)
               asv       =1.0/as
               ax        =asv*sj(i,j,1)
               ay        =asv*sj(i,j,2)
c left state
               rl        =qcl(i,j,1)
               rlv       =1.0/rl
               ul        =qcl(i,j,2)
               vl        =qcl(i,j,3)
               pl        =qcl(i,j,4)
               unl       =ax*ul+ay*vl
               x1l       =rlv*qsl(i,j,1)
               x2l       =rlv*qsl(i,j,2)
               x3l       =rlv*qsl(i,j,3)
               qsl4      =bet1*qcl(i,j,1)+bet2*qsl(i,j,1)
     &                                   +bet3*qsl(i,j,3)
               qsl5      =bet4*qcl(i,j,1)+bet5*qsl(i,j,2)
     &                                   +bet6*qsl(i,j,3)
               ano       =amas1(1)*qsl(i,j,1)
               ann       =amas1(2)*qsl(i,j,2)
               anno      =amas1(3)*qsl(i,j,3)
               ano2      =amas1(4)*qsl4
               ann2      =amas1(5)*qsl5
               ansum     =ano+ann+anno+ano2+ann2
               tl        =pl/(8.314*ansum)
               cvav      =ano*cvi(1)+ann*cvi(2)+anno*cvi(3)+ano2*cvi(4)
     &                                                     +ann2*cvi(5)
               gml       =1.0+8.314*ansum/cvav
               cl        =sqrt(gml*pl*rlv)
               evw       =qvl(i,j,1)
               x6l       =rlv*evw
               eh        =h0i(1)*ano+h0i(2)*ann+h0i(3)*anno
               ek        =0.5*rl*(ul*ul+vl*vl)
               el        =cvav*tl+evw+eh+ek
               hl        =(el+pl)*rlv
c right state
               rr        =qcr(i,j,1)
               rrv       =1.0/rr
               ur        =qcr(i,j,2)
               vr        =qcr(i,j,3)
               pr        =qcr(i,j,4)
               unr       =ax*ur+ay*vr
               x1r       =rrv*qsr(i,j,1)
               x2r       =rrv*qsr(i,j,2)
               x3r       =rrv*qsr(i,j,3)
               qsr4      =bet1*qcr(i,j,1)+bet2*qsr(i,j,1)
     &                                   +bet3*qsr(i,j,3)
               qsr5      =bet4*qcr(i,j,1)+bet5*qsr(i,j,2)
     &                                   +bet6*qsr(i,j,3)
               ano       =amas1(1)*qsr(i,j,1)
               ann       =amas1(2)*qsr(i,j,2)
               anno      =amas1(3)*qsr(i,j,3)
               ano2      =amas1(4)*qsr4
               ann2      =amas1(5)*qsr5
               ansum     =ano+ann+anno+ano2+ann2
               tr        =pr/(8.314*ansum)
               cvav      =ano*cvi(1)+ann*cvi(2)+anno*cvi(3)+ano2*cvi(4)
     &                                                     +ann2*cvi(5)
               gmr       =1.0+8.314*ansum/cvav
               cr        =sqrt(gmr*pr*rrv)
               evw       =qvr(i,j,1)
               x6r       =rrv*evw
               eh        =h0i(1)*ano+h0i(2)*ann+h0i(3)*anno
               ek        =0.5*rr*(ur*ur+vr*vr)
               er        =cvav*tr+evw+eh+ek
               hr        =(er+pr)*rrv
c u_half(+-) and p_half(+-) using a common speed of sound on L/R state
            swvv      =(1.d0+min(1.d0,10.d0*abs(pr-pl)/min(pr,pl)))*0.5
c            swdps     =min(float(ish(i,j)+ish(i,j-1)),1.d0)
            swdps     =min(dble(ish(i,j)+ish(i,j-1)),1.d0)
               if(swdps.ne.0.0) swvv=1.0
c  sc(ro,p) = 1      iff strong shock -> Hanel 87
c           = p/ro   otherwise
               prav      =min(pl*rlv,pr*rrv)
               scl       =(1.0-swdps)*pl*rlv+swdps*prav
               scr       =(1.0-swdps)*pr*rrv+swdps*prav
               cmax      =max(cl,cr)
               cml       =(1.0-swdps)*cmax+swdps*cl
               cmr       =(1.0-swdps)*cmax+swdps*cr
               rml       =unl/cml
               rmr       =unr/cmr
               al        =2.0*scl/(scl+scr)
               ar        =2.0*scr/(scl+scr)
               if(abs(rml).le.1.0)then
                  plp    =pl*(rml+1.0)**2*(2.0-rml)*0.25
                  ulp0   =(unl+abs(unl))*0.5
                  ulp1   =(unl+cml)**2/(4.0*cml)
                  ulp    =al*(ulp1-ulp0)+ulp0
               elseif(rml.gt.1.0)then
                  plp    =pl
                  ulp    =unl
               else
                  plp    =0.0
                  ulp    =0.0
               endif
               if(abs(rmr).le.1.0)then
                  prm    = pr*(rmr-1.0)**2*(2.0+rmr)*0.25
                  urm0   = (unr-abs(unr))*0.5
                  urm1   =-(unr-cmr)**2/(4.0*cmr)
                  urm    = ar*(urm1-urm0)+urm0
               elseif(rmr.lt.-1.0)then
                  prm    =pr
                  urm    =unr
               else
                  prm    =0.0
                  urm    =0.0
               endif
c numerical flux
               rul       =ulp*rl
               rur       =urm*rr
               fm        =rul+rur
               f1d       =0.5*(fm*(x1r+x1l)-abs(fm)*(x1r-x1l))
               f1v       =rul*x1l+rur*x1r
               f1        =swdps*f1v+(1.0-swdps)*f1d
               f2d       =0.5*(fm*(x2r+x2l)-abs(fm)*(x2r-x2l))
               f2v       =rul*x2l+rur*x2r
               f2        =swdps*f2v+(1.0-swdps)*f2d
               f3d       =0.5*(fm*(x3r+x3l)-abs(fm)*(x3r-x3l))
               f3v       =rul*x3l+rur*x3r
               f3        =swdps*f3v+(1.0-swdps)*f3d
               f6d       =0.5*(fm*(x6r+x6l)-abs(fm)*(x6r-x6l))
               f6v       =rul*x6l+rur*x6r
               f6        =swdps*f6v+(1.0-swdps)*f6d
               f7d       =0.5*(fm*(unr+unl)-abs(fm)*(unr-unl))
               f7v       =rul*unl+rur*unr
               f7        =swvv*f7v+(1.0-swvv)*f7d+(plp+prm)
               f9d       =0.5*(fm*(hr+hl)-abs(fm)*(hr-hl))
               f9v       =rul*hl+rur*hr
               f9        =swdps*f9v+(1.0-swdps)*f9d
c tangential velocity components
               ult       =ul-ax*unl
               vlt       =vl-ay*unl
               urt       =ur-ax*unr
               vrt       =vr-ay*unr
c numerical fluxes
               qsl(i,j,1) =as*f1
               qsl(i,j,2) =as*f2
               qsl(i,j,3) =as*f3
               qvl(i,j,1) =as*f6
               qcl(i,j,1) =as*fm
               qcl(i,j,2) =as*(ax*f7+swdps*(rul*ult+rur*urt)
     &                              +(1.0-swdps)*0.5*(fm*(ult+urt)
     &                                          -abs(fm)*(urt-ult)))
               qcl(i,j,3) =as*(ay*f7+swdps*(rul*vlt+rur*vrt)
     &                              +(1.0-swdps)*0.5*(fm*(vlt+vrt)
     &                                          -abs(fm)*(vrt-vlt)))
               qcl(i,j,4) =as*f9
            end do
         end do
c increment
         do j=1,jmax-1
            do i=1,imax-1
               dsq(i,j,1) =dsq(i,j,1)-qsl(i,j+1,1)+qsl(i,j,1)
               dsq(i,j,2) =dsq(i,j,2)-qsl(i,j+1,2)+qsl(i,j,2)
               dsq(i,j,3) =dsq(i,j,3)-qsl(i,j+1,3)+qsl(i,j,3)
               dvq(i,j,1) =dvq(i,j,1)-qvl(i,j+1,1)+qvl(i,j,1)
               dcq(i,j,1) =dcq(i,j,1)-qcl(i,j+1,1)+qcl(i,j,1)
               dcq(i,j,2) =dcq(i,j,2)-qcl(i,j+1,2)+qcl(i,j,2)
               dcq(i,j,3) =dcq(i,j,3)-qcl(i,j+1,3)+qcl(i,j,3)
               dcq(i,j,4) =dcq(i,j,4)-qcl(i,j+1,4)+qcl(i,j,4)
            end do
         end do
      endif
c* termination
      return
      end
c
      subroutine vflx(mode)
c     *********************************************************
c     *     viscous flux function                             *
c     *********************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /aray3/ qcl(nx,ny,4),qcr(nx,ny,4),
     &               qsl(nx,ny,3),qsr(nx,ny,3),
     &               qvl(nx,ny,1),qvr(nx,ny,1)
      common /aray4/ dcq(nx,ny,4),dsq(nx,ny,3),dvq(0:nx,0:ny,1)
      common /aray7/ scs(nx,ny,3),scv(nx,ny,1)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /cntl5/ nprc,iprc
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
c* entry
      if(nprc.eq.0) return
c* (i)
      if(mode.eq.1) then
         do i=1,imax
            do j=1,jmax-1
               vmv       =2.0/(vol(i,j)+vol(i-1,j))
               dsx       =vmv*si(i,j,1)
               dsy       =vmv*si(i,j,2)
c derivatives
               drdx      =dsx*(qc(i,j,1)-qc(i-1,j,1))
               drdy      =dsy*(qc(i,j,1)-qc(i-1,j,1))
               dudx      =dsx*(uc(i,j)-uc(i-1,j))
               dudy      =dsy*(uc(i,j)-uc(i-1,j))
               dvdx      =dsx*(vc(i,j)-vc(i-1,j))
               dvdy      =dsy*(vc(i,j)-vc(i-1,j))
               dttdx     =dsx*(tt(i,j)-tt(i-1,j))
               dttdy     =dsy*(tt(i,j)-tt(i-1,j))
               dtvdx     =dsx*(tv(i,j)-tv(i-1,j))
               dtvdy     =dsy*(tv(i,j)-tv(i-1,j))
               dx1dx     =dsx*(xc(i,j,1)-xc(i-1,j,1))
               dx1dy     =dsy*(xc(i,j,1)-xc(i-1,j,1))
               dx2dx     =dsx*(xc(i,j,2)-xc(i-1,j,2))
               dx2dy     =dsy*(xc(i,j,2)-xc(i-1,j,2))
               dx3dx     =dsx*(xc(i,j,3)-xc(i-1,j,3))
               dx3dy     =dsy*(xc(i,j,3)-xc(i-1,j,3))
               dx4dx     =dsx*(xc(i,j,4)-xc(i-1,j,4))
               dx4dy     =dsy*(xc(i,j,4)-xc(i-1,j,4))
               dx5dx     =dsx*(xc(i,j,5)-xc(i-1,j,5))
               dx5dy     =dsy*(xc(i,j,5)-xc(i-1,j,5))
c interface values
               rf        =0.5*(qc(i,j,1)+qc(i-1,j,1))
               uf        =0.5*(uc(i,j)+uc(i-1,j))
               vf        =0.5*(vc(i,j)+vc(i-1,j))
               vmuf      =0.5*(vmue(i,j)+vmue(i-1,j))
               vkaf      =0.5*(vkap(i,j)+vkap(i-1,j))
               vkvf      =0.5*(vkav(i,j)+vkav(i-1,j))
               ds1f      =0.5*(ds(i,j,1)+ds(i-1,j,1))
               ds2f      =0.5*(ds(i,j,2)+ds(i-1,j,2))
               ds3f      =0.5*(ds(i,j,3)+ds(i-1,j,3))
               ds4f      =0.5*(ds(i,j,4)+ds(i-1,j,4))
               ds5f      =0.5*(ds(i,j,5)+ds(i-1,j,5))
               hs1f      =0.5*(hs(i,j,1)+hs(i-1,j,1))
               hs2f      =0.5*(hs(i,j,2)+hs(i-1,j,2))
               hs3f      =0.5*(hs(i,j,3)+hs(i-1,j,3))
               hs4f      =0.5*(hs(i,j,4)+hs(i-1,j,4))
               hs5f      =0.5*(hs(i,j,5)+hs(i-1,j,5))
               es1f      =0.5*(evs(i,j,1)+evs(i-1,j,1))
               es2f      =0.5*(evs(i,j,2)+evs(i-1,j,2))
               es3f      =0.5*(evs(i,j,3)+evs(i-1,j,3))
               es4f      =0.5*(evs(i,j,4)+evs(i-1,j,4))
               es5f      =0.5*(evs(i,j,5)+evs(i-1,j,5))
               if(i.eq.1) then
                  vfyv   =dvdy
! i = imax の時は微分ができないから修正
               else if(i.eq.imax) then
                  vfyv   =dvdy
               else
                  vfyv   =2.0*vf/(y(i,j)+y(i,j+1))
               endif
               txx       =0.6666667*vmuf*(2.0*dudx-dvdy-vfyv)
               tyy       =0.6666667*vmuf*(2.0*dvdy-dudx-vfyv)
               txy       =vmuf*(dudy+dvdx)
               r1u       =(rf*ds1f)*dx1dx
               r1v       =(rf*ds1f)*dx1dy
               r2u       =(rf*ds2f)*dx2dx
               r2v       =(rf*ds2f)*dx2dy
               r3u       =(rf*ds3f)*dx3dx
               r3v       =(rf*ds3f)*dx3dy
               r4u       =(rf*ds4f)*dx4dx
               r4v       =(rf*ds4f)*dx4dy
               r5u       =(rf*ds5f)*dx5dx
               r5v       =(rf*ds5f)*dx5dy
               qhx       =r1u*hs1f+r2u*hs2f+r3u*hs3f+r4u*hs4f+r5u*hs5f
               qhy       =r1v*hs1f+r2v*hs2f+r3v*hs3f+r4v*hs4f+r5v*hs5f
               qex       =r1u*es1f+r2u*es2f+r3u*es3f+r4u*es4f+r5u*es5f
               qey       =r1v*es1f+r2v*es2f+r3v*es3f+r4v*es4f+r5v*es5f
               dqtdx     =vkaf*dttdx
               dqtdy     =vkaf*dttdy
               dqvdx     =vkvf*dtvdx
               dqvdy     =vkvf*dtvdy
               dqdx      =dqtdx+dqvdx+qhx
               dqdy      =dqtdy+dqvdy+qhy
c viscous flux function
               qsl(i,j,1) =si(i,j,1)*r1u+si(i,j,2)*r1v
               qsl(i,j,2) =si(i,j,1)*r2u+si(i,j,2)*r2v
               qsl(i,j,3) =si(i,j,1)*r3u+si(i,j,2)*r3v
               qvl(i,j,1) =si(i,j,1)*(dqvdx+qex)+si(i,j,2)*(dqvdy+qey)
               qcl(i,j,2) =si(i,j,1)*txx+si(i,j,2)*txy
               qcl(i,j,3) =si(i,j,1)*txy+si(i,j,2)*tyy
               qcl(i,j,4) =uf*qcl(i,j,2)+vf*qcl(i,j,3)+si(i,j,1)*dqdx
     &                                                +si(i,j,2)*dqdy
            end do
         end do
c increment
         do i=1,imax-1
            do j=1,jmax-1
               dsq(i,j,1) =dsq(i,j,1)+qsl(i+1,j,1)-qsl(i,j,1)
               dsq(i,j,2) =dsq(i,j,2)+qsl(i+1,j,2)-qsl(i,j,2)
               dsq(i,j,3) =dsq(i,j,3)+qsl(i+1,j,3)-qsl(i,j,3)
               dvq(i,j,1) =dvq(i,j,1)+qvl(i+1,j,1)-qvl(i,j,1)
               dcq(i,j,2) =dcq(i,j,2)+qcl(i+1,j,2)-qcl(i,j,2)
               dcq(i,j,3) =dcq(i,j,3)+qcl(i+1,j,3)-qcl(i,j,3)
               dcq(i,j,4) =dcq(i,j,4)+qcl(i+1,j,4)-qcl(i,j,4)
            end do
         end do
c* in-normal
      elseif(mode.eq.2) then
         do j=1,jmax
            do i=1,imax-1
               vmv       =2.0/(vol(i,j)+vol(i,j-1))
               dsx       =vmv*sj(i,j,1)
               dsy       =vmv*sj(i,j,2)
c derivatives
               drdx      =dsx*(qc(i,j,1)-qc(i,j-1,1))
               drdy      =dsy*(qc(i,j,1)-qc(i,j-1,1))
               dudx      =dsx*(uc(i,j)-uc(i,j-1))
               dudy      =dsy*(uc(i,j)-uc(i,j-1))
               dvdx      =dsx*(vc(i,j)-vc(i,j-1))
               dvdy      =dsy*(vc(i,j)-vc(i,j-1))
               dttdx     =dsx*(tt(i,j)-tt(i,j-1))
               dttdy     =dsy*(tt(i,j)-tt(i,j-1))
               dtvdx     =dsx*(tv(i,j)-tv(i,j-1))
               dtvdy     =dsy*(tv(i,j)-tv(i,j-1))
               dx1dx     =dsx*(xc(i,j,1)-xc(i,j-1,1))
               dx1dy     =dsy*(xc(i,j,1)-xc(i,j-1,1))
               dx2dx     =dsx*(xc(i,j,2)-xc(i,j-1,2))
               dx2dy     =dsy*(xc(i,j,2)-xc(i,j-1,2))
               dx3dx     =dsx*(xc(i,j,3)-xc(i,j-1,3))
               dx3dy     =dsy*(xc(i,j,3)-xc(i,j-1,3))
               dx4dx     =dsx*(xc(i,j,4)-xc(i,j-1,4))
               dx4dy     =dsy*(xc(i,j,4)-xc(i,j-1,4))
               dx5dx     =dsx*(xc(i,j,5)-xc(i,j-1,5))
               dx5dy     =dsy*(xc(i,j,5)-xc(i,j-1,5))
c interface values
               rf        =0.5*(qc(i,j,1)+qc(i,j-1,1))
               uf        =0.5*(uc(i,j)+uc(i,j-1))
               vf        =0.5*(vc(i,j)+vc(i,j-1))
               vmuf      =0.5*(vmue(i,j)+vmue(i,j-1))
               vkaf      =0.5*(vkap(i,j)+vkap(i,j-1))
               vkvf      =0.5*(vkav(i,j)+vkav(i,j-1))
               ds1f      =0.5*(ds(i,j,1)+ds(i,j-1,1))
               ds2f      =0.5*(ds(i,j,2)+ds(i,j-1,2))
               ds3f      =0.5*(ds(i,j,3)+ds(i,j-1,3))
               ds4f      =0.5*(ds(i,j,4)+ds(i,j-1,4))
               ds5f      =0.5*(ds(i,j,5)+ds(i,j-1,5))
               hs1f      =0.5*(hs(i,j,1)+hs(i,j-1,1))
               hs2f      =0.5*(hs(i,j,2)+hs(i,j-1,2))
               hs3f      =0.5*(hs(i,j,3)+hs(i,j-1,3))
               hs4f      =0.5*(hs(i,j,4)+hs(i,j-1,4))
               hs5f      =0.5*(hs(i,j,5)+hs(i,j-1,5))
               es1f      =0.5*(evs(i,j,1)+evs(i,j-1,1))
               es2f      =0.5*(evs(i,j,2)+evs(i,j-1,2))
               es3f      =0.5*(evs(i,j,3)+evs(i,j-1,3))
               es4f      =0.5*(evs(i,j,4)+evs(i,j-1,4))
               es5f      =0.5*(evs(i,j,5)+evs(i,j-1,5))
               vfyv      =2.0*vf/(y(i,j)+y(i+1,j))
               txx       =0.6666667*vmuf*(2.0*dudx-dvdy-vfyv)
               tyy       =0.6666667*vmuf*(2.0*dvdy-dudx-vfyv)
               txy       =vmuf*(dudy+dvdx)
               r1u       =(rf*ds1f)*dx1dx
               r1v       =(rf*ds1f)*dx1dy
               r2u       =(rf*ds2f)*dx2dx
               r2v       =(rf*ds2f)*dx2dy
               r3u       =(rf*ds3f)*dx3dx
               r3v       =(rf*ds3f)*dx3dy
               r4u       =(rf*ds4f)*dx4dx
               r4v       =(rf*ds4f)*dx4dy
               r5u       =(rf*ds5f)*dx5dx
               r5v       =(rf*ds5f)*dx5dy
               qhx       =r1u*hs1f+r2u*hs2f+r3u*hs3f+r4u*hs4f+r5u*hs5f
               qhy       =r1v*hs1f+r2v*hs2f+r3v*hs3f+r4v*hs4f+r5v*hs5f
               qex       =r1u*es1f+r2u*es2f+r3u*es3f+r4u*es4f+r5u*es5f
               qey       =r1v*es1f+r2v*es2f+r3v*es3f+r4v*es4f+r5v*es5f
               dqtdx     =vkaf*dttdx
               dqtdy     =vkaf*dttdy
               dqvdx     =vkvf*dtvdx
               dqvdy     =vkvf*dtvdy
               dqdx      =dqtdx+dqvdx+qhx
               dqdy      =dqtdy+dqvdy+qhy
c viscous flux function
               qsl(i,j,1) =sj(i,j,1)*r1u+sj(i,j,2)*r1v
               qsl(i,j,2) =sj(i,j,1)*r2u+sj(i,j,2)*r2v
               qsl(i,j,3) =sj(i,j,1)*r3u+sj(i,j,2)*r3v
               qvl(i,j,1) =sj(i,j,1)*(dqvdx+qex)+sj(i,j,2)*(dqvdy+qey)
               qcl(i,j,2) =sj(i,j,1)*txx+sj(i,j,2)*txy
               qcl(i,j,3) =sj(i,j,1)*txy+sj(i,j,2)*tyy
               qcl(i,j,4) =uf*qcl(i,j,2)+vf*qcl(i,j,3)+sj(i,j,1)*dqdx
     &                                                +sj(i,j,2)*dqdy
            end do
         end do
c increment
         do j=1,jmax-1
            do i=1,imax-1
               dsq(i,j,1) =dsq(i,j,1)+qsl(i,j+1,1)-qsl(i,j,1)
               dsq(i,j,2) =dsq(i,j,2)+qsl(i,j+1,2)-qsl(i,j,2)
               dsq(i,j,3) =dsq(i,j,3)+qsl(i,j+1,3)-qsl(i,j,3)
               dvq(i,j,1) =dvq(i,j,1)+qvl(i,j+1,1)-qvl(i,j,1)
               dcq(i,j,2) =dcq(i,j,2)+qcl(i,j+1,2)-qcl(i,j,2)
               dcq(i,j,3) =dcq(i,j,3)+qcl(i,j+1,3)-qcl(i,j,3)
               dcq(i,j,4) =dcq(i,j,4)+qcl(i,j+1,4)-qcl(i,j,4)
            end do
         end do
      endif
c* termination
      return
      end
c
      subroutine saxi
c     *********************************************************
c     *     axisymmteric term                                 *
c     *********************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /aray3/ qcl(nx,ny,4),qcr(nx,ny,4),
     &               qsl(nx,ny,3),qsr(nx,ny,3),
     &               qvl(nx,ny,1),qvr(nx,ny,1)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /coor3/ xg(0:nx,0:ny),yg(0:nx,0:ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /cntl5/ nprc,iprc
      common /cntl6/ nacu
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /wbcd1/ twall,cefn,cefo
      common /fscd1/ fsv,alpha
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /fscd3/ vmueinf,vkapinf,vkavinf,
     &               ds1inf,ds2inf,ds3inf,ds4inf,ds5inf,
     &               evs1inf,evs2inf,evs3inf,evs4inf,evs5inf,
     &               hs1inf,hs2inf,hs3inf,hs4inf,hs5inf
c* axisymmetric term-1 (inviscid term)
      do j=1,jmax-1
         do i=1,imax-1
            ygvc       =vc(i,j)/yg(i,j)
            qsr(i,j,1) =-ygvc*qs(i,j,1)
            qsr(i,j,2) =-ygvc*qs(i,j,2)
            qsr(i,j,3) =-ygvc*qs(i,j,3)
            qvr(i,j,1) =-ygvc*qv(i,j,1)
            qcr(i,j,1) =-ygvc*qc(i,j,1)
            qcr(i,j,2) =uc(i,j)*qcr(i,j,1)
            qcr(i,j,3) =vc(i,j)*qcr(i,j,1)
            qcr(i,j,4) =-ygvc*(qc(i,j,4)+pc(i,j))
         end do
      end do
      if(nprc.eq.0) return
c* axisymmetric term-2 (diffusion and viscous terms)
      do j=1,jmax-1
         do i=1,imax-1
            volv       =1.0/vol(i,j)
            sxx        =(0.25*volv)*(si(i,j,1)+si(i+1,j,1))
            sxy        =(0.25*volv)*(si(i,j,2)+si(i+1,j,2))
            sex        =(0.25*volv)*(sj(i,j,1)+sj(i,j+1,1))
            sey        =(0.25*volv)*(sj(i,j,2)+sj(i,j+1,2))
c            sex        =(0.25*volv)*(sj(i,j,1)+sj(i+1,j,1))
c            sey        =(0.25*volv)*(sj(i,j,2)+sj(i+1,j,2))
            dudx       =sxx*(uc(i+1,j  )-uc(i-1,j  ))
     &                 +sex*(uc(i  ,j+1)-uc(i  ,j-1))
            dudy       =sxy*(uc(i+1,j  )-uc(i-1,j  ))
     &                 +sey*(uc(i  ,j+1)-uc(i  ,j-1))
            dvdx       =sxx*(vc(i+1,j  )-vc(i-1,j  ))
     &                 +sex*(vc(i  ,j+1)-vc(i  ,j-1))
            dvdy       =sxy*(vc(i+1,j  )-vc(i-1,j  ))
     &                 +sey*(vc(i  ,j+1)-vc(i  ,j-1))
            dx1dy      =sxy*(xc(i+1,j  ,1)-xc(i-1,j  ,1))
     &                 +sey*(xc(i  ,j+1,1)-xc(i  ,j-1,1))
            dx2dy      =sxy*(xc(i+1,j  ,2)-xc(i-1,j  ,2))
     &                 +sey*(xc(i  ,j+1,2)-xc(i  ,j-1,2))
            dx3dy      =sxy*(xc(i+1,j  ,3)-xc(i-1,j  ,3))
     &                 +sey*(xc(i  ,j+1,3)-xc(i  ,j-1,3))
            dx4dy      =sxy*(xc(i+1,j  ,4)-xc(i-1,j  ,4))
     &                 +sey*(xc(i  ,j+1,4)-xc(i  ,j-1,4))
            dx5dy      =sxy*(xc(i+1,j  ,5)-xc(i-1,j  ,5))
     &                 +sey*(xc(i  ,j+1,5)-xc(i  ,j-1,5))
            dttdy      =sxy*(tt(i+1,j  )-tt(i-1,j  ))
     &                 +sey*(tt(i  ,j+1)-tt(i  ,j-1))
            dtvdy      =sxy*(tv(i+1,j  )-tv(i-1,j  ))
     &                 +sey*(tv(i  ,j+1)-tv(i  ,j-1))
            ygv        =1.0/yg(i,j)
            rygv       =ygv*qc(i,j,1)
            txy        =vmue(i,j)*(dudy+dvdx)
            tyy        =0.6667*vmue(i,j)*(2.0*dvdy-dudx-ygv*vc(i,j))
            thh        =0.6667*vmue(i,j)*(2.0*ygv*vc(i,j)-dudx-dvdy)
            ds1cy      =ds(i,j,1)*dx1dy
            ds2cy      =ds(i,j,2)*dx2dy
            ds3cy      =ds(i,j,3)*dx3dy
            ds4cy      =ds(i,j,4)*dx4dy
            ds5cy      =ds(i,j,5)*dx5dy
            qsr(i,j,1) =qsr(i,j,1)+rygv*ds1cy
            qsr(i,j,2) =qsr(i,j,2)+rygv*ds2cy
            qsr(i,j,3) =qsr(i,j,3)+rygv*ds3cy
            qcr(i,j,2) =qcr(i,j,2)+ygv*txy
            qcr(i,j,3) =qcr(i,j,3)+ygv*(tyy-thh)
            qcr(i,j,4) =qcr(i,j,4)+ygv*(uc(i,j)*txy+vc(i,j)*tyy
     &                            +vkap(i,j)*dttdy+vkav(i,j)*dtvdy)
     &                      +rygv*(ds1cy*hs(i,j,1)+ds2cy*hs(i,j,2)
     &                            +ds3cy*hs(i,j,3)+ds4cy*hs(i,j,4)
     &                            +ds5cy*hs(i,j,5))
            qvr(i,j,1) =qvr(i,j,1)+ygv*vkav(i,j)*dtvdy
     &                      +rygv*(ds1cy*evs(i,j,1)+ds2cy*evs(i,j,2)
     &                            +ds3cy*evs(i,j,3)+ds4cy*evs(i,j,4)
     &                            +ds5cy*evs(i,j,5))
         end do
      end do
c* termination
      return
      end
c
      subroutine tsol
c     **************************************************************
c     *     calculates vibrational-electronic temperature tv and   *
c     *     translational-rotational temperature tt                *
c     *      inputs: qc=conservative variables                     *
c     *      outputs: tt=translational-rotational temperature [K]  *
c     *               tv=vibrational-electroic temperature [K]     *
c     *      uses 5-species, 2-temperature model                   *
c     *      species index: 1=o,2=n,3=no,4=o2,5=n2                 *
c     *      accuracy verified against numerical jacobians,        *
c     *               april 19, 1988, chul park                    *
c     **************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /parm3/ mdd
      common /stat1/ loop,locl
c* set primitive variables
      ttmax =-1.0d30
      tvmax =-1.0d30
      ttmin = 1.0d30
      tvmin = 1.0d30
      do j=1,jmax-1
         do i=1,imax-1
            rv        =1.0/qc(i,j,1)
            uc(i,j)   =rv*qc(i,j,2)
            vc(i,j)   =rv*qc(i,j,3)
            qs4       =bet1*qc(i,j,1)+bet2*qs(i,j,1)+bet3*qs(i,j,3)
            qs5       =bet4*qc(i,j,1)+bet5*qs(i,j,2)+bet6*qs(i,j,3)
            xc(i,j,1) =rv*qs(i,j,1)
            xc(i,j,2) =rv*qs(i,j,2)
            xc(i,j,3) =rv*qs(i,j,3)
            xc(i,j,4) =rv*qs4
            xc(i,j,5) =rv*qs5
c* ano, ann, anno, ano2, and ann2 are number densities in mole/m3
            rhou      =qc(i,j,2)
            rhov      =qc(i,j,3)
            e         =qc(i,j,4)
            ano       =qs(i,j,1)*amas1(1)
            ann       =qs(i,j,2)*amas1(2)
            anno      =qs(i,j,3)*amas1(3)
            ano2      =qs4      *amas1(4)
            ann2      =qs5      *amas1(5)
            ev        =qv(i,j,1)
c* normalized vibrational-electronic energy evp=ev(j/m3)/(avogadro*k)
            evp       =ev/8.314
c* determines vibrational-electronic temperature tv
c      uses approximate energy levels (everything is in mks units)
c      inputs
c        evp=ev/(av*k), av=6.022e23, k=1.38e-23
c        tv=first estimate of tv
c      output
c        tv=calculated value of tv
c      use newton's iteration method.
c newton iteration
            do m=1,10 ! Newton 法の最大反復回数
               tvPre     =tv(i,j) ! 前回の Tv
               tv1       =1.0/tv(i,j)
               xa        =exp(226.0*tv1)
               xb        =exp(-3785.5*tv1)
               xa2       =xa*xa
               xa4       =xa2*xa2
               xa8       =xa4*xa4
               xa9       =xa8*xa
               xa10      =xa8*xa2
               xa11      =xa10*xa
               xa12      =xa10*xa2
               xa14      =xa12*xa2
               xa15      =xa12*xa2*xa
               xb2       =xb*xb
               xb3       =xb2*xb
               xb4       =xb2*xb2
               xb5       =xb4*xb
               xb6       =xb5*xb
               xb7       =xb6*xb
               xb11      =xb5*xb6
               xb13      =xb6*xb7
               ev1       =2260.0/(xa10-1.0)
               ev2       =3390.0/(xa15-1.0)
               ev3       =2712.0/(xa12-1.0)
               ev4       =(22713.0*xb3+18927.5*xb5)/(3.0+2.0*xb3+xb5)
               ev5       =113565.0*xb6/(9.0+5.0*xb6)
               ev6       =264985.0*xb7/(4.0+10.0*xb7)
               eva       =ano2*ev1+ann2*ev2+anno*ev3+ano2*ev4+ano*ev5
     &                                                       +ann*ev6
               dev1      =ano2*22600.0*xa9/((xa10-1.0)**2)
               dev2      =ann2*50850.0*xa14/((xa15-1.0)**2)
               dev3      =anno*32544.0*xa11/((xa12-1.0)**2)
               dev4      =(ano2/((3.0+2.0*xb3+xb5)**2))*((68139.0*xb2
     &                   +94637.5*xb4)*(3.0+2.0*xb3+xb5)
     &                   -(6.0*xb2+5.0*xb4)*(22713.0*xb3+18927.5*xb5))
               dev5      =ano*113565.0*(6.0*xb5*(9.0+5.0*xb6)
     &                                -30.0*xb11)/((9.0+5.0*xb6)**2)
               dev6      =ann*264985.0*(7.0*xb6*(4.0+10.0*xb7)
     &                                -70.0*xb13)/((4.0+10.0*xb7)**2)
               devdta    =tv1*tv1*(226.0*xa*(dev1+dev2+dev3)
     &                           +3785.5*xb*(dev4+dev5+dev6))
               tv(i,j)   =tv(i,j)+(evp-eva)/devdta
c              tv(i,j)   =dmax1(tv(i,j),0.5d0*tvinf)  !!!
c              tv(i,j)   =dmax1(tv(i,j),tvinf)  !!!
               tvRate    =dabs((tv(i,j)-tvPre)/tv(i,j)) ! 変化率
               if (tvRate.lt.1.0d-05) exit ! 収束半径 : 桁数指定
            end do
c* determine translational-rotational temperautre
            rv        =1.0/qc(i,j,1)
            cvt       =e-ev-0.5*(rhou**2+rhov**2)*rv
     &                -(h0i(1)*ano+h0i(2)*ann+h0i(3)*anno)
            cvav      =ano*cvi(1)+ann*cvi(2)+anno*cvi(3)+ano2*cvi(4)
     &                                                  +ann2*cvi(5)
            tt(i,j)   =cvt/cvav
            ! for Robust
            if(tt(i,j).lt.0.0d0) then
               tt(i,j)= 1.0d-30 ! restrict
            end if
c            tt(i,j)   =dmax1(tt(i,j),0.5d0*ttinf) !!!
c            tt(i,j)   =dmax1(tt(i,j),ttinf) !!!
c* set primitive varibales(II)
            ansum     =ano+ann+anno+ano2+ann2
            gm        =1.0+8.314*ansum/cvav
            pc(i,j)   =8.314*ansum*tt(i,j)
            cc(i,j)   =sqrt(gm*pc(i,j)*rv)
c
            if(tt(i,j).gt.ttmax) then
               ittmax= i
               jttmax= j
            end if
            if(tt(i,j).lt.ttmin) then
               ittmin= i
               jttmin= j
            end if
            if(tv(i,j).gt.tvmax) then
               itvmax= i
               jtvmax= j
            end if
            if(tv(i,j).lt.tvmin) then
               itvmin= i
               jtvmin= j
            end if
            ttmax     =dmax1(tt(i,j),ttmax)
            tvmax     =dmax1(tv(i,j),tvmax)
            ttmin     =dmin1(tt(i,j),ttmin)
            tvmin     =dmin1(tv(i,j),tvmin)
         end do
      end do
      ! Tv_min のステップ変化を出力する
      write(31,103) loop,tvmin,x(itvmin,jtvmin),y(itvmin,jtvmin)
103   format(i6,e11.3,e11.3,e11.3)
c*print control
      mw  =mod(loop,mdd)
      if(loop.le.100) mw =0
      if(mw.eq.0) then
         write(6,*) ''
         write(6,*) 'tsol) ttmax=',ttmax,' tvmax=',tvmax
         write(6,*) 'tsol) ttmin=',ttmin,' tvmin=',tvmin
         write(6,*) 'tsol) ittmax=',ittmax,' jttmax=',jttmax
         write(6,*) 'tsol) itvmax=',itvmax,' jtvmax=',jtvmax
         write(6,*) 'tsol) ittmin=',ittmin,' jttmin=',jttmin
         write(6,*) 'tsol) itvmin=',itvmin,' jtvmin=',jtvmin
      endif
c* termination
      return
      end
c
      subroutine htbs
c     **************************************************************
c     *     reference temperature behind shock wave                *
c     **************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /chem2/ tss(nx)
c* set
      do i=1,imax-1
         tss(i) =3991.2
         do j=1,jmax-1
            tss(i)=max(tss(i),tt(i,j))
         end do
      end do
c* termination
      return
      end
c
      subroutine sorc
c     ************************************************************
c     *     chemical source terms                                *
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /aray7/ scs(nx,ny,3),scv(nx,ny,1)
      common /aray8/ ajs(nx,ny,3),ajv(nx,ny,1)
      common /aray9/ tvev(nx,0:ny)
      common /coor1/ imax,jmax
      common /cntl10/nvis,nrlx
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /dtst1/ dt(nx,ny),dt0,dti
      common /cntl5/ nprc,iprc
      common /cntl8/ nchm
      common /cntl9/ qq,fmul3,fmul4
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem2/ tss(nx)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /chem6/ telim
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
c----------------------------------------------------------------------
c     >chemistry package to calculate source and its jacobian
c     >the original version in aiaa paper 89-0685
c     >modification made in feb 1989 to expand vibrational relaxation
c     >model to include all 5 colliding species, and to include the
c      thermal dissociation of n2
c
c inputs:
c     qc =conservative variables
c     tv =vibrational-electronic temperature, K.
c     tt =translational-rotational temperature, K.
c
c outputs:
c     sc =source
c
c     >uses 5-species, 4-reaction, 2-temperature model
c     >species index: 1=o,2=n,3=no,4=o2,5=n2
c     >reaction index: 1.o2+m=o+o+m, 2.n2+o=no+n, 3.no+o=o2+n,
c                      4.n2+m=n+n+m, 5.no+m=n+o+m
c     >accuracy verified against numerical jacobians, aug 3, 1989,
c      chul park
c----------------------------------------------------------------------
      dimension crat1(5),crat4(5),crat5(5),crate(3),
     &          crat1a(5),crat4a(5),crat5a(5),cratea(3)
      dimension aka1(5),aka2(5),aka3(5),aka4(5),aka5(5)
c----------------------------------------------------------------------
c     >equilibrium constants aka, are consistent with mole/kg unit,
c      taken from park, "nonequilibrium hypersonic aerothermodynamics",
c      at number density of 1.0e17 cm-3.
c     >for reaction 1,3 and 5, a2 is increased by 13.8155 for mole/m3
c      unit
c----------------------------------------------------------------------
      data aka1/0.553880, 16.275511,  1.776300,  -6.57200,  0.031445/
      data aka2/0.976460,  0.890430,  0.745720,  -3.96420,  0.007123/
      data aka3/0.004815, -1.744300, -1.222700,  -0.95824, -0.045545/
      data aka4/1.535100, 15.421600,  1.299300, -11.49400, -0.006980/
      data aka5/0.558890, 14.531080,  0.553960,  -7.53040, -0.014089/
c----------------------------------------------------------------------
c     >crat1a(i) is the o2 dissociation rate coef in m3/(mole-s) with
c      species i as the colliding partner
c     >cratea(2) and cratea(3) are the no exchange rate coef in
c      m3/(mole-s)
c     >crate4a(i) is the n2 dissociation rate coef in m3/(mole-s) with
c      species i as the colliding partner
c     >crate5a(i) is the no dissociation rate coef in m3/(mole-s) with
c      species i as the colliding partner
c----------------------------------------------------------------------
      data crat1a/2*1.0e16, 3*2.0e15/,
c     &     cratea/0., 6.4e11, 8.4e6/,
c modified by Ishihara
     &     cratea/0.0d0,5.7d6,8.4d6/, ! Bose and Candler
     &     crat4a/2*3.0e16, 3*7.0e15/,
     &     crat5a/3*1.1e11, 2*5.0e9/
c----------------------------------------------------------------------
c     >ebnoka,ebo2ka,ebn2ka are preferential removal energies in deg k
c      taken presently top be 100% of dissociation energy
c----------------------------------------------------------------------
c      data ebnoka/75500./,ebo2ka/59500./,ebn2ka/113200./
c modified by Ishihara
      data ebnoka/75500./,ebo2ka/59360./,ebn2ka/113200./
c----------------------------------------------------------------------
c     >aij, bij are a- and b-values of millikan and white for molecule
c      i in collision with species j. i=1,2,3:no,o2,n2.
c                                     j=1,2,3,4,5:o,n,no,o2,n2
c     >approximate aij, bij model used
c     >exbij is exp(bij)
c----------------------------------------------------------------------
c      data a11/40.0/, a12/120.0/, a13/ 40.0/, a14/120.0/, a15/120.0/,
c     &     a21/40.0/, a22/ 60.0/, a23/120.0/, a24/120.0/, a25/160.0/,
c     &     a31/60.0/, a32/160.0/, a33/160.0/, a34/160.0/, a35/220.0/
c      data exb11/1.319e9/, exb12/2.174e9/, exb13/4.852e8/,
c     &     exb14/2.174e9/, exb15/2.174e9/,
c     &     exb21/1.319e9/, exb22/1.785e8/, exb23/2.174e9/,
c     &     exb24/2.174e9/, exb25/9.745e9/,
c     &     exb31/1.785e8/, exb32/9.745e9/, exb33/9.745e9/,
c     &     exb34/9.745e9/, exb35/7.200e10/
c aij, bij is used from park, JTHT Vol.7, No.3, 1993
c modified by Ishihara
      data a11/49.5d0/,a12/49.5d0/,a13/ 49.5d0/,a14/49.5d0/,a15/49.5d0/,
     &     a21/47.7d0/,a22/72.4d0/,a23/136.d0/,a24/138.d0/,a25/134.d0/,
     &     a31/72.4d0/,a32/180.d0/,a33/225.d0/,a34/229.d0/,a35/221.d0/
c exbij is exp(aij*bij+18.42)
c modified by Ishihara
      data exb11/7.9910268d8/,exb12/7.9910268d8/,exb13/7.9910268d8/,
     &     exb14/7.9910268d8/,exb15/7.9910268d8/,
     &     exb21/1.6670143d9/,exb22/2.9603848d8/,exb23/5.7519226d9/,
     &     exb24/6.2760083d9/,exb25/5.2055958d9/,
     &     exb31/2.9603848d8/,exb32/1.1164445d10/,exb33/7.2910609d10/,
     &     exb34/8.5818462d10/,exb35/6.0687283d10/
c----------------------------------------------------------------------
c     >avmas is average molecular weight, assumed to be 22 g/mole
c     >avmas1 is its reciprocal
c----------------------------------------------------------------------
      data avmas/0.022/
      avmas1 =1.0/avmas
c* entry
      if(nchm.eq.0) then
         do j=1,jmax-1
           do i=1,imax-1
               scs(i,j,1) =0.0
               scs(i,j,2) =0.0
               scs(i,j,3) =0.0
               scv(i,j,1) =0.0
               ajs(i,j,1) =0.0
               ajs(i,j,2) =0.0
               ajs(i,j,3) =0.0
               ajv(i,j,1) =0.0
               tvev(i,j)  =0.0
            end do
         end do
         return
      endif
c ************************
c section 0. set constants
c ************************
      rmx   =1.0d110
c----------------------------------------------------------------------
c     >fmul3 is a multiplicative factor for collision limiting cross
c      section. the larger the fmul4, the faster is relaxation rate
c     >fmul4 is a dividing factor for the approximate millikan-white
c      p-tau in atm-sec. the larger the fmul4, the faster is
c      relaxation rate
c----------------------------------------------------------------------
      pd    =0.3
c      fmul3 =1.0
c      fmul4 =1.0
c
      ebnok =ebnoka*pd
      ebo2k =ebo2ka*pd
      ebn2k =ebn2ka*pd
      do k=1,5
         crat1(k) =fmul3*crat1a(k)
         crat4(k) =fmul3*crat4a(k)
         crat5(k) =fmul3*crat5a(k)
      end do
      crate(2) =fmul3*cratea(2)
      crate(3) =fmul3*cratea(3)
      g214     =2.14e-15/fmul4
      h85      =fmul4*8.2e-5
c *********************************************************
c section 1. convert computational variables into mks units
c *********************************************************
c ano, ann, anno, ano2, and ann2 are number densities in mole/m3
      scmn = 1.0e30
      scmx =-1.0e30
      do j=1,jmax-1
         do i=1,imax-1
            rhoo   =qs(i,j,1)
            rhon   =qs(i,j,2)
            rhono  =qs(i,j,3)
            rhoo2  =xc(i,j,4)*qc(i,j,1)
            rhon2  =xc(i,j,5)*qc(i,j,1)
            ev     =qv(i,j,1)
            rho    =qc(i,j,1)
            rhou   =qc(i,j,2)
            rhov   =qc(i,j,3)
            e      =qc(i,j,4)
            ano    = rhoo*amas1(1)
            ann    = rhon*amas1(2)
            anno   =rhono*amas1(3)
            ano2   =rhoo2*amas1(4)
            ann2   =rhon2*amas1(5)
            ansum  =ano+ann+anno+ano2+ann2
            ansum1 =1.0/ansum
c ***********************************************************
c section 2. derivatives of tv with respect to flow variables
c ***********************************************************
c derivative of ev with respect to tv dev/dtv
            tv1    =1.0/tv(i,j)
            xa     =exp(226.0*tv1)
            xb     =exp(-3785.5*tv1)
            xa2    =xa*xa
            xa4    =xa2*xa2
            xa8    =xa4*xa4
            xa9    =xa8*xa
            xa10   =xa8*xa2
            xa11   =xa10*xa
            xa12   =xa10*xa2
            xa14   =xa12*xa2
            xa15   =xa12*xa2*xa
            xb2    =xb*xb
            xb3    =xb2*xb
            xb4    =xb2*xb2
            xb5    =xb4*xb
            xb6    =xb5*xb
            xb7    =xb6*xb
            xb11   =xb5*xb6
            xb13   =xb6*xb7
            ev1    =2260.0/(xa10-1.0)
            ev2    =3390.0/(xa15-1.0)
            ev3    =2712.0/(xa12-1.0)
            ev4    =(22713.0*xb3+18927.5*xb5)/(3.0+2.0*xb3+xb5)
            ev5    =113565.0*xb6/(9.0 +5.0*xb6)
            ev6    =264985.0*xb7/(4.0+10.0*xb7)
            evb    =ano2*ev1+ann2*ev2+anno*ev3+ano2*ev4
     &             +ano*ev5 +ann*ev6
            dev1   =ano2*22600.0*xa9/((xa10-1.0)**2)
            dev2   =ann2*50850.0*xa14/((xa15-1.0)**2)
            dev3   =anno*32544.0*xa11/((xa12-1.0)**2)
            dev4p  =(1.0/((3.0+2.0*xb3+xb5)**2))*((68139.0*xb2
     &             +94637.5*xb4)*(3.0+2.0*xb3+xb5)
     &             -(6.0*xb2+5.0*xb4)*(22713.0*xb3+18927.5*xb5))
            dev4   =dev4p*ano2
            dev5p  =113565.0*(6.0*xb5*(9.0+5.0*xb6)-30.0*xb11)
     &             /((9.0+5.0*xb6)**2)
            dev5   =dev5p*ano
            dev6p  =264985.0*(7.0*xb6*(4.0+10.0*xb7)-70.0*xb13)
     &             /((4.0+10.0*xb7)**2)
            dev6   =dev6p*ann
            devdtv =tv1*tv1*(226.0*xa*(dev1+dev2+dev3)
     &             +3785.5*xb*(dev4+dev5+dev6))
c derivatives of ev with respect to species mole fractions, in k
            dedo   =113565.0*xb6/(9.0 +5.0*xb6)
            dedn   =264985.0*xb7/(4.0+10.0*xb7)
            dedno  =2712.0/(xa12-1.0)
            dedo2  =(2260.0/(xa10-1.0)+(22713.0*xb3+18927.5*xb5)
     &                                /(3.0+2.0*xb3+xb5))
            dedn2  =3390.0/(xa15-1.0)
c
            devdr  =alp1*dedo2+alp4*dedn2
            devdro =amas1(1)*dedo+alp2*dedo2
            devdrn =amas1(2)*dedn+alp5*dedn2
            devrno =amas1(3)*dedno+alp3*dedo2+alp6*dedn2
c
c derivatives of tv by flow variables
            dedtv1 =1.0/(devdtv+1.0e-30)
c            dtvdr  =0.0
c modified by Ishihara
            dtvdr  =dedtv1*devdr
c
            dtvdru =0.0
            dtvdrv =0.0
            dtvde  =0.0
c            dtvdro =-devdro*dedtv1
c            dtvdrn =-devdrn*dedtv1
c            dtvrno =-devrno*dedtv1
c modified by Ishihara
            dtvdro = devdro*dedtv1
            dtvdrn = devdrn*dedtv1
            dtvrno = devrno*dedtv1
c
            dtvdev = dedtv1/8.314
            tvev(i,j) =dtvdev
c ***********************************
c section 3. partial derivatives of t
c ***********************************
            rho1   =1.0/qc(i,j,1)
            ddm    =0.5*rho1*(rhou**2+rhov**2)
c            cvt    =e-ev-ddm-(h0i(1)*ano+h0i(2)*ann+h0i(3)*anno)
c modified by Ishihara
            cvt    =e-ev-0.5d0*(rhou**2+rhov**2)
     &              -(h0i(1)*ano+h0i(2)*ann+h0i(3)*anno)
            cvav   =ano*cvi(1)+ann*cvi(2)+anno*cvi(3)+ano2*cvi(4)
     &                                               +ann2*cvi(5)
            cvav1  =1.0/cvav
            cvav2  =cvav1*cvav1
            t1     =1.0/tt(i,j)
c derivatives of t by the flow variables
c original  dtdr   = cvav1*0.5*(rhou**2+rhov**2)*rho1*rho1
c            dtdr   = cvav2*( 0.5*ddm*rho1*cvav
c     &                      -(alp1*cvi(4)+alp4*cvi(5))*cvt)
c modified by Ishihara
            dtdr   = cvav2*(-ddm*rho1*cvav
     &                      -(alp1*cvi(4)+alp4*cvi(5))*cvt)
            dtdru  =-rhou*cvav1*rho1
            dtdrv  =-rhov*cvav1*rho1
            dtde   = cvav1
            dtdro  = cvav2*(-h0i(1)*amas1(1)*cvav
     &                    -(cvi(1)*amas1(1)+alp2*cvi(4))*cvt)
            dtdrn  = cvav2*(-h0i(2)*amas1(2)*cvav
     &                    -(cvi(2)*amas1(2)+alp5*cvi(5))*cvt)
            dtdrno = cvav2*(-h0i(3)*amas1(3)*cvav
     &                    -(cvi(3)*amas1(3)
     &                     +alp3*cvi(4)+alp6*cvi(5))*cvt)
            dtdro2 =-cvav2*cvt*cvi(4)*amas1(4)
            dtdrn2 =-cvav2*cvt*cvi(5)*amas1(5)
            dtdev  =-cvav1
c *************************************************************
c section 4.
c     >chemical source terms for o, n, and no
c     >reaction rates w1,w2,w3 are for o,n,no, in mole/(m3-sec)
c     >use two-temperature model ta=sqrt(t*tv)
c     >reaction index (1) o2+m=o+o+m
c                     (2) n2+o=no+n
c                     (3) no+o=o2+n
c                     (4) n2+m=n+n+m
c                     (5) no+m=n+o+m
c     >number densities ano,ann,anno,ano2,ann2 are in moles/m3
c *************************************************************
            ta     =(tt(i,j)**qq)*(tv(i,j)**(1.-qq))
            tasqr  =sqrt(ta)
            ta1    =1.0/ta
            tasqr1 =1.0/tasqr
            ta2    =ta1*ta1
            ta3    =ta1*tasqr1
            ta16   =ta**(-1.6)
            z      =ta*1.0e-4
            z1     =10000.0*ta1
            zlog   =log(z1)
c modified by Ishihara
            z      =dmax1(ta*1.0d-4,telim*1.0d-4)
            z1     =dmin1(10000.d0*ta1,10000.d0/telim)
c
            te     =tt(i,j)
            tesqr  =sqrt(te)
            te1    =1.0/te
            te2    =te1*te1
            tesqr1 =1./tesqr
            te3    =te1*tesqr1
            te16   =te**(-1.6)
            z2     =te*1.0e-4
            z21    =10000.0*te1
            z2log  =log(z21)
c modified by Ishihara
            z2     =dmax1(te*1.0d-4,telim*1.0d-4)
            z21    =dmin1(10000.d0*te1,10000.d0/telim)
c modified by Ishihara
c for Bose and Candler NO exchange reaction
            te042  =te**0.42d0
cc reaction (1). reaction rate coefs are from park, 1988 jan aiaa paper
c            cini   =crat1(1)*ano +crat1(2)*ann+crat1(3)*anno
c     &             +crat1(4)*ano2+crat1(5)*ann2
c            exp1f  =exp(-59500.0*ta1)
c            exp1r  =exp(-59500.0*te1-aka1(1)*z2   -aka1(2)
c     &                              -aka1(3)*z2log-aka1(4)*z21
c     &                              -aka1(5)*z21*z21)
c            cta1f  =cini*ta3*exp1f*ano2
c            cte1r  =cini*te3*exp1r*ano*ano
c            w1     =cta1f-cte1r
c            w1f    =cta1f
c            w1r    =cte1r
c reaction (1). reaction rate coefs arefrom Park, JTHT, Vol.15, No1, 2001
c modified by Ishihara
            cini     =crat1(1)*ano +crat1(2)*ann +crat1(3)*anno
     &               +crat1(4)*ano2+crat1(5)*ann2
            exp1f    =dexp(-59360.0d0*ta1)
            exp1r    =dexp(-59360.0d0*te1-aka1(1)*z2 -aka1(2)
     &                -aka1(3)*z2log-aka1(4)*z21-aka1(5)*z21*z21)
            cta1f    =cini*ta3*exp1f*ano2
            cte1r    =cini*te3*exp1r*ano*ano
            w1       =cta1f-cte1r
            w1f    =cta1f
            w1r    =cte1r
cc reaction (2). reaction rate coef is taken from hanson
c            exp2f  =exp(-38400.0*te1)
c            exp2r  =exp(-38400.0*te1-aka2(1)*z2   -aka2(2)
c     &                              -aka2(3)*z2log-aka2(4)*z21
c     &                              -aka2(5)*z21*z21)
c            w2     =crate(2)*te1*(ann2*ano*exp2f-anno*ann*exp2r)
c            w2f    =crate(2)*te1*ann2*ano*exp2f
c            w2r    =crate(2)*te1*anno*ann*exp2r
c reaction (2). reaction rate coef is taken from Bose and Candler
c modified by Ishihara
            exp2f    =dexp(-42938.0d0*te1)
            exp2r    =dexp(-42938.0d0*te1-aka2(1)*z2 -aka2(2)
     &              -aka2(3)*z2log-aka2(4)*z21-aka2(5)*z21*z21)
            w2       =crate(2)*te042*(ann2*ano*exp2f-anno*ann*exp2r)
            w2f      =crate(2)*te042*ann2*ano*exp2f
            w2r      =crate(2)*te042*anno*ann*exp2r
cc reaction (3). reaction rate coef is taken from hanson
c            exp3f  =exp(-19450.0*te1)
c            exp3r  =exp(-19450.0*te1-aka3(1)*z2   -aka3(2)
c     &                              -aka3(3)*z2log-aka3(4)*z21
c     &                              -aka3(5)*z21*z21)
c            w3     =crate(3)*(anno*ano*exp3f-ano2*ann*exp3r)
c            w3f    =crate(3)*anno*ano*exp3f
c            w3r    =crate(3)*ano2*ann*exp3r
c reaction (3). reaction rate coef is taken from Bose and Candler
c modified by Ishihara
            exp3f    =dexp(-19400.0d0*te1)
            exp3r    =dexp(-19400.0d0*te1-aka3(1)*z2 -aka3(2)
     &               -aka3(3)*z2log-aka3(4)*z21-aka3(5)*z21*z21)
            w3       =crate(3)*(anno*ano*exp3f-ano2*ann*exp3r)
            w3f      =crate(3)*anno*ano*exp3f
            w3r      =crate(3)*ano2*ann*exp3r
c reaction (4). reaction rate coefs are from park, 1988 jan aiaa paper
            cini4  =crat4(1)*ano +crat4(2)*ann+crat4(3)*anno
     &             +crat4(4)*ano2+crat4(5)*ann2
            exp4f  =exp(-113200.0*ta1)
            exp4r  =exp(-113200.0*te1-aka4(1)*z2   -aka4(2)
     &                               -aka4(3)*z2log-aka4(4)*z21
     &                               -aka4(5)*z21*z21)
            cta4f  =cini4*ta16*exp4f*ann2
            cte4r  =cini4*te16*exp4r*ann*ann
            w4     =cta4f-cte4r
            w4f    =cta4f
            w4r    =cte4r
c reaction (5). reaction rate coefs are from park, 1988 jan aiaa paper
c modified by K.NIIZUMA
            cini5  =crat5(1)*ano +crat5(2)*ann+crat5(3)*anno
     &             +crat5(4)*ano2+crat5(5)*ann2
            exp5f  =exp(-75500.0*ta1)
            exp5r  =exp(-75500.0*te1-aka5(1)*z2   -aka5(2)
     &                              -aka5(3)*z2log-aka5(4)*z21
     &                              -aka5(5)*z21*z21)
            cta5f  =cini5*exp5f*anno
            cte5r  =cini5*exp5r*ann*ano
            w5     =cta5f-cte5r
            w5f    =cta5f
            w5r    =cte5r
c rate of change of number density of o,n,no (mole/m3-s)
            dodt   = 2.0*w1-w2-w3+w5
            dndt   = w2+w3+2.0*w4+w5
            dnodt  = w2-w3-w5
            do2dt  =-w1+w3
            dn2dt  =-w2-w4
c
c *****************************************
c section 5. derivaties of reaction sources
c *****************************************
c derivatives of reaction sources w1,w2,w3 with respect to ta
c modified by K.NIIZUMA
c            dw1dta =ta1*(-1.5 +59500.0*ta1)*cta1f
c modified by Ishihara
            dw1dta   =ta1*(-1.5d0+59360.0d0*ta1)*cta1f
            dw2dta =0.0
            dw3dta =0.0
            dw4dta =ta1*(-1.6+113200.0*ta1)*cta4f
            dw5dta =ta1*       75500.0*ta1 *cta5f
c            dw1dt  =-te1*(-1.5+59500.0*te1)*cte1r
c     &             +cini*te3*exp1r*ano*ano*(-aka1(3)*te1
c     &             +aka1(1)*1.0e-4-aka1(4)*te1*z21
c     &             -2.0*te1*aka1(5)*z21*z21)
c modified by Ishihara
            dw1dt    =-te1*(-1.5d0+59360.0d0*te1)*cte1r
     &               +cini*te3*exp1r*ano*ano
     &               *(-aka1(3)*te1    +    aka1(1)*1.0d-4
     &               -aka1(4)*te1*z21-2.0d0*aka1(5)*te1*z21*z21)
c            dw2dt  = te1*(-1.0+38400.0*te1)*w2
c     &             +crate(2)*te1*exp2r*anno*ann*(-aka2(3)*te1
c     &             +aka2(1)*1.e-4-te1*aka2(4)*z21
c     &             -2.0*te1*aka2(5)*z21*z21)
c Bose and Candler modified by Ishihara
            dw2dt    = te1*(0.42d0+42938.0d0*te1)*w2
     &               +crate(2)*te1*exp2r*anno*ann
     &               *(-aka2(3)*te1      +    aka2(1)*1.0d-4
     &                 -aka2(4)*te1*z21-2.0d0*aka2(5)*te1*z21*z21)
c            dw3dt  = te1*      19450.0*te1 *w3
c     &             +crate(3)    *exp3r*ano2*ann*(-aka3(3)*te1
c     &             +aka3(1)*1.e-4-te1*aka3(4)*z21
c     &             -2.0*te1*aka3(5)*z21*z21)
c Bose and Candler modified by Ishihara
            dw3dt    = te1*19400.0d0*te1*w3
     &               +crate(3)*exp3r*ano2*ann
     &               *(-aka3(3)*te1      +    aka3(1)*1.0d-4
     &                 -aka3(4)*te1*z21-2.0d0*aka3(5)*te1*z21*z21)
            dw4dt  =-te1*(-1.6+113200.0*te1)*cte4r
     &             +cini4*te16*exp4r*ann*ann*(-aka4(3)*te1
     &             +aka4(1)*1.e-4-aka4(4)*te1*z21
     &             -2.0*te1*aka4(5)*z21*z21)
            dw5dt  =-te1*       75500.0*te1 *cte5r
     &             +cini5     *exp5r*ann*ano*(-aka5(3)*te1
     &             +aka5(1)*1.e-4-aka5(4)*te1*z21
     &             -2.0*te1*aka5(5)*z21*z21)
c modification ends
c derivatives of average temperature ta by flow variables
            ttinfqr =ta*tv1
            tvtsqr  =ta*t1
c dtadx's are 0.5*(sqrt(tv/t)*(dt/dx)+sqrt(t/tv)*(dtv/dx)), etc
            dtadr  =qq*tvtsqr*dtdr  +(1.-qq)*ttinfqr*dtvdr
            dtadru =qq*tvtsqr*dtdru +(1.-qq)*ttinfqr*dtvdru
            dtadrv =qq*tvtsqr*dtdrv +(1.-qq)*ttinfqr*dtvdrv
            dtade  =qq*tvtsqr*dtde  +(1.-qq)*ttinfqr*dtvde
            dtadro =qq*tvtsqr*dtdro +(1.-qq)*ttinfqr*dtvdro
            dtadrn =qq*tvtsqr*dtdrn +(1.-qq)*ttinfqr*dtvdrn
            dtarno =qq*tvtsqr*dtdrno+(1.-qq)*ttinfqr*dtvrno
            dtadev =qq*tvtsqr*dtdev +(1.-qq)*ttinfqr*dtvdev
            cini1  =1.0/cini
            cini41 =1.0/cini4
            cini51 =1.0/cini5
c derivatives of w1
      dw1dro =(crat1(1)*amas1(1)+alp2*crat1(4))*cini1*w1
     &       +cini*(ta3*exp1f*alp2-te3*2.0*ano*exp1r*amas1(1))
     &       +dw1dta*dtadro+dw1dt*dtdro
      dw1drn =(crat1(2)*amas1(2)+alp5*crat1(5))*cini1*w1
     &       +dw1dta*dtadrn+dw1dt*dtdrn
      dw1drno=(crat1(3)*amas1(3)+alp3*crat1(4)+alp6*crat1(5))*cini1*w1
     &       +cini*ta3*exp1f*alp3
     &       +dw1dta*dtarno+dw1dt*dtdrno
      dw1dev =dw1dta*dtadev+dw1dt*dtdev
c derivatives of w2
c      dw2dro  = crate(2)*te1*exp2f*ann2*amas1(1)+dw2dt*dtdro
c      dw2drn  = crate(2)*te1*(exp2f*ano*alp5-exp2r*anno*amas1(2))
c     &          +dw2dt*dtdrn
c      dw2drno = crate(2)*te1*(exp2f*ano*alp6-exp2r*ann*amas1(3))
c     &          +dw2dt*dtdrno
c      dw2dev  =dw2dta*dtadev+dw2dt*dtdev
c Bose and Candler modified by Ishihara
      dw2dro   =crate(2)*te042*exp2f*ann2*amas1(1)+dw2dt*dtdro
      dw2drn   =crate(2)*te042*(exp2f*ano*alp5
     &                         -exp2r*anno*amas1(2))+dw2dt*dtdrn
      dw2drno  =crate(2)*te042*(exp2f*ano*alp6
     &                         -exp2r*ann*amas1(3))+dw2dt*dtdrno
c      dw2dro2  =dw2dt*dtdro2
c      dw2drn2  =crate(2)*te042*exp2f*ano *amas1(5)+dw2dt*dtdrn2
      dw2dev   =dw2dta*dtadev+dw2dt*dtdev
c derivatives of w3
      dw3dro  = crate(3)*(exp3f*anno*amas1(1)-exp3r*ann*alp2)
     &         +dw3dt*dtdro
      dw3drn  =-crate(3)*exp3r*ano2*amas1(2)+dw3dt*dtdrn
      dw3drno = crate(3)*(exp3f*ano*amas1(3)-exp3r*ann*alp3)
     &         +dw3dt*dtdrno
      dw3dev  =dw3dta*dtadev+dw3dt*dtdev
c derivatives of w4
      dw4dro =(crat4(1)*amas1(1)+alp2*crat1(4))*cini41*w4
     &        +dw4dta*dtadro+dw4dt*dtdro
      dw4drn =(crat4(2)*amas1(2)+alp5*crat1(5))*cini41*w4
     &        +cini4*(ta16*exp4f*alp5-te16*2.0*ann*exp4r*amas1(2))
     &        +dw4dta*dtadrn+dw4dt*dtdrn
      dw4drno=(crat4(3)*amas1(3)+alp3*crat4(4)+alp6*crat4(5))*cini41*w4
     &        +cini4*ta16*exp4f*alp6
     &        +dw4dta*dtarno+dw4dt*dtdrno
      dw4dev  =dw4dta*dtadev+dw4dt*dtdev
c derivatives of w5
c modified by K.NIIZUMA
      dw5dro =(crat5(1)*amas1(1)+alp2*crat5(4))*cini51*w5
     &       -cini5*ann*exp5r*amas1(1)
     &       +dw5dta*dtadro+dw5dt*dtdro
      dw5drn =(crat5(2)*amas1(2)+alp5*crat5(5))*cini51*w5
     &       -cini5*ano*exp5r*amas1(2)
     &       +dw5dta*dtadrn+dw5dt*dtdrn
      dw5drno=(crat5(3)*amas1(3)+alp3*crat5(4)+alp6*crat5(5))*cini51*w5
     &       +cini5*    exp5f*amas1(3)
     &       +dw5dta*dtarno+dw5dt*dtdrno
      dw5dev  =dw5dta*dtadev+dw5dt*dtdev
c modification ends
c section 6.  >fill in source sc and jacobian sj elements
c              for species variables
c             >derivatives of number densities of o,n,no
c              (mole/m3-sec) by flow variables
            scs(i,j,1)  =amass(1)* dodt
            scs(i,j,2)  =amass(2)* dndt
            scs(i,j,3)  =amass(3)*dnodt
c
            sjac16      =amass(1)*(2.0*dw1dro -dw2dro -dw3dro +dw5dro )
            sjac17      =amass(1)*(2.0*dw1drn -dw2drn -dw3drn +dw5drn )
            sjac18      =amass(1)*(2.0*dw1drno-dw2drno-dw3drno+dw5drno)
            sjac111     =amass(1)*(2.0*dw1dev -dw2dev -dw3dev +dw5dev )
c
            sjac26      =amass(2)*(dw2dro +dw3dro +2.0*dw4dro +dw5dro )
            sjac27      =amass(2)*(dw2drn +dw3drn +2.0*dw4drn +dw5drn )
            sjac28      =amass(2)*(dw2drno+dw3drno+2.0*dw4drno+dw5drno)
            sjac211     =amass(2)*(dw2dev +dw3dev +2.0*dw4dev +dw5dev )
c
            sjac36      =amass(3)*(dw2dro -dw3dro -dw5dro )
            sjac37      =amass(3)*(dw2drn -dw3drn -dw5drn )
            sjac38      =amass(3)*(dw2drno-dw3drno-dw5drno)
            sjac311     =amass(3)*(dw2dev -dw3dev -dw5dev )
c
            sjac46  = amass(4)*(-dw1dro +dw3dro )
            sjac47  = amass(4)*(-dw1drn +dw3drn )
            sjac48  = amass(4)*(-dw1drno+dw3drno)
            sjac411 = amass(4)*(-dw1dev+dw3dev)
c
            sjac56  = amass(5)*(-dw2dro -dw4dro )
            sjac57  = amass(5)*(-dw2drn -dw4drn )
            sjac58  = amass(5)*(-dw2drno-dw4drno)
            sjac511 = amass(5)*(-dw2dev-dw4dev)
c
c diagonal terms for implicit scheme of Eberhardt & Imlay
c
            aa11         =min(abs(sjac16 ),rmx)
            aa12         =min(abs(sjac17 ),rmx)
            aa13         =min(abs(sjac18 ),rmx)
            aa16         =min(abs(sjac111),rmx)
            ajs(i,j,1)  =sqrt(aa11*aa11+aa12*aa12+aa13*aa13
     &                       +aa16*aa16)
c
            aa21         =min(abs(sjac26 ),rmx)
            aa22         =min(abs(sjac27 ),rmx)
            aa23         =min(abs(sjac28 ),rmx)
            aa26         =min(abs(sjac211),rmx)
            ajs(i,j,2)  =sqrt(aa21*aa21+aa22*aa22+aa23*aa23
     &                       +aa26*aa26)
c
            aa31         =min(abs(sjac36 ),rmx)
            aa32         =min(abs(sjac37 ),rmx)
            aa33         =min(abs(sjac38 ),rmx)
            aa36         =min(abs(sjac311),rmx)
            ajs(i,j,3)  =sqrt(aa31*aa31+aa32*aa32+aa33*aa33
     &                       +aa36*aa36)
c
c *************************************************************
c section 7. vibrational relaxation times and their derivatives
c *************************************************************
            t13         =tt(i,j)**(-0.3333333)
            h85t        =h85*tt(i,j)
c            exp20       =exp(-20.0*t13)
c            exp11       =exp20*exp20
c            exp13       =exp11
c            exp21       =exp11
c            exp22       =exp11*exp20
c            exp31       =exp22
c            exp12       =exp22*exp22
c            exp14       =exp12
c            exp15       =exp12
c            exp23       =exp12
c            exp24       =exp12
c            exp25       =exp12*exp11
c            exp32       =exp25
c            exp33       =exp25
c            exp34       =exp25
c            exp35       =exp25*exp22
c modified by Ishihara
            exp11    =dexp(-a11*t13)
            exp12    =dexp(-a12*t13)
            exp13    =dexp(-a13*t13)
            exp14    =dexp(-a14*t13)
            exp15    =dexp(-a15*t13)
            exp21    =dexp(-a21*t13)
            exp22    =dexp(-a22*t13)
            exp23    =dexp(-a23*t13)
            exp24    =dexp(-a24*t13)
            exp25    =dexp(-a25*t13)
            exp31    =dexp(-a31*t13)
            exp32    =dexp(-a32*t13)
            exp33    =dexp(-a33*t13)
            exp34    =dexp(-a34*t13)
            exp35    =dexp(-a35*t13)
c
            ant         =qc(i,j,1)*avmas1
            ant2        =ant*ant
            f1          =h85t*(ano*exb11*exp11+ann*exb12*exp12
     &                        +anno*exb13*exp13+ano2*exb14*exp14
     &                        +ann2*exb15*exp15)
            f2          =h85t*(ano*exb21*exp21+ann*exb22*exp22
     &                        +anno*exb23*exp23+ano2*exb24*exp24
     &                        +ann2*exb25*exp25)
            f3          =h85t*(ano*exb31*exp31+ann*exb32*exp32
     &                        +anno*exb33*exp33+ano2*exb34*exp34
     &                        +ann2*exb35*exp35)
            t2          =sqrt(tt(i,j))
            t3          =tt(i,j)*t2
c taua's are reciprocals of relaxation time tau
            x1          =g214*t3*f1+ant
            x2          =g214*t3*f2+ant
            x3          =g214*t3*f3+ant
            xx1         =1.0/x1
            xx2         =1.0/x2
            xx3         =1.0/x3
            xx12        =xx1*xx1
            xx22        =xx2*xx2
            xx32        =xx3*xx3
            taua1       =ant*f1*xx1
            taua2       =ant*f2*xx2
            taua3       =ant*f3*xx3
            tau1        =1.0/taua1
            tau2        =1.0/taua2
            tau3        =1.0/taua3
            p           =8.205e-5*ansum*tt(i,j)
c derivatives of f and taua by t
            df1dt       =f1*t1+0.33333*h85*t13*(ano*a11*exb11*exp11
     &                    +ann*a12*exb12*exp12+anno*a13*exb13*exp13
     &                   +ano2*a14*exb14*exp14+ann2*a15*exb15*exp15)
            df2dt       =f2*t1+0.33333*h85*t13*(ano*a21*exb21*exp21
     &                    +ann*a22*exb22*exp22+anno*a23*exb23*exp23
     &                   +ano2*a24*exb24*exp24+ann2*a25*exb25*exp25)
            df3dt       =f3*t1+0.33333*h85*t13*(ano*a31*exb31*exp31
     &                    +ann*a32*exb32*exp32+anno*a33*exb33*exp33
     &                   +ano2*a34*exb34*exp34+ann2*a35*exb35*exp35)
            dtu1dt      =ant*df1dt*xx1
     &                  -g214*t2*(1.5*f1+tt(i,j)*df1dt)*ant*f1*xx12
            dtu2dt      =ant*df2dt*xx2
     &                  -g214*t2*(1.5*f2+tt(i,j)*df2dt)*ant*f2*xx22
            dtu3dt      =ant*df3dt*xx3
     &                  -g214*t2*(1.5*f3+tt(i,j)*df3dt)*ant*f3*xx32
c derivatives of taua by number density of colliding particles
c (mole/m3)
            df1dn1      =h85t*exb11*exp11
            df1dn2      =h85t*exb12*exp12
            df1dn3      =h85t*exb13*exp13
            df1dn4      =h85t*exb14*exp14
            df1dn5      =h85t*exb15*exp15
            df2dn1      =h85t*exb21*exp21
            df2dn2      =h85t*exb22*exp22
            df2dn3      =h85t*exb23*exp23
            df2dn4      =h85t*exb24*exp24
            df2dn5      =h85t*exb25*exp25
            df3dn1      =h85t*exb31*exp31
            df3dn2      =h85t*exb32*exp32
            df3dn3      =h85t*exb33*exp33
            df3dn4      =h85t*exb34*exp34
            df3dn5      =h85t*exb35*exp35
c
            dtu1do  =ant2*df1dn1*xx12
            dtu1dn  =ant2*df1dn2*xx12
            dtu1dno =ant2*df1dn3*xx12
            dtu1do2 =ant2*df1dn4*xx12
            dtu1dn2 =ant2*df1dn5*xx12
            dtu2do  =ant2*df2dn1*xx22
            dtu2dn  =ant2*df2dn2*xx22
            dtu2dno =ant2*df2dn3*xx22
            dtu2do2 =ant2*df2dn4*xx22
            dtu2dn2 =ant2*df2dn5*xx22
            dtu3do  =ant2*df3dn1*xx32
            dtu3dn  =ant2*df3dn2*xx32
            dtu3dno =ant2*df3dn3*xx32
            dtu3do2 =ant2*df3dn4*xx32
            dtu3dn2 =ant2*df3dn5*xx32
c
            dt1do      =dtu1dt*dtdro +dtu1do *amas1(1)+dtu1do2*alp2
            dt1dn      =dtu1dt*dtdrn +dtu1dn *amas1(2)+dtu1dn2*alp4
            dt1dno     =dtu1dt*dtdrno+dtu1dno*amas1(3)
     &                               +dtu1do2*alp3+dtu1dn2*alp6
            dt2do      =dtu2dt*dtdro +dtu2do *amas1(1)+dtu2do2*alp2
            dt2dn      =dtu2dt*dtdrn +dtu2dn *amas1(2)+dtu2dn2*alp4
            dt2dno     =dtu2dt*dtdrno+dtu2dno*amas1(3)
     &                               +dtu2do2*alp3+dtu1dn2*alp6
            dt3do      =dtu3dt*dtdro +dtu3do *amas1(1)+dtu3do2*alp2
            dt3dn      =dtu3dt*dtdrn +dtu3dn *amas1(2)+dtu3dn2*alp4
            dt3dno     =dtu3dt*dtdrno+dtu3dno*amas1(3)
     &                               +dtu3do2*alp3+dtu3dn2*alp6
c derivatives of taua by flow variables
            dt1dev      =dtu1dt*dtdev
            dt2dev      =dtu2dt*dtdev
            dt3dev      =dtu3dt*dtdev
c ******************************************************
c section 8. vibrational source term and its derivatives
c            species order 1=no,2=o2,3=n2
c ******************************************************
            xae         =exp(226.0*t1)
            xae2        =xae*xae
            xae4        =xae2*xae2
            xae8        =xae4*xae4
            xae10       =xae8*xae2
            xae12       =xae10*xae2
            xae15       =xae12*xae2*xae
            eve1        =2712.0/(xae12-1.0)
            eve2        =2260.0/(xae10-1.0)
            eve3        =3390.0/(xae15-1.0)
            eva1        =ev3
            eva2        =ev1
            eva3        =ev2
c modified by K.NIIZUMA
            ts          =tss(i)
            s           =3.5*exp(-5000.0/ts)
c modification ends
c alternative method
c           tsvs        =(t-tv(i,j))/(ts-tv(i,j))
c           ft          =abs(tsvs)
c           fts         =ft**(s-1.)
c           fts1        =ft**(s-2.)
            tsvs        =(ts-tv(i,j))/(ts-tvinf)
            ft          =abs(tsvs)
            fts         =ft**(s-1.0)
            fts1        =ft**(s-2.0)
c modified by Ishihara
            if(nrlx.eq.0) then
               fts      =1.d0
               fts1     =0.d0
cd
            end if
            wv1         =8.314*taua1*anno*(eve1-eva1)*fts
            wv2         =8.314*taua2*ano2*(eve2-eva2)*fts
            wv3         =8.314*taua3*ann2*(eve3-eva3)*fts
            wv          =wv1+wv2+wv3
c derivatives of evek=eve/k by t, and evk=ev/k by tv
            devek1      =(eve1*t1)**2*xae12
            devek2      =(eve2*t1)**2*xae10
            devek3      =(eve3*t1)**2*xae15
            devk1       =(eva1*tv1)**2*xa12
            devk2       =(eva2*tv1)**2*xa10
            devk3       =(eva3*tv1)**2*xa15
            xb          =exp(-3785.5*tv1)
            yy          =3785.5*tv1*tv1*xb
            devk4       =dev4p*yy
            devk5       =dev5p*yy
            devk6       =dev6p*yy
c alternative method
c           d1          = (s-1.0)*fts1/(ts-tv(i,j))
c           d2          =-(s-1.0)*fts1*(ts-t)/((ts-tv(i,j)**2)
            d1          =0.0
            d2          =-(s-1.0)*fts1/(ts-tvinf)
c parameters bit and cit
            b1t         =devek1*fts+(eve1-eva1)*d1
            b2t         =devek2*fts+(eve2-eva2)*d1
            b3t         =devek3*fts+(eve3-eva3)*d1
            c1tv        =-devk1*fts+(eve1-eva1)*d2
            c2tv        =-devk2*fts+(eve2-eva2)*d2
            c3tv        =-devk3*fts+(eve3-eva3)*d2
c
c derivatives of wvi with respect to flow variables
c
            dv1ev       =wv1*dt1dev*tau1
     &                  +8.314*taua1*anno*(b1t*dtdev+c1tv*dtvdev)
            dv2ev       =wv2*dt2dev*tau2
     &                  +8.314*taua2*ano2*(b2t*dtdev+c2tv*dtvdev)
            dv3ev       =wv3*dt3dev*tau3
     &                  +8.314*taua3*ann2*(b3t*dtdev+c3tv*dtvdev)
            dv1dro      =wv1*dt1do*tau1
     &                  +8.314*taua1*anno*(b1t*dtdro+c1tv*dtvdro)
            dv1drn      =wv1*dt1dn*tau1
     &                  +8.314*taua1*anno*(b1t*dtdrn+c1tv*dtvdrn)
            dv1rno      =wv1*dt1dno*tau1
     &                  +8.314*taua1*(amas1(3)*(eve1-eva1)*fts1
     &                              +anno*(b1t*dtdrno+c1tv*dtvrno))
            dv2dro      =wv2*dt2do*tau2
     &                  +8.314*taua2*(alp2*amas1(4)*(eve2-eva2)*fts1
     &                               +ano2*(b2t*dtdro+c2tv*dtvdro))
            dv2drn      =wv2*dt1dn*tau2
     &                  +8.314*taua2*ano2*(b2t*dtdrn+c2tv*dtvdrn)
            dv2rno      =wv2*dt2dno*tau2
     &                  +8.314*taua2*(alp3*amas1(4)*(eve2-eva2)*fts1
     &                              +ano2*(b2t*dtdrno+c2tv*dtvrno))
            dv3dro      =wv3*dt3do*tau3
     &                  +8.314*taua3*ann2*(b3t*dtdro+c3tv*dtvdro)
            dv3drn      =wv3*dt1dn*tau3
     &                  +8.314*taua3*(alp5*amas1(5)*(eve3-eva3)*fts1
     &                              +ann2*(b3t*dtdrn+c3tv*dtvdrn))
            dv3rno      =wv3*dt3dno*tau3
     &                  +8.314*taua3*(alp6*amas1(5)*(eve3-eva3)*fts1
     &                              +ann2*(b3t*dtdrno+c3tv*dtvrno))
c
            dwvdev      =dv1ev+dv2ev+dv3ev
            dwvdro      =dv1dro+dv2dro+dv3dro
            dwvdrn      =dv1drn+dv2drn+dv3drn
            dwvrno      =dv1rno+dv2rno+dv3rno
c
c ************************************************
c section 9. assemble the vibrational source terms
c ************************************************
c wvt is the total source for vibronic energy in j/(m3-sec)
c           wvt=wv+8.314*((ebnok)*dnodt+(ebo2k+ev4)
c     &        *do2dt+(ebn2k)*dn2dt+ev5*dodt+ev6*dndt)
        wvt     =wv+8.314*(ebnok*dnodt+ebo2k*do2dt+ebn2k*dn2dt)
c
c derivatives of wvt
c
        dpd1dt =  ebnok     *(          dw2dt-dw3dt          -dw5dt)
        dpd2dt = (ebo2k+ev4)*(   -dw1dt      +dw3dt                )
     &                 +ev5 *(2.0*dw1dt-dw2dt-dw3dt          +dw5dt)
        dpd3dt =  ebn2k     *(         -dw2dt          -dw4dt      )
     &                 +ev6 *(          dw2dt+dw3dt+2.0*dw4dt+dw5dt)
c
        dpd1dro = dpd1dt*dtdro +ebnok*sjac36 *amas1(3)
        dpd1drn = dpd1dt*dtdrn +ebnok*sjac37 *amas1(3)
        dpd1rno = dpd1dt*dtdrno+ebnok*sjac38 *amas1(3)
        dpd1dev = dpd1dt*dtdev +ebnok*sjac311*amas1(3)
c
        dpd2dro = dpd2dt*dtdro +(ebo2k+ev4)*sjac46 *amas1(4)
     &                                +ev5 *sjac16 *amas1(1)
        dpd2drn = dpd2dt*dtdrn +(ebo2k+ev4)*sjac47 *amas1(4)
     &                                +ev5 *sjac17 *amas1(1)
        dpd2rno = dpd2dt*dtdrno+(ebo2k+ev4)*sjac48 *amas1(4)
     &                                +ev5 *sjac18 *amas1(1)
        dpd2dev = dpd2dt*dtdev +(ebo2k+ev4)*sjac411*amas1(4)
     &                                +ev5 *sjac111*amas1(1)
c
        dpd3dro = dpd3dt*dtdro +ebn2k*sjac56 *amas1(5)
     &                         +ev6  *sjac26 *amas1(2)
        dpd3drn = dpd3dt*dtdrn +ebn2k*sjac57 *amas1(5)
     &                         +ev6  *sjac27 *amas1(2)
        dpd3rno = dpd3dt*dtdrno+ebn2k*sjac58 *amas1(5)
     &                         +ev6  *sjac28 *amas1(2)
        dpd3dev = dpd3dt*dtdev +ebn2k*sjac511*amas1(5)
     &                         +ev6  *sjac211*amas1(2)
c
        dpddro  = dpd1dro+dpd2dro+dpd3dro
        dpddrn  = dpd1drn+dpd2drn+dpd3drn
        dpdrno  = dpd1rno+dpd2rno+dpd3rno
        dpddev  = dpd1dev+dpd2dev+dpd3dev
c
        dwtdev  = dwvdev+dpddev
        dwtdro  = dwvdro+dpddro
        dwtdrn  = dwvdrn+dpddrn
        dwtrno  = dwvrno+dpdrno
c
c *********************************************************
c section 10. >fill in source sr and jacobian sjac elements
c              for vibrational variables
c             >index i=1:vibrational energy, 2:o, 3:n, 4:no
c *********************************************************
            scv(i,j,1)  =wvt
            sjac66      =dwtdro
            sjac67      =dwtdrn
            sjac68      =dwtrno
            sjac611     =dwtdev
            a61         =min(abs(sjac66 ),rmx)
            a62         =min(abs(sjac67 ),rmx)
            a63         =min(abs(sjac68 ),rmx)
            a66         =min(abs(sjac611),rmx)
            ajv(i,j,1)  =abs(a66)
        end do
      end do
c
c* termination
      return
      end
c
      subroutine heat
c     ************************************************************
c     *     set boundary values                                  *
c     ************************************************************
      implicit real*8 (a-h,o-z)
c      parameter(nx=101,ny=101)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /cntl5/ nprc,iprc
      common /stat1/ loop,locl
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /wbcd1/ twall,cefn,cefo
c
      dimension pw(0:nx,2),p1(0:nx,2),s(0:nx)
c
      open(99,file='./output/qw', status='unknown')
c      open(98,file='./output/qwalongy',status='unknown')
!      open(97,file='./output/surfvariable',status='unknown')
!      open(96,file='./output/qh',status='unknown')
!      open(95,file='./output/xcwall',status='unknown')
      write(99,101)"#","  s"," Q_tot [J/m2]"," Q_T [J/m2]"
     &             ," Q_Tv [J/m2]"," Q_hs [J/m2]"
c      write(98,101)"#","  y"," Q_tot [J/m2]"," Q_T [J/m2]"
c     &             ," Q_Tv [J/m2]"," Q_hs [J/m2]"
!      write(97,101)"#","  s"," Tt [K]"," Tv [K]", "velo [m2]"
!     $             , "xc1", "xc2", "xc3", "xc4", "xc5"
!      write(96,101)"#","  s", "qh", "qh1", "qh2", "qh3"
!     &                , "qh4", "qh5"
!      write(95,101)"#","  s", "xc1w", "xc2w", "xc3w", "xc4w"
!     &                , "xc5w"
c
      rwall   = sqrt(x(1,1)**2+y(1,1)**2)
      pw(0,1) = -rwall
      pw(0,2) = 0.d0
      s(0) = 0.d0
c
      do i=1,imax-1
         pw(i,1) = (x(i,1)+x(i+1,1))/2.
         pw(i,2) = (y(i,1)+y(i+1,1))/2.
         p1(i,1) = (x(i,1)+x(i+1,1)+x(i,2)+x(i+1,2))/4.
         p1(i,2) = (y(i,1)+y(i+1,1)+y(i,2)+y(i+1,2))/4.
         dx      = sqrt((pw(i,1)-pw(i-1,1))**2+(pw(i,2)-pw(i-1,2))**2)
         s(i)    = s(i-1)+dx
         dl1     = sqrt((p1(i,1)-pw(i,1))**2+(p1(i,2)-pw(i,2))**2)
c
         rhow    = (qc(i,0,1)+qc(i,1,1))/2.
         ds1w    = (ds(i,0,1)+ds(i,1,1))/2.
         ds2w    = (ds(i,0,2)+ds(i,1,2))/2.
         ds3w    = (ds(i,0,3)+ds(i,1,3))/2.
         ds4w    = (ds(i,0,4)+ds(i,1,4))/2.
         ds5w    = (ds(i,0,5)+ds(i,1,5))/2.
         hs1w    = (hs(i,0,1)+hs(i,1,1))/2.
         hs2w    = (hs(i,0,2)+hs(i,1,2))/2.
         hs3w    = (hs(i,0,3)+hs(i,1,3))/2.
         hs4w    = (hs(i,0,4)+hs(i,1,4))/2.
         hs5w    = (hs(i,0,5)+hs(i,1,5))/2.
         xc1w    = (xc(i,0,1)+xc(i,1,1))/2.
         xc2w    = (xc(i,0,2)+xc(i,1,2))/2.
         xc3w    = (xc(i,0,3)+xc(i,1,3))/2.
         xc4w    = (xc(i,0,4)+xc(i,1,4))/2.
         xc5w    = (xc(i,0,5)+xc(i,1,5))/2.
         vkpw    = (vkap(i,0)+vkap(i,1))/2.
         vkvw    = (vkav(i,0)+vkav(i,1))/2.
c
         dtt     = tt(i,1)  -twall
         dtv     = tv(i,1)  -twall
         dxc1    = xc(i,1,1)-xc1w
         dxc2    = xc(i,1,2)-xc2w
         dxc3    = xc(i,1,3)-xc3w
         dxc4    = xc(i,1,4)-xc4w
         dxc5    = xc(i,1,5)-xc5w
         dxh1    = ds1w*dxc1*hs1w
         dxh2    = ds2w*dxc2*hs2w
         dxh3    = ds3w*dxc3*hs3w
         dxh4    = ds4w*dxc4*hs4w
         dxh5    = ds5w*dxc5*hs5w
c
         qwtt    = vkpw*dtt
         qwtv    = vkvw*dtv
         qwh     = rhow*(dxh1+dxh2+dxh3+dxh4+dxh5)
c
         qw      = (qwtt+qwtv+qwh)/dl1
         qwtt    = qwtt/dl1
         qwtv    = qwtv/dl1
         qwh     = qwh /dl1
c
!         velo1 = dsqrt(uc(i,1)**2+vc(i,1)**2)
c
         write(99,100) s(i),qw,qwtt,qwtv,qwh,dtt,dtv
         yg  =0.5d0*(y(i,1)+y(i+1,1))
c         write(98,100) yg,qw,qwtt,qwtv,qwh
!         write(97,100) s(i),tt(i,1),tv(i,1),velo1,xc(i,1,1),xc(i,1,2)
!     &                 ,xc(i,1,3),xc(i,1,4),xc(i,1,5)
!         write(96,100) s(i),qwh,rhow*dxh1/dl1,rhow*dxh2/dl1
!     &                 ,rhow*dxh3/dl1,rhow*dxh4/dl1,rhow*dxh5/dl1
!         write(95,100) s(i),xc1w,xc2w,xc3w,xc4w,xc5w
      enddo
c
      close (99)
c      close (98)
!      close (97)
!      close (96)
!      close (95)
c
 100  format(1p9e14.6)
 101  format(a,a13,8a14)
c
      return
      end

      subroutine heatrms
c     ************************************************************
c     *     calculation heat rms  each time step                 *
c     ************************************************************
      implicit real*8 (a-h,o-z)
c      parameter(nx=101,ny=101)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /cntl5/ nprc,iprc
      common /stat1/ loop,locl
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /wbcd1/ twall,cefn,cefo
      common /rms/ rmshflx !for heat residual
c
      dimension pw(0:nx,2),p1(0:nx,2),s(0:nx)
      dimension qw(0:nx),preqw(0:nx)
c
c
! boundary
      rwall   = sqrt(x(1,1)**2+y(1,1)**2)
      pw(0,1) = -rwall  !x-cordinate at center
      pw(0,2) = 0.d0    !y-cordinate at center
      s(0) = 0.d0
! reset rmshflx
      rmshflx=0.0d0
!
      do i=1,imax-1
         pw(i,1) = (x(i,1)+x(i+1,1))/2.
         pw(i,2) = (y(i,1)+y(i+1,1))/2.
         p1(i,1) = (x(i,1)+x(i+1,1)+x(i,2)+x(i+1,2))/4.
         p1(i,2) = (y(i,1)+y(i+1,1)+y(i,2)+y(i+1,2))/4.
         dx      = sqrt((pw(i,1)-pw(i-1,1))**2+(pw(i,2)-pw(i-1,2))**2)
         s(i)    = s(i-1)+dx
         dl1     = sqrt((p1(i,1)-pw(i,1))**2+(p1(i,2)-pw(i,2))**2)
c
         rhow    = (qc(i,0,1)+qc(i,1,1))/2.
         ds1w    = (ds(i,0,1)+ds(i,1,1))/2.
         ds2w    = (ds(i,0,2)+ds(i,1,2))/2.
         ds3w    = (ds(i,0,3)+ds(i,1,3))/2.
         ds4w    = (ds(i,0,4)+ds(i,1,4))/2.
         ds5w    = (ds(i,0,5)+ds(i,1,5))/2.
         hs1w    = (hs(i,0,1)+hs(i,1,1))/2.
         hs2w    = (hs(i,0,2)+hs(i,1,2))/2.
         hs3w    = (hs(i,0,3)+hs(i,1,3))/2.
         hs4w    = (hs(i,0,4)+hs(i,1,4))/2.
         hs5w    = (hs(i,0,5)+hs(i,1,5))/2.
         xc1w    = (xc(i,0,1)+xc(i,1,1))/2.
         xc2w    = (xc(i,0,2)+xc(i,1,2))/2.
         xc3w    = (xc(i,0,3)+xc(i,1,3))/2.
         xc4w    = (xc(i,0,4)+xc(i,1,4))/2.
         xc5w    = (xc(i,0,5)+xc(i,1,5))/2.
         vkpw    = (vkap(i,0)+vkap(i,1))/2.
         vkvw    = (vkav(i,0)+vkav(i,1))/2.
c
         dtt     = tt(i,1)  -twall
         dtv     = tv(i,1)  -twall
         dxc1    = xc(i,1,1)-xc1w
         dxc2    = xc(i,1,2)-xc2w
         dxc3    = xc(i,1,3)-xc3w
         dxc4    = xc(i,1,4)-xc4w
         dxc5    = xc(i,1,5)-xc5w
         dxh1    = ds1w*dxc1*hs1w
         dxh2    = ds2w*dxc2*hs2w
         dxh3    = ds3w*dxc3*hs3w
         dxh4    = ds4w*dxc4*hs4w
         dxh5    = ds5w*dxc5*hs5w
c
         qwtt    = vkpw*dtt
         qwtv    = vkvw*dtv
         qwh     = rhow*(dxh1+dxh2+dxh3+dxh4+dxh5)
c
! for rms hflx
         preqw(i) = qw(i)
         qw(i)      = (qwtt+qwtv+qwh)/dl1
!         qwtt    = qwtt/dl1
!         qwtv    = qwtv/dl1
!         qwh     = qwh /dl1
! clac rmshflx
         rmthflx = qw(i)-preqw(i)
         rmshflx = rmshflx+rmthflx*rmthflx
      enddo
c
      return
      end


      subroutine tecp
c     ************************************************************
c     *     tecplot
c     ************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nx=901,ny=901)
      common /aray1/ qc(0:nx,0:ny,4),qs(0:nx,0:ny,3),qv(0:nx,0:ny,1)
      common /aray2/ pc(0:nx,0:ny),cc(0:nx,0:ny),tt(0:nx,0:ny),
     &               tv(0:nx,0:ny),uc(0:nx,0:ny),vc(0:nx,0:ny),
     &               xc(0:nx,0:ny,5)
      common /coor1/ imax,jmax
      common /coor2/ x(nx,ny),y(nx,ny)
      common /metr1/ si(nx,ny,3),sj(nx,ny,3)
      common /metr2/ vol(0:nx,0:ny)
      common /cntl5/ nprc,iprc
      common /stat1/ loop,locl
      common /parm1/ nstp,nout
      common /chem1/ amass(5),amas1(5),cvi(5),h0i(3)
      common /chem3/ vmue(0:nx,0:ny),vkap(0:nx,0:ny),vkav(0:nx,0:ny)
      common /chem4/ ds(0:nx,0:ny,5),evs(0:nx,0:ny,5),hs(0:nx,0:ny,5)
      common /chem5/ fr,alp1,alp2,alp3,alp4,alp5,alp6,
     &                  bet1,bet2,bet3,bet4,bet5,bet6
      common /fscd2/ rinf,ruinf,rvinf,einf,evinf,r1inf,r2inf,r3inf,
     &               r4inf,r5inf,uinf,vinf,pinf,ttinf,tvinf,cinf,
     &               gaminf,xminf,x1inf,x2inf,x3inf,x4inf,x5inf
      common /wbcd1/ twall,cefn,cefo
      common /knud1/ aLocalKnud(0:nx,0:ny),frPth(0:nx,0:ny),
     &               thermVelo(0:nx,0:ny)
c
      dimension fv(0:nx,0:ny,5),fq(0:nx,0:ny)
      dimension flag(0:nx,0:ny),qn(0:nx,0:ny,20) ! qn 13 to 18
      dimension icom(20),pcom(20)
      character fname1*20
      character file_tec*50
c
      do j=0,jmax
         do i=0,imax
            flag(i,j)=0.d0
         end do
      end do
      do j=1,jmax-1
         do i=1,imax-1
            flag(i,j)=1.d0
         end do
      end do
c
      do j=1,jmax
         do i=1,imax
            volsum1  =1.d0/(vol(i-1,j-1)*flag(i-1,j-1)
     &                  +vol(i-1,j  )*flag(i-1,j  )
     &                  +vol(i  ,j-1)*flag(i  ,j-1)
     &                  +vol(i  ,j  )*flag(i  ,j  ))
            do k=1,4
               qn(i,j,k)=(qc(i-1,j-1,k)*vol(i-1,j-1)*flag(i-1,j-1)
     &              +qc(i-1,j  ,k)*vol(i-1,j  )*flag(i-1,j  )
     &              +qc(i  ,j-1,k)*vol(i  ,j-1)*flag(i  ,j-1)
     &              +qc(i  ,j  ,k)*vol(i  ,j  )*flag(i  ,j  ))
     &             *volsum1
            enddo
            qn(i,j,5)=(pc(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +pc(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &           +pc(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &             +pc(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
            qn(i,j,6)=(tt(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +tt(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &           +tt(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +tt(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
            qn(i,j,7)=(tv(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +tv(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &           +tv(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +tv(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
            qn(i,j,8)=(vmue(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +vmue(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &           +vmue(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +vmue(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
            qn(i,j,9)=(vkap(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +vkap(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &           +vkap(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +vkap(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
            qn(i,j,10)=(vkav(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +vkav(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &           +vkav(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +vkav(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
            qn(i,j,11)=(cc(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +cc(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &           +cc(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +cc(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
            qn(i,j,12)=(uc(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +uc(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &           +uc(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +uc(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
            qn(i,j,13)=(vc(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +vc(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &           +vc(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +vc(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
! output mass fraction
! O
            qn(i,j,14)=(xc(i-1,j-1,1)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +xc(i-1,j  ,1)*vol(i-1,j  )*flag(i-1,j  )
     &           +xc(i  ,j-1,1)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +xc(i  ,j  ,1)*vol(i  ,j  )*flag(i  ,j  ))*volsum1
!N
            qn(i,j,15)=(xc(i-1,j-1,2)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +xc(i-1,j  ,2)*vol(i-1,j  )*flag(i-1,j  )
     &           +xc(i  ,j-1,2)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +xc(i  ,j  ,2)*vol(i  ,j  )*flag(i  ,j  ))*volsum1
!NO
            qn(i,j,16)=(xc(i-1,j-1,3)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +xc(i-1,j  ,3)*vol(i-1,j  )*flag(i-1,j  )
     &           +xc(i  ,j-1,3)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +xc(i  ,j  ,3)*vol(i  ,j  )*flag(i  ,j  ))*volsum1
!O2
            qn(i,j,17)=(xc(i-1,j-1,4)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +xc(i-1,j  ,4)*vol(i-1,j  )*flag(i-1,j  )
     &           +xc(i  ,j-1,4)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +xc(i  ,j  ,4)*vol(i  ,j  )*flag(i  ,j  ))*volsum1
!N2
            qn(i,j,18)=(xc(i-1,j-1,5)*vol(i-1,j-1)*flag(i-1,j-1)
     &           +xc(i-1,j  ,5)*vol(i-1,j  )*flag(i-1,j  )
     &           +xc(i  ,j-1,5)*vol(i  ,j-1)*flag(i  ,j-1)
     &           +xc(i  ,j  ,5)*vol(i  ,j  )*flag(i  ,j  ))*volsum1
! Knudsen Number
! Knudsen
            qn(i,j,19)=(aLocalKnud(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &        +aLocalKnud(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &        +aLocalKnud(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &        +aLocalKnud(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
            qn(i,j,20)=(thermVelo(i-1,j-1)*vol(i-1,j-1)*flag(i-1,j-1)
     &        +thermVelo(i-1,j  )*vol(i-1,j  )*flag(i-1,j  )
     &        +thermVelo(i  ,j-1)*vol(i  ,j-1)*flag(i  ,j-1)
     &        +thermVelo(i  ,j  )*vol(i  ,j  )*flag(i  ,j  ))*volsum1
         enddo
      enddo
c
      do i=1,imax
         do j=1,jmax
            rho1       = 1.d0/qn(i,j,1)
            uu         = qn(i,j,12)
            vv         = qn(i,j,13)
            velo       = sqrt(vv*vv+uu*uu)
            fv(i,j,1)  = uu
            fv(i,j,2)  = vv
            fv(i,j,3)  = velo
            fv(i,j,4)  = velo/qn(i,j,11)
            fv(i,j,5)  = qn(i,j,11)
         enddo
      enddo
c
! 201  format('VARIABLES = "X","Y"'
!     &                  ,',"P","T","Tv","Rho"'
!     &                 ,',"U","V","Velo","Mach","a"')
 201  format('VARIABLES = "X","Y"'
     &                  ,',"P","T","Tv","Rho"'
     &                  ,',"U","V","Velo","Mach","a"'
     &                  ,',"O","N","NO","O2","N2"
     &                  ,"Knud","mu","thermVelo"')
 202  format('ZONE I=',i4,',J=',i4,',F=POINT')
c
      open(10,file='./output/prop.plt',form='formatted')
!      k = loop/nout
!      write (file_tec,'(a,i6.6,a)') "./output/prp",k,".plt"
!      open(10,file=file_tec,form='formatted')
      write(10,201)
      write(10,202) imax,jmax
      do j=1,jmax
         do i=1,imax
            write(10,*) x(i,j),y(i,j)
     &           ,qn(i,j,5),qn(i,j,6),qn(i,j,7),qn(i,j,1)
     &           ,(fv(i,j,isp),isp=1,5)
     &           ,qn(i,j,14),qn(i,j,15),qn(i,j,16)
     &           ,qn(i,j,17),qn(i,j,18)
     &           ,qn(i,j,19),qn(i,j,8),qn(i,j,20)
         enddo
         write(10,*)
      enddo
      close(10)
c
      open(11,file='./output/TEMP')
      open(12,file='./output/PDVM')
      open(13,file='./output/Mass')
      open(14,file='./output/Energy')
!      s   = 0.
      xcen = 0.5d0*(x(1,1)+x(2,1))
      ycen = 0.5d0*(y(1,1)+y(2,1))
      s = -sqrt(xcen**2+ycen**2)
      u0  = sqrt(uc(1,0)**2+vc(1,0)**2)
      u1  = sqrt(uc(1,1)**2+vc(1,1)**2)
      u   = 0.5d0*(u0+u1)
      p   = 0.5d0*(pc(1,0)+pc(1,1))
      c   = 0.5d0*(cc(1,0)+cc(1,1))
      ttw = 0.5d0*(tt(1,0)+tt(1,1))
      tvw = 0.5d0*(tv(1,0)+tv(1,1))
!      um  = u/cc(1,j)
      um  = u/cc(1,1)
      ro  = 0.5d0*(qc(1,0,1)+qc(1,1,1))
      e   = 0.5d0*(qc(1,0,4)+qc(1,1,4))
      ev  = 0.5d0*(qv(1,0,1)+qv(1,1,1))
      x1  = 0.5d0*(xc(1,0,1)+xc(1,1,1))
      x2  = 0.5d0*(xc(1,0,2)+xc(1,1,2))
      x3  = 0.5d0*(xc(1,0,3)+xc(1,1,3))
      x4  = 0.5d0*(xc(1,0,4)+xc(1,1,4))
      x5  = 0.5d0*(xc(1,0,5)+xc(1,1,5))
      write(11,101)"#","  s","T [K] ","Tv [K]"
      write(11,100)s,ttw,tvw
      write(12,101)"#","  s","P [MPa]","Rho [kg/m3]"
     &              ,"Velo [m/s]","  Mach"
      write(12,100)s,p*1.e-6,ro,u,um
      write(13,101)"#","  s","    O","    N",
     &             "   NO","   O2","   N2"
      write(13,100)s,x1,x2,x3,x4,x5
      write(14,101)"#","  s"," E [J/m3]"," Ev [J/m3]"
      write(14,100)s,e,ev
      do j=1,jmax-1
         xcen = 0.25*(x(1,j+1)+x(1,j)+x(2,j+1)+x(2,j))
         ycen = 0.25*(y(1,j+1)+y(1,j)+y(2,j+1)+y(2,j))
!         s = sqrt(xcen**2+ycen**2)/abs(x(1,1))-1.d0
!         s = sqrt(xcen**2+ycen**2)
         s = -sqrt(xcen**2+ycen**2)
         u = sqrt(uc(1,j)**2+vc(1,j)**2)
         um= u/cc(1,j)
         write(11,100)s,tt(1,j),tv(1,j)
         write(12,100)s,pc(1,j)*1.e-6,qc(1,j,1),u,um
         write(13,100)s,(xc(1,j,isp),isp=1,5)
         write(14,100)s,qc(1,j,4),qv(1,j,1)
      enddo
      close(11)
      close(12)
      close(13)
      close(14)
c
!      open(10,file='prop.gnu')
!      write(10,203) "#","  X","Y"
!     &             ,"P","T","Tv","Rho"
!     &             ,"U","V","Velo","Mach","a"
!      do j=1,jmax
!      do i=1,imax
!         write(10,204) x(i,j),y(i,j)
!     &   ,qn(i,j,5),qn(i,j,6),qn(i,j,7),qn(i,j,1)
!     &   ,(fv(i,j,isp),isp=1,5)
!      enddo
!      write(10,*)
!      enddo
!      do i=1,imax
!      do j=1,jmax
!         write(10,204) x(i,j),y(i,j)
!     &   ,qn(i,j,5),qn(i,j,6),qn(i,j,7),qn(i,j,1)
!     &   ,(fv(i,j,isp),isp=1,5)
!      enddo
!      write(10,*)
!      enddo
!      close(10)
c
 100  format(6e14.6)
 101  format(a,a13,5a14)
! 204  format(1p12e14.6)
 203  format(a,a13,11a14)
c
      return
      end

