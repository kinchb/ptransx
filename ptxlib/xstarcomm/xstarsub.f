      subroutine xstarsub(x_energies, x_mi, x_intmi, density, Fe_abund,
     $                    heat, temp, niter, new_xee, new_temp,
     $                    x_c_emis, x_absorp, x_elines, x_linemis,
     $                    heat_vals)
c
      implicit none
c
      save
c
      include './PARAM'
c
      real*8 x_energies(nxx)
      real*8 x_mi(nxx)
      real*8 x_intmi(nxx)
      real*8 density
      real*8 Fe_abund
      real*8 heat
      real*8 temp
      integer niter
c
      real*8 new_xee
      real*8 new_temp
      real*8 x_c_emis(nxx)
      real*8 x_absorp(nxx)
      real*8 x_elines(nnnl)
      real*8 x_linemis(nnnl)
      real*8 heat_vals(9)
c
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character*1 kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     line luminosities
      real*8 elum(3,nnnl),elumo(3,nnnl)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
      real*8 fline(2,nnnl),flinel(ncn)
c     line optical depths
      real*8 tau0(2,nnnl)
c     energy bins
      real*8 epi(ncn)
c      continuum lum
      real*8 zrems(4,ncn),zremso(4,ncn),
     $          zremsz(ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn),brcems(ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     level populations
      real*8 xilev(nnml),bilev(nnml),rniss(nnml)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
      real*8 elumab(2,nnml),elumabo(2,nnml)
      real*8 tauc(2,nnml)
c     ion abundances
      real*8 xii(nni)
c     heating and cooling
      real*8 htt(nni),cll(nni)
      real*8 rrrt(nni),pirt(nni)
      integer nlevs(nni)
c     element abundances
      real*8 abel(nl),ababs(nl)
c     the saved rates
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
c     compton heating data
      real*8 decomp(ncomp,ncomp),ecomp(ncomp),sxcomp(ncomp)
      character(63) atcredate
c
c     local variables
c     state variables
      real*8 p,rdel,r19,xi,xcol,zeta,rdelo,r,t,xpx,delr,xpx0,r0
c     heating-cooling variables
      real*8 httot,cltot,htcomp,clcomp,clbrems,elcter,cllines,
     $     clcont,hmctot,fh,fe
c     input parameters
      character*16 knam,knam2
      character*80 kmodelname,specfile,spectype
      real*8 enlum,emult,taumax,xeemin,xlum,rmax,xpxcol,trad
      real*8 cfrac,critf,vturbi,xlum2,xee,radexp
      integer lcdd,ncn2
c     variables associated with thermal equilibrium solution
      integer nmat,ntotit,lnerrd
c     switches
      integer lprid,lpril,lpri,lwri,lpriu,nlimdt,lpris,lfpi
      integer lprisv,lforce,lpri2
      integer  nnmax,nlimd,lunlog,nloopctl,numrec,npass,lfix
c     temporary for xwrite
      character*133 tmpst
      real*8 ectt
c     temporary for spectrum
      real*8 eptmp(ncn),zrtmp(ncn)
      integer nlprnt(18),nlnprnt
c     temporary integers
      integer ll,mm,kk,ldir,jk,jkp,jk2
      integer nlsvn,ncsvn
      real*8 eliml,elimh
      integer istatus, iunit,iunit2,numcon2,ncut,ncutmx,iunit3
c     times
      real*8 tinf,ttot,t1s,t2s
      integer luna,lunb,lun30,lun31,ntmp,np2
      integer lpril2
C     storing info for parameters
      character*20 parname(55)
      character*10 partype(55)
      real*8 parms(55)
      character*30 parcomm(55)
      integer nparms, specunit,irecl
c     temporary line pointers
      integer nlin(nnnl)
      real*8 elin(nnnl),rctmp(2,nnnl),fpr2
      real*8 cemtmp(2,nnml)
c
      real*8 ct1,ct2,time
c
c     -----------------------------------------------------
c     variables we need for our PTRANSX interface
c
      real*8 cemis(ncn)
      real*8 scatt(ncn)
      real*8 tmpsum
      integer ii, jj
      character*50 infile, outfile, abundfile
      real*8 t_min
      real*8 mhd_heat
      real*8 clsup
      integer tmp

      integer test
      integer test2
c
c     print *, test
      test = test + 1
c
c     print *, test2
c
c     -----------------------------------------------------
c
c     local definitions
c
C     Parameter Names
C
      data parname/'cfrac','temperature',
     $   'lcpres','pressure','density','spectrum',
     $   'spectrum_file','spectun','trad',
     $   'rlrad38','column','rlogxi',
     $   'nsteps','niter','lwrite',
     $   'lprint','lstep',
     $   'habund','heabund',
     $   'liabund','beabund','babund','cabund',
     $   'nabund','oabund','fabund','neabund',
     $   'naabund','mgabund','alabund',
     $   'siabund','pabund','sabund',
     $   'clabund','arabund','kabund',
     $   'caabund','scabund','tiabund',
     $   'vabund','crabund','mnabund ',
     $   'feabund','coabund','niabund',
     $   'cuabund','znabund','emult','taumax','xeemin',
     $   'critf','vturbi','npass','modelname',
     $   'loopcontrol'/
      data partype/'real','real',
     $    'integer','real','real','string',
     $    'string','integer','real',
     $    'real','real','real',
     $    'integer','integer','integer',
     $    'integer','integer',
     $    'real','real','real',
     $    'real','real','real',
     $    'real','real','real',
     $    'real','real','real',
     $    'real','real',
     $    'real','real','real',
     $    'real','real','real',
     $    'real','real','real',
     $    'real','real','real',
     $    'real','real','real',
     $    'real','real','real',
     $    'real','real','real','integer','string',
     $    'integer'/
      data parcomm/' ','Units of 10**4 K',
     $     '1=yes, 0=no','dynes/cm**2','cm**(-3)',' ',
     $     ' ','0=energy, 1=photons','or alpha',
     $     '/10**38 erg/sec','cm**(-2)',' ',
     $     ' ',' ','1=yes, 0=no',
     $     '1=yes, 0=no',' ',
     $     ' ',' ',' ',
     $     ' ',' ',' ',
     $     ' ',' ',' ',
     $     ' ',' ',' ',
     $     ' ',' ',' ',
     $     ' ',' ',
     $     ' ',' ',' ',
     $     ' ',' ',' ',
     $     ' ',' ',' ',
     $     ' ',' ',' ',
     $     ' ',' ',' ',
     $     ' ',' ',' ',
     $     ' ',' ',
     $     ' '/
      data nlnprnt/14/,nlprnt/11,22,1,19,23,24,5,10,16,14,4,6,15,
     $            21,7,3*0/
      nparms=55
c
      call remtms(ct1)
c
c     -----------------------------------------------------
c     get command-line arguments to determine what file to
c     read in and what file to write out to
c
c     -----------------------------------------------------
c     if we are going to modify epi, it must be done here,
c     and we must set ncn2 appropriately
c
      ncn2 = nxx
      do ii = 1, ncn2
          epi(ii) = x_energies(ii)
      enddo
c
c     -----------------------------------------------------
c     basically, this section exists so that we don't
c     accidentally break xstarsetup
c
      lpri=0
      lunlog=6
      nlimd=99
      numrec=3
      lwri=0
      lfix=0
      npass=1
      t=100.
      r=0.
      xpxcol=1.e+10
      xpx=1.e+10
      lcdd=1
      zeta=2.
      emult=0.5
      taumax=5.
      xeemin=0.01
      trad=-1.
      xlum=1.
      cfrac=0.
      emult=0.75
      taumax=5.
      xeemin=0.1
      nmat=nds
      nloopctl=0
      critf=1.d-8
      vturbi=1.
      xcol=0.
      radexp=0.
      do mm=1,ncn
        zremsz(mm)=0.
        enddo
      lpris=lpri
      lpri=0
      spectype='pow     '
      xlum2=1.
      if (spectype.eq.'bbody   ')
     $   call starf(trad,xlum2,epi,ncn2,zremsz,lpri,lunlog)
      if (spectype.eq.'pow     ')
     $   call ispec4(trad,xlum,epi,ncn2,zremsz,lpri,lunlog)
      if (spectype.eq.'brems   ')
     $   call ispec(trad,xlum,epi,ncn2,zremsz,lpri,lunlog)
      call ispecgg(xlum,epi,ncn2,zremsz,
     $               lpri,lunlog)
      call ispcg2(zremsz,epi,ncn2,enlum,lpri,lunlog)
      xee=1.21
      abel(1)=1.
      do mm=3,30
          abel(mm)=1.
       enddo
      abel(2)=1.
      xi=10.**zeta
      r19=sqrt(xlum/amax1(1.e-36,xi*xpx))
      r=r19*1.e+19
c
c     -----------------------------------------------------
c     set abundances here:
c
c     *all solar*
      do ii = 1, nl
          if (Fe_abund.lt.0.0) then
              abel(ii) = 0.0
          else
              abel(ii) = 1.0
          endif
      enddo
      abel(1) = 1.0
      abel(2) = 1.0
c
c     except Fe
      abel(26) = abs(Fe_abund)
c
      if (Fe_abund.gt.0.0) then
          test2 = test2 + 1
      endif
c
c     H/He/Fe only*
c     do ii = 1, nl
c         abel(ii) = 0.0
c     enddo
c     abel(1)  = 1.0
c     abel(2)  = 1.0
c     abel(26) = 1.0
c
c     xee = 1.0
c
c     -----------------------------------------------------
c     we do not want to modify this call to xstarsetup
c
      if (test.eq.1) then
          call xstarsetup(lnerrd,nlimd,
     $       lpri,lprid,lunlog,tinf,critf,
     $       t,trad,r,delr,xee,xpx,ababs,abel,cfrac,xlum,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,atcredate,
     $       decomp,ecomp,sxcomp,
     $       zrems,zremsz,
     $       tau0,dpthc,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,nlevs,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,ncsvn,rates,vsav,idrates,
     $       ntotit,
     $       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,
     $       xilev,bilev,rniss,nmat,elum,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,nlin,elin)
      elseif (test2.eq.1) then
          call xstarsetup(lnerrd,nlimd,
     $       lpri,lprid,lunlog,tinf,critf,
     $       t,trad,r,delr,xee,xpx,ababs,abel,cfrac,xlum,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,atcredate,
     $       decomp,ecomp,sxcomp,
     $       zrems,zremsz,
     $       tau0,dpthc,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,nlevs,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,ncsvn,rates,vsav,idrates,
     $       ntotit,
     $       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,
     $       xilev,bilev,rniss,nmat,elum,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,nlin,elin)
      endif
c
c    ------------------------------------------------------
c    read in some data, then call xstarcalc
c
      xpx = density
      mhd_heat = heat
      t = temp
      t_min = 0.0
      nlimdt = niter
c     nlimdt = 0
      do ii = 1, ncn2
          bremsa(ii) = x_mi(ii)
      enddo
      do ii = 1, ncn2
          bremsint(ii) = x_intmi(ii)
      enddo
c
      vturbi = 0.0
      critf  = 1.d-8
      clsup  = 0.0
c
c     print *, "got here"
c     print *, niter
c     print *, 'in xstar: t = ', t
c
          call xstarcalc(lpri2,lnerrd,nlimdt,
     $       lpri,lprid,lunlog,tinf,vturbi,critf,
     $       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       decomp,ecomp,sxcomp,
     $       tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,nlevs,ncsvn,rates,vsav,idrates,
     $       ntotit,
     $       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,
     $         cllines,clcont,htcomp,clcomp,clbrems,
     $       xilev,bilev,rniss,nmat,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,fline,flinel, t_min, mhd_heat, clsup)
c
c     -----------------------------------------------------
c
c     print *, "got here"
c
      new_xee = xee
      new_temp = t
c
c     print *, 'in xstar: new_temp = ', new_temp
c     print *, 'in xstar: new_xee = ', new_xee
c
      if (new_xee.lt.1.0) then
          print *, 'in xstar: xee is small; supplied T is ', temp 
      endif
c
c     multiply by 4pi (a relic from XSTAR's past)
c     and then copy rccemis(2) into cemis for later on
      rccemis(2, 1) = (1.0 / (4.0 * 3.14159)) * brcems(1)
      rccemis(2, 2) = (1.0 / (4.0 * 3.14159)) * brcems(2)
      rccemis(2, 3) = (1.0 / (4.0 * 3.14159)) * brcems(3)
      cemis(1)      = rccemis(2, 1)
      cemis(2)      = rccemis(2, 2)
      cemis(3)      = rccemis(2, 3)
      do ii = 4, ncn2
c         if (isnan(rccemis(2, ii))) then
c             print *, 'in xstar: nan in rccemis at ', ii
c         endif
c         if (isnan(brcems(ii))) then
c             print *, 'in xstar: nan in brcems at ', ii
c         endif
          rccemis(2, ii) = (2.0 * rccemis(2, ii))
     $    + ((1.0 / (4.0 * 3.14159)) * brcems(ii))
          cemis(ii) = rccemis(2, ii)
      enddo
c
      do ii = 1, ncn2
          x_c_emis(ii) = cemis(ii)
          if (isnan(x_c_emis(ii))) then
c             print *, 'in xstar: nan in x_c_emis at ', ii
              if (ii.eq.1) then
                  x_c_emis(ii) = x_c_emis(ii+1)
              elseif (ii.eq.ncn2) then
                  x_c_emis(ii) = x_c_emis(ii-1)
              else
                  x_c_emis(ii) = 0.5 * (x_c_emis(ii-1) + x_c_emis(ii+1))
              endif
          endif
      enddo
c
      do ii = 1, ncn2
          x_absorp(ii) = opakc(ii)
c         if (isnan(x_absorp(ii))) then
c             print *, 'in xstar: nan in x_absorp at ', ii
c         endif
      enddo
c
      do ii = 1, nnnl
          x_elines(ii) = 12398.4/elin(ii)
      enddo
c
      do ii = 1, nnnl
          x_linemis(ii) = rcem(2, ii)
      enddo
c
c     now bin the line emissivity into the continuum emissivity
c
c     do ii = 1, ncn2-1
c         tmpsum = 0.0
c         do jj = 1, nnnl
c             if ((12398.4/elin(jj) >= epi(ii)).and.
c    $        (12398.4/elin(jj) < epi(ii+1))) then
c                 tmpsum = tmpsum + rcem(2, jj)
c             endif
c         enddo
c         if (tmpsum.gt.0.0) then
c             tmpsum =
c    $        tmpsum/(1.602177e-12 * (epi(ii+1) - epi(ii)))
c             rccemis(2, ii) = rccemis(2, ii)
c    $        + (tmpsum / (2.0 * 3.14159))
c         endif
c     enddo
c
c     write output here ... first line is xee, then the
c     left column is emissivity and the right is opacity
c
      heat_vals(1) = httot
      heat_vals(2) = htcomp
      heat_vals(3) = httot-htcomp-mhd_heat
      heat_vals(4) = cltot
      heat_vals(5) = clcomp
      heat_vals(6) = clbrems
      heat_vals(7) = cllines
      heat_vals(8) = clcont
      heat_vals(9) = clsup
c
      tmp = -1
      call getlun(tmp)
c
c     open(unit = 54, file = "xout_abund_0.dat")
c     write(54, *) nni
c     do ii = 1, nni
c         write(54, *) xii(ii)
c     enddo
c     close(54)
c
c     -----------------------------------------------------
c
      return
      end
c
c     -----------------------------------------------------
c     -           SUBROUTINE DEFINITIONS BELOW            -
c     -----------------------------------------------------
c
      subroutine amcol(n,l,temp,ic,z1,rm,ne,sum,ecm,cn)
c
c subroutine amcol determines the rate of angular momentum changing
c collisions in hydrogenic atoms due to collisions with ions.
c the codes are based on the method of hummer & storey (1987)
c      author: M. Bautista
c
c        z1=charge of colliding ion
c        ic=ionic charge of hydrogenic atom
c        rm=reduced mass of colliding system in units of electron mass
c        ne=electron number density
c        sum = sum of spontaneous transitions out of level n,l
c        cn = transition rate for nl -> nl-1
c        cn(nl -> nl-1) = cn
c        cn(nl -> nl+1) = cn*(2.*l+1)/(2.*l-1)
c
      implicit none
c
      real*8 ne,temp,z1,rm,sum,ecm,cn
      real*8 pc1,pc2,pc3,dnl,pc,qnl,den
      integer n,l,ic
c
      pc1=1.181+log10(temp/ne)
      pc2=pc1
      pc3=pc1
      dnl=6.*z1/ic*z1/ic*n*n*(n*n-l*l-l-1)
      if(sum.ne.0.) pc2=10.95+log10(temp/rm/sum/sum)
      if(ecm.ne.0.) pc3=log10(temp/rm/ecm/ecm)-11.22
      pc=min(pc1,pc2,pc3)
      qnl=9.933e-6*sqrt(rm/temp)*dnl*(11.538+log10(temp/dnl/rm)+pc)
c
      den=l*(n*n-l*l)+(l+1)*(n*n-(l+1)*(l+1))
      cn=qnl*l*(n*n-l*l)/den
c
      return
      end
      subroutine amcrs(n,l,temp,ic,z1,rm,ne,sum,ecm,psi,il,cn,
     $                 lpri,lun11)
c
c the angular momentum changing collision rates are calculated using
c either the pengelly & seaton (1964) formula (amcol) or the impact
c parameter method of seaton (1962) (impact) if the energy levels are
c non-degenerate.  the ps routine is used if the ratio of amcol/impact
c is greater than 0.94 since this code is faster.  ** beware - there
c may be problems if ne is too large ( > 1.e+7).  pc1 will be used in
c amcol rather than pc3 and the change will not occur.
c      author: M. Bautista
c
c     n = principal quantum number of initial state
c     l = orbital quantum number of initial state
c     temp = temperature in kelvin
c     ic = ionic charge of target particle
c     z1 = charge of incident particle
c     rm = mass of incident particle in units of electron mass me
c     ne = electron number density
c     sum = total spontaneous transition rate out of n,l
c     cn = transition rate for nl -> nl-1
c     ecm = energy difference between nl and nl-1
c     psi = see notes for defn
c     il = flag to decide whether to use impact or amcol
c     cn = transition rate for nl -> nl-1
c
      implicit none
c
      real*8 ne
      real*8 temp,z1,rm,sum,ecm,psi,cr
      real*8 en,dnl,rho1,rhom,cn,rat
      integer n,l,ic,il,lun11,lpri
c
      en=real(n)
      dnl=6.*z1/ic*z1/ic*n*n*(n*n-l*l-l-1)
      rho1=0.72/sum
      if(ecm.ne.0.) rho1=min(rho1,5.946e-12/ecm)
      rhom=3.929e11*rho1*temp/sqrt(dnl)/rm
      if(rhom.lt.10.) go to 30
      call amcol(n,l,temp,ic,z1,rm,ne,sum,ecm,cn)
      cn=0
c mab il=0
        il=0
        if(ecm.ne.0.) then
          if(il.eq.0) then
          if (lpri.gt.1)
     $    write (lun11,*)'call impact 1',en,l,temp,ic,z1,rm,ecm,psi
          call impact(en,l,temp,ic,z1,rm,ecm,psi,cr)
          rat=cn/cr
          cn=cr
          if(rat.gt. 0.94) il=1
          endif
        endif
c     go to 40
c
 30    if(ecm.eq.0.) then
      call velimp(n,l,temp,ic,z1,rm,ne,sum,cn)
      else
      if (lpri.gt.1)
     $write (lun11,*)'call impact 2',en,l,temp,ic,z1,rm,ecm,psi
      call impact(en,l,temp,ic,z1,rm,ecm,psi,cn)
      endif

c      if(ne.gt.1.e14) then
c     call impact(en,l,temp,ic,z1,rm,ecm,psi,cn)
c      endif
c
c 40    continue
c
      return
      end
      subroutine anl1(ni,nf,lf,iq,alm,alp,lpri,lun11)
c
c this subroutine is used to calculate the values of the
c spontaneous transition rates for electric dipole transitions from
c level ni,lf+1 and ni,lf-1 to level nf,lf.
c the transition probabilities (a values) are calculated
c using the gordon (1929) formula.
c        iq=ionic charge
c      author: M. Bautista
c
      implicit none
c
      real*8 y1,y2,x1,x2,x3,x4,x5,t
      real*8 alm,alp,an,dum
      integer ni,nf,lf,iq,lpri,lun11
      integer li,n,l,np,ia1,ia2,ib,ic
      real*8 x,rev,rn
c
      alm=0.
c
c **** for case a set lower limit of nf=1, for case b nf starts at 2
c
      do 40 li=lf-1,lf+1,2
      if(li.lt.0) go to 40
      if(lf.gt.li) go to 100
      n=ni
      np=nf
      l=li
      go to 101
 100   n=nf
      np=ni
      l=lf
 101   continue
c
       call dfact(n+l,x1)
       call dfact(np+l-1,x2)
       call dfact(2*l-1,x3)
       call dfact(n-l-1,x4)
       call dfact(np-l,x5)
      ia1=-n+l+1
      ia2=ia1-2
      ib=-np+l
      ic=2*l
      x=-4.*n*np/((n-np)*(n-np))
       call hgf(ia1,ib,ic,x,y1)
       call hgf(ia2,ib,ic,x,y2)
      rev=abs(n-np)
      rn=float(n+np)
      t=(l+1)*log((4.*float(n*np)))+(rn-2*l-2)*log(rev)
      t=t-log(4.e0)-rn*log(rn)
      y1=abs((y1-y2*(rev/rn)**2))
      y1=log(y1)+t
      t=2.*y1+x1+x2-2.*x3-x4-x5
      t=exp(t)
      an=2.6761e09*iq**4*max(li,lf)*t/(2.*li+1)
         dum=(1./nf/nf-1./ni/ni)**3.
      an=dum*an
      if(li.lt.lf) alm=an
      if(li.gt.lf) alp=an
c
 40    continue
c
      if (lpri.gt.1) then
        write (lun11,*)'in anl1:',li,lf,t,ni,nf,iq
        write (lun11,*) rn,n,np,rev,y1,y2,x
        write (lun11,*) ia1,ia2,ib,ic,x1,x2,x3,x4,x5,l,an
        endif
c
      return
      end
       subroutine binemis(lun11,lpri,xlum,
     $       t,vturbi,epi,ncn2,dpthc,
     $       idat1,rdat1,kdat1,nptrs,
     $       npar,npnxt,npfi,
     $       nplin,nlsvn,
     $       eliml,elimh,elum,zrems,zremsz,ilsv,
     $       zrtmp,ewsv,elsv,nlsv)
c
C
      implicit none
c
      include './PARAM'
      integer nbtpp
      parameter (nbtpp=ncn)
c
c
c     passed parameters
      integer lun11
c     master data
      integer idat1(nidat1),nptrs(nptt,ndat2)
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl)
c     line luminosities
      real*8 elum(3,nnnl)
c     energy bins
      real*8 epi(ncn)
c     continuum lum
      real*8 zrems(4,ncn),zremsz(ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
      real*8 zrtmp(6,ncn)
      integer kl, nilin, nkdt,nidt,lcon,lrtyp,ltyp,ml
      integer nlsvn, ln, ll, numcon, lnn, nbtmp, nrdt
      integer  verbose
      integer lup,ndtmp,mllz,nkdt2,iion,nitmp,
     $        ltyp2,lrtyp2,lcon2,nrdt2,nidt2,iltmp,mlpar
      real*8 eliml, elimh, elmmtpp, dele, etmp, elin, aasmall
      real*8 dele2,egam,profile,deler,delea,rdat2(20000),vturbi,aatmp
      real*8 optpp,e00,deleepi,etptst,tst,sume,zrsum1,zrsum2,deletpp,
     $        deleused,tmpe,zrtp2,zrtp1,bbb,deleused2,xlum
      integer ml1,mlmin,mlmax,ij,mlm,ldir,mloff,ml1m,ml1mn,mlc,ncut,
     $        ml1mx,ml2,np1k,np1i,np1r,np1i2,ml1max,ml1min
c     arrays containing printout info
      integer ilsv(nnnl)
      real*8 ewsv(nnnl),elsv(nnnl)
      real*8 etpp(nbtpp),zrtpp2(2,nbtpp),zrtmps(2,nbtpp)
      integer ldon(2)
c
      integer lpri,ncn2,nlsv,nelin
      real*8 t,delet,deleturb,deleth,e0,vth,vturb,dpcrit
      real*8 zro
c
c     externally defined functions
      integer nbinc
      real*8 voigte
c
      data dpcrit/1.e-6/
c
      verbose=lpri
c
c     open and prepare the fits file for spectral data
      if(verbose.gt.0) write (lun11,*)'in binemis:'

c     build spectra data tables
      numcon=ncn2
      bbb=vturbi
      do ll=1,ncn2
        zrtmp(4,ll)=0.
        zrtmp(5,ll)=0.
        zrtmps(1,ll)=0.
        zrtmps(2,ll)=0.
        zrtpp2(1,ll)=0.
        zrtpp2(2,ll)=0.
        enddo
      do  lnn=1,nlsvn
        ln=lnn
        ml=nplin(ln)
        call drd(ltyp,lrtyp,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,ml-1,
     $          nptrs,0,lun11)
        elin=abs(rdat1(np1r))
        egam=rdat1(np1r+2)
        lup=idat1(np1i+1)
        nilin=npar(ml)
        call drd(ltyp,lrtyp,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,nilin-1,
     $          nptrs,0,lun11)
        nelin=npar(nilin)
        nilin=idat1(np1i+2)
c       get nuclear mass
        ml=nelin
        call drd(ltyp,lrtyp,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,ml-1,
     $          nptrs,0,lun11)
        aatmp=rdat1(np1r+1)
        elmmtpp=(elum(2,ln)+elum(1,ln))/2.
        if (lpri.ne.0)
     $     write (lun11,*)ln,elin,elmmtpp,nilin,nelin,egam,
     $              lup,aatmp
        if (((ln.gt.0).and.(ln.le.nnnl))
     $    .and.((elin.gt.eliml).and.(elin.lt.elimh))
     $    .and.(elmmtpp.gt.1.e-12*xlum).and.(aatmp.gt.1.e-24)
     $    .and.((nilin.gt.0).and.(nilin.le.nni)))
     $       then

c         line parameters
          etmp=12398.54/elin
          nbtmp=nbinc(etmp,epi,ncn2)
c
c         find associated type 86 data
          iion=1
          nitmp=npfi(13,iion)
          call drd(ltyp,lrtyp,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,nitmp-1,
     $          nptrs,0,lun11)
          if (lpri.ne.0)
     $      write (lun11,*)'searching for ion'
          do while ((idat1(np1i-1+nidt).ne.nilin).and.(iion.lt.nni))
            iion=iion+1
            nitmp=npfi(13,iion)
            call drd(ltyp,lrtyp,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,nitmp-1,
     $          nptrs,0,lun11)
            if (lpri.ne.0)
     $        write (lun11,*)iion,idat1(np1i-1+nidt),nilin,nitmp
            enddo
          ndtmp=npfi(41,iion)
          delea=0.
          if (ndtmp.gt.0) then
            if (lpri.ne.0)
     $        write (lun11,*)'  found ion',lup,ndtmp
            mllz=npar(ndtmp)
            call drd(ltyp2,lrtyp2,lcon2,
     $          nrdt,np1r,nidt,np1i2,nkdt,np1k,ndtmp-1,
     $         nptrs,0,lun11)
            iltmp=idat1(np1i2+1)
            mlpar=mllz
            do while ((ndtmp.ne.0).and.(lup.ne.iltmp)
     $         .and.(mlpar.eq.mllz))
               call drd(ltyp2,lrtyp2,lcon2,
     $           nrdt,np1r,nidt,np1i2,nkdt,np1k,ndtmp-1,
     $          nptrs,0,lun11)
              iltmp=idat1(np1i2+1)
              if (lpri.ne.0)
     $           write (lun11,*)'   ',nidt2,iltmp,ndtmp
              ndtmp=npnxt(ndtmp)
              mlpar=0
              if (ndtmp.ne.0) mlpar=npar(ndtmp)
              enddo
            endif
          if (lup.eq.iltmp) then
            delea=rdat2(3)*(4.14e-15)
            egam=rdat2(4)
            endif
c
c         cheat for narrow line plot
c         delea=0.
c
c         a list of all the deles
c           delea=auger natural width in eV
c           deleturb=turbulent width
c           deleth=thermal Doppler width
c           dele=thermal+turbulent width
c           deler=radiative natural width
c           dele2=total width, natural + Doppler
c           deletpp=goal of resolution of internal grid=dele/8
c           deleepi=xstar grid spacing
c           deleused=spacing of internal grid=deleepi/int(deleepi/depetpp)
c           delet=energy offset from line center in units of dele (local)
c
c         thermal width quantities
          vth=(1.2e+1)*sqrt(t/aatmp)
          vturb=max(bbb,vth)
          e0=(12398.54)/max(elin,1.e-24)
          deleturb=e0*(vturb/3.e+5)
          deleth=e0*(vth/3.e+5)
c         old expression
c          dele=deleth+deleturb
c         new expression
          dele=sqrt(deleth*deleth+deleturb*deleturb)
          deler=egam*(4.14e-15)
          dele2=delea+deler+dele
          aasmall=(delea+deler)/(1.e-36+dele)/12.56
c
          ml1=nbtmp
          if (lpri.ge.1) write (lun11,*)
     &   'e0,elin,elum1,elum2,ml1,deleth,delea:',
     &    e0,elin,elum(1,ln),elum(2,ln),ml1,deleth,delea
c
c         calculate profile on temporary grid
          e00=epi(ml1)
          etmp=e0
c         deleepi is the grid spacing of the epi grid
c         deletpp is the physical energy spacing needed
c           for an accurate integration of the voigt profile
c         ncut is the ratio of these two quantities,
c           used for rebinning the calculated voigt profile
          deleepi=epi(ml1+1)-epi(ml1)
          deletpp=dele
          ncut=int(deleepi/deletpp)
          ncut=max(ncut,1)
          ncut=min(ncut,nbtpp/10)
          deleused=deleepi/float(ncut)
          mlc=0
          ldir=1
          ldon(1)=0
          ldon(2)=0
          mlmin=nbtpp
          mlmax=1
          ml1min=ncn+1
          ml1max=0
          ml2=nbtpp/2
          if (lpri.ne.0) write (lun11,*)'ncut=',ncut,deleused,deletpp,
     $                                    deleepi
c
c         calculate profile at continuum bin closest to line center
          delet=(e00-etmp)/dele
          if (aasmall.gt.1.e-9) then
              profile=voigte(abs(delet),aasmall)/1.772
            else
              profile=exp(-delet*delet)/1.772
            endif
          profile=profile/dele/(1.602197e-12)
          etpp(ml2)=e00
          zrtpp2(1,ml2)=elum(1,ln)*profile
          zrtpp2(2,ml2)=elum(2,ln)*profile
          tst=1.
c
c         now put profile on temporary grid
c         work outward in both directions from line center
          do while ((ldon(1)*ldon(2).eq.0).and.(mlc.lt.nbtpp/2))
c
            mlc=mlc+1
c
c           alternate directions
            do ij=1,2
              ldir=-ldir
c
c             test to see if done in this direction
              if (ldon(ij).ne.1) then
c
c               index into temporary grid
                mlm=ml2+ldir*mlc

                etptst=e00+float(ldir*mlc)*deleused
c
c               test to see if within allowed range
                if ((mlm.lt.nbtpp).and.(mlm.gt.1)
     $            .and.(etptst.gt.0.).and.(etptst.lt.epi(ncn2))) then
c
c                 calculate index extremes for later use
c                 ml1m is index into epi grid
c                 ml1min and ml1max are extremes of ml1m
c                 mlmin and mlmax are extremes of mlm
                  mlmin=min(mlm,mlmin)
                  mlmax=max(mlm,mlmax)
c
c                 store energy binc
                  etpp(mlm)=e00+float(ldir*mlc)*deleused
c
c                 calculate profile
                  delet=(etpp(mlm)-etmp)/dele
                  if (aasmall.gt.1.e-9) then
                      profile=voigte(abs(delet),aasmall)/1.772
                    else
                      profile=exp(-delet*delet)/1.772
                    endif
                  profile=profile/dele/(1.602197e-12)
c
                  zrtpp2(1,mlm)=elum(1,ln)*profile
                  zrtpp2(2,mlm)=elum(2,ln)*profile
                  tst=profile
c
c                 print
                  if (lpri.ne.0) write (lun11,*) 'first write',
     $               mlm,etpp(mlm),ij,
     $               deleused,delet,mlmin,mlmax,ml1,
     $               mlc,mloff,mod(mloff,ncut),profile,zrtpp2(2,mlm)
c
c                 end of test for within range
                  endif
c
c               test to see if done in this direction:
c                 profile not too small
c                 index within range
c                 energy within range
c                 within specified number of doppler widths (50)
                if (((tst.lt.dpcrit)
     $               .or.(mlm.le.1).or.(mlm.ge.nbtpp)
     $               .or.(etptst.le.0.).or.(etptst.ge.epi(ncn2))
     $               .or.(mlc.gt.nbtpp)
     $               .or.(abs(delet).gt.max(50.,200.*aasmall)))
     $               .and.(ml1min.lt.ml1-2).and.(ml1max.gt.ml1+2)
     $               .and.(ml1min.ge.1).and.(ml1max.le.ncn))
     $                ldon(ij)=1
c
c               end of test for done in this direction
                endif
c
c             end of loop over directions
              enddo
c
c           end of loop over energies
            enddo
c
c         store into continuum bins
          sume=0.
          zrsum1=0.
          zrsum2=0.
          ml1min=nbinc(etpp(mlmin),epi,ncn2)
          ml1max=nbinc(etpp(mlmax),epi,ncn2)
          ml1m=ml1min
          if (lpri.ne.0) write (lun11,*)'renormalizing profile',
     $       ml2,mlmin,mlmax,ml1,ml1min,ml1max
          mlmin=max(mlmin,2)
          mlmax=min(mlmax,nbtpp)
c          ml1m=nbinc(etpp(mlmin),epi,ncn2)
          ml1mx=1
          ml1mn=ncn2
c
c         step through temp grid bins
c         and  sum over intervals
          do mlm=mlmin+1,mlmax
c
            tmpe=abs(etpp(mlm)-etpp(mlm-1))
            sume=sume+tmpe
            zrsum1=zrsum1+(zrtpp2(1,mlm)+zrtpp2(1,mlm-1))*tmpe/2.
            zrsum2=zrsum2+(zrtpp2(2,mlm)+zrtpp2(2,mlm-1))*tmpe/2.
c            if (lpri.ne.0) write (lun11,*)mlm,etpp(mlm),ml1m,epi(ml1m),
c     $         sume,zrsum1,zrsum2
c
c           test to see if you have reached epi grid boundary
            if (etpp(mlm).gt.epi(ml1m)) then
c
c             store current sum
              if (mlm.eq.mlmax) ml1m=max(1,ml1m-1)
              if (sume.gt.1.d-24) then
                zrtp2=zrsum2/sume
                zrtp1=zrsum1/sume
                do while ((etpp(mlm).gt.epi(ml1m)).and.(ml1m.lt.ncn2))
                  zrtmps(1,ml1m)=zrtp1
                  zrtmps(2,ml1m)=zrtp2
                  if (lpri.ne.0) write (lun11,*)mlm,ml1m,
     $               epi(ml1m),epi(ml1m+1),etpp(mlm),zrtmps(2,ml1m)
                  ml1m=ml1m+1
                  enddo
                endif
c
c             reset interval sums
              zrsum2=0.
              zrsum1=0.
              sume=0.
c
c             end of test for epi bin boundary
              endif
c
c           end of rebinning loop
            enddo
c
          do ml1m=ml1min,ml1max
            zrtmp(5,ml1m)=zrtmp(5,ml1m)+zrtmps(2,ml1m)
            zrtmp(4,ml1m)=zrtmp(4,ml1m)+zrtmps(1,ml1m)
            if (lpri.ne.0) write (lun11,*)ml1m,
     $           epi(ml1m),zrtmps(1,ml1m),zrtmp(4,ml1m)
            enddo
          endif
        enddo


      if (lpri.ne.0) write (lun11,*)'after first binemis loop'
      do kl=1,numcon
         if (lpri.ne.0) write (lun11,*)kl,epi(kl),zrems(2,kl),
     $          zrems(3,kl),zrtmp(4,kl),zrtmp(5,kl),dpthc(1,kl),
     $          zremsz(kl)
         zrtmp(4,kl)=zrtmp(4,kl)+zrems(2,kl)

         zrtmp(5,kl)=zrtmp(5,kl)+zrems(3,kl)
         zrems(2,kl)=zrtmp(4,kl)
         zrems(3,kl)=zrtmp(5,kl)
         zrtmp(3,kl)=zremsz(kl)*exp(-dpthc(1,kl))
         zrtmp(2,kl)=zremsz(kl)
         zrtmp(1,kl)=epi(kl)
         zrtmp(6,kl)=zrems(4,kl)
         enddo
c

      return
      end
      subroutine bkhsgo(sg,et,d,b,na,a,epi,ncn2,t,lpri,lfast,lun11)
c
c     this routine does the work in computing cross sections by the
c     method of barfield, et. al.
c     author:  T. Kallman (from xstar1)
c
      implicit none
c
      include './PARAM'

      integer na
c
      real*8 sg(ncn),b(na),a(11,na),epi(ncn)
      integer lpri,lfast,lun11,ncn2,nbinc
      real*8 d,t
      integer lprisv,jj,nb1,j,i,nphint,lk,kk,nskp,lrcalc
      real*8 tmp,xx,yy,sgtmp,et,epii
c
      if (lpri.gt.1) write (lun11,*)'in bkhsgo:'
     $      ,na,b,t
      lprisv=lpri
c
      jj = 1
      yy=0.
      tmp=0.
      nb1=max(1,nbinc(et,epi,ncn2))
      do j=1,nb1
         sg(j)=0.
         enddo
      i=nb1
      do while ((i.le.nphint).and.(jj.le.na))
        epii = epi(i)
        if (lpri.gt.1) write (lun11,*)i,epii,et
        xx = epii*(1.e-3) - d

        if ( xx.gt.0. ) then
          if (lpri.gt.1) write (lun11,*)d,xx,jj
          if ( xx.ge.b(jj) ) jj = jj + 1
          xx = max(xx,0.)
          yy = log10(xx)
          tmp = 0.
          do  lk = 1,11
             kk = 12 - lk
             tmp = a(kk,jj) + yy*tmp
             if (lpri.gt.1)
     $                 write (lun11,*)lk,kk,yy,tmp,a(kk,jj)
             enddo
          tmp = min(max(-50.,tmp),24.)
          sgtmp = 10.**(tmp-24.)
          sg(i)=sgtmp
          if (lpri.gt.1)
     $             write (lun11,*)i,epii,xx,
     $               tmp,sgtmp
         endif
       call enxt(et,nb1,lpri,epi,ncn2,t,lfast,lun11,
     $                  i,nskp,nphint,lrcalc)
       i=i+nskp
       enddo
      if (i.lt.nphint) then
        do j=i,nphint
          sg(j)=0.
          enddo
        endif
c
      if (lpri.gt.1) write (lun11,*)'leaving bkhsgo'
      lpri=lprisv
c
      return
      end
      subroutine bremem(lpri,lun11,xee,xpx,t,epi,ncn2,brcems,opakc)
c
c     this routine computes emissivities due to thermal bremsstrahlung.
c     author:  T. Kallman (from xstar1)
c
      implicit none
c
      include './PARAM'
c
      real*8 epi(ncn),brcems(ncn),xpx,t,xee,opakc(ncn)
      integer lpri,lun11,ncn2,numcon,kl,kk
      real*8 cc,xnx,enz2,zz,temp,gam,gau,brtmp,fbg,ekt,t6
      integer lskp,lprisv
      real*8 bbee
c
c      data cc/8.223e-15/
      data cc/1.032e-13/
c
      lskp=1
c
      ekt = t*(0.861707)
      t6 = t/100.
c
      lprisv=lpri
      if (lpri.gt.0) write (lun11,*)'in bremem',t
c
      numcon=ncn2
      do kl = 1,numcon,lskp
         brcems(kl) = 0.
         enddo
c
      xnx=xpx*xee
      enz2=(1.4)*xnx
      zz=1.
      do kk = 1,numcon,lskp
         temp = epi(kk)/ekt
         gam = zz*zz*(0.158)/t6
         gau = 1.
         if ( temp.lt.100. ) gau = fbg(temp,gam)
         brtmp = cc*xnx*enz2*gau*exp(-temp)/sqrt(t)
         brcems(kk) = brcems(kk) + brtmp
         bbee=0.
c         if ((brtmp.gt.1.e-34).and.(temp.lt.50.)) then
c           bbee=2.*(epi(kk)/3.99836e-8)**3/(exp(temp)-1.+1.e-24)
c           opakc(kk)=opakc(kk)+brtmp/(1.e-24+bbee)
c           endif
         if ( lpri.gt.0 ) write (lun11,99001) kk,
     &                zz,enz2,gam,temp,gau,brtmp,bbee,opakc(kk)
         enddo
c
      lpri=lprisv
c
      return
99001 format (' ',i6,8e12.4)
      end
      subroutine calc_maxwell_rates(lun11,lpri,coll_type,min_T,max_T,
     $  tarr, om, dE,  T,
     $ z,  degl,  degu,
     $ exc_rate, dex_rate)
c
c /*! \file calc_maxwell_rates.c
c *  \brief Calculates Maxwellian-averaged collision strengths
c */
c  /*! \brief Calculates Maxwellian-averaged collision strengths.
c    Calculates Maxwellian-averaged collision strengths for any
c    transition, given input values taken from the APED (collision type,
c    temperature minimum/maximum, and data vectors Tarr and om.  Also
c    requires the energy separating the two levels, the electron (or proton)
c    temperature, the proton number, and the upper and lower level
c    degeneracies.  Returns both the excitation and the deexcitation rate.
c    \param coll_type The collision type, as defined in calc_maxwell_rates.h
c    \param min_T The minimum temperature (in Kelvin) for which the
c    transition is defined.
c    \param max_T The maximum temperature (in Kelvin) for which the
c    transition is defined.
c    \param Tarr
c    \param om
c    \param dE The delta Energy of the transition, in keV
c    \param T  Temperature, in Kelvin, of either the electrons or protons
c    involved in the transition.
c    \param Z  The # of protons in the target ion.
c    \param degl The degeneracy of the lower level in the transition
c    \param degu The degeneracy of the upper level in the transition
c    \param exc_rate The calculated excitation rate coefficient, in cm^3/s.
c    \param dex_rate The calculated deexcitation rate coefficient, in cm^3/s.
c  */
c_______________________________________________________________
c
      real*8 calc_spline
      real*8 interpol_huntd
      real*8 calc_kato
      real*8 calc_sampson_s
      real*8 calc_sampson_p
      real*8 calc_sampson_h
      real*8 exint1
      integer coll_type
      real*8  min_T
      real*8 max_T
      real*8 Tarr(100)
      real*8 om(100)
      real*8 dE
      real*8 T
      integer lun11,lpri
      integer Z
      real*8 degl
      real*8 degu
      real*8 exc_rate
      real*8 dex_rate,pow
c
      integer calc_type
      integer N_interp, num_spline
      real*8 a,b,c,d,e
      real*8 chi, chi_inv, upsilon, rate_coeff,  st, logT, expint_2
      real*8 xs(5)
      real*8 xs9(9)
      data xs/0.00, 0.25, 0.50, 0.75, 1.00/
      data xs9/0.00, 0.125, 0.25, 0.375, 0.50, 0.675, 0.80, 0.925, 1.00/
      real*8 yp(5), yp9(9), spline_dat(9)
c /*! \file calc_maxwell_rates.h
c  * \brief Collisional excitation data types
c  * Data are stored in a number of different formats;
c  * this header lists all the currently-used versions
c  * which are found in calc_maxwell_rates.c
c  */
      integer MAX_UPS
      data MAX_UPS/20 /
c /*      !< Length of temp/om arrays in collision structure */
      integer MAX_CHI
      data MAX_CHI/200. /
c /*      !< Do not calculate rate if kT/dE greater than this.*/
      real*8 KBOLTZ
      data KBOLTZ/8.617385e-8  /
c /*      !< in units of keV/K */
      real*8 M_E
      data M_E/2.7182818284590452354  /
c /*      !< Euler e */
      real*8 UPSILON_COLL_COEFF
      data UPSILON_COLL_COEFF/8.629e-6  /
c /*      !< sqrt{2 pi / kB} hbar^2/m_e^{1.5} */
      integer E_UPSILON
      data E_UPSILON/1    /
c /*      !< Electron upsilon values (unitless) */
      integer E_RATE_COEFF
      data E_RATE_COEFF/2 /
c /*      !< Electron rate coefficient (cm^3/s) */
      integer P_UPSILON
      data P_UPSILON/3    /
c /*      !< Proton upsilon values (unitless) */
      integer P_RATE_COEFF
      data P_RATE_COEFF/4 /
c /*      !< Proton rate coefficient (cm^3/s) */

      integer BURGESSTULLY
      data BURGESSTULLY/1  /
c     /*      !< Burgess-Tully-type data*/
      integer CHIANTI_1
      data CHIANTI_1/11   /
c /*      !< CHIANTI pre-4.0 type 1 data (5 pt spline) */
      integer CHIANTI_2
      data CHIANTI_2/12   /
c /*      !< CHIANTI pre-4.0 type 2 data (5 pt spline) */
      integer CHIANTI_3
      data CHIANTI_3/13   /
c /*      !< CHIANTI pre-4.0 type 3 data (5 pt spline) */
      integer CHIANTI_4
      data CHIANTI_4/14   /
c /*      !< CHIANTI pre-4.0 type 4 data (5 pt spline) */
      integer CHIANTI_5
      data CHIANTI_5/15   /
c /*      !< CHIANTI pre-4.0 type 5 data (5 pt spline) */
      integer CHIANTI_6
      data CHIANTI_6/16   /
c /*      !< CHIANTI pre-4.0 type 6 data (5 pt spline) */

      integer CHIANTI4_1
      data CHIANTI4_1/21   /
c /*      !< CHIANTI 4.0 type 1 data (9 pt spline) */
      integer CHIANTI4_2
      data CHIANTI4_2/22   /
c /*      !< CHIANTI 4.0 type 2 data (9 pt spline) */
      integer CHIANTI4_3
      data CHIANTI4_3/23   /
c /*      !< CHIANTI 4.0 type 3 data (9 pt spline) */
      integer CHIANTI4_4
      data CHIANTI4_4/24   /
c /*      !< CHIANTI 4.0 type 4 data (9 pt spline) */
      integer CHIANTI4_5
      data CHIANTI4_5/25   /
c /*      !< CHIANTI 4.0 type 5 data (9 pt spline) */
      integer CHIANTI4_6
      data CHIANTI4_6/26   /
c /*      !< CHIANTI 4.0 type 6 data (9 pt spline) */

      integer SGC_1
      data SGC_1/31   /
c /*      !< Sampson, Goett and Clark (1983) S-type He-like data */
      integer SGC_2
      data SGC_2/32   /
c /*      !< Sampson, Goett and Clark (1983) P-type He-like data */
      integer SGC_3
      data SGC_3/33   /
c /*      !< Sampson, Goett and Clark (1983) S-type H-like data */
      integer KATO_NAKAZAKI_1
      data KATO_NAKAZAKI_1/41  /
c /*      !< Kato and Nakazaki (1989), ADNDT 42, 313 */
      integer KATO_NAKAZAKI_2
      data KATO_NAKAZAKI_2/42  /
c /*      !< Kato and Nakazaki (1989), ADNDT 42, 313 */
c /* These must be spaced by at least MAX_UPS */
      integer INTERP_E_UPSILON
      data INTERP_E_UPSILON/100     /
c /*      !< Include both left & right boundaries */
      integer INTERP_P_UPSILON
      data INTERP_P_UPSILON/200     /
c /*      !< Include both left & right boundaries */
      integer INTERP_E_RATE_COEFF
      data INTERP_E_RATE_COEFF/300  /
c /*      !< Include both left & right boundaries */
      integer INTERP_P_RATE_COEFF
      data INTERP_P_RATE_COEFF/400  /
c /*      !< Include both left & right boundaries */

      integer INTERP_E_UPS_OPEN
      data INTERP_E_UPS_OPEN/150    /
c /*      !< Include neither boundary */
      integer INTERP_P_UPS_OPEN
      data INTERP_P_UPS_OPEN/250    /
c /*      !< Include neither boundary */
      integer INTERP_E_RATE_OPEN
      data INTERP_E_RATE_OPEN/350   /
c /*      !< Include neither boundary */
      integer INTERP_P_RATE_OPEN
      data INTERP_P_RATE_OPEN/450   /
c /*      !< Include neither boundary */

      integer INTERP_E_UPS_INC_MIN
      data INTERP_E_UPS_INC_MIN/500   /
c /*      !< Include only minimum; max is out */
      integer INTERP_P_UPS_INC_MIN
      data INTERP_P_UPS_INC_MIN/600   /
c /*      !< Include only minimum; max is out */
      integer INTERP_E_RATE_INC_MIN
      data INTERP_E_RATE_INC_MIN/700  /
c /*      !< Include only minimum; max is out */
      integer INTERP_P_RATE_INC_MIN
      data INTERP_P_RATE_INC_MIN/800  /
c /*      !< Include only minimum; max is out */

      integer INTERP_E_UPS_INC_MAX
      data INTERP_E_UPS_INC_MAX/550   /
c /*      !< Include only maximum; min is out */
      integer INTERP_P_UPS_INC_MAX
      data INTERP_P_UPS_INC_MAX/650   /
c /*      !< Include only maximum; min is out */
      integer INTERP_E_RATE_INC_MAX
      data INTERP_E_RATE_INC_MAX/750  /
c /*      !< Include only maximum; min is out */
      integer INTERP_P_RATE_INC_MAX
      data INTERP_P_RATE_INC_MAX/850  /
c /*      !< Include only maximum; min is out */
      integer PROTON_BT
      data PROTON_BT/1001 /
c /*      !< For Burgess-Tully Proton excitation rates */

      calc_type = -1
      chi  = dE / (KBOLTZ*T)
      chi_inv = (KBOLTZ*T) / dE

      logT = log10(T)
      st = 0
      upsilon = 0.0
      rate_coeff = 0.0
      if (lpri.ge.2)
     $     write (lun11,*)'in calc_maxwell_rates,',t,dE,z,degl,degu,
     $                   kboltz,chi,calc_type,coll_type

      exc_rate = 0
      dex_rate = 0

c  Check to make sure temperature is inside valid range for this */
c  transition otherwise return defaults set above.  */
      if (lpri.ge.2) write (lun11,*)T,min_T,max_T
      if ((T .lt. min_T).or.(T .gt. max_T)) then
c         call errmess(lun11,35,"calc_maxwell_rates, Bad transition ")
        return
        endif

c     burgesstully=1
      if (coll_type .eq. BURGESSTULLY) then
        expint_2 = exint1(chi,2)
        a = om(1)
        b = om(2)
        c = om(3)
        d = om(4)
        e = om(5)
        upsilon = a + b*chi*expint_2 + c*chi*(1-chi*expint_2)
     $     + d*(chi/2)*(1-chi*(1-chi*expint_2)) + e*expint_2
        calc_type = E_UPSILON
        if (lpri.ge.2)
     $   write (lun11,*)'coll_type=BURGESSTULLY',chi,upsilon,calc_type
        endif

c     chianti_1=1
c     chianti_2=12
c     chianti_3=13
c     chianti_4=14
c     chianti_5=15
c     chianti_6=16
      if ((coll_type.eq.CHIANTI_1   ).or.
     $     (coll_type.eq.CHIANTI_2   ).or.
     $     (coll_type.eq.CHIANTI_3   ).or.
     $     (coll_type.eq.CHIANTI_4   ).or.
     $     (coll_type.eq.CHIANTI_5   ).or.
     $     (coll_type.eq.CHIANTI_6   )) then
        c = om(1)

        if (coll_type.eq.CHIANTI_1  ) st = 1 - (log(c)/log(chi_inv+c))
        if (coll_type.eq.CHIANTI_2  ) st = chi_inv/(chi_inv+c)
        if (coll_type.eq.CHIANTI_3  ) st = chi_inv/(chi_inv+c)
        if (coll_type.eq.CHIANTI_4  ) st = 1 - (log(c)/log(chi_inv+c))
        if (coll_type.eq.CHIANTI_5  ) st = chi_inv/(chi_inv+c)
        if (coll_type.eq.CHIANTI_6  ) st = chi_inv/(chi_inv+c)

        do mm=1,5
          spline_dat(mm) = (om(mm))
          enddo
        call prep_spline(xs, spline_dat, 5, yp)
        upsilon = calc_spline(xs, spline_dat, yp, 5, st)

        if (coll_type.eq.CHIANTI_1  ) upsilon=upsilon*log(chi_inv+M_E)
c        if (coll_type .eq. CHIANTI_2)   /* Do nothing */
        if (coll_type.eq.CHIANTI_3  ) upsilon = upsilon/(chi_inv+1)
        if (coll_type.eq.CHIANTI_4  ) upsilon = upsilon*log(chi_inv+c)
        if (coll_type.eq.CHIANTI_5  ) upsilon = upsilon/chi_inv
        if (coll_type.eq.CHIANTI_6  ) upsilon = pow(10.d0,upsilon)
        if (coll_type.eq.CHIANTI4_6 ) then
          calc_type = P_UPSILON
          else
          calc_type = E_UPSILON
          endif
        if (lpri.ge.2)
     $   write (lun11,*)'coll_type=CHIANTI 1-6',chi,upsilon,calc_type
        endif

c     chianti4_1=21
c     chianti4_2=22
c     chianti4_3=23
c     chianti4_4=24
c     chianti4_5=25
c     chianti4_6=26
      if ((coll_type.eq.CHIANTI4_1  ).or.
     $     (coll_type.eq.CHIANTI4_2  ).or.
     $     (coll_type.eq.CHIANTI4_3  ).or.
     $     (coll_type.eq.CHIANTI4_4  ).or.
     $     (coll_type.eq.CHIANTI4_5  ).or.
     $     (coll_type.eq.CHIANTI4_6  )) then
        num_spline = om(1)
        c = om(2)
        if (lpri.gt.1)
     $    write (lun11,*)'in chianti4:',om(1),num_spline,om(2),c

        if (coll_type.eq.CHIANTI4_1  )st=1 - (log(c)/log(chi_inv+c))
        if (coll_type.eq.CHIANTI4_2  )st=chi_inv/(chi_inv+c)
        if (coll_type.eq.CHIANTI4_3  )st=chi_inv/(chi_inv+c)
        if (coll_type.eq.CHIANTI4_4  )st=1 - (log(c)/log(chi_inv+c))
        if (coll_type.eq.CHIANTI4_5  )st=chi_inv/(chi_inv+c)
        if (coll_type.eq.CHIANTI4_6  )st=chi_inv/(chi_inv+c)

        do mm=1,num_spline
           spline_dat(mm) = (om(2+mm))
           if (lpri.ne.0) write (lun11,*)mm,spline_dat(mm)
           enddo
        if (num_spline .eq. 5) then
            call prep_spline(xs, spline_dat, num_spline, yp)
            upsilon = calc_spline(xs, spline_dat, yp, num_spline, st)
          else
            if (num_spline .eq. 9) then
                if (lpri.ne.0) write (lun11,*)'before prep_spline:',
     $            xs9,spline_dat
                call prep_spline(xs9, spline_dat, num_spline, yp9)
                if (lpri.ne.0) write (lun11,*)'before calc_spline:',
     $               st,yp9
                upsilon = calc_spline(xs9,spline_dat,yp9,num_spline,st)
              else
                call errmess(lun11,71,"calc_maxwell_rates, Chianti data
     $with %d values, not 5 or 9 as expected")
                stop
c     $           num_spline)
                return
              endif
          endif

        if (coll_type.eq.CHIANTI4_1) upsilon=upsilon*log(chi_inv+M_E)
c        if (coll_type .eq. CHIANTI4_2)   /* Do nothing */
        if (coll_type.eq.CHIANTI4_3) upsilon = upsilon/(chi_inv+1)
        if (coll_type.eq.CHIANTI4_4) upsilon=upsilon*log(chi_inv+c)
        if (coll_type.eq.CHIANTI4_5) upsilon = upsilon/chi_inv
        if (coll_type.eq.CHIANTI4_6) upsilon = pow(10.d0,upsilon)
        if (coll_type.eq.CHIANTI4_6) then
            calc_type =P_UPSILON
          else
            calc_type=E_UPSILON
          endif
        if (lpri.ge.2)
     $   write (lun11,*)'coll_type=CHIANTI4 1-6',chi,upsilon,calc_type
        endif
c

c     SGC_1=31
      if (coll_type .eq.SGC_1       ) then
c /* S-type transitions, He-like */
        upsilon = calc_sampson_s(om, Z, T)
        calc_type =E_UPSILON
        endif

c     SGC_2=32
      if (coll_type .eq. SGC_2) then
c /* P-type transitions, He-like */
        upsilon = calc_sampson_p(om, Z, T)
        calc_type =E_UPSILON
        endif

c     SGC_3=33
      if (coll_type .eq. SGC_3) then
c/* S-type transitions, H-like */
        upsilon = calc_sampson_h(om, Z, T)
        calc_type =E_UPSILON
        endif

c     KATO_NAKAZAKI_1=41
      if (coll_type .eq.KATO_NAKAZAKI_1) then
        upsilon = calc_kato(1, om, Z, T)
        calc_type =E_UPSILON
        endif

c     KATO_NAKAZAKI_2=42
      if (coll_type .eq.KATO_NAKAZAKI_2) then
        upsilon = calc_kato(2, om, Z, T)
        calc_type =E_UPSILON
        endif

c     PROTON_BT=1001
      if (coll_type .eq.PROTON_BT   ) then
        a = om(1)
        b = om(2)
        c = om(3)
c        /* 0.8 for the np/ne ratio, */
        if ((logT .gt. om(3)).and.( logT .lt. om(4))) then
          rate_coeff=0.8*pow(10.d0, (a + b*(logT) + c*logT*logT) )
          endif
        calc_type =P_RATE_COEFF
        if (lpri.gt.1) write (lun11,*)'coll_type=PROTON_BT'
        endif

c     INTERP_E_UPSILON=100
c     MAX_UPS=20
      if (lpri.gt.1) write (lun11,*)'coll_type:',coll_type,
     $   INTERP_E_UPSILON,INTERP_E_UPSILON+MAX_UPS,MAX_UPS
      if ((coll_type .ge.INTERP_E_UPSILON) .and.
     $    (coll_type .le.INTERP_E_UPSILON + MAX_UPS)) then
        N_interp = coll_type - INTERP_E_UPSILON
        upsilon = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = E_UPSILON
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_UPSILON'
        endif

c     INTERP_E_UPS_OPEN=150
      if ((coll_type .ge. INTERP_E_UPS_OPEN) .and.
     $     (coll_type .le. INTERP_E_UPS_OPEN + MAX_UPS)) then
        N_interp = coll_type - INTERP_E_UPS_OPEN
        if ((T .gt. min_T).and.(T .lt. max_T))
     $      upsilon = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = E_UPSILON
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_UPS_OPEN'
        endif

c     INTERP_E_UPS_INC_MIN=500
      if ((coll_type .ge. INTERP_E_UPS_INC_MIN) .and.
     $     (coll_type .le. INTERP_E_UPS_INC_MIN + MAX_UPS)) then
        N_interp = coll_type - INTERP_E_UPS_INC_MIN
        if (T .lt. max_T)
     $     upsilon = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = E_UPSILON
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_UPS_INC_MIN'
        endif

c     INTERP_E_UPS_INC_MAX=550
      if ((coll_type .ge. INTERP_E_UPS_INC_MAX) .and.
     $     (coll_type .le. INTERP_E_UPS_INC_MAX + MAX_UPS)) then
        N_interp = coll_type - INTERP_E_UPS_INC_MAX
        if (T .gt. min_T)
     $     upsilon = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = E_UPSILON
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_UPS_INC_MAX'
        endif

c     INTERP_P_UPSILON=200
      if ((coll_type .ge. INTERP_P_UPSILON) .and.
     $     (coll_type .le. INTERP_P_UPSILON + MAX_UPS)) then
        N_interp = coll_type - INTERP_P_UPSILON
        upsilon = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = P_UPSILON
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_UPSILON'
        endif

c     INTERP_P_UPS_OPEN=250
      if ((coll_type .ge. INTERP_P_UPS_OPEN) .and.
     $     (coll_type .le. INTERP_P_UPS_OPEN + MAX_UPS)) then
        N_interp = coll_type - INTERP_P_UPS_OPEN
        if ((T .gt. min_T).and.(T .lt. max_T))
     $     upsilon = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = P_UPSILON
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_UPS_OPEN'
        endif

c     INTERP_P_UPS_INC_MIN=600
      if ((coll_type .ge. INTERP_P_UPS_INC_MIN) .and.
     $     (coll_type .le. INTERP_P_UPS_INC_MIN + MAX_UPS)) then
        N_interp = coll_type - INTERP_P_UPS_INC_MIN
        if (T .lt. max_T)
     $     upsilon = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = P_UPSILON
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_UPS_INC_MIN'
        endif

c     INTERP_P_UPS_INC_MIN=650
      if ((coll_type .ge. INTERP_P_UPS_INC_MAX) .and.
     $     (coll_type .le. INTERP_P_UPS_INC_MAX + MAX_UPS)) then
        N_interp = coll_type - INTERP_P_UPS_INC_MAX
        if (T .gt. min_T)
     $     upsilon = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = P_UPSILON
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_UPS_INC_MAX'
        endif

c     INTERP_E_RATE_COEFF=300
      if ((coll_type .ge. INTERP_E_RATE_COEFF) .and.
     $     (coll_type .le. INTERP_E_RATE_COEFF + MAX_UPS)) then
        N_interp = coll_type - INTERP_E_RATE_COEFF
        rate_coeff = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = E_RATE_COEFF
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_RATE_COEFF'
        endif

c     INTERP_E_RATE_OPEN=350
      if ((coll_type .ge. INTERP_E_RATE_OPEN) .and.
     $     (coll_type .le. INTERP_E_RATE_OPEN + MAX_UPS)) then
        N_interp = coll_type - INTERP_E_RATE_OPEN
        if ((T .gt. min_T).and.(T .lt. max_T))
     $     rate_coeff = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = E_RATE_COEFF
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_RATE_OPEN'
        endif

c     INTERP_E_RATE_INC_MIN=700
      if ((coll_type .ge. INTERP_E_RATE_INC_MIN) .and.
     $     (coll_type .le. INTERP_E_RATE_INC_MIN + MAX_UPS)) then
        N_interp = coll_type - INTERP_E_RATE_INC_MIN
        if (T .lt. max_T)
     $     rate_coeff = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = E_RATE_COEFF
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_RATE_INC_MIN'
        endif

c     INTERP_E_RATE_INC_MAX=750
      if ((coll_type .ge. INTERP_E_RATE_INC_MAX) .and.
     $     (coll_type .le. INTERP_E_RATE_INC_MAX + MAX_UPS)) then
        N_interp = coll_type - INTERP_E_RATE_INC_MAX
        if (T .gt. min_T)
     $     rate_coeff = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = E_RATE_COEFF
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_E_RATE_INC_MAX'
        endif

c     INTERP_P_RATE_COEFF=400
      if ((coll_type .ge. INTERP_P_RATE_COEFF) .and.
     $     (coll_type .le. INTERP_P_RATE_COEFF + MAX_UPS)) then
        N_interp = coll_type - INTERP_P_RATE_COEFF
        rate_coeff = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = P_RATE_COEFF
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_RATE_COEFF'
        endif

c     INTERP_P_RATE_OPEN=450
      if ((coll_type .ge. INTERP_P_RATE_OPEN) .and.
     $     (coll_type .le. INTERP_P_RATE_OPEN + MAX_UPS)) then
        N_interp = coll_type - INTERP_P_RATE_OPEN
        if ((T .gt. min_T).and.(T .lt. max_T))
     $     rate_coeff = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = P_RATE_COEFF
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_RATE_OPEN'
        endif

c     INTERP_P_RATE_INC_MIN<800
      if ((coll_type .ge. INTERP_P_RATE_INC_MIN) .and.
     $     (coll_type .le. INTERP_P_RATE_INC_MIN + MAX_UPS)) then
        N_interp = coll_type - INTERP_P_RATE_INC_MIN
        if (T .lt. max_T)
     $     rate_coeff = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = P_RATE_COEFF
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_RATE_INC_MIN'
        endif

c     INTERP_P_RATE_INC_MAX=850
      if ((coll_type .ge. INTERP_P_RATE_INC_MAX) .and.
     $     (coll_type .le. INTERP_P_RATE_INC_MAX + MAX_UPS)) then
        N_interp = coll_type - INTERP_P_RATE_INC_MAX
        if (T .gt. min_T)
     $     rate_coeff = interpol_huntd(N_interp, Tarr, om, T)
        calc_type = P_RATE_COEFF
        if (lpri.gt.1) write (lun11,*)'coll_type=INTERP_P_RATE_INC_MAX'
        endif

      if (calc_type .eq. -1) then
        call errmess(lun11,43,"calc_trans_rates, Undefined collision ty
     $pe.")
        return
        endif

      if (calc_type .eq. E_UPSILON) then

        if (upsilon .le. 0.)  then
          upsilon = 0.
c /* Negative upsilon is unphysical. */
          endif

        if (chi .lt. MAX_CHI) then
          exc_rate=UPSILON_COLL_COEFF*upsilon*exp(-chi)/(sqrt(T)*degl)
          dex_rate=UPSILON_COLL_COEFF * upsilon / (sqrt(T)*degu)
          endif
        endif

      if (calc_type .eq. P_UPSILON) then

        if (upsilon .le. 0.)  then
          upsilon = 0.
c /* Negative upsilon is unphysical. */
          endif

        exc_rate = 0.
        dex_rate = 0.
        call errmess(lun11,59,"calc_rate,  Can't calculate collision str
     $ength for protons.")
        return
        endif

      if ((calc_type .eq. P_RATE_COEFF)
     $  .or.(calc_type .eq. E_RATE_COEFF)) then
        if (rate_coeff .lt. 0) then
          rate_coeff= 0.
c /* Negative rate coefficient is unphysical. */
          endif
        exc_rate = rate_coeff
        dex_rate = rate_coeff*exp(chi)*(degl/degu)
        endif
c
      return
      end
c_______________________________________________________________
      real*8 function calc_sampson_s(om, Z, Te)

      real*8 om(14)
      integer Z
      real*8 Te

c  /* These routines come from Sampson, Goett, & Clark, ADNDT 29, 467 */

      real*8 result
      real*8 kT, y
      real*8 dE, z2s_h, a1_gam, a1e_gam, a2, c0, c1, c2
      real*8 a2e, cere, cere1, s, se
      real*8 a1, a1y, E1y, Ery, Er1y
      real*8 Zeff, Zeff_e, term
      integer re
      real*8 Z2gamma, Z2gamma_e,exint_n,pow
      real*8 KBOLTZ,fre,frem
      data KBOLTZ/8.617385e-8/
c      /* kB = (keV/K) */

      Z2gamma=0.0
      Z2gamma_e=0.0
      dE     = om(1)
      a1_gam = om(2)
      a1e_gam= om(3)
      z2s_h  = om(4)
      a2     = om(5)
      c0     = om(6)
      c1     = om(7)
      c2     = om(8)
      a2e    = om(9)
      cere   = om(10)
      cere1  = om(11)
      re     = om(12)
      s      = om(13)
      se     = om(14)

      kT = KBOLTZ*Te
      y = dE/kT

      a1 = a2+1.0
      a1y = a1*y
      E1y = exint_n(y,-1.d0,1)
      Ery = exint_n(a1y,-1.d0,1)
      Er1y= exint_n(a1y,Ery,2)

      term = (c1*Ery + c2*Er1y/a1)
      if ((a1_gam .ne.0.0).and.(term.gt.0)) then
          Z2gamma = c0 + 1.333*z2s_h*E1y*exp(y) + y*exp(a1y)*term
       else
c          { /* Avoid doing exponential */
          Z2gamma = c0 + 1.333*z2s_h*E1y*exp(y)
       endif

      a1 = a2e+1
      a1y = a1*y
      Ery = exint_n(a1y,-1.d0,re)
      Er1y = exint_n(a1y,-1.d0,re+1)

      fre=re
      frem=re-1
      term = (cere*Ery/(pow(a1,frem)) + cere1*Er1y/(pow(a1,fre)))
      if ((a1e_gam .ne. 0.0).and. (term .gt. 0)) then
          Z2gamma_e = y*exp(a1y)*term
         else
          Z2gamma_e = 0.0
         endif

      Zeff = Z - s
      Zeff_e = Z - se

      result = a1_gam*Z2gamma/(Zeff*Zeff)
     $  + a1e_gam*Z2gamma_e/(Zeff_e*Zeff_e)

      calc_sampson_s=result
      return
      end
c_______________________________________________________________
      real*8 function calc_sampson_p(om, Z,Te)

      real*8 om(8)
      real*8 result
      real*8 kT, y
      real*8 dE, a, z2s, c0, cr, cr1, r, s
      real*8 Z2gamma, Zeff, a1, a1y, E1y, Ery, Er1y
      real*8 term,te,exint_n,pow
      integer z

      dE  = om(1)
      a   = om(2)
      z2s = om(3)
      c0  = om(4)
      cr  = om(5)
      cr1 = om(6)
      r   =  om(7)
      s   = om(8)

      kT = KBOLTZ*Te
      y  = dE/kT

      a1 = a+1
      a1y = a1*y

      E1y = exint_n(y,-1.d0,1)
      Ery = exint_n(a1y,-1.d0,int(r))
      Er1y= exint_n(a1y,-1.d0,int(r+1))

      term = (cr*Ery/pow(a1,r-1) + cr1*Er1y/pow(a1,r))
      if (term.gt. 0) then
          Z2gamma = c0 + 1.333*z2s*E1y*exp(y) + y*exp(a1y)*term
        else
          Z2gamma = c0 + 1.333*z2s*E1y*exp(y)
        endif

      Zeff = Z - s

      result = Z2gamma / (Zeff*Zeff)

      calc_sampson_p=result
c
      return
      end
c_______________________________________________________________
      real*8 function calc_sampson_h(om,Z,Te)

      real*8 om(7)
      integer Z
      real*8 Te

      real*8 result
      real*8 kT, y,exint_n,pow
      real*8 dE, z2s, a, c0,c1,c2,c_sw
      real*8 a1, a1y, E1y, Ery, Er1y, Zeff2
      real*8 term

      dE   = om(1)
      z2s  = om(2)
      a    = om(3)
      c0   = om(4)
      c1   = om(5)
      c2   = om(6)
      c_sw = om(7)
      Zeff2 =  Z*Z

      kT  = KBOLTZ*Te
      y   = dE/kT
      a1  = a+1
      a1y = a1*y
      E1y = exint_n(y,-1.d0,1)
      Ery = exint_n(a1y,-1.d0,1)
      Er1y = exint_n(a1y,Ery,2)
c      /* This is E_2(a1y) */

      term = (c1*Ery + c2*Er1y/a1)
      if (term.gt.0) then
          result = c0 + 1.333*z2s*E1y*exp(y) + y*exp(a1y)*term
        else
          result = c0 + 1.333*z2s*E1y*exp(y)
        endif

      result = result*2*c_sw/Zeff2

      calc_sampson_h=result
      return
      end
c_______________________________________________________________
      real*8 function calc_kato(coll_type,par,z,te)

c /* This fit comes from Kato & Nakazaki, 1989, Atomic Data and Nuclear
c           Data Tables 42, 2, 313 */

      real*8 result
      real*8 dE, kT
      real*8 A, B, C, D, E, P, Q, X1
      real*8 y,exint_n,pow
      real*8 E1y, E1Xy
      real*8 ups_nr, ups_r
      real*8 term1, term2, term3
      real*8 par(10),te
      integer z,coll_type

      dE = par(1)
c       /* in keV */
      A  = par(2)
      B  = par(3)
      C  = par(4)
      D  = par(5)
      E  = par(6)
      P  = par(7)
      Q  = par(8)
      X1 = par(9)
c
      kT = KBOLTZ*Te
      y = dE/kT

      go to (1,2,3)coll_type
 1      continue
c        case (1):  /* Simple case (eq 6, in above reference) */
          E1y = exint_n(y,-1.d0,1)
          term1 = A/y + C + (D/2)*(1-y)
          term2 = B - C*y + (D/2)*y*y + E/y
          result = y*(term1 + exp(y)*E1y*term2)
          go to 9000

 2        continue
c        case (2):  /* More complex case (eq 10-12 in above reference) */
          E1Xy = exint_n(X1*y,-1.d0,1)
          term3 = exp(y*(1-X1))

          term1 = A/y + C/X1 + (D/2)*(1/(X1*X1) - y/X1) + (E/y)*log(X1)
          term2 = exp(y*X1)*E1Xy*(B - C*y + D*y*y/2 + E/y)

          ups_nr = y*term3*( term1 + term2 )
          ups_r = P*((1 + 1/y) - term3*(X1 + 1/y)) + Q*(1-term3)

          result = ups_nr + ups_r
          go to 9000
 3        continue
 9000     continue

      calc_kato=result
c
      return
      end
c  /* Subroutine exint1(x,jump),
c   translated from Fortran.
c   This subroutine can be called by an external main program, such
c   as call_exint1.c.
c
c   jump = 1: exint1 = E1(x);
c   jump = 2: exint1 = exp(x) * E1(x);
c   jump = 3: exint1 = x * exp(x) * E1(x);
c   Returns a real*8 precision floating pointeger value.
c*/

      real*8 function exint1(x,jump)
      real*8 x
      integer jump
      real*8 EI_1
      real*8 EI_2
      real*8 EI_3
      real*8 x2
      real*8 x3
      real*8 x4,retval,pow
      real*8 a(10)
      data
     $a/7.122452e-07,1.766345e-06,2.928433e-05,.0002335379,.001664156,
     $.01041576, .05555682, .2500001, .9999999, .57721566490153/
      real*8 b(8)
      data
     $b/8.5733287401,18.059016973,8.6347608925,.2677737343,9.5733223454,
     $     25.6329561486, 21.0996530827, 3.9584969228/
c
      if (x.eq.0.d0) then
         exint1=0.
         return
      else
        if (x .le. 1.0d0) then
          EI_1 = ((((((((a(1) * x - a(2)) * x + a(3)) * x
     $       - a(4)) * x + a(5)) * x - a(6)) * x
     $       + a(7)) * x - a(8)) * x + a(9)) * x
     $       - log(x) - a(10)
          EI_2 = exp(x) * EI_1
          EI_3 = x * EI_2
        else
          x2 = pow(x,2.d0)
          x3 = pow(x,3.d0)
          x4 = pow(x,4.d0)
          EI_3 = (x4 + b(1) * x3 + b(2) * x2
     $     + b(3) * x + b(4)) /
     $      (x4 + b(5) * x3 + b(6) *x2
     $       + b(7) * x + b(8))
          EI_1 = EI_3 / (x * exp(x))
          EI_2 = EI_3 / x
        endif
      endif
c  /* In K&R C, the switch argument has to be of type
c     int.Can terminate with return or break.
c     */
      go to (1,2,3)jump
 1    continue
c     case 1:
      retval=EI_1
      go to 9000
 2    continue
c     case 2:
      retval=EI_2
      go to 9000
 3    continue
c     case 3:
      retval=EI_3
 9000 continue
c      default:
      exint1=retval
      return
      end
c_______________________________________________________________
      real*8 function exint_n(x,E1x_in,n)
      real*8 x
      real*8 E1x_in
      integer n

c/* Calculates E_n(x), for n = 1,2,3,4. If n.gt.1, and E1x_in .ge. 0, then
c   routine uses E1x for E_1(x) in calculating recurrance.
c   Returns a real*8 precision floating pointeger value. */
c
      real*8 E1x
      real*8 x2, x3, x4

      real*8 a(10)
      data
     $a/7.122452e-07,1.766345e-06,2.928433e-05,.0002335379,.001664156,
     $.01041576, .05555682, .2500001, .9999999, .57721566490153/
      real*8 b(8)
      data
     $b/8.5733287401,18.059016973,8.6347608925,.2677737343,9.5733223454,
     $25.6329561486, 21.0996530827, 3.9584969228/
c
      E1x = E1x_in
      if ((n.eq.1).or.(E1x .lt. 0)) then
        if (x.eq.0.) then
           exint_n=0.
           return
        else
          if (x .le. 1.0) then
            E1x = ((((((((a(1)*x-a(2))*x+a(3))*x-a(4))*x
     $         +a(5))*x-a(6))*x+a(7))*x
     $         -a(8))*x+a(9))*x - log(x) - a(10)
          else
             x2 = x*x
             x3 = x2*x
             x4 = x2*x2
             E1x = (x4+b(1)*x3+b(2)*x2+b(3)*x+b(4))/
     $         (x4+b(5)*x3+b(6)*x2+b(7)*x+b(8))
             E1x = E1x / (x * exp(x))
          endif
        endif
      endif
c
      go to (1,2,3,4)n
 1      continue
c       case (1):
        retval= (E1x)
        go to 9000
 2      continue
c       case (2):
        retval=(exp(-x) - x*E1x)
        go to 9000
 3      continue
c       case (3):
        retval= (x*x*E1x + exp(-x)*(1-x))/2.
        go to 9000
 4      continue
c       case (4):
        retval= (exp(-x)*(2-x+x*x) - x*x*x*E1x)/6.
        go to 9000
        call errmess(6,27,"exint_n. Bad value for n = ")
 9000   continue
        exint_n=retval
        return
        end
      subroutine errmess(lun11,nlen,str1)
      character*(*) str1
      integer lun11
      write (lun11,*) 'in errmess:',str1
      return
      end
c_______________________________________________________________
      real*8 function interpol_huntd(n,x,y,z)
      integer n
      real*8 x(n)
      real*8 y(n)
      real*8 z
      real*8 grad,d1,df,f,pow
      integer jl, jm, ju
      logical inc
c
      inc = .false.
      if (x(n-1) .gt. x(1)) inc=.true.
      if (( inc .and.((z .lt. x(1)) .or. (z .gt. x(n-1)))) .or.
     $ (.not.(inc).and.(z .gt. x(1) .or. z .lt. x(n-1)))) then
c        write (6,*)"interpol_huntd: Asking for, min is, max is",z,
c     $        x(1),x(n-1)
c        write (6,*)"interpol_huntd: Cannot extrapolate"
        return
        endif
      jl = 0
      ju = n
      do while (ju - jl .gt. 1)
        jm = (ju + jl) / 2
        if ((z .gt. x(jm)) .eqv. inc) then
          jl = jm
        else
          ju = jm
        endif
      enddo
c     /* ------	Z is sandwiched between JL and JU ------ */
      if ((x(jl) .gt. 0.).and.(x(ju).gt.0.).and.
     $    (y(jl) .gt. 0.).and.(y(ju) .gt. 0.)) then
        grad = (log10(y(ju)) - log10(y(jl))) /
     $      (log10(x(ju)) - log10(x(jl)))
        df = grad * (log10(z) - log10(x(jl)))
        d1 = log10(y(jl)) + df
        f = pow(10.d0, d1)
      else
        f = y(jl)+(y(ju)-y(jl))*((z-x(jl))/(x(ju)-x(jl)))
      endif
      interpol_huntd=f
      return
      end
c_______________________________________________________________
c /* Rewritten from Numerical recipes.  Free version must be obtained */
c /* before release! */

      subroutine  prep_spline(x,y,n,y2)
c
      real*8 x(n)
      real*8 y(n)
      integer n
      real*8 y2(n)

      integer i,k
      real*8 p, qn, sig, un
      real*8 u(500)

      y2(1)=0.0
      u(1)=0.0

c      write (6,*)'in prep_spline',n
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0
        y2(i)=(sig-1.0)/p
        u(i)=(y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1))
        u(i)=(6.0*u(i)/(x(i+1)-x(i-1))-sig*u(i-1))/p
c        write (6,*)i,x(i),y(i),sig,p,y2(i),u(i)
      enddo

      qn=0.0
      un=0.0

      y2(n-1)=(un-qn*u(n-2))/(qn*y2(n-2)+1.0)
c      write (6,*)un,qn,u(n-2),y2(n-2),y2(n-1)

      do k=n-2,0,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
         enddo

c      free(u)
      return
      end
c_______________________________________________________________

      real*8 function calc_spline(xa, ya,y2a,n,x)
c
      integer n
      real*8 xa(n)
      real*8 ya(n)
      real*8 y2a(n)
      real*8 x

      integer klo,khi,k
      real*8 h,b,a
      real*8 result

      klo=0
      khi=n-1

c      write (6,*)'in calc_spline',x,n
      do while (khi-klo .gt. 1)
        k=(khi+klo)/2
c        write (6,*)k,khi,klo,xa(k),x
        if (xa(k) .gt. x) then
            khi=k
          else
            klo=k
          endif
        enddo
      h=xa(khi)-xa(klo)
      if (h.eq.0.)then
        write (6,*)"calc_spline: Bad x array input"
        stop
        endif

      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      result=a*ya(klo)+b*ya(khi)
     $ +((a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi))*(h*h)/6.0

      calc_spline=result
      return
      end
      real*8 function pow(x,y)
      real*8 x,y
      pow=x**y
      return
      end
      subroutine calt57(te,den2,e,ep,n,cion,crec,lun11,lpri)
c
c  this routine calculates rates for type 57 data, collisional ionization
c
c  temp   temperature in K
c  den    electron density in cm-3
c  e      level energy in eV (first real*8 in type 6)
c  ep     ionization potential in eV (forth real*8 in type 6)
c  n      level's principal quatum number (first integer in type 6)
c  cion   ionization rate in s-1.cm+3
c  crec   3-body recombination rate in s-1.cm+6. THIS VALUE MUST BE
c         MUTIPLIED BY (stat. weigth level/stat. weigth recombining
c          level)
c      author: M. Bautista
c
       implicit none
c
       real*8 te,den,den2,e,ep,cion,crec
       integer n,lun11,lpri
       real*8 rk,cb,rio,rc,temp,tmin,rno,rno2,ciono,
     $      beta,wte,wtm,e2,e3,ete,etm
       real*8 rn
c
       data e3/0./,e2/0./

       cion=0.
       crec=0.
       rn=float(n)
       rk=1.16058e+4
       cb=13.605692*1.6021e-19/1.3805e-23
       if (lpri.gt.1)
     $  write (lun11,*)'entering calt57:',den2,te,
     $   ep,e,n
       if (ep.lt.e) return
       rio=(ep-e)/13.6
       rc=sqrt(rio)*rn
c       if (den.gt.1.e+18) then
c        print*,'density too high for cics'
c        return
c       endif
       den=min(den2,1.e+18)
       tmin=3.8e+4*rc*sqrt(rc)
       temp=te
       if (te.lt.tmin) then
        temp=tmin
       endif
       rno=SQRT(1.8887E+8*rc/den**0.3333)
       rno2=(1.814e26*(rc**6)/2./den)**0.13333
       rno=min(rno,rno2)
       if (int(rno).gt.n) then
        call irc(n,temp,rc,rno,ciono,lpri,lun11)
       if (lpri.gt.1)
     $  write (lun11,*)'in calt57:',rno,den,rc,tmin,te,
     $   ciono
c
c extrapolates to actual temperature below Tmin
c
        cion=0.
        crec=0.
        if (te.lt.tmin) then
         beta=.25*(sqrt((100.*rc+91.)/(4.*rc+3.))-5.)
         wte=(log(1.+te/cb/rio))**(beta/(1.+te/cb*rio))
         wtm=(log(1.+tmin/cb/rio))**(beta/(1.+tmin/cb*rio))
         call eint(rio/te*cb,ete,e2,e3)
         if (ete.lt.1.e-20) return
         call eint(rio/tmin*cb,etm,e2,e3)
         if (lpri.gt.1)
     $    write (lun11,*)'in calt57:',
     $     n,temp,tmin,rc,rno,cion,te,ete,etm,wte,wtm,cion
         cion=ciono*sqrt(tmin/te)*ete/(etm+1.e-30)*wte/(wtm+1.e-30)
        else
         cion=ciono
        endif
c
        if (cion.le.1.e-24) return
c
c        se=log(cion)
c        sr=se-36.1136-1.5*log(te)+(13.6*
c     c                 rc*rc*(1./float(n*n)-1./rno/rno)*rk/te)
        if (lpri.gt.1)
     $   write (lun11,*)te,n,rno,rk,te,cion
        cion=cion/float(n*n)
c        crec=expo(sr)/float(n*n)
        crec=cion*(2.0779e-16)*
     $    exp((ep-e)*rk/te)
c     $    exp(13.6*rc*rc*(1./float(n*n)-1./rno/rno)*rk/te)
     $              /te**(1.5)
       endif
        return
        end
       subroutine calt60_62(temp,m,idata,np1r,np1i,rdat1,idat1,Upsilon)
c
c  This rutine takes the coefficients in data type 60 and 62 (reals
c  as dtype and integers as itype) and return effective collision
c  strengths according to fits by Callaway (1994).
c  "temp" is the temperature in Kelvin and "m" is the number of
c  reals in dtype. "idata" is the data type either 60 or 62.
c      author: M. Bautista
c
        implicit none
        include './PARAM'
c
        integer m
        real*8 rdat1(nrdat1)
        integer idat1(nidat1)
c        real*8 dtype(m)
        real*8 temp,Upsilon
c        integer itype(m)
        integer i,idata,np1i,np1r
        real*8 t1,de,tmax,tt,rat
c
        t1=temp*6.33652e-6
        if (temp.gt.1.e9) t1=6.33652e+3
        de=1./float(idat1(np1i))**2-1./float(idat1(np1i-1+2))**2
        tmax=4.*de
        tmax=1.
        tt=t1
        if (t1.gt.tmax) tt=tmax
        if (idata.eq.60) then
         rat=0.
         do i=1,m-2
          rat=rat+rdat1(np1r-1+i+2)*(tt**(i-1))
         enddo
         Upsilon=rat
        else
         rat=0.
         do i=1,m-5
          rat=rat+rdat1(np1r-1+i+2)*(tt**(i-1))
         enddo
         Upsilon=rat+rdat1(np1r-1+m-2)*log(rdat1(np1r-1+m-1)*tt)
     $               *exp(-rdat1(np1r-1+m)*tt)
        endif
c
         if (t1.gt.tt) then
          upsilon=Upsilon*(1.+log(t1/tmax)/(log(t1/tmax)+1.))
         endif
c
      return
      end
      subroutine calt66(temp,np1r,rdat1,gamma)
c
c   Takes coefficients in data type 66 and returns effective collision
c    strenghts for He-like ions according to Kato & Nakazaki (1989)
c    eq. (6).
c      author: M. Bautista
c
       implicit none
       include './PARAM'
c
       real*8 rdat1(nrdat1)
       real*8 temp,gamma
       real*8 eboltz,y,a,b,c,d,e,em1,gam1,dele,gam2,gam3
       integer np1r
c
       eboltz=1.160443e+04
       dele=rdat1(np1r)
       y=dele/temp*eboltz
c
       if (y.lt.1.e-20) then
        print*,'error in calt66. y too low. y=',y
        return
       endif
       if (y.gt.1.e+20) then
        gamma=0.
        return
       endif
c
       if (y.gt.77.) y=77.
       call expint(y,em1)
       a=rdat1(np1r-1+2)
       b=rdat1(np1r-1+3)
       c=rdat1(np1r-1+4)
       d=rdat1(np1r-1+5)
       e=rdat1(np1r-1+6)
       gam1=y*((a/y+c)+d*.5*(1.-y))+em1*(b-c*y+d*y*y*.5+e/y)
       dele=rdat1(np1r-1+7)
       y=dele/temp*eboltz
       if (y.gt.77.) y=77.
       call expint(y,em1)
       a=rdat1(np1r-1+8)
       b=rdat1(np1r-1+9)
       c=rdat1(np1r-1+10)
       d=rdat1(np1r-1+11)
       e=rdat1(np1r-1+12)
       gam2=y*((a/y+c)+d*.5*(1.-y))+em1*(b-c*y+d*y*y*.5+e/y)
       dele=rdat1(np1r-1+13)
       y=dele/temp*eboltz
       if (y.gt.77.) y=77.
       call expint(y,em1)
       a=rdat1(np1r-1+14)
       b=rdat1(np1r-1+15)
       c=rdat1(np1r-1+16)
       d=rdat1(np1r-1+17)
       e=rdat1(np1r-1+18)
       gam3=y*((a/y+c)+d*.5*(1.-y))+em1*(b-c*y+d*y*y*.5+e/y)
       gamma=gam1+gam2+gam3
       return
       end
       subroutine calt67(temp,np1r,rdat1,gamma)
c
c   Takes coefficients in data type 67 and returns effective collision
c    strenghts for He-like ions according to Keenan, McCann, & Kingston
c    (1987) eq. (2)
c      author: M. Bautista
c
       implicit none
c
      integer nptmpdim
      parameter (nptmpdim=200000)
      include './PARAM'
c
       real*8 rdat1(nrdat1)
       real*8 gamma,temp,tp
       integer np1r
c
       tp=log10(temp)
       gamma=rdat1(np1r-1+1)+rdat1(np1r-1+2)*tp+rdat1(np1r-1+3)*tp*tp
c
       return
       end
       subroutine calt68(temp,np1r,np1i,rdat1,idat1,gamma)
c
c    Takes coefficients in data type 68 and returns effective collision
c    strenghts for He-like ions according to Sanpson & Zhang.
c      author: M. Bautista
c
       implicit none
       include './PARAM'
c
       real*8 rdat1(nrdat1)
       integer idat1(nidat1)
       real*8 z,tt,gamma,temp
       integer np1r,np1i
c
       z=float(idat1(np1i-1+3))
       tt=log10(temp/z/z/z)
       gamma=rdat1(np1r-1+1)+rdat1(np1r-1+2)*tt+rdat1(np1r-1+3)*tt*tt
c
       return
       end
       subroutine calt69(temp,m,np1r,rdat1,gamma)
c
c   Takes coefficients in data type 69 and returns effective collision
c    strenghts for He-like ions according to Kato & Nakazaki (1989)
c    eq. (6). m is the dimension of dtype69
c      author: M. Bautista
c
       implicit none
       include './PARAM'
c
       integer m
       real*8 rdat1(nrdat1)
       real*8 gamma,temp,eboltz,dele,y
       real*8 a,b,c,d,e,p,q,x1,em1,gnr,gr
       integer np1r
c
       eboltz=1.160443e+04
       dele=rdat1(np1r)
       y=dele/temp*eboltz
c
       if (y.lt.1.e-20) then
        print*,'error in calt69. y too low. y=',y
        return
       endif
       if (y.gt.1.e+20) then
        gamma=0.
        return
       endif
c
       if (y.gt.77.) y=77.
       if (y.lt.5.e-2) y=5.e-2
       call expint(y,em1)
       a=rdat1(np1r-1+2)
       b=rdat1(np1r-1+3)
       c=rdat1(np1r-1+4)
       d=rdat1(np1r-1+5)
       e=rdat1(np1r-1+6)
       if (m.eq.6) then
        gamma=y*((a/y+c)+d*.5*(1.-y))+em1*(b-c*y+d*y*y*.5+e/y)
c        write(2,*)temp,y,em1,gamma
       else
        p=rdat1(np1r-1+7)
        q=rdat1(np1r-1+8)
        x1=rdat1(np1r-1+9)
        call expint(y*x1,em1)
        gnr=a/y+c/x1+d*.5*(1./x1/x1-y/x1)+e/y*log(x1)+
     #   em1/y/x1*(b-c*y+d*y*y*.5+e/y)
        gnr=gnr*y*exp(y*(1.-x1))
        gr=p*(1.+1./y)*(1.-exp(y*(1.-x1))*(x1+1/y)/(1.+1./y)) +
     #     q*(1.-exp(y*(1.-x1)))
        gamma=gnr+gr
       endif
c
      return
      end
      subroutine calt70(temp,den,eth,ic,m,np1r,np1i,rdat1,idat1,
     1                  nx,xe,xs,rec,al,lun11,lpri)
c
c  This routine takes the coefficients in data type 70 (dtype70 reals
c  in itype70 integers) and returns the recombination rate (in s-1cm-3)
c  and the correstpondent phot. x-section for the superlevel. m is the
c  dimension of dtype70. nx is the number of points in the x-section
c  xe() contains the photon energy in Ry and xx() is the x-section
c  in Mb.
c  temp, den, and ic are the temperature, electron density
c  and effective charge of the ion respectively.
c  eth is the threshold energy for the superlevel in Ry.
c      author: M. Bautista
c

       implicit none
      integer nptmpdim
c     parameter (nptmpdim=200000)
      parameter (nptmpdim=10000)
      include './PARAM'
c
      integer m
      real*8 rdat1(nrdat1)
      integer idat1(nidat1)
      real*8 xe(nptmpdim),xs(nptmpdim),rne,rte,rm,
     $      rec1,rec2,rec,scale,al,crit,temp,den,eth,dt
      integer nden,ntem,nxs,in,it,kt1,i1,imax,nx,
     $      lun11,lpri,lprim,ic,i,np1r,np1i
c
c     alpf: fitting coef. for hydrogenic recombination n=12,l=0
c      dimension alpf(3)
c      data alpf/-7.1094841E-02,-9.0274535E-02,-14.26129/
c
      if (lpri.gt.1)
     $  write (lun11,*)'in calt70:',temp,den,eth,ic,m,rdat1(np1r),
     $                  idat1(np1i-1+1)
      rne=log10(den)
      rte=log10(temp)
      nden=idat1(np1i-1+1)
      ntem=idat1(np1i-1+2)
      nxs=idat1(np1i-1+3)
      if (nden.gt.1) then
      if (rne.gt.rdat1(np1r-1+nden)) then
c       write (lun11,*)'DENSITY TOO HIGH AT SUPREC'
c       write (lun11,*)'z=',ic,' temp=',temp,' Ne=',den
c       return
       rne=min(rne,rdat1(np1r-1+nden))
      endif
      if (rne.le.rdat1(np1r-1+1)) then
       in=1
      else
       in=int(rne/rdat1(np1r-1+nden)*nden)-1
       if (in.ge.nden) in=in-1
 5     in=in+1
       if (in.lt.nden .and. rne.ge.rdat1(np1r-1+in+1)) goto 5
       if (rne.lt.rdat1(np1r-1+in)) then
        in=in-2
        goto 5
       endif
      endif
      else
       in=1
      endif
      if (rte.lt.rdat1(np1r-1+nden+1)) then
       it=1
      else
       dt=(rdat1(np1r-1+nden+ntem)-rdat1(np1r-1+nden+1))/float(ntem)
       it=int((rte-rdat1(np1r-1+nden+1))/dt)
 6     it=it+1
       if (it.ge.ntem) then
        it=ntem-1
       else
        if (rte.ge.rdat1(np1r-1+nden+it+1)) goto 6
        if (rte.lt.rdat1(np1r-1+nden+it)) then
         it=it-2
         goto 6
        endif
       endif
      endif
      kt1=nden+ntem+(in-1)*ntem+it
      rm=(rdat1(np1r-1+kt1+1)-rdat1(np1r-1+kt1))
     $    /(rdat1(np1r-1+nden+it+1)-
     #    rdat1(np1r-1+nden+it))
      rec1=rdat1(np1r-1+kt1)+rm*(rte-rdat1(np1r-1+nden+it))
      if (nden.gt.1) then
       kt1=kt1+ntem
       rm=(rdat1(np1r-1+kt1+1)-rdat1(np1r-1+kt1))
     $    /(rdat1(np1r-1+nden+it+1)-
     #    rdat1(np1r-1+nden+it))
       rec2=rdat1(np1r-1+kt1)+rm*(rte-rdat1(np1r-1+nden+it))
       rm=(rec2-rec1)/(rdat1(np1r-1+in+1)-rdat1(np1r-1+in))
       rec=rec1+rm*(rne-rdat1(np1r-1+in))
      else
       rec=rec1
      endif
      if (lpri.gt.1) write (lun11,*)nden,ntem,it,in,kt1,
     $             rdat1(np1r-1+kt1),rm,rec2,rec1,rec
      rec=10.**rec
c
      i1=ntem*nden+ntem+nden
      do i=1,nxs
       xe(i)=rdat1(np1r-1+i1+(i-1)*2+1)
       xs(i)=rdat1(np1r-1+i1+(i-1)*2+2)
      if (lpri.gt.1)
     $  write (lun11,*)i,xe(i),xs(i)
      enddo
      lprim=0
      call milne(temp,nxs,xe,xs,eth,al,lun11,lprim)
      scale=rec/(1.e-24+al)
      if (lpri.gt.1)
     $ write (lun11,*)'in calt70:',rec,al,scale,xs(1),nxs,eth
      crit=1.e-6
      imax=nxs
      do i=1,nxs
       xs(i)=xs(i)*scale
       xs(i)=min(xs(i),1000000.)
       if (xs(i).gt.xs(1)*crit) imax=i
      if (lpri.gt.1)
     $  write (lun11,*)i,xe(i),xs(i)
      enddo
      nxs=imax
      nx=nxs
c
      return
      end
      subroutine calt71(temp,den,ic,m,np1r,np1i,rdat1,idat1,wav,aij,
     $                  lun11,lpri)
c
c  This rutine takes the coefficients in data type 71 (dtype71 reals
c  in itype71 integers) and returns the radiative transition prbability
c  (in s-1) from the superlevels to the spectroscopic level given by
c  itype71(3).
c  The wavelength for the transition is also given in wav
c  temp, den, and ic are the temperature, electron density
c  and effective charge of the ion respectively.
c      author: M. Bautista
c
       implicit none
      include './PARAM'
c
      integer m
      real*8 rdat1(nrdat1)
      integer idat1(nidat1)
      real*8 wav,aij,temp,den,rne,rte,dtmp,rm,rec1,
     $     rec,rec2
      integer lun11,lpri,nden,ntem,in,it,kt1,ic
      integer javi,np1r,np1i
c
      javi=m
      m=javi
c
      rne=log10(den)
      rte=log10(temp)
      nden=idat1(np1i-1+1)
      ntem=idat1(np1i-1+2)
      if (lpri.ne.0) write (lun11,*)'in calt71:',nden,ntem

      if (nden.eq.1 .and. ntem.eq.1) then
        if (rdat1(np1r-1+3).gt.30.) then
          dtmp=log10(rdat1(np1r-1+3))
        else
          dtmp=rdat1(np1r-1+3)
        endif
       wav=rdat1(np1r-1+4)
       aij=10.**dtmp
c       aij=min(aij,1.e+12)
       if (lpri.ne.0) write (lun11,*)'early return',aij,wav
       return
      endif
      if (rne.gt.rdat1(np1r-1+nden)) then
        if (lpri.ne.0) then
          write (lun11,*)'DENSITY TOO HIGH AT CALT71'
          write (lun11,*)'z=',ic,' temp=',temp,' Ne=',den,nden,
     $                   rdat1(np1r-1+nden)
          endif
        rne=min(rne,rdat1(np1r-1+nden))
      endif
      if (rte.gt.(rdat1(np1r-1+nden+ntem)+1.)) then
       rte=rdat1(np1r-1+nden+ntem)+1.
      endif
      if (rte.lt.(rdat1(np1r-1+nden+1)-1.)) then
       rte=rdat1(np1r-1+nden+1)-1.
      endif
c
      wav=rdat1(np1r-1+nden*ntem+nden+ntem+1)
      if (rne.le.rdat1(np1r-1+1)) then
       in=1
      else
       in=0
 5     in=in+1
       if (rne.ge.rdat1(np1r-1+in+1).and.in.lt.nden) goto 5
      endif
      if (rte.lt.rdat1(np1r-1+nden+1)) then
       it=1
      else
       it=0
 6     it=it+1
       if (it.ge.ntem) then
        it=ntem-1
       else
        if (rte.ge.rdat1(np1r-1+nden+it+1)) goto 6
       endif
      endif
      kt1=nden+ntem+(in-1)*ntem+it
      rm=(rdat1(np1r-1+kt1+1)-rdat1(np1r-1+kt1))
     $   /(rdat1(np1r-1+nden+it+1)-
     #    rdat1(np1r-1+nden+it))
      rec1=rdat1(np1r-1+kt1)+rm*(rte-rdat1(np1r-1+nden+it))
      kt1=kt1+ntem
      rm=(rdat1(np1r-1+kt1+1)-rdat1(np1r-1+kt1))
     $    /(rdat1(np1r-1+nden+it+1)-
     #    rdat1(np1r-1+nden+it))
      rec2=rdat1(np1r-1+kt1)+rm*(rte-rdat1(np1r-1+nden+it))
c
      rm=(rec2-rec1)/(rdat1(np1r-1+in+1)-rdat1(np1r-1+in))
      rec=rec1+rm*(rne-rdat1(np1r-1+in))
      aij=10.**rec
c      aij=min(aij,1.e+12)
      if (lpri.ne.0) write (lun11,*)'late return',rm,rec2,
     $       rec1,rec,aij,wav
c
      return
      end
      subroutine calt72(temp,np1r,rdat1,nrdt,rate)
c
c   Takes coefficients in data type 72 and returns capture rates
c   (in s^-1) for DR through satellite levels considered explicitly.
c      author: M. Bautista
c
       implicit none
       include './PARAM'
c
c
       real*8 rdat1(nrdat1),temp,rate,dele,s,rtmp
       integer np1r,nrdt
c
       dele=rdat1(np1r-1+2)
       s=4.141292e-22/(temp**1.5)
       rtmp=rdat1(np1r+2)
       if (nrdt.lt.3) rtmp=1.
       rate=s*exp(-dele/temp)*rdat1(np1r-1+1)*rtmp/2.*.5
c
       return
       end
       subroutine calt73(temp,np1r,np1i,rdat1,idat1,crate)
c
c   Takes coefficients in data type 73 and returns excitation rate times
c   the statistical weight of the lower level (w_i C(i,j) in s^-1
c   cm^3).
c      author: M. Bautista
c
       implicit none
       include './PARAM'
c
       real*8 rdat1(nrdat1)
       integer idat1(nidat1)
       real*8 temp,crate,boltzk,const,z,y,gam,zeff,z2s
       real*8 a,co,cr,cr1,r,e1,em1,ee1,ee2,ee3,er,er1,qij
       integer np1r,np1i
c
       boltzk=1.578876e+05
       const=5.46538e-11
       z=float(idat1(np1i-1+3))
       y=z*z*rdat1(np1r)*boltzk/temp
c
       gam=0.
       if (rdat1(np1r-1+2).ge. 0.1) gam=-.2
       if (rdat1(np1r-1+2).gt. 0.01 .and. rdat1(np1r-1+2).lt. 0.1)
     $    gam=0.
       if (rdat1(np1r-1+2).le. 0.01) gam=0.2
       zeff=float(idat1(np1i+2))-gam
       z2s=rdat1(np1r-1+2)
       a=rdat1(np1r-1+3)
       co=rdat1(np1r-1+4)
       cr=rdat1(np1r-1+5)
       cr1=rdat1(np1r-1+6)
       r=rdat1(np1r-1+7)
       if (y.gt.40.) then
        crate=0.
        return
       endif
        call expint(y,em1)
        e1=em1/y*exp(-y)
       if (y*a+y.le.80) then
        call eint(y*a+y,ee1,ee2,ee3)
       else
        ee1=0.
        ee2=0.
        ee3=0.
       endif
       er=0.
       er1=0.
       if (r.eq.1.) then
        er=ee1
        er1=ee2
       endif
       if (r.eq.2.) then
        er=ee2
        er1=ee3
       endif
       if (y*a+y.le.40) then
        qij=co*exp(-y)+1.55*z2s*e1+y*exp(y*a)*(cr*er/(a+1.)**(r-1.)
     #      +cr1*er1/(a+1.)**r)
       else
        qij=co*exp(-y)+1.55*z2s*e1
       endif
       crate=qij*boltzk/temp*sqrt(temp)/zeff/zeff*const
       if (crate.lt.0.) crate=0.
c
       return
       end
       subroutine calt74(temp,np,xse,xss,nddd,np1r,rdat1,rate,alpha)
c
c   Takes coefficients in data type 74 and any given radiation flux
c   array given in xse(i) (energiesin eV) and xss(i) (flux values)
c   and return the the analytic integral over resonances in the
c   cross sections represented by delta functions.
c   The routine also returns the DR recombination coefficient (in
c   s-1cm-3) for the given value of temp (in Kelvins). alpha MUST
c   be mutiplied by the stadistical of the recombined state and
c   divided by that of the recombining state.
c
c   np is the number of points xse() and nd is the number of real
c   values in dtype74()
c      author: M. Bautista
c
       implicit none
       include './PARAM'
c
       real*8 rdat1(nrdat1)
       integer nddd
       real*8 temp,rate,alpha,te,ryk,ry,xt,x,hgh,
     $      rm,xsec,factor
       integer np,m,i,ipos,ip,np1r
       real*8 xse(np),xss(np)
c
       te=temp*1.38066e-16
       ryk=4.589343e+10
       ry=13.605692
       m=(nddd-1)/2
c
       xt=rdat1(np1r-1+1)
       x=rdat1(np1r-1+2)
       hgh=rdat1(np1r-1+2+m)
       alpha=0.
       if (x/ryk/te.lt.40.) then
        alpha=exp(-x/ryk/te)*(x+xt)*(x+xt)*hgh
       endif
       do i=2,m
        x=rdat1(np1r-1+1+i)*ry
        hgh=rdat1(np1r-1+1+i+m)
        x=x/ry
        if (x/ryk/te.lt.40.) then
         alpha=alpha+exp(-x/ryk/te)*(x+xt)*(x+xt)*hgh
        endif
       enddo
        factor=213.9577e-9
        alpha=alpha*factor/(te**1.5)/ryk/ryk
c
       if (xse(np).lt.(x+xt)*ry) then
        rate=0.e0
        return
       endif
c
       x=(rdat1(np1r-1+2)+rdat1(np1r-1+1))*ry
       i=np/2
  5    if (xse(i).ge.x) then
        i=i-1
        goto 5
       endif
       i=i-1
  10   i=i+1
       if(xse(i).lt.x.and.xse(i+1).ge.x)then
        ipos=i
       else
        goto 10
       endif
c
       rm=(xss(ipos+1)-xss(ipos))/(xse(ipos+1)-xse(ipos))
       xsec=xss(ipos)+rm*(x-xse(ipos))
       rate=xsec*rdat1(np1r-1+2+m)
       do  i=2,m
        x=(rdat1(np1r-1+1+i)+rdat1(np1r-1+1))*ry
        hgh=rdat1(np1r-1+1+i+m)
        ip=ipos
 20     if (xse(ip).lt.x) then
         ip=ip+1
         goto 20
        endif
        ip=ip-2
        ip=ip+1
        if (ip.le.np) then
          ipos=ip
          rm=(xss(ipos+1)-xss(ipos))/(xse(ipos+1)-xse(ipos))
          xsec=xss(ipos)+rm*(x-xse(ipos))
          rate=rate+xsec*hgh
          endif
        enddo
c
        rate=rate*4.752e-22
c
      return
      end
      subroutine calt77(lpri,lun11,temp,den,m,
     $                       np1r,np1i,rdat1,idat1,cul,clu)
c
c  This rutine takes the coefficients in data type 77 (dtype77 reals
c  and itype77 integers) and returns the collisional transition rates
c  (in s-1) from the superlevel (cul) and to the superlevel (clu) from
c  the spectroscopic level given by itype77(3).
c  The wavelength for the transition is also given in wav
c  temp, den, and ic are the temperature, electron density
c  and effective charge of the ion respectively.
c      author: M. Bautista
c
      implicit none
      include './PARAM'
c
c
      real*8 rdat1(nrdat1)
      integer idat1(nidat1)
      integer m
      real*8 temp,den,cul,clu,rne,rte,div,rm,rec1,rec2,
     $         gg,wav,rec,xt
      integer nden,ntem,in,it,kt1,nll,k,nl1,nl2,il,
     $        lpri,lun11,np1r,np1i
c
      rne=log10(den)
      rte=log10(temp)
      nden=idat1(np1i-1+1)
      ntem=idat1(np1i-1+2)
      if (rne.gt.rdat1(np1r-1+nden)) then
c       print*,'DENSITY TOO HIGH AT CALT77'
c       print*,'z=',ic,' temp=',temp,' Ne=',den,nden,rdat1(np1r-1+nden)
       rne=rdat1(np1r-1+nden)
      endif
      if (rte.gt.(rdat1(np1r-1+nden+ntem)+1.)) then
c       print*,'TEMPERATURE TOO HIGH AT CALT77'
c       print*,'z=',ic,' temp=',temp,' Ne=',den
         rte=(rdat1(np1r-1+nden+ntem)+1.)
      endif
      if (rte.lt.(rdat1(np1r-1+nden+1)-1.)) then
       rte=rdat1(np1r-1+nden+1)-1.
      endif
c
      wav=rdat1(np1r-1+nden*ntem+nden+ntem+1)
      if (rne.le.rdat1(np1r-1+1)) then
       in=1
      else
       in=0
 5     in=in+1
       if (rne.ge.rdat1(np1r-1+in+1).and.in.lt.nden) goto 5
      endif
      if (rte.lt.rdat1(np1r-1+nden+1)) then
       it=1
      else
       it=0
 6     it=it+1
       if (it.ge.ntem) then
        it=ntem-1
       else
        if (rte.ge.rdat1(np1r-1+nden+it+1)) goto 6
       endif
      endif
c
      kt1=nden+ntem+(in-1)*ntem+it
      div=rdat1(np1r-1+nden+it+1)-rdat1(np1r-1+nden+it)
      rm=(rdat1(np1r-1+kt1+1)-rdat1(np1r-1+kt1))/(div+1.e-36)
      rec1=rdat1(np1r-1+kt1)+rm*(rte-rdat1(np1r-1+nden+it))
      kt1=kt1+ntem
      rm=(rdat1(np1r-1+kt1+1)-rdat1(np1r-1+kt1))
     $    /(rdat1(np1r-1+nden+it+1)-
     #    rdat1(np1r-1+nden+it)+1.e-36)
      rec2=rdat1(np1r-1+kt1)+rm*(rte-rdat1(np1r-1+nden+it))
c
      rm=(rec2-rec1)/(rdat1(np1r-1+in+1)-rdat1(np1r-1+in)+1.e-36)
      rec=rec1+rm*(rne-rdat1(np1r-1+in))
      if (lpri.ne.0) write (lun11,*)'in calt77:',
     $ temp,den,nden,ntem,rte,rne,wav,in,it,div,rm,
     $ rec1,rec2,rec
      cul=10.**rec
c
      nll=idat1(np1i-1+3)
      k=0
  7   k=k+1
      nl1=k*(k-1)/2+1
      nl2=(k+1)*k/2+1
      if (nll.ge.nl2) goto 7
      il=nll-nl1
      gg=float(2*il+1)*2.
      xt=1.43817e+8/wav/temp
      if (xt.lt.100) then
       clu=cul*exp(-xt)/gg
      else
       clu=0.e0
      endif
c
c
      return
      end
      subroutine chisq(ajisb,cjisb,indb,nindb,
     $   ipmat,x,lun11,lpri)
c
      include './PARAM'
c
      real*8 ajisb(2,ndb),cjisb(ndb)
      integer indb(2,ndb),ll1mx(nd),ll2mx(nd)
      real*8 x(nd),bmat(nd),xo(nd),sum(nd),
     $       term1mx(nd),term2mx(nd),term1,term2
c
c        check the solution
         write (lun11,*)'checking the solution'
         err=0.
        do mm=1,ipmat
          sum(mm)=0.
          term1mx(mm)=0.
          term2mx(mm)=0.
          ll1mx(mm)=0
          ll2mx(mm)=0
          enddo
        write (lun11,*)'rate equations'
        do ll=1,nindb
          mm=min(ipmat,indb(1,ll))
          nn=min(ipmat,indb(2,ll))
          if (mm.ne.nn) then
            term1=ajisb(1,ll)*x(nn)
            term2=ajisb(2,ll)*x(mm)
            sum(mm)=sum(mm)+term1-term2
            if (term1.gt.term1mx(mm)) then
               term1mx(mm)=term1
               ll1mx(mm)=ll
               endif
            if (term2.gt.term2mx(mm)) then
               term2mx(mm)=term2
               ll2mx(mm)=ll
               endif
            write (lun11,*)ll,mm,nn,ajisb(1,ll),ajisb(2,ll),
     $         x(nn),x(mm),sum(mm),term1,term2
 9246       format (1h ,i4,4e12.4,2i4)
            endif
          enddo
        write (lun11,*)'sums'
        do mm=1,ipmat
          write (lun11,*)mm,x(mm),sum(mm),term1mx(mm),
     $       indb(2,ll1mx(mm)),term2mx(mm),indb(2,ll2mx(mm))
          enddo
c
         if (lpri.gt.2)
     $    write (lun11,*)'leaving leqt'
c
      return
      end
      real*8 function cmpfnc(decomp,ecomp,sxcomp,ee,sxx,lun11,lpri)
c
c     this routine is used in the relativistic compton calculation
c
      implicit none
c
      include './PARAM'
c
      real*8 decomp(ncomp,ncomp),ecomp(ncomp),sxcomp(ncomp)
      real*8 ee,eetp,sxx,ddedsx,ddede,dele,delsx,sxtp
      integer lun11,lpri,mm,ll,mmm1,llm1,nc2
c
c     Not used
      integer javi

      javi=lpri
      lpri=javi
c
      eetp=ee
      sxtp=sxx
      nc2=ncomp
      cmpfnc=0.
      if (eetp.gt.1.e-4) then
        call hunt3(ecomp,nc2,eetp,mm,0,lun11)
        call hunt3(sxcomp,nc2,sxtp,ll,0,lun11)
c       if ((mm.gt.1).and.(ll.gt.1)) then
               mm=max(2,min(ncomp,mm))
               ll=max(2,min(ncomp,ll))
               mmm1 = mm - 1
               llm1 = ll - 1
               ddedsx = (decomp(ll,mm)-decomp(llm1,mm)
     $                  +decomp(ll,mmm1)-decomp(llm1,mmm1))
     $                     /(2.*(sxcomp(ll)-sxcomp(llm1)))
               ddede = (decomp(ll,mm)-decomp(ll,mmm1)
     $                  +decomp(llm1,mm)-decomp(llm1,mmm1))
     &                 /(2.*(ecomp(mm)-ecomp(mmm1)))
               dele = ee - ecomp(mmm1)
               delsx = sxx - sxcomp(llm1)
               cmpfnc = ddedsx*delsx + ddede*dele + decomp(llm1,mmm1)
        else
               cmpfnc = 4.*sxx - ee
        endif
c      if (lpri.ne.0)
c     $ write (lun11,*)'in cmpfnc',ee,sxx,ll,mm,cmpfnc
c
      return
      end
      subroutine comp(lpri,lun11,epi,ncn2,bremsa,r,cmp1,cmp2)
c
c
c
c
c     this subroutine computes the heating - cooling due to compton
c     scattering.  the rate is returned in the common block coheat.
c
c
      parameter (ncn=10000)
c
      real*8 bremsa(ncn),epi(ncn),r,cmp1,cmp2
c
c
      data c2/8.219e-06/
c
      lprisv=lpri
c      lpri=3
      if (lpri.ge.1) write (lun11,*)'in comp'
c
      sigth = 6.65e-25
      c1=1.95639e-6
      tmp1 = 0.
      tmp2 = 0.
      c2 = 0.
      r19=r/1.e+19
      fpr2=r19*r19
c
c     due to continuum.
      fac1 = sigth*bremsa(1)*epi(1)*(1.-c2*epi(1))
      fac3 = sigth*bremsa(1)*4.
      numcon=ncn2
      do 100 i = 2,numcon
         delt = epi(i) - epi(i-1)
         fac2 = sigth*bremsa(i)*epi(i)*(1.-c2*epi(i))
         tmp1 = tmp1 + (fac1+fac2)*delt/2.
         fac1 = fac2
         fac4 = sigth*bremsa(i)*4.
         tmp2 = tmp2 + (fac3+fac4)*delt/2.
         fac3 = fac4
         if ( lpri.gt.2 ) write (lun11,99001) i,epi(i),bremsa(i),
     &                           fac1,fac3,tmp1,tmp2
 100  continue
c
      ebar = tmp1*4./(1.e-30+tmp2)
      if ( lpri.gt.2 ) write (lun11,*) 'ebar=',ebar
c
c
      if (lpri.gt.2)  write (lun11,*)c1,tmp1,tmp2
      cmp1 = c1*tmp1
      cmp2 = c1*tmp2
      if (lpri.gt.2) write (lun11,*)cmp1,cmp2
      lpri=lprisv
c
      return
99001 format (' ',i4,6e12.4)
      end
      subroutine comp3(lpri,lun11,epi,ncn2,bremsa,r,cmp1,cmp2)
c
c
c     this subroutine computes the heating - cooling due to compton
c     scattering.  the rate is returned in the common block coheat.
c      using ferland expression
c
      parameter (ncn=10000)
c
      real*8 bremsa(ncn),epi(ncn),r,cmp1,cmp2
c
c
      data c2/8.219e-06/
c
      lprisv=lpri
c      lpri=3
      if (lpri.ge.1) write (lun11,*)'in comp'
c
      sigth = 6.65e-25
      c1=1.95639e-6
      tmp1 = 0.
      tmp2 = 0.
      c2 = 0.
c     trying the old expression
c      c2 = 21./5./5.11e+5
      r19=r/1.e+19
      fpr2=r19*r19
c
c     due to continuum.
      ery=epi(1)/13.6
      alpha=1./(1.+ery*(1.1792e-4+7.084e-10*ery))
      beta=1.-alpha*ery*(1.1792e-4+2.*7.084e-10*ery)/4.
      fac1 = sigth*bremsa(1)*epi(1)*alpha
      fac3 = sigth*bremsa(1)*4.*alpha*beta
      numcon=ncn2
      do 100 i = 2,numcon
         ery=epi(i)/13.6
         alpha=1./(1.+ery*(1.1792e-4+7.084e-10*ery))
         beta=1.-alpha*ery*(1.1792e-4+2.*7.084e-10*ery)/4.
         delt = epi(i) - epi(i-1)
         fac2 = sigth*bremsa(i)*epi(i)*alpha
         tmp1 = tmp1 + (fac1+fac2)*delt/2.
         fac1 = fac2
         fac4 = 4.*sigth*bremsa(i)*alpha*beta
         tmp2 = tmp2 + (fac3+fac4)*delt/2.
         fac3 = fac4
         if ( lpri.gt.2 ) write (lun11,99001) i,epi(i),bremsa(i),
     &                           fac1,fac3,tmp1,tmp2
 100  continue
c
      ebar = tmp1*4./(1.e-30+tmp2)
      if ( lpri.gt.2 ) write (lun11,*) 'ebar=',ebar
c
c
      if (lpri.gt.2)  write (lun11,*)c1,tmp1,tmp2
      cmp1 = c1*tmp1
      cmp2 = c1*tmp2
      if (lpri.gt.2) write (lun11,*)cmp1,cmp2
      lpri=lprisv
c
      return
99001 format (' ',i4,6e12.4)
      end
      subroutine comp2(lpri,lun11,epi,ncn2,bremsa,t,
     $  r,decomp,ecomp,sxcomp,cmp1,cmp2)
c
c     this sub-routine computes the heating - cooling due to compton
c     scattering.  the rate is returned in the common block coheat.
c     relativistic version, using rates from Guilbert (1986)c
c     author:  T. Kallman (from xstar1)
c
      implicit none
c

      include './PARAM'
c
      real*8 decomp(ncomp,ncomp),ecomp(ncomp),sxcomp(ncomp)
      real*8 bremsa(ncn),epi(ncn)
      real*8 t,r,cmp1,cmp2,emc2,sigth0,tmp1,
     $     ekt,xx,sxx,zrmstp,eee,ee,sum1,sum2,sum3,tmp1o,eeeo,
     $     eeo,ans,cfake,hfake,cohc,cmpfnc
      integer lpri,lun11,ncn2,lprisv,numcon,kl
c
c     Not used
      real*8 javir
c
      data emc2/5.11e+5/,sigth0/6.65e-25/
c
      javir=r
c      r=javir
c
      lprisv=lpri
c      lpri=2
      if (lpri.ge.1) write (lun11,*)'in comp2'
c
      sigth0 = 6.65e-25
      tmp1 = 0.
c
      ekt = t*0.861707
      xx = emc2/(ekt+1.e-10)
      sxx = 1./xx
      zrmstp = bremsa(1)
      eee = epi(1)
      ee = eee/emc2
      tmp1 = zrmstp*cmpfnc(decomp,ecomp,sxcomp,ee,sxx,lun11,lpri)
      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      numcon=ncn2
      do kl = 2,numcon
         tmp1o = tmp1
         eeeo = eee
         eeo = ee
         eee = epi(kl)
         ee = eee/emc2
         zrmstp = bremsa(kl)
         tmp1 = zrmstp*cmpfnc(decomp,ecomp,sxcomp,ee,sxx,lun11,lpri)
         sum1 = sum1 + (tmp1+tmp1o)*(eee-eeeo)/2.
         sum2 = sum2 + (bremsa(kl)+bremsa(kl-1))*(eee-eeeo)/2.
         sum3 = sum3 + (bremsa(kl)*ee+bremsa(kl-1)*eeo)*(eee-eeeo)/2.
         if (lpri.ne.0) write (lun11,*)kl,eee,ee,zrmstp,tmp1,sum1
         enddo
      ans = sum1
      cfake=sum2*sigth0
      hfake=sum3*sigth0
      cohc = -ans*sigth0
      cmp1=hfake
      cmp2=(-cohc+hfake)/ekt


      if (lpri.ne.0)
     $ write (lun11,*)'cmp1,cmp2:',cmp1,cmp2,cfake,hfake
c
      lpri=lprisv
c
      return
      end
      subroutine deleafnd(jkk,lup,ml,
     $   nrdti,np1ri,nidti,np1ii,nkdti,np1ki,
     $   idat1,rdat1,kdat1,nptrs,np2,
     $   npar,npnxt,npfi,npfirst,
     $   nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $   npconi2,ncsvn,delea,lfnd,lpri,lun11)
c
      implicit none
c
      include './PARAM'
c
      character(49) kdesc(ntyp),kdesc2
      character(29) krdesc(ntyp)
c     master data
      integer idat1(nidat1),nptrs(nptt,ndat2)
      integer np2
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
      integer jkk,mllz,lup,lfnd,iltmp,lcon,ltyp,lrtyp,
     $   nrdt,np1r,nidt,np1i,nkdt,np1k,ml,nilin,
     $   nrdti,np1ri,nidti,np1ii,nkdti,np1ki,
     $   nlsvn,ncsvn,lpri,lun11,nitmp,iion,mlm,ndtmp
      real*8 delea
c
c     find associated type 86 data
c
c     this is not needed
      go to 9092
      iion=0
c      nilin=npar(ml)
      nilin=npfirst(13)
      mlm=nilin-1
      call drd(ltyp,lrtyp,lcon,
     $           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $           nptrs,0,lun11)
      if (lpri.ge.1)
     $   write (lun11,*)'searching for ion',jkk,mlm,
     $   idat1(np1i+nidt-1)
       do while ((idat1(np1i+nidt-1).ne.jkk)
     $     .and.(iion.lt.nni))
        iion=iion+1
        nitmp=npfi(13,iion)
        mlm=nitmp-1
        call drd(ltyp,lrtyp,lcon,
     $           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $           nptrs,0,lun11)
        if (lpri.ge.1)
     $    write (lun11,*)iion,idat1(np1i+nidt-1),
     $           nitmp
        enddo
 9092   continue
      iion=jkk
c     get damping parameter for iron auger damped lines
      ndtmp=npfi(41,iion)
      if (lpri.gt.1) write (lun11,*)'ndtmp=',iion,ndtmp
      if (ndtmp.eq.0) then
          lfnd=0
          return
        else
          if (lpri.gt.1)
     $    write (lun11,*)'  found ion',lup,ndtmp
          mllz=npar(ndtmp)
          mlm=ndtmp-1
          call drd(ltyp,lrtyp,lcon,
     $           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $           nptrs,0,lun11)
          iltmp=idat1(np1i+1)
          do while ((ndtmp.ne.0).and.(lup.ne.iltmp)
     $            .and.(npar(ndtmp).eq.mllz))
            mlm=ndtmp-1
            call drd(ltyp,lrtyp,lcon,
     $             nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $             nptrs,0,lun11)
            if (lpri.ge.2)
     $        call dprinto(ltyp,lrtyp,lcon,
     $                   nrdt,np1r,nidt,np1i,nkdt,np1k,
     $                   rdat1,idat1,kdat1,lun11)
            iltmp=idat1(np1i+1)
            if (lpri.gt.1)
     $       write (lun11,*)'   ',nidt,
     $       iltmp,ndtmp,lup
            ndtmp=npnxt(ndtmp)
            enddo
        endif
      if (lpri.gt.1) write (lun11,*)'lup,iltmp',
     $                     lup,iltmp
      if (lup.eq.iltmp) then
          lfnd=1
          delea=rdat1(np1r+2)*(4.136e-15)
          if (lpri.gt.1) write (lun11,*)rdat1(np1r+2),delea
        else
          lfnd=0
        endif

c
      return
      end
      subroutine deletefile(filename,status)
c
c     a simple little routine to delete a fits file
c     author:  T. Bridgman
c
      implicit none
      integer status,unit,blocksize
      character*(*) filename
      character*50 kcom
      character*1 ktmp
      integer mm,ll,lenact
c
      data kcom/'rm -f                                             '/
c
      mm=lenact(filename)

      do ll=1,mm
        read (filename(ll:ll),'(a1)')ktmp
        write (kcom(6+ll:6+ll),'(a1)')ktmp
        enddo
      do ll=mm+7,50
        write (kcom(ll:ll),'(a1)')' '
        enddo
c      write (6,*)'executing:',kcom
c      this is slow but it works
c      call system(kcom)
c      return
c
c     simply return if status is greater than zero
c      write (6,*)'in deletefile',unit,filename,status
      if (status .gt. 0)return
c
c     get an unused logical unit number to use to open the fits file
c      call ftgiou(unit,status)
      call getlun(unit)
c      write (6,*)'after ftgiou',unit,status
c
c     try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)
c      write (6,*)'after ftopen',unit,status
c
      if (status .eq. 0)then
c         file was opened;  so now delete it
          call ftdelt(unit,status)
c      write (6,*)'after ftdelt 1',unit,status
      else if (status .eq. 103)then
c         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
      else
c         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
c      write (6,*)'after ftdelt 2',unit,status
      end if


c     free the unit number for later reuse
c      call ftfiou(unit, status)
      close(unit)
c      write (6,*)'after ftfiou',unit,status
c
      end
      subroutine dfact(n,x)
c
c to calculate the factorial of an integer n.  the output x is
c the natural log of n factorial, in double precision.
c     author:  T. Kallman
c
      implicit none
      integer n,i
      real*8 x
c
      x=0.
      if(n.eq.0) return
      do i=1,n
       x=x+log(float(i))
       enddo
c
      return
      end
      subroutine dprint(ltyp,lrtyp,lcon,
     $  lrdat,rdat,lidat,idat,lkdat,kdat,
     $  np1r,np1i,np1k,np2,
     $  idat1,rdat1,kdat1,nptrs,lpri,lun11)
c
c     this  routine prints one element of the database
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
      integer nptmpdim
      parameter (nptmpdim=200000)
c
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
      integer nptrs(nptt,ndat2)
      integer ltyp,lrtyp,lcon,lrdat,lidat,lkdat,
     $  np1r,np1i,np1k,np2,ml,lpri,lun11,lprisv
      real*8 rdat(nptmpdim)
      integer idat(nptmpdim)
      character(1) kdat(nptmpdim)
cc
      lprisv=lpri
      if (lpri.ne.0)
     $ write (lun11,*)'in dprint, np2=',np2
      if (np2.ge.ndat2) then
          write (6,*) 'data index error'
          return
          endif
      nptrs(1,np2)=1
      nptrs(2,np2)=ltyp
      nptrs(3,np2)=lrtyp
      nptrs(4,np2)=lcon
      nptrs(5,np2)=lrdat
      nptrs(6,np2)=lidat
      nptrs(7,np2)=lkdat
      nptrs(8,np2)=np1r
      nptrs(9,np2)=np1i
      nptrs(10,np2)=np1k
      if (lpri.ne.0) then
        write (lun11,*)'in dprint:',np2,ltyp,lrtyp,lrdat,lidat,lkdat
        write (lun11,*)'          ',lcon,np1r,np1i,np1k
        endif
      np2=np2+1
      if (lrdat.gt.0) then
        do ml=1,lrdat
           rdat1(np1r)=rdat(ml)
           np1r=np1r+1
           enddo
        endif
      if (lidat.gt.0) then
        do  ml=1,lidat
           idat1(np1i)=idat(ml)
           np1i=np1i+1
           enddo
        endif
      if (lkdat.eq.0) then
        do ml=1,lkdat
            kdat1(np1k)=kdat(ml)
            np1k=np1k+1
            enddo
        endif
c
      if ((np1k.gt.nkdat1).or.(np1i.gt.nidat1).or.(np1r.gt.nrdat1)
     $   .or.(np2.gt.ndat2)) then
        write (lun11,*)'dprint index error,',np1k,np1i,np1r,np2
        return
        endif
c
      lpri=lprisv
c
      return
      end
      subroutine dprinto(ltyp,lrtyp,lcon,
     $  nrdt,np1r,nidt,np1i,nkdt,np1k,rdat1,idat1,kdat1,lun11)
c
c
      implicit none
      include './PARAM'
      integer ltyp,lrtyp,lcon,nrdt,nidt,nkdt,lun11,lsp,ml2,nkd,
     $        nkd2,ll2,ml,mm,ll
      integer nptmpdim
      parameter (nptmpdim=200000)
c
      character(1) kblnk,kperc
      character(400000) kdtt
      character(1) ktst,kdtt2(400000)
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
      integer np1r,np1i,np1k
c
      data kblnk/' '/,kperc/'%'/
c
c      write (lun11,*)ltyp,lrtyp,lcon,nrdt,nidt,nkdt,np1r,np1i,np1k
c     $  (rdat(mm),mm=1,nrdt),(idat1(np1i-1+mm),mm=1,nidt),
c     $  kblnk,(kdat1(np1k-1+mm),mm=1,nkdt),kblnk,kperc
c      return
c
      write (kdtt(1:18),9823)ltyp,lrtyp,lcon
      write (kdtt(19:37),9823)nrdt,nidt,nkdt
c
c
      lsp=0
      ml2=0
      if (lsp.eq.1) go to 9009
      nkd=38
      nkd2=39
      if (nrdt.gt.0) then
c        nkd2=nkd+nrdt*13
        nkd2=nkd
        do 303 ml=1,nrdt
           ml2=nkd+(ml-1)*15
             write (kdtt(ml2:ml2+14),'(1pe15.7)')rdat1(np1r-1+ml)
             nkd2=nkd2+15
 303       continue
        endif
      nkd=nkd2
      write (kdtt(nkd:nkd),'(a1)')kblnk
      nkd=nkd2+1
      if (nidt.gt.0) then
        nkd2=nkd+nidt*8
        do 302 ml=1,nidt
           ml2=nkd+(ml-1)*8
           write (kdtt(ml2:ml2+7),'(i8)')idat1(np1i-1+ml)
           ml2=ml2+8
 302       continue
        endif
      nkd=nkd2
      if (nkdt.gt.0) then
        write (kdtt(nkd:nkd),'(a1)')kblnk
        nkd=nkd+1
        nkd2=nkd+nkdt
c        write (lun11,*)nkd,nkdt,nkd2,(kdat1(np1k-1+mm),mm=1,nkdt)
        do 301 ml=1,nkdt
          ml2=nkd+ml-1
          write (kdtt(ml2:ml2),'(a1)')kdat1(np1k-1+ml)
           ml2=ml2+1
 301      continue
         endif
 9823    format (3i6)
c       write (lun11,*)'before write:'
c       write (lun11,*)kdtt
      ml2=ml2-1
c
c
      ll2=0
c     remove spaces
      ktst=kperc
      do 3301 ll=1,ml2
c         ktsto=ktst
         read(kdtt(ll:ll),'(a1)')ktst
c         if ((ktst.eq.kblnk).and.(ktsto.eq.kblnk)) go to 3301
         ll2=ll2+1
         kdtt2(ll2)=ktst
 3301    continue
c
      write (lun11,911)(kdtt2(mm),mm=1,ll2),kblnk,kperc
 911  format (400000a1)
c
      return
c
 9009 continue
c      call dprints(ltyp,lrtyp,lcon,
c     $  nrdt,rdat,nidt,idat,nkdt,kdat,lun11)
c
c
      return
      end
      subroutine dprinto2(ltyp,lrtyp,lcon,
     $  nrdt,rdat,nidt,idat,nkdt,kdat,lun11)
c
c     this  routine prints one element of the database
c     author:  T. Kallman
c
c
      implicit none
c
      integer nptmpdim
      parameter (nptmpdim=200000)
c
      real*8 rdat(nptmpdim)
      integer idat(nptmpdim)
      character(1) kblnk,kperc,kdat(nptmpdim)
      integer ltyp,lrtyp,lcon,nrdt,nidt,nkdt,lun11,lsp,ml2,nkd,
     $        nkd2,ll2,ml,mm,ll
c
      character(400000) kdtt
      character(1) ktst,kdtt2(400000)
c
      data kblnk/' '/,kperc/'%'/
c
c      write (lun11,*)ltyp,lrtyp,lcon,nrdt,nidt,nkdt,np1r,np1i,np1k
c     $  (rdat(mm),mm=1,nrdt),(idat1(np1i-1+mm),mm=1,nidt),
c     $  kblnk,(kdat1(np1k-1+mm),mm=1,nkdt),kblnk,kperc
c      return
c
      write (kdtt(1:18),9823)ltyp,lrtyp,lcon
      write (kdtt(19:37),9823)nrdt,nidt,nkdt
c
c
      lsp=0
      ml2=0
      if (lsp.eq.1) go to 9009
      nkd=38
      nkd2=39
      if (nrdt.gt.0) then
c        nkd2=nkd+nrdt*13
        nkd2=nkd
        do 303 ml=1,nrdt
           ml2=nkd+(ml-1)*15
             write (kdtt(ml2:ml2+14),'(1pe15.7)')rdat(ml)
             nkd2=nkd2+15
 303       continue
        endif
      nkd=nkd2
      write (kdtt(nkd:nkd),'(a1)')kblnk
      nkd=nkd2+1
      if (nidt.gt.0) then
        nkd2=nkd+nidt*8
        do 302 ml=1,nidt
           ml2=nkd+(ml-1)*8
           write (kdtt(ml2:ml2+7),'(i8)')idat(ml)
           ml2=ml2+8
 302       continue
        endif
      nkd=nkd2
      if (nkdt.gt.0) then
        write (kdtt(nkd:nkd),'(a1)')kblnk
        nkd=nkd+1
        nkd2=nkd+nkdt
c        write (lun11,*)nkd,nkdt,nkd2,(kdat1(np1k-1+mm),mm=1,nkdt)
        do 301 ml=1,nkdt
          ml2=nkd+ml-1
          write (kdtt(ml2:ml2),'(a1)')kdat(ml)
           ml2=ml2+1
 301      continue
         endif
 9823    format (3i6)
c       write (lun11,*)'before write:'
c       write (lun11,*)kdtt
      ml2=ml2-1
c
c
      ll2=0
c     remove spaces
      ktst=kperc
      do 3301 ll=1,ml2
c         ktsto=ktst
         read(kdtt(ll:ll),'(a1)')ktst
c         if ((ktst.eq.kblnk).and.(ktsto.eq.kblnk)) go to 3301
         ll2=ll2+1
         kdtt2(ll2)=ktst
 3301    continue
c
      write (lun11,911)(kdtt2(mm),mm=1,ll2),kblnk,kperc
 911  format (400000a1)
c
      return
c
 9009 continue
c      call dprints(ltyp,lrtyp,lcon,
c     $  nrdt,rdat,nidt,idat,nkdt,kdat,lun11)
c
c
      return
      end
      subroutine dprints(ltyp,lrtyp,lcon,
     $  nrdt,rdat,nidt,idat,nkdt,kdat,lun11)
c
c     this  routine prints one element of the database
c     author:  T. Kallman
c
      implicit none
      integer nptmpdim
      parameter (nptmpdim=200000)
      real*8 rdat(nptmpdim)
      integer idat(nptmpdim)
      character(1) kdat(nptmpdim)
      character(20000) kdtt
      character(1) kblnk,ktst,kperc,kdtt2(nptmpdim)
      integer ltyp,lrtyp,lcon,nrdt,nidt,nkdt,lun11,
     $        nkd,nkd2,ml2,itmp,ml,ll2,ll,mm
      real*8 rtmp
c
      data kblnk/' '/,kperc/'%'/
c
c      do ll=1,20000
c         write (kdtt(ll:ll),'(a1)')kblnk
c         enddo
      write (kdtt(1:18),9823)ltyp,lrtyp,lcon
      write (kdtt(19:37),9823)nrdt,nidt,nkdt
c
      nkd=38
      nkd2=37
      rtmp=rdat(1)
      if (1.gt.nrdt) rtmp=0.
      ml2=nkd2
      write (kdtt(ml2:ml2+12),'(1pe13.5)')rtmp
      ml2=nkd+13
      write (kdtt(ml2-1:ml2-1),'(a1)')kblnk
      rtmp=rdat(2)
c      rtmp=rdat(nrdt)
      if (2.gt.nrdt) rtmp=0.
      write (kdtt(ml2:ml2+12),'(1pe13.5)')rtmp
      nkd2=nkd+2*13
      nkd=nkd2
c
      write (kdtt(nkd:nkd),'(a1)')kblnk
      nkd=nkd2+1
      ml2=nkd
      itmp=idat(1)
      if (1.gt.nidt) itmp=0
      ml2=nkd2
      write (kdtt(ml2:ml2+5),'(i6)') itmp
      ml2=ml2+6
      itmp=0
      if (nidt.gt.1) itmp=idat(nidt-1)
      write (kdtt(ml2:ml2+5),'(i6)') itmp
      ml2=ml2+6
      itmp=idat(nidt)
      if (1.gt.nidt) itmp=0
      write (kdtt(ml2:ml2+5),'(i6)') itmp
      nkd2=nkd+3*6-1
        write (kdtt(nkd:nkd),'(a1)')kblnk
      nkd=nkd2
c
      nkd=nkd2
      ml2=nkd
      if (nkdt.ne.0) then
        write (kdtt(nkd:nkd),'(a1)')kblnk
        nkd=nkd+1
        nkd2=nkd+nkdt
c        write (lun11,*)nkd,nkdt,nkd2,(kdat(mm),mm=1,nkdt)
        do  ml=1,nkdt
          ml2=nkd+ml-1
          write (kdtt(ml2:ml2),'(a1)')kdat(ml)
           ml2=ml2+1
          enddo
        endif
 9823    format (4i6)
c       write (lun11,*)'before write:'
c       write (lun11,*)kdtt
      ml2=ml2-1
c
      ll2=0
c     remove spaces
      ktst=kperc
      do ll=1,ml2
         read(kdtt(ll:ll),'(a1)')ktst
         ll2=ll2+1
         kdtt2(ll2)=ktst
         enddo
c
      write (lun11,911)(kdtt2(mm),mm=1,ll2),kblnk,kperc
 911  format (20000a1)
c
c
c
      return
      end
      subroutine dprints2(ltyp,lrtyp,lcon,
     $  nrdt,rdat,nidt,idat,nkdt,kdat,lun11)
c
c     this  routine prints one element of the database
c     author:  T. Kallman
c
      implicit none
      integer nptmpdim
      parameter (nptmpdim=200000)
      real*8 rdat(nptmpdim)
      integer idat(nptmpdim)
      character(1) kdat(nptmpdim)
      character(20000) kdtt
      character(1) kblnk,ktst,kperc,kdtt2(nptmpdim)
      integer ltyp,lrtyp,lcon,nrdt,nidt,nkdt,lun11,
     $        nkd,nkd2,ml2,itmp,ml,ll2,ll,mm
      real*8 rtmp
c
      data kblnk/' '/,kperc/'%'/
c
c      do ll=1,20000
c         write (kdtt(ll:ll),'(a1)')kblnk
c         endif
      write (kdtt(1:18),9823)ltyp,lrtyp,lcon
      write (kdtt(19:37),9823)nrdt,nidt,nkdt
c
      nkd=38
      nkd2=37
      rtmp=rdat(1)
      if (1.gt.nrdt) rtmp=0.
      ml2=nkd2
      write (kdtt(ml2:ml2+12),'(1pe13.5)')rtmp
      ml2=nkd+13
      write (kdtt(ml2-1:ml2-1),'(a1)')kblnk
      rtmp=0.
      if (2.le.nrdt) rtmp=rdat(nrdt)
      write (kdtt(ml2:ml2+12),'(1pe13.5)')rtmp
      nkd2=nkd+2*13
      nkd=nkd2
c
      write (kdtt(nkd:nkd),'(a1)')kblnk
      nkd=nkd2+1
      ml2=nkd
      itmp=idat(1)
      if (1.gt.nidt) itmp=0
      ml2=nkd2
      write (kdtt(ml2:ml2+5),'(i6)') itmp
      ml2=ml2+6
      itmp=idat(nidt)
      if (1.gt.nidt) itmp=0
      write (kdtt(ml2:ml2+5),'(i6)') itmp
      nkd2=nkd+2*6-1
        write (kdtt(nkd:nkd),'(a1)')kblnk
      nkd=nkd2
c
      nkd=nkd2
      ml2=nkd
      if (nkdt.ne.0) then
        write (kdtt(nkd:nkd),'(a1)')kblnk
        nkd=nkd+1
        nkd2=nkd+nkdt
c        write (lun11,*)nkd,nkdt,nkd2,(kdat(mm),mm=1,nkdt)
        do  ml=1,nkdt
          ml2=nkd+ml-1
          write (kdtt(ml2:ml2),'(a1)')kdat(ml)
           ml2=ml2+1
          enddo
        endif
 9823    format (4i6)
c       write (lun11,*)'before write:'
c       write (lun11,*)kdtt
      ml2=ml2-1
c
      ll2=0
c     remove spaces
      ktst=kperc
      do ll=1,ml2
         read(kdtt(ll:ll),'(a1)')ktst
         ll2=ll2+1
         kdtt2(ll2)=ktst
         enddo
c
      write (lun11,911)(kdtt2(mm),mm=1,ll2),kblnk,kperc
 911  format (20000a1)
c
      return
      end
      subroutine drd(ltyp,lrtyp,lcon,lrdat,np1r,lidat,np1i,lkdat,np1k,
     &                 np2,nptrs,lpri,lun11)
c
c     this routine reads one element from the database
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
c
c     master data
      integer  nptrs(nptt,ndat2)
      integer nrd,lcon,ltyp,lrtyp,lrdat,lidat,
     $        lkdat,np2,lpri,lun11,np1,np1r,np1i,np1k
c
      if ( lpri.ne.0 ) write (lun11,*) 'in drd, np2=' , np2,
     &                                  ntyp
c      if ((ltyp.le.0).or.(ltyp.gt.ntyp))
c     $    stop 'data typing error'
      nrd = 0
      lcon=1
        np2 = np2 + 1
        nrd = nrd + 1
        np1 = nptrs(1,np2)
        ltyp = nptrs(2,np2)
        lrdat = nptrs(5,np2)
        lidat = nptrs(6,np2)
        lkdat = nptrs(7,np2)
        lrtyp = nptrs(3,np2)
        lcon = nptrs(4,np2)
        np1r = nptrs(8,np2)
        np1i = nptrs(9,np2)
        np1k = nptrs(10,np2)
        if ( lpri.ne.0 ) write (lun11,*) 'in dread:' , np2 , np1 ,ltyp,
     &                                 lrtyp , lrdat , lidat
        if ( lpri.ne.0 ) write (lun11,99001) lkdat , lcon , np1r ,np1i,
     &                        np1k
      lcon = 0
      if ( lpri.ne.0 ) write (lun11,*) 'leaving drd' , np2
c
      return
99001 format (8x,5i8)
      end
      subroutine dread(ltyp,lrtyp,lcon,lrdat,rdat,lidat,idat,lkdat,kdat,
     &                 np2,idat1,rdat1,kdat1,nptrs,lpri,lun11)
c
c     this routine reads one element from the database
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=200000)
c
c     master data
      integer idat1(nidat1) , nptrs(nptt,ndat2)
      real*8  rdat1(nrdat1)
      character(1) kdat1(nkdat1)
      real*8 rdat(nptmpdim)
      integer idat(nptmpdim)
      character(1) kdat(nptmpdim)
      integer mlr,mli,mlk,nrd,lcon,ltyp,lrtyp,lrdat,lidat,
     $        lkdat,np2,lpri,lun11,np1,np1r,np1i,np1k,ml,mm
c
      if ( lpri.ne.0 ) write (lun11,*) 'in dread, np2=' , np2,
     &                                  ntyp
c      if ((ltyp.le.0).or.(ltyp.gt.ntyp))
c     $    stop 'data typing error'
      mlr = 0
      mli = 0
      mlk = 0
      nrd = 0
      lcon=1
      do while (lcon.ne.0)
        np2 = np2 + 1
        nrd = nrd + 1
        np1 = nptrs(1,np2)
        ltyp = nptrs(2,np2)
        lrdat = nptrs(5,np2)
        lidat = nptrs(6,np2)
        lkdat = nptrs(7,np2)
        lrtyp = nptrs(3,np2)
        lcon = nptrs(4,np2)
        np1r = nptrs(8,np2)
        np1i = nptrs(9,np2)
        np1k = nptrs(10,np2)
        if ( lpri.ne.0 ) write (lun11,*) 'in dread:' , np2 , np1 ,ltyp,
     &                                 lrtyp , lrdat , lidat
        if ( lpri.ne.0 ) write (lun11,99001) lkdat , lcon , np1r ,np1i,
     &                        np1k
        if ( lrdat.ne.0 ) then
          do ml = 1 , lrdat
            rdat(mlr+ml) = rdat1(np1r+ml-1)
            if ( lpri.ne.0 ) write (lun11,*) mlr , np1r , rdat1(np1r) ,
     &                              rdat(mlr+ml)
            enddo
          np1r=np1r+lrdat-1
          mlr=mlr+lrdat
          if ( lpri.ne.0 ) write (lun11,*) 'rdat=' ,
     &                           (rdat(mm),mm=1,lrdat) , np2
          endif
        if ( lidat.ne.0 ) then
          do ml = 1 , lidat
            idat1(np1i-1+mli+ml) = idat1(np1i+ml-1)
            if ( lpri.ne.0 ) write (lun11,*) mli , np1i , idat1(np1i) ,
     &                              idat1(np1i-1+mli+ml) , np2
            enddo
          mli = mli + lidat
          np1i = np1i + lidat-1
          if ( lpri.ne.0 ) write (lun11,*) 'idat=' ,
     &                           (idat1(np1i-1+mm),mm=1,lidat)
          endif
        if ( lkdat.ne.0 ) then
          do ml = 1 , lkdat
            kdat(mlk+ml) = kdat1(np1k+ml-1)
c           write (lun11,*)mlk,np1k,kdat1(np1k),kdat(mlk),np2
            enddo
          mlk = mlk + lkdat
          np1k = np1k + lkdat-1
          if ( lpri.ne.0 ) write (lun11,*) 'kdat=' ,
     &                           (kdat(mm),mm=1,lkdat)
          endif
        enddo
c      np2=np2-nrd+1
      lidat = mli
      lrdat = mlr
      lkdat = mlk
      lcon = 0
c
c     the last pointer has to point to the next empty space.
c      nptr1(np2+1)=np1+1
c
c      call remtms(tt1)
c      tread = tread + abs(tt1-tt0)
c
      if ( lpri.ne.0 ) write (lun11,*) 'leaving dread' , np2
c
      return
99001 format (8x,5i8)
      end
      subroutine dsec(lnerr,nlim,
     $       lpri,lppri,lun11,tinf,vturbi,critf,
     $       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,
     $       zrems,zremso,
     $       elumab,elumabo,elum,elumo,
     $       decomp,ecomp,sxcomp,
     $       tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,nlevs,ncsvn,rates,vsav,idrates,
     $       ntotit,
     $       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,
     $         cllines,clcont,htcomp,clcomp,clbrems,
     $       xilev,bilev,rniss,nmat,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,fline,flinel, t_min, mhd_heat)
c
c     this routine solves for temperature and electron density by the
c     double secant method
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
c
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     line luminosities
      real*8 elum(3,nnnl),elumo(3,nnnl)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
      real*8 fline(2,nnnl),flinel(ncn)
c     line optical depths
      real*8 tau0(2,nnnl)
c     energy bins
      real*8 epi(ncn)
c      continuum lum
      real*8 zrems(4,ncn),zremso(4,ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn),brcems(ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     level populations
      real*8 xilev(nnml),bilev(nnml),rniss(nnml)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
      real*8 elumab(2,nnml),elumabo(2,nnml)
      real*8 tauc(2,nnml)
c     ion abundances
      real*8 xii(nni)
c     heating and cooling
      real*8 htt(nni),cll(nni)
      real*8 rrrt(nni),pirt(nni)
      integer nlevs(nni)
c     element abundances
      real*8 abel(nl)
c     the saved rates
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
c     compton heating data
      real*8 decomp(ncomp,ncomp),ecomp(ncomp),sxcomp(ncomp)
c     state variables
      real*8 p,r,t,xpx,delr
c     heating-cooling variables
      real*8 httot,cltot,htcomp,clcomp,clbrems,elcter,cllines,
     $     clcont,hmctot
c     input parameters
      real*8 trad,tinf
      real*8 cfrac,critf,vturbi,xee
      integer lcdd,ncn2,lpri,lun11,np2,nlim
c     variables associated with thermal equilibrium solution
      integer nmat,ntotit
c     temporary for xwrite
      character(133) tmpst
      integer nlsvn,ncsvn
c
c     local variables
      integer nnt,nntt,lnerr,lppri0,lppri,nlimt,nlimx,nnxx,
     $        nlimtt,nlimxx,iht,ilt,iuht,iult,ihx,ilx,nnx
      real*8 crite,crith,critt,fact,facx,epst,epsx,epstt,to,
     $     tl,th,xeel,xeeh,elctrl,elctrh,hmctth,hmcttl,tst,
     $     testt
c
c     -----------------
      real*8 t_min
      integer been_here
      real*8 mhd_heat
c
      been_here = 0
c     -----------------
c
      if (lpri.ne.0) write (lun11,*)'in dsec'
c
      crite=1.e-05
c     crite=1.e-06
c     crith=1.e-02
      crith=1.e-03
c     crith=1.e-06
      critt=1.e-03
c     critt=1.e-09
c
      ntotit=0
      nnt = 0
      nntt=0
      lnerr = 0
      lppri0 = lppri
      nlimt =max(nlim,0)
      nlimx=abs(nlim)
      nlimtt=max0(nlimt,1)
      nlimxx=max0(nlimx,1)
      if (lpri.ne.0)
     $ write (lun11,*)'nlimtt,nlimxx,lppri--',nlimtt,nlimxx,lppri
      fact = 1.2
      facx = 1.2
      epst = crith
      epsx = crite
      epstt = critt
      to = 1.e+30
      tl = 0.
      th = 0.
      xeel = 0.
      xeeh = 1.
      elctrl = 1.
      elctrh = -1.
      hmctth = 0.
      hmcttl = 0.
c
      iht = 0
      ilt = 0
      iuht = 0
      iult = 0
c
 100  nnx = 0
      t=max(t,tinf)
      if (t.lt.tinf*1.01) then
          nlimt=0
          nlimtt=0
          nlimx=0
          nlimxx=0
        else
          nlimt =max(nlim,0)
          nlimx=abs(nlim)
          nlimxx=nlimx
        endif
c      if (t.lt.tinf) return
      nnxx=0
      ihx = 0
      ilx = 0
 200  continue
      if ( lppri.ne.0 ) then
        write (lun11,99001)
     $   nnx,xee,xeel,xeeh,elcter,elctrl,elctrh,
     $   nnt,t,tl,th,hmctot,hmcttl,hmctth
        write (tmpst,99001)
     $   nnx,xee,xeel,xeeh,elcter,elctrl,elctrh,
     $   nnt,t,tl,th,hmctot,hmcttl,hmctth
        call xwrite(tmpst,10)
        endif
      call func(lpri,lun11,vturbi,critf,
     $       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       decomp,ecomp,sxcomp,
     $       tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,nlevs,ncsvn,rates,vsav,idrates,
     $       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,
     $         cllines,clcont,htcomp,clcomp,clbrems,
     $       xilev,bilev,rniss,nmat,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,fline,flinel, mhd_heat)
      if ( lppri.ne.0 ) then
        write (lun11,99001)
     $   nnx,xee,xeel,xeeh,elcter,elctrl,elctrh,
     $   nnt,t,tl,th,hmctot,hmcttl,hmctth
        write (tmpst,99001)
     $   nnx,xee,xeel,xeeh,elcter,elctrl,elctrh,
     $   nnt,t,tl,th,hmctot,hmcttl,hmctth
        call xwrite(tmpst,10)
99001   format (' in dsec -- ',i4,6(1pe9.2),i4,6(1pe9.2))
        endif
      ntotit=ntotit+1
      nnx = nnx + 1
      nnxx=nnxx+1
      if (nnxx.ge.nlimxx) go to 300
      tst=abs(elcter)/max(1.e-10,xee)
      if (tst.lt.epsx) go to 300
      if ( elcter.lt.0 ) then
            ihx = 1
            xeeh = xee
            elctrh = elcter
            if ( ilx.ne.1 ) then
               xee = xee*facx
c              print *, xee
               goto 200
               endif
         else
            ilx = 1
            xeel = xee
            elctrl = elcter
            if ( ihx.ne.1 ) then
               xee = xee/facx
c              print *, xee
               goto 200
               endif
         endif
         xee = (xeel*elctrh-xeeh*elctrl)/(elctrh-elctrl)
c        print *, xee
         goto 200
c
c
 300  continue
c
c    ------------------------------------------------------------------
c
      if (been_here.gt.4) then
            t = t_min
            goto 500
      endif
      if (t.lt.t_min) then
            print *, "^^^"
            print *, "t = ", t 
            print *, " ... is less than t_min = ", t_min
            t = t_min
            print *, "... now t = ", t
            print *, "^^^"
            been_here = been_here + 1
      endif
      if (t.gt.t_min) been_here = been_here - 1
c
c    ------------------------------------------------------------------
c
      nntt=nntt+1
      nnt = nnt + 1
      if ( abs(hmctot).le.epst )  goto 500
      if (nntt.ge.nlimtt) go to 500
      if ( nnt.lt.nlimt ) then
         if ( hmctot.lt.0 ) then
            iht = 1
            th = t
            hmctth = hmctot
            iuht = 1
            if ( iult.eq.0 ) hmcttl = hmcttl/2.
            iult = 0
            if ( ilt.ne.1 ) then
               t = t/fact
               goto 100
            endif
         else
            ilt = 1
            tl = t
            hmcttl = hmctot
            iult = 1
            if ( iuht.eq.0 ) hmctth = hmctth/2.
            iuht = 0
            if ( iht.ne.1 ) then
               t = t*fact
               goto 100
            endif
         endif
         testt = abs(1.-t/to)
         if ( testt.lt.epstt ) then
            lnerr = -2
            if ( lppri.ne.0 ) then
               write (lun11,99004)
               write (lun11,99006) nnt,t,tl,th,hmctot,hmcttl,
     &                         hmctth
            endif
            goto 500
         else
            to = t
            t = (tl*hmctth-th*hmcttl)/(hmctth-hmcttl)
            goto 100
         endif
      endif
c
      lnerr = 2
      write (lun11,99002)
      write (lun11,99006) nnt,t,tl,th,hmctot,hmcttl,hmctth
c
 500  if ( lppri.ne.0 ) write (lun11,99007) testt,epst,hmctot
      lppri = lppri0
c
      return
99002 format (' ','**** note: in dsec --  too many iterations **** ')
99004 format (' ',' warrning -- dsec not converging ')
99006 format (' ',' temperature ',i4,6(1pe16.8))
99007 format (' ',' finishing dsec -- test,epst,hmctot',3(1pe16.8))
      end
      real*8 function ee1exp(x)
c
c     this routine computes the first exponential integral.
c
      implicit none
c
      real*8 x
c
      if ( x.ge.1. ) then
        ee1exp=(1./x)*(0.250621+x*(2.334733+x))/(1.68153+x*(3.330657+x))
        return
      endif
c
      ee1exp = (-log(x)-0.57721566+
     &      x*(0.99999193+x*(-0.24991055+x*(0.05519968+
     &      x*(-0.00976004+x*0.0010707857)))))*exp(x)
c
      return
      end
      real*8 function ee1expo(x)
c
c     this routine computes the first exponential integral.
c
      implicit none
c
      real*8 x,expo
c
      if ( x.ge.1. ) then
      ee1expo=(1./x)*(0.250621+x*(2.334733+x))/(1.68153+x*(3.330657+x))
      return
      endif
c
      ee1expo = (-log(x)-0.57721566+
     &      x*(0.99999193+x*(-0.24991055+x*(0.05519968+
     &      x*(-0.00976004+x*0.0010707857)))))*expo(x)
c
      return
      end
       subroutine eint(t,e1,e2,e3)
c
c  returns the values of the exponential integral function of order
c  1, 2, and 3
c     author:  T. Kallman
c
       implicit none
       real*8 t,e1,e2,e3,ss,expo
c
       e1=0.
       e2=0.
       e3=0.
c       if (t.gt.50.) return
       call expint(t,ss)
       e1=ss/t/exp(t)
       e2=exp(-t)-t*e1
       e3=0.5*(expo(-t)-t*e2)
       return
       end
      subroutine ener(epi,ncn2)
c
c     sets up energy grin
c     author: T. Kallman
c
      implicit none
c
      include './PARAM'
c
      real*8 epi(ncn)
      integer numcon,numcon2,numcon3,ncn2,ll,ll2
      real*8 ebnd1,ebnd2,ebnd2o,dele
c
      numcon = ncn2
      if (numcon.lt.4) write (6,*) 'in ener: numcon error'
      numcon2=max(2,ncn2/50)
      numcon3=numcon-numcon2
      ebnd1=0.1
c     nb changed energy grid for H only
      ebnd2=4.e+4
c      ebnd2=4.e+1
      ebnd2o=ebnd2
      dele=(ebnd2/ebnd1)**(1./float(numcon3-1))
      epi(1)=ebnd1
      do ll=2,numcon3
        epi(ll)=epi(ll-1)*dele
        enddo
      ebnd2=1.e+6
      ebnd1=ebnd2o
      dele=(ebnd2/ebnd1)**(1./float(numcon2-1))
      do ll2=1,numcon2
        ll=ll2+numcon3
        epi(ll)=epi(ll-1)*dele
        enddo
c
      return
      end
      subroutine enxt(eth,nb1,lpri,epi,ncn2,t,lfast,lun11,
     $                  jk,nskp,nphint,lrcalc)
c
c     finds next energy bin for photoionizaion rate integrations
c     author: T. Kallman
c
      implicit none
c
      include './PARAM'
c
      real*8 epi(ncn)
      real*8 ergsev,bk,tm,t,eth,bktm,exptst,epii
      integer nb1,ncn2,lfast,lun11,jk,nphint,lrcalc,lpri
      integer nskp,numcon2,nbinc,numcon3,nskp1,numcon,nskp2
c
      data ergsev/1.602197e-12/
      data bk/1.38062e-16/
c
       if (lpri.gt.2)
     $  write (lun11,*)'in enxt:',eth,nb1,t,lfast,jk,lpri,
     $                    epi(1),epi(ncn2),ncn2
      tm=t*1.e4
      bktm=bk*tm/ergsev
      if (lfast.le.2) then
         numcon2=max(2,ncn2/50)
         nphint=ncn2-numcon2
         nskp=1
         nskp2=1
      elseif (lfast.eq.3) then
         nphint=nbinc(max(3.*eth,eth+3.*bktm),epi,ncn2)
         nphint=max(nphint,nb1+1)
         nskp=max(1,int((nphint-nb1)/16))
         nskp2=nskp
      else
        nphint=nbinc(1.d+4,epi,ncn2)
        nskp=1
        nskp2=1
        endif
      nskp1=nskp
      epii=epi(jk)
      exptst=(epii-eth)/bktm
      if (exptst.lt.3.) then
         lrcalc=1
         nskp=nskp1
       else
         lrcalc=0
         nskp=nskp2
       endif
       nphint=max(nphint,nb1+nskp)
       numcon=ncn2
       numcon2=max(2,ncn2/50)
       numcon3=numcon-numcon2
       nphint=min(nphint,numcon3)
       if (lpri.gt.2)
     $  write (lun11,*)'in enxt:',eth,nb1,t,lfast,jk,nskp,
     $   nphint,lrcalc
c
      return
      end
      subroutine erc(n,m,t,ic,se,sd,a,lun11,lpri)
c
c erc calculates the excitation rate, se [cm**3/s],  for atomic
c transitions between lower state n and upper state m in hydrogen
c due to electron collisions.  the energy loss rate, sl [ev*cm**3/s],
c from the electron gas is also determined.  (cf. johnson,1972)
c sd = deexcitation rate;   sg = energy gained by electron gas
c sm is a quantity symmetrical in n and m, used in models
c ***  the quantity em1 is required from subr. expint in this program
c     author:  M. Bautista
c
      implicit none
c
      integer n,m,ic,lun11,lpri
      real*8 t,a,ym,xn,s,sd,se,f,z,yn,dif,e1y,e1z,e2,rn,rm,
     $     ann,bnn,sm,expo,bn,ric
c
c
      rn=float(n)
      rm=float(m)
      ric=float(ic)
      if (lpri.gt.1)
     $ write(lun11,*)'erc',n,m,t,ic,a
      sm=0.
      if(ic.ne.1) then
       if (ic.lt.10) then
        ym=157803.*ic*ic/t/m/m
        if (ym.gt.40.) then
         sd=0.
         se=0.
         return
        endif
        if (lpri.gt.1)
     $   write (lun11,*)'before impactn:',
     $       n,m,t,ic,a,sm
        call impactn(n,m,t,ic,a,sm,lun11,lpri)
        if (lpri.gt.1)
     $   write (lun11,*)'after impactn:',
     $       n,m,t,ic,a,sm
        ym=157803.*ric*ric/t/(rm*rm)
        xn=(1./(rn*rn)-1./(rm*rm))
        yn=157803.*ric*ric*xn/t
        s=sm/(rn*rn)/exp(ym)
        sd=s*rn*rn/(rm*rm)
        if (yn.lt.40.) then
         se=s*exp(-yn)
        else
         se=0.e0
        endif
        if (lpri.gt.1)
     $    write (lun11,*)ym,xn,yn,s,sd,se
       else
        if (lpri.gt.1)
     $   write (lun11,*)'calling szcoll:',
     $       n,m,t,se,ic
        call szcoll(n,m,t,se,ic)
        ym=157803.*ric*ric/t/(rm*rm)
        xn=(1./(rn*rn)-1./(rm*rm))
        yn=157803.*ric*ric*xn/t
        sd=se*exp(min(50.,yn))*(rn*rn)/(rm*rm)
        if (lpri.gt.1)
     $   write (lun11,*)'after szcoll:',
     $       ym,xn,yn,sd,ric,rm,xn,rn,t
       endif
      else
      xn=(1./(rn*rn)-1./(rm*rm))
      f=-1.2456e-10*a/xn/xn
      yn=157803.*xn/t
      ym=157803./t/(rm*rm)
      z=1.94*xn*rn**0.43+yn
      if(n.eq.1) z=yn+0.45*xn
      dif=z-yn
c
      if (lpri.gt.1)
     $ write(lun11,*)'before expint:',yn,z
       call expint(yn,e1y)
       call expint(z,e1z)
      if (lpri.gt.1)
     $ write(lun11,*)'after expint:',yn,e1y,z,e1z
      e2=(1.-e1y)/yn-expo(-dif)*(1.-e1z)/z
      ann=-2.*f*rm*rm/xn/rn/rn
      bn=(4.-18.63/rn+36.24/rn/rn-28.09/(rn**3))/rn
      if(n.eq.1) bn=-0.603
      bnn=(1.+4./(xn*rn*rn*3)+bn/(rn**4*xn*xn))*4./(rm**3*xn*xn)
      if (lpri.gt.1)
     $ write(lun11,*)e2,rn,rm,ann,bn,bnn
c
      s=ann*((1./yn+0.5)*e1y/yn-(1./z+0.5)*e1z*expo(-dif)/z)
      s=s+e2*(bnn-ann*log(2./xn))
      s=1.095e-10*yn*yn*sqrt(t)*s/xn
      sm=s*rn*rn*expo(ym)
      if (lpri.gt.1)
     $  write(lun11,*)s,sm,yn,xn,t,dif
       sd=s*rn/rm*rn/rm
       se=s*expo(-yn)
      if (lpri.gt.1)
     $  write(lun11,*)'erc=',se,sd
c      sg=13.60*sq*xn
c      sl=13.60*s*xn
        endif
      if (lpri.gt.1)
     $  write(lun11,*)'erc return:',se,sd
c
      return
      end
      function exp10(x)
c
      implicit none
      real*8 exp10,x
c
      exp10=exp(2.30259*x)
c
      return
      end
      subroutine expint(x,em1)
c
c expint is a subroutine to calculate the value of e1, the exponential
c integral or em1=x*expo(x)*e1 at the point x.  the polynomial
c expressions that are used come from abromowitz and stegen
c     author:  T. Kallman
c
      implicit none
      real*8 x,em1,b1,b2,b3,b4,c1,c2,c3,c4,a0,a1,a2,a3,a4,a5,e1,expo
c
      if(x.le.1.) go to 100
c
      b1=9.5733223454
      b2=25.6329561486
      b3=21.0996530827
      b4=3.9584969228
      c1=8.5733287401
      c2=18.0590169730
      c3=8.6347608925
      c4=0.2677737343
      em1=x**4+c1*x**3+c2*x*x+c3*x+c4
      em1=em1/(x**4+b1*x*x*x+b2*x*x+b3*x+b4)
c      e1=em1/x/expo(x)
      go to 200
c
 100   continue
      a0=-0.57721566
      a1=0.99999193
      a2=-0.24991055
      a3=0.05519968
      a4=-0.00976004
      a5=0.00107857
      if (x.gt.0)then
      e1= a0+a1*x+a2*x*x+a3*x**3+a4*x**4+a5*x**5-log(x)
      else
      e1=-a0+a1*x+a2*x*x+a3*x**3+a4*x**4+a5*x**5-log(-x)
      endif
      em1=e1*x*expo(x)
c
 200   continue
c
      return
      end
      real*8 function expo(x)
c
      implicit none
c
      real*8 x,crit
c
      crit=600.
c      expo=exp(x)
c      return
c      crit=60.
c      if (x.lt.-crit) then
c        expo=1.e-24
c      else
c        xtmp=min(x,crit)
        expo=exp(min(max(x,-crit),crit))
c      endif
c
      return
      end
      subroutine fact(n,x)
c
c to calculate the factorial of an integer n.  the output x is
c the natural log of n factorial.
c
      implicit none
      real*8 x
      integer i,n
c
      x=0.
      if(n.ne.0) then
      do  i=1,n
        x=x+log(float(i))
        enddo
      endif
c
      return
      end
      real*8 function fbg(u,gam)
c
c     this function computes the free-free gaunt factor
c      u=h nu/kt
c      gam=z**2 ry/kt
c         z=charge of scattering ion
c         ry=rydberg constant
c         kt=kt, etc.
c
      implicit none
      save
c
      real*8 a,a1,a2,a3,ai,ak,born,g1,g2,gam,
     &     gam1,gam2,gam3,p,power,t,u,u1,u2,expo
      real*8 u4
      integer m,m1,n
c      real*8  t,ai,ak,u4
      dimension a(6,7,3),gam2(6),gam3(6)
      dimension a1(6,7),a2(6,7),a3(6,7)
c
      equivalence (a1(1,1),a(1,1,1)),(a2(1,1),a(1,1,2)),
     &             (a3(1,1),a(1,1,3))
c
      data gam2/.7783,1.2217,2.6234,4.3766,20.,70./
      data gam3/1.,1.7783,3.,5.6234,10.,30./
      data a1/1.001,1.004,1.017,1.036,1.056,1.121,1.001,
     &     1.005,1.017,1.046,1.073,1.115,.9991,1.005,
     &     1.030,1.055,1.102,1.176,.9970,1.005,1.035,
     &     1.069,1.134,1.186,.9962,1.004,1.042,1.100,
     &     1.193,1.306,.9874,.9962,1.047,1.156,1.327,
     &     1.485,.9681,.9755,.8363,1.208,1.525,1.955/
      data a2/.30290,.16160,.04757,.01300,.00490,-.00320,
     &     .49050,.21550,.08357,.02041,.00739,.00029,
     &     .65400,.28330,.08057,.03257,.00759,-.00151,
     &     1.0290,.39100,.12660,.05149,.01274,.00324,
     &     .95690,.48910,.17640,.05914,.01407,-.00024,
     &     1.2690,.75790,.32600,.10770,.02800,.00548,
     &     1.3270,1.0170,1.3980,.20500,.06050,.00187/
      data a3/ - 1.3230,-.25400,-.01571,-.001000,-.000184,
     &     .00008,-4.7620,-.33860,-.03571,-.001786,-.000300,
     &     .00001,-8.3490,-.42060,-.02571,-.003429,-.000234,
     &     .00005,-13.231,-.59000,-.04571,-.005714,-.000445,
     &     -.00004,-7.6720,-.68520,-.06430,-.005857,-.000420,
     &     .00004,-7.1430,-.99470,-.12000,-.010070,-.000851,
     &     -.00004,-3.1750,-1.1160,-.84140,-.018210,-.001729,
     &     .00023/
c
      gam1 = gam*1000.
      if ( gam1.gt.100. ) then
         power = -.134/(gam**.2097)
         fbg = 1.5*(3.*u)**power
         return
      else
         u2 = u**2
c
c*****compute born approximation gaunt factor
c
         u1 = u/2.
         t = u1/3.75
         u4 = u1/2.
         if ( u1.gt.2. ) then
c
            ak = 1.2533141 - .07832358/u4 + .02189568/u4**2 -
     &           .01062446/u4**3 + .00587872/u4**4 - .00251540/u4**5 +
     &           .00053208/u4**6
            ak = ak/(expo(u1)*sqrt(u1))
         else
            ai = 1.0 + 3.5156229*t**2 + 3.0899424*t**4 +
     &           1.2067492*t**6 + 0.2659732*t**8 + 0.0360768*t**10 +
     &           0.0045813*t**12
            ak = -1.*log(u4)*ai - .57721566 + .42278420*u4**2 +
     &           .23069758*u4**4 + .0348859*u4**6 + .00262698*u4**8 +
     &           .00010750*u4**10 + .0000074*u4**12
         endif
         born = .5513*expo(u1)*ak
c
c*****compute polymonial factor to multiply born approximation
c
         m=0
         n=0
         if ( gam1.ge.1. ) then
            if ( u.ge..003 ) then
               if ( u.le..03 ) n = 1
               if ( (u.le..3) .and. (u.gt..03) ) n = 2
               if ( (u.le.1.) .and. (u.gt..3) ) n = 3
               if ( (u.le.5.) .and. (u.gt.1.) ) n = 4
               if ( (u.le.15.) .and. (u.gt.5.) ) n = 5
               if ( u.gt.15. ) n = 6
               if ( gam1.le.1.7783 ) m = 1
               if ( (gam1.le.3.) .and. (gam1.gt.1.7783) ) m = 2
               if ( (gam1.le.5.6234) .and. (gam1.gt.3.) ) m = 3
               if ( (gam1.le.10.) .and. (gam1.gt.5.6234) ) m = 4
               if ( (gam1.le.30.) .and. (gam1.gt.10.) ) m = 5
               if ( (gam1.le.100.) .and. (gam1.gt.30.) ) m = 6
               m1 = m + 1
               g1 = (a(n,m,1)+a(n,m,2)*u+a(n,m,3)*u2)*born
               g2 = (a(n,m1,1)+a(n,m1,2)*u+a(n,m1,3)*u2)*born
               p = (gam1-gam3(m))/gam2(m)
               fbg = (1.0-p)*g1 + p*g2
               return
            endif
         endif
      endif
      fbg = born
c
      return
      end
      subroutine fbg2(u,gam, fbg)
c
c     this function computes the free-free gaunt factor
c      u=h nu/kt
c      gam=z**2 ry/kt
c         z=charge of scattering ion
c         ry=rydberg constant
c         kt=kt, etc.
c
      implicit none
      save
c
      real*8 fbg
      real*8 a,a1,a2,a3,ai,ak,born,g1,g2,gam,
     &     gam1,gam2,gam3,p,power,t,u,u1,u2,expo
      real*8 u4
      integer m,m1,n
c      real*8  t,ai,ak,u4
      dimension a(6,7,3),gam2(6),gam3(6)
      dimension a1(6,7),a2(6,7),a3(6,7)
c
      equivalence (a1(1,1),a(1,1,1)),(a2(1,1),a(1,1,2)),
     &             (a3(1,1),a(1,1,3))
c
      data gam2/.7783,1.2217,2.6234,4.3766,20.,70./
      data gam3/1.,1.7783,3.,5.6234,10.,30./
      data a1/1.001,1.004,1.017,1.036,1.056,1.121,1.001,
     &     1.005,1.017,1.046,1.073,1.115,.9991,1.005,
     &     1.030,1.055,1.102,1.176,.9970,1.005,1.035,
     &     1.069,1.134,1.186,.9962,1.004,1.042,1.100,
     &     1.193,1.306,.9874,.9962,1.047,1.156,1.327,
     &     1.485,.9681,.9755,.8363,1.208,1.525,1.955/
      data a2/.30290,.16160,.04757,.01300,.00490,-.00320,
     &     .49050,.21550,.08357,.02041,.00739,.00029,
     &     .65400,.28330,.08057,.03257,.00759,-.00151,
     &     1.0290,.39100,.12660,.05149,.01274,.00324,
     &     .95690,.48910,.17640,.05914,.01407,-.00024,
     &     1.2690,.75790,.32600,.10770,.02800,.00548,
     &     1.3270,1.0170,1.3980,.20500,.06050,.00187/
      data a3/ - 1.3230,-.25400,-.01571,-.001000,-.000184,
     &     .00008,-4.7620,-.33860,-.03571,-.001786,-.000300,
     &     .00001,-8.3490,-.42060,-.02571,-.003429,-.000234,
     &     .00005,-13.231,-.59000,-.04571,-.005714,-.000445,
     &     -.00004,-7.6720,-.68520,-.06430,-.005857,-.000420,
     &     .00004,-7.1430,-.99470,-.12000,-.010070,-.000851,
     &     -.00004,-3.1750,-1.1160,-.84140,-.018210,-.001729,
     &     .00023/
c
      gam1 = gam*1000.
      if ( gam1.gt.100. ) then
         power = -.134/(gam**.2097)
         fbg = 1.5*(3.*u)**power
         return
      else
         u2 = u**2
c
c*****compute born approximation gaunt factor
c
         u1 = u/2.
         t = u1/3.75
         u4 = u1/2.
         if ( u1.gt.2. ) then
c
            ak = 1.2533141 - .07832358/u4 + .02189568/u4**2 -
     &           .01062446/u4**3 + .00587872/u4**4 - .00251540/u4**5 +
     &           .00053208/u4**6
            ak = ak/(expo(u1)*sqrt(u1))
         else
            ai = 1.0 + 3.5156229*t**2 + 3.0899424*t**4 +
     &           1.2067492*t**6 + 0.2659732*t**8 + 0.0360768*t**10 +
     &           0.0045813*t**12
            ak = -1.*log(u4)*ai - .57721566 + .42278420*u4**2 +
     &           .23069758*u4**4 + .0348859*u4**6 + .00262698*u4**8 +
     &           .00010750*u4**10 + .0000074*u4**12
         endif
         born = .5513*expo(u1)*ak
c
c*****compute polymonial factor to multiply born approximation
c
         m=0
         n=0
         if ( gam1.ge.1. ) then
            if ( u.ge..003 ) then
               if ( u.le..03 ) n = 1
               if ( (u.le..3) .and. (u.gt..03) ) n = 2
               if ( (u.le.1.) .and. (u.gt..3) ) n = 3
               if ( (u.le.5.) .and. (u.gt.1.) ) n = 4
               if ( (u.le.15.) .and. (u.gt.5.) ) n = 5
               if ( u.gt.15. ) n = 6
               if ( gam1.le.1.7783 ) m = 1
               if ( (gam1.le.3.) .and. (gam1.gt.1.7783) ) m = 2
               if ( (gam1.le.5.6234) .and. (gam1.gt.3.) ) m = 3
               if ( (gam1.le.10.) .and. (gam1.gt.5.6234) ) m = 4
               if ( (gam1.le.30.) .and. (gam1.gt.10.) ) m = 5
               if ( (gam1.le.100.) .and. (gam1.gt.30.) ) m = 6
               m1 = m + 1
               g1 = (a(n,m,1)+a(n,m,2)*u+a(n,m,3)*u2)*born
               g2 = (a(n,m1,1)+a(n,m1,2)*u+a(n,m1,3)*u2)*born
               p = (gam1-gam3(m))/gam2(m)
               fbg = (1.0-p)*g1 + p*g2
               return
            endif
         endif
      endif
      fbg = born
c
      return
      end
      real*8 function ff2(x,lpri,lun11)
c
c     expression usind in collisional ionization rate coefficient following
c     arnaud and raymond (1992)
c     author:  T. Kallman
c
      implicit none
c
      real*8 x
      real*8 q(15),p(15)
      real*8 pp,qq,ptst,xprod
      integer j,lun11,lpri
c
      data q/1.,2.1958e+2,2.0984e+4,1.1517e+6,4.0349e+7,
     $      9.4900e+8,1.5345e+10,1.7182e+11,1.3249e+12,
     $      6.9071e+12,2.3531e+13,4.9432e+13,5.7760e+13,
     $      3.0225e+13,3.3641e+12/
      data p/1.,2.1658e+2,2.0336e+4,1.0911e+6,3.7114e+7,
     $       8.3963e+8,1.2889e+10,1.3449e+11,9.4002e+11,
     $       4.2571e+12,1.1743e+13,1.7549e+13,1.0806e+13,
     $       4.9776e+11,0./
c
      xprod=1.
      if (lpri.ne.0)
     $ write (lun11,*)'in ff2:',x
      pp=0.
      qq=0.
      do j=1,15
        ptst=1./xprod
        if ((ptst.lt.1.e+20).and.(xprod.lt.1.e+24/x)) then
          pp=pp+p(j)/xprod
          qq=qq+q(j)/xprod
          xprod=xprod*x
          if (lpri.ne.0)
     $     write (lun11,*)j,xprod,pp,qq,ptst
          endif
        enddo
      ff2=pp/(1.e-20+qq)/x/x
c
      return
      end
      subroutine fheader(unit,knam,atcredate,mdlname,status)

C
C   File Name:    fheader.f
C   Author:       W.T. Bridgman
C   Date:         January 1999
C   Abstract:     Routines for writing a stardard primary FITS
C                 header for XSTAR output files.
C
C     Create a FITS file with an empty primary header
C     nrhs columns and nrows row
C
C     Parameters:
C        unit    integer            File unit number
C        knam    char*16            File name to create
C        mdlname char*30            Model name for this run
C        status  integer            Returned status code
C
      implicit none
c
      character(16) knam, filename
      character(30) mdlname
c     the atomic data creation date
      character(63) atcredate
      integer unit, status
      integer lun11

      integer bitpix,naxis,naxes(2),group
      logical simple,extend
      integer blocksize

      status=0
c
      filename=knam

C     Delete the file if it already exists, so we can recreate it
      call deletefile(filename,status)
c      write (6,*)'in fheader after deletefile',filename,status

C     Get an unused Logical Unit Number to use to open the FITS file
c      call ftgiou(unit,status)
      call getlun(unit)
c      write (6,*)'in fheader after ftgiou',unit,status

C     open the FITS file, with write access
      lun11=6
c      write (lun11,*)'opening fits unit ',unit, filename
      blocksize=1
      call ftinit(unit,filename,blocksize,status)
      if (status .gt. 0)call printerror(lun11,status)
c
c      if (status .gt. 0)stop


C     initialize parameters for primary array
      simple=.true.
      bitpix=16
      naxis=0
      naxes(1)=0
      naxes(2)=0
      extend=.true.
      group=1
C     write the required primary header keywords
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,group,extend,status)
      if (status .gt. 0)call printerror(lun11,status)

C     now add additional keywords
      call ftpcom(unit,'***********************************',status)
      call ftpkys(unit,'CREATOR','XSTAR version 2.2.1bn26',
     $ 'Program which generated this file',status)
      if (status .gt. 0)call printerror(lun11,status)

C     Extract the system date
      call ftpdat(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

C     Save run-specific information
      call ftpkys(unit,'MODEL',mdlname,'Model name for this run',status)
      if (status .gt. 0)call printerror(lun11,status)

      call ftpkys(unit,'ATDATA',atcredate,
     $  'Atomic data creation date',status)
      if (status .gt. 0)call printerror(lun11,status)


      return
      end
      subroutine find53(stmpp,etmpp,ntmp,efnd,sg,jlo,lun11,lpri)
c
c     this routine finds the continuum bin index.
c     author T. Kallman
c
      implicit none
c
      integer ntmp
      real*8 stmpp(ntmp),etmpp(ntmp)
      integer jlo,lun11,lpri
      real*8 efnd,sg
      integer ml2,mlp
      real*8 del1,del2,alg1,alg2,algtmp
c
c      lpri=0
c
c      if ((efnd.ge.etmpp(1)).and.(efnd.le.etmpp(ntmp))) then
      if (lpri.ne.0) write (lun11,*)'in find53:',efnd,ntmp,
     $    etmpp(1),etmpp(ntmp),stmpp(1),stmpp(ntmp)
      if ((efnd.ge.0.).and.(efnd.le.etmpp(ntmp))) then
        call hunt3(etmpp,ntmp,efnd,jlo,0,lun11)
        ml2=max(jlo,1)
        ml2=min(ml2,ntmp-1)
        mlp=ml2+1
        if (mlp.eq.ntmp) then
            alg1=log(max(stmpp(mlp),1.e-26)/max(stmpp(ml2),1.e-26))
            alg2=log(max(etmpp(mlp),1.e-26)/max(etmpp(ml2),1.e-26))
            algtmp=alg1/alg2
            sg=stmpp(ml2)*(efnd/etmpp(ml2))**algtmp
          else
            del1=(efnd-etmpp(ml2))/(etmpp(mlp)-etmpp(ml2))
            del2=(efnd-etmpp(mlp))/(etmpp(mlp)-etmpp(ml2))
            sg=-stmpp(ml2)*del2+stmpp(mlp)*del1
          endif
         sg=max(0.,sg)
         if (lpri.ne.0)
     $     write (lun11,*)sg,ml2,stmpp(ml2),stmpp(mlp),
     $           del1,del2,efnd,etmpp(ml2)
         else
              sg=0.
         endif
c
      return
      end
      subroutine fitsclose(lun11,unit,status)
C
C     Close the file & release the unit number
c     author: T. Bridgman
C
C     Parameters:
C        unit    integer            File unit number
C        status  integer            Returned status code
c
      implicit none
      integer unit, status,lun11
c
c      write (6,*)'closing fits unit ',unit
      call ftclos(unit, status)
c      call ftfiou(unit, status)
      close(unit)
      if (status .gt. 0)call printerror(lun11,status)
c
      return
      end
      function flinabs(ptmp)
c
c     this function will treat the absorption of incident
c     continuum in the line
c
      implicit none
c
      real*8 flinabs, ptmp
c
c      flinabs=0.
      flinabs=1.
c
      return
      end
      subroutine fnappend(knam,nint)
      character*16 knam
      if (nint.gt.9) then
          write (knam(3:4),'(i2)')nint
        else
          write (knam(3:3),'(a1)')'0'
          write (knam(4:4),'(i1)')nint
        endif
      return
      end
      subroutine fparmlist(unit,hdunum,mdlname,npar,parname,partype,
     $                    parval,parcomm,nloopctl,status,lun11)
c
c     author: T. Bridgman
c     parameters:
c        unit    integer            file unit number
c        hdunum  integer            number of last hdu written
c        mdlname char*30            model name for this run
c        npar    integer            number of parameters passed
c        parname char*20(999)       parameter name
c        partype char*10(999)       parameter type
c        parval  real(999)          parameter values converted to reals
c        parcomm char*30(999)       parameter comments & string values
c        nloopctl integer           loop control parameter
c        status  integer            returned status code
c
      implicit none
c     passed parameters
      character(30) mdlname
      integer unit, status, hdunum, npar, nloopctl
      character(20) parname(55)
      character(10) partype(55)
      real*8 parval(55)
      real parval4(55)
      character(30) parcomm(55)
c     parameter info
      integer idat1(6000000)   !jg
      integer lun11
      integer tfields,nrows,varidat
      character(16) ttype(5),tform(5),tunit(5)
      integer colnum,frow,felem,hdutype
      integer ll,mm
      character(30) extname
      character(4) ktmp2
c
      data tform/'1I','20A','1E','10A','30A'/
      data ttype/'index','parameter','value','type','comment'/
      data tunit/' ',' ',' ',' ',' '/
c
      nrows=npar
      varidat=0
c
      do mm=1,55
        parval4(mm)=parval(mm)
        enddo
c
c     move to the last hdu (hdunum) in the file
      call ftmahd(unit,hdunum,hdutype,status)
      if (status .gt. 0)call printerror(lun11,status)
c
c     append a new empty extension after the last hdu
      call ftcrhd(unit,status)
      if (status .gt. 0)call printerror(lun11,status)
c
c     define parameters for the binary table (see the above data statements)
      tfields=5
c
c     build extension name
      extname='PARAMETERS'
      if(nloopctl.gt.0) then
          write(ktmp2,'(i4.4)')nloopctl
c          extname='parameters_' // ktmp2
          endif
c
c     write the required header parameters for the binary table
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,
     $              varidat,status)
      if (status .gt. 0)call printerror(lun11,status)
c
c     save run-specific information
      call ftpkys(unit,'MODEL',mdlname,'model name for this run',status)
      if (status .gt. 0)call printerror(lun11,status)
c
c     set 'global' parameters for writing fits columns
      frow=1
      felem=1
c
c     column  1  (index)
      colnum=1
      do ll=1,nrows
         idat1(ll)=ll
         enddo
      call ftpclj(unit,colnum,frow,felem,nrows,idat1,status)
      if (status .gt. 0)call printerror(lun11,status)
c
c     column  2  (parameter name)
      colnum=2
      call ftpcls(unit,colnum,frow,felem,nrows,parname,status)
      if (status .gt. 0)call printerror(lun11,status)

c     column  3  (parameter value)
      colnum=3
      call ftpcle(unit,colnum,frow,felem,nrows,parval4,status)
      if (status .gt. 0)call printerror(lun11,status)

c     column  4 (parameter type)
      colnum=4
      call ftpcls(unit,colnum,frow,felem,nrows,partype,status)
      if (status .gt. 0)call printerror(lun11,status)

c     column  5 (parameter comment)
      colnum=5
      call ftpcls(unit,colnum,frow,felem,nrows,parcomm,status)
      if (status .gt. 0)call printerror(lun11,status)

c----------------------------------------------------------------
c     compute checksums
      call ftpcks(unit,status)
      if (status .gt. 0)call printerror(lun11,status)
c
      return
      end
      subroutine freef(lpri,lun11,epi,ncn2,t,xpx,xee,opakc)
c
c     this sub-routine computes the free-free opacity and
c     include it into the total one (opakc)
c     author:  J. Garcia (July 2008)
c
      implicit none
c
      include './PARAM'
c
      real*8 opakc(ncn),epi(ncn)
      real*8 t, opaff, ekt, t6, temp
      real*8 xpx, xee, xnx, enz2, cc
      real*8 gam, gau, fbg, zz
      integer numcon,lpri,lun11,ncn2,kk
c
c
      data cc/2.614e-37/
c
      if (lpri.gt.0) write (lun11,*)'in freef',t
c
      numcon=ncn2
      xnx=xpx*xee
      ekt = t*(0.861707)
      t6 = t/100.
      enz2=(1.4)*xnx
      zz=1.
      do kk=1,numcon
         temp = epi(kk)/ekt
         gam = zz*zz*(0.158)/t6
         gau = 1.
!         if ( temp.lt.100. ) gau = fbg(temp,gam)
!         if ( temp.lt.100. )
!     1    gau = 10.**(0.2258*epi(kk)**(0.08)*(4.094-log10(epi(kk)))    !JG
!     2      +log10(t6)*(0.133*(4.094-log10(epi(kk)))-0.2)-0.538)   !JG

         opaff = cc*xnx*enz2*gau/sqrt(t)/epi(kk)**3.
     1           *(1. - exp(-temp))

         opakc(kk) = opakc(kk) + opaff
c
!!! THIS IS A TEST
!         if(opakc(kk).gt.6.65e-7)opakc(kk)=6.65e-7
c
      enddo
c
      return
      end
      subroutine fstepr(unit,hdunum,radin,radout,rdel,t,pres,
     $                xcol,xee,xpx,xi,
     $                idat1,rdat1,kdat1,nptrs,npnxt,npfi,
     $                npfirst,npar,npilev,
     $                xilev,bilev,rniss,nloopctl,
     $                lun11,status)
C
C
C     Write Data for each radial zone to an individual extension
c     author: T. Bridgman
C
C     Append a FITS extension binary table containing
C     nrhs columns and at most nrhdimj rows
C
C     Parameters:
C        unit    integer            File unit number
C        hdunum  integer            Number of last HDU written
C        radin   real*8               inner radius of shell
C        radout  real*8               outer radius of shell
C                                   nb but now it is delr in the call
C        temp    real*8               temperature of shell
C        pres    real*8               pressure in shell
C        nrhdimj  integer            Maximum number of rows
C        idat1   integer(nidat1)    Needed by the atomic database
C        rdat1   real(nidat1)       Needed by the atomic database
C        kdat1   char*nidat1        Needed by the atomic database
C        nptrs                      Needed by the atomic database
C        npnxt                      Needed by the atomic database
C        npfi                       Needed by the atomic database
C        npfirst                    Needed by the atomic database
C        npcon                      Needed by the atomic database
C        npconi                     Needed by the atomic database
C        npcon2                     Needed by the atomic database
C        xilev   real(nrhdimj)       Fractional level population array
C        cemab   real(2,nrhdimj)     Recombination emission
C        opakab  real(nrhdimj)       Opacity
C        tauc    real(2,nrhdimj)     Optical depth
C        poptol  real*8               Tolerance for population level
C        nloopctl integer           Loop control variable
C        nzone   integer            Pass number through iteration process
C        status  integer            Returned status code
C
      implicit none
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=200000)
c
C     Allocation for passed parameters
      real*8 xilev(nnml),bilev(nnml),rniss(nnml)
      real*8 rdat1(nrdat1)
      real*8 radin, radout,rdel, t, pres, xcol,xee,xpx,xi
      real rtmp
      integer unit,hdunum, nrows, status, nloopctl
      integer idat1(nidat1),nptrs(nptt,ndat2)

c     pointers to master data
      integer npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni),npar(ndat2)
      integer npilev(nd,nni)
C     Internal work areas
      real rwrk1(nnml),rwrk2(nnml), elev(nnml)
      integer ntptr(nnml)
      integer natomic(nnml), mllev(nnml),nupper(nnml)
      character(10) kion(nnml)
      character(20) klevt(nnml)
      integer tfields,varidat
      character(16) ttype(9),tform(9),tunit(9)
      integer colnum,frow,felem,hdutype, klel, mlel, jk, ltyp
      integer lrtyp, lcon, nrdt, nidt, mmlv, mm, lun11, lpril,lpri
      integer mllel, klion, mlion, jkk, kl
      integer mt2, mlleltp, nnz, nions
      character(43) extname
      character(30) ktmp2
      character(1) kdat1(nkdat1)
C     Database manipulation quantities
      real*8  xeltp
      integer  nkdt
      character(1) klev(100,nd)
      real*8 rnissl(nnml)
      integer j,nkdti,np1ki
      real*8 rlev(10,nd)
      integer ilev(10,nd),nlpt(nd),iltp(nd),nlev
      integer mm2,mmtmp,kkkl,lk,mlm
      integer np1i,np1r,np1k
      real*8 eth
      character(10) kdtmp

      data tform/'1J','1I','1E','8A','1I','20A','1E','1E','1I'/
      data ttype/'index','ion_index','e_excitation','ion',
     $ 'atomic_number','ion_level','population','lte',
     $  'upper index'/
      data tunit/' ',' ','eV',' ',' ',' ',' ',' ',' '/
c
      varidat=0
      lpril=0
      lpri=lpril
c
      status=0
c
C     Move to the last HDU (hdunum) in the file
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Moving to end-of-FITS file'
      call ftmahd(unit,hdunum,hdutype,status)
      if (status .gt. 0)call printerror(lun11,status)
c
c
C     append a new empty extension after the last HDU
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Create the new extension'
      call ftcrhd(unit,status)
      if (status .gt. 0)call printerror(lun11,status)
C
C     Extracting data from the Atomic Database here
c
C
C     lpril is flag for printing debug information
       nions=0
      if (lpril.ne.0) then
        write (lun11,*)'raw data'
        do j=1,nnml
          if (xilev(j).gt.1.e-37)
     $     write (lun11,*)j,xilev(j)
          enddo
        endif
c
C     initialize line counter
      mmlv=0
C     First look for element data (jk is element index)
      klel=11
      mlel=npfirst(klel)
      jkk=0
      jk=0
      do while (mlel.ne.0)
c
c       get element data
       jk=jk+1
        mt2=mlel-1
        call drd(ltyp,lrtyp,lcon,
     $     nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,
     $     nptrs,0,lun11)
        mllel=idat1(np1i+nidt-1)
        xeltp=rdat1(np1r)
        nnz=idat1(np1i)
        if (lpril.ne.0)
     $    write (lun11,*)'element:',jk,mlel,mllel,nnz,
     $                  (kdat1(np1k-1+mm),mm=1,nkdt),xeltp
C       ignore if the abundance is small
        if (xeltp.lt.1.e-10) then
            jkk=jkk+nnz
          else
c
c           now step thru ions (jkk is ion index)
            klion=12
            mlion=npfirst(klion)
            jkk=0
            kl=0
            do while ((mlion.ne.0).and.(kl.lt.nnz))
c
              jkk=jkk+1
C             retrieve ion name from kdati
              mlm=mlion-1
              call drd(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidt,np1i,nkdti,np1ki,mlm,
     $            nptrs,0,lun11)

C             if not accessing the same element, skip to the next element
              mlleltp=idat1(np1i+nidt-2)
              if (mlleltp.eq.mllel) then
c
                kl=kl+1
                if (lpril.ne.0)
     $              write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
     $                          (kdat1(np1ki-1+mm),mm=1,nkdti)
c
c               get level data
                call func2l(jkk,lpril,lun11,t,xee,xpx,
     $              idat1,rdat1,kdat1,nptrs,
     $              npar,npnxt,npfi,
     $              rnissl,rlev,ilev,
     $              nlpt,iltp,nlev,klev)
c
c               step thru levels
                do mm2=1,nlev
c
c                 get level pointer
                  mmtmp=npilev(mm2,jkk)
                  if (mmtmp.ne.0) then
                    kkkl=mmtmp
                    mmlv=mmtmp
c
c                   test for level pop
                    if (xilev(kkkl).gt.1.d-34) then
c
c                     get data
                      eth=rlev(1,mm2)
                      nions=nions+1
                      mllev(nions)=idat1(np1i+nidt-2)
C                     Note that rwrk1 must be written to the file before
C                     it is overwritten in subsequent columns
                      rwrk1(nions)=xilev(mmlv)
                      rwrk2(nions)=rniss(mmlv)
                      elev(nions)=eth
                      ntptr(nions)=kkkl
                      natomic(nions)=nnz
                      nupper(nions)=mm2
                      do mm=1,nkdti
                        write (kdtmp(mm:mm),'(a1)')kdat1(np1ki-1+mm)
                        enddo
                      do mm=nkdti+1,9
                        write (kdtmp(mm:mm),'(a1)')' '
                        enddo
                      kion(nions)=kdtmp
                      write(klevt(nions),'(20a1)')
     $                        (klev(mm,mm2),mm=1,20)
                      if (lpri.ne.0) then
                        write (lun11,*)nions,xilev(mmlv),
     $                         rdat1(np1r),nnz,mmlv,kkkl
                        write (lun11,9296)kkkl,
     $                      (kdat1(np1i-1+mm),mm=1,20),
     $                      (klev(lk,mm2),lk=1,20),eth,xilev(kkkl),
     $                      rniss(kkkl),bilev(kkkl)
 9296                   format (1x,i6,1x,(40a1),7(1pe13.5))
                        endif
c
c                     end of test for level pop
                      endif
c
c                   end of test for level pointer
                    endif
c
c                 end of step thru levels
                  enddo
c
c               end of test for element
                endif
c
C             Go to next ion
              mlion=npnxt(mlion)
              enddo
c
C           end of test for abundance
            endif
c
        mlel=npnxt(mlel)
C       Go to next element
        enddo

c

C     End of atomic database extraction
C----------------------------------------------------------------
C     define parameters for the binary table (see the above data statements)
      nrows=nions
      tfields=9
C     Build extension name
      extname='XSTAR_RADIAL'
      if(nloopctl.gt.0) then
          write(ktmp2,'(I4.4)')nloopctl
          extname='XSTAR_RADIAL_' // ktmp2
          endif

      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Write table headers'
C     write the required header parameters for the binary table
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,
     $              varidat,status)
      if (status .gt. 0)call printerror(lun11,status)

      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Add some more keywords'

C     Write some model parameters in the extension header
      call ftpcom(unit,'***********************************',status)
      if (status .gt. 0)call printerror(lun11,status)

      call ftpcom(unit,'Model Keywords',status)
      if (status .gt. 0)call printerror(lun11,status)

C     Write values to 3 decimal places
      rtmp=radin
      call ftpkye(unit,'RINNER',rtmp,3,'[cm] Inner shell radius',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=radout
      call ftpkye(unit,'ROUTER',rtmp,3,'[cm] Outer shell radius',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=rdel
      call ftpkye(unit,'RDEL',rtmp,3,'[cm] distance from face',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=t
      call ftpkye(unit,'TEMPERAT',rtmp,3,'[10**4K] Shell Temperature',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=pres
      call ftpkye(unit,'PRESSURE',rtmp,3,'[dynes/cm**2] Shell Pressure',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xcol
      call ftpkye(unit,'COLUMN',rtmp,3,'[/cm**2] Column ',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=xee
      call ftpkye(unit,'XEE',rtmp,3,'electron fraction',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xpx
      call ftpkye(unit,'DENSITY',rtmp,3,'[/cm**3] Density',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xi
      call ftpkye(unit,'LOGXI',rtmp,3,
     $ '[erg cm/s] log(ionization parameter)',status)
      if (status .gt. 0)call printerror(lun11,status)

C-------------------------------------------------------------------
C     Step through the columns and write them to the file
C
C     set 'global' parameters for writing FITS columns
      frow=1
      felem=1

C     column  1  (Line number)
      colnum=1
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Writing Column ',colnum
      call ftpclj(unit,colnum,frow,felem,nions,ntptr,status)
      if (status .gt. 0)call printerror(lun11,status)

C     column  2 (Level number of this ion)
      colnum=2
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Writing Column ',colnum
      call ftpclj(unit,colnum,frow,felem,nions,mllev,status)
      if (status .gt. 0)call printerror(lun11,status)

C     column  3  (Energy)
      colnum=3
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nions,elev,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  4  (Ion)
      colnum=4
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Writing Column ',colnum
      call ftpcls(unit,colnum,frow,felem,nions,kion,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  5  (Atomic Number)
      colnum=5
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Writing Column ',colnum
      call ftpclj(unit,colnum,frow,felem,nions,natomic,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  6 (Level Designation)
      colnum=6
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Writing Column ',colnum
      call ftpcls(unit,colnum,frow,felem,nions,klevt,status)
      if (status .gt. 0)call printerror(lun11,status)

C----------------------------------------------------------------
C     column 7 (Level population)
C     rwrk1 can be safely overwritten after this step

      colnum=7
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nions,rwrk1,status)
      if (status .gt. 0)call printerror(lun11,status)


      colnum=8
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nions,rwrk2,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  9 (upper level index)
      colnum=9
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr: Writing Column ',colnum
      call ftpclj(unit,colnum,frow,felem,nions,nupper,status)
      if (status .gt. 0)call printerror(lun11,status)

c----------------------------------------------------------------
C     Compute checksums
      call ftpcks(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

      return
      end
      subroutine fstepr2(unit,hdunum,radin,radout,rdel,temp,pres,
     $                xcol,xee,xpx,xi,
     $                idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $                nplin,nlsvn,rcem,oplin,tau0,nloopctl,
     $                lun11,status)
C
C     Append a FITS extension binary table containing
C     nrhs columns and at most nrhdimj rows
c     author: T. Bridgman
C
C     Parameters:
C        unit    integer            File unit number
C        hdunum  integer            Number of last HDU written
C        radin   real*8               inner radius of shell
C        radout  real*8               outer radius of shell
C                                   nb but now it is delr in the call
C        temp    real*8               temperature of shell
C        pres    real*8               pressure in shell
C        nrhdimj  integer            Maximum number of rows
C        idat1   integer(nidat1)    Needed by the atomic database
C        rdat1   real(nidat1)       Needed by the atomic database
C        kdat1   char*nidat1        Needed by the atomic database
C        nptrs                      Needed by the atomic database
C        npnxt                      Needed by the atomic database
C        npfi                       Needed by the atomic database
C        npfirst                    Needed by the atomic database
C        npcon                      Needed by the atomic database
C        npconi                     Needed by the atomic database
C        npcon2                     Needed by the atomic database
C        xilev   real(nrhdimj)       Fractional level population array
C        cemab   real(2,nrhdimj)     Recombination emission
C        opakab  real(nrhdimj)       Opacity
C        tauc    real(2,nrhdimj)     Optical depth
C        poptol  real*8               Tolerance for population level
C        nloopctl integer           Loop control variable
C        nzone   integer            Pass number through iteration process
C        status  integer            Returned status code
C
      implicit none
      include './PARAM'
      integer mllz

      integer nptmpdim
      parameter (nptmpdim=400000)
c
C     Allocation for passed parameters
      real*8 tau0(2,nnnl), rcem(2,nnnl)
      real*8 rdat1(nrdat1)
      real rtmp
      real*8 radin, radout,rdel, temp, pres,xcol,xee,xpx,xi
      integer unit,hdunum, nrows, status, nloopctl
      integer idat1(nidat1),nptrs(nptt,ndat2)
c     line opacities
      real*8 oplin(nnnl)

c     pointers to master data
      integer npnxt(ndat2)
      integer npfi(ntyp,nni),npar(ndat2)
      integer nplin(nnnl)

C     Internal work areas
      real rwrk1(nptmpdim)
      integer ntptr(nptmpdim)
      character(10) kion(nptmpdim)
      character(20) klevl(nptmpdim),klevu(nptmpdim),kblnk20
      integer tfields,varidat
      character(16) ttype(10),tform(10),tunit(10)
      integer colnum,frow,felem,hdutype,ll, ltyp
      integer lrtyp, lcon, nrdt, nidt, mm, lun11, lpril,lpri
      integer jkk, nlev, mlpar
      integer nlplmx,ln,lnn,ml,nilin2,nlpl,lmm,kltmpn,kltmpo,
     $         llo,lup,llofnd,lupfnd,nlines,nlsvn,
     $         k,kl2,kk,nilin,lm,mlm
      integer np1i,np1r,np1k
      real*8 eliml,elimh,elin,elmmtpp,elcomp
      real*8 rniss(nnml)
      real*8 aij,ergsev,etst,gglo,ggup,flin,ener
      integer idest1,idest2,ilevlo,ilevup,j,kl,ktt,lk
      real*8  rlev(10,nd)
      integer ilv(10,nd),nlpt(nd),iltp(nd)
      character(1) klev(100,nd)
      character(33) extname
      character(20) ktmp2,klablo,klabup
      character(1) kdat1(nkdat1)
      character(9) kinam1

C     Database manipulation quantities
      integer nkdt
      character(1) kblnk,kdtmp(200)
      integer kltmp(nptmpdim)
      real elsv(nptmpdim)
      logical done

      data kblnk/' '/
      data kblnk20/'                    '/
c
      data tform/'1J','1E','8A','20A','20A','1E','1E','1E',
     $ '1E','1E'/

      data ttype/'index','wavelength','ion',
     $ 'lower_level','upper_level','emis_inward',
     $ 'emis_outward','opacity','tau_in','tau_out'/

      data tunit/' ','A',' ',' ',' ','erg/cm^3/s',
     $ 'erg/cm^3/s','/cm',' ',' '/

      varidat=0
c
      lpri=0
      lpril=lpri
c

C     Move to the last HDU (hdunum) in the file
c     if (lpri.ne.0)
c    $ write(lun11,*)'fstepr2: Moving to end-of-FITS file'
      call ftmahd(unit,hdunum,hdutype,status)
      if (status .gt. 0)call printerror(lun11,status)

C     append a new empty extension after the last HDU
c     if (lpri.ne.0)
c    $ write (lun11,*)'fstepr2: Create the new extension'
      call ftcrhd(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

C----------------------------------------------------------------
C
C     Extracting data from the Atomic Database here
C
c     if (lpri.ne.0)
c    $ write (lun11,*)' '
c
c     print important lines
c     if (lpri.ne.0)
c    $ write (lun11,*)'emission line luminosities (erg/sec/10**38))',
c    $                 nlsvn
         nlplmx=200000
c
c     step through lines
      nlpl=0
      do lnn=1,nlsvn
c
        if ((rcem(1,lnn).gt.1.d-64).or.
     $      (rcem(2,lnn).gt.1.d-64).or.
     $      (oplin(lnn).gt.1.d-64)) then
c
c         get line data
          ln=lnn
          ml=nplin(ln)
          mlm=ml-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          if (lpri.ne.0) write (lun11,*)ln,ml,rdat1(np1r)
c
c         exclude rate type 14
          elin=abs(rdat1(np1r))
          if ((lrtyp.ne.14).and.(abs(elin).gt.0.1)
     $       .and.(abs(elin).lt.9.e+9)) then
c
            ergsev=1.602197e-12
            ener=ergsev*(12398.41)/max(elin,1.e-24)
            etst=ener/ergsev
            idest1=idat1(np1i)
            idest2=idat1(np1i+1)
            aij=rdat1(np1r+2)
            if (lpri.ne.0) write (lun11,*)'line data',elin,ener,etst,
     $                       idest1,idest2,aij
c
c           get ion data
            nilin=npar(ml)
            mlm=nilin-1
            call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $        nptrs,0,lun11)
            do ktt=1,min(8,nkdt)
              write (kinam1(ktt:ktt),'(a1)')kdat1(np1k-1+ktt)
              enddo
            do ktt=nkdt+1,9
              write (kinam1(ktt:ktt),'(a1)')kblnk
              enddo
c
c           now find level data
            jkk=idat1(np1i+nidt-1)
            if (lpri.ne.0) write (lun11,*)'ion',kinam1,jkk
            call func2l(jkk,lpri,lun11,temp,xee,xpx,
     $              idat1,rdat1,kdat1,nptrs,
     $              npar,npnxt,npfi,
     $              rniss,rlev,ilv,
     $              nlpt,iltp,nlev,klev)

            ggup=rlev(2,idest1)
            gglo=rlev(2,idest2)
            do lk=1,20
              write (klablo(lk:lk),'(a1)')klev(lk,idest1)
              write (klabup(lk:lk),'(a1)')klev(lk,idest2)
              enddo
            flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo)
            ilevlo=idest1
            ilevup=idest2
c
            nlpl=nlpl+1
            nlpl=min(nlpl,nlplmx)
            ntptr(nlpl)=lnn
            elsv(nlpl)=elin
            kion(nlpl)=kinam1
            klevl(nlpl)=klablo
            klevu(nlpl)=klabup
            j=ln
            if (lpri.ne.0)
     $        write (lun11,9929)j,elin,kinam1,
     $        (klev(mm,ilevlo),mm=1,20),(klev(mm,ilevup),mm=1,20),
     $        rlev(1,ilevlo),rlev(1,ilevup),rlev(2,ilevlo),
     $        rlev(2,ilevup),rlev(3,ilevlo),rlev(3,ilevup),
     $        ilv(1,ilevlo),ilv(1,ilevup),ilv(2,ilevlo),ilv(2,ilevup),
     $        ilv(3,ilevlo),ilv(3,ilevup)
 9929       format (1h ,i9,1pe13.5,1x,a9,1x,2(20a1,1x),6(1pe13.5),
     $          6i6)
            if (lpri.ne.0)
     $       write (lun11,*)j,elin,oplin(j),rcem(1,j),
     $                      rcem(2,j)
c
            endif
c
          endif
c
        enddo

c      if (nlpl.le.0) return
      nlpl=max(nlpl,1)
c

C     End of atomic database extraction
C----------------------------------------------------------------
C     define parameters for the binary table (see the above data statements)
      nrows=nlpl
      if (lpri.ne.0)
     $ write (lun11,*)'before header write'
      tfields=10
C     Build extension name
      extname='XSTAR_RADIAL'
      if(nloopctl.gt.0) then
          write(ktmp2,'(I4.4)')nloopctl
          extname='XSTAR_RADIAL_' // ktmp2
          endif

      if (lpri.ne.0)
     $ write (lun11,*)'fstepr2: Write table headers'
C     write the required header parameters for the binary table
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,
     $              varidat,status)
      if (status .gt. 0)call printerror(lun11,status)

      if (lpri.ne.0)
     $ write (lun11,*)'fstepr2: Add some more keywords'

C     Write some model parameters in the extension header
      call ftpcom(unit,'***********************************',status)
      if (status .gt. 0)call printerror(lun11,status)

      call ftpcom(unit,'Model Keywords',status)
      if (status .gt. 0)call printerror(lun11,status)

C     Write values to 3 decimal places
      rtmp=radin
      call ftpkye(unit,'RINNER',rtmp,3,'[cm] Inner shell radius',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=radout
      call ftpkye(unit,'ROUTER',rtmp,3,'[cm] Outer shell radius',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=rdel
      call ftpkye(unit,'RDEL',rtmp,3,'[cm] distance from face',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=temp
      call ftpkye(unit,'TEMPERAT',rtmp,3,'[10**4K] Shell Temperature',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=pres
      call ftpkye(unit,'PRESSURE',rtmp,3,'[dynes/cm**2] Shell Pressure',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xcol
      call ftpkye(unit,'COLUMN',rtmp,3,'[/cm**2] Column ',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=xee
      call ftpkye(unit,'XEE',rtmp,3,'electron fraction',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xpx
      call ftpkye(unit,'DENSITY',rtmp,3,'[/cm**3] Density',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xi
      call ftpkye(unit,'LOGXI',rtmp,3,
     $ '[erg cm/s] log(ionization parameter)',status)
      if (status .gt. 0)call printerror(lun11,status)

      if (lpri.ne.0)
     $ write (lun11,*)'after header write'
C-------------------------------------------------------------------
C     Step through the columns and write them to the file
C
C     set 'global' parameters for writing FITS columns
      frow=1
      felem=1


C     column  1  (Line number)
      colnum=1
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr2: Writing Column ',colnum,nlpl
      nlines=nlpl
      call ftpclj(unit,colnum,frow,felem,nlines,ntptr,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  2  (wavelength)
      colnum=2
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr2: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,elsv,status)
      if (status .gt. 0)call printerror(lun11,status)
      if (status .gt. 0) return


C     column  3  (Ion)
      colnum=3
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr2: Writing Column ',colnum
      call ftpcls(unit,colnum,frow,felem,nlines,kion,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  4 (lower Level Designation)
      colnum=4
      if (lpri.ne.0)
     $ write (lun11,*)'fstepr2: Writing Column ',colnum
      call ftpcls(unit,colnum,frow,felem,nlines,klevl,status)
      if (status .gt. 0)call printerror(lun11,status)

C     column  5 (Level Designation)
      colnum=5
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr2: Writing Column ',colnum
      call ftpcls(unit,colnum,frow,felem,nlines,klevu,status)
      if (status .gt. 0)call printerror(lun11,status)

C----------------------------------------------------------------

C     column  6
      colnum=6
      do ll=1,nlines
         rwrk1(ll)=0.
         if (ntptr(ll).ne.0)
     $      rwrk1(ll)=rcem(1,ntptr(ll))
         if (lpri.ne.0) write (lun11,*)ll,ntptr(ll),rwrk1(ll)
         enddo
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr2: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status)
      if (status .gt. 0)call printerror(lun11,status)

C     column  7
      colnum=7
      do ll=1,nlines
         rwrk1(ll)=0.
         if (ntptr(ll).ne.0)
     $    rwrk1(ll)=rcem(2,ntptr(ll))
         enddo
      if (lpri.ne.0)
     $ write (lun11,*)'fstepr2: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  8
      colnum=8
      do ll=1,nlines
         rwrk1(ll)=0.
         if (ntptr(ll).ne.0)
     $    rwrk1(ll)=oplin(ntptr(ll))
         enddo
      if (lpri.ne.0)
     $ write (lun11,*)'fstepr2: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  9
      colnum=9
      do ll=1,nlines
         rwrk1(ll)=0.
         if (ntptr(ll).ne.0)
     $    rwrk1(ll)=tau0(1,ntptr(ll))
         enddo
      if (lpri.ne.0)
     $ write (lun11,*)'fstepr2: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status)
      if (status .gt. 0)call printerror(lun11,status)

C     column  10
      colnum=10
      do ll=1,nlines
         rwrk1(ll)=0.
         if (ntptr(ll).ne.0)
     $    rwrk1(ll)=tau0(2,ntptr(ll))
         enddo
      if (lpri.ne.0)
     $ write (lun11,*)'fstepr2: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status)
      if (status .gt. 0)call printerror(lun11,status)


c----------------------------------------------------------------
C     Compute checksums
      call ftpcks(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

c
c
      return
      end
      subroutine fstepr3(unit,hdunum,radin,radout,rdel,t,prs,abel,
     $             xcol,xee,xpx,xi,
     $             idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $             npfirst,npilev,npconi2,ncsvn,
     $             rniss,cemab,cabab,opakab,tauc,nloopctl,
     $             lun11,status)
C
C     Append a FITS extension binary table containing
C     nrhs columns and at most nrhdimj rows
c     author: T. Bridgman
C
C     Parameters:
C        unit    integer            File unit number
C        hdunum  integer            Number of last HDU written
C        radin   real*8               inner radius of shell
C        radout  real*8               outer radius of shell
C                                   nb but now it is delr in the call
C        t    real*8               temperature of shell
C        prs    real*8               pressure in shell
C        nrhdimj  integer            Maximum number of rows
C        idat1   integer(nidat1)    Needed by the atomic database
C        rdat1   real(nidat1)       Needed by the atomic database
C        kdat1   char*nidat1        Needed by the atomic database
C        nptrs                      Needed by the atomic database
C        npnxt                      Needed by the atomic database
C        npfi                       Needed by the atomic database
C        npfirst                    Needed by the atomic database
C        npcon                      Needed by the atomic database
C        npconi                     Needed by the atomic database
C        npcon2                     Needed by the atomic database
C        xilev   real(nrhdimj)       Fractional level population array
C        cemab   real(2,nrhdimj)     Recombination emission
C        opakab  real(nrhdimj)       Opacity
C        tauc    real(2,nrhdimj)     Optical depth
C        poptol  real*8               Tolerance for population level
C        nloopctl integer           Loop control variable
C        nzone   integer            Pass number through iteration process
C        status  integer            Returned status code
C
      implicit none
      include './PARAM'
      integer mllz,mllz2

      integer nptmpdim
      parameter (nptmpdim=400000)
c
C     Allocation for passed parameters
      real*8 rdat1(nrdat1)
      real rtmp
      real*8 radin, radout,rdel, t, prs, xcol,xee,xpx,xi
      integer unit,hdunum, nrows, status, nloopctl
      integer idat1(nidat1),nptrs(nptt,ndat2)
c     line opacities
      real*8 oplin(nnnl)
      real*8 abel(nl)
      real*8 tauc(2,nnml)
      real*8 cemab(2,nnml),opakab(nnml),cabab(nnml)

c     pointers to master data
      integer npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni),npar(ndat2)
      integer npilev(nd,nni)
      integer nplin(nnnl)
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)

C     Internal work areas
      real*8 rniss(nnml)
      real*8 rlev(10,nd)
      integer ilev(10,nd),nlpt(nd),iltp(nd),nlev
      real rwrk1(nptmpdim),rwrk2(nptmpdim),rwrk3(nptmpdim),
     $  rwrk4(nptmpdim),rwrk5(nptmpdim),
     $  rwrk6(nptmpdim),rwrk7(nptmpdim),rwrk8(nptmpdim)
      integer ntptr(nptmpdim),ntptr2(nptmpdim)
      character(20) klevl(nptmpdim),klevu(nptmpdim),kblnk20
      character(1) klev(100,nd)
      integer tfields,varidat
      character(16) ttype(12),tform(12),tunit(12)
      integer colnum,frow,felem,hdutype,ll, ltyp
      integer lrtyp, lcon, nrdt, nidt, mm, lun11, lpril,lpri
      integer jkk,nidti
      integer nlplmx,ln,lnn,ml,nilin2,nlpl,lmm,kltmpn,kltmpo,
     $         llo,lup,llofnd,lupfnd,nlines,nlsvn,
     $         k,kl2,kk,nilin,lm,mlm,kkkl,idest1,idest2,lk,
     $         mlel,mlion,mllel,mlleltp,mlpar,mt2,mltype,nkdti,
     $         jk,kl,klel,klion,ncsvn,nlevmx,nnz,np1ki,mmlv
      integer np1i,np1r,np1k
      real*8 eliml,elimh,elin,elmmtpp,elcomp,eth,xeltp
      character(33) extname
      character(20) ktmp20
      character(8) ktmp8,kion(nptmpdim)
      character(20) ktmp2
      character(1) kdat1(nkdat1)
c     needed for upper level search
      integer jkk3,nlevp,ndtmp,iltmp,lcon2,lrtyp2,ltyp2,
     $         np1r2,nrdt2,np1i2,nidt2,np1k2,nkdt2
      real*8 ett

C     Database manipulation quantities
      integer nkdt
      character(1) kblnk,kdtmp(200)
      integer kltmp(50000),ilevlo(50000),ilevup(50000)
      real elsv(50000)
      logical done

      data kblnk/' '/
      data kblnk20/'                    '/
c
      data tform/'1J','1J','1E','8A','20A','20A','1E','1E','1E',
     $ '1E','1E','1E'/

      data ttype/'rrc index','level index','energy','ion',
     $ 'lower_level','upper_level','emis_inward',
     $ 'emis_outward','integrated absn','opacity',
     $  'tau_in','tau_out'/

      data tunit/' ',' ','ev',' ',' ',' ','erg/cm^3/s',
     $ 'erg/cm^3/s','erg/cm^3/s','/cm',' ',' '/

      varidat=0
c
      lpri=0
      lpril=lpri
c

C     Move to the last HDU (hdunum) in the file
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Moving to end-of-FITS file'
      call ftmahd(unit,hdunum,hdutype,status)
      if (status .gt. 0)call printerror(lun11,status)

C     append a new empty extension after the last HDU
      if (lpri.ne.0)
     $ write (lun11,*)'fstepr3: Create the new extension'
      call ftcrhd(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

C----------------------------------------------------------------
C
C     Extracting data from the Atomic Database here
C
      if (lpri.ne.0)
     $ write (lun11,*)'in fstepr3 '
c      lpril=0
      kltmpo=0
c
C     First look for element data (jk is element index)
      klel=11
      mlel=npfirst(klel)
      jk=0
      kk=0
      jkk=0
c
c     step through elements
      do while (mlel.ne.0)
c
c       get element data
        jk=jk+1
        mt2=mlel-1
        call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,
     $        nptrs,0,lun11)
        mllel=idat1(np1i+nidt-1)
        xeltp=rdat1(np1r)
        xeltp=abel(mllel)
        nnz=idat1(np1i)
        if (lpri.ge.1)
     $        write (lun11,*)'element:',jk,mlel,mllel,nnz,
     $                  (kdat1(np1k-1+mm),mm=1,nkdt)
c
C       ignore if the abundance is small
        if (xeltp.lt.1.e-10) then
            jkk=jkk+nnz
          else
c
c           now step thru ions (jkk is ion index)
            klion=12
            mlion=npfirst(klion)
            jkk=0
            kl=0
            do while ((mlion.ne.0).and.(kl.lt.nnz))
              jkk=jkk+1
c
C             retrieve ion name from kdati
              mlm=mlion-1
              call drd(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidti,np1i,nkdti,np1ki,mlm,
     $            nptrs,0,lun11)
c
C             if not accessing the same element, skip to the next element
              mlleltp=idat1(np1i+nidti-2)
              if (mlleltp.eq.mllel) then
c
                kl=kl+1
                if (lpri.ge.1)
     $            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
     $                        (kdat1(np1ki+mm-1),mm=1,nkdti)
c
c               now find level data
                call func2l(jkk,lpri,lun11,t,xee,xpx,
     $              idat1,rdat1,kdat1,nptrs,
     $              npar,npnxt,npfi,
     $              rniss,rlev,ilev,
     $              nlpt,iltp,nlev,klev)
c
c               now step through rate type 7 data
                mltype=7
                ml=npfi(mltype,jkk)
                mllz=0
                if (ml.ne.0) mllz=npar(ml)
                mlpar=0
                if (ml.ne.0) mlpar=npar(ml)
                do while ((ml.ne.0).and.(mlpar.eq.mllz))
c
c                 get rrc data
                  kkkl=npconi2(ml)
                  if (lpri.ne.0) write (lun11,*)kkkl,ml,idest1,
     $                    cemab(1,kkkl),cemab(2,kkkl)
c
c                 test for non-zero rrc data
                  if ((kkkl.gt.0).and.(kkkl.le.ndat2)
     $                .and.((cemab(1,kkkl).gt.1.e-36)
     $                .or.(cemab(2,kkkl).gt.1.e-36)
     $                .or.(cabab(kkkl).gt.1.e-36)
     $                .or.(opakab(kkkl).gt.1.e-36))) then
c
c
c                   increment buffer counter
                    kk=kk+1
                    if (kk.ge.nptmpdim) then
                      write (lun11,*)'buffer overflow in fstepr33'
                      kk=nptmpdim
                      endif
c
c                   get rrc  data
                    mlm=ml-1
                    call drd(ltyp,lrtyp,lcon,
     $                nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $                nptrs,0,lun11)
                    idest1=idat1(np1i+nidt-2)
                    nlevp=nlev
                    idest2=nlevp+idat1(np1i-1+nidt-3)-1
c
c                   label for lower level
                    do lk=1,20
                      write (ktmp20(lk:lk),'(a1)')klev(lk,idest1)
                      enddo
                    klevl(kk)=ktmp20
c
c                   label for upper level
                    write (ktmp20(1:20),'(a20)')'continuum           '
                    klevu(kk)=ktmp20
c
c                   ion label
                    do lk=1,nkdti
                      write (ktmp8(lk:lk),'(a1)')kdat1(np1ki+lk-1)
                      enddo
                    do lk=nkdti+1,8
                      write (ktmp8(lk:lk),'(a1)')kblnk
                      enddo
c
                    eth=rlev(4,idest1)-rlev(1,idest1)
                    ett=eth
c
c                   get upper level data
                    if (idest2.gt.nlevp) then
                      jkk3=jkk+1
                      if (lpri.gt.1)
     $                  write (lun11,*)jkk3,ndtmp,nlevp,idest2
                      ndtmp=npfi(13,jkk3)
                      if (lpri.gt.1)
     $                  write (lun11,*)jkk3,ndtmp,nlevp,idest2
                      if (ndtmp.le.0) stop 'ndtmp error'
                      mllz2=npar(ndtmp)
                      iltmp=0
                      do while ((ndtmp.ne.0).and.
     $                    (iltmp.ne.(idest2-nlevp+1)).and.
     $                    (npar(ndtmp).eq.mllz2))
                        mlm=ndtmp-1
                        call drd(ltyp2,lrtyp2,lcon2,
     $                    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $                    nptrs,0,lun11)
                        iltmp=idat1(np1i2+nidt2-2)
                        if (lpri.gt.1) write (lun11,*)nidt2,iltmp,ndtmp
                        ndtmp=npnxt(ndtmp)
                        if (ndtmp.le.0) stop 'ndtmp error'
                        enddo
c                     NB fix to excited level PI and rec
                      ett=ett+rdat1(np1r2)
                      eth=ett
                      if (lpri.gt.1)
     $                  write (lun11,*) ndtmp,iltmp,idest2,ett
c                     label for lower level
                      ktmp20=kblnk20
                      do lk=1,nkdt2
                        write (ktmp20(lk:lk),'(a1)')kdat1(np1k2+lk-1)
                        enddo
                      klevu(kk)=ktmp20
                      endif
c
c                   other data
                    kion(kk)=ktmp8
                    elsv(kk)=eth
                    ilevlo(kk)=idest1
                    ilevup(kk)=idest2
                    ntptr(kk)=kkkl
                    mmlv=npilev(idest1,jkk)
                    ntptr2(kk)=mmlv
                    if (lpri.ge.1)
     $                  write (lun11,981)kkkl,eth,idest1,
     $                    cemab(1,kkkl),cemab(2,kkkl)
 981                  format (1x,i6,1pe11.3,i6,6(1pe11.3))
c
c                   done with this rrc
                    endif
c
c                 end of loop over rrcs
                  ml=npnxt(ml)
                  if (ml.ne.0) mlpar=npar(ml)
                  enddo
c
c               end of test for element
                endif
c
C             Go to next ion
              mlion=npnxt(mlion)
              enddo
c
c         end of test for non-zero element abund
          endif
c
        mlel=npnxt(mlel)
C       Go to next element
        enddo
c
      nlpl=max(nlpl,1)
c

C     End of atomic database extraction
C----------------------------------------------------------------
C     define parameters for the binary table (see the above data statements)
      nrows=kk
      if (lpri.ne.0)
     $ write (lun11,*)'before header write'
      tfields=12
C     Build extension name
      extname='XSTAR_RADIAL'
      if(nloopctl.gt.0) then
          write(ktmp2,'(I4.4)')nloopctl
          extname='XSTAR_RADIAL_' // ktmp2
          endif

      if (lpri.ne.0)
     $ write (lun11,*)'fstepr3: Write table headers'
C     write the required header parameters for the binary table
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,
     $              varidat,status)
      if (status .gt. 0)call printerror(lun11,status)

      if (lpri.ne.0)
     $ write (lun11,*)'fstepr3: Add some more keywords'

C     Write some model parameters in the extension header
      call ftpcom(unit,'***********************************',status)
      if (status .gt. 0)call printerror(lun11,status)

      call ftpcom(unit,'Model Keywords',status)
      if (status .gt. 0)call printerror(lun11,status)

C     Write values to 3 decimal places
      rtmp=radin
      call ftpkye(unit,'RINNER',rtmp,3,'[cm] Inner shell radius',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=radout
      call ftpkye(unit,'ROUTER',rtmp,3,'[cm] Outer shell radius',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=rdel
      call ftpkye(unit,'RDEL',rtmp,3,'[cm] distance from face',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=t
      call ftpkye(unit,'TEMPERAT',rtmp,3,'[10**4K] Shell Temperature',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=prs
      call ftpkye(unit,'PRESSURE',rtmp,3,'[dynes/cm**2] Shell Pressure',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xcol
      call ftpkye(unit,'COLUMN',rtmp,3,'[/cm**2] Column ',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=xee
      call ftpkye(unit,'XEE',rtmp,3,'electron fraction',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xpx
      call ftpkye(unit,'DENSITY',rtmp,3,'[/cm**3] Density',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xi
      call ftpkye(unit,'LOGXI',rtmp,3,
     $ '[erg cm/s] log(ionization parameter)',status)
      if (status .gt. 0)call printerror(lun11,status)

      if (lpri.ne.0)
     $ write (lun11,*)'after header write'
C-------------------------------------------------------------------
C     Step through the columns and write them to the file
C
C     set 'global' parameters for writing FITS columns
      frow=1
      felem=1


C     column  1  (continuum index)
      colnum=1
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Writing Column ',colnum,nlpl
      nlines=kk
      call ftpclj(unit,colnum,frow,felem,nlines,ntptr,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  2  (level index)
      colnum=2
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Writing Column ',colnum,nlpl
      nlines=kk
      call ftpclj(unit,colnum,frow,felem,nlines,ntptr2,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  3  (wavelength)
      colnum=3
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,elsv,status)
      if (status .gt. 0)call printerror(lun11,status)
      if (status .gt. 0) return


C     column  4  (Ion)
      colnum=4
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Writing Column ',colnum
      call ftpcls(unit,colnum,frow,felem,nlines,kion,status)
      if (status .gt. 0)call printerror(lun11,status)


C     column  5 (lower Level Designation)
      colnum=5
      if (lpri.ne.0)
     $ write (lun11,*)'fstepr3: Writing Column ',colnum
      call ftpcls(unit,colnum,frow,felem,nlines,klevl,status)
      if (status .gt. 0)call printerror(lun11,status)

C     column  6 (Level Designation)
      colnum=6
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Writing Column ',colnum
      call ftpcls(unit,colnum,frow,felem,nlines,klevu,status)
      if (status .gt. 0)call printerror(lun11,status)

C----------------------------------------------------------------

      do ll=1,nlines
         rwrk1(ll)=0.
         if (ntptr(ll).ne.0) then
           rwrk3(ll)=cemab(1,ntptr(ll))
           rwrk4(ll)=cemab(2,ntptr(ll))
           rwrk5(ll)=cabab(ntptr(ll))
           rwrk6(ll)=opakab(ntptr(ll))
           rwrk7(ll)=tauc(1,ntptr(ll))
           rwrk8(ll)=tauc(2,ntptr(ll))
           endif
         enddo
c
c
C     column  7
      colnum=7
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk3,status)
      if (status .gt. 0)call printerror(lun11,status)

c
C     column  8
      colnum=8
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk4,status)
      if (status .gt. 0)call printerror(lun11,status)

c
C     column  9
      colnum=9
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk5,status)
      if (status .gt. 0)call printerror(lun11,status)

c
C     column  10
      colnum=10
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk6,status)
      if (status .gt. 0)call printerror(lun11,status)

c
C     column  11
      colnum=11
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr3: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk7,status)
      if (status .gt. 0)call printerror(lun11,status)

C     column  12
      colnum=12
      if (lpri.ne.0)
     $ write (lun11,*)'fstepr3: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk8,status)
      if (status .gt. 0)call printerror(lun11,status)



c----------------------------------------------------------------
C     Compute checksums
      call ftpcks(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

      return
      end
      subroutine fstepr4(unit,hdunum,radin,radout,rdel,t,prs,
     $             xcol,xee,xpx,xi,
     $             idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $             epi,ncn2,dpthc,opakc,rccemis,nloopctl,
     $             lun11,status)
C
C     Append a FITS extension binary table containing
C     nrhs columns and at most nrhdimj rows
c     author: T. Kallman
C
C     Parameters:
C        unit    integer            File unit number
C        hdunum  integer            Number of last HDU written
C        radin   real*8               inner radius of shell
C        radout  real*8               outer radius of shell
C                                   nb but now it is delr in the call
C        t    real*8               temperature of shell
C        prs    real*8               pressure in shell
C        nrhdimj  integer            Maximum number of rows
C        idat1   integer(nidat1)    Needed by the atomic database
C        rdat1   real(nidat1)       Needed by the atomic database
C        kdat1   char*nidat1        Needed by the atomic database
C        nptrs                      Needed by the atomic database
C        npnxt                      Needed by the atomic database
C        npfi                       Needed by the atomic database
C        npfirst                    Needed by the atomic database
C        npcon                      Needed by the atomic database
C        npconi                     Needed by the atomic database
C        npcon2                     Needed by the atomic database
C        xilev   real(nrhdimj)       Fractional level population array
C        cemab   real(2,nrhdimj)     Recombination emission
C        opakab  real(nrhdimj)       Opacity
C        tauc    real(2,nrhdimj)     Optical depth
C        poptol  real*8               Tolerance for population level
C        nloopctl integer           Loop control variable
C        nzone   integer            Pass number through iteration process
C        status  integer            Returned status code
C
      implicit none
      include './PARAM'
      integer mllz

      integer nptmpdim
      parameter (nptmpdim=ncn)
c
C     Allocation for passed parameters
      real*8 rdat1(nrdat1)
      real rtmp
      real*8 radin, radout,rdel, t, prs, xcol,xee,xpx,xi
      integer unit,hdunum, nrows, status, nloopctl
      integer idat1(nidat1),nptrs(nptt,ndat2)
c     energy bins
      real*8 epi(ncn)
c     continuum opacities
      real*8 opakc(ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn)
      integer ncn2

c     pointers to master data
      integer npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni),npar(ndat2)

C     Internal work areas
      real ntptr(nptmpdim)
      real rwrk1(nptmpdim)
      integer tfields,varidat
      character(16) ttype(5),tform(5),tunit(5)
      integer colnum,frow,felem,hdutype,ll, ltyp
      integer lrtyp, lcon, nrdt, nidt, mm, lun11, lpril,lpri
      integer jkk,nidti,nlines
      real*8 eliml,elimh,elin,elmmtpp,elcomp,eth,xeltp
      character(33) extname
      character(20) ktmp20
      character(20) ktmp2
      character(1) kdat1(nkdat1)

C     Database manipulation quantities
      integer nkdt
      character(1) kblnk,kdtmp(200)
      logical done

      data kblnk/' '/
c
      data tform/'1J','1E','1E','1E','1E'/

      data ttype/'index','energy','opacity','fwd dpth',
     $ 'bck dpth'/

      data tunit/' ','ev','/cm',' ',' '/

      varidat=0
c
      lpri=0
      lpril=lpri
c
      if (lpri.ne.0)
     $ write (lun11,*)'in fstepr4 input hdu',hdunum
c
C     Move to the last HDU (hdunum) in the file
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr4: Moving to end-of-FITS file'
      call ftmahd(unit,hdunum,hdutype,status)
      if (status .gt. 0)call printerror(lun11,status)

C     append a new empty extension after the last HDU
      if (lpri.ne.0)
     $ write (lun11,*)'fstepr4: Create the new extension'
      call ftcrhd(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

C----------------------------------------------------------------
C
C     Extracting data from the Atomic Database here
C
c      lpril=0
c

C     End of atomic database extraction
C----------------------------------------------------------------
C     define parameters for the binary table (see the above data statements)
      nrows=ncn2
      if (lpri.ne.0)
     $ write (lun11,*)'before header write'
      tfields=5
C     Build extension name
      extname='XSTAR_RADIAL'
      if(nloopctl.gt.0) then
          write(ktmp2,'(I4.4)')nloopctl
          extname='XSTAR_RADIAL_' // ktmp2
          endif

      if (lpri.ne.0)
     $ write (lun11,*)'fstepr4: Write table headers'
C     write the required header parameters for the binary table
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,
     $              varidat,status)
      if (status .gt. 0)call printerror(lun11,status)

      if (lpri.ne.0)
     $ write (lun11,*)'fstepr4: Add some more keywords'

C     Write some model parameters in the extension header
      call ftpcom(unit,'***********************************',status)
      if (status .gt. 0)call printerror(lun11,status)

      call ftpcom(unit,'Model Keywords',status)
      if (status .gt. 0)call printerror(lun11,status)

C     Write values to 3 decimal places
      rtmp=radin
      call ftpkye(unit,'RINNER',rtmp,3,'[cm] Inner shell radius',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=radout
      call ftpkye(unit,'ROUTER',rtmp,3,'[cm] Outer shell radius',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=rdel
      call ftpkye(unit,'RDEL',rtmp,3,'[cm] distance from face',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=t
      call ftpkye(unit,'TEMPERAT',rtmp,3,'[10**4K] Shell Temperature',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=prs
      call ftpkye(unit,'PRESSURE',rtmp,3,'[dynes/cm**2] Shell Pressure',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xcol
      call ftpkye(unit,'COLUMN',rtmp,3,'[/cm**2] Column ',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)

      rtmp=xee
      call ftpkye(unit,'XEE',rtmp,3,'electron fraction',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xpx
      call ftpkye(unit,'DENSITY',rtmp,3,'[/cm**3] Density',
     $ status)
      if (status .gt. 0)call printerror(lun11,status)
c
      rtmp=xi
      call ftpkye(unit,'LOGXI',rtmp,3,
     $ '[erg cm/s] log(ionization parameter)',status)
      if (status .gt. 0)call printerror(lun11,status)

      if (lpri.ne.0)
     $ write (lun11,*)'after header write'
C-------------------------------------------------------------------
C     Step through the columns and write them to the file
C
C     set 'global' parameters for writing FITS columns
      frow=1
      felem=1

      do mm=1,ncn2
        ntptr(mm)=mm
        enddo

C     column  1  (continuum index)
      colnum=1
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr4: Writing Column ',colnum
      nlines=ncn2
      call ftpclj(unit,colnum,frow,felem,nlines,ntptr,status)
      if (status .gt. 0)call printerror(lun11,status)
c
      do mm=1,ncn2
        rwrk1(mm)=epi(mm)
        enddo
c
C     column  2 energy
      colnum=2
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr4: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status)
      if (status .gt. 0)call printerror(lun11,status)
      if (status .gt. 0) return

c
      do mm=1,ncn2
        rwrk1(mm)=opakc(mm)
        enddo
c
C     column  3 opacity
      colnum=3
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr4: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status)
      if (status .gt. 0)call printerror(lun11,status)
      if (status .gt. 0) return

c
      do mm=1,ncn2
        rwrk1(mm)=dpthc(1,mm)
        enddo
c
C     column  4 depth forward
      colnum=4
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr4: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status)
      if (status .gt. 0)call printerror(lun11,status)
      if (status .gt. 0) return

c
      do mm=1,ncn2
        rwrk1(mm)=dpthc(2,mm)
        enddo
c
C     column  5 depth backward
      colnum=5
      if (lpri.ne.0)
     $ write(lun11,*)'fstepr4: Writing Column ',colnum
      call ftpcle(unit,colnum,frow,felem,nlines,rwrk1,status)
      if (status .gt. 0)call printerror(lun11,status)
      if (status .gt. 0) return




c----------------------------------------------------------------
C     Compute checksums
      call ftpcks(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

      return
      end
      subroutine func(lpri,lun11,vturbi,critf,
     $       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       decomp,ecomp,sxcomp,
     $       tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,nlevs,ncsvn,rates,vsav,idrates,
     $       xiin,rrrts,pirts,htt,cll,httot,cltot,hmctot,elcter,
     $         cllines,clcont,htcomp,clcomp,clbrems,
     $       xilevt,bilevt,rnist,nmat,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,fline,flinel, mhd_heat)

c
c     this routine steps through data and calculates
c     new version attempts to avoid rates for unabundant ions
c     author: T. Kallman
c
c     with data structures designed for Lucy's iterative method
c       nsup is a pointer from level n to superlevel N
c
c     no longer calls full func3 in main loop.
c     func3 calls moved to funcsyn
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
      real*8 zrems(4,ncn),zremso(4,ncn)
      real*8 elum(3,nnnl),elumo(3,nnnl)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
c     line optical depths
      real*8 tau0(2,nnnl)
c     energy bins
      real*8 epi(ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn),brcems(ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
      real*8 fline(2,nnnl),flinel(ncn)
c     level populations
      real*8 xilevt(nnml),bilevt(nnml),rnist(nnml)
c     ion abundances
      real*8 xiin(nni)
      real*8 rrrts(nni),pirts(nni)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
      real*8 elumab(2,nnml),elumabo(2,nnml)
      real*8 tauc(2,nnml)
      real*8 rlev(10,nd)
      integer ilev(10,nd),nlpt(nd),iltp(nd)
      real*8 htt(nni),cll(nni)
      integer nlevs(nni)
c     element abundances
      real*8 abel(nl)
c     the saved rates
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
      real*8 decomp(ncomp,ncomp),ecomp(ncomp),sxcomp(ncomp)
c
      character(1) klev(100,nd)
c
      real*8 rrrt(31),pirt(31),xin(31),xitmp(31)
      real*8 ajisb(2,ndb),cjisb(ndb)
      integer indb(2,ndb)
      real*8 xilev(nd),rniss(nnml)
      real*8 rrcor(nni),pirtt(31)
      real*8 bmat(nd),bmatl(nd)
      real*8 rnisl(nd)
      real*8 x(nd)
      integer ipsv(31),nsup(nd)
      integer lpri,lun11,lcdd,ncn2,np2,ncsvn,nmat,nlsvn
      real*8 vturbi,critf,t,trad,r,delr,xee,xpx,cfrac,p,
     $     hmctot,elcter,cllines,clcont,htcomp,clcomp,clbrems
      real*8 xh1,xh0,httot,cltot
      real*8 rnisum,crith,cltmp,cmp1,cmp2,
     $     cltot2,enelec,httot2,httmp,pirtsum,rniss2,rrrtt,
     $     rtdm,tt1,tt2,xintp,xeltp,ximax,xilast,xintp2,
     $     xisum,xipp,cl,ht
      real*8 httotd,cltotd,hmctotd,elcterd,
     $     htcompd,clcompd,clbremsd
      integer nlev,nindb,lprisv,
     $     jkk,ipmat,ltyp,ldir,llp,imax,ilimh,
     $     lrtyp,lcon,nrdt,nidt,nkdt,ll,jkkl,ipmatsv,
     $     iliml,jk,kl1,mm,kl,kl2,klion,klel,klp,llm,lp,lm,
     $     lprim,lprif,lpril,lpritp,lsum,lsumt,ndtmp,mlel,
     $     ml1,mmt,mllel,mlion,mleltp,mmtmp,nit,nit2,nit3,nitmx,
     $     nitmx2,nlevm,nnz,nnzp,nsp,mlm,np1i,np1r,np1k
c
      real*8 mhd_heat
c
      lprisv=lpri
      lprif=lpri
      lpritp=0
c      if (lpri.ge.1) lpritp=2
      if (lprif.ne.0)
     $  write (lun11,*)'in func, inputs:',t,
     $         xee,xpx,lcdd,p,abel(1),delr
       if (lcdd.ne.1)
     $   xpx = p/1.38e-12/max(t,1.e-24)
c
      xh0=xpx*xiin(1)*abel(1)
      xh1=xpx*(1.-xiin(1))*abel(1)
c
c      zero emissivitiesd and opacities
       do ll=1,nni
         htt(ll)=0.
         cll(ll)=0.
         xiin(ll)=0.
         enddo
       do ll=1,29
         rrrt(ll)=0.
         pirt(ll)=0.
         enddo
       do ll=1,nnml
         xilevt(ll)=0.
         bilevt(ll)=0.
         cemab(1,ll)=0.
         cemab(2,ll)=0.
         opakab(ll)=0.
         enddo
       elcter=0.
       httot=mhd_heat
       cltot=0.
       clcont=0.
       cllines=0.
       do ll=1,nnnl
         rcem(2,ll)=0.
         rcem(1,ll)=0.
         oplin(ll)=0.
         fline(1,ll)=0.
         fline(2,ll)=0.
         enddo
       do ll=1,ncn2
         rccemis(1,ll)=0.
         rccemis(2,ll)=0.
         opakc(ll)=0.
         opakscatt(ll)=0.
         enddo
       do ll=1,ncn2
         flinel(ll)=0.
         enddo
c
c
c      now calculate.  first step thru elements
       jkk=0
       jkkl=0
       klel=11
       mlel=npfirst(klel)
       jk=0
       xilast=0.
       httot=mhd_heat
       cltot=0.
       do while (mlel.ne.0)
         mlm=mlel-1
         call drd(ltyp,lrtyp,lcon,
     $     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $     nptrs,0,lun11)
         mllel=0
         if (nidt.gt.0) then
           if (lprif.ne.0)
     $         write (lun11,9339)(kdat1(np1k-1+mm),mm=1,nkdt)
 9339      format (1x, ' element:',12a1)
           mllel=idat1(np1i)
           jk=mllel
           nnz=idat1(np1i)
           nnzp=nnz+1
           xeltp=0.
           if (jk.gt.0) xeltp=abel(jk)
           if (xeltp.gt.1.e-24) then
c
c            now step thru ions first pass: func1
             if (lprif.ne.0) write (lun11,*)' first pass'
             klion=12
             mlion=npfirst(klion)
             jkk=0
             kl=0
             do while ((mlion.ne.0).and.(kl.lt.nnz))
               jkk=jkk+1
               mlm=mlion-1
               call drd(ltyp,lrtyp,lcon,
     $           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $           nptrs,0,lun11)
               mleltp=npar(mlion)
               if (mleltp.eq.mlel) then
                 kl=kl+1
                 call func2l(jkk,lpritp,lun11,t,xee,xpx,
     $                  idat1,rdat1,kdat1,nptrs,
     $                  npar,npnxt,npfi,
     $                  rniss,rlev,ilev,
     $                  nlpt,iltp,nlev,klev)
                 if (lprif.ne.0)
     $            write (lun11,9338)(kdat1(np1k-1+mm),mm=1,nkdt)
9338             format (1x, ' ion:',8a1)
                 if (lprif.ne.0) write (lun11,9328)
                 call func1(jkk,kl,nnz,
     $               lpri,lun11,vturbi,
     $               t,trad,r,xee,xpx,xh1,xh0,
     $               epi,ncn2,bremsa,bremsint,
     $               idat1,rdat1,kdat1,nptrs,np2,
     $               npar,npnxt,npfi,npfirst,
     $               nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $               npconi2,ncsvn,rates,vsav,idrates,
     $               rniss,rlev,ilev,
     $               nlpt,iltp,nlev,klev,
     $               pirtt,rrrtt)
                 rrcor(jkk)=1.
                 pirtsum=0.
                 do mm=kl,nnzp
                   pirtsum=pirtsum+pirtt(mm)
                   enddo
                 pirt(kl)=pirtsum
                 rrrt(kl)=rrrtt
                 endif
               mlion=npnxt(mlion)
               enddo
c
c            do ion balance
             kl2=kl
             call istruc(pirt,rrrt,xin,kl2,lpritp,lun11)
             xisum=0.
             ldir=+1
             do mm=1,kl
                 kl1=mm
                 xiin(jkk-kl+kl1)=xin(kl1)
                 xisum=xisum+xin(kl1)
                 enddo
             xipp=1.-xisum
             klp=kl+1
             xin(klp)=xipp
             iliml=0
             ilimh=0
             imax=0
             ximax=0.
             do mm=1,klp
               if (xin(mm).gt.ximax) then
                 ximax=xin(mm)
                 imax=mm
                 endif
               enddo
             imax=max(min(nnz,imax),1)
             llp=imax
             llm=imax
             iliml=imax
             ilimh=imax
             lp=0
             lm=0
             if (imax.ne.klp) then
                 lsumt=nlevs(jkk-kl+imax)
               else
                 lsumt=0
               endif
             ndtmp=nd
             mmt=0
             do while ((lsumt.lt.ndtmp).and.(ldir.ne.0))
               mmt=mmt+1
               lsum=lsumt
               iliml=min(iliml,llm)
               ilimh=max(ilimh,llp)
               if ((llp.ge.klp).or.(xin(llp).lt.critf)) then
                 ldir=-1
                 lp=1
                 endif
               if ((llm.le.1).or.(xin(llm).lt.critf)) then
                 ldir=+1
                 lm=1
                 endif
               if ((lm.ne.1).and.(lp.ne.1)) then
                 if (xin(llp+1).gt.xin(llm-1)) then
                     ldir=+1
                   else
                     ldir=-1
                   endif
                 endif
               if ((lp.eq.1).and.(lm.eq.1)) ldir=0
               if (ldir.eq.+1) then
                   llp=llp+1
                   if (llp.ne.klp) then
                       lsumt=lsum+nlevs(jkk-kl+llp)
                     else
                       lsumt=lsum
                     endif
                   endif
               if (ldir.eq.-1) then
                   llm=llm-1
                   lsumt=lsum+nlevs(jkk-kl+llm)
                   endif
               ilimh=max(ilimh-1,iliml+1)
               enddo
             if (lprif.ne.0) then
                 write (lun11,*)'ion fractions:',iliml,ilimh,lsum
                 write (lun11,*)'ion, pi rate,    rec rate,   fraction'
                 do mm=1,kl
                   write (lun11,9023)mm,pirt(mm),rrrt(mm),xin(mm)
 9023              format (1x,i4,3(1pe10.2))
                   enddo
                 endif
c
c
c            now step thru ions for second pass
             if (lprif.ne.0) write (lun11,*)' second pass'
 9328        format ('     ion    process     d1    d2 ',
     $         '    rec use    ans1      ans2  ',
     $     '   ionization  recombination')
 9329        format ('     ion    process     d1    d2 ',
     $  '             rec use     ans1      ans2     ans3       ans4  ',
     $         '  aji(lo,up) aji(up,lo)',
     $         ' aji(lo,lo) aji(up,up)')
 9330        format ('     ion    process     d1    d2 ',
     $   '             rec use     ans1      ans2     ans3       ans4 ',
     $      '  emiss.    opac.     energy   rec   heat      cool    ')
c
             if (lprif.ne.0) write (lun11,*)'zeroing:',nmat
             do ml1=1,nmat
               bmat(ml1)=0.
               nsup(ml1)=0
               enddo
             klion=12
             mlion=npfirst(klion)
             jkk=0
             kl=0
             nindb=0
             ipmat=0
             nsp=1
             do while ((mlion.ne.0).and.(kl.lt.nnz))
               jkk=jkk+1
               mlm=mlion-1
               call drd(ltyp,lrtyp,lcon,
     $           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $           nptrs,0,lun11)
               mleltp=npar(mlion)
               if (mleltp.eq.mlel) then
                 kl=kl+1
                 ipsv(kl)=-1
                 if (lprif.ne.0)
     $             write (lun11,9338)(kdat1(np1k-1+mm),mm=1,nkdt)
                 if (lprif.ne.0) write (lun11,9329)
                 if ((kl.ge.iliml).and.(kl.le.ilimh)) then
                   call func2l(jkk,lpritp,lun11,t,xee,xpx,
     $                  idat1,rdat1,kdat1,nptrs,
     $                  npar,npnxt,npfi,
     $                  rniss,rlev,ilev,
     $                  nlpt,iltp,nlev,klev)
                   if (lprif.ne.0) write (lun11,*)'ipmat=',ipmat
                   if (lprif.ne.0) write (lun11,*)'before func2',nindb
                   ipsv(kl)=ipmat
                   rrrtt=rrrt(kl)
                   call func2(jkk,kl,ilimh,
     $                   lpri,lun11,vturbi,
     $                   t,trad,r,xee,xpx,xh1,xh0,cfrac,
     $                   epi,ncn2,bremsa,bremsint,tau0,tauc,
     $                   idat1,rdat1,kdat1,nptrs,np2,
     $                   npar,npnxt,npfi,npfirst,
     $                   nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $                   npconi2,ncsvn,rates,vsav,idrates,
     $                   rniss,rlev,ilev,nmat,
     $                   nlpt,iltp,nlev,klev,ajisb,cjisb,indb,
     $                   rrrtt,ipmat,nindb,
     $                   rcem,oplin,opakc,opakscatt,
     $                   cemab,cabab,opakab,fline,flinel)
c                   call func2a(jkk,kl,ilimh,
c       $                 lpri,lun11,lfpi,vturbi,
c       $                 t,trad,r,xee,xpx,xh1,xh0,
c       $                 epi,ncn2,bremsa,bremsint,
c       $                 idat1,rdat1,kdat1,nptrs,np2,
c       $                 npar,npnxt,npfi,npfirst,
c       $                 nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
c       $                 npconi2,ncsvn,rates,vsav,idrates,
c       $                 rlev,ilev,
c       $                 nlpt,iltp,nlev,klev,
c       $                 ajis,cjis,ipmat,ipsv)
c
c                  condense the reaction matrix to omit levels without type 7 data
c                    some definition:
c                       nlev is the number of levels in the ion (including cont.)
c                       nmat is the maximum index of levels linked by this ion
c                       ipmat is the base index into the master arrays for this ion
c                          i.e. the array elements for  this ion start at ipmat+1
c                          and go to ipmat+nmat
c
                   if (lpri.ge.2) write (lun11,*)'level mapping:'
                   do mm=1,nmat
                     if (mm.le.nlev) then
                       mmtmp=npilev(mm,jkk)
                       if (mmtmp.gt.0) then
                         x(mm+ipmat)=xilevt(mmtmp)
                         endif
                       endif
                     enddo
c
c                  set up superlevel pointers
                   nsup(1+ipmat)=nsp
                   nsp=nsp+1
                   do mm=2,nlev-1
                     nsup(mm+ipmat)=nsp
                     enddo
                   nsp=nsp+1
                   ipmat=ipmat+nlev-1
                   if (ipmat.gt.nd) stop 'ipmat too large.
     $                                    Increase critf value.'
                   endif
                 endif
               mlion=npnxt(mlion)
               enddo
c
c
c
             nsup(ipmat+1)=nsp
             x(ipmat+1)=0.
             ipmat=ipmat+1
             nmat=ipmat
c
c            now calculate populations:  second pass, full list
             call remtms(tt1)
             nitmx=200
             nitmx2=200
             lprim=0
c             if (lpri.ne.0) lprim=4
             if (lpri.ne.0)
     $         write (lun11,*)'before msolvelucy',ipmat
             call msolvelucy(ajisb,cjisb,indb,nindb,nsup,nsp,ipmat,
     $          bmat,x,ht,cl,nit,nit2,nit3,nitmx,nitmx2,lun11,lprim)
             if (lprim.ge.2)
     $         call chisq(ajisb,cjisb,indb,nindb,
     $                    ipmat,x,lun11,lpri)
             call remtms(tt2)
             if (lpri.gt.0)
     $        write (lun11,981)abs(tt2-tt1),nit,nit2,nit3,ht,cl
 981         format (1x,'after msolvelucy',(1pe11.3),3i4,2(1pe11.3))
             do mm=1,ipmat
               bmatl(mm)=x(mm)
               enddo
             cltot=cltot+cl*xeltp
             httot=httot+ht*xeltp
c
c            now calculate lte abundances
c            step thru ions
             if (lprif.ne.0) write (lun11,*)'calculating lte abundances'
             ipmatsv=ipmat
             ipmat=0
             mlion=npfirst(klion)
             jkk=jkk-nnz
             kl=0
             do while ((mlion.ne.0).and.(kl.lt.nnz))
               ltyp=klion
               mlm=mlion-1
               call drd(ltyp,lrtyp,lcon,
     $           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $           nptrs,0,lun11)
               mleltp=npar(mlion)
               if (mleltp.eq.mlel) then
                 kl=kl+1
                 jkk=jkk+1
                 if ((kl.ge.iliml).and.(kl.le.ilimh)) then
                   call func2l(jkk,lpri,lun11,t,xee,xpx,
     $                  idat1,rdat1,kdat1,nptrs,
     $                  npar,npnxt,npfi,
     $                  rniss,rlev,ilev,
     $                  nlpt,iltp,nlev,klev)
                   nlevm=nlev-1
                   rnisum=0.
                   do mm=1,nlevm
                     rnisum=rnisum+rniss(mm)
                     enddo
                   rniss2=rniss(nlev)
                   rrrt(kl)=pirt(kl)*rnisum/(rniss2+1.e-28)
                   endif
                 endif
               mlion=npnxt(mlion)
               enddo
             call istruc(pirt,rrrt,xitmp,kl2,lpritp,lun11)
             mlion=npfirst(klion)
             jkk=jkk-nnz
             kl=0
             do while ((mlion.ne.0).and.(kl.lt.nnz))
               ltyp=klion
               mlm=mlion-1
               call drd(ltyp,lrtyp,lcon,
     $           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $           nptrs,0,lun11)
               mleltp=npar(mlion)
               if (mleltp.eq.mlel) then
                 kl=kl+1
                 jkk=jkk+1
                 if ((kl.ge.iliml).and.(kl.le.ilimh)) then
                   call func2l(jkk,lpri,lun11,t,xee,xpx,
     $                  idat1,rdat1,kdat1,nptrs,
     $                  npar,npnxt,npfi,
     $                  rniss,rlev,ilev,
     $                  nlpt,iltp,nlev,klev)
c                  func2l returns the saha levels relative to the continuum
c                  in rniss.  we calculate the lte populations..
                   nlevm=nlev-1
                   if (kl.eq.nnz) nlevm=nlev
                   do mm=1,nlevm
                     rnisl(mm+ipmat)=rniss(mm)
     $                    *(xitmp(kl)+xitmp(kl+1))
                     enddo
                   ipmat=ipmat+nlev-1
                   endif
                 endif
               mlion=npnxt(mlion)
               enddo
c
c            step thru ions
             if (lprif.ne.0) write (lun11,*)' third pass',ipmat
             mlion=npfirst(klion)
             ipmat=ipmat+1
             ipmatsv=ipmat
             ipmat=0
             jkk=jkk-nnz
             kl=0
             do while ((mlion.ne.0).and.(kl.lt.nnz))
               ltyp=klion
               mlm=mlion-1
               call drd(ltyp,lrtyp,lcon,
     $           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $           nptrs,0,lun11)
               mleltp=npar(mlion)
               if (mleltp.eq.mlel) then
                 kl=kl+1
                 jkk=jkk+1
                 call func2i(jkk,
     $             idat1,rdat1,kdat1,nptrs,
     $             npfi,npar,npnxt,nlev)
c                this loop avoids a messy error
                 do mm=1,nlev
                   mmtmp=npilev(mm,jkk)
                   if (mmtmp.gt.0) then
                     xilevt(mmtmp)=0.
                     if (lpri.gt.1)
     $                write (lun11,*)mm,mmtmp,xilevt(mmtmp),ipmat,
     $                               mm+ipmat
                     endif
                   enddo
                 pirts(jkk)=pirt(kl)
                 rrrts(jkk)=rrrt(kl)
                 if ((kl.ge.iliml).and.(kl.le.ilimh)) then
                   if (lprif.ne.0)
     $              write (lun11,9338)(kdat1(np1k-1+mm),mm=1,nkdt)
                   if (lprif.ne.0) write (lun11,9330)
                   call func2l(jkk,lpri,lun11,t,xee,xpx,
     $                  idat1,rdat1,kdat1,nptrs,
     $                  npar,npnxt,npfi,
     $                  rniss,rlev,ilev,
     $                  nlpt,iltp,nlev,klev)
                   nlevm=nlev-1
c                  retrieve saved abundances
                   xintp=0.
                   if (lpri.gt.1) write (lun11,*)'saving populations'
c
c                  nb this code makes H and He fully ionized
c                   if ((jkk.eq.1).or.(jkk.eq.3))
c     $               bmatl(nlev+ipmat)=1.
c
                   do mm=1,nlevm
c
c                    nb this code makes H and He fully ionized
c                     if (jkk.le.3) then
c                       bmatl(mm+ipmat)=0.
c                       endif
c
                     xilev(mm)=bmatl(mm+ipmat)
                     xintp=xintp+xilev(mm)
                     mmtmp=npilev(mm,jkk)
                     if (mmtmp.gt.0) then
                       rnist(mmtmp)=rnisl(mm+ipmat)
                       if (mmtmp.gt.nnml) stop 'mmtmp error'
                       xilevt(mmtmp)=bmatl(mm+ipmat)
                       if (lpri.gt.1)
     $                  write (lun11,*)mm,mmtmp,xilevt(mmtmp)
                       endif
                     enddo
                   xiin(jkk)=xintp
                   xin(kl)=xintp
                   mmtmp=npilev(nlev,jkk)
                   if (mmtmp.gt.nnml) stop 'mmtmp error'
                   xilevt(mmtmp)=bmatl(ipmat+nlev)
                   rnist(mmtmp)=rnisl(ipmat+nlev)
c                  this is an ungraceful solution to this problem
                   xintp2=bmatl(ipmat+nlev)
c                  this is a bad approxmiation...
c                   if (kl.lt.nnz)
c       $              xintp2=xintp2
c       $                 +bmatl(ipmat+nlev+1)
c       $                 +bmatl(ipmat+nlev+2)
                   xilev(nlev)=xintp2
c
                   rrrts(jkk)=pirts(jkk)*xintp/(xintp2+1.e-28)
                   if (lpri.ne.0) then
                     write (lun11,*)'ipmat=',ipmat
                   call func3p(jkk,jkkl,lpri,lun11,vturbi,
     $                 t,trad,r,xee,xpx,xh1,xh0,cfrac,
     $                 epi,ncn2,bremsa,bremsint,tau0,tauc,
     $                 idat1,rdat1,kdat1,nptrs,np2,
     $                 npar,npnxt,npfi,npfirst,
     $                 nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $                 npconi2,ncsvn,rates,vsav,idrates,
     $                 rniss,rlev,ilev,
     $                    nlpt,iltp,nlev,klev,
     $                 xeltp,rrcor(jkk),httmp,cltmp,cllines,clcont,rtdm,
     $                 bmatl,ipmat,ipmatsv,
     $                 rcem,oplin,rccemis,opakc,opakscatt,
     $                 cemab,cabab,opakab,fline,flinel)
                     endif
c
                   ipmat=ipmat+nlev-1
                   endif
                 endif
c
               mlion=npnxt(mlion)
               enddo
c
             xisum=0.
             do kl1=1,kl
               xisum=xisum+xin(kl1)
               enelec=float(kl1-1)
               elcter=elcter+xin(kl1)*enelec*xeltp
               enddo
             enelec=float(kl)
             elcter=elcter+max(0.,1.-xisum)*enelec*xeltp
c
             endif

           endif
c
         if  (mlel.ne.0) mlel=npnxt(mlel)
c
         enddo
c
      lpril=0
c      call comp2(lpril,lun11,epi,ncn2,bremsa,t,
c     $   r,decomp,ecomp,sxcomp,cmp1,cmp2)
c      call comp(lpri,lun11,epi,ncn2,bremsa,r,cmp1,cmp2)
      call comp3(lpri,lun11,epi,ncn2,bremsa,r,cmp1,cmp2)
c     call freef(lpri,lun11,epi,ncn2,t,xpx,xee,opakc)
c     call bremem(lpril,lun11,xee,xpx,t,epi,ncn2,brcems,opakc)
      call heatf(jkk,lpri,lun11,
     $       t,r,cfrac,delr,xee,xpx,abel,
     $       epi,ncn2,bremsa,
     $       idat1,rdat1,kdat1,nptrs,
     $       npar,npnxt,npfi,npfirst,nplin,nlsvn,
     $       npconi2,ncsvn,rlev,ilev,nlev,klev,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       rcem,oplin,rccemis,opakc,opakscatt,cemab,fline,flinel,
     $       brcems,cmp1,cmp2,httot,cltot,hmctot,
     $             cllines,clcont,htcomp,clcomp,clbrems, mhd_heat)
c
       elcter=xee-elcter
c
      if (lprif.ne.0) write (lun11,*)'leaving func'
c
      lprisv=lpri
c
c
      return
      end
      subroutine func1(jkk,kl,nnz,lpri,lun11,vturbi,
     $       t,trad,r,xee,xpx,xh1,xh0,
     $       epi,ncn2,bremsa,bremsint,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,ncsvn,rates,vsav,idrates,
     $       rniss,rlev,ilev,
     $          nlpt,iltp,nlev,klev,
     $       pirt,rrrt2)
c
c     this routine calculates rates affecting ion balance
c     author: T. Kallman
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     energy bins
      real*8 epi(ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rcdum(2,ncn)
      real*8 fline(2,nnnl)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     element abundances
c     the saved rates
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
c     state variables
      real*8 r,t,xpx
c     heating-cooling variables
c     input parameters
      real*8 trad
      real*8 vturbi,xee
      integer ncn2,lpri,lun11,np2
      integer nlsvn,ncsvn
      real*8 rniss(nd)
      real*8 rlev(10,nd)
      integer ilev(10,nd),nlpt(nd),iltp(nd)
      real*8 pirt(31)
      character(1) klev(100,nd)
      character(49) kdesc2
      real*8 tsq,ans1,ans2,xh1,xh0,rrrt2
      real*8 abund1,abund2,ptmp1,ptmp2,ans3,ans4,opakb1
      integer idest1,idest2,idest3,idest4
      integer np1i,np1r,np1k
      integer nnzp,nlevmx,mltype,ml,mllz,nlev,mlpar,
     $  ltyp,lrtyp,lcon,nrdt,nidt,nkdt,mlrdesc,llo,lup,
     $  nnz,mm,jkk,lk,lpriu,kl,mlm,lfpi
c
c      if (lpri.ne.0)
c     $  write (lun11,*)'in func1, inputs:',t,
c     $         xee,xpx,xnx
c
c     lfpi value:  calculate photoionization rates only
      lfpi=1
c
c     zero temporaries
      tsq=sqrt(t)
      nnzp=nnz+1
      do mm=1,nnzp
        pirt(mm)=0.
        enddo
      rrrt2=0.
c     now find all the rates affecting this ion
c     step thru types
      mltype=1
      do while (mltype.lt.ntyp)
        mlrdesc=mltype
        ml=npfi(mltype,jkk)
        if (((mlrdesc.eq.6).or.(mlrdesc.eq.1).or.(mlrdesc.eq.8)
     $    .or.(mlrdesc.eq.15).or.(mlrdesc.eq.7).or.(mlrdesc.eq.42))
     $    .and.(ml.ne.0)) then
          mllz=npar(ml)
          mlpar=npar(ml)
          do while ((ml.ne.0).and.(mlpar.eq.mllz))
c           step thru records of this type
            if (nptrs(3,ml).eq.mlrdesc) then
              mlm=ml-1
              call drd(ltyp,lrtyp,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $          nptrs,0,lun11)
c              if (lpri.ne.0)  call dprinto(ltyp,lrtyp,lcon,
c     $          nrdt,np1r,nidt,np1i,nkdt,np1k,
c     $          rdat1,idat1,kdat1,lun11)
c             calculate rates
              lpriu=lpri
              abund1=0.
              abund2=0.
              ptmp1=0.
              ptmp2=0.
              idest1=0
              if (lrtyp.eq.7) idest1=idat1(np1i+nidt-2)
              if (((lrtyp.eq.1).and.(ltyp.ne.53))
     $             .or.((lrtyp.eq.7).and.(idest1.eq.1))
     $             .or.(lrtyp.eq.15).or.(lrtyp.eq.42)) then
                 call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $             nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $             ans3,ans4,idest1,idest2,idest3,idest4,
     $             abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $             opakc,opakscatt,rcdum,fline,lpriu,kdesc2,
     $             r,t,trad,tsq,xee,xh1,xh0,
     $             epi,ncn2,bremsa,bremsint,
     $             rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $             idat1,rdat1,kdat1,nptrs,np2,
     $             npar,npnxt,npfi,npfirst,
     $             nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $             npconi2,ncsvn,rates,vsav,idrates)
                 if (idest1.eq.1) then
                   llo=idest3
                   lup=idest4-idest3
                   pirt(kl+lup)=pirt(kl+lup)+ans1
                   if (lpri.ge.1)
     $              write (lun11,9001)jkk,lrtyp,ltyp,llo,lup+idest3,
     $                      ml,ans1,pirt(kl+lup)
 9001               format (1x,6i6,' ion ',1pe10.3,14x,1pe10.3)
                   endif
                 endif
              if ((lrtyp.eq.6).or.(lrtyp.eq.8)) then
                 call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $             nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $             ans3,ans4,idest1,idest2,idest3,idest4,
     $             abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $             opakc,opakscatt,rcdum,fline,lpriu,kdesc2,
     $             r,t,trad,tsq,xee,xh1,xh0,
     $             epi,ncn2,bremsa,bremsint,
     $             rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $             idat1,rdat1,kdat1,nptrs,np2,
     $             npar,npnxt,npfi,npfirst,
     $             nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $             npconi2,ncsvn,rates,vsav,idrates)
                 rrrt2=rrrt2+ans1
c                Commented... llo has a funny value...    jg
                 llo=1
                 if (lpri.ge.1)
     $              write (lun11,9002)jkk,lrtyp,ltyp,llo,ml,
     $                    ans1,rrrt2
 9002            format (1x,4i6,6x,i6,' ion ',14x,
     $                              1pe10.3,14x,1pe10.4)
                 endif
              endif
            ml=npnxt(ml)
            if (ml.ne.0) mlpar=npar(ml)
            enddo
          endif
        mltype=mltype+1
        enddo
c
c
      return
      end
      subroutine func2(jkk,kl,ilimh,lpriz,lun11,vturbi,
     $       t,trad,r,xee,xpx,xh1,xh0,cfrac,
     $       epi,ncn2,bremsa,bremsint,tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,ncsvn,rates,vsav,idrates,
     $       rniss,rlev,ilev,nmat,
     $          nlpt,iltp,nlev,klev,ajisb,cjisb,indb,
     $         rrrtot,ipmat,nindb,
     $       rcem,oplin,opakc,opakscatt,
     $       cemab,cabab,opakab,fline,flinel)
c
c     this routine calculates rates affecting level populations
c     author: T. Kallman
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     energy bins
      real*8 epi(ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rcdum(2,ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line optical depths
      real*8 tau0(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
      real*8 fline(2,nnnl),flinel(ncn)
c     element abundances
c     the saved rates
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
c     state variables
      real*8 r,t,xpx
c     heating-cooling variables
c     input parameters
      real*8 trad
      real*8 vturbi,xee
      integer ncn2,lpri,lun11,lfpi,np2
      integer nlsvn,ncsvn
      real*8 rlev(10,nd)
      integer ilev(10,nd),nlpt(nd),iltp(nd)
      character(1) klev(100,nd)
      character(49) kdesc2
      real*8 tsq,ans1,ans2,xh1,xh0,cfrac,rrrtot
      real*8 abund1,abund2,ptmp1,ptmp2,ans3,ans4,opakb1
      integer idest1,idest2,idest3,idest4,ipmat,nindb
c     continuum flux
      real*8 tauc(2,nnml)
c
      character(1) kblnk
      real*8 rniss(nnml)
      real*8 ajisb(2,ndb),cjisb(ndb)
      integer indb(2,ndb)
      real*8 rrrt3,tau1,tau2,airtmp,rrrtot2,e1,e2,pescl,pescv,ptmp,eth
      integer np1i,np1r,np1k
      integer nlev,lpriz,nlevmx,mltype,ml,mllz,mlrdesc,lpriu,
     $        llo,lup,ilimh,jkk,kl,nmat,ml1,ltyp,
     $        lrtyp,lcon,nrdt,nidt,nkdt,lk,ll,kkkl,jkkl,mlpar,mlm,mm
c
      data kblnk/' '/
c
c
c      if (lpri.gt.1)
c     $  write (lun11,*)'in func2, inputs:',t,
c     $         xee,xpx
c
c     lfpi value: photoionization and recombination, no opacities
      lfpi=2
c
      lpri=lpriz
      tsq=sqrt(t)
c     zero temporaries
      rrrt3=rrrtot
      rrrtot=0.
      rrrtot2=0.
c
c     now find the rates affecting this ion
      mltype=0
      do while (mltype.lt.ntyp)
        mltype=mltype+1
        mlrdesc=mltype
        if (((mlrdesc.lt.10).or.(mlrdesc.gt.13))
     $    .and.(mlrdesc.gt.0)) then
          ml=npfi(mltype,jkk)
          if (ml.gt.0) then
            mllz=npar(ml)
c           step thru records of this type
            mlpar=npar(ml)
            do while ((ml.ne.0).and.(mlpar.eq.mllz))
              if (nptrs(3,ml).eq.mlrdesc) then
                mlm=ml-1
                call drd(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $            nptrs,0,lun11)
c               calculate rates
                lpriu=lpri
                abund1=0.
                abund2=0.
                ptmp2=1.
                ptmp1=0.
                if ((lrtyp.eq.7).or.((lrtyp.eq.1).and.(ltyp.ne.53)))
     $            then
                  kkkl=npconi2(ml)
                  if ((kkkl.gt.0).and.(kkkl.le.ndat2)) then
                    tau1=tauc(1,kkkl)
                    tau2=tauc(2,kkkl)
                    ptmp1=pescv(tau1)*(1.-cfrac)
                    ptmp2=pescv(tau2)*(1.-cfrac)
     $                  +2.*pescv(tau1+tau2)*cfrac
                    ptmp=(ptmp1+ptmp2)
c                   note that radiative rates and emissivities will
c                     have the escape probabilities in them
c                     when calculated in ucalc
                    call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $                nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $                ans3,ans4,idest1,idest2,idest3,idest4,
     $                abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $                opakc,opakscatt,rcdum,fline,lpriu,kdesc2,
     $                r,t,trad,tsq,xee,xh1,xh0,
     $                epi,ncn2,bremsa,bremsint,
     $                rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $                idat1,rdat1,kdat1,nptrs,np2,
     $                npar,npnxt,npfi,npfirst,
     $                nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $                npconi2,ncsvn,rates,vsav,idrates)
c                   this statement prevents double counting of PI from
c                     excited levels.  All rate type 1 data should not have
c                     lower level associated with it, but some does.
                    if ((lrtyp.eq.1).and.(idest1.ne.1)) ans1=0.
                    cabab(kkkl)=ans3
                    opakab(kkkl)=opakb1
                    cemab(1,kkkl)=ptmp1*ans4/(ptmp1+ptmp2)
                    cemab(2,kkkl)=ptmp2*ans4/(ptmp1+ptmp2)
                    if (lrtyp.eq.1) then
                      ans2=0.
                      ans4=0.
                      endif
                    llo=idest1
                    lup=idest2
                    if (kl.eq.ilimh) lup=min(nlev,lup)
                    if ((idest1.gt.0).and.(idest2.gt.0).and.
     $                (llo.le.nd).and.(lup.le.nd)) then
                      eth=rlev(4,idest1)-rlev(1,idest1)
                      airtmp=ans2
                      rrrt3=rrrt3-airtmp*max(0.,1.-ptmp)
                      nindb=nindb+1
                      if (nindb.gt.ndb) stop 'array indexing error'
                      ajisb(1,nindb)=ans1
                      ajisb(2,nindb)=airtmp
                      cjisb(nindb)=0.
                      indb(1,nindb)=lup+ipmat
                      indb(2,nindb)=llo+ipmat
                      if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                      nindb=nindb+1
                      if (nindb.gt.ndb) stop 'array indexing error'
                      ajisb(1,nindb)=airtmp
                      ajisb(2,nindb)=ans1
                      cjisb(nindb)=0.
                      indb(1,nindb)=llo+ipmat
                      indb(2,nindb)=lup+ipmat
                      nindb=nindb+1
                      if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                      if (nindb.gt.ndb) stop 'array indexing error'
                      ajisb(1,nindb)=-ans1
                      ajisb(2,nindb)=-ans1
                      cjisb(nindb)=-ans3*xpx
                      indb(1,nindb)=llo+ipmat
                      indb(2,nindb)=llo+ipmat
                      if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                      nindb=nindb+1
                      if (nindb.gt.ndb) stop 'array indexing error'
                      ajisb(1,nindb)=-airtmp
                      ajisb(2,nindb)=-airtmp
                      cjisb(nindb)=ans4*xpx
                      indb(1,nindb)=lup+ipmat
                      indb(2,nindb)=lup+ipmat
                      if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                      nmat=max(nmat,llo)
                      nmat=max(nmat,lup)
                      rrrtot=rrrtot+airtmp
                      if (idest2.le.nlev)
     $                  rrrtot2=rrrtot2+airtmp
                      if ((lpri.ge.1))
     $                  write (lun11,9004)jkk,lrtyp,ltyp,idest1,
     $                    idest2,llo,lup,ml,ans1,ans2,ans3,ans4,
     $                    cemab(1,kkkl)+cemab(2,kkkl),opakab(kkkl),eth,
     $                    kkkl,ptmp1,ptmp2
9004                    format(1x,8i6,' level',7(1pe10.3),i6,
     $                     7(1pe10.3))
                      endif
                     endif
                   endif
                 if ((lrtyp.eq.5).or.(lrtyp.eq.40)) then
                   call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $               nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $               ans3,ans4,idest1,idest2,idest3,idest4,
     $               abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $               opakc,opakscatt,rcdum,fline,lpriu,kdesc2,
     $               r,t,trad,tsq,xee,xh1,xh0,
     $               epi,ncn2,bremsa,bremsint,
     $               rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $               idat1,rdat1,kdat1,nptrs,np2,
     $               npar,npnxt,npfi,npfirst,
     $               nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $               npconi2,ncsvn,rates,vsav,idrates)
                   llo=idest1
                   lup=idest2
                   if (kl.eq.ilimh) lup=min(nlev,lup)
                   if ((llo.ne.0).and.(lup.ne.0).and.
     $               (llo.le.nd).and.(lup.le.nd)) then
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=ans1
                     ajisb(2,nindb)=ans2
                     cjisb(nindb)=0.
                     indb(1,nindb)=lup+ipmat
                     indb(2,nindb)=llo+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=ans2
                     ajisb(2,nindb)=ans1
                     cjisb(nindb)=0.
                     indb(1,nindb)=llo+ipmat
                     indb(2,nindb)=lup+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=-ans1
                     ajisb(2,nindb)=-ans1
                     cjisb(nindb)=0.
                     indb(1,nindb)=llo+ipmat
                     indb(2,nindb)=llo+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=-ans2
                     ajisb(2,nindb)=-ans2
                     cjisb(nindb)=0.
                     indb(1,nindb)=lup+ipmat
                     indb(2,nindb)=lup+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                     nmat=max(nmat,llo)
                     nmat=max(nmat,lup)
                     if (lpri.ge.1)
     $                 write (lun11,9004)jkk,lrtyp,ltyp,idest1,idest2,
     $                   llo,lup,ml,ans1,ans2,ans3,ans4
                     endif
                   endif
                 if ((lrtyp.eq.3).or.(lrtyp.eq.23)) then
                   call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $               nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $               ans3,ans4,idest1,idest2,idest3,idest4,
     $               abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $               opakc,opakscatt,rcdum,fline,lpriu,kdesc2,
     $               r,t,trad,tsq,xee,xh1,xh0,
     $               epi,ncn2,bremsa,bremsint,
     $               rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $               idat1,rdat1,kdat1,nptrs,np2,
     $               npar,npnxt,npfi,npfirst,
     $               nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $               npconi2,ncsvn,rates,vsav,idrates)
                   if ((idest1.gt.0).and.(idest2.gt.0).and.
     $               (idest1.le.nlev).and.(idest2.le.nlev).and.
     $               (idest1.le.nd).and.(idest2.le.nd))
     $               then
                     e1=rlev(1,idest1)
                     e2=rlev(1,idest2)
                     if ((e1/(1.e-24+e2)-1.).lt.0.01) then
                         lup=idest2
                         llo=idest1
                       else
                         lup=idest1
                         llo=idest2
                       endif
                     if (kl.eq.ilimh) lup=min(nlev,lup)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=ans2
                     ajisb(2,nindb)=ans1
                     cjisb(nindb)=0.
                     indb(1,nindb)=llo+ipmat
                     indb(2,nindb)=lup+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=ans1
                     ajisb(2,nindb)=ans2
                     cjisb(nindb)=0.
                     indb(1,nindb)=lup+ipmat
                     indb(2,nindb)=llo+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=-ans1
                     ajisb(2,nindb)=-ans1
                     cjisb(nindb)=0.
                     indb(1,nindb)=llo+ipmat
                     indb(2,nindb)=llo+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=-ans2
                     ajisb(2,nindb)=-ans2
                     cjisb(nindb)=0.
                     indb(1,nindb)=lup+ipmat
                     indb(2,nindb)=lup+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                     nmat=max(nmat,llo)
                     nmat=max(nmat,lup)
                     if (lpri.ge.1)
     $                 write (lun11,9004)jkk,lrtyp,ltyp,idest1,
     $                   idest2,llo,lup,ml,ans1,ans2,ans3,ans4
                     endif
                   endif
                 if ((lrtyp.eq.4).or.(lrtyp.eq.14).or.(lrtyp.eq.9))
     $             then
                   jkkl=nplini(ml)
                   if ((jkkl.lt.nnnl).and.(jkkl.gt.0)) then
                     tau1=tau0(1,jkkl)
                     tau2=tau0(2,jkkl)
                     ptmp1=pescl(tau1)*(1.-cfrac)
                     ptmp2=pescl(tau2)*(1.-cfrac)
     $                      +2.*pescl(tau1+tau2)*cfrac
c                     ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau2)*cfrac
c                    note that radiative rates and emissivities will
c                      have the escape probabilities in them
c                      when calculated in ucalc
                     call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $                 nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $                 ans3,ans4,idest1,idest2,idest3,idest4,
     $                 abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $                 opakc,opakscatt,rcdum,fline,lpriu,kdesc2,
     $                 r,t,trad,tsq,xee,xh1,xh0,
     $                 epi,ncn2,bremsa,bremsint,
     $                 rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $                 idat1,rdat1,kdat1,nptrs,np2,
     $                 npar,npnxt,npfi,npfirst,
     $                 nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $                 npconi2,ncsvn,rates,vsav,idrates)
                     if ((idest1.gt.0).and.(idest2.gt.0).and.
     $                 (idest1.le.nlev).and.(idest2.le.nlev).and.
     $                 (idest2.le.nd).and.(idest1.le.nd)) then
                       e1=rlev(1,idest1)
                       e2=rlev(1,idest2)
                       if (e1.lt.e2) then
                           lup=idest2
                           llo=idest1
                         else
                           lup=idest1
                           llo=idest2
                         endif
                       if (kl.eq.ilimh) lup=min(nlev,lup)
                       rcem(1,jkkl)=ans4*ptmp1/(ptmp1+ptmp2)
                       rcem(2,jkkl)=ans4*ptmp2/(ptmp1+ptmp2)
                       oplin(jkkl)=opakb1
c                      these terms represent radiative exctiation.
c                       must turn them off to attain LTE.
                       nindb=nindb+1
                       if (nindb.gt.ndb) stop 'array indexing error'
                       ajisb(1,nindb)=ans2
                       ajisb(2,nindb)=ans1
                       cjisb(nindb)=0.
                       indb(1,nindb)=lup+ipmat
                       indb(2,nindb)=llo+ipmat
                       if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                               ajisb(1,nindb),ajisb(2,nindb)
                       nindb=nindb+1
                       if (nindb.gt.ndb) stop 'array indexing error'
                       ajisb(1,nindb)=ans1
                       ajisb(2,nindb)=ans2
                       cjisb(nindb)=0.
                       indb(1,nindb)=llo+ipmat
                       indb(2,nindb)=lup+ipmat
                       if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                               ajisb(1,nindb),ajisb(2,nindb)
                       nindb=nindb+1
                       if (nindb.gt.ndb) stop 'array indexing error'
                       ajisb(1,nindb)=-ans2
                       ajisb(2,nindb)=-ans2
                       indb(1,nindb)=llo+ipmat
                       indb(2,nindb)=llo+ipmat
                       cjisb(nindb)=-ans3*xpx
                       if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                       nindb=nindb+1
                       if (nindb.gt.ndb) stop 'array indexing error'
                       ajisb(1,nindb)=-ans1
                       ajisb(2,nindb)=-ans1
                       indb(1,nindb)=lup+ipmat
                       indb(2,nindb)=lup+ipmat
                       cjisb(nindb)=ans4*xpx
                       if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                       nmat=max(nmat,llo)
                       nmat=max(nmat,lup)
                       if ((lpri.ge.1))
     $                 write (lun11,9009)jkk,lrtyp,ltyp,idest1,
     $                   idest2,llo,lup,ml,ans1,ans2,ans3,ans4,
     $                   tau1,tau2,jkkl,ptmp1,ptmp2
 9009                  format (1x,8i6,' level',6(1pe10.3),
     $                       i6,2(1pe10.3))
                       endif
                     endif
                   endif
                 if (lrtyp.eq.42) then
                   call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $               nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $               ans3,ans4,idest1,idest2,idest3,idest4,
     $               abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $               opakc,opakscatt,rcdum,fline,lpriu,kdesc2,
     $               r,t,trad,tsq,xee,xh1,xh0,
     $               epi,ncn2,bremsa,bremsint,
     $               rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $               idat1,rdat1,kdat1,nptrs,np2,
     $               npar,npnxt,npfi,npfirst,
     $               nplin,nplini,nlsvn,npcon,npconi,npilev,
     $               npilevi,npconi2,ncsvn,rates,vsav,idrates)
                   if ((idest1.gt.0).and.(idest2.gt.0).and.
     $               (idest1.le.nlev).and.(idest2.le.nlev).and.
     $               (idest2.le.nd).and.(idest1.le.nd))
     $               then
                     e1=rlev(1,idest1)
                     e2=rlev(1,idest2)
                     if (e1.lt.e2) then
                         lup=idest2
                         llo=idest1
                       else
                         lup=idest1
                         llo=idest2
                       endif
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=ans1
                     ajisb(2,nindb)=0.
                     cjisb(nindb)=0.
                     indb(1,nindb)=lup+ipmat
                     indb(2,nindb)=llo+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=0.
                     ajisb(2,nindb)=ans1
                     cjisb(nindb)=0.
                     indb(1,nindb)=llo+ipmat
                     indb(2,nindb)=lup+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=-ans1
                     ajisb(2,nindb)=-ans1
                     cjisb(nindb)=-ans3*xpx
                     indb(1,nindb)=llo+ipmat
                     indb(2,nindb)=llo+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                     nmat=max(nmat,llo)
                     nmat=max(nmat,lup)
                     if ((lpri.ge.1))
     $                 write (lun11,9004)jkk,lrtyp,ltyp,idest1,idest2,
     $                 llo,lup,ml,ans1,ans2,ans3,ans4
                     endif
                   endif
                 if (lrtyp.eq.41) then
                   call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $               nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $               ans3,ans4,idest1,idest2,idest3,idest4,
     $               abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $               opakc,opakscatt,rcdum,fline,lpriu,kdesc2,
     $               r,t,trad,tsq,xee,xh1,xh0,
     $               epi,ncn2,bremsa,bremsint,
     $               rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $               idat1,rdat1,kdat1,nptrs,np2,
     $               npar,npnxt,npfi,npfirst,
     $               nplin,nplini,nlsvn,npcon,npconi,npilev,
     $               npilevi,npconi2,ncsvn,rates,vsav,idrates)
                   if ((idest1.gt.0).and.(idest2.gt.0).and.
     $               (idest1.le.nlev).and.
     $               (idest2.le.nd).and.(idest1.le.nd))
     $               then
                     lup=idest2
                     llo=idest1
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=ans1
                     ajisb(2,nindb)=0.
                     cjisb(nindb)=0.
                     indb(1,nindb)=lup+ipmat
                     indb(2,nindb)=llo+ipmat
                    if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                               ajisb(1,nindb),ajisb(2,nindb)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=0.
                     ajisb(2,nindb)=ans1
                     cjisb(nindb)=0.
                     indb(1,nindb)=llo+ipmat
                     indb(2,nindb)=lup+ipmat
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                               ajisb(1,nindb),ajisb(2,nindb)
                     nindb=nindb+1
                     if (nindb.gt.ndb) stop 'array indexing error'
                     ajisb(1,nindb)=-ans1
                     ajisb(2,nindb)=-ans1
                     indb(1,nindb)=llo+ipmat
                     indb(2,nindb)=llo+ipmat
                     cjisb(nindb)=-ans3*xpx
                     if (lpri.gt.1)
     $                 write (lun11,*)nindb,indb(1,nindb),indb(2,nindb),
     $                        ajisb(1,nindb),ajisb(2,nindb),cjisb(nindb)
                     nmat=max(nmat,llo)
                     nmat=max(nmat,lup)
                     if ((lpri.ge.1))
     $                 write (lun11,9004)jkk,lrtyp,ltyp,idest1,idest2,
     $                   llo,lup,ml,ans1,ans2,ans3,ans4
                     endif
                   endif
                 endif
               ml=npnxt(ml)
               if (ml.ne.0) mlpar=npar(ml)
               enddo
            endif
          endif
        enddo
c
      if (lpri.ne.0)
     $   write (lun11,*)'rrrtt=',rrrtot,rrrtot2
c
c
      return
      end
      subroutine func2i(jkk,
     $       idat1,rdat1,kdat1,nptrs,
     $       npfi,npar,npnxt,nlev)
c
c     this routine counts the levels for each ion
c     author: T. Kallman
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2)
      integer npfi(ntyp,nni)
      integer np1i,np1r,np1k
c
      integer nlevmx,mltype,ml,mllz,nlev,nidt,nrdt,nkdt,
     $        jkk,ltyp,lrtyp,lcon,mlpar,lun11,mlm
c
c          now find level data
c          step thru types
           nlevmx=0
           mltype=13
           ml=npfi(mltype,jkk)
           mllz=npar(ml)
c          step thru records of this type
           mlpar=npar(ml)
           do while ((ml.ne.0).and.(mlpar.eq.mllz))
              mlm=ml-1
              call drd(ltyp,lrtyp,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $          nptrs,0,lun11)
              nlevmx=nlevmx+1
              nlev=idat1(np1i+nidt-2)
              ml=npnxt(ml)
              if (ml.ne.0) mlpar=npar(ml)
              enddo
           nlev=nlevmx
c
      return
      end
      subroutine func2l(jkk,lpri,lun11,t,xee,xpx,
     $       idat1,rdat1,kdat1,nptrs,
     $       npar,npnxt,npfi,
     $       rniss,rlev,ilev,
     $       nlpt,iltp,nlev,klev)
c
c     this routine calculates rates affecting level populations
c     author: T. Kallman
c

      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2)
      integer npfi(ntyp,nni)
      real*8 rlev(10,nd)
      integer ilev(10,nd),nlpt(nd),iltp(nd)
      integer np1i,np1r,np1k
      character(1) klev(100,nd),kblnk
c
      integer nlevmx,mltype,ml,mllz,nlev,nidt,nrdt,nkdt,lun11,
     $        lpri,ltyp,lrtyp,lcon,jkk,mm,lk,mlpar,mlm
      real*8 xee,xpx,t,bb
      real*8 rniss(nnml)
c
      data kblnk/' '/
c
      if (lpri.gt.1)
     $  write (lun11,*)'in func2l, inputs:',t,
     $          xee,xpx
c
c     now find level data
c     step thru types
      nlevmx=0
      mltype=13
      ml=npfi(mltype,jkk)
      if (lpri.gt.1) write (lun11,*)'jkk=',jkk,ml,npar(ml)
      mllz=npar(ml)
c     step thru records of this type
      mlpar=npar(ml)
      do while ((ml.ne.0).and.(mlpar.eq.mllz))
         mlm=ml-1
         call drd(ltyp,lrtyp,lcon,
     $     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $     nptrs,0,lun11)
         nlev=idat1(np1i+nidt-2)
         nlevmx=max(nlevmx,nlev)
         if ((nlev.gt.0).and.(nlev.le.nd)) then
           nlpt(nlev)=ml
           iltp(nlev)=ltyp
 9101      format (1x,'level quantities:',4i6,4(1pe12.5),3i6,8a1)
           if (lpri.gt.1) write (lun11,9101)
     $       ml,nlev,ltyp,lrtyp,(rdat1(np1r+mm-1),mm=1,4),
     $       idat1(np1i),idat1(np1i+1),
     $       idat1(np1i+2),(kdat1(np1k+mm-1),mm=1,8)
           do  lk=1,nrdt
             rlev(lk,nlev)=rdat1(np1r+lk-1)
             enddo
           do lk=1,nidt
             ilev(lk,nlev)=idat1(np1i+lk-1)
             enddo
           do lk=1,nkdt
             klev(lk,nlev)=kdat1(np1k+lk-1)
             enddo
           do lk=nkdt+1,100
             klev(lk,nlev)=kblnk
             enddo
           endif
         ml=npnxt(ml)
         if (ml.ne.0) mlpar=npar(ml)
         enddo
      nlev=nlevmx
      call levwk(rniss,bb,lpri,rlev,ilev,
     $   nlpt,iltp,nlev,klev,t,xee,xpx,lun11)
c
      return
      end
      subroutine func3(jkk,jkkl,lpri,lun11,vturbi,
     $       t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,
     $       epi,ncn2,bremsa,bremsint,tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,ncsvn,rates,vsav,idrates,
     $       rniss,rlev,ilev,
     $          nlpt,iltp,nlev,klev,
     $       xeltp,rrcor,htt,cll,cllines,clcont,rrrt,
     $       xilev,ipmat,ipmatsv,
     $       rcem,oplin,rccemis,opakc,opakscatt,
     $       cemab,cabab,opakab,fline,flinel, clsup)
c
c     this routine calculates rates affecting emission and
c        absorption
c     author: T. Kallman
c
c      note that the abundances are passed in in the array xilev
c      this array is indexed for the element as a whole
c      and for each ion the offset is the index ipmat.
c      so that for each ion the levels begin at ipmat+1 ...
c      the same does not go for rniss, the lte populations,
c      which are numbered from 1 ...
c
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line optical depths
      real*8 tau0(2,nnnl)
c     energy bins
      real*8 epi(ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn)
c     line opacities
      real*8 oplin(nnnl)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     element abundances
c     the saved rates
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
c     state variables
      real*8 r,t,xpx,delr
c     heating-cooling variables
c     input parameters
      real*8 trad
      real*8 vturbi,xee
      integer ncn2,lpri,lun11,np2
      integer nlsvn,ncsvn
      real*8 rlev(10,nd)
      integer ilev(10,nd),nlpt(nd),iltp(nd)
      character(1) klev(100,nd)
      character(49) kdesc2
      real*8 fline(2,nnnl),flinel(ncn)
c     level populations
      real*8 xilev(nd)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
      real*8 tauc(2,nnml)
      character(1) kblnk
      real*8 tsq,ans1,ans2,xh1,xh0,cfrac
      real*8 abund1,abund2,ptmp1,ptmp2,ans3,ans4,opakb1,
     $     xeltp,rrcor,cllines,clcont,htt,cll
      integer idest1,idest2,idest3,idest4
      real*8 rniss(nd)
      real*8 abundtot,rrrt
      real*8 tau1,tau2,e1,e2,pescl,pescv,
     $     cemtmp1,cemtmp2,czzz,elin,ener,htsum,eth,opakbb,
     $     rcemm,rcsum,ergsev
      integer nlev,nlevmx,mltype,ml,mllz,mlrdesc,lpriu,
     $     llo,lup,jkk,ipmat,ltyp,
     $     lrtyp,lcon,nrdt,nidt,nkdt,lk,kkkl,jkkl,ipmatsv,
     $     lprisv,ml3,mm,nb1,nbinc,mlpar,mlm,lfpi
      integer np1i,np1r,np1k
      real*8 clsup
c
c
      data kblnk/' '/
c
      ergsev=1.602197e-12
      lprisv=lpri
c
      if (lpri.gt.0)
     $  write (lun11,*)'in func3, inputs:',t,xee,xpx,delr,ipmat,ipmatsv
c

c
c     lfpi mode:  opacities only
      lfpi=3
c
      abundtot=0.
      rrrt=0.
      tsq=sqrt(t)
c
      lprisv=lpri
c
      htt=0.
      cll=0.
      if (lpri.ne.0) then
        write (lun11,*)'level populations:'
        do mm=1,nlev
          write (lun11,9022)mm,(klev(ml,mm),ml=1,20),
     $      rlev(1,mm),rlev(2,mm),
     $      xilev(mm+ipmat),rniss(mm)
 9022     format (i4,20a1,4(1pe10.3))
          enddo
        endif
c
c
c     now do other  rates
      lpriu=lpri
      mltype=9
      mlrdesc=mltype
      ml=npfi(mltype,jkk)
      if (ml.ne.0) then
        mllz=npar(ml)
        mlpar=npar(ml)
        do while ((ml.ne.0).and.(mlpar.eq.mllz))
          mlm=ml-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          idest1=idat1(np1i)
          idest2=idat1(np1i+1)
          jkkl=nplini(ml)
          if ((rdat1(np1r).gt.0.01).and.(jkkl.ne.0)
     $      .and.(idest1.gt.0).and.(idest2.gt.0).and.
     $      (idest1.lt.nlev).and.(idest2.lt.nlev)) then
            e1=rlev(1,idest1)
            e2=rlev(1,idest2)
            eth=abs(e2-e1)
            if (e1.lt.e2) then
                lup=idest2+ipmat
                llo=idest1+ipmat
              else
                lup=idest1+ipmat
                llo=idest2+ipmat
              endif
            abund1=xilev(llo)*xpx*xeltp
            abund2=xilev(lup)*xpx*xeltp
            tau1=tau0(1,jkkl)
            tau2=tau0(2,jkkl)
            ptmp1=
     $        pescl(tau1)*(1.-cfrac)
            ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau1+tau2)*cfrac
c            ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau2)*cfrac
            call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $        ans3,ans4,idest1,idest2,idest3,idest4,
     $        abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $        opakc,opakscatt,rccemis,fline,lpriu,kdesc2,
     $        r,t,trad,tsq,xee,xh1,xh0,
     $        epi,ncn2,bremsa,bremsint,
     $        rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $        idat1,rdat1,kdat1,nptrs,np2,
     $        npar,npnxt,npfi,npfirst,
     $        nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $        npconi2,ncsvn,rates,vsav,idrates)
c
            cemtmp1=abund2*ptmp1*ans4
            cemtmp2=abund2*ptmp2*ans4
            rcsum=ans4*abund2
            cll=cll+rcsum
            clcont=clcont+rcsum
            if ((lpri.ge.1))
     $        write (lun11,9002)jkk,lrtyp,ltyp,idest1,idest2,
     $        llo,lup,ml,ans1,ans2,ans3,ans4,
     $        cemtmp1+cemtmp2,opakb1,eth,
     $        jkkl,cll,htt
            endif
          ml=npnxt(ml)
          mlpar=0
          if (ml.ne.0) mlpar=npar(ml)
          enddo
        endif
c

c     now do other  rates
      mltype=7
      mlrdesc=mltype
      ml=npfi(mltype,jkk)
      mllz=0
      if (ml.ne.0) mllz=npar(ml)
      mlpar=mllz
      do while ((ml.ne.0).and.(mlpar.eq.mllz))
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        idest1=idat1(np1i+nidt-2)
        idest2=nlev+idat1(np1i-1+nidt-3)-1
        kkkl=npconi2(ml)
        if ((kkkl.ne.0).and.(kkkl.le.ndat2)
     $     .and.(idest1.gt.0)) then
          llo=idest1+ipmat
          lup=idest2+ipmat
          eth=rlev(4,idest1)-rlev(1,idest1)
          abund1=xilev(llo)*xeltp
          abund2=xilev(lup)*xeltp
          if (lup.gt.ipmatsv) then
            lup=min(lup,ipmatsv)
            abund2=0.
            endif
          nb1=nbinc(eth,epi,ncn2)
          if ((lup.le.ipmatsv).and.
     $      ((cabab(kkkl).gt.1.e-34).or.
     $        ((cemab(1,kkkl)+cemab(2,kkkl)).gt.1.e-34))) then
            tau1=tauc(1,kkkl)
            tau2=tauc(2,kkkl)
            ptmp1=
     $        pescv(tau1)*(1.-cfrac)
            ptmp2=pescv(tau2)*(1.-cfrac)+2.*pescv(tau1+tau2)*cfrac
            call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $        ans3,ans4,idest1,idest2,idest3,idest4,
     $        abund1,abund2,ptmp1,ptmp2,xpx,opakab(kkkl),
     $        opakc,opakscatt,rccemis,fline,lpriu,kdesc2,
     $        r,t,trad,tsq,xee,xh1,xh0,
     $        epi,ncn2,bremsa,bremsint,
     $        rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $        idat1,rdat1,kdat1,nptrs,np2,
     $        npar,npnxt,npfi,npfirst,
     $        nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $        npconi2,ncsvn,rates,vsav,idrates)
            if (lpri.ge.1)
     $        write (lun11,9002)jkk,lrtyp,ltyp,
     $          idest1,idest2,llo,lup,ml,ans1,ans2,
     $          ans3,ans4,cemab(1,kkkl)+cemab(2,kkkl),
     $          opakab(kkkl),eth,
     $          kkkl,cll,htt
 9002         format (1x,8i6,' h-c ',
     $          7(1pe10.3),i6,2(1pe10.3),4(1pe10.3))
            endif
          endif
        ml=npnxt(ml)
        mlpar=0
        if (ml.ne.0) mlpar=npar(ml)
        enddo
c
      mltype=42
      mlrdesc=mltype
      ml=npfi(mltype,jkk)
      if (ml.ne.0) then
        mllz=npar(ml)
        mlpar=npar(ml)
        do while ((ml.ne.0).and.(mlpar.eq.mllz))
          mlm=ml-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          idest1=idat1(np1i+nidt-2)
          idest2=idat1(np1i+nidt-3)
          if (idest1.gt.0) then
            llo=idest1+ipmat
            lup=idest2+ipmat
            eth=rlev(4,idest1)-rlev(1,idest1)
            abund1=xilev(llo)*xeltp
            abund2=xilev(lup)*xeltp
            if (lup.gt.ipmatsv) then
              abund2=0.
              lup=ipmatsv
              endif
            if ((lup.le.ipmatsv).and.
     $      (xilev(llo)/(1.e-36+xilev(1+ipmat)).gt.1.e-24)) then
              call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $          ans3,ans4,idest1,idest2,idest3,idest4,
     $          abund1,abund2,ptmp1,ptmp2,xpx,opakbb,
     $          opakc,opakscatt,rccemis,fline,lpriu,kdesc2,
     $          r,t,trad,tsq,xee,xh1,xh0,
     $          epi,ncn2,bremsa,bremsint,
     $          rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $          idat1,rdat1,kdat1,nptrs,np2,
     $          npar,npnxt,npfi,npfirst,
     $          nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $          npconi2,ncsvn,rates,vsav,idrates)
              htsum=ans3*xpx*abund1
              htt=htt+htsum
              czzz=0.
              if (lpri.ge.1)
     $        write (lun11,9002)jkk,lrtyp,ltyp,
     $          idest1,idest2,llo,lup,ml,ans1,ans2,
     $          ans3,ans4,czzz,opakbb,eth,kkkl,cll,htt
              endif
            endif
          ml=npnxt(ml)
          if (ml.ne.0) mlpar=npar(ml)
          enddo
        endif
c
      mltype=1
      mlrdesc=mltype
      ml=npfi(mltype,jkk)
      mllz=0
      if (ml.ne.0) mllz=npar(ml)
      mlpar=mllz
      do while ((ml.ne.0).and.(mlpar.eq.mllz))
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        if (.not.((mlrdesc.eq.1).and.((ltyp.eq.93).or.(ltyp.eq.53)))
     $      .and.(.not.((mlrdesc.eq.7).and.(idat1(np1i+nidt-2).ne.1))))
     $      then
          if (nidt.gt.3) then
            idest1=idat1(np1i+nidt-2)
            idest2=nlev+idat1(np1i+nidt-4)-1
            kkkl=npconi2(ml)
c            if (lpri.ne.0) write (lun11,*)'kkkl=',kkkl,idest1,ltyp
            if ((kkkl.ne.0).and.(kkkl.le.ndat2)
     $        .and.(idest1.gt.0)) then
              llo=idest1+ipmat
              lup=idest2+ipmat
c              eth=rlev(4,idest1)-rlev(1,idest1)
              eth=rdat1(np1r)
              abund1=xilev(llo)*xeltp
              abund2=xilev(lup)*xeltp
              if (lup.gt.ipmatsv) then
                abund2=0.
                lup=ipmatsv
                endif
c              if (lpri.ne.0) write (lun11,*)lup,ipmatsv,llo,xilev(llo),
c     $               xilev(1+ipmat)
              if ((lup.le.ipmatsv).and.
     $          (xilev(llo)/(1.e-36+xilev(1+ipmat)).gt.1.e-24)) then
                tau1=tauc(1,kkkl)
                tau2=tauc(2,kkkl)
                ptmp1=
     $          pescv(tau1)*(1.-cfrac)
                ptmp2=pescv(tau2)*(1.-cfrac)+2.*pescv(tau1+tau2)*cfrac
                call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $            nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $            ans3,ans4,idest1,idest2,idest3,idest4,
     $            abund1,abund2,ptmp1,ptmp2,xpx,opakab(kkkl),
     $            opakc,opakscatt,rccemis,fline,lpriu,kdesc2,
     $            r,t,trad,tsq,xee,xh1,xh0,
     $            epi,ncn2,bremsa,bremsint,
     $            rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $            idat1,rdat1,kdat1,nptrs,np2,
     $            npar,npnxt,npfi,npfirst,
     $            nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $            npconi2,ncsvn,rates,vsav,idrates)
                htsum=ans3*xpx*abund1
                rrrt=rrrt+xilev(lup)*(ptmp1+ptmp2)*ans2
                abundtot=abundtot+xilev(lup)
                htt=htt+htsum
                if (lpri.ge.1)
     $            write (lun11,9002)jkk,lrtyp,ltyp,
     $            idest1,idest2,llo,lup,ml,ans1,ans2,
     $            ans3,ans4,cemab(1,kkkl)+cemab(2,kkkl),
     $            opakab(kkkl),eth,
     $            kkkl,cll,htt
                endif
              endif
            endif
          endif
        ml=npnxt(ml)
        mlpar=0
        if (ml.ne.0) mlpar=npar(ml)
        enddo
c
      mltype=4
      mlrdesc=mltype
      ml=npfi(mltype,jkk)
      if (ml.ne.0) then
        mllz=npar(ml)
        mlpar=npar(ml)
        do while ((ml.ne.0).and.(mlpar.eq.mllz))
             mlm=ml-1
             call drd(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $            nptrs,0,lun11)
             idest1=idat1(np1i)
             idest2=idat1(np1i+1)
             jkkl=nplini(ml)
             if ((rdat1(np1r).gt.0.01).and.(jkkl.ne.0)
     $          .and.(idest1.gt.0).and.(idest2.gt.0).and.
     $          (idest1.lt.nlev).and.(idest2.lt.nlev)) then
               e1=rlev(1,idest1)
               e2=rlev(1,idest2)
               if (e1.lt.e2) then
                   lup=idest2+ipmat
                   llo=idest1+ipmat
                 else
                   lup=idest1+ipmat
                   llo=idest2+ipmat
                 endif
c               if (lpri.ne.0) write (lun11,*)idest1,idest2
               if ((xilev(llo)/(1.e-36+xilev(1+ipmat)).gt.1.e-24).or.
     $             (xilev(lup)/(1.e-36+xilev(1+ipmat)).gt.1.e-24)) then
                 abund1=xilev(llo)*xpx*xeltp
                 abund2=xilev(lup)*xpx*xeltp
                 ml3=nplin(jkkl)
                 tau1=tau0(1,jkkl)
                 tau2=tau0(2,jkkl)
                 ptmp1=
     $                pescl(tau1)*(1.-cfrac)
                 ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau1+tau2)*cfrac
c                 ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau2)*cfrac
                 lpriu=lpri
c                we need to call ucalc again because rcem
c                already has the abundance in from func3p
                 call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $             nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $             ans3,ans4,idest1,idest2,idest3,idest4,
     $             abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $             opakc,opakscatt,rccemis,fline,lpriu,kdesc2,
     $             r,t,trad,tsq,xee,xh1,xh0,
     $             epi,ncn2,bremsa,bremsint,
     $             rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $             idat1,rdat1,kdat1,nptrs,np2,
     $             npar,npnxt,npfi,npfirst,
     $             nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $             npconi2,ncsvn,rates,vsav,idrates)
                 rcem(1,jkkl)=abund2*ans4*ptmp1/(ptmp1+ptmp2)
                 rcem(2,jkkl)=abund2*ans4*ptmp2/(ptmp1+ptmp2)
                 ml3=nplin(jkkl)
                 oplin(jkkl)=opakb1*abund1
                 if (ml3.ne.0) then
                   mlm=ml3-1
                   call drd(ltyp,lrtyp,lcon,
     $               nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $               nptrs,0,lun11)
                   elin=abs(rdat1(np1r))
                   ener=12398.41/elin
                   nb1=nbinc(ener,epi,ncn2)
                   nb1=max(2,min(ncn2-1,nb1))
c                   fline(1,jkkl)=(ans1*xilev(lup)-ans2*xilev(llo))
c     $                *xpx*xeltp*ener*ergsev*ptmp1
c                   fline(2,jkkl)=(ans1*xilev(lup)-ans2*xilev(llo))
c     $                *xpx*xeltp*ener*ergsev*ptmp2
c                   flinel(nb1)=flinel(nb1)+(fline(1,jkkl)+fline(2,jkkl))
c     $               *2./(epi(nb1+1)-epi(nb1-1))/ergsev
                   cll=cll+rcem(1,jkkl)+rcem(2,jkkl)
                   htt=htt+abund1*ans3
                   cllines=cllines+rcem(1,jkkl)+rcem(2,jkkl)
                   if ((lpri.ge.1))
     $               write (lun11,9002)jkk,lrtyp,
     $                 ltyp,idest1,idest2,
     $                 llo,lup,ml,ans1,ans2,ans3,ans4,
     $                 rcem(1,jkkl)+rcem(2,jkkl),oplin(jkkl),
     $                 rdat1(np1r),jkkl,cll,htt
     $                  ,ptmp1,ptmp2
c     $                 ,rnrb,rnrb*xpx*xeltp*ener*ergsev
                   endif
                 endif
               endif
             ml=npnxt(ml)
             if (ml.ne.0) mlpar=npar(ml)
             enddo
        endif
c
           mltype=14
           mlrdesc=mltype
           ml=npfi(mltype,jkk)
           mllz=0
           if (ml.ne.0)  mllz=npar(ml)
           mlpar=mllz
           do while ((ml.ne.0).and.(mlpar.eq.mllz))
             mlm=ml-1
             call drd(ltyp,lrtyp,lcon,
     $         nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $         nptrs,0,lun11)
             idest1=idat1(np1i-1+nidt-3)
             idest2=idat1(np1i+nidt-3)
             if ((idest1.gt.0).and.(idest2.gt.0).and.
     $          (idest1.lt.nlev).and.(idest2.lt.nlev)) then
               e1=rlev(1,idest1)
               e2=rlev(1,idest2)
               eth=abs(e2-e1)
               if (e1.lt.e2) then
                   lup=idest2+ipmat
                   llo=idest1+ipmat
                 else
                   lup=idest1+ipmat
                   llo=idest2+ipmat
                 endif
               if ((xilev(llo)/(1.e-36+xilev(1+ipmat)).gt.1.e-24).or.
     $             (xilev(lup)/(1.e-36+xilev(1+ipmat)).gt.1.e-24)) then
                 abund1=xilev(llo)*xpx*xeltp
                 abund2=xilev(lup)*xpx*xeltp
                 tau1=tau0(1,jkkl)
                 tau2=tau0(2,jkkl)
                 ptmp1=
     $                pescl(tau1)*(1.-cfrac)
                 ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau1+tau2)*cfrac
                 lpriu=lpri
                 call ucalc(ltyp,lrtyp,ml,lcon,jkk,vturbi,
     $             nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $             ans3,ans4,idest1,idest2,idest3,idest4,
     $             abund1,abund2,ptmp1,ptmp2,xpx,opakb1,
     $             opakc,opakscatt,rccemis,fline,lpriu,kdesc2,
     $             r,t,trad,tsq,xee,xh1,xh0,
     $             epi,ncn2,bremsa,bremsint,
     $             rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfpi,lun11,
     $             idat1,rdat1,kdat1,nptrs,np2,
     $             npar,npnxt,npfi,npfirst,
     $             nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $             npconi2,ncsvn,rates,vsav,idrates)
                 rcemm=abund2*ans4
                 rccemis(2,3)=rccemis(2,3)+
     $                  rcemm/(epi(4)-epi(3)+1.e-24)/ergsev/12.56
c                if (lpri.ne.0) write (lun11,*)jkkl,tau0(1,jkkl),
                 cll=cll+rcemm
                 clcont=clcont+rcemm
                 clsup=clsup+rcemm
                 if ((lpri.ge.1))
     $              write (lun11,9002)jkk,lrtyp,ltyp,
     $               idest1,idest2,llo,lup,ml,ans1,ans2,
     $               ans3,ans4,rcemm,opakb1,eth,
     $               kkkl,cll,htt
                 endif
               endif
             ml=npnxt(ml)
             mlpar=0
             if (ml.ne.0) mlpar=npar(ml)
             enddo
c
      rrrt=rrrt/max(1.e-24,abundtot)
      lpri=lprisv
c
c
      return
      end
      subroutine func3p(jkk,jkkl,lpri,lun11,vturbi,
     $       t,trad,r,xee,xpx,xh1,xh0,cfrac,
     $       epi,ncn2,bremsa,bremsint,tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,ncsvn,rates,vsav,idrates,
     $       rniss,rlev,ilev,
     $          nlpt,iltp,nlev,klev,
     $       xeltp,rrcor,htt,cll,cllines,clcont,rrrt,
     $       xilev,ipmat,ipmatsv,
     $       rcem,oplin,rccemis,opakc,opakscatt,
     $       cemab,cabab,opakab,fline,flinel)
c
c     this routine prints level populations
c     author: T. Kallman
c
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line optical depths
      real*8 tau0(2,nnnl)
c     energy bins
      real*8 epi(ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn)
c     line opacities
      real*8 oplin(nnnl)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     element abundances
c     the saved rates
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
c     state variables
      real*8 r,t,xpx
c     heating-cooling variables
c     input parameters
      real*8 trad
      real*8 vturbi,xee
      integer ncn2,lpri,lun11,np2
      integer nlsvn,ncsvn
      real*8 rlev(10,nd)
      integer ilev(10,nd),nlpt(nd),iltp(nd)
      character(1) klev(100,nd)
      character(49) kdesc2
      real*8 fline(2,nnnl),flinel(ncn)
c     level populations
      real*8 xilev(nd)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
      real*8 tauc(2,nnml)
      character(1) kblnk
      real*8 tsq,ans1,ans2,xh1,xh0,cfrac
      real*8 abund1,abund2,ptmp1,ptmp2,ans3,ans4,opakb1,
     $     xeltp,rrcor,cllines,clcont,htt,cll
      integer idest1,idest2,idest3,idest4
      real*8 rniss(nd)
      real*8 abundtot,rrrt
      real*8 tau1,tau2,e1,e2,pescl,pescv,
     $     cemtmp1,cemtmp2,czzz,elin,ener,htsum,eth,opakbb,
     $     rcemm,rcsum,ergsev
      integer nlev,nlevmx,mltype,ml,mllz,mlrdesc,lpriu,
     $     llo,lup,jkk,ipmat,ltyp,
     $     lrtyp,lcon,nrdt,nidt,nkdt,lk,kkkl,jkkl,ipmatsv,
     $     lprisv,ml3,mm,nb1,nbinc,mlpar,mlm
      integer np1i,np1r,np1k
c
c
      data kblnk/' '/
c
      ergsev=1.602197e-12
c
      abundtot=0.
      rrrt=0.
      tsq=sqrt(t)
c
      lprisv=lpri
c
      if (lpri.gt.0)
     $  write (lun11,*)'in func3p, inputs:',t,
     $         xee,xpx
c
c
      htt=0.
      cll=0.
      if (lpri.ne.0) then
        write (lun11,*)'level populations:'
        do mm=1,nlev
          write (lun11,9022)mm,(klev(ml,mm),ml=1,20),
     $      rlev(1,mm),rlev(2,mm),
     $      xilev(mm+ipmat),rniss(mm)
 9022     format (i4,20a1,4(1pe10.3))
          enddo
        endif
c
c     now do other  rates
      mltype=7
      mlrdesc=mltype
      ml=npfi(mltype,jkk)
      mllz=0
      if (ml.ne.0) mllz=npar(ml)
      mlpar=mllz
      do while ((ml.ne.0).and.(mlpar.eq.mllz))
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        idest1=idat1(np1i+nidt-2)
        idest2=nlev+idat1(np1i-1+nidt-3)-1
        kkkl=npconi2(ml)
        if ((kkkl.ne.0).and.(kkkl.le.ndat2)
     $     .and.(idest1.gt.0)) then
          llo=idest1+ipmat
          lup=idest2+ipmat
          eth=rlev(4,idest1)-rlev(1,idest1)
          abund1=xilev(llo)*xeltp
          abund2=xilev(lup)*xeltp
          cabab(kkkl)=cabab(kkkl)*abund1*xpx
          cemab(1,kkkl)=cemab(1,kkkl)*abund2*xpx
          cemab(2,kkkl)=cemab(2,kkkl)*abund2*xpx
          if ((lpri.ge.1))
     $        write (lun11,9002)jkk,lrtyp,ltyp,idest1,idest2,
     $        llo,lup,ml,
     $        cemab(1,kkkl)+cemab(2,kkkl),opakab(kkkl),eth,
     $        kkkl,cll,htt
 9002         format (1x,8i6,' h-cp',
     $          40x,3(1pe10.3),i6,2(1pe10.3),4(1pe10.3))
          endif
        ml=npnxt(ml)
        mlpar=0
        if (ml.ne.0) mlpar=npar(ml)
        enddo
c
      mltype=4
      mlrdesc=mltype
      ml=npfi(mltype,jkk)
      if (ml.ne.0) then
        mllz=npar(ml)
        mlpar=npar(ml)
        do while ((ml.ne.0).and.(mlpar.eq.mllz))
             mlm=ml-1
             call drd(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $            nptrs,0,lun11)
             idest1=idat1(np1i)
             idest2=idat1(np1i+1)
             jkkl=nplini(ml)
             if ((rdat1(np1r).gt.0.01).and.(jkkl.ne.0)
     $          .and.(idest1.gt.0).and.(idest2.gt.0).and.
     $          (idest1.lt.nlev).and.(idest2.lt.nlev)) then
               e1=rlev(1,idest1)
               e2=rlev(1,idest2)
               if (e1.lt.e2) then
                   lup=idest2+ipmat
                   llo=idest1+ipmat
                 else
                   lup=idest1+ipmat
                   llo=idest2+ipmat
                 endif
c               if (lpri.ne.0) write (lun11,*)idest1,idest2
               if ((xilev(llo)/(1.e-36+xilev(1+ipmat)).gt.1.e-24).or.
     $             (xilev(lup)/(1.e-36+xilev(1+ipmat)).gt.1.e-24)) then
                 abund1=xilev(llo)*xpx*xeltp
                 abund2=xilev(lup)*xpx*xeltp
                 tau1=tau0(1,jkkl)
                 tau2=tau0(2,jkkl)
                 ptmp1=
     $                pescl(tau1)*(1.-cfrac)
                 ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau1+tau2)*cfrac
c                 ptmp2=pescl(tau2)*(1.-cfrac)+2.*pescl(tau2)*cfrac
                 lpriu=lpri
                 rcem(1,jkkl)=rcem(1,jkkl)*abund2
                 rcem(2,jkkl)=rcem(2,jkkl)*abund2
                 oplin(jkkl)=oplin(jkkl)*abund1
                 if ((lpri.ge.1))
     $               write (lun11,9002)jkk,lrtyp,
     $                 ltyp,idest1,idest2,
     $                 llo,lup,ml,
     $                 rcem(1,jkkl)+rcem(2,jkkl),oplin(jkkl),
     $                 rdat1(np1r),jkkl,cll,htt
     $                  ,ptmp1,ptmp2
                 endif
               endif
             ml=npnxt(ml)
             if (ml.ne.0) mlpar=npar(ml)
             enddo
        endif
c
      lpri=lprisv
c
c
      return
      end
      subroutine funcsyn(lpri,lun11,vturbi,critf,
     $       t,trad,r,delr,xee,xpx,abel,cfrac,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       decomp,ecomp,sxcomp,
     $       tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,nlevs,ncsvn,rates,vsav,idrates,
     $       xiin,rrrts,pirts,htt,cll,httot,cltot,hmctot,elcter,
     $         cllines,clcont,htcomp,clcomp,clbrems,
     $       xilevt,bilevt,rnist,nmat,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,fline,flinel, clsup)

c
c     calculates opacities and emissivities and does transfer
c     level populations and integrates continuum emissivities
c     and opacities are assumed as input
c
c     author: T. Kallman
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
      real*8 zrems(4,ncn),zremso(4,ncn)
      real*8 elum(3,nnnl),elumo(3,nnnl)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
c     line optical depths
      real*8 tau0(2,nnnl)
c     energy bins
      real*8 epi(ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn),brcems(ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
      real*8 fline(2,nnnl),flinel(ncn)
c     level populations
      real*8 xilevt(nnml),bilevt(nnml),rnist(nnml)
c     ion abundances
      real*8 xiin(nni)
      real*8 rrrts(nni),pirts(nni)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
      real*8 elumab(2,nnml),elumabo(2,nnml)
      real*8 tauc(2,nnml)
      real*8 rlev(10,nd)
      integer ilev(10,nd),nlpt(nd),iltp(nd)
      real*8 htt(nni),cll(nni)
      integer nlevs(nni)
c     element abundances
      real*8 abel(nl)
c     the saved rates
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
      real*8 decomp(ncomp,ncomp),ecomp(ncomp),sxcomp(ncomp)
c
      character(1) klev(100,nd)
c
      real*8 rrrt(31),pirt(31),xin(31),xitmp(31)
      real*8 ajisb(2,ndb),cjisb(ndb)
      integer indb(2,ndb)
      real*8 xilev(nd),rniss(nnml)
      real*8 rrcor(nni),pirtt(31)
      real*8 bmat(nd),bmatl(nd)
      real*8 rnisl(nd)
      real*8 x(nd)
      integer ipsv(31),nsup(nd)
      integer lpri,lun11,lcdd,ncn2,np2,ncsvn,nmat,nlsvn
      real*8 httotd,cltotd,hmctotd,elcterd,
     $     htcompd,clcompd,clbremsd
      real*8 vturbi,critf,t,trad,r,delr,xee,xpx,cfrac,p,
     $     hmctot,elcter,cllines,clcont,htcomp,clcomp,clbrems
      real*8 xh1,xh0,httot,cltot
      real*8 rnisum,crith,cltmp,cmp1,cmp2,
     $     cltot2,enelec,httot2,httmp,pirtsum,rniss2,rrrtt,
     $     rtdm,tt1,tt2,xintp,xeltp,ximax,xilast,xintp2,
     $     xisum,xipp,cl,ht
      integer nlev,nindb,lprisv,
     $     jkk,ipmat,ltyp,ldir,llp,imax,ilimh,
     $     lrtyp,lcon,nrdt,nidt,nkdt,ll,jkkl,ipmatsv,
     $     iliml,jk,kl1,mm,kl,kl2,klion,klel,klp,llm,lp,lm,
     $     lprim,lprif,lpril,lpritp,lsum,lsumt,ndtmp,mlel,
     $     ml1,mmt,mllel,mlion,mleltp,mmtmp,nit,nit2,nit3,nitmx,
     $     nitmx2,nlevm,nnz,nnzp,nsp,mlm,np1i,np1r,np1k
      real*8 clsup
c
      lprisv=lpri
      lprif=lpri
      lpritp=0
      if (lprif.ne.0)
     $  write (lun11,*)'in funcsyn, inputs:',t,
     $         xee,xpx,lcdd,p,abel(1),delr
       if (lcdd.ne.1)
     $   xpx = p/1.38e-12/max(t,1.e-24)
c
      xh0=xpx*xiin(1)*abel(1)
      xh1=xpx*(1.-xiin(1))*abel(1)
c
c      zero emissivitiesd and opacities
c      note that here variables
c      on the level grid  are
c      already calculated in func3p
       do ll=1,nnnl
         fline(1,ll)=0.
         fline(2,ll)=0.
         rcem(1,ll)=0.
         rcem(2,ll)=0.
         oplin(ll)=0.
         enddo
       do ll=1,ncn2
         rccemis(1,ll)=0.
         rccemis(2,ll)=0.
         opakc(ll)=0.
         opakscatt(ll)=0.
         enddo
       do ll=1,ncn2
         flinel(ll)=0.
         enddo
c
c
c      now calculate.  first step thru elements
       jkk=0
       jkkl=0
       klel=11
       mlel=npfirst(klel)
       jk=0
       xilast=0.
       do while (mlel.ne.0)
         mlm=mlel-1
         call drd(ltyp,lrtyp,lcon,
     $     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $     nptrs,0,lun11)
         mllel=0
         if (nidt.gt.0) then
           if (lprif.ne.0)
     $         write (lun11,9339)(kdat1(np1k-1+mm),mm=1,nkdt)
 9339      format (1x, ' element:',12a1)
           mllel=idat1(np1i)
           jk=mllel
           nnz=idat1(np1i)
           nnzp=nnz+1
           xeltp=0.
           if (jk.gt.0) xeltp=abel(jk)
           if (xeltp.gt.1.e-24) then
c
c            find ion indeces
             if (lprif.ne.0) write (lun11,*)' finding ion indeces'
             klion=12
             mlion=npfirst(klion)
             jkk=0
             kl=0
             do while ((mlion.ne.0).and.(kl.lt.nnz))
               jkk=jkk+1
               mlm=mlion-1
               call drd(ltyp,lrtyp,lcon,
     $           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $           nptrs,0,lun11)
               mleltp=npar(mlion)
               if (mleltp.eq.mlel) then
                 kl=kl+1
                 endif
               mlion=npnxt(mlion)
               enddo
c
c            unpack ion fractions
             xisum=0.
             ldir=+1
             do mm=1,kl
                 kl1=mm
                 xin(kl1)=xiin(jkk-kl+kl1)
                 xisum=xisum+xin(kl1)
                 enddo
             xipp=1.-xisum
             klp=kl+1
             xin(klp)=xipp
c
c            find iliml, ilimh
             iliml=0
             ilimh=0
             imax=0
             ximax=0.
             do mm=1,klp
               if (xin(mm).gt.ximax) then
                 ximax=xin(mm)
                 imax=mm
                 endif
               enddo
             imax=max(min(nnz,imax),1)
             llp=imax
             llm=imax
             iliml=imax
             ilimh=imax
             lp=0
             lm=0
             if (imax.ne.klp) then
                 lsumt=nlevs(jkk-kl+imax)
               else
                 lsumt=0
               endif
             ndtmp=nd
             mmt=0
             do while ((lsumt.lt.ndtmp).and.(ldir.ne.0))
               mmt=mmt+1
               lsum=lsumt
               iliml=min(iliml,llm)
               ilimh=max(ilimh,llp)
               if ((llp.ge.klp).or.(xin(llp).lt.critf)) then
                 ldir=-1
                 lp=1
                 endif
               if ((llm.le.1).or.(xin(llm).lt.critf)) then
                 ldir=+1
                 lm=1
                 endif
               if ((lm.ne.1).and.(lp.ne.1)) then
                 if (xin(llp+1).gt.xin(llm-1)) then
                     ldir=+1
                   else
                     ldir=-1
                   endif
                 endif
               if ((lp.eq.1).and.(lm.eq.1)) ldir=0
               if (ldir.eq.+1) then
                   llp=llp+1
                   if (llp.ne.klp) then
                       lsumt=lsum+nlevs(jkk-kl+llp)
                     else
                       lsumt=lsum
                     endif
                   endif
               if (ldir.eq.-1) then
                   llm=llm-1
                   lsumt=lsum+nlevs(jkk-kl+llm)
                   endif
               ilimh=max(ilimh-1,iliml+1)
               enddo
             if (lpri.ne.0) write (lun11,*)'iliml,ilimh:',iliml,ilimh
c
c            step thru ions
             if (lprif.ne.0) write (lun11,*)' third pass',ipmat
             mlion=npfirst(klion)
             ipmat=ipmat+1
             ipmat=0
             jkk=jkk-nnz
             kl=0
             do while ((mlion.ne.0).and.(kl.lt.nnz))
               ltyp=klion
               mlm=mlion-1
               call drd(ltyp,lrtyp,lcon,
     $           nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $           nptrs,0,lun11)
               mleltp=npar(mlion)
               if (mleltp.eq.mlel) then
                 kl=kl+1
                 jkk=jkk+1
                 call func2i(jkk,
     $             idat1,rdat1,kdat1,nptrs,
     $             npfi,npar,npnxt,nlev)
                 if ((kl.ge.iliml).and.(kl.le.ilimh)) then
                   if (lprif.ne.0)
     $              write (lun11,9338)(kdat1(np1k-1+mm),mm=1,nkdt)
9338               format (1x, ' ion:',8a1)
                   call func2l(jkk,lpri,lun11,t,xee,xpx,
     $                  idat1,rdat1,kdat1,nptrs,
     $                  npar,npnxt,npfi,
     $                  rniss,rlev,ilev,
     $                  nlpt,iltp,nlev,klev)
                   nlevm=nlev-1
c                  retrieve saved abundances
                   do mm=1,nlevm
                     bmatl(mm+ipmat)=xilev(mm)
                     mmtmp=npilev(mm,jkk)
                     if (mmtmp.gt.0) then
                       if (mmtmp.gt.nnml) stop 'mmtmp error'
                       bmatl(mm+ipmat)=xilevt(mmtmp)
                       if (lpri.gt.1)
     $                  write (lun11,*)mm,mmtmp,xilevt(mmtmp)
                       endif
c
c                    nb this code makes H and He fully ionized
c                     if (jkk.le.3) then
c                       bmatl(mm+ipmat)=0.
c                       endif
c
                     enddo
                   mmtmp=npilev(nlev,jkk)
                   if (mmtmp.gt.nnml) stop 'mmtmp error'
                   bmatl(ipmat+nlev)=xilevt(mmtmp)
c
c                  nb this code makes H and He fully ionized
c                   if ((jkk.eq.1).or.(jkk.eq.3))
c     $               bmatl(nlev+ipmat)=1.
c
                   rrcor(jkk)=1.
                   ipmatsv=ipmat+nlev
                   call func3(jkk,jkkl,lpri,lun11,vturbi,
     $                 t,trad,r,delr,xee,xpx,xh1,xh0,cfrac,
     $                 epi,ncn2,bremsa,bremsint,tau0,tauc,
     $                 idat1,rdat1,kdat1,nptrs,np2,
     $                 npar,npnxt,npfi,npfirst,
     $                 nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $                 npconi2,ncsvn,rates,vsav,idrates,
     $                 rniss,rlev,ilev,
     $                    nlpt,iltp,nlev,klev,
     $                 xeltp,rrcor(jkk),httmp,cltmp,cllines,clcont,rtdm,
     $                 bmatl,ipmat,ipmatsv,
     $                 rcem,oplin,rccemis,opakc,opakscatt,
     $                 cemab,cabab,opakab,fline,flinel, clsup)
c
                   ipmat=ipmat+nlev-1
                   endif
                 endif
c
               mlion=npnxt(mlion)
               enddo
c
             endif

           endif
c
         if  (mlel.ne.0) mlel=npnxt(mlel)
c
         enddo
c
      lpril=0
c     do tranfer.  assumes comp2 and brems have been called
c     already
      call heatt(jkk,lpri,lun11,
     $       t,r,cfrac,delr,xee,xpx,abel,
     $       epi,ncn2,bremsa,
     $       idat1,rdat1,kdat1,nptrs,
     $       npar,npnxt,npfi,npfirst,nplin,nlsvn,
     $       npconi2,ncsvn,rlev,ilev,nlev,klev,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       rcem,oplin,rccemis,opakc,opakscatt,cemab,fline,flinel,
     $       brcems,cmp1,cmp2,httotd,cltotd,hmctotd,
     $             cllines,clcont,htcompd,clcompd,clbremsd)
c
      if (lprif.ne.0) write (lun11,*)'leaving funcsyn'
c
      lprisv=lpri
c
      return
      end
      subroutine fwrtascii(unit,extname,rdati,ncol,
     $                      nidat1j,klabs, kform, kunits,lun11)

c     write an ascii table extension containing
c     ncol columns and nidat1 rows
c     author: T. Bridgman
c
c     parameters:
c        unit    integer            file unit number
c        extname char*30            name of the ascii extension
c        rdati   real(ncol*nidat1j)  data array
c        ncol    integer            number of columns
c        nrhdim  integer            maximum number of rows & columns
c        nidat1j  integer            actual number of rows
c        klabs   char*16(ncol)      column labels
c        kform   char*16(ncol)      column numeric format
c        kunits  char*15(ncol)      column units
c
c     modifications:
c        1998/12/17, wtb: fix fits keyword format problem.  enhanced
c                    parameter list for more flexibility.
c        1999/01/04, wtb: added file creation date, model name, creator
c                    code & checksum
c        1999/01/25, wtb: convert this routine so it just writes an
c                    ascii table extension.
c
      implicit none
      include './PARAM'
      integer nrhmx,nrhmx1

      parameter (nrhmx1=999)
      parameter (nrhmx=3999)

c     passed parameters
      real*8 rdati(nrhmx1,nrhmx)
      real rdat(nrhmx)
      character(16) klabs(nrhmx), kform(nrhmx), kunits(nrhmx)
      integer ncol, nidat1j
!      character(30) extname  !jg
      character(10) extname

      integer unit, status, tfields, nrows, rowlen, verbose,lun11
      integer tbcol(nrhmx),felem,frow,colnum,kk,ll
c
      status=0
      verbose=0
      tfields=ncol
      nrows=nidat1j
      rowlen=0
      tbcol(1)=0
c     append a new empty extension onto the end of the primary array
      call ftcrhd(unit,status)

      if(verbose.gt.0) write(6,*)'fwrtascii: writing header table'
c     write the required header parameters for the ascii table
      call ftphtb(unit,rowlen,nrows,tfields,klabs,tbcol,kform,kunits,
     &            extname,status)
      if (status .gt. 0)call printerror(lun11,status)
c
c     map each column to a 1-d array before writing to the file
      do kk=1,tfields
        if(verbose.gt.0) write(6,*)'fwrtascii: building column ',kk
        frow=1
        felem=1
        colnum=kk
        do ll=1,nidat1j
          rdat(ll)=sngl(rdati(kk,ll))
          enddo
        if(verbose.gt.0) write(6,*)'fwrtascii: writing column ',kk
        call ftpcle(unit,colnum,frow,felem,nrows,rdat,status)
        enddo
      if (status .gt. 0)call printerror(lun11,status)

c     compute checksums
      if(verbose.gt.0) write(6,*)'fwrtascii: writing checksum'
      call ftpcks(unit,status)
c     check for any error, and if so print out error messages
      if (status .gt. 0)call printerror(lun11,status)
c
      end
      subroutine getlun(lun)
      implicit none
      integer lunu, lun
      save lunu   !jg
c
!      common /lpass/lunu
c
      if(lunu.lt.20)lunu=20  !jg
      lunu=lunu+1
      if (lun.eq.-1) then
          lunu=0
      else
          lun=lunu
      endif
c
      return
      end
      subroutine gull1(n,rs,gus,gls,lpri,lun11)
c
c this subroutine calculates the value of |g(n,l;r,l')|**2
c given n the principal qn and r for all l=[o,n-1] and
c l'=l+1 or l'=l-1.  ref burgess (1964), brockelhurst (1971)
c      author: m. bautista
c
      implicit none
c
      integer n,lpri,lun11
      real*8  cn,clu,cll,g0,gu(100),gl(100),fn,pi,s,r
      real*8  dn,dl
      real*8 f1,rs
      real*8 gus(100),gls(100)
      integer n1,l
c
      dn=dfloat(n)
      r=dble(rs)
      pi=2.d0*dacos(0.d+0)
      n1=2*n-1
      call fact(n1,f1)
      g0=0.5d0*dlog(pi/2.d0)+dlog(8.d0*dn)
     $   +dn*dlog(4.d0*dn)-dble(f1)
        if(r.eq.0.d0) then
        gu(n)=g0-2.d0*dn
        else
      s=dsqrt(r)
      gu(n)=g0-2.d0*datan(dfloat(n)*s)/s
     $  -0.5d0*dlog(1.d0-dexp(-2.d0*pi/s))
        end if
      gu(n)=dexp(gu(n))
c
      fn=1.d-300/gu(n)
      gu(n)=gu(n)*fn
c
      if(n.eq.1) go to 40
      gu(n-1)=(2.d0*dn-1.d0)*(1.d0+dn*dn*r)*dn*gu(n)
      gl(n)=(1.d0+dn*dn*r)*gu(n)/(2.d0*dn)
      gl(n-1)=(2.d0*dn-1.d0)*(4.d0+(dn-1.d0)*(1.d0+dn*dn*r))*gl(n)
c
      do 10 l=n-1,3,-1
      dl=dfloat(l)
      gu(l-1)=(4.d0*dn*dn
     $  -4.d0*dl*dl+dl*(2.d0*dl-1.d0)*(1.d0+dn*dn*r))*gu(l)
      gu(l-1)=gu(l-1)-4*dn*dn*(dn-dl)
     $        *(dn+dl)*(1+(dl+1)*(dl+1)*r)*gu(l+1)
      gl(l-1)=(4*dn*dn-4*(dl-1)*(dl-1)
     $        +(dl-1)*(2*dl-1)*(1+dn*dn*r))*gl(l)
      gl(l-1)=gl(l-1)-4*dn*dn*(dn-dl)
     $        *(dn+dl)*(1+(dl-1)*(dl-1)*r)*gl(l+1)
 10    continue
      gl(1)=0.d0
      gu(1)=(4.d0*dn*dn-16.d0+6.d0*(1.d0+dn*dn*r))*gu(2)
      gu(1)=gu(1)-4.d0*dn*dn*(dn-2.d0)*(dn+2.d0)*(1.d0+9.d0*r)*gu(3)
c
      cn=dlog(dn)-dn*dlog(4.d0*dn*dn)
     $    -(2.d0*dn+4.d0)*dlog(1.d0+dn*dn*r)
      gu(1)=cn+dlog(1.d0+r)+2.d0*dlog(gu(1))-2.d0*dlog(fn)
      clu=cn+dlog(1.d0+r)
      cll=cn
c
      do 30 l=1,n-1
      dl=dfloat(l)
      clu=clu+dlog(4.d0*dn*dn*(dn-dl)
     $    *(dn+dl)*(1.d0+(dl+1.d0)*(dl+1.d0)*r))
      cll=cll+dlog(4.d0*dn*dn*(dn-dl)
     $    *(dn+dl)*(1.d0+(dl-1.d0)*(dl-1.d0)*r))
      gu(l+1)=clu+2.d0*dlog(gu(l+1))-2.d0*dlog(fn)
      gl(l+1)=cll+2.d0*dlog(gl(l+1))-2.d0*dlog(fn)
 30    continue
      go to 60
c
 40    gl(1)=0.d0
      gu(1)=2.d0*dlog(gu(1))-dlog(4.d0)
     $      -5.d0*dlog(1.d0+r)-2.d0*dlog(fn)
      if (lpri.ne.0)
     $ write (lun11,*)'in gull1:',r,fn,gu(1),gu(2),gl(1),gl(2)
c converts results to single precision to give in retudn
       do l=1,100
        gus(l)=sngl(gu(l))
        gls(l)=sngl(gl(l))
       enddo

 60    return
      end
      subroutine heatf(jkk,lpri,lun11,
     $       t,r,cfrac,delr,xee,xpx,abel,
     $       epi,ncn2,bremsa,
     $       idat1,rdat1,kdat1,nptrs,
     $       npar,npnxt,npfi,npfirst,nplin,nlsvn,
     $       npconi2,ncsvn,rlev,ilev,nlev,klev,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       rcem,oplin,rccemis,opakc,opakscatt,cemab,fline,flinel,
     $       brcems,cmp1,cmp2,httot,cltot,hmctot,
     $             cllines,clcont,htcomp,clcomp,clbrems, mhd_heat)
c
c     this routine calculates heating and cooling.
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl)
c     pointers to line data
      integer npconi2(ndat2)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
c     energy bins
      real*8 epi(ncn)
c     continuum flux
      real*8 bremsa(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn),brcems(ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     state variables
      real*8 r,t,xpx,delr,delrl
c     heating-cooling variables
c     input parameters
      real*8 xee
      integer ncn2,lpri,lun11
      integer nlsvn,ncsvn
      real*8 fline(2,nnnl),flinel(ncn)
c     level populations
      real*8 abel(nl)
      real*8 cemab(2,nnml)
      character(1) kblnk
      real*8 zrems(4,ncn),zremso(4,ncn)
      real*8 elum(3,nnnl),elumo(3,nnnl)
      real*8 elumab(2,nnml),elumabo(2,nnml)
      real*8 xeltp,cllines,clcont,cmp1,cmp2,cltot,
     $     hmctot,htcomp,clcomp,clbrems,etst,ekt,epiio,
     $     epii,fac,fpr2,optpp,optp2,tmpc1,tautmp,cfrac,
     $     tmpc2,tmpho,xnx,httot,tmpc,tmph,tmpco,hmctmp
      real*8 elin,ener,eth,r19,ergsev,tmpscat,hmctmpo
      real*8 hpctot,hpctmp,hpctmpo
      character(1) klev(100,nd)
      real*8 tmp2,tmp2o
      integer lskp
      real*8 rlev(10,nd)
      integer ilev(10,nd)
      integer idest1,jk,klel,kl,
     $     klion,mlleltp,mllel,mlel,mlion,mt2,nilin,nnz,numcon,
     $     nblin
      integer np1i,np1r,np1k,np1ki
      integer nlev,nlevmx,mltype,ml,mllz,jkk,ltyp,
     $     lrtyp,lcon,nrdt,nidt,nkdt,lk,kkkl,
     $     lprisv,mm,nbinc,mlpar,mlm
c
      real*8 mhd_heat
c
      data kblnk/' '/
      data ergsev/1.602197e-12/
c
      lprisv=lpri
c
      xnx=xpx*xee
      if (lpri.ge.1) lpri=2
      if (lpri.gt.1) write (lun11,*)'in heatt',httot,cltot,delr,r
      if (lpri.gt.1) write (lun11,*)ncsvn
      numcon=ncn2
      r19=r*(1.e-19)
      fpr2=12.56*r19*r19
c
c     comment these out to implement scattering
      clbrems=0.
      lskp=1
      tmp2=0.
      do kl=1,numcon
        tmp2o=tmp2
        tmp2 =brcems(kl)
        if ( kl.ge.2 ) clbrems=clbrems+(tmp2+tmp2o)
     &                *(epi(kl)-epi(kl-lskp))*ergsev/2.
        enddo
c
      delrl=delr
c
      fac=delrl
      ekt = t*(0.861707)
      htcomp = cmp1*xnx*ergsev
      clcomp = ekt*cmp2*xnx*ergsev
      httot=httot+htcomp
      cltot=cltot+clcomp+clbrems
      if (lpri.ge.1) write (lun11,9953)htcomp,clcomp,cmp1,cmp2,
     $   clbrems,httot,cltot
      hmctot=2.*(httot-cltot)/(1.e-37+httot+cltot)
c     print *, hmctot, t
      if (lpri.ge.1) write (lun11,*)hmctot
 9953 format (1h , ' compton heating, cooling=',8e12.4)
      lpri=lprisv
c

      return
      end
      subroutine heatt(jkk,lpri,lun11,
     $       t,r,cfrac,delr,xee,xpx,abel,
     $       epi,ncn2,bremsa,
     $       idat1,rdat1,kdat1,nptrs,
     $       npar,npnxt,npfi,npfirst,nplin,nlsvn,
     $       npconi2,ncsvn,rlev,ilev,nlev,klev,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       rcem,oplin,rccemis,opakc,opakscatt,cemab,fline,flinel,
     $       brcems,cmp1,cmp2,httot,cltot,hmctot,
     $             cllines,clcont,htcomp,clcomp,clbrems)
c
c     this routine calculates heating and cooling.
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl)
c     pointers to line data
      integer npconi2(ndat2)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
c     energy bins
      real*8 epi(ncn)
c     continuum flux
      real*8 bremsa(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn),brcems(ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     state variables
      real*8 r,t,xpx,delr,delrl
c     heating-cooling variables
c     input parameters
      real*8 xee
      integer ncn2,lpri,lun11
      integer nlsvn,ncsvn
      real*8 fline(2,nnnl),flinel(ncn)
c     level populations
      real*8 abel(nl)
      real*8 cemab(2,nnml)
      character(1) kblnk
      real*8 zrems(4,ncn),zremso(4,ncn)
      real*8 elum(3,nnnl),elumo(3,nnnl)
      real*8 elumab(2,nnml),elumabo(2,nnml)
      real*8 xeltp,cllines,clcont,cmp1,cmp2,cltot,
     $     hmctot,htcomp,clcomp,clbrems,etst,ekt,epiio,
     $     epii,fac,fpr2,optpp,optp2,tmpc1,tautmp,cfrac,
     $     tmpc2,tmpho,xnx,httot,tmpc,tmph,tmpco,hmctmp
      real*8 elin,ener,eth,r19,ergsev,tmpscat,hmctmpo
      real*8 hpctot,hpctmp,hpctmpo
      character(1) klev(100,nd)
      real*8 tmp2,tmp2o
      integer lskp
      real*8 rlev(10,nd)
      integer ilev(10,nd)
      integer idest1,jk,klel,kl,
     $     klion,mlleltp,mllel,mlel,mlion,mt2,nilin,nnz,numcon,
     $     nblin
      integer np1i,np1r,np1k,np1ki
      integer nlev,nlevmx,mltype,ml,mllz,jkk,ltyp,
     $     lrtyp,lcon,nrdt,nidt,nkdt,lk,kkkl,
     $     lprisv,mm,nbinc,mlpar,mlm
c
c
      data kblnk/' '/
      data ergsev/1.602197e-12/
c
c
      lprisv=lpri
c
      xnx=xpx*xee
      if (lpri.ge.1) lpri=2
      if (lpri.gt.1) write (lun11,*)'in heatt',httot,cltot,delr,r
      if (lpri.gt.1) write (lun11,*)ncsvn
      numcon=ncn2
      r19=r*(1.e-19)
      fpr2=12.56*r19*r19
c
c     comment these out to implement scattering
      clbrems=0.
      lskp=1
      tmp2=0.
      do kl=1,numcon
        tmp2o=tmp2
        tmp2 =brcems(kl)
        if ( kl.ge.2 ) clbrems=clbrems+(tmp2+tmp2o)
     &                *(epi(kl)-epi(kl-lskp))*ergsev/2.
        enddo
c
c      delrl=max(delr,1.e-24*r)
        delrl=delr
        httot=0.
        cltot=0.
        hmctot=0.
        hpctot=0.
        epii=epi(1)
        hmctmp=0.
        hpctmp=0.
        tmpc=0.
        tmph=0.
        do kl=1,numcon
          optpp=opakc(kl)
          optp2=max(1.e-34,optpp)
c         for outward only
          epiio=epii
          epii=epi(kl)
          tmpho=tmph
          tautmp=optp2*delrl
          fac=1.
          if (tautmp.gt.0.01)
     $      fac=(1.-exp(-tautmp))/tautmp
          tmph=bremsa(kl)*optp2
          tmpco=tmpc
          tmpc1=rccemis(1,kl)
          tmpc2=rccemis(2,kl)
c          tmpc1=rccemis(1,kl)+brcems(kl)*(1.-cfrac)/2.
c          tmpc2=rccemis(2,kl)+brcems(kl)*(1.+cfrac)/2.
          tmpc=(tmpc1+tmpc2)*12.56
          hmctmpo=hmctmp
          hpctmpo=hpctmp
          hmctmp=(tmph-tmpc)*fac
          hpctmp=(tmph+tmpc)*fac
c         testing lte
c          zrems(1,kl)=12.56*(tmpc1+tmpc2)*fpr2/(1.e-34+optp2)
c         this is the good expression
          zrems(1,kl)=max(0.,zremso(1,kl)
     $         -(tmph-12.56*(tmpc1+tmpc2)-flinel(kl))*fac*delrl*fpr2
     $        )
          zrems(2,kl)=zremso(2,kl)+12.56*tmpc1*fac*delrl*fpr2
          zrems(3,kl)=zremso(3,kl)+12.56*tmpc2*fac*delrl*fpr2
          tmpscat=opakscatt(kl)*bremsa(kl)*fac
          zrems(4,kl)=zremso(4,kl)+tmpscat*delrl*fpr2
          if (kl.gt.1) then
            hmctot=hmctot+(hmctmp+hmctmpo)
     $       *(epii-epiio)*ergsev/2.
            hpctot=hpctot+(hpctmp+hpctmpo)
     $       *(epii-epiio)*ergsev/2.
            endif
          httot=(hmctot+hpctot)/2.
          cltot=(-hmctot+hpctot)/2.
          if (lpri.ge.1) write (lun11,9009)kl,epii,optpp,
     $     bremsa(kl),tmph,tmpc,flinel(kl),httot,cltot
     $       ,rccemis(1,kl)+rccemis(2,kl),hmctot,tautmp,fac
c     $       ,zrems(1,kl),zrems(2,kl),zrems(3,kl)
 9009     format (1x,i6,15(1pe12.5))
          enddo
        if (lpri.ge.1)
     $   write (lun11,*)'continuum heating, cooling:',
     $       hmctot,hpctot,(hmctot+hpctot)/2.,-(hmctot-hpctot)/2.
        clcont=cltot
        do jkk=1,nlsvn
          jk=jkk
          ml=nplin(jk)
          mlm=ml-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          elin=abs(rdat1(np1r))
          if ((elin.lt.1.e+8).and.(elin.gt.1.)
     $      .and.(ml.ne.0).and.(lrtyp.eq.4)) then
            nilin=npar(ml)
            ener=ergsev*(12398.41)/max(elin,1.e-24)
            etst=ener/ergsev
            nblin=nbinc(etst,epi,ncn2)
c           this statement is in v221boxx
            optp2=0.
c            optp2=opakc(nblin)
            tautmp=optp2*delrl
            fac=1.
c            if (tautmp.gt.1.e-4)
c     $        fac=(1.-exp(-tautmp))/tautmp
            tmph=elumo(1,jk)*optp2/fpr2
            tmpc1=rcem(1,jk)
            tmpc2=rcem(2,jk)
            tmpc=tmpc1+tmpc2
c           must remove these statements when radiative excitation is included
            hmctot=hmctot-tmpc*fac
            hpctot=hpctot+tmpc*fac
            cltot=cltot+tmpc
            elum(1,jk)=max(0.,elumo(1,jk)+(-tmph+tmpc1)*fac*delrl*fpr2)
            tmph=elumo(2,jk)*optp2/fpr2
            elum(2,jk)=max(0.,elumo(2,jk)+(-tmph+tmpc2)*fac*delrl*fpr2)
            if (lpri.ge.1)
     $       write (lun11,9019)jk,nilin,nblin,etst,
     $        rcem(1,jk),rcem(2,jk),delrl,fpr2,cltot,
     $        elumo(2,jk),elum(2,jk),optp2,tmph
c            write (lun11,*)'line sum',jk,nilin,nblin,etst,
c     $        rcem(2,jk),optp2,fac,delrl,fpr2,tmph,tmpc2,
c     $        elumo(2,jk),elum(2,jk),(elum(2,jk)-elumo(2,jk))/delrl/fpr2
 9019     format (1x,3i6,14(1pe12.4))
            endif
          enddo
c
c       calculate rrc luminosities
C       First look for element data (jk is element index)
          if (lpri.ge.2)
     $     write (lun11,*)'rrc print:'
          klel=11
          mlel=npfirst(klel)
          jk=0
          jkk=0
          do while (mlel.ne.0)
            jk=jk+1
            mt2=mlel-1
            call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,
     $        nptrs,0,lun11)
            if (nidt.gt.0) then
              mllel=idat1(np1i+nidt-1)
              xeltp=rdat1(np1r)
              xeltp=abel(mllel)
              nnz=idat1(np1i)
              if (lpri.ge.2)
     $        write (lun11,*)'element:',jk,mlel,mllel,nnz,
     $                  (kdat1(np1k-1+mm),mm=1,nkdt)
C           ignore if the abundance is small
            if (xeltp.lt.1.e-10) then
                jkk=jkk+nnz
              else
c               now step thru ions (jkk is ion index)
                klion=12
                mlion=npfirst(klion)
                jkk=0
                kl=0
                do while ((mlion.ne.0).and.(kl.lt.nnz))
                  jkk=jkk+1
C                 retrieve ion name from kdati
                  mlm=mlion-1
                  call drd(ltyp,lrtyp,lcon,
     $              nrdt,np1r,nidt,np1i,nkdt,np1ki,mlm,
     $              nptrs,0,lun11)
C                 if not accessing the same element, skip to the next element
                  mlleltp=idat1(np1i+nidt-2)
                  if (mlleltp.eq.mllel) then
                    kl=kl+1
                    if (lpri.ge.2)
     $              write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
     $                          (kdat1(np1ki+mm-1),mm=1,nkdt)
c                   now find level data
c                   step thru types
                    nlevmx=0
                    mltype=13
                    ml=npfi(mltype,jkk)
                    mllz=0
                    if (ml.ne.0) mllz=npar(ml)
c                   step thru records of this type
                    mlpar=0
                    if (ml.ne.0) mlpar=npar(ml)
                    do while ((ml.ne.0).and.(mlpar.eq.mllz))
                      mlm=ml-1
                      call drd(ltyp,lrtyp,lcon,
     $                  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $                  nptrs,0,lun11)
                      nlev=idat1(np1i+nidt-2)
                      nlevmx=max(nlevmx,nlev)
                      if ((nlev.gt.0).and.(nlev.le.nd)) then
                        do  lk=1,nrdt
                          rlev(lk,nlev)=rdat1(np1r+lk-1)
                          enddo
                        do lk=1,nidt
                          ilev(lk,nlev)=idat1(np1i+lk-1)
                          enddo
                        do lk=1,nkdt
                          klev(lk,nlev)=kdat1(np1k+lk-1)
                          enddo
                        do lk=nkdt+1,20
                          klev(lk,nlev)=kblnk
                          enddo
                        endif
                      ml=npnxt(ml)
                      if (ml.ne.0) mlpar=npar(ml)
                      enddo
                    nlev=nlevmx
                    mltype=7
                    ml=npfi(mltype,jkk)
                    mllz=0
                    if (ml.ne.0) mllz=npar(ml)
                    mlpar=0
                    if (ml.ne.0) mlpar=npar(ml)
                    do while ((ml.ne.0).and.(mlpar.eq.mllz))
c                     step thru records of this type
                      mlm=ml-1
                      call drd(ltyp,lrtyp,lcon,
     $                  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $                  nptrs,0,lun11)
                      kkkl=npconi2(ml)
                      idest1=idat1(np1i+nidt-2)
                      if ((kkkl.gt.0).and.(kkkl.le.ndat2)
     $                .and.((cemab(1,kkkl).gt.1.e-36)
     $                .or.(cemab(2,kkkl).gt.1.e-36))) then
                        eth=rlev(4,idest1)-rlev(1,idest1)
                        tmpc=(cemab(1,kkkl)+cemab(2,kkkl))
                        optp2=0.
                        fac=delrl
                        tmph=elumabo(1,kkkl)*optp2/fpr2
                        elumab(1,kkkl)=max(0.,elumabo(1,kkkl)
     $                            +(-tmph+tmpc)*fac*fpr2/2.)
                        elumab(2,kkkl)=max(0.,elumabo(2,kkkl)
     $                            +(-tmph+tmpc)*fac*fpr2/2.)
                        if (lpri.ge.2)
     $                  write (lun11,981)kkkl,eth,idest1,
     $                    cemab(1,kkkl),cemab(2,kkkl),
     $                    elumabo(1,kkkl),elumab(1,kkkl)
 981                    format (1x,i6,1pe11.3,i6,6(1pe11.3))
                        endif
                      ml=npnxt(ml)
                      if (ml.ne.0) mlpar=npar(ml)
                      enddo
                    endif
C                 Go to next ion
                  mlion=npnxt(mlion)
                  enddo
                endif
              endif
            mlel=npnxt(mlel)
C           Go to next element
            enddo
        cllines=cltot-clcont
        if (lpri.ge.1)
     $   write (lun11,*)'line cooling',cltot,cllines
c
c
      fac=delrl
      ekt = t*(0.861707)
      htcomp = cmp1*xnx*ergsev
      clcomp = ekt*cmp2*xnx*ergsev
c      these statements needed to implement scattering
c      httot=htcomp*fac*fpr2+httot
c      cltot=clcomp*fac*fpr2+cltot
c      hmctot=hmctot+(htcomp-clcomp)
c      hpctot=hpctot+(htcomp+clcomp)
      hmctot=hmctot+(htcomp-clcomp-clbrems)
      hpctot=hpctot+(htcomp+clcomp+clbrems)
      httot=httot+htcomp
      cltot=cltot+clcomp+clbrems
      if (lpri.ge.1) write (lun11,9953)htcomp,clcomp,cmp1,cmp2,
     $   clbrems,(hmctot+hpctot)/2.,-(hmctot-hpctot)/2.
      hmctot=2.*hmctot/(1.e-37+hpctot)
      if (lpri.ge.1) write (lun11,*)hmctot
 9953 format (1h , ' compton heating, cooling=',8e12.4)
      lpri=lprisv
c

      return
      end
      subroutine hgf(ia,ib,ic,x,hyp)
c
c     subroutine hgf calculates the value, hyp, of the
c     hypergeometric fn at x for constants ia,ib,ic
c     real*8  ser,hyp
c     author:  M. Bautista
c
      implicit none
      integer ia,ib,ic,i,j,n
      real*8 x,hyp,ser
c
      ser=1.
      hyp=1.
      i=-ia
      j=-ib
      i=min(i,j)
      do 10 n=0,i
      ser=ser*(ia+n)*(ib+n)*x/((n+1.)*(ic+n))
      hyp=hyp+ser
 10    continue
c
      return
      end
      subroutine hphotx(ener,ic,nq,xsec,lun11,lpri)
c
c     ener   is the photon energy in ryds with respect to the ionization
c         threshold
c     xsec   is an array containing the cross section in mb (18^{-18} cm^2)
c         for all l=[0,nq-1]
c     ic     ion charge
c     np     principal quantum number
c     ll     angular momentum number
c        real*8  en,cons,r,rk,theta1,theta2,gu,gl
c     author:  M. Bautista
c
c
      implicit none
c
      real*8 gu(100),gl(100),xsec(100)
      real*8 cons,ener,en,r,rk,theta1,theta2
      integer ic,nq,lun11,lpri,lm
c
      cons=.54492*acos(0.)
c        write (lun11,*)'in hphotx:',ener,ic,nq,cons
        en=ener
c        r=dsqrt(en)
        r=sqrt(en)
        rk=r/float(ic*ic)
        if (lpri.ne.0)
     $   write (lun11,*)'before call gull1:',nq,rk,r,en
        call gull1(nq,rk*rk,gu,gl,lpri,lun11)
c        call gull1(nq,rk,gu,gl,lun11)
       do lm=0,nq-1
        theta1=(1.+nq*nq*rk*rk)*exp(gu(lm+1))
        theta2=(1.+nq*nq*rk*rk)*exp(gl(lm+1))
        if (lpri.ne.0)
     $   write (lun11,*)'after call gull1:',lm,gu(lm+1),gl(lm+1),
     $           theta1,theta2
        xsec(lm+1)=cons*((lm+1)*theta1+lm*theta2)/(2.*lm+1.)
        xsec(lm+1)=float(nq*nq)/float(ic*ic)*xsec(lm+1)
       enddo
c
      return
      end
      subroutine hunt3(xx,n,x,jlo,lpri,lun11)
c
      implicit none
c
      integer n,jlo,lpri,lun11,nint,jhi,inc,jm,nintmx
      real*8 xx(n),x
      logical ascnd
c
      data nintmx/1000/
c
      nint=0
      if (lpri.gt.1) write (lun11,*)'in hunt',n,x
c      do m=1,n
c        write (lun11,*)m,xx(m)
c        enddo
      jlo=1
      ascnd=.false.
      if (xx(n).gt.xx(1)) ascnd=.true.
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=1
        jhi=n+1
        if (lpri.gt.1) write (lun11,*)'initializing',jlo,jhi
        go to 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if (lpri.gt.1) write (lun11,*)'hunt up ',jlo,jhi,xx(jlo),xx(jhi)
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          if (lpri.gt.1) write (lun11,*)'double the increment',inc
          go to 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if (lpri.gt.1) write (lun11,*)'hunt down ',jlo,jhi,xx(jlo),
     $                        xx(jhi)
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          if (lpri.gt.1) write (lun11,*)'double the increment',inc
          go to 2
        endif
      endif
      jm=n/2
3     continue
      if (lpri.gt.1) write (lun11,*)'bisection phase',jlo,jhi,jm,xx(jm)
      jlo=min(jlo,n)
      jlo=max(jlo,1)
      nint=nint+1
      if (lpri.gt.1) write (lun11,*)'bisection phase',jlo,jhi,jm,xx(jm)
      if ((jhi-jlo.eq.1).or.(nint.gt.nintmx)) return
      if (lpri.gt.1) write (lun11,*)'bisection phase',jlo,jhi,jm,xx(jm)
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      go to 3
c
      end
      subroutine huntf(xx,n,x,jlo,lpri,lun11)
c
c     this version of hunt assumes equally spaced data in log
c     author:  T. Kallman
c
      implicit none
c
      integer n,jlo,lpri,lun11
      real*8 xx(n),x,xtmp,tst,tst2
c
      xtmp=max(x,xx(2))
      jlo=int((n-1)*log(xtmp/xx(1))/log(xx(n)/xx(1)))+1
      if (jlo.lt.n) then
        tst=abs(log(x/xx(jlo)))
        tst2=abs(log(x/xx(jlo+1)))
        if (tst2.lt.tst) jlo=jlo+1
        endif
      jlo=max(1,jlo)
      jlo=min(n,jlo)
      if (lpri.ne.0)
     $  write (lun11,*)'in huntf',n,xx(1),xx(n),jlo,xx(jlo),x
c
      return
      end
      subroutine impact(en,l,temp,ic,z1,rm,ecm,psi,cr)
c
c impact parameter collision cross-sections using the method of seaton.
c     author:  M. Bautista
c
      implicit none
c
      real*8 en,temp,z1,rm,ecm,psi,cr
      integer l,ic
c
      real*8  b,xsi,phi,bo,xsw,phw,del
      real*8 tk,fi,wo,ev,po,w,wi,ff,crinc,ric
      integer inc,jm,j
c
      tk=8.617e-5*temp
      inc=1
      jm=90*inc
      cr=0.
      fi=0.
      wo=0.
      b=10.d0
      ric=float(ic)
      ev=abs(ecm)/8065.48
c
      po=(3.*en*en-l*(l+1))/2./ic
c
c strong coupling
c
 21    del=b/100.d0/inc
      do 20 j=1,jm
      b=b-del
      call impcfn(b,xsi,phi)
c     write (lun11,*)b,xsi,phi
      w=ric*rm*ev/sngl(b)*sqrt(2.*sngl(xsi)*psi)
      wi=w+ecm/8065.48/2.
c     write (lun11,*)wi/tk
      if(wi/tk.ge.100.) go to 13
      if(wi.le.0.) go to 20
c
c weak coupling
c
      bo=dble(po*ev/2./w*sqrt(wi*rm/13.60))
      call impcfn(bo,xsw,phw)
c
c the minimum of the weak and strong coupling x-sections is used
      ff=min(sngl(xsi/2.d0+phi),sngl(phw))
      ff=ff*exp(-wi/tk)
c
      crinc=(fi+ff)/2.*(wi-wo)
      cr=crinc+cr
      if(cr.lt.1.e-20) go to 20
      fi=ff
      wo=wi
      if(crinc/cr.lt.1.e-5) go to 13
c
 20    continue
      go to 21
 13       cr=6.900e-5*z1*z1*sqrt(rm/temp)*psi*cr/tk
c
      return
      end
      subroutine impactn(n,m,temp,ic,amn,cmm,lun11,lpri)
c
c impact parameter collision cross-sections using the method of seaton.
c impactn.ftn calculates the electron collisional excitation rate for
c transitions between principal quantum number n and m in hydrogenic
c atoms with ionic charge ic.  it is assumed that rm=1 and z1=1.
c cmm is the symmetrical quantity used in the models.
c     author:  M. Bautista
c
      implicit none
c
      real*8 temp,ecm,psi,cr,amn,cmm
      integer n,m,lun11,lpri,ic
c
      real*8  b,xsi,phi,bo,xsw,phw,del
      real*8 xm,rm,z1,tk,ecm3,po,fi,wo,ev,wi,w,ff,crinc
      integer inc,jm,j
c
c
      if (lpri.ne.0)
     $ write (lun11,*)'in impactn:',n,m,temp,ic,amn
      cmm=0.
      xm=157888.*float(ic*ic)/temp/float(m*m)
      if(xm.gt.60) return
      rm=1.
      z1=1.
      tk=8.617e-5*temp
      inc=1
      jm=90*inc
c
      ecm=109737.*float(ic*ic)*(1./float(n*n)-1./float(m*m))
      ecm3=ecm**3
      ecm=-ecm
      psi=1.644e+5*amn/ecm3
       po=(5.*float(n*n)+1)/4./float(ic)
c
      cr=0.
      fi=0.
      wo=0.
      b=10.d0
      ev=abs(ecm)/8065.48
c
      if (lpri.ne.0)
     $ write (lun11,*)xm,tk,ecm,ecm3,psi,po,ev
c
c strong coupling
c
 21    del=b/100.d0/dfloat(inc)
      do 20 j=1,jm
      b=b-del
      call impcfn(b,xsi,phi)
      w=float(ic)*rm*ev/sngl(b)*sqrt(2.*sngl(xsi)*psi)
      wi=w+ecm/8065.48/2.
      if (lpri.ne.0)
     $ write (lun11,*)'in 20 loop:',j,b,xsi,phi,w,wi
      if(wi/tk.ge.100.) go to 13
      if(wi.le.0.) go to 20
c
cc weak coupling
cc
       bo=dble(po*ev/2./w*sqrt(wi*rm/13.60))
       call impcfn(bo,xsw,phw)
cc
cc the minimum of the weak and strong coupling x-sections is used
       ff=min(sngl(xsi/2.d0+phi),sngl(phw))
c
c only the strong coupling calculation is used
      ff=sngl(xsi/2.d0+phi)
      ff=ff*exp(-wi/tk)
c
      crinc=(fi+ff)/2.*(wi-wo)
      cr=crinc+cr
      if(cr.lt.1.e-20) go to 20
      fi=ff
      wo=wi
      if (lpri.ne.0)
     $ write (lun11,*)'weak coupling:',bo,xsw,phw,ff,crinc,cr
      if((crinc/cr.lt.1.e-5).and.(crinc.gt.1.e-7)) go to 13
c
 20    continue
      go to 21
 13       cr=6.900e-5*z1*z1*sqrt(rm/temp)*psi*cr/tk
      cmm=cr*m*m*exp(xm)
c
      if (lpri.ne.0)
     $ write (lun11,*)'done with impactn:',cr,cmm
c
      return
      end
      subroutine impcfn(x,xsi,phi)
c
c data for functions used in the impact parameter method are generated
c using polynomials fitted to seaton's (1962) values using least square
c     author:  M. Bautista
c
      implicit none
c
      real*8  a(6),b(6),x,xsi,phi,pi,y
      integer n
c
      pi=2.d0*dacos(0.d0)
      a(1)=0.9947187d0
      a(2)=0.6030883d0
      a(3)=-2.372843d0
      a(4)=1.864266d0
      a(5)=-0.6305845d0
      a(6)=8.1104480d-02
      b(1)=0.2551543d0
      b(2)=-0.5455462d0
      b(3)=0.3096816d0
      b(4)=4.2568920d-02
      b(5)=-2.0123060d-02
      b(6)=-4.9607030d-03
c
      if(x.gt.2.d0) go to 25
      xsi=0.d0
      phi=0.d0
      do 20 n=1,6
      xsi=xsi+a(n)*x**(n-1)
      y=dlog(x)
      phi=phi+b(n)*y**(n-1)
 20    continue
      if(x.eq.1.d0) phi=b(1)
      if(x.lt.0.05d0) then
      xsi=1.0d0+0.01917d0/0.05d0*x
      y=dlog(1.1229d0/x)
      phi=y+x*x/4.d0*(1.d0-2.d0*y*y)
      endif
      go to 30
c
 25    xsi=pi*x*dexp(-2.d0*x)*(1.d0+0.25d0/x+1.d0/32.d0/x/x)
      phi=pi/2.d0*dexp(-2.d0*x)*(1.d0+0.25d0/x-3.d0/32.d0/x/x)
c
 30    continue
c
      return
      end
      subroutine init(lunlog,bremsa,bremsint,tau0,dpthc,tauc,
     $   xii,rrrt,pirt,htt,cll,httot,cltot,
     $   cllines,clcont,htcomp,clcomp,clbrems,
     $   xilev,rcem,oplin,rccemis,brcems,opakc,opakscatt,
     $   cemab,cabab,opakab,elumab,elumabo,elum,elumo,
     $   zrems,zremso,rates,vsav,idrates,fline,flinel)
c
c     this routine initializes everything
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
c     line luminosities
      real*8 elum(3,nnnl),elumo(3,nnnl)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
c     line optical depths
      real*8 tau0(2,nnnl)
      real*8 fline(2,nnnl),flinel(ncn)
c     continuum lum
      real*8 zrems(4,ncn),zremso(4,ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn),brcems(ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     level populations
      real*8 xilev(nnml)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
      real*8 elumab(2,nnml),elumabo(2,nnml)
      real*8 tauc(2,nnml)
c     ion abundances
      real*8 xii(nni)
c     heating/cooling
      real*8 htt(nni),cll(nni)
      real*8 rrrt(nni),pirt(nni)
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
      real*8 httot,cltot,cllines,clcont,htcomp,clcomp,clbrems
      integer i,j,lunlog
      character(133) tmpst
c
      httot=0.
      cltot=0.
      cllines=0.
      clcont=0.
      htcomp=0.
      clcomp=0.
      clbrems=0.
c
      do i = 1,ncn
         rccemis(1,i)=0.
         rccemis(2,i)=0.
         brcems(i)=0.
         flinel(i)=0.
         zrems(1,i)=0.
         zrems(2,i)=0.
         zrems(3,i)=0.
         zrems(4,i)=0.  !jg
         zremso(1,i)=0.
         zremso(2,i)=0.
         zremso(3,i)=0.
         zremso(4,i)=0. !jg
         bremsint(i)=0.
         bremsint(i)=0.
         bremsa(i)=0.
         dpthc(1,i) = 0.
         dpthc(2,i)=0.
c         dpthc(2,i)=1.e+10
         opakc(i)=0.
         opakscatt(i)=0.
         enddo
       do  i = 1,nnnl
         fline(1,i)=0.
         fline(2,i)=0.
         rcem(1,i)=0.
         rcem(2,i)=0.
         elum(1,i)=0.
         elum(2,i)=0.
         elum(3,i)=0.
         elumo(1,i)=0.
         elumo(2,i)=0.
         elumo(3,i)=0.
c         tau0(2,i)=1.e+20
c         tau0(1,i) = 1.e+20
         tau0(2,i)=0.
         tau0(1,i) = 0.
         oplin(i)=0.
         enddo
c      write (lunlog,*)'NB no backward escape'
c      write (tmpst,*)'NB no backward escape'
c      call xwrite(tmpst,10)
      do i=1,nni
         xii(i)=0.
         htt(i)=0.
         cll(i)=0.
         rrrt(i)=0.
         pirt(i)=0.
         enddo
       do i = 1,nnml
         elumab(1,i)=0.
         elumab(2,i)=0.
         elumabo(1,i)=0.
         elumabo(2,i)=0.
         cabab(i)=0.
         cemab(1,i)=0.
         cemab(2,i)=0.
         opakab(i)=0.
c         xilev(i)=1.
         xilev(i)=0.
         tauc(1,i) = 0.
         tauc(2,i) =0.
         enddo
      do i=1,ndat2
         do j=1,4
           rates(j,i)=0.
           vsav(j,i)=0.
           enddo
         idrates(1,i)=0
         idrates(2,i)=0
         enddo
c
      return
      end
      subroutine intin(x1,x2,x0,t,ri2,ri3,lpri,lun11)
c
c     this routine does the integrals needed by milne
c     author:  M. Bautista
c
      implicit none
c
      real*8  x1,x2,x0,t,ri2,ri3
      real*8  ryk,s1,s2,s0,del,rr
      integer lpri,lun11
c
       ryk=7.2438d+15
       s1=x1*ryk/t
       s2=x2*ryk/t
       s0=x0*ryk/t
       del=ryk/t
       if (lpri.gt.1)
     $  write (lun11,*)'in intin:',s1,s2,s0,del
       if ((s1-s0).lt.90.d0) then
c           ri2=dexpo(s0-s1)*(s1*s1+2.*s1+2.)-dexpo(s0-s2)*(s2*s2+2.*s2+2.)
           ri2=dexp(s0-s1)*((s1*s1+2.d0*s1+2.d0)
     $        -dexp(s1-s2)*(s2*s2+2.d0*s2+2.d0))/del/dsqrt(del)
           if (lpri.gt.1)
     $     write (lun11,*)'ri2=',ri2
           if ((s0.lt.1.d-3).and.(s2.lt.1.d-3).and.(s1.lt.1.d-3))
     $       ri2=0.d0
         else
           ri2=0.d0
         endif
c       ri2=ri2/(del**1.5)
       if (lpri.gt.1)
     $     write (lun11,*)'ri2=',ri2
c       rr=dexpo(s0-s1)*(s1**3)-dexpo(s0-s2)*(s2**3)
       rr=dexp(s0-s1)*((s1**3)-dexp(s1-s2)*(s2**3))
       if (lpri.gt.1)
     $     write (lun11,*)'rr=',rr
       ri3=(rr/del/dsqrt(del)+3.d0*ri2)/del
       if (lpri.gt.1)
     $  write (lun11,*)'in intin:',s1,s2,s0,del,ri2,rr,ri3
       return
       end
      subroutine ioneqm(z,a,s,n,m,l,lpri,lun11)
      implicit none
c
c     this routine computes ionization equilibrium
c     solves a system of ionization equations, attempting
c     to avoid overflow problems.
c     author:  T. Kallman
c
c
      integer m, n
      real*8  z(m),a(m),s(n),q(31)
      real*8  eps,delt,pl,suml,tst
      real*8  sumg, pg
      integer i,j,jk,jmax,k,l,ll,lpri,mmn
      integer mmx,lun11
      integer lprisv
c
      data eps/1.e-6/
      data delt/1.e-28/
c
      lprisv=lpri
      if (lpri.ge.1) lpri=2
      if (lpri.ge.2) write (lun11,*)'in ioneqm'
c
c     initialize
      do jk = 1,n
         s(jk) = 0.d0
         enddo
c
      if ( lpri.ge.2 ) write (lun11,99001)
c
c     form naive ratio
      do j = 1,m
         q(j) = a(j)/(z(j)+delt)
         enddo
c
c     step thru and search for max. q value
      jk = l
 300  jk = jk + 1
      if ( (jk.lt.n) .and. (q(jk-1).lt.1.d0) ) goto 300
      jmax = jk
c
      if ( lpri.ge.2 ) write (lun11,99002) n,m,l,jmax
c
c     step forwards
      suml = 0.d0
      if ( jmax.ne.n ) then
         pl = 1.d0
         mmx = jmax - 1
 350     mmx = mmx + 1
         pl = pl/(q(mmx)+delt)
         suml = suml + pl
         tst = pl/(suml+delt)
         if ( lpri.ge.2 ) write (lun11,99003) mmx,q(mmx),suml,pl
         if ( (tst.gt.eps) .and. (mmx.lt.m) ) goto 350
      endif
c
c     step backwards
      sumg = 0.d0
      if ( jmax.ne.l ) then
         pg = 1.d0
         mmn = jmax
 400     mmn = mmn - 1
         pg = pg*q(mmn)
         sumg = sumg + pg
         tst = pg/(sumg+delt)
         if ( lpri.ge.2 ) write (lun11,99004) mmn,q(mmn),sumg,pg
         if ( (tst.gt.eps) .and. (mmn.gt.l) ) goto 400
      endif
c
c
      s(jmax) = 1.d0/(1.d0+suml+sumg)
      if ( jmax.ne.n ) then
         do j = jmax,mmx
            s(j+1) = s(j)/(q(j)+delt)
            enddo
        endif
c
      if ( jmax.ne.l ) then
         k = jmax - mmn
         do i = 1,k
            j = jmax - i
            s(j) = s(j+1)*q(j)
            enddo
        endif
c
      if ( lpri.ge.2 ) write (lun11,99005) (ll,q(ll),s(ll),ll=l,n)
c
      lpri=lprisv
c
      return
99001 format (' ',' in ioneqm ')
99002 format (' ',' n,m,l,jmax --',4i4)
99003 format (' ','in greater than loop, j,q,sum,p --',i4,3e12.4)
99004 format (' ','in less than loop, j,q,sum,p --',i4,3e12.4)
99005 format (' ',i4,2e12.4)
      end
      subroutine irc(n,t,rc,rno,se,lpri,lun11)
c
c irc calculates the excitation rate, se [cm**3/s], for ionization
c of hydrogen atoms from state n due to electron collisions, assuming
c the continuum starts at level rno.  the
c energy loss rate, elost [ev*cm**3/s], is also determined.
c cin is the 3-body recombination rate, determined from cni by
c detailed balance.
c ref. johnson (1972)
c     author:  m. bautista
c
      implicit none
c
      real*8 t,rc,rno,se
      integer n
      real*8 xo,yn,an,bn,rn,g0,g1,g2,zn,ey
      real*8 ez
      integer lpri,lun11
c
      if (lpri.ne.0) write (lun11,*)'in irc',
     $ n,t,rc,rno
      if(rc.ne.1.) then                          ! mab
       call szirc(n,t,rc,rno,se,lpri,lun11)
       return
      endif
c
      xo=1.-n*n/rno/rno
      yn=xo*157803./(t*n*n)
      if(n-2) 100,200,300
 100   an=1.9603*n*(1.133/3./xo**3-0.4059/4./xo**4+0.07014/5./xo**5)
      bn=2./3.*n*n/xo*(3.+2./xo-0.603/xo/xo)
      rn=0.45
      go to 400
c
 200   an=1.9603*n*(1.0785/3./xo**3-0.2319/4./xo**4+0.02947/5./xo**5)
      bn=(4.-18.63/n+36.24/(n*n)-28.09/(n*n*n))/n
      bn=2./3.*n*n/xo*(3.+2./xo+bn/xo/xo)
      rn=0.653
      go to 400
c
 300   g0=(0.9935+0.2328/n-0.1296/(n*n))/3./xo**3
      g1=-(0.6282-0.5598/n+0.5299/(n*n))/(n*4.)/xo**4
      g2=(0.3887-1.181/n+1.470/(n*n))/(n*n*5.)/xo**5
      an=1.9603*n*(g0+g1+g2)
      bn=(4.-18.63/n+36.24/(n*n)-28.09/(n*n*n))/n
      bn=(3.+2./xo+bn/xo/xo)*2.*n*n/3./xo
      rn=1.94*n**(-1.57)
c
 400   continue
      rn=rn*xo
      zn=rn+yn
      call expint(yn,ey)
      call expint(zn,ez)
      se=an*(ey/yn/yn-exp(-rn)*ez/zn/zn)
      ey=1.+1./yn-ey*(2./yn+1.)
      ez=exp(-rn)*(1.+1./zn-ez*(2./zn+1.))
      se=se+(bn-an*log(2.*n*n/xo))*(ey-ez)
      se=se*sqrt(t)*yn*yn*n*n*1.095e-10/xo
      if (lpri.ne.0) write (lun11,*)'in irc',
     $ xo,yn,an,bn,rn,zn,ey,ez,se
c      cii=se*n*n
c mab
c     cii=cii/100.
c      cni=se*exp(-yn)
c      cin=4.144219e-16*n*n/t/sqrt(t)*se
c      elost=se*13.60/(n*n)
c
      return
      end
      subroutine ispcg2(zremsz,epi,ncn2,enlum,lpri,lun11)
c
c     this subroutine calculates number luminosity
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
      real*8 zremsz(ncn),epi(ncn)
      integer ncn2,lpri,lun11
      real*8 enlum
      real*8 sum2,sum3,sum4,sum5
      integer jk
      integer numcon
c
      if (lpri.ge.1) write (lun11,*)'in ispec2'
      numcon=ncn2
      sum2 = 0.
      sum3 = 0.
      sum4 = 0.
      sum5 = 0.
      do jk = 1,numcon
         if (jk.gt.1)
     $     sum5 = sum5+(zremsz(jk)+zremsz(jk-1))
     &             *(epi(jk)-epi(jk-1))/2.
         if ( epi(jk).ge.13.6 ) then
            sum2 = sum2+(zremsz(jk)/epi(jk)+zremsz(jk-1)/epi(jk-1))
     &             *(epi(jk)-epi(jk-1))/2.
           if ( epi(jk).le.24.48 )
     $        sum3 = sum3+(zremsz(jk)/epi(jk)+zremsz(jk-1)/epi(jk-1))
     &             *(epi(jk)-epi(jk-1))/2.

         endif
         if ((epi(jk).ge.24.48).and.(epi(jk).le.54.4))
     $     sum4 = sum4+(zremsz(jk)/epi(jk)+zremsz(jk-1)/epi(jk-1))
     &             *(epi(jk)-epi(jk-1))/2.
          if (lpri.ge.1)
     $     write (lun11,*)jk,epi(jk),zremsz(jk),sum2
          enddo
      enlum = sum2
c     write (lun11,*)'U(1-1.8),U(1.8-4):',sum3,sum4
c     write (lun11,*)'Lbol=',sum5*1.602197e-12
c
      return
      end
      subroutine ispec(tp,xlum,epi,ncn2,zremsz,lpri,lun11)
c
c
c     this subroutine generates the initial spectrum.
c     brems stores the flux to be used
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
      real*8 zremsz(ncn),epi(ncn)
      integer ncn2,lpri,lun11
      real*8 zremsi(ncn),ergsev,const,xlum
      real*8 sum,ekt,tp
      integer i,numcon,lprisv
      real*8 expo
c
      data ergsev/1.602197e-12/
c
      numcon=ncn2
      ekt=1000.*(0.861707)*tp
      sum=0.
      lprisv=lpri
      if (lpri.ge.1) write (lun11,*)'in ispec',tp,xlum
      do i=1,numcon
         zremsi(i)=expo(-epi(i)/ekt)
         if (lpri.gt.1) write (lun11,*)i,epi(i),zremsi(i)
         if (.not.((epi(i).lt.13.6).or.(epi(i).gt.1.36e+4)
     $        .or.(i.le.1)))
     $    sum=sum+(zremsi(i)+zremsi(i-1))*(epi(i)-epi(i-1))/2.
         enddo
c
      const=xlum/sum/ergsev
      do i=1,numcon
         zremsz(i)=zremsi(i)*const
         if (lpri.ge.1)
     $        write (lun11,*)i,epi(i),zremsi(i),const,zremsz(i)
         enddo
      lpri=lprisv
c
      return
      end
      subroutine ispec4(tp,xlum,epi,ncn2,zremsz,lpri,lun11)
c
c     this subroutine generates the initial spectrum.
c     power law spectrum
c     brems stores the flux to be used
c     author:  T. Kallman
c
c
      implicit none
c
      include './PARAM'
c
      real*8 zremsz(ncn),epi(ncn)
      integer ncn2,lpri,lun11
      real*8 zremsi(ncn),ergsev,const,xlum
      real*8 sum,ecut
      integer i,numcon,nb1,nb2,nbinc,lprisv
      real*8 tp
c
      data ergsev/1.602197e-12/
c
      numcon=ncn2
      ecut=0.01
      sum=0.
      lprisv=lpri
      if (lpri.ge.1) write (lun11,*)'in ispec4',tp,xlum
      nb1=nbinc(13.6d0,epi,ncn2)
      nb2=nbinc(1.36d+4,epi,ncn2)
      do i=1,numcon
         zremsi(i)=1.e-24
         if (epi(i).gt.ecut)
     $    zremsi(i)=epi(i)**tp
         if (lpri.gt.1) write (lun11,*)i,epi(i),zremsi(i)
         if ((i.ge.nb1).and.(i.le.nb2))
     $    sum=sum+(zremsi(i)+zremsi(i-1))*(epi(i)-epi(i-1))/2.
         enddo
c
      const=xlum/sum/ergsev
      do i=1,numcon
         zremsz(i)=zremsz(i)+zremsi(i)*const
         if (lpri.ge.1)
     $        write (lun11,*)i,epi(i),zremsi(i),const,zremsz(i)
         enddo
      lpri=lprisv
c
      return
      end
      subroutine ispecg(eptmp,zrtmp,nret,epi,ncn2,zremsz,xlum,
     $                  lpri,lun11)
c
c     this subroutine generates the initial spectrum.
c     brems stores the flux to be used
c     generic renormalization
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
      integer nret
      real*8 zremsz(ncn),epi(ncn)
      integer ncn2,lpri,lun11
      real*8 ergsev,const,xlum
      real*8 sum,tmp,tmpo, exp10
      integer numcon
      integer jlo,kl
      real*8 zremsi(ncn),eptmp(nret),zrtmp(nret)
      real*8 x,epmx,epmn,zr1,zr2,ep1,ep2,alx,aly,y
      integer jk
c
      data ergsev/1.602197e-12/
c
c        linear interpolation in log
      jlo = 0
c     if (lpri.ge.1) write (lun11,*)'in ispecg:',nret
c     if ( lpri.gt.2 ) write (lun11,*) (ll,eptmp(ll),zrtmp(ll),ll=1,nret)
      numcon=ncn2
      do kl = 1,numcon
         x = epi(kl)
         zremsi(kl) = 0.
         epmx = max(eptmp(1),eptmp(nret))
         epmn = min(eptmp(1),eptmp(nret))
         if ( lpri.gt.2 ) write (lun11,*) kl,x,epmx,epmn
         if ( (x.le.epmx) .and. (x.ge.epmn) ) then
            call hunt3(eptmp,nret,x,jlo,lpri,lun11)
            jlo = max0(jlo,1)
            zr1 = log10(max(zrtmp(jlo+1),1.e-24))
            zr2 = log10(max(zrtmp(jlo),1.e-24))
            ep1 = log10(max(eptmp(jlo+1),1.e-24))
            ep2 = log10(max(eptmp(jlo),1.e-24))
            alx = log10(x)
            alx = max(alx,ep2)
            alx = min(alx,ep1)
            aly = (zr1-zr2)*(alx-ep2)/(ep1-ep2+1.e-24) + zr2
            y = exp10(aly)
            zremsi(kl) = y
            if ( lpri.gt.2 ) write (lun11,*) kl,x,jlo,zr1,zr2,
     &                              ep1,ep2,y
         endif
         enddo
c
      sum = 0.
      tmp = zremsi(1)
      if ( lpri.gt.2 ) write (lun11,*) ' in ispecg'
      do jk = 2,ncn2
         tmpo = tmp
         tmp = zremsi(jk)
         if ( lpri.gt.2 ) write (lun11,*) jk,epi(jk),tmp,tmpo,sum
         if ( (epi(jk).ge.13.6) .and. (epi(jk).le.1.36e+4) ) then
            sum = sum + (tmp+tmpo)*(epi(jk)-epi(jk-1))/2.
            endif
         enddo
      sum = sum*ergsev
      const = xlum/sum
      do jk = 1,ncn2
         if ( lpri.gt.2 ) write (lun11,*) jk,epi(jk),zremsz(jk),
     &                                zremsi(jk)
         zremsz(jk) = zremsz(jk) + zremsi(jk)*const
         enddo
c
      return
      end
      subroutine ispecgg(xlum,epi,ncn2,zremsz,
     $               lpri,lun11)
c
c     this subroutine generates the initial spectrum.
c     renormalization
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
      real*8 epi(ncn),zremsz(ncn)
      integer numcon,ncn2,i
      real*8 ergsev,sum,const,xlum
      integer lpri, lun11
c
      data ergsev/1.602197e-12/
c
      numcon=ncn2
      sum=0.
      if (lpri.gt.1) write (lun11,*)'in ispec',xlum
      do i=1,numcon
         if (lpri.gt.1) write (lun11,*)i,epi(i),zremsz(i)
         if ((epi(i).ge.13.6).and.(epi(i).le.1.36e+4)
     $        .and.(i.gt.1))
     $    sum=sum+(zremsz(i)+zremsz(i-1))*(epi(i)-epi(i-1))/2.
         enddo
c
      const=xlum/sum/ergsev
      do i=1,numcon
         zremsz(i)=zremsz(i)*const
         if (lpri.gt.1)
     $        write (lun11,*)i,epi(i),const,zremsz(i)
         enddo
c
      return
      end
      subroutine istruc(zeff,alpha,xitp,nnz,lpri,lun11)
c
      implicit none
c
      integer nnz,lpri,lun11
c
      real*8 zeff(31),alpha(31),xitp(31),xisum
      real*8  z8(31),a8(31),x8(31)
      integer mm,nnzp1,ill
c
      if (lpri.ne.0)
     $ write (lun11,*)'ion rates:',nnz
      do mm=1,nnz
          z8(mm)=dble(zeff(mm))
          a8(mm)=dble(alpha(mm))
          if (lpri.ne.0)
     $     write (lun11,9901)zeff(mm),alpha(mm)
9901      format (1x,2(1pe11.3))
          enddo
c
      nnzp1 = nnz + 1
      ill=1
      call ioneqm(z8,a8,x8,nnzp1,nnz,ill,lpri,lun11)
c
      xisum=0.
      do mm=1,nnz
          xitp(mm)=sngl(x8(mm))
          xisum=xisum+xitp(mm)
          enddo
      xitp(nnz+1)=max(0.,1.-xisum)
c
      return
      end

      integer function lenact(cbuf)
      implicit none
      character cbuf*(*)
c---
c function to return the active length of a character string, not
c counting any trailing blanks.  n.b. an all blank string will
c return zero as the length.
c---
c cbuf    i    string whose length is to be measured.
c---
c 1988-jun-13 - standard fortran version [aft]
c---
      integer   i
c---
      do 190 i=len(cbuf),1,-1
         if(cbuf(i:i).ne.' ') then
            lenact=i
            return
         end if
  190 continue
      lenact=0
      return
      end

      subroutine leqt2f(a,m,n,np,b,idgt,wkarea,ier,lun11,lpri)
c
      implicit none
c
      include './PARAM'
c
      integer indx(nd),ier,lun11,lpri,n,np,m,npp
      real*8 a(np,np),b(np),wkarea(1)
      real*8  ao(ndss,ndss),bo(ndss),btmp,tmp,sum,errmx,err
      real*8  an(ndss,ndss),bn(ndss),d,tmpmx
      integer mm,ll2,mmmx,mmmxo,jk,idgt,kl
c
c     n had better be less than nd
c
c     Not used
      integer javi
      real*8 javir
      javi=m
c      m=javi
      javi=ier
      javir=wkarea(1)
c      wkarea(1)=javir
      javi=idgt
c
      do jk=1,n
        bo(jk)=dble(b(jk))
        bn(jk)=dble(b(jk))
        do kl=1,n
           an(jk,kl)=dble(a(jk,kl))
           ao(jk,kl)=dble(a(jk,kl))
           enddo
        enddo
c
      npp=ndss
      if (lpri.gt.1)
     $ write (lun11,*)'before ludcmp',n,npp,np
      call ludcmp(an,n,npp,indx,d,lun11,lpri)
      npp=ndss
      if (lpri.gt.1)
     $ write (lun11,*)'after ludcmp',n,npp
      call lubksb(an,n,npp,indx,bn,lun11,lpri)
      if (lpri.gt.1)
     $ write (lun11,*)'after lubksb'
      npp=ndss
      call mprove(ao,an,n,npp,indx,bo,bn,lun11,lpri)
      if (lpri.gt.2)
     $ write (lun11,*)'after mprove',n,npp,np
c
c        check the solution
         if (lpri.gt.2) write (lun11,*)'checking the solution'
         errmx=0.d0
         do  ll2=1,n
          sum=0.d0
          tmpmx=0.d0
          mmmx=0
          mmmxo=0
          do  mm=1,n
            btmp=bn(mm)
            tmp=dble(a(ll2,mm))*max(0.d0,btmp)
            if (abs(tmp).ge.tmpmx) then
              mmmxo=mmmx
              mmmx=mm
              tmpmx=max(tmpmx,abs(tmp))
              endif
            sum=sum+tmp
            enddo
          sum=sum-dble(b(ll2))
          err=sum/max(1.d-24,tmpmx)
          errmx=max(errmx,abs(err))
          if (lpri.gt.2) write (lun11,9246)ll2,bn(ll2),tmpmx,sum,err,
     $                                     mmmx,mmmxo
 9246     format (1h ,i4,4e12.4,2i4)
          enddo
c
      do jk=1,n
         if (lpri.gt.2)
     $    write (lun11,*)jk,b(jk)
         if (bn(jk).lt.1.d-36) bn(jk)=0.d0
         if (bn(jk).gt.1.d+36) bn(jk)=1.d+36
         b(jk)=sngl(bn(jk))
         enddo
c
         if (lpri.gt.2)
     $    write (lun11,*)'leaving leqt'
c      ier=0
c      wkarea(1)=0.
c      idgt=0
c
      return
      end
      subroutine levwk(rniss,bb,lpri,rlev,ilev,
     $          nlpt,iltp,nlev,klev,t,xee,xpx,lun11)
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
      character(1) klev(100,nd)
      real*8 rniss(nnml)
      real*8 rlev(10,nd)
      integer ilev(10,nd),nlpt(nd),iltp(nd)
      real*8 ergsev,bk,t,bktm,q2,rs,ethion,emltlv,
     $     eexlv,ethsht,explev2,bb,expo
      integer lpri,lprisv,nlev,lun11,ll
c
      real*8 xnx, xpx, xee, tm
      integer mm

      data ergsev/1.602197e-12/
      data bk/1.38062e-16/
c
      lprisv=lpri
c      lpri=0
      xnx=xpx*xee
      bb=1.
      tm=t*1.e4
      bktm=bk*tm/ergsev
      q2=2.07e-16*xnx*(tm**(-1.5))
      emltlv=rlev(2,nlev)
      rs=q2/emltlv
      ethion=rlev(1,nlev)
      if (lpri.gt.1)
     $ write (lun11,9902)tm,bktm,q2,
     $    emltlv,rs,ethion,xnx
 9902 format (1x,'in levwk',8(1pe11.3))
      rniss(nlev)=1.
      do ll=1,nlev-1
        eexlv=rlev(1,ll)
        emltlv=rlev(2,ll)
        ethsht=(ethion-eexlv)/bktm
        ethsht=max(ethsht,0.)
        explev2=expo(-ethsht)
        rniss(ll)=emltlv/(explev2/rs)
        bb=bb+rniss(ll)
        if (lpri.gt.1)
     $   write (lun11,9901)ll,eexlv,emltlv,ethsht,explev2,
     $     rniss(ll),rs,bb,ilev(1,ll),iltp(ll),nlpt(ll),
     $     (klev(mm,ll),mm=1,8)
 9901   format (1x,i4,7(1pe11.3),3i6,8a1)
        enddo
        do ll=1,nlev
          rniss(ll)=rniss(ll)/bb
          enddo
       bb=1.
       lpri=lprisv
c
      return
      end
      subroutine linopac(lprie,lun11,optpp,ans2,sigvtherm,vtherm,bremsa,
     $                   rcem1,rcem2,elin,vturbi,t,aatmp,delea,epi,ncn2,
     $                   opakc,opakscatt,rccemis,fline,lfast)
c
      implicit none
c
c     this routine puts line opacity into continuum bins
c     author:  T. Kallman
c
c
      include './PARAM'
      integer nbtpp
      parameter (nbtpp=20000)
c
      real*8 epi(ncn),opakc(ncn),dpthc(2,ncn),opakscatt(ncn)
     $       ,bremsa(ncn)
      integer ldon(2)
      real*8  rccemis(2,ncn),fline(2,nnnl)
c
      real*8  prftmp,sum,rcem1,rcem2,ans2,sigvtherm,vtherm
c
      real*8 vturbi,delr,t,eliml,elimh,
     $  dpcrit,bbb,optpp,delea,aatmp,elin,etmp,vth,
     $  vturb,deleturb,deleth,dele,aasmall,
     $  deleused,deleepi,delet,deletpp,e00,dpthmx,
     $  dpthtmp,e0,tst,opsum,optmpo,profile,optp2,
     $  tmpopmx,tmpopo,etptst,opsv4,sum2,optmp2,optmp2o
      integer ldirt,lpri,lun11,
     $  ncn2,llk,lnn,ln,ml,lup,nilin,nelin,
     $  nbtmp,iion,nitmp,ndtmp,mllz,
     $  iltmp,i,ij,lind,
     $  lcon,ldir,mlm,ml1m,ltyp,lrtyp,ml1,ml2,ml1min,
     $  ml1max,mlc,mloff,mlmin,mlmax,ncut,nidt,
     $  nkdt,nrdt,np2,ncsvn,lprie
      integer nbinc,mlpar,lfast
      real*8 voigte
      integer np1i,np1r,np1k
c
c     temporary grid for use in calculating profile
      real*8 etpp(nbtpp),optpp2(nbtpp),optmp(ncn)
c
c      real*8  tmpew,tmpewo,tmpop,tmpe,sum,sume
      real*8 tmpew,tmpewo,tmpop,tmpe,sume,rnormchk,ergsev
c
      data dpcrit/1.e-6/,ergsev/1.602197e-12/
c
c
      lpri=lprie
c      lpri=0
c
c     test whether line is in range
      if ((elin.gt.1.e+8).or.(elin.lt.1.)) return
c
c     for scattering model, add in line opacity
      bbb=vturbi
      elin=abs(elin)
c     thermal width quantities
      vth=(1.29E+1)*sqrt(t/aatmp)
      vturb=bbb
c      e0=(12398.41)/max(elin,1.E-24)
      e0=(12398.41)/max(elin,1.E-24)
      if (e0.le.epi(1)) return
      deleturb=e0*(vturb/3.E+5)
      deleth=e0*(vth/3.E+5)
c     old expression
c     dele=deleth+deleturb
c     new expression
      dele=sqrt(deleth*deleth+deleturb*deleturb)
      aasmall=delea/(1.E-24+dele)/12.56
c
c     continuum bin for line
      ml1=nbinc(e0,epi,ncn2)
      ml1=max(min(ncn-1,ml1),2)
c
c     here is what we do to get the heating right
      prftmp=2./(epi(ml1+1)-epi(ml1-1))
      opsv4=optpp*dele
c
c     print line quantities
      if (lpri.ge.1) write (lun11,*)
     &   'e0,optpp,dpcrit*opakc(ml1),ml1,deleth,delea:',
     &    e0,optpp, dpcrit*opakc(ml1),ml1,deleth,delea
      if (lpri.ge.1) write (lun11,*)optpp,prftmp,opsv4,dele,
     $       opsv4*prftmp,rcem1,rcem2
c
c     test for simple calculation
      if (lfast.gt.2) then
c
c         single bin calculation
          opakc(ml1)=opakc(ml1)+opsv4*prftmp
          opakscatt(ml1)=opakscatt(ml1)+opsv4*prftmp
          rccemis(1,ml1)=rccemis(1,ml1)+rcem1*prftmp/ergsev/12.56
          rccemis(2,ml1)=rccemis(2,ml1)+rcem2*prftmp/ergsev/12.56
          return
c
c       full profile calculation
        else
c
c         calculate profile on temporary grid
c         set up temporary grid
          e00=epi(ml1)
          etmp=e0
c         deleepi is the grid spacing of the epi grid
c         deletpp is the physical energy spacing needed
c           for an accurate integration of the voigt profile
c         ncut is the ratio of these two quantities,
c           used for rebinning the calculated voigt profile
          deleepi=epi(ml1+1)-epi(ml1)
c         expanding step to make broader lines
          deletpp=dele
          ncut=int(deleepi/deletpp)
          ncut=max(ncut,1)
          ncut=min(ncut,nbtpp/10)
          deleused=deleepi/float(ncut)
          mlc=0
          ldir=1
          ldon(1)=0
          ldon(2)=0
          mlmin=nbtpp
          mlmax=1
          ml1min=ncn+1
          ml1max=0
          ml2=nbtpp/2
          if (lpri.ge.1) write (lun11,*)'ncut=',ncut,deleused,deletpp,
     $                                  deleepi
c
c         calculate profile at continuum bin closest to line center
          delet=(e00-etmp)/dele
          if (aasmall.gt.1.e-6) then
              profile=voigte(abs(delet),aasmall)/1.772
            else
              profile=exp(-delet*delet)/1.772
            endif
          etpp(ml2)=e00
          optpp2(ml2)=optpp*profile
          tst=1.
c
c         now put profile on temporary grid
c         work outward in both directions from line center
          do while ((ldon(1)*ldon(2).eq.0).and.(mlc.lt.nbtpp/2))
c
            mlc=mlc+1
c
c           alternate directions
            do ij=1,2
              ldir=-ldir
c
c             test to see if done in this direction
              if (ldon(ij).ne.1) then
c
c               index into temporary grid
                mlm=ml2+ldir*mlc
c
c               energy of temporary grid point
                etptst=e00+float(ldir*mlc)*deleused
c
c               test to see if within allowed range
                if ((mlm.le.nbtpp).and.(mlm.ge.1)
     $           .and.(etptst.gt.0.).and.(etptst.lt.epi(ncn2))) then
c
c                 calculate index extremes for later use
c                 ml1m is index into epi grid
c                 ml1min and ml1max are extremes of ml1m
c                 mlmin and mlmax are extremes of mlm
                  mlmin=min(mlm,mlmin)
                  mlmax=max(mlm,mlmax)
c
c                 store energy binc
                  etpp(mlm)=e00+float(ldir*mlc)*deleused

c                 calculate profile
                  delet=(etpp(mlm)-etmp)/dele
                  if (aasmall.gt.1.e-9) then
                      profile=voigte(abs(delet),aasmall)/1.772
                    else
                      profile=exp(-delet*delet)/1.772
                    endif
c
c                 calculate opacity
                  optpp2(mlm)=optpp*profile
c                  tst=optpp2(mlm)*delr
                  tst=profile
c
c                 print
                  if (lpri.ge.1) write (lun11,*) 'first write',
     $             mlm,etpp(mlm),ij,
     $             deleused,delet,mlmin,mlmax,ml1,
     $             mlc,mloff,mod(mloff,ncut),profile,optpp2(mlm),
     $             tst
c
c                 end of test for within range
                  endif
c
c               test to see if done in this direction:
c                 profile not too small
c                 index within range
c                 energy within range
c                 within specified number of doppler widths (50)
                if (((tst.lt.dpcrit)
     $               .or.(mlm.le.1).or.(mlm.ge.nbtpp)
     $               .or.(etptst.le.0.).or.(etptst.ge.epi(ncn2))
     $               .or.(mlc.gt.nbtpp)
     $               .or.(abs(delet).gt.max(50.,200.*aasmall)))
     $               .and.(ml1min.lt.ml1-2).and.(ml1max.gt.ml1+2)
     $               .and.(ml1min.ge.1).and.(ml1max.le.ncn))
     $                ldon(ij)=1
c
c               end of test for done in this direction
                endif
c
c             end of loop over directions
              enddo
c
c           end of loop over energies
            enddo
c
c         store into continuum bins
          sum=0.
          opsum=0.
          tmpop=0.
          tmpopmx=0.
          sume=0.
          ml1min=nbinc(etpp(mlmin),epi,ncn2)
          ml1max=nbinc(etpp(mlmax),epi,ncn2)
          ml1m=ml1min
          if (lpri.ge.1) write (lun11,*)'renormalizing profile',
     $       ml2,mlmin,mlmax,ml1m,ml1min,ml1max
          tmpew=0.
          mlmin=max(mlmin,2)
          mlmax=min(mlmax,nbtpp)
c
c         step through temp grid bins
c         and  sum over intervals
          do mlm=mlmin+1,mlmax
c
            tmpopo=tmpop
            tmpop=optpp2(mlm)
            tmpopmx=max(tmpopmx,tmpop)
            tmpe=abs(etpp(mlm)-etpp(mlm-1))
c
c           update interval sum
            sume=sume+tmpe
            opsum=opsum+(tmpop+tmpopo)*tmpe/2.
c
c           test to see if you have reached epi grid boundary
            if (etpp(mlm).gt.epi(ml1m)) then
c
c             store current sum
              optmpo=opakc(ml1m)
              if (sume.gt.1.d-34) then
                optp2=opsum/sume
               do while ((etpp(mlm).gt.epi(ml1m)).and.(ml1m.lt.ncn2))
                  opakc(ml1m)=opakc(ml1m)+optp2
c                 print
                  if (lpri.ge.1) write (lun11,*)mlm,ml1m,
     $             epi(ml1m),epi(ml1m+1),etpp(mlm),opakc(ml1m),
     $               optmpo,optpp2(mlm),opsum,sume
                  ml1m=ml1m+1
                  enddo
                endif
c
c             reset interval sums
              tmpopmx=0.
              opsum=0.
              sume=0.
c
c             end of test for epi bin boundary
              endif
c
c           end of rebinning loop
            enddo
 9000     continue
c
c         norm check
c          rnormchk=ewsv(nlsv)/optpp/dele/(1.e-34+delr)
c          if (lpri.ne.0) write (lun11,*)'norm check',nilin,elin,optpp,
c     $          dele,aasmall,ewsv(nlsv),rnormchk
c
c       end of test for fast calculation
        endif
c
c
      return
      end
      subroutine lubksb(a,n,np,indx,b,lun11,lpri)
      implicit none
c
      integer i, ii, j, ll, n, np
      integer lpri, lun11
c
      integer indx(np)
      real*8  a(np,np), b(np), sum
c
      if (lpri.gt.1)
     $ write (lun11,*)'in lubksb',n,np
      ii = 0
      do 100 i = 1 , n
         ll = indx(i)
c         write (lun11,*)'i,ll:',i,ll
         sum = b(ll)
         b(ll) = b(i)
         if ( ii.ne.0 ) then
            do 20 j = ii , i - 1
               if (lpri.gt.1)
     $          write (lun11,*)'i,j,ii:',i,j,ii,sum
               sum = sum - a(i,j)*b(j)
 20         continue
         elseif ( sum.ne.0.d0 ) then
            ii = i
         endif
         if (lpri.gt.1)
     $    write (lun11,*)'i,sum:',i,sum
         b(i) = sum
 100  continue
      do 200 i = n , 1 , -1
         sum = b(i)
         if ( i.lt.n ) then
            do 120 j = i + 1 , n
               sum = sum - a(i,j)*b(j)
 120        continue
         endif
         b(i) = sum/a(i,i)
 200  continue
      return
      end
      subroutine ludcmp(a,n,np,indx,d,lun11,lpri)
c
c
      implicit none
c
      real*8  tiny
      parameter (tiny=1.0d-20)
      include './PARAM'
      integer np, n
c
c
      real*8  a(np,np),vv(nd),d,aamax,sum,dum
      integer indx(n),lun11,lpri,j,i,imax,k
c
      if (lpri.gt.1)
     $ write (lun11,*)'in ludcmp:'
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
          if (lpri.gt.1) write (lun11,*)i,j,a(i,j),aamax
11      continue
        if (aamax.eq.0.d0) then
           if (lpri.gt.1)
     $      write (lun11,9902)
           return
9902       format (1h ,'singular matrix.' )
         end if
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.d0
        imax=0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        imax=max(imax,1)
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        if (lpri.gt.1)
     $   write (lun11,*)'j,imax:',j,imax
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.d0)a(j,j)=tiny
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.0.d0)a(n,n)=tiny
c
      return
      end
      subroutine milne(temp4,nt,x4,y4,eth4,alpha4,lun11,lpri)
C
c     this routine calculates the Milne relation for type 53 data
c       x    = array of energies in Ry with respect to the threshold
c       y    = array of cross sections in Mb
c       nt   = number of points in topbase arrays
c       eth  = threshold energy in Ry
c       alpha= recombination rate for level n,lo
c     author:  M. Bautista
c
      implicit none
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
      integer lun11,lpri,nt,i
      real*8 x4(nt),y4(nt),temp4,alpha4,eth4
      real*8   temp,x(nptmpdim),y(nptmpdim),eth,alpha,st,ry,
     $    sum,s1,s2,v1,v2,rb,ra,ri2,ri3,sumo,crit
c
c      lpri=2
c
      ry=2.17896d-11
c
       temp=(temp4)
       eth=(eth4)
       do i=1,nt
         x(i)=(x4(i))
         y(i)=(y4(i))
         enddo
c
       st=(x(1)+eth)*ry
       sumo=1.d0
       sum=0.d0
       if (lpri.gt.1)
     $   write (lun11,*)'in milne:',temp,nt,eth,x(1),y(1)
       i=1
       crit=0.01d0
       do while ((abs(sum-sumo).gt.crit*sum).and.(i.lt.nt))
         i=i+1
         s1=(x(i-1)+eth)*ry
         s2=(x(i)+eth)*ry
         if (s2.lt.s1) return
         v1=y(i-1)
         v2=y(i)
         if (lpri.gt.1) write (lun11,*)'i=',i,x(i),y(i),
     $     s1,s2,v1,v2
         if ((v1.ne.0.d0).or.(v2.ne.0.d0)) then
           rb=(v2-v1)/(s2-s1+1.d-24)
           ra=v2-rb*s2
           call intin(s1,s2,st,temp,ri2,ri3,lpri,lun11)
           sumo=sum
           sum = sum + (ra*ri2 + rb*ri3)
           if (lpri.gt.1) write (lun11,*)i,x(i),y(i),
     $     s1,s2,v1,v2,ra,rb,ri2,ri3,sum
           endif
         enddo
c      alpha=sum*.79788*2.4917e+25/(temp**1.5)
       alpha=sum*.79788d0*40.4153d0
       if (lpri.gt.1)
     $    write (lun11,*)'alpha=',alpha,sum
       alpha4=sngl(alpha)
c
      return
      end
      subroutine mprove(a,alud,n,np,indx,b,x,lun11,lpri)
      implicit none
      integer np, n,nl
      parameter (nl=10000)
      real*8  a(np,np), alud(np,np), b(n), x(n), r(nl), sdp
      integer indx(n), lun11, lpri, i, j
      if (lpri.gt.1)
     $ write (lun11,*)'in mprove:'
      do 12 i=1,n
        sdp=-b(i)
        do 11 j=1,n
          sdp=sdp+a(i,j)*x(j)
11      continue
        r(i)=sdp
12    continue
      call lubksb(alud,n,np,indx,r,lun11,lpri)
      do 13 i=1,n
        x(i)=x(i)-r(i)
13    continue
c
      return
      end
      subroutine msolvelucy(ajisb,cjisb,indb,nindb,nsup,nspmx,
     $   ipmat,bmat,x,ht,cl,niter,nit2,nit3,nitmx,nitmx2,lun11,lpri)
c
c     solves lucy iteration
c     author:  T. Kallman
c
      implicit none
c
c     solves lucy iteration
c     author:  T. Kallman
c
      include './PARAM'
c
c
      real*8 ajisb(2,ndb),cjisb(ndb)
      integer indb(2,ndb)
      real*8 x(nd),bmat(nd),xo(nd),xoo(nd)
      real*8 bmatsup(ndss),ajissup(ndss,ndss),xsup(ndss),rr(nd)
      real*8 cjissup(ndss,ndss)
      real*8 wkarea(1)
      real*8 tt1, crit, crit2, diff, tt2, diff2
      real*8 riu(nds),ril(nds),rui(nds),rli(nds),xsum,tst,cl
      real*8 ht, clp, htp
      integer nsup(nd), mm, ipmat, ll
      integer lpril, lpri, lprisv, idgt, ier
      integer lun11, niter, nit3, nitmx, nspmx
      integer nn, nsp, ngood, nspm, nspn, nspcon
      integer nit2, nitmx2, m2, nindb
c
c
c       step thru levels, and form calculate superlevel quantities
      lprisv=lpri
      call remtms(tt1)
c      lpri=0
      if (lpri.gt.1)
     $ write (lun11,*)'in msolvelucy',lpri,nindb
c      crit=1.e-2
c      crit2=1.e-2
      crit=1.e-4
      crit2=1.e-4
      diff=1.
      niter=0
      nit3=0
      do while ((diff.gt.crit).and.(niter.lt.nitmx))
        niter=niter+1
        if (lpri.gt.1) write (lun11,*)'iteration=',niter
        if (lpri.gt.1) write (lun11,*)'initial populations:'
        do mm=1,ipmat
          xo(mm)=x(mm)
          if (lpri.gt.3) write (lun11,*)mm,x(mm),nsup(mm)
          enddo
        do mm=1,nspmx
          xsup(mm)=0.
          bmatsup(mm)=0.
          do nn=1,nspmx
            ajissup(mm,nn)=0.
            cjissup(mm,nn)=0.
            enddo
          enddo
        do mm=1,ipmat
          nsp=nsup(mm)
          xsup(nsp)=xsup(nsp)+x(mm)
          enddo
        call remtms(tt2)
        if (lpri.gt.1)
     $    write (lun11,*)'before constucting matrix',abs(tt2-tt1)
        tt1=tt2
        if (lpri.gt.1)
     $      write (lun11,*)'constucting the condensed matirx:'
        ngood=0
        do mm=1,ipmat
          nspm=nsup(mm)
          rr(mm)=x(mm)/(1.e-36+xsup(nspm))
          if (xsup(nspm).le.1.e-36) rr(mm)=1.
          enddo
        do ll=1,nindb
          mm=min(ipmat,indb(1,ll))
          nn=min(ipmat,indb(2,ll))
          nspm=nsup(mm)
          nspn=nsup(nn)
          if ((nspn.ne.nspm)
     $         .and.(nspn.ne.0).and.(nspm.ne.0)
     $         .and.((abs(rr(mm)).gt.1.e-36)
     $           .or.(abs(rr(nn)).gt.1.e-36))
     $         .and.((abs(ajisb(1,ll)).gt.1.e-36)
     $           .or.(abs(ajisb(2,ll)).gt.1.e-36))) then
              ajissup(nspm,nspn)=ajissup(nspm,nspn)
     $            +(ajisb(1,ll))*rr(nn)
              ajissup(nspm,nspm)=ajissup(nspm,nspm)
     $            -(ajisb(2,ll))*rr(mm)
c              ajissup(nspn,nspm)=ajissup(nspn,nspm)
c     $            +(ajisb(2,ll))*rr(nn)
c              ajissup(nspn,nspn)=ajissup(nspn,nspn)
c     $            -(ajisb(1,ll))*rr(mm)
              cjissup(nspm,nspn)=cjissup(nspm,nspn)
     $            +(cjisb(ll))*rr(mm)
              ngood=ngood+1
              if (lpri.gt.3) write (lun11,91)ll,mm,nn,nspm,nspn,
     $              rr(mm),rr(nn),ajisb(1,ll),ajisb(2,ll),
     $              ajissup(nspm,nspn),ajissup(nspm,nspm)
91            format (1x,'used ',5i6,6(1pe13.5))
            endif
          enddo
        call remtms(tt2)
        if (lpri.gt.1)
     $     write (lun11,*)'after constucting matrix',abs(tt2-tt1),
     $                       ngood
        if (lpri.gt.1) then
          if (lpri.gt.1) write (lun11,*)'the condensed populations:'
          do nsp=1,nspmx
            write (lun11,*)nsp,xsup(nsp)
            enddo
          if (lpri.gt.1) write (lun11,*)'the condensed matrix:'
          do nspm=1,nspmx
            do nspn=1,nspmx
              if (abs(ajissup(nspm,nspn)).gt.1.e-37)
     $         write (lun11,*)nspm,nspn,ajissup(nspm,nspn)
              enddo
            enddo
          endif
c        put in number conservation
c         nspcon=1
         nspcon=nspmx
         do mm=1,nspmx
           ajissup(nspcon,mm)=1.
           bmatsup(mm)=0.
           enddo
        bmatsup(nspcon)=1.
        lpril=0
        call remtms(tt1)
        if (lpri.gt.2)
     $    write (lun11,*)'before leqt',abs(tt2-tt1)
        call leqt2f(ajissup,1,nspmx,ndss,bmatsup,idgt,wkarea,ier,
     $                      lun11,lpril)
         call remtms(tt2)
         if (lpri.gt.2)
     $    write (lun11,*)'after leqt',abs(tt2-tt1)
        if (lpri.gt.2) write (lun11,*)'the new condensed populations:'
        do mm=1,nspmx
          xsup(mm)=bmatsup(mm)
          if (lpri.gt.2) write (lun11,*)mm,xsup(mm)
          enddo
        if (lpri.gt.3) write (lun11,*)'new populations'
        do mm=1,ipmat
          nsp=nsup(mm)
          x(mm)=rr(mm)*xsup(nsp)
          if (lpri.gt.3) write (lun11,*)mm,nsp,rr(mm),x(mm)
          enddo
        nit2=0
        diff2=10.
        do while ((nit2.lt.nitmx2).and.(diff2.ge.crit2))
          nit2=nit2+1
          nit3=nit3+1
          if (lpri.gt.2) write (lun11,*)'before calculate new x(mm)',
     $                                   nit2,nit3
          call remtms(tt2)
          if (lpri.gt.2)
     $    write (lun11,*)'in diff2 loop',abs(tt2-tt1)
          tt1=tt2
          do mm=1,ipmat
            riu(mm)=0.
            rui(mm)=0.
            ril(mm)=0.
            rli(mm)=0.
            enddo
          if (lpri.gt.3) write (lun11,*)'the riu calculation'
          do ll=1,nindb
            mm=indb(1,ll)
            nn=min(ipmat,indb(2,ll))
            if (nn.gt.mm) then
                riu(mm)=riu(mm)+abs(ajisb(2,ll))
                rui(mm)=rui(mm)+abs(ajisb(1,ll))*x(nn)
                if (lpri.gt.3) write (lun11,*)mm,nn,ajisb(2,ll),
     $           ajisb(1,ll),x(nn),riu(mm),rui(mm)
              endif
            enddo
          if (lpri.gt.3) write (lun11,*)'the ril calculation'
          do ll=1,nindb
            mm=indb(1,ll)
            nn=min(ipmat,indb(2,ll))
            if (nn.lt.mm) then
c               I hope the indeces are in the right order here
                ril(mm)=ril(mm)+abs(ajisb(2,ll))
                rli(mm)=rli(mm)+abs(ajisb(1,ll))*x(nn)
                if (lpri.gt.3) write (lun11,*)mm,nn,ajisb(2,ll),
     $           ajisb(1,ll),x(nn),ril(mm),rli(mm)
              endif
            enddo
          do mm=1,ipmat
            xoo(mm)=x(mm)
            x(mm)=(rli(mm)+rui(mm))/(ril(mm)+riu(mm)+1.e-24)
            if (lpri.gt.3) write (lun11,*)mm,riu(mm),rui(mm),ril(mm),
     $                                     rli(mm),x(mm),xoo(mm)
            enddo
          xsum=0.
          do mm=1,ipmat
            xsum=xsum+x(mm)
            enddo
          if (lpri.gt.3) write (lun11,*)'new and old populations',
     $                      xsum
          do mm=1,ipmat
            x(mm)=x(mm)/(1.e-24+xsum)
            enddo
          m2=1
          diff2=0.
          tst=0.
          do while ((diff2.lt.1.e+3)
     $           .and.(m2.le.ipmat).and.(tst.lt.1.e+3))
            if (lpri.gt.3) write (lun11,*)m2,x(m2),xoo(m2),xo(m2),
     $            diff2
            tst=1.
            if (x(m2).gt.1.e-22) tst=xoo(m2)/x(m2)
            diff2=diff2+(tst-1.)*(tst-1.)
            m2=m2+1
            enddo
          if (lpri.gt.2) write (lun11,*) 'diff2=',diff2,nit2
          enddo
        diff=0.
        m2=1
        do while ((m2.le.ipmat).and.(diff.lt.1.e+3))
          tst=xo(m2)*(1.e-30)
          if ((x(m2).gt.tst).and.(diff.lt.1.e+10).and.(x(m2).gt.1.e-35))
     $     diff=diff+(min(1.e+10,(xo(m2)-x(m2))/(xo(m2)+x(m2))))**2
          if (lpri.gt.3) write (lun11,*)m2,x(m2),xo(m2),
     $            diff
          m2=m2+1
          enddo
        if (lpri.gt.2) write (lun11,*) 'diff=',diff
      enddo
c
      if (lpri.ge.1) write (lun11,*)'heating-cooling in msolvelucy:'
      cl=0.
      ht=0.
      do ll=1,nindb
        mm=min(indb(1,ll),ipmat)
        nn=min(indb(2,ll),ipmat)
        if (cjisb(ll).gt.0.) then
              cl=cl+x(mm)*cjisb(ll)
            else
              ht=ht-x(mm)*cjisb(ll)
            endif
          if ((lpri.ge.1).and.(abs(cjisb(ll)).gt.1.e-24))
     $         write (lun11,981)ll,mm,nn,x(mm),cjisb(ll),ht,cl
 981           format (1x,3i6,4(1pe11.3))
        enddo
      go to 9090
      if (lpri.gt.2) write (lun11,*)'heating-cooling superlevels:'
      clp=0.
      htp=0.
      do mm=1,nspmx
        do nn=1,nspmx
          if (cjissup(mm,nn).gt.0.) then
              clp=clp+xsup(mm)*cjissup(mm,nn)
            else
              htp=htp-xsup(mm)*cjissup(mm,nn)
            endif
          if ((lpri.gt.2).and.(abs(cjissup(mm,nn)).gt.1.e-24))
     $         write (lun11,*)mm,nn,xsup(mm),cjissup(mm,nn),htp,clp
          enddo
        enddo
      ht=htp
      cl=clp
 9090 continue
c
      lpri=lprisv
c
      return
      end
      function nbinc(e,epi,ncn2)
      implicit none
c
c     this function bins the continuum
c     lines between   epi(i) and   epi(i+1) are put in bin number i.
c     energies between 0 and   epi(1) are put in bin number 50.
c     author:  T. Kallman
c
      include './PARAM'
c
      real*8 e
      integer jlo, lun11,nbinc, ncn2, numcon, numcon2,
     &        numcon3
      real*8 epi(ncn)
c
      lun11=6
      numcon=ncn2
      numcon2=max0(2,ncn2/50)
      numcon3=numcon-numcon2
      call huntf(epi,numcon3,e,jlo,0,lun11)
c     call hunt3(epi,numcon3,e,jlo,0,lun11)
c      if (abs(e-epi(jlo+1)).lt.abs(e-epi(jlo))) jlo=jlo+1
      nbinc=jlo
c
      return
      end
      real*8 function pescl(tau)
c
c     this routine calculates escape probability for a line transition
c     inputs: optical depths-- tau for a line transition
c
      implicit none
c
      real*8 tau
      real*8 pi,tauw,aa,bb
c
      data pi/3.1415927/
c
      tauw=1.e5
c     *** need to determine tauw from line profiles?***
      if(tau.lt.1.0) then
        if(tau.lt.1.e-5) then
          pescl=1.0
          go to 10
        end if
        aa=2.0*tau
        pescl=(1.0-exp(-aa))/aa
        go to 10
      end if
      bb=0.5*sqrt(log(tau))/(1.0+tau/tauw)
      pescl=1./(tau*sqrt(pi)*(1.2+bb))
10    continue
      pescl=pescl/2.0
c
      return
      end
      real*8 function pescv(tau)
c
c     this routine calculates escape probability for
c       continuum
c     inputs: optical depths-- tau, energy
c     author:  T. Kallman
c
      implicit none
c
      real*8 tau
      real*8 taubar,eps
c
c     NB optically thin rrcs
c     test for lte
      pescv=0.5
      return
c      if (e.lt.13.6) return
c     fudge because of too much case b
c      taubar=tau/200.
c      taubar=tau/5.
      taubar=tau*100.
      pescv=1./(taubar+1.)
c      pescv=expo(-taubar)
c     fudge because of numerical problem with large tau.
      eps=1.e-6
      pescv=max(pescv,eps)
      pescv=pescv/2.
c
      return
      end

      subroutine pexs(nmin,kdim,zc,eion,far,gam,scal,
     +                e,axs,ierr,lpri,lun11)
c
c     Compute photoexcitation cross-section assuming one
c     Rydberg serie converging to a threshold.
c
c     nmin = starting princ. quant. num. of the serie
c     zc = effective charge, Z-Ne+1
c     eion = threshold energy in Ry
c     far = oscillator strength of the serie member with
c            n=nmin
c     gam = resonance width in Ry
c     e = external energy grid in Ry
c     kdim = dimension of the external energy grid
c     axs = cross section
c     scal = scaling factor
c     ierr = error indicator (=0: OK; .ne.0:error)
c     author:  P. Palmeri
c
      implicit none
c
      integer nmax,nmin
      real*8  pi
c
      parameter(pi=3.14159d0,nmax=30)
c
      real*8  x(nmax),a(nmax),e(*),axs(*)
      real*8  zc,eion,far,gam,scal
      integer lpri,ierr,lun11,kdim,im
      integer nres,jmin,jmax,jres,jj,i,ii,ij,kk,n
      real*8  del,xmin,xres,axtp,res
c
      data x,a/nmax*0.,nmax*0./
c
      if (lpri.ne.0) write (lun11,*)'in pexs',nmin,kdim,zc,
     $        eion,far,gam,scal
c
      ierr=0
      if(nmin.ge.nmax) then
       ierr=1
       return
      endif
      nres=nmax
      do 10 n=nmin,nmax
c
c   energy of the resonance ...
c
       x(n)=-(zc/dble(n))**2.d0
c
c   area of the resonance ...
c
       a(n)=8.06725d0*far*dble(nmin**3)/dble(n**3)
c
c   search for unresolved limit in term of member ...
c
       if(n.gt.nmin) then
        del=x(n)-x(n-1)
         res=gam/2.d0
        if(del.gt.res) nres=n
       endif

        if (lpri.ne.0) write (lun11,*)n,x(n),a(n),del,res

   10 continue
c
c   define shifted energy range ...
c
      xmin=x(nmin)-30.d0*gam
      xres=x(nres)
      jmin=1
      jmax=kdim
      jres=jmax
      do 20 i=1,kdim
        axs(i)=0.d0
       e(i)=e(i)-eion
       im=max(1,i-1)
       if(i.gt.1.and.e(im).le.xmin.and.e(i).gt.xmin)
     +    jmin=i-1
       if(i.gt.1.and.e(im).le.xres.and.e(i).gt.xres)
     +    jres=i-1
       if(i.gt.1.and.e(im).lt.0.d0.and.e(i).ge.0.d0)
     +    jmax=i-1
   20 continue
      if(jmin.eq.jmax) jmax=jmin+1
      if (lpri.ne.0) write (lun11,*)'jmin,jres:',jmin,jres,
     $              jmax,nmin,nmax
      do 30 ii=nmin,nmax
       do 30 jj=jmin,jres
c
c   constant-width Lorentzian resonances ...
c
        axtp=a(ii)/pi*gam/2.d0
     +  /((e(jj)-x(ii))**2+(gam/2.d0)**2)
     +  +a(ii)/pi*gam/2.d0/((abs(e(jj))-x(ii))**2
     +   +(gam/2.d0)**2)
c
c   near-threshold pill-up (oscill. strength conservation)
c
        axs(jj)=axs(jj)+axtp
      if ((lpri.ne.0).and.(axs(jj).gt.1.d-24))
     $  write (lun11,*)'30 loop',jj,e(jj),axtp,axs(jj)
   30 continue
c
c   near-threshold extrapolation ...
c
      do ij=jres+1,jmax
        axs(ij)=axs(jres)
        if (lpri.ne.0) write (lun11,*)ij,axs(ij)
      enddo
c
c   scaling of the xs ...
c
      do 40 kk=1,kdim
       axs(kk)=scal*axs(kk)
c
c   return to the "usual" energy grid ...
c
       e(kk)=e(kk)+eion
       if ((lpri.ne.0).and.(axs(kk).gt.1.d-24))
     $  write (lun11,*)kk,e(kk),axs(kk)
   40 continue
c
c
      return
      end
      subroutine phextrap(etmp,stmp,ntmp,ntmp2,ett,ncn2,lpri,lun11)
c
c     this routine does the extrapolation of the ip photoionization cross
c     sections
c     called by ucalc, for data types 53, 49.
c     author: T. Kallman
c
      implicit none
c
      integer ntmp,ntmp2
c
      real*8 etmp(ntmp2),stmp(ntmp2)
      real*8 ett
      integer ncn2,lpri,lun11
      real*8 dele,dels,s1,e1,e2,s2
      integer nadd
c
      dele=1.3
      nadd=0
      dels=dele**3
      s1=stmp(ntmp-1)
      e1=etmp(ntmp-1)*13.6+ett
      if (lpri.gt.1) write (lun11,*)'in phextrap:',ntmp,e1,s1
      do while ((s1.gt.1.e-27).and.(nadd+ntmp.lt.ncn2)
     $         .and.(e1.lt.2.e+4))
        e2=e1*dele
        s2=s1/dels
        nadd=nadd+1
        stmp(nadd+ntmp-1)=s2
        etmp(nadd+ntmp-1)=(e2-ett)/13.6
        e1=e2
        s1=s2
        enddo
      ntmp=nadd+ntmp-1
c
      return
      end
      subroutine phint53(stmpp,etmpp,ntmp,ethi,pirt,rrrt,piht,rrcl,
     $ abund1,abund2,ptmp1,ptmp2,xpx,opakab,rnist,
     $ opakc,rccemis,lpri,epi,ncn2,bremsa,t,swrat,xnx,
     $ lfast,lun11)
c
c
c     this routine does the integration over the spectrum as required by
c     photo.
c     uses power law piecewise analytic integrals, assuming \sigma~e^{-3}
c        for recombination integrals.
c     special version for photemis
c
      implicit none
c
      integer ntmp
c
      include './PARAM'
c
      real*8  bremsa(ncn),epi(ncn),stmpp(ntmp),etmpp(ntmp)
      real*8 rccemis(2,ncn),opakc(ncn)
      real*8 sgbar(ncn),bremsint(ncn)
      integer lpri,ncn2,lfast,lun11
      real*8 ethi,pirt,rrrt,piht,rrcl,abund1,abund2,ptmp1,ptmp2,xpx,
     $     opakab,t,swrat,xnx,crit
      real*8 eth,ergsev,bk,tm,bktm,ener,sgtmp,epii,sgtp,optmp,
     $     sumr,sumh,sumi,sumc,tempi,
     $     bremtmp,tempr,sumho,exptst,
     $     atmp2,rnist,tsq,ethsht,optmp2,
     $     wwir,wwih,bbnurjp,e2t,e2,e1,e1o,bremtmpp,
     $     epiip,e2o,tempc,tempcp,exptmpp,rctmp2,rctmp1,
     $     s2,s2to,s2o,s2t,sgtpp,sum,tempip,temphp,temprp,temph,
     $     wwirp,wwicp,wwic,wwii,wwihp,wwiip,enermx,exptsto
      integer lprisv,numcon2,nphint,nb1,kl,jk,nbn,klmax,nbinc
c
      logical done
c
      data ergsev/1.602197e-12/
      data bk/1.38062e-16/
c
c
c     internal to this routine we use the threshold binned
c      eth=ethi
      nb1=nbinc(ethi,epi,ncn2)
      eth=epi(nb1)
      lprisv=lpri
c
      tm=t*1.e4
      bktm=bk*tm/ergsev
      tsq = sqrt(t)
c      rnist=(2.61e-21)*swrat/t/tsq
c
      if (lpri.ge.1)
     $  write (lun11,*)'in phint53:',
     $      eth,xnx,swrat,t,
     $      stmpp(1),etmpp(1),ntmp,lfast,ptmp1,ptmp2
     $      ,abund1,abund2,rnist,crit
      ethsht=eth/bktm
      ethsht=max(ethsht,0.)
c
      numcon2=max(2,ncn2/50)
      nphint=ncn2-numcon2
c
      sumr = 0.
      sumh = 0.
      sumho = 0.
      sumc = 0.
      sumi=0.
      ener=eth+etmpp(1)*(13.605692)
      nb1=nbinc(ener,epi,ncn2)
      do while ((epi(nb1).lt.ener).and.(nb1.lt.nphint))
        nb1=nb1+1
        enddo
      nb1=nb1-1
c      if (lpri.ge.1) write (lun11,*)'nb1=',nb1,ener,nphint
      if (nb1.ge.nphint) return
c      if (epi(kl+1).lt.ethi) return
      tempr=0.

c     step through cross section, map onto kl grid
      jk=1
      enermx=eth+etmpp(ntmp)*(13.605692)
      nbn=nbinc(enermx,epi,ncn2)
      nbn=max(nbn,min(nb1+1,ncn2-1))
      ener=eth+etmpp(jk)*(13.605692)
      sgtmp=stmpp(1)
      sgbar(max(1,nb1-1))=0.
      kl=nb1
      jk=1
      e1=epi(kl)
      e2=eth+etmpp(jk)*(13.605692)
      s2=stmpp(jk)
      if (e1.lt.e2) then
        kl=kl+1
        e1=epi(kl)
        endif
      e1o=e2
      sum=0.
      done=.false.
      do while (.not.done)
c       step through where jk grid is finer
c        if (lpri.ne.0) write (lun11,*)'mapping loop:',jk,kl,e2,e1
        do while ((e2.lt.e1).and.(jk.lt.(ntmp-1)))
          jk=jk+1
          e2o=e2
          s2o=s2
          e2=eth+etmpp(jk)*(13.605692)
          s2=stmpp(jk)
          sum=sum+(s2+s2o)*(e2-e2o)/2.
c          if (lpri.ne.0) write (lun11,*)'jk loop',jk,e2,s2,e2o,s2o,sum
          enddo
c       kl bin exceeds jk bin, subtract off extra
        sum=sum-(s2+s2o)*(e2-e2o)/2.
c       now interpolate to find value at kl bin
        e2t=e1
        if (e2-e2o.gt.1.e-8) then
             s2t=s2o+(s2-s2o)*(e2t-e2o)/(e2-e2o+1.e-24)
           else
             s2t=s2o
           endif
        s2to=s2t
c       now update sum at kl bin
        sum=sum+(s2t+s2o)*(e2t-e2o)/2.
c       save
        if (abs(e1-e1o).gt.1.e-34) then
            sgbar(kl)=sum/(e1-e1o)
          else
            sgbar(kl)=0.
          endif
c        if (lpri.ne.0) write (lun11,*)'saving:',kl,e1,sgbar(kl),sum,s2t
        e1o=e1
c       increment kl
        kl=kl+1
        e1=epi(kl)
c       step through where kl grid is finer
        do while (e1.lt.e2)
          e2t=e1
          if (e2-e2o.gt.1.e-8) then
              s2t=s2o+(s2-s2o)*(e2t-e2o)/(e2-e2o)
            else
              s2t=s2o
            endif
          s2to=s2t
          sum=(s2t+s2to)*(e1-e1o)/2.
          sgbar(kl)=sum/(e1-e1o)
c          if (lpri.ne.0) write (lun11,*)'kl loop',kl,e1,
c     $                                   sgbar(kl),sum,s2t
          e1o=e1
          kl=kl+1
          e1=epi(kl)
          enddo
c       update sum for remaining bit
        sum=(s2+s2t)*(e2-e2t)/2.
c        if (lpri.ne.0) write (lun11,*)'testing for done:',kl,nphint,
c     $                                                    jk,ntmp
        if ((kl.gt.nphint-1).or.(jk.ge.ntmp-1))
     $       done=.true.
        enddo
      klmax=kl-1
c
c
c     preliminary setup
      sgtpp=sgbar(nb1)
      bremtmpp=bremsa(nb1)/(12.56)
      epiip=epi(nb1)
      temprp=(12.56)*sgtpp*bremtmpp/epiip
      temphp=temprp*epiip
      exptst=(epiip-eth)/bktm
      exptmpp=exp(-exptst)
      bbnurjp=epiip*epiip*epiip*(1.571e+22)*2.
      tempip=rnist*(bremtmpp+bbnurjp)
     $  *sgtpp*exptmpp/epiip*(ptmp1+ptmp2)
      tempcp=tempip*epiip
c
      kl=nb1
      epii=epi(kl)
c      if (lpri.ne.0) write (lun11,*)'kl=',kl,klmax,sumh,sumho
      rctmp1=0.
      rctmp2=0.
      do while (kl.lt.klmax)
c
c       the basics
        sgtmp=max(0.,sgbar(kl))
        sgtp=sgtmp
        sgtpp=sgbar(kl+1)
        bremtmp=bremsa(kl)/(12.56)
        bremtmpp=bremsa(kl+1)/(12.56)
        epii=epi(kl)
        epiip=epi(kl+1)
c
c       pi rate
        tempr=temprp
        temprp=(12.56)*sgtpp*bremtmpp/epiip
        wwir=(epiip-epii)/2.
        wwirp=wwir
        sumr = sumr + (tempr*wwir+temprp*wwirp)
c
c       heat
        temph=temphp
        temphp=temprp*epiip
        wwih=wwir
        wwihp=wwih
        sumho=sumh
        sumh = sumh + (temph*wwih+temphp*wwihp)
c
c       rec
        exptsto=exptst
        exptst=(epiip-eth)/bktm
        if ((exptsto.lt.200.).and.(lfast.ge.2)) then
          bremtmpp=bremsa(kl+1)/(12.56)
          exptmpp=exp(-exptst)
          bbnurjp=epiip*epiip*epiip*(1.571e+22)*2.
          tempi=tempip
          tempip=rnist*(bremtmpp+bbnurjp)
     $        *sgtpp*exptmpp*12.56/epiip
          atmp2=tempip*epiip
          tempip=tempip*(ptmp1+ptmp2)
          wwii=wwir
          wwiip=wwir
          sumi = sumi + (tempi*wwii+tempip*wwiip)
c
c         cool
          tempc=tempcp
          tempcp=tempip*epiip
          wwic=wwir
          wwicp=wwir
          sumc = sumc+tempc*wwic+tempcp*wwicp
c
          rctmp1=abund2*atmp2*ptmp1*xpx/12.56
          rccemis(1,kl)=rccemis(1,kl)+rctmp1
          rctmp2=abund2*atmp2*ptmp2*xpx/12.56
          rccemis(2,kl)=rccemis(2,kl)+rctmp2
c
          endif
c
c       emiss and opac
c       the emission must be fudged to get the right cooling with a
c         trapezoid integration.
        optmp=abund1*xpx*sgtp
        opakc(kl)=opakc(kl)+optmp
c
        if (kl.le.(nb1+1)) then
          optmp2=rnist*exptmpp*sgtp*abund2*(ptmp1+ptmp2)*xpx
          opakab=optmp-optmp2
          endif
c
c       print
        if (lpri.ge.1) then
          write (lun11,*)jk,kl,epi(kl),kl,nphint
          write (lun11,901)jk,kl,epi(kl),sgtp,bremtmp,tempr,sumr
c     $                     ,exptsto,tempi,sumi,tempip,wwir
 901      format(1x,2i6,10(1pe11.3))
          endif
c     $      write (lun11,901)kl,epii,sgtp,bremtmp,
c     $         tempr,temprp,wwir,wwirp,sumr,
c     $         temph,temphp,wwih,wwihp,sumh,
c     $         tempi,tempip,wwii,wwiip,sumi,
c     $         tempc,tempcp,wwic,wwicp,sumc
c     $         ,rctmp1,rctmp2,exptst
c 901       format(1x,'found something',i6,25(1pe11.3))
        if (sgtp.lt.0.) then
          write (6,*) 'phint error'
          return
          endif
c
        kl=kl+1
        enddo
c
      pirt = pirt + sumr
      rrrt = rrrt + sumi
      piht = piht + sumh*ergsev
      rrcl = rrcl + sumc*ergsev

      if (lpri.ge.1) write (lun11,*)'in phint53:',eth,pirt,rrrt
     $         ,piht,rrcl
      lpri=lprisv
c
      do kl=nb1,nbn
        sgbar(kl)=0.
        enddo
c
c
      return
      end
      subroutine phint53hunt(stmpp,etmpp,ntmp,ethi,pirt,rrrt,piht,rrcl,
     $ lpri,epi,ncn2,bremsa,t,swrat,xnx,crit,lfast,lun11)
c
c
c     this routine does the integration over the spectrum as required by
c     photo.
c     this is my version from 11/1/99 which performs successive bisections
c      g stands for good
c     author:  T. Kallman
c
      implicit none
c
      integer ntmp
c
      include './PARAM'
c
      real*8  bremsa(ncn),epi(ncn),stmpp(ntmp),etmpp(ntmp)
      real*8 ansar1(ncn),ansar2(ncn)
      integer luse(ncn)
      integer lpri,ncn2,lfast,lun11
      real*8 ethi,pirt,rrrt,piht,rrcl,
     $     t,swrat,xnx,crit
      real*8 eth,ergsev,bk,tm,bktm,ener,epii,
     $     sumr,sumh,sumi,sumc,tempi,enero,
     $     bremtmp,tempr,tempro,deld,sumho,exptst,
     $     tempi1,tempi2,tempio,atmp2,rnist,atmp2o,bbnurj,delt,
     $     exptmp,emaxx,efnd,ethsht,etst,sumro,sumio,tst3,tst2,
     $     tst,sumco,tst4,tst1,tsq,tsti,sgtmp
      integer lprisv,numcon2,nphint,nb1,kl,itmp,jlo,lprif,
     $     nskp,ndelt,npass,numcon,numcon3,nbinc
c
c
      data ergsev/1.602197e-12/
      data bk/1.38062e-16/
      data delt/1.e-28/
c
c     initialize.
       pirt =0.
       rrrt =0.
       piht=0.
       rrcl=0.
c
      nb1=nbinc(ethi,epi,ncn2)+1
      eth=ethi
      numcon=ncn2
      numcon2=max(2,ncn2/50)
      numcon3=numcon-numcon2
      if (lpri.ge.1) write (lun11,*)'in phint53:',
     $      eth,xnx,swrat,t,nb1,
     $      stmpp(1),etmpp(1),ntmp,lfast
      if (nb1.ge.numcon3) return
c
      tm=t*1.e4
      bktm=bk*tm/ergsev
      tsq = sqrt(t)
      rnist=(5.216e-21)*swrat/t/tsq
c
c     from levwk
      ethsht=eth/bktm
      ethsht=max(ethsht,0.)
c
      lprisv=lpri
c
c     first find range
      kl=numcon3/2
      emaxx=etmpp(ntmp)*(13.605692)+eth
      nphint=nbinc(emaxx,epi,ncn2)
      ndelt=nphint-nb1
      if (lpri.ge.1) write (lun11,*)'in phint53:',
     $      eth,rnist,xnx,swrat,t,ndelt,nphint,nb1,
     $      stmpp(1),etmpp(1),etmpp(ntmp),ntmp,lfast
      ndelt=max(ndelt,1)
      delt=float(ndelt)
      itmp=int(log(delt)/(0.69315)+0.5)
 1011 continue
      ndelt=2**itmp
      nphint=nb1+ndelt
      etst=0.
      if (nphint.le.numcon3)
     $ etst=(epi(nphint)-eth)/13.605692
c      write (lun11,*)ndelt,nphint,etst,numcon3,etmpp(ntmp),ntmp,itmp
      if ((nphint.gt.numcon3).or.(etst.gt.etmpp(ntmp))) then
        itmp=itmp-1
        if (itmp.gt.1) go to 1011
        endif
      if (lpri.ge.1) write (lun11,*)'in phint53:',
     $      eth,rnist,xnx,swrat,t,ndelt,nphint,nb1,
     $      stmpp(1),etmpp(1),ntmp,lfast
c
      do kl=max(1,nb1-1),nphint
        luse(kl)=0
        enddo
c
c     now step thru successive approximations
      nskp=ndelt
      sumr = 0.
      sumh = 0.
      sumc = 0.
      sumi=0.
      npass=0
      do while (((tst3.gt.crit).or.(tst1.gt.crit)
     $    .or.(tst2.gt.crit).or.(tst4.gt.crit).or.(sumi.le.1.e-24))
     $    .and.(nskp.gt.1))
        npass=npass+1
        nskp=max(1,nskp/2)
c
        sumro=sumr
        sumho=sumh
        sumio=sumi
        sumco=sumc
        sumr = 0.
        sumh = 0.
        sumc = 0.
        sumi=0.
        tempr = 0.
        tempi=0.
        atmp2=0.
        jlo=1
        ener = epi(nb1)
        kl=nb1-1
        kl=max(kl,1)
        do while (kl.le.nphint)
          enero=ener
          epii=epi(kl)
          ener=epii
          bremtmp=bremsa(kl)/(25.3)
          tempio=tempi
          atmp2o=atmp2
          sgtmp=0.
          if (ener.ge.eth) then
            if (luse(kl).eq.0) then
                efnd=(ener-eth)/13.605692
                lprif=0
                if (lpri.gt.1) lprif=1
                call find53(stmpp,etmpp,ntmp,efnd,sgtmp,jlo,lun11,lprif)
                ansar1(kl)=sgtmp
                if (lprif.ge.1) write (lun11,*)'after find53:',
     $                  jlo,efnd,sgtmp
                exptst=(epii-eth)/bktm
                exptmp=exp(-exptst)
                bbnurj=epii*epii*epii
                tempi1=rnist*bbnurj*sgtmp*exptmp*(1.571e+22)/epii
                tempi2=rnist*bremtmp*sgtmp*exptmp/epii
c                tempi2=0.
                tempi=tempi1+tempi2
                atmp2=tempi*epii
                ansar2(kl)=atmp2
              else
                sgtmp=ansar1(kl)
                atmp2=ansar2(kl)
                tempi=atmp2/epii
              endif
            endif
          tempro=tempr
          tempr=(25.3)*sgtmp*bremtmp/epii
          deld = ener - enero
          tst=(tempr+tempro)*deld/2.
          sumr = sumr + tst
          sumh=sumh+(tempr*ener+tempro*enero)*deld/2.
          tsti = (tempi+tempio)*deld/2.
          sumi = sumi + tsti
          sumc = sumc+(atmp2+atmp2o)*deld/2.
          if ((lpri.ge.1).and.(npass.le.1)) then
              write (lun11,*)kl,ener,luse(kl),sgtmp,bremtmp,atmp2
              write (lun11,*)kl,ener,bremtmp,bbnurj
              write (lun11,*) tempr,tempi,rnist,exptmp,bbnurj
              write (lun11,*) sumr,sumi,sumh,sumc,sgtmp
              endif
          luse(kl)=1
          kl=kl+nskp
          enddo
c
c
        tst3=abs((sumio-sumi)/(sumio+sumi+1.e-24))
        tst1=abs((sumro-sumr)/(sumro+sumr+1.e-24))
        tst2=abs((sumho-sumh)/(sumho+sumh+1.e-24))
        tst4=abs((sumco-sumc)/(sumco+sumc+1.e-24))
        if (lpri.ge.1) write (lun11,*)'after pass:',npass,
     $     sumr,sumh,sumi,sumc,tst1,tst2,tst3,tst4,nskp
        enddo
c
         pirt = pirt + sumr
         rrrt = rrrt + xnx*sumi
         piht = piht + sumh*ergsev
         rrcl = rrcl + xnx*sumc*ergsev
c
c
         if (lpri.ge.1) write (lun11,*)'in phint53:',eth,pirt,rrrt
     $         ,piht,rrcl,npass,sumc,ergsev
         lpri=lprisv
c
      return
      end
      subroutine phint53pl(sth,e1,alph,ethi,pirt,rrrt,piht,rrcl,
     $ abund1,abund2,ptmp1,ptmp2,xpx,opakab,
     $ opakc,rccemis,lpri,epi,ncn2,bremsa,t,swrat,xnx,
     $ lfast,lun11)
c
c
c     this routine does the integration over the spectrum as required by
c     photo.
c     uses power law piecewise analytic integrals, assuming \sigma~e^{-3}
c        for recombination integrals.
c
      implicit none
c
      include './PARAM'
c
      real*8  bremsa(ncn),epi(ncn)
      real*8 rccemis(2,ncn),opakc(ncn)
      real*8 sth,e1,alph
      integer lpri,ncn2,lfast,lun11
      real*8 ethi,pirt,rrrt,piht,rrcl,abund1,abund2,ptmp1,ptmp2,xpx,
     $     opakab,t,swrat,xnx
      real*8 eth,ergsev,bk,tm,bktm,ener,sgtmp,epii,sgtp,optmp,
     $     sumr,sumh,sumi,sumc,tempi,atmp2,
     $     bremtmp,tempr,sumho,exptst,rnist,tsq,ethsht,temphp,
     $     bremtmpp,bbnurjp,exptmpp,rctmp2,rctmp1,wwih,wwic,wwir,
     $     sgtpp,tempcp,tempc,temph,tempip,temprp,wwicp,wwiip,wwii,
     $     wwihp,epiip,wwirp
      integer lprisv,numcon2,nphint,nb1,kl,klmax,nbinc
c
      data ergsev/1.602197e-12/
      data bk/1.38062e-16/
c
      eth=ethi
      lprisv=lpri
c
      tm=t*1.e4
      bktm=bk*tm/ergsev
      tsq = sqrt(t)
      rnist=(2.61e-21)*swrat/t/tsq
c
      if (lpri.ge.1) write (lun11,*)'in phint53pl:',
     $      eth,xnx,swrat,t,
     $      lfast,ptmp1,ptmp2
     $      ,abund1,abund2
      ethsht=eth/bktm
      ethsht=max(ethsht,0.)
c
      numcon2=max(2,ncn2/50)
      nphint=ncn2-numcon2
c
      sumr = 0.
      sumh = 0.
      sumho = 0.
      sumc = 0.
      sumi=0.
      ener=ethi+13.605692*e1
      nb1=nbinc(ener,epi,ncn2)
      do while ((epi(nb1).lt.ener).and.(nb1.lt.nphint))
        nb1=nb1+1
        enddo
      if (lpri.ge.1) write (lun11,*)'nb1=',nb1,ener,nphint
      if (nb1.ge.nphint) return
c      if (epi(kl+1).lt.ethi) return
      tempr=0.
      rctmp1=0.
      rctmp2=0.
c
c     preliminary setup
      sgtpp=sth
      bremtmpp=bremsa(nb1)/(12.56)
      epiip=epi(nb1)
      temprp=(12.56)*sgtpp*bremtmpp/epiip
      temphp=temprp*epiip
      exptst=(epiip-eth)/bktm
      exptmpp=exp(-exptst)
      bbnurjp=epiip*epiip*epiip*(1.571e+22)*2.
      tempip=rnist*(bremtmpp+bbnurjp*(ptmp1+ptmp2))
     $  *sgtpp*exptmpp/epiip
      tempcp=tempip*epiip
      klmax=ncn2
c
      kl=nb1
      epii=epi(kl)
      if (lpri.ne.0) write (lun11,*)'kl=',kl,klmax,sumh,sumho
      do while ((kl.lt.klmax)
     $    .and.(abs(sumh/(sumho+1.e-24)-1.).gt.1.e-6))
c
c       the basics
        sgtmp=max(0.,sth*(epi(kl)/ener)**alph)
        sgtp=sgtmp
        sgtpp=sth*(epi(kl+1)/ener)**alph
        bremtmp=bremsa(kl)/(12.56)
        bremtmpp=bremsa(kl+1)/(12.56)
        epii=epi(kl)
        epiip=epi(kl+1)
c
c       pi rate
        tempr=temprp
        temprp=(12.56)*sgtpp*bremtmpp/epiip
        wwir=(epiip-epii)/2.
        wwirp=wwir
        sumr = sumr + (tempr*wwir+temprp*wwirp)
c
c       heat
        temph=temphp
        temphp=temprp*epiip
        wwih=wwir
        wwihp=wwih
        sumho=sumh
        sumh = sumh + (temph*wwih+temphp*wwihp)
c
c       rec
        exptst=(epiip-eth)/bktm
        if (exptst.lt.30.) then
          bremtmpp=bremsa(kl+1)/(25.3)
          exptmpp=exp(-exptst)
          bbnurjp=epiip*epiip*epiip*(1.571e+22)*2.
          tempi=tempip
          tempip=rnist*(bremtmpp+bbnurjp*(ptmp1+ptmp2))
     $        *sgtpp*exptmpp/epiip
          wwii=wwir
          wwiip=wwir
          sumi = sumi + (tempi*wwii+tempip*wwiip)
c
c         cool
          tempc=tempcp
          tempcp=tempip*epiip
          wwic=wwir
          wwicp=wwir
          sumc = sumc+tempc*wwic+tempcp*wwicp
c
          atmp2=tempc
          rctmp1=abund2*atmp2*ptmp1*xpx*xnx
          rctmp2=abund2*atmp2*ptmp2*xpx*xnx
          rccemis(1,kl)=rccemis(1,kl)+rctmp1
          rccemis(2,kl)=rccemis(2,kl)+rctmp2
c         print *, "Do we make it here (in phint53pl)?"
c
          endif
c
c       emiss and opac
c       the emission must be fudged to get the right cooling with a
c         trapezoid integration.
        optmp=abund1*xpx*sgtp
        if (kl.le.(nb1+1)) opakab=optmp
        opakc(kl)=opakc(kl)+optmp
c
c       print
         if ((lpri.ge.1).and.(abs(rctmp1+rctmp2).gt.1.e-24)
     $         .and.(ener.gt.6000.))
     $      write (lun11,901)kl,epii,sgtp,bremtmp,
     $         tempr,temprp,wwir,wwirp,sumr,
     $         temph,temphp,wwih,wwihp,sumh,
     $         tempi,tempip,wwii,wwiip,sumi,
     $         tempc,tempcp,wwic,wwicp,sumc
     $         ,rctmp1,rctmp2
 901       format(1x,'found something',i6,25(1pe11.3))
        if (sgtp.lt.0.) then
          write (6,*) 'phint error'
          return
          endif
c
        kl=kl+1
        enddo
c
      pirt = pirt + sumr
      rrrt = rrrt + xnx*sumi
      piht = piht + sumh*ergsev
      rrcl = rrcl + xnx*sumc*ergsev

      if (lpri.ge.1) write (lun11,*)'in phint53:',eth,pirt,rrrt
     $         ,piht,rrcl
      lpri=lprisv
c
c
      return
      end
      subroutine phint53old(stmpp,etmpp,ntmp,ethi,pirt,rrrt,piht,rrcl,
     $ abund1,abund2,ptmp1,ptmp2,xpx,opakab,rnist,
     $ opakc,rccemis,lpri,epi,ncn2,bremsa,t,swrat,xnx,
     $ lfast,lun11)
c
c
c     this routine does the integration over the spectrum as required by
c     photo.
c     this is my version from 11/1/99 which performs successive bisections
c      g stands for good
c
      implicit none
c
      integer ntmp
c
      include './PARAM'
c
      real*8  bremsa(ncn),epi(ncn),stmpp(ntmp),etmpp(ntmp)
      real*8 rccemis(2,ncn),opakc(ncn)
      integer lpri,ncn2,lfast,lun11
      real*8 ethi,pirt,rrrt,piht,rrcl,abund1,abund2,ptmp1,ptmp2,xpx,
     $     opakab,t,swrat,xnx
      real*8 eth,ergsev,bk,tm,bktm,ener,sgtmp,epii,sgtp,ansar1,optmp,
     $     sumr,sumh,sumi,sumc,tempi,atmp2,enero,sgtmpo,dels,
     $     sgmn,sgmx,bremtmp,tempr,tempro,epiio,deld,sumho,exptst,
     $     tempi1,tempi2,tempio,rnist,tsq,tempi1o,tempi2o,optmp2,
     $     atmp2o,bbnurj,exptmp,atmp1,atmp1o,bbnu,sss,xee,rcctot,
     $     rccemiso,expo
      integer lprisv,numcon2,nphint,nb1,kl,jk,
     $     nbinc,jkp
c
      data ergsev/1.602197e-12/
      data bk/1.38062e-16/
c
c     some constants
      eth=ethi
      lprisv=lpri
      tm=t*1.e4
      bktm=bk*tm/ergsev
      tsq = sqrt(t)
c     this constant has 4 pi built in
c      rnist=(2.61e-21)*swrat/t/tsq
c
c      print input parameters
      if (lpri.ge.1) write (lun11,*)'in phint53:',
     $      eth,rnist,xnx,swrat,t,
     $      stmpp(1),etmpp(1),ntmp,lfast,abund1,abund2
     $      ,ptmp1,ptmp2
c
c     limits of energy indeces
      numcon2=max(2,ncn2/50)
      nphint=ncn2-numcon2
c
c     zero temporaries
      sumr = 0.
      sumh = 0.
      sumho=0.
      tempr = 0.
      sumc = 0.
      sumi=0.
      tempi=0.
      atmp1=0.
      atmp2=0.
c
c     starting values for integration
      ener=ethi+etmpp(1)*(13.605692)
      nb1=nbinc(ener,epi,ncn2)
      kl=nb1
      if (kl.ge.nphint) return
      if (epi(kl+1).lt.ethi) return
      sgtmp=stmpp(1)
      dels=0.
      epii=epi(kl)
      jk=0
      rcctot=0.
      jkp=min(jk+1,ntmp)
c
c     step through cross section indeces
      do while ((jk.lt.ntmp).and.
     $          (ethi+etmpp(jkp)*(13.605692).lt.epi(nphint-1)))
c     $    .and.(abs(sumh/(sumho+1.e-24)-1.).gt.1.e-6))
c
c       get cross section dependent quantities
        jk=jk+1
        enero=ener
        ener=ethi+etmpp(jk)*(13.605692)
        sgtmpo=sgtmp
        sgtmp=stmpp(jk)
c
c       test for whether cross section grid energy
c         is greater than master energy grid energy
        if (epi(kl+1).lt.ener) then
c
c         calculate derivitive for interpolation
          dels=(sgtmp-sgtmpo)/(ener-enero+1.e-24)
          sgmn=min(sgtmp,sgtmpo)
          sgmx=max(sgtmp,sgtmpo)
c
c         step through photon energies
          do while ((kl.lt.nphint).and.(epi(kl+1).lt.ener))
c
c           the cross section
            kl=kl+1
            sgtp=max(sgmn,min(sgmx,sgtmpo+dels*(epi(kl)-enero)))
            ansar1=sgtp
c
c           the photoionization rate
            bremtmp=bremsa(kl)/(12.56)
            tempro=tempr
            epiio=epii
            epii=epi(kl)
            tempr=(12.56)*sgtp*bremtmp/epii
            deld = epii - epiio
            sumr = sumr + (tempr+tempro)*deld/2.
c
c           the photoionization heating rate
            sumho=sumh
            sumh=sumh+(tempr*epii+tempro*epiio)*deld/2.
c
c           skip if lfast=1
            if (lfast.gt.1) then
c
c             quantities associated with recombination
c             this expression causes trouble when the threshold is not
c              at eth
c              exptst=(-eth+epii)/bktm
c             this is a possible fix
              exptst=(-(eth+max(0.,etmpp(1)*(13.605692)))+epii)/bktm
              exptmp=0.
c
c             test for exponential in rrc
              if (exptst.lt.200.) then
c
                exptmp=expo(-exptst)
                bbnurj=epii*epii*epii*(1.571e+22)*2.
                tempi1=rnist*bbnurj*exptmp*sgtp*12.56/epii
                tempi2=rnist*bremtmp*exptmp*sgtp*12.56/epii
                tempi1=tempi1*(ptmp1+ptmp2)
                tempi2=tempi2*(ptmp1+ptmp2)
                tempi1o=tempi1
                tempi2o=tempi2
                Atmp1o=atmp1
                Atmp2o=atmp2
c               the old way with stimulated recombination treated as emission
c                tempi=(tempi1+tempi2)*(ptmp1+ptmp2)
                atmp1=tempi1*epii
                atmp2=tempi2*epii
c
c               the recombination rate
                sumi = sumi + (tempi1+tempi1o+tempi2+tempi2o)*deld/2.
c
c               the recombination cooling rate
                sumc = sumc+(atmp1+atmp1o+atmp2+atmp2o)*deld/2.
c
c               the rrc emissivity
                rccemis(1,kl)=rccemis(1,kl)+
     $                 abund2*atmp1*ptmp1*xpx/12.56
                rccemis(2,kl)=rccemis(2,kl)+
     $                abund2*atmp1*ptmp2*xpx/12.56
c
c               total rrc emissivity
                if (kl.gt.1)
     $           rcctot=rcctot+(rccemis(1,kl)+rccemis(1,kl-1)+
     $                       rccemis(2,kl)+rccemis(2,kl-1))
     $              *(epi(kl)-epi(kl-1))*(1.602197e-12)*12.56/2.
c
c               end of test for exponential in rrc
                endif
c
              if (kl.le.(nb1+1)) then
                optmp=abund1*ansar1*xpx
c               need to take out the 4 pi
                optmp2=rnist*exptmp*sgtp*abund2*(ptmp1+ptmp2)*xpx
                opakab=optmp-optmp2
                endif

c             test for lfast=3
              if (lfast.gt.2) then
c
c               opacity
                optmp=abund1*ansar1*xpx
c               need to take out the 4 pi
                optmp2=rnist*exptmp*sgtp*abund2*(ptmp1+ptmp2)*xpx
                if (kl.le.(nb1+1)) opakab=optmp-optmp2
c               new way treating stimulated recombination as negative absorption
                opakc(kl)=opakc(kl)+max(0.,optmp-optmp2)
c
c               the source function and photon occupation number
                bbnu=bbnurj/(expo(epii/bktm)-1.)
c               for a successful comparison must write s in sr^-1
                sss=(rccemis(1,kl)+rccemis(2,kl))/(1.e-34+opakc(kl))
c
c               end of test for lfast=3
                endif
c
c             end of test for lfast=1
              endif
c
c           diagnostic print
            if (lpri.ge.1)
     $        write (lun11,902)kl,jk,epi(kl),sgtp,bbnurj,bremtmp,bbnu,
     $          atmp2,tempr,tempi1,tempi2,sumi,optmp,optmp2,
     $          opakc(kl),rccemis(1,kl),rccemis(2,kl),
     $          sss,sss/(1.e-34+bbnu),
     $          rcctot,rcctot/abund2/xpx,sumc*ergsev,exptst,exptmp
 902        format (1x,2i6,23(1pe12.4))
c
c           end of loop over master energy grid
            enddo
c
c         end of test for energy grid points
          endif
c
c       end of loop over cross section grid
        enddo
c
      pirt = pirt + sumr
      rrrt = rrrt + sumi
      piht = piht + sumh*ergsev
      rrcl = rrcl + sumc*ergsev

      if (lpri.ge.1) write (lun11,*)'in phint53:',eth,pirt,rrrt
     $         ,piht,rrcl
      lpri=lprisv
c
c
      return
      end
      subroutine phintfo(sigc,ethi,pirt,rrrt,piht,rrcl,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lpri,epi,ncn2,bremsa,t,swrat,xnx,lfast,lun11)
c
c     this routine does the integration over the spectrum as required by
c     type 59 data
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
      real*8  bremsa(ncn),epi(ncn)
      real*8 opakc(ncn)
      real*8 sigc(ncn)
      integer lpri,ncn2,lfast,lun11
      real*8 ethi,pirt,rrrt,piht,rrcl,abund1,abund2,xpx,
     $     opakab,t,swrat,xnx
      real*8 eth,ergsev,bk,tm,bktm,ener,sgtmp,epii,ansar1,optmp,
     $     sumr,sumh,sumi,sumc,tempi,enero,
     $     bremtmp,tempr,tempro,exptst,
     $     tempi1,tempi2,tempio,atmp2,rnist,tsq,ethsht,
     $     atmp2o,bbnurj,exptmp,deld,tsti,
     $     bbnu,htmpo,htmp,tmpp,sigth,htsum,eps,tst
      integer lprisv,nphint,nb1,kl,
     $     lprie,nskp,nbinc,lrcalc,numcon2
c
c
      data ergsev/1.602197e-12/
      data bk/1.38062e-16/
c
      data eps/1.e-2/
c
      nb1=nbinc(ethi,epi,ncn2)
      eth=ethi
c
      tm=t*1.e4
      bktm=bk*tm/ergsev
      tsq = sqrt(t)
      rnist=(5.216e-21)*swrat/t/tsq
c
c     from levwk
      ethsht=eth/bktm
      ethsht=max(ethsht,0.)
c
c     initialize.
      pirt =0.
      rrrt =0.
c
      lprisv=lpri
c
c     continuum
      tempro = 0.
      tempr = 0.
      sumr = 0.
      sumh = 0.
      sumc = 0.
      tempio = 0.
      tempi=0.
      sumi=0.
      htsum=0.
      htmp=0.
      sgtmp = 0.
      numcon2=max(2,ncn2/50)
      nphint=ncn2-numcon2
      nphint=max(nphint,nb1+1)
      if (lpri.ge.1) write (lun11,*)'in phintf:',
     $      eth,rnist,xnx,swrat,t,abund1,abund2,
     $      nphint,nb1
      ener = epi(nb1)
      sgtmp=0.
      bbnu=0.
      kl=nb1
      atmp2=0.
      tst=1.e+10
      do while ((kl.le.nphint).and.
     $            ((lfast.le.2).or.(tst.gt.eps)))
            enero=ener
            ener = epi(kl)
            epii=ener
            sgtmp = sigc(kl)
            if (kl.eq.nb1) sigth=sgtmp
            bremtmp=bremsa(kl)/(25.3)
            tempro=tempr
            tempr=(25.3)*sgtmp*bremtmp/epii
            deld = ener - enero
            tst=(tempr+tempro)*deld/2.
            ansar1=sgtmp
            sumr = sumr + tst
            sumh=sumh+(tempr*ener+tempro*enero)*deld*ergsev/2.
            exptst=(epii-eth)/bktm
            exptmp=exp(-exptst)
            bbnurj=epii*epii*epii*(1.571e+22)
            tempi1=rnist*bbnurj*exptmp*sgtmp/epii
            tempi2=rnist*bremtmp*exptmp*sgtmp/epii
            tempi=tempi1+tempi2
            atmp2o=atmp2
            atmp2=tempi1*epii
            tsti = (tempi+tempio)*deld/2.
            sumi = sumi + tsti
            sumc = sumc+(atmp2+atmp2o)*deld*ergsev/2.
c            rctmp1=abund2*ansar2*ptmp1*xpx*xnx
c            rctmp2=abund2*ansar2*ptmp2*xpx*xnx
c            rccemis(1,kl)=rccemis(1,kl)+rctmp1
c            rccemis(2,kl)=rccemis(2,kl)+rctmp2
            optmp=abund1*ansar1*xpx
            if (kl.le.nb1+1) opakab=optmp
            htmpo=htmp
            tmpp=optmp
            htmp=bremsa(kl)*tmpp
            htsum=htsum+(htmp+htmpo)
     $               *(epi(kl)-epi(kl-1))*(1.602197e-12)/2.
            opakc(kl)=opakc(kl)+optmp
            if (lpri.gt.1) then
              write (lun11,*)kl,ener,bremtmp,bbnu,bbnurj
              write (lun11,*) tempr,tempi,rnist,exptmp,bbnurj
              write (lun11,*) sumr,sumi,sumh,sumc,sgtmp
              endif
            tempro = tempr
            tempio = tempi
            enero = ener
            tst=1.
            lprie=0
            call enxt(ethi,nb1,lprie,epi,ncn2,t,lfast,lun11,
     $                  kl,nskp,nphint,lrcalc)
            kl=kl+nskp
            enddo
c
         pirt = pirt + sumr
         rrrt = rrrt + xnx*sumi
         piht = piht + sumh
         rrcl = rrcl + xnx*sumc

         if (lpri.ge.1)
     $    write (lun11,*)'in phintf:',eth,sigth,
     $     pirt,rrrt,piht,rrcl
         lpri=lprisv
c
      return
      end
      subroutine pprint(jj,jkstep,
     $ tp,xlum,lwri,lpri,r,t,xpx,p,lcdd,numrec,npass,
     $ nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,zremsz,epi,ncn2,
     $ abel,cfrac,emult,taumax,xeemin,spectype,specfile,specunit,
     $ kmodelname,nloopctl,nparms,parname,partype,parms,parcomm,
     $ atcredate,lun11,tinf,xcol,vturbi,critf,radexp,
     $ delr,rdel,enlum,xee,ababs,
     $ bremsa,bremsint,tau0,dpthc,tauc,
     $ idat1,rdat1,kdat1,nptrs,np2,
     $ npar,npnxt,npfi,npfirst,
     $ nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $ npconi2,ncsvn,
     $ ntotit,lnerrd,
     $ xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,
     $         cllines,clcont,htcomp,clcomp,clbrems,
     $ xilev,bilev,rniss,
     $ rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,opakab,
     $ cabab,elumab,elum,zrems)
c
c     this routine prints
c     author:  T. Kallman
c
c     variable categories:
c     step indeces:
c      jj,jkstep,
c     input parameters:
c      tp,xlum,lwri,lpri,xpx,p,lcdd,numrec,npass,
c      nnmax,nlimd,rmax,xpxcol,xi,zeta,lfix,
c      abel,cfrac,emult,taumax,xeemin,
c      tinf,xcol,vturbi,critf,ababs,
c     other input
c      spectype,specfile,specunit,kmodelname,nloopctl,
c      nparms,parname,partype,parms,parcomm,lun11,
c     input spectrum
c      zremsz,epi,ncn2,bremsa,
c     database quantities
c      idat1,rdat1,kdat1,nptrs,np2,
c      npar,npnxt,npfi,npfirst,nplin,nplini,
c      npcon,npconi,npilev,npilevi,npconi2,
c      nlsvn,ncsvn,
c     step diagnostics
c      ntotit,lnerrd,
c     state variables
c      r,t,xee,xii,xilev,bilev,rniss,
c     derived state variables
c      delr,rdel,enlum,
c     rates
c      rrrt,pirt,htt,cll,httot,cltot,hmctot,
c      cllines,clcont,htcomp,clcomp,clbrems,
c     emissivities and opacities
c      rcem,oplin,rccemis,brcems,opakc,cemab,opakab,
c     optical depths and luminosities
c      tau0,dpthc,tauc,elumab,elum,zrems
c
c
C     A plethora of printing options...
C
C        jj
C         1 - 500 strongest emissions lines, sorted by strength
C         19 - RRC luminosities
C         23 - 500 strongest absorption lines, sorted by strength
C         24 - absorption edges
C         2 - print input parameter listC
c         4 - continuum opacity and emissivity
c         5 - energy sums
c         6 - continuum luminosities and depths
c         8 - line list
c         10 - ion abundances and thermal rates (erg/sec)
C        11 - Write FITS file with summary of ion abundances
C        12 - append abundance values to the data array for xout_abund1.fits
C             Doesn't actually write the file, just accumulates values.
c        13 - blank space
c        14 - line opacities and emissivities
c        15 - line luminosities
C        18 - line wavelengths and levels
c        20 - line finding list
c        21 - level opacities and emissivities
c         7 - level populations
C        17 - print column headings for this pass
C         9 - print short summary line of the radial zone
c        16 - times
c        22 - ionization parameter etc.
c        25 - outputting to common block
c        26 - ferland print
C
C     Modifications:
C        1998/12/17, WTB: Fix FITS keyword format problem for
C                       writeascii routine.  Removed dependence on
C                       writeimage routine.
C        1999/01/04, WTB: Added model name to writeascii parameter list
C        1999/01/05, WTB: Removed log(Xi)& log(U1) columns from calls #11
C                       & #12
C                       Inserted '_' in spaces for ion names in call #11
C
      implicit none
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
      integer idat1(nidat1),
     $      nptrs(nptt,ndat2)
c     $      ,np1r,np1i,np1k,np2
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     line luminosities
      real*8 elum(3,nnnl)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
c     line optical depths
      real*8 tau0(2,nnnl)
c     energy bins
      real*8 epi(ncn)
c     continuum lum
      real*8 zrems(4,ncn),
     $          zremsz(ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
c     continuum flux
      real*8 bremsa(ncn)
      real*8 bremsint(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn),brcems(ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     level populations
      real*8 xilev(nnml),bilev(nnml),rniss(nnml)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
      real*8 elumab(2,nnml)
      real*8 tauc(2,nnml)
c     ion abundances
      real*8 xii(nni)
c     heating/cooling
      real*8 htt(nni),cll(nni)
c     the atomic data creation date
      character(63) atcredate
      real*8 rrrt(nni),pirt(nni)
      real*8 abel(nl),ababs(nl)
      real*8 xcoltmp(nni)
      integer kltmp(5000)
      real*8 zrtmp(999,3999),zrtmpcol(999,3999),
     $   zrtmpc(999,3999),zrtmph(999,3999)
c      real*8 epi2(ncn),zrems2(3,ncn)
c
c      common /ewout/newout,lnewo(nnnl),kdewo(8,nnnl),
c     $  kdewol(20,nnnl),kdewou(20,nnnl),aijewo(nnnl),flinewo(nnnl),
c     $  ggloewo(nnnl),ggupewo(nnnl),
c     $  elewo(nnnl),tau0ewo(nnnl),elout(2,nnnl),zrtmp,epi2,zrems2
c     A feature added November 2007 is output of the strongest lines,
c     sorted by element and ion into a common block called 'ewout'
c     The contents of the common block are:
c       newout:  number of lines in the list.
c       lnewo:   array conatining line indexes.
c       kdewo:   character array containing the name of the ion
c       kdewol:  character array containing the name of the lower level
c       kdewou:  character array containing the name of the upper level
c       aijewo:  array containing A values for the lines
c       flinewo: array containing f values for the lines
c       ggloewo: array containing statistical weights for the lower levels
c       ggupewo: array containing statistical weights for the upper levels
c       elewo:   array containing the line wavelengths
c       tau0ewo: array containing the line center depths
c       elout:   array containing line luminosities in xstar units (erg/s/10^38)
c       zrtmp: a 2d array 999x3999, containing a zone-by-zone summary of the
c                state of the gas.  The second index is the zone number
c                the first index is as follows:
c                1: radius (cm)
c                2: log(xi)
c                3: electron fraction (relative to nuclei)
c                4: nucleus number density
c                5: pressure (dynes/cm^2)
c                6: tempeatre/10^4 K
c                7:  heating-cooling/(heating+cooling)
c                8-..: ion fractions for all the ions in the model
c                    a model with non-zero abundance for all elements
c                    will have 168 ions (excluding bare nuclei)
c                    numbered such that 1=H0, 2=He0, 3=He+, 4=C0,
c                    ... 168=Ni27+.  In this case the upper limit
c                    for these columns will be 168+8=176
c       epi: energy grid in eV, length=99999
c       zrems: spectrum in erg/s/erg/10**38.  This is a 2d array 3x99999,
c           where column 1=transmitted outward flux, including diffuse
c           emission, column 2=diffuse inward emission (no direct flux)
c           column 3=diffuse outward emission (no direct flux)
c
c
c       character(1) kdewo,kdewol,kdewou
c       character(1) klevl(20),klevu(20)
c       real*8  flinewo, aijewo,
c      real*8 ggloewo, ggupewo, elewo, tau0ewo
c      real*8 elout
c      integer nilino, jkktmp, lup, lnewo, lupfnd, llofnd
c      integer newout
c
      character(20) parname(55)
      character(10) partype(55)
      real*8 parms(55)
      character(30) parcomm(55),kmodelname
      character(8) spectype, specfile
      character(1) kdtmp(100),kblnk,klablo(20),klabup(20)
      integer nparms, specunit, nloopctl
      character(8) kabstring(30)
      character(8) ktmp8
      character(20) ktmp20,klevu,klevl,kblnk20
      character(9) kinam1
      character(133) tmpst
      character(16) knam,klabs(3999),kunits(3999),kform(3999),ktmp,
     $              kblnk16
      character(1) klev(100,nd),kdat(nptmpdim)
      integer klen, unit, status
      real*8  rlev(10,nd)
      integer ilv(10,nd),nlpt(nd),iltp(nd)
      real*8  elsv(nnnl)
      integer jpnt(nnnl)
c jg
      real*8 xnx, xpx, xee, eliml, elimh
      real*8 elmmtpp, elcomp, xeltp, cabcompare
      real*8 abund1, cfrac, t, p, tp, xlum
      real*8 xpxcol, zeta, taumax, xeemin, critf, radexp
      real*8 vturbi, opsum, tstar, fstr, rsum1
      real*8 rsum2, sgtmp, tmp, crayj, fstro
      real*8 emult, ekkr, delte, rssmn, elsum, ergsev
      real*8 sumtmp1, sumtmp2, r19, tmp1, tmp2, tmp1o
      real*8 tmp2o, err, sum1, sum2, sum3, sum4
      real*8 r, ener, etst, aij, gglo, ggup, flin
      real*8 httot, cltot, htcomp, clcont, cllines
      real*8 clcomp, clbrems, uu1, enlum, alguu1
      real*8 skse, ecc, ekt, sksec, zetac, enn0
      real*8 egam, rdel, hmctot, vvthermsc
      real*8 elmtp, elmtpb, ethi, ethc, terr
      real*8 ett, optpp, optppo, tmp2c, xcol,fpr2
      real*8 tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
      real*8 ttot, enlumx
      real*8 uux, alguux, eth, abundel
      real*8  xi, delr, elin, flux1, flux2
      real*8 rmax, tinf, rdum, dep, rss, bbe, rocc

      integer jj, lun11, lpril, kltmpo, nlplmx, lm
      integer nlpl, lnn, nlsvn, ln, ml, ltyp, lrtyp
      integer lcon, nrdt, nidt, nkdt, nilin
      integer lmm, kl2, k
      integer kltmpn, mm, kk, j, ipmat, klel, mlel
      integer jkk, nnz, klion,mlion, jk, mt2, mllel
      integer kl, nkdti, mlleltp, nlevmx, mltype
      integer mlpar, lk, nilin2, mllz, nlev, kkkl, idest1
      integer llo, np2, nlevp, idest2, jkk3, ndtmp
      integer mllz2, iltmp, ltyp2, lrtyp2, lcon2, nrdt2
      integer nidt2, nkdt2, lcdd, numrec, nlimd, lwri, lpri
      integer lfix, npass, numcon, ncn2, i, jlk, lun11sv
      integer ktt, lfnd, nell, lkk, nelin, jkl, mmlv, nidti
      integer nlyc, nry, jkstep, ilevup, ilevlo, jkko
      integer niter, jjj, jp1, npi, npc, mll, ltypc
      integer ilevt, npio, mllev, mlcu, mm2, mmtmp, ll
      integer ntotit, nb1, nb10
      integer lnerrd, nbinc
      integer nnmax, ncsvn,mlm
      integer np1i,np1r,np1k,np1i2,np1r2,np1k2,np1ki
c
      logical done
c
c     Not used
      integer javi
      real*8 javir
c      character(80) javik
c
      save zrtmpc,zrtmph,zrtmp
c
      data kblnk20/'                    '/
      data kblnk/' '/,kblnk16/'                '/
      data kabstring/'H abund=','Heabund=','Liabund=',
     $               'Beabund=','B abund=','C abund=','N abund=',
     $               'O abund=','F abund=','Neabund=','Naabund=',
     $               'Mgabund=','Alabund=','Siabund=','P abund=',
     $               'S abund=','Clabund=','Arabund=','K abund=',
     $               'Caabund=','Scabund=','Tiabund=','V abund=',
     $               'Crabund=','Mnabund=','Feabund=','Coabund=',
     $               'Niabund=','Cuabund=','Znabund='/
c

      javi=nnmax
      javir=rmax
      javi=nparms
c      javik=parname(1)
c      javik=partype(1)
      javir=tinf
      javir=bremsa(1)
      bremsa(1)=javir
      javi=nplini(1)
      javi=npilevi(1)
      javi=npconi(1)
      javi=ncsvn
      javir=bilev(1)
      javir=xi
      javi=lnerrd
c      lnerrd=javi
c

      xnx=xpx*xee
c
      if ((jj.ne.9).and.(jj.ne.12))
     $ write (lun11,9211)jj
 9211 format (1x, 'print option:',i2)
      if ((jj.le.0).or.(jj.gt.27)) return
c
      goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     $23,24,25,26,27),
     $  jj
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
1     continue
c
      lpril=0
c     500 strongest emission lines, sorted by strength
      write (lun11,*)'emission line luminosities (erg/sec/10**38))'
      kltmpo=0
      nlplmx=500
      eliml=0.1
      elimh=1.0e10
c     find the strongest lines.
      do  lm=1,nlplmx
        kltmp(lm)=0
        enddo
c
c     step through lines
      nlpl=1
      do lnn=1,nlsvn
c
c       get line data
        ln=lnn
        ml=nplin(ln)
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        elin=abs(rdat1(np1r))
c
c       exclude rate type 14
        if (lrtyp.ne.14) then
c
c         get ion data
          nilin=npar(ml)
          mlm=nilin-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          nilin2=idat1(np1i+nidt-1)
c
c         get lum and test for strength, wavelength
          elmmtpp=(elum(2,ln)+elum(1,ln))/2.
          if (lpril.ne.0)
     $       write (lun11,*)lnn,elin,nilin,elmmtpp,ln,ml
          if ((ln.gt.0).and.(ln.lt.nnnl)
     $       .and.(elin.ge.eliml).and.(elin.le.elimh)
     $       .and.(elin.le.8.9e+6)
     $       .and.(elmmtpp.gt.1.e-36)
     $       .and.(nilin2.gt.0).and.(nilin2.le.nni))
     $        then
c
c           insertion sort
            lmm=0
            elcomp=1.e+10
            do while ((lmm.lt.nlpl).and.(elmmtpp.lt.elcomp))
              lmm=lmm+1
              kl2=kltmp(lmm)
              elcomp=0.
              if (kl2.gt.0)
     $          elcomp=(elum(2,kl2)+elum(1,kl2))/2.
              enddo
            if (lpril.ne.0)
     $       write (lun11,8516)ln,elin,elmmtpp,lmm,nlpl,kl2,elcomp
 8516       format (1h ,i4,2e12.4,3i4,e12.4)
            kltmpo=ln
            do  k=lmm,min(nlplmx,nlpl)
              if ((lpril.ne.0).and.(kltmp(k).ne.0))
     $          write (lun11,*)'in 557 loop',k,kltmp(k),kltmpo
              kltmpn=kltmp(k)
              kltmp(k)=kltmpo
              kltmpo=kltmpn
              enddo
           nlpl=min(nlplmx,nlpl+1)
           if (lpril.ne.0)
     $       write (lun11,*)'done with 557 loop',lm
            endif
c           end of insertion
c
          endif
c
        enddo
c
      if (nlpl.gt.0) kltmp(nlpl)=kltmpo
c
c     printing loop
      write (lun11,959)
c
c     step through lines
      do  kk=1,nlpl
        if (lpril.ne.0)
     $    write (lun11,*)'kk=',kk
        ln=kltmp(kk)
        if (ln.ne.0) then
c
c         get line data
          ml=nplin(ln)
          if (ml.ne.0) then
            if (lpril.ne.0)
     $      write (lun11,*)'   ',ln,ml
            mlm=ml-1
            call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $        nptrs,0,lun11)
            elin=abs(rdat1(np1r))
c
c           get ion data
            nilin=npar(ml)
            mlm=nilin-1
            call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $        nptrs,0,lun11)
            do mm=1,nkdt
              kdtmp(mm)=kdat1(np1k-1+mm)
              enddo
            do mm=nkdt+1,9
              kdtmp(mm)=kblnk
              enddo
c            nilin=idat1(np1i+2)
            if (lpril.ne.0)
     $      write (lun11,*)ml,nilin,npar(ml)
c
c           print
            write (lun11,9955)kk,ln,(kdtmp(mm),mm=1,9),elin,
     $      elum(1,ln),elum(2,ln)
 9955       format (1x,2i8,1x,9a1,3(1pe13.5))
c
            endif
c
          endif
c
        enddo
c
      write (lun11,993)
c
      return
c
c
 19   continue
c
c     print 500 strongest recombination continua
      lpril=0
      write (lun11,*)'recombination continuum luminosities',
     $  '(erg/sec/10**38))'
      write (lun11,*)'index, ion, level, energy (eV), RRC luminosity '
c
C     lpril is flag for printing debug information
      if (lpril.ne.0) then
        write (lun11,*)'raw data'
        do j=1,nnml
          if (xilev(j).gt.1.e-37)
     $     write (lun11,*)j,xilev(j),elumab(1,j)
          enddo
        endif
c
C     First look for element data (jk is element index)
      klel=11
      mlel=npfirst(klel)
      jk=0
      kk=0
      jkk=0
c
c     step through elements
      do while (mlel.ne.0)
c
c       get element data
        jk=jk+1
        mt2=mlel-1
        call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,
     $        nptrs,0,lun11)
        mllel=idat1(np1i+nidt-1)
        xeltp=rdat1(np1r)
        xeltp=abel(mllel)
        nnz=idat1(np1i)
        if (lpril.ge.1)
     $        write (lun11,*)'element:',jk,mlel,mllel,nnz,
     $                  (kdat1(np1k-1+mm),mm=1,nkdt)
c
C       ignore if the abundance is small
        if (xeltp.lt.1.e-10) then
            jkk=jkk+nnz
          else
c
c           now step thru ions (jkk is ion index)
            klion=12
            mlion=npfirst(klion)
            jkk=0
            kl=0
            do while ((mlion.ne.0).and.(kl.lt.nnz))
              jkk=jkk+1
c
C             retrieve ion name from kdati
              mlm=mlion-1
              call drd(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidti,np1i,nkdti,np1ki,mlm,
     $            nptrs,0,lun11)
c
C             if not accessing the same element, skip to the next element
              mlleltp=idat1(np1i+nidti-2)
              if (mlleltp.eq.mllel) then
c
                kl=kl+1
                if (lpril.ge.1)
     $            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
     $                        (kdat1(np1ki+mm-1),mm=1,nkdti)
c
c               now find level data
                call func2l(jkk,lpril,lun11,t,xee,xpx,
     $              idat1,rdat1,kdat1,nptrs,
     $              npar,npnxt,npfi,
     $              rniss,rlev,ilv,
     $              nlpt,iltp,nlev,klev)
c
c               now step through rate type 7 data
                mltype=7
                ml=npfi(mltype,jkk)
                mllz=0
                if (ml.ne.0) mllz=npar(ml)
                mlpar=0
                if (ml.ne.0) mlpar=npar(ml)
                do while ((ml.ne.0).and.(mlpar.eq.mllz))
c
c                 get rrc data
                  kkkl=npconi2(ml)
                  if (lpril.ne.0) write (lun11,*)kkkl,ml,idest1,
     $                    elumab(1,kkkl),elumab(2,kkkl)
c
c                 test for non-zero rrc data
                  if ((kkkl.gt.0).and.(kkkl.le.ndat2)
     $                .and.((elumab(1,kkkl).gt.1.e-36)
     $                .or.(elumab(2,kkkl).gt.1.e-36))) then
c
c                   get rrc  data
                    mlm=ml-1
                    call drd(ltyp,lrtyp,lcon,
     $                nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $                nptrs,0,lun11)
                    idest1=idat1(np1i+nidt-2)
                    nlevp=nlev
                    idest2=nlevp+idat1(np1i-1+nidt-3)-1
c
c                   label for lower level
                    do lk=1,20
                      write (ktmp20(lk:lk),'(a1)')klev(lk,idest1)
                      enddo
                    klevl=ktmp20
c
c                   label for upper level
                    write (ktmp20(1:20),'(a20)')'continuum           '
                    klevu=ktmp20
c
c                   ion label
                    do lk=1,nkdti
                      write (ktmp8(lk:lk),'(a1)')kdat1(np1ki+lk-1)
                      enddo
                    do lk=nkdti+1,8
                      write (ktmp8(lk:lk),'(a1)')kblnk
                      enddo
c
                    eth=rlev(4,idest1)-rlev(1,idest1)
                    ett=eth
c
c                   get upper level data
                    if (idest2.gt.nlevp) then
                      jkk3=jkk+1
                      if (lpril.gt.1)
     $                  write (lun11,*)jkk3,ndtmp,nlevp,idest2
                      ndtmp=npfi(13,jkk3)
                      if (lpril.gt.1)
     $                  write (lun11,*)jkk3,ndtmp,nlevp,idest2
                      if (ndtmp.le.0) stop 'ndtmp error'
                      mllz=npar(ndtmp)
                      iltmp=0
                      do while ((ndtmp.ne.0).and.
     $                    (iltmp.ne.(idest2-nlevp+1)).and.
     $                    (npar(ndtmp).eq.mllz))
                        mlm=ndtmp-1
                        call drd(ltyp2,lrtyp2,lcon2,
     $                    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $                    nptrs,0,lun11)
                        iltmp=idat1(np1i2+nidt2-2)
                        if (lpril.gt.1) write (lun11,*)nidt2,iltmp,ndtmp
                        ndtmp=npnxt(ndtmp)
                        if (ndtmp.le.0) stop 'ndtmp error'
                        enddo
c                     NB fix to excited level PI and rec
                      ett=ett+rdat1(np1r2)
                      eth=ett
                      if (lpril.gt.1)
     $                  write (lun11,*) ndtmp,iltmp,idest2,ett
c                     label for lower level
                      ktmp20=kblnk20
                      do lk=1,nkdt2
                        write (ktmp20(lk:lk),'(a1)')kdat1(np1k2+lk-1)
                        enddo
                      klevu=ktmp20
                      endif
c
c                   other data
                    mmlv=npilev(idest1,jkk)
                    write (lun11,9293)kkkl,mmlv,ktmp8,idest1,idest2,
     $                  klevl,klevu,eth,elumab(1,kkkl),elumab(2,kkkl)
c
c                   done with this rrc
                    endif
c
c                 end of loop over rrcs
                  ml=npnxt(ml)
                  if (ml.ne.0) mlpar=npar(ml)
                  enddo
c
c               end of test for element
                endif
c
C             Go to next ion
              mlion=npnxt(mlion)
              enddo
c
c         end of test for non-zero element abund
          endif
c
        mlel=npnxt(mlel)
C       Go to next element
        enddo
c
      write (lun11,993)
c
      return
c
c
 23   continue
c
      lpril=0
c     print 500 strongest absoprtion lines
      write (lun11,*)'line depths'
      kltmpo=0
      nlplmx=500
      eliml=0.1
      elimh=1.0e10
c     find the strongest lines.
      do  lm=1,nlplmx
        kltmp(lm)=0
        enddo
c
c     step through lines
      nlpl=1
      do lnn=1,nlsvn
c
c       get line data
        ln=lnn
        ml=nplin(ln)
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        elin=abs(rdat1(np1r))
c
c       exclude rate type 14
        if (lrtyp.ne.14) then
c
c         get ion data
          nilin=npar(ml)
          mlm=nilin-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          nilin2=idat1(np1i+nidt-1)
c
c         get lum and test for strength, wavelength
          elmmtpp=tau0(1,ln)
          if (lpril.ne.0)
     $       write (lun11,*)lnn,elin,nilin,elmmtpp,ln,ml
          if ((ln.gt.0).and.(ln.lt.nnnl)
     $       .and.(elin.ge.eliml).and.(elin.le.elimh)
     $       .and.(elin.le.8.9e+6)
     $       .and.(elmmtpp.gt.1.e-36)
     $       .and.(nilin2.gt.0).and.(nilin2.le.nni))
     $        then
c
c           insertion sort
            lmm=0
            elcomp=1.e+10
            do while ((lmm.lt.nlpl).and.(elmmtpp.lt.elcomp))
              lmm=lmm+1
              kl2=kltmp(lmm)
              elcomp=0.
              if (kl2.gt.0)
     $          elcomp=tau0(1,kl2)
              enddo
            if (lpril.ne.0)
     $       write (lun11,8516)ln,elin,elmmtpp,lmm,nlpl,kl2,elcomp
            kltmpo=ln
            do  k=lmm,min(nlplmx,nlpl)
              if ((lpril.ne.0).and.(kltmp(k).ne.0))
     $          write (lun11,*)'in 557 loop',k,kltmp(k),kltmpo
              kltmpn=kltmp(k)
              kltmp(k)=kltmpo
              kltmpo=kltmpn
              enddo
           nlpl=min(nlplmx,nlpl+1)
           if (lpril.ne.0)
     $       write (lun11,*)'done with 557 loop',lm
            endif
c           end of insertion
c
          endif
c
        enddo
c
      if (nlpl.gt.0) kltmp(nlpl)=kltmpo
c
c     printing loop
      write (lun11,959)
 959  format (1x,'index, ion, wavelength, transmitted, reflected')
c
c     step through lines
      do  kk=1,nlpl
        if (lpril.ne.0)
     $    write (lun11,*)'kk=',kk
        ln=kltmp(kk)
        if (ln.ne.0) then
c
c         get line data
          ml=nplin(ln)
          if (ml.ne.0) then
            if (lpril.ne.0)
     $      write (lun11,*)'   ',ln,ml
            mlm=ml-1
            call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $        nptrs,0,lun11)
            elin=abs(rdat1(np1r))
c
c           get ion data
            nilin=npar(ml)
            mlm=nilin-1
            call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $        nptrs,0,lun11)
            do mm=1,nkdt
              kdtmp(mm)=kdat1(np1k-1+mm)
              enddo
            do mm=nkdt+1,9
              kdtmp(mm)=kblnk
              enddo
c            nilin=idat1(np1i+2)
            if (lpril.ne.0)
     $      write (lun11,*)ml,nilin,npar(ml)
c
c           print
            write (lun11,9955)kk,ln,(kdtmp(mm),mm=1,9),elin,
     $      tau0(1,ln),tau0(2,ln)
c
            endif
c
          endif
c
        enddo
c
      write (lun11,993)
c
      return
c
c
 24   continue
c
      lpril=0
c     print 500 strongest absorption edges
      write (lun11,*)'absorption edge depths'
      write (lun11,*)'index, ion, level, energy (eV), depth '
c
C     lpril is flag for printing debug information
      if (lpril.ne.0) then
        write (lun11,*)'raw data'
        do j=1,nnml
          if (xilev(j).gt.1.e-37)
     $     write (lun11,*)j,xilev(j),tauc(1,j)
          enddo
        endif
c
C     First look for element data (jk is element index)
      klel=11
      mlel=npfirst(klel)
      jk=0
      kk=0
      jkk=0
c
c     step through elements
      do while (mlel.ne.0)
c
c       get element data
        jk=jk+1
        mt2=mlel-1
        call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,
     $        nptrs,0,lun11)
        mllel=idat1(np1i+nidt-1)
        xeltp=rdat1(np1r)
        xeltp=abel(mllel)
        nnz=idat1(np1i)
        if (lpril.ge.1)
     $        write (lun11,*)'element:',jk,mlel,mllel,nnz,
     $                  (kdat1(np1k-1+mm),mm=1,nkdt)
c
C       ignore if the abundance is small
        if (xeltp.lt.1.e-10) then
            jkk=jkk+nnz
          else
c
c           now step thru ions (jkk is ion index)
            klion=12
            mlion=npfirst(klion)
            jkk=0
            kl=0
            do while ((mlion.ne.0).and.(kl.lt.nnz))
              jkk=jkk+1
c
C             retrieve ion name from kdati
              mlm=mlion-1
              call drd(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidti,np1i,nkdti,np1ki,mlm,
     $            nptrs,0,lun11)
c
C             if not accessing the same element, skip to the next element
              mlleltp=idat1(np1i+nidti-2)
              if (mlleltp.eq.mllel) then
c
                kl=kl+1
                if (lpril.ge.1)
     $            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
     $                        (kdat1(np1ki+mm-1),mm=1,nkdti)
c
c               now find level data
                call func2l(jkk,lpril,lun11,t,xee,xpx,
     $              idat1,rdat1,kdat1,nptrs,
     $              npar,npnxt,npfi,
     $              rniss,rlev,ilv,
     $              nlpt,iltp,nlev,klev)
c
c               now step through rate type 7 data
                mltype=7
                ml=npfi(mltype,jkk)
                mllz=0
                if (ml.ne.0) mllz=npar(ml)
                mlpar=0
                if (ml.ne.0) mlpar=npar(ml)
                do while ((ml.ne.0).and.(mlpar.eq.mllz))
c
c                 get rrc data
                  kkkl=npconi2(ml)
                  if (lpril.ne.0) write (lun11,*)kkkl,ml,idest1,
     $                    elumab(1,kkkl),elumab(2,kkkl)
c
c                 test for non-zero rrc data
                  if ((kkkl.gt.0).and.(kkkl.le.ndat2)
     $                .and.((elumab(1,kkkl).gt.1.e-36)
     $                .or.(elumab(2,kkkl).gt.1.e-36))) then
c
c                   get rrc  data
                    mlm=ml-1
                    call drd(ltyp,lrtyp,lcon,
     $                nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $                nptrs,0,lun11)
                    idest1=idat1(np1i+nidt-2)
                    nlevp=nlev
                    idest2=nlevp+idat1(np1i-1+nidt-3)-1
c
c                   label for lower level
                    do lk=1,20
                      write (ktmp20(lk:lk),'(a1)')klev(lk,idest1)
                      enddo
                    klevl=ktmp20
c
c                   label for upper level
                    write (ktmp20(1:20),'(a20)')'continuum           '
                    klevu=ktmp20
c
c                   ion label
                    do lk=1,nkdti
                      write (ktmp8(lk:lk),'(a1)')kdat1(np1ki+lk-1)
                      enddo
                    do lk=nkdti+1,8
                      write (ktmp8(lk:lk),'(a1)')kblnk
                      enddo
c
                    eth=rlev(4,idest1)-rlev(1,idest1)
                    ett=eth
c
c                   get upper level data
                    if (idest2.gt.nlevp) then
                      jkk3=jkk+1
                      if (lpril.gt.1)
     $                  write (lun11,*)jkk3,ndtmp,nlevp,idest2
                      ndtmp=npfi(13,jkk3)
                      if (lpril.gt.1)
     $                  write (lun11,*)jkk3,ndtmp,nlevp,idest2
                      if (ndtmp.le.0) stop 'ndtmp error'
                      mllz=npar(ndtmp)
                      iltmp=0
                      do while ((ndtmp.ne.0).and.
     $                    (iltmp.ne.(idest2-nlevp+1)).and.
     $                    (npar(ndtmp).eq.mllz))
                        mlm=ndtmp-1
                        call drd(ltyp2,lrtyp2,lcon2,
     $                    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $                    nptrs,0,lun11)
                        iltmp=idat1(np1i2+nidt2-2)
                        if (lpril.gt.1) write (lun11,*)nidt2,iltmp,ndtmp
                        ndtmp=npnxt(ndtmp)
                        if (ndtmp.le.0) stop 'ndtmp error'
                        enddo
c                     NB fix to excited level PI and rec
                      ett=ett+rdat1(np1r2)
                      eth=ett
                      if (lpril.gt.1)
     $                  write (lun11,*) ndtmp,iltmp,idest2,ett
c                     label for lower level
                      ktmp20=kblnk20
                      do lk=1,nkdt2
                        write (ktmp20(lk:lk),'(a1)')kdat1(np1k2+lk-1)
                        enddo
                      klevu=ktmp20
                      endif
c
c                   other data
                    mmlv=npilev(idest1,jkk)
                    write (lun11,9293)kkkl,mmlv,ktmp8,idest1,idest2,
     $                  klevl,klevu,eth,tauc(1,kkkl),tauc(2,kkkl)
 9293               format(1x,2i6,1x,a8,2i6,1x,2(a20,1x),3(1pe11.3))
c
c                   done with this rrc
                    endif
c
c                 end of loop over rrcs
                  ml=npnxt(ml)
                  if (ml.ne.0) mlpar=npar(ml)
                  enddo
c
c               end of test for element
                endif
c
C             Go to next ion
              mlion=npnxt(mlion)
              enddo
c
c         end of test for non-zero element abund
          endif
c
        mlel=npnxt(mlel)
C       Go to next element
        enddo
c
      return
c
c
 2    continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C     Print list of input parameters
c
c
c
      write (lun11,'(1x)')
      write (lun11,*)'input parameters:'
      write (lun11,*)'covering fraction=',cfrac
      write (lun11,*)'temperature (/10**4K)=',t
      write (lun11,*)'constant pressure switch (1=yes, 0=no)=',1-lcdd
      write (lun11,*)'pressure (dyne/cm**2)=',p
      write (lun11,*)'density (cm**-3)=',xpx
      write (lun11,*)'spectrum type=',spectype
      write (lun11,*)'spectrum file=',specfile
      write (lun11,*)'spectrum units? (0=energy, 1=photons)',specunit
      write (lun11,*)'radiation temperature or alpha=',tp
      write (lun11,*)'luminosity (/10**38 erg/s)=',xlum
      write (lun11,*)'column density (cm**-2)=',xpxcol
      write (lun11,*)'log(ionization parameter)=',zeta
      do j=1,nl
         write (lun11,*)kabstring(j),abel(j)
         enddo
      write (lun11,*)'model name=',kmodelname
      write (lun11,*)'number of steps=',numrec
      write (lun11,*)'number of iterations=',nlimd
      write (lun11,*)'write switch (1=yes, 0=no)=',lwri
      write (lun11,*)'print switch (1=yes, 0=no)=',lpri
      write (lun11,*)'step size choice switch=',lfix
      write (lun11,*)'loop control (0=standalone)=',nloopctl
      write (lun11,*)'number of passes=',npass
      write (lun11,*)'emult=',emult
      write (lun11,*)'taumax=',taumax
      write (lun11,*)'xeemin=',xeemin
      write (lun11,*)'critf=',critf
      write (lun11,*)'vturbi=',vturbi
      write (lun11,*)'ncn2=',ncn2
      write (lun11,*)'radexp=',radexp
      write (lun11,'(1x)')
      write (lun11,993)
      return

c      r19=r*(1.e-19)
c      uu1=enlum/(12.56*xpx*r19*r19)/3.e+10
c      alguu1=log10(max(1.e-24,uu1))
c      skse=xlum/(xpx*r19*r19)
c      zeta=log10(max(1.e-24,skse))
c      ecc=2.998e+10
c      ekt=t*(0.861707)*ergsev
c      sksec=skse/12.56/((1.+xee)*ekt*ecc)
c      zetac=log10(max(1.e-24,sksec))
c      enn0=xpx
c      nlyc=nbinc(13.7,epi,ncn2)
c      nry=nlyc+1
c      egam=zremsz(nry)/(2.*12.56*enn0*ecc*r19*r19+1.e-24)
c      nry=nbinc(13.6,epi,ncn2)+1
c      write (lun11,993)
c      write (lun11,*)'input parameters:'
c      write (lun11,*)'continuum luminosity=',xlum
c      write (lun11,*)'pressure or density=',xpx
c      write (lun11,*)'radius=',r
c      write (lun11,*)'ionization parameter=',skse
c      write (lun11,993)
c      write (tmpst,993)
c      call xwrite(tmpst,10)
c      write (tmpst,*)'input parameters:'
c      call xwrite(tmpst,10)
c      write (tmpst,*)'continuum luminosity=',xlum
c      call xwrite(tmpst,10)
c      write (tmpst,*)'pressure or density=',xpx
c      call xwrite(tmpst,10)
c      write (tmpst,*)'radius=',r
c      call xwrite(tmpst,10)
c      write (tmpst,*)'ionization parameter=',skse
c      call xwrite(tmpst,10)
c      write (tmpst,993)
c      call xwrite(tmpst,10)


3     continue
c
C     Write the parameter list to the FITS file
C     When changing this list, make sure nparms, parname, partype,
C     and parcomm are also properly updated.  Watch out for the
C     model name which is currently parcomm(37)
      parms(1)=cfrac
      parms(2)=t
      parms(3)=1-lcdd
      parms(4)=p
      parms(5)=xpx
      parms(6)=0.0
      parcomm(6)=spectype
      parms(7)=0.0
      parcomm(7)=specfile
      parms(8)=specunit
      parms(9)=tp
      parms(10)=xlum
      parms(11)=xpxcol
      parms(12)=zeta
      parms(13)=numrec
      parms(14)=nlimd
      parms(15)=lwri
      parms(16)=lpri
      parms(17)=lfix
      do j=1,nl
         parms(17+j)=abel(j)
         enddo
      parms(17+nl+1)=emult
      parms(17+nl+2)=taumax
      parms(17+nl+3)=xeemin
      parms(17+nl+4)=critf
      parms(17+nl+5)=vturbi
      parms(17+nl+6)=npass
      parms(17+nl+7)=0.0
      parcomm(17+nl+7)=kmodelname
      parms(17+nl+8)=nloopctl
c
      return
c
c
4     continue
c
c     print continuum opacity and emissivity
      write (lun11,*)
     $ 'continuum opacity and emissivities (/cm**3/sec/10**38)'
      write (lun11,*)'channel, energy,      opacity,    sigma*e**3,'
     $, 'scattered,  rec. in,   rec. out,  brem. em., source, bbe,',
     $'photon occ'
      opsum=0.
      tstar=t
      ekkr=xnx*(6.65e-25)
      ekkr=max(1.e-20,ekkr)
      optpp=max(opakc(1),ekkr)
      fstr=0.
      rsum1=0.
      rsum2=0.
      numcon=ncn2
c
c     step thru continuum bins
      do 135 kl=2,numcon
c
c        sigma*e**3
         sgtmp=(opakc(kl)*(epi(kl)/1000.)**3)/max(1.e-24,xpx)
c
c        calculate sum and rosseland mean
         if ((kl.gt.1).and.(epi(kl).gt.100.))
     $    opsum=opsum+(opakc(kl)+opakc(kl-1))*(epi(kl)-epi(kl-1))/2.
         i=kl
         tmp = epi(i)*1.16/tstar
         crayj = 1./tmp
         fstro = fstr
         if ( tmp.le.50. ) then
            if ( tmp.gt.1.e-4 ) crayj = 1./(exp(tmp)-1.)
            crayj=crayj*crayj
c            fstr= cconst*tmp*crayj*epi(i)**3/tstar
            fstr= tmp*crayj*epi(i)**3/tstar
         endif
         optppo=optpp
         optpp=max(opakc(kl),ekkr)
         delte=epi(kl)-epi(kl-1)
         rsum1=min(1.e+20,rsum1+(fstr/optpp+fstro/optppo)*delte/2.)
         rsum2=min(1.e+20,rsum2+(fstr+fstro)*delte/2.)
c
c        source function
         rss=(rccemis(1,kl)+rccemis(2,kl)+brcems(kl)/12.56)/
     $           (1.e-34+opakc(kl))
c         rss=(rccemis(1,kl)+rccemis(2,kl))/(1.e-34+opakc(kl))
c
c        planck function
         bbe=2.*epi(kl)**3*(1.5642e+22)
     $       /(exp(epi(kl)/(0.861707*t))-1.+1.e-34)
c
c        photon occupation number
         rocc=rss/(bbe+1.e-34)
c
c        print
         write (lun11,967)kl,epi(kl),opakc(kl),sgtmp,opakscatt(kl),
     $            rccemis(1,kl),rccemis(2,kl), brcems(kl),rss,bbe,rocc
967      format (1h ,i6,10(1pe13.5))
c
135      continue
c
c     print summed opacities
      write (lun11,*)'opsum cont=',opsum
      rssmn=rsum2/rsum1
c      ens1=rssmn/(t*1.e+4)**(-3.5)/xpx/xpx/1.66e-24
      write (lun11,*)'rosseland mean opacity=',t,rssmn
      write (lun11,993)
c
      return
C
5     continue
c
c     print energy sums
      elsum=0.
      ergsev=1.602197e-12
      do jlk=1,nlsvn
         ln=jlk
         ml=nplin(ln)
         mlm=ml-1
         call drd(ltyp,lrtyp,lcon,
     $     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $     nptrs,0,lun11)
         elin=abs(rdat1(np1r))
         if ((elin.lt.1.e+8).and.(elin.gt.1.)) then
           elsum=elsum+elum(1,ln)+elum(2,ln)
           endif
         enddo
      sumtmp1=0.
      sumtmp2=0.
      ergsev=1.602197e-12
      r19=r*1.e-19
      tmp1=zremsz(1)*(1.-exp(-dpthc(1,1)))
      tmp2=zrems(3,1)+zrems(2,1)
      do jk=2,ncn2
         tmp1o=tmp1
         tmp1=zremsz(jk)*(1.-exp(-dpthc(1,jk)))
         sumtmp1=sumtmp1+(tmp1+tmp1o)*(epi(jk)-epi(jk-1))*ergsev/2.
         tmp2o=tmp2
         tmp2=zrems(2,jk)+zrems(3,jk)
         sumtmp2=sumtmp2+(tmp2+tmp2o)*(epi(jk)-epi(jk-1))*ergsev/2.
         enddo
      err=(sumtmp1-sumtmp2-elsum)/(sumtmp1+1.e-24)
      write (lun11,9981)sumtmp1,sumtmp2,elsum,err
 9981 format (1x,'energy sums: abs, cont, line, err:',4(1pe13.5))
      write (lun11,993)
c      write (lun11,*),httot,cltot,fpr2dr,httot/(sumtmp1+1.e-24),
c     $  cltot/(elsum+sumtmp2+1.e-24)
c
      return
c
c
6     continue
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     print continuum luminosities and depths
      write (lun11,*)'continuum luminosities (/sec/10**38) and depths'
      write (lun11,*)'channel,energy,inc.,trn. lum.,ref. lum.,'
     $,'scattered,backward depth,forward depth'
c
      numcon=ncn2
      sum1=0.
      sum2=0.
      sum3=0.
      sum4=0.
      ergsev=1.602197e-12
      r19=r*1.e-19
      fpr2=12.56*r19*r19
c
c     step thru continuum bins
      do kl=1,numcon
c
c       planck function
        bbe=2.*epi(kl)**3*(1.5642e+22)
     $       /(exp(epi(kl)/(0.861707*t))-1.+1.e-34)
c
c       photon occupation number
        rocc=zrems(1,kl)/(bbe+1.e-34)/fpr2/12.56
c
        write (lun11,968)kl,epi(kl),zremsz(kl),
     $    zrems(1,kl),zrems(2,kl),zrems(3,kl),zrems(4,kl),
     $    dpthc(1,kl),dpthc(2,kl),bbe,rocc
c
c       sums
        if (kl.gt.1) then
           sum1=sum1+(zremsz(kl)+zremsz(kl-1))
     $         *(epi(kl)-epi(kl-1))*ergsev/2.
           sum2=sum2+(zrems(1,kl)+zrems(1,kl-1))
     $         *(epi(kl)-epi(kl-1))*ergsev/2.
           sum3=sum3+(zrems(2,kl)+zrems(2,kl-1))
     $         *(epi(kl)-epi(kl-1))*ergsev/2.
           sum4=sum4+(zrems(3,kl)+zrems(3,kl-1))
     $         *(epi(kl)-epi(kl-1))*ergsev/2.
           endif
968      format (1h ,i6,10(1pe13.5))
c
         enddo
c
      write (lun11,*)'norms:'
      write (lun11,9698)sum1,sum2,sum3,sum4
 9698 format (20x,4(1pe13.5))
      write (lun11,993)
c
      return
c
c
8     continue
c

      lpril=0
      write (lun11,*)'line list'
      write (lun11,9943)
9943  format (1x,'     wave          element  ion  glo          gup   '
     $ ,'      fij')
c
c     step through lines
      nlpl=1
      do lnn=1,nlsvn
c
c       get line data
        ln=lnn
        ml=nplin(ln)
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        elin=abs(rdat1(np1r))
c
c
c       exclude rate type 14
        if ((lrtyp.ne.14).and.(abs(elin).gt.0.1)
     $       .and.(abs(elin).lt.9.e+9)) then
c
          ergsev=1.602197e-12
          ener=ergsev*(12398.41)/max(elin,1.e-24)
          etst=ener/ergsev
          idest1=idat1(np1i)
          idest2=idat1(np1i+1)
          aij=rdat1(np1r+2)
          if (ltyp.eq.82) aij=rdat1(np1r+3)
c
c         get ion data
          nilin=npar(ml)
          mlm=nilin-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          do ktt=1,min(8,nkdt)
            write (kinam1(ktt:ktt),'(a1)')kdat1(np1k-1+ktt)
            enddo
          do ktt=nkdt+1,9
            write (kinam1(ktt:ktt),'(a1)')kblnk
            enddo
c
c          if (lpri.ge.1)
c     $      write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
c     $          (kdat1(np1ki+mm-1),mm=1,nkdti)
c
c         now find level data
          jkk=idat1(np1i+nidt-1)
          call func2l(jkk,lpril,lun11,t,xee,xpx,
     $              idat1,rdat1,kdat1,nptrs,
     $              npar,npnxt,npfi,
     $              rniss,rlev,ilv,
     $              nlpt,iltp,nlev,klev)

          ggup=rlev(2,idest1)
          gglo=rlev(2,idest2)
          do lk=1,20
            klablo(lk)=klev(lk,idest1)
            klabup(lk)=klev(lk,idest2)
            enddo
          flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo)
          write (lun11,9944)lnn,elin,kinam1,aij,flin,gglo,ggup,
     $             klablo,klabup
9944      format (1h ,i9,e12.4,1x,a9,4(1pe12.4),1x,20a1,1x,20a1)
c
          endif
        enddo
      write (lun11,993)
c
      return
c
c
10    continue
c
      write (lun11,*)'ion abundances and thermal rates (erg/sec)'
      write (lun11,947)
947   format (1x,'index, ion, abundance, recombination, ionization,',
     $' heating, cooling: ')
c
c     step thru ions
      klion=12
      mlion=npfirst(klion)
      lk=0
      do while (mlion.ne.0)
c
c        get ion data
         lk=lk+1
         ltyp=klion
         mlm=mlion-1
         call drd(ltyp,lrtyp,lcon,
     $     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $     nptrs,0,lun11)
         do mm=1,nkdt
           kdtmp(mm)=kdat1(np1k-1+mm)
           enddo
         do mm=nkdt+1,9
           kdtmp(mm)=kblnk
           enddo
c
c        get element data
         nell=npar(mlion)
         mlm=nell-1
         call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
c         write (lun11,*)mlion,lk,np1i,nidt,np1i+nidt-1,
c     $      idat1(np1i+nidt-1),ababs(idat1(np1i+nidt-1)),mlm
         if ((idat1(np1i+nidt-1).gt.0)
     $     .and.(idat1(np1i+nidt-1).le.nl)) then
           abundel=ababs(idat1(np1i+nidt-1))
c
c          print out
           if (abundel.gt.1.e-15)
     $      write (lun11,9046)lk,(kdtmp(mm),mm=1,9),
     $      xii(lk),rrrt(lk),pirt(lk),htt(lk),cll(lk)
9046       format (1x,i4,1x,9a1,5(1pe16.8))
c
           endif
c
         mlion=npnxt(mlion)
         enddo
c
      write (lun11,*)'total heating, cooling:',
     $            httot,cltot
      write (lun11,*)'partial heating rates: photo,compton',
     $            httot-htcomp,htcomp
      write (lun11,*)'partial cooling rates: rec,lines,brems,compton',
     $            clcont,cllines,clcomp,clbrems
      write (lun11,993)
c
c
      return
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Write FITS file with summary of ion abundances
C     by radial zone
C
11    continue
c
      knam='xout_abund1.fits'
      call fheader(unit,knam,atcredate,kmodelname,status)

c      write (lun11,*)'in pprint, 11: numrec=',numrec
c      call writeimage(knam)
      do mm=1,999
        kunits(mm)=kblnk16
        klabs(mm)=kblnk16
        kform(mm)=kblnk16
        enddo
      klabs(1)='radius          '
      kform(1)='E11.3'
      kunits(1)='cm'
      klabs(2)='delta_r         '
      kform(2)='E11.3'
      kunits(2)='cm'
      klabs(3)='ion_parameter   '
      kform(3)='E11.3'
      kunits(3)='erg*cm/s'
      klabs(4)='x_e             '
      kform(4)='E11.3'
      kunits(4)=' '
      klabs(5)='n_p             '
      kform(5)='E11.3'
      kunits(6)='cm**(-3)'
      klabs(6)='pressure        '
      kform(6)='E11.3'
      kunits(6)='dynes/cm**2'
      klabs(7)='temperature     '
      kform(7)='E11.3'
      kunits(7)='10**4 K'
      klabs(8)='frac_heat_error'
      kform(8)='E11.3'
      kunits(8)=' '
C     Search for the ion names in the database
      klion=12
      mlion=npfirst(klion)
      do lkk=1,nni
        xcoltmp(lkk)=0.
        enddo
      lk=0
      do while (mlion.ne.0)
           lk=lk+1
           ltyp=klion
           mlm=mlion-1
           call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
c
c          get element abundance
           nelin=npar(mlion)
           ml=nelin
           mlm=ml-1
           call drd(ltyp,lrtyp,lcon,
     $       nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $       nptrs,0,lun11)
           mllel=idat1(np1i+nidt-1)
           xeltp=ababs(mllel)
c
c          go back to ion data
           mlm=mlion-1
           call drd(ltyp,lrtyp,lcon,
     $       nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $       nptrs,0,lun11)

C          Compute string length from character array by searching backwards
C          for first non-blank character
           klen=0
           do mm=1,nkdt
             kdat(mm)=kdat1(np1k-1+mm)
             enddo
           do mm=nkdt+1,9
             kdat(mm)=kblnk
             enddo
           do mm=1,9
             if(kdat(10-mm).ne.' '.and.klen.eq.0) then
                klen=10-mm
             endif
             enddo
c           write (lun11,*)'kdat:',(kdat(mm),mm=1,9)
c           write (lun11,*)'klen:',klen
C          Replace ' ' in ion names to '_' to match FITS standard
           do mm=1,9
             if(kdat(mm).eq.' '.and.mm.lt.klen) then
                write (ktmp(mm:mm),'(a1)')'_'
             else
                write (ktmp(mm:mm),'(a1)')kdat(mm)
             endif
             enddo
           do mm=10,16
             write (ktmp(mm:mm),'(a1)')' '
             enddo
           do jkl=2,numrec
c             write (lun11,*)jkl,lk,zrtmp(2,jkl),zrtmp(8+lk,jkl),
c     $                             zrtmp(5,jkl),xeltp,xcoltmp(lk)
             xcoltmp(lk)=xcoltmp(lk)
     $         +(zrtmp(8+lk,jkl)*zrtmp(5,jkl)
     $             +zrtmp(8+lk,jkl-1)*zrtmp(5,jkl-1))
     $         *(zrtmp(2,jkl)-zrtmp(2,jkl-1))*xeltp/2.
             enddo
           klabs(8+lk)=ktmp
           kform(8+lk)='E11.3'
           kunits(8+lk)=' '
           mlion=npnxt(mlion)
           enddo
c
      call fwrtascii(unit,'ABUNDANCES',zrtmp,8+lk,
     $                  numrec,klabs,kform,kunits,lun11)
c
c     calculate columns
      numrec=1
      do lkk=1,lk
        zrtmpcol(8+lkk,numrec)=xcoltmp(lkk)
        enddo
      do lkk=1,8
        zrtmpcol(lkk,numrec)=0.
        enddo
c
      call fwrtascii(unit,'COLUMNS   ',zrtmpcol,8+lk,
     $                  numrec,klabs,kform,kunits,lun11)
c
      klabs(8+lk+1)='compton'
      kform(8+lk+1)='E11.3'
      kunits(8+lk+1)=' '
      klabs(8+lk+2)='total'
      kform(8+lk+2)='E11.3'
      kunits(8+lk+2)=' '
      call fwrtascii(unit,'HEATING                                   ',
     $ zrtmph,8+lk+2,numrec,klabs,kform,kunits,lun11)
c
      klabs(8+lk+1)='compton'
      kform(8+lk+1)='E11.3'
      kunits(8+lk+1)=' '
      klabs(8+lk+2)='brems'
      kform(8+lk+2)='E11.3'
      kunits(8+lk+2)=' '
      klabs(8+lk+3)='total'
      kform(8+lk+3)='E11.3'
      kunits(8+lk+3)=' '
      call fwrtascii(unit,'COOLING                                    ',
     $ zrtmpc,8+lk+3,numrec,klabs,kform,kunits,lun11)
c
      call fitsclose(lun11,unit,status)
c
      return
c
c
12    continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Add ionic abundances info in this radial zone to array for
C     eventual inclusion in xout_abund1.fits
C     Modifies zrtmp
C
      ergsev=1.602197e-12
      r19=r*(1.e-19)
c      write (lun11,*)enlum,xpx,r,xlum,jkstep
      uu1=enlum/(12.56*xpx*r19*r19)/3.e+10
c      write (lun11,*)uu1
      alguu1=log10(max(1.e-24,uu1))
      skse=xlum/(xpx*r19*r19)
      zeta=log10(max(1.e-24,skse))
      ecc=2.998e+10
      ekt=t*(0.861707)*ergsev
c      sksec=skse/(12.56*((1.+xee)*ekt+pradl/(1.e-24+xpx))*ecc)
      sksec=skse/12.56/((1.+xee)*ekt*ecc)
      zetac=log10(max(1.e-24,sksec))
      enn0=xpx
      nlyc=nbinc(13.7d0,epi,ncn2)
      nry=nlyc+1
      egam=zremsz(nry)/(2.*12.56*enn0*ecc*r19*r19+1.e-24)
      nry=nbinc(13.6d0,epi,ncn2)+1
C     Copy the values for radial zone jkstep
      if (jkstep.gt.3999) return
      zrtmp(1,jkstep)=r
      zrtmp(2,jkstep)=rdel
      zrtmp(3,jkstep)=zeta
      zrtmp(4,jkstep)=xee
      zrtmp(5,jkstep)=xpx
      zrtmp(6,jkstep)=p
      zrtmp(7,jkstep)=t
      zrtmp(8,jkstep)=hmctot
      do lk=1,8
        zrtmpc(lk,jkstep)=zrtmp(lk,jkstep)
        zrtmph(lk,jkstep)=zrtmp(lk,jkstep)
        enddo
      klion=12
      mlion=npfirst(klion)
      lk=0
      do while (mlion.ne.0)
        lk=lk+1
        zrtmp(8+lk,jkstep)=xii(lk)
        zrtmpc(8+lk,jkstep)=htt(lk)
        zrtmph(8+lk,jkstep)=cll(lk)
        mlion=npnxt(mlion)
        enddo
      zrtmph(8+lk+1,jkstep)=htcomp
      zrtmph(8+lk+2,jkstep)=httot
      zrtmpc(8+lk+1,jkstep)=clcomp
      zrtmpc(8+lk+2,jkstep)=clbrems
      zrtmpc(8+lk+3,jkstep)=cltot
c
      return
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c
13    continue
c
993   format (1h )
c
      return
c
c
14    continue
c
      write (lun11,900)
900   format ('line opacities and emissivities',
     $ ' (erg/cm**3/sec/10**38)')
      write (lun11,915)
915   format (1x,'index,wavelength,energy,ion,opacity,rec. em.,',
     $'coll. em.,fl. em.,di. em.,cx. em.')
c
c     step through lines
      nlpl=1
c     write (lun11,*)nlsvn
      do lnn=1,nlsvn
c
c       get line data
        ln=lnn
        ml=nplin(ln)
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        elin=abs(rdat1(np1r))
        if ((lrtyp.ne.14).and.(abs(elin).gt.0.1)
     $       .and.(abs(elin).lt.9.e+9)) then

          ergsev=1.602197e-12
          ener=ergsev*(12398.41)/max(elin,1.e-24)
          etst=ener/ergsev
          idest1=idat1(np1i)
          idest2=idat1(np1i+1)
          aij=rdat1(np1r+2)
c
c         get ion data
          nilin=npar(ml)
          mlm=nilin-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
           do ktt=1,min(8,nkdt)
            write (kinam1(ktt:ktt),'(a1)')kdat1(np1k-1+ktt)
            enddo
          do ktt=nkdt+1,9
            write (kinam1(ktt:ktt),'(a1)')kblnk
            enddo
c
          j=ln
          write (lun11,904)j,elin,etst,kinam1,oplin(j),rcem(1,j),
     $                      rcem(2,j)
904       format (1h ,i9,2(1pe13.5),1x,a9,6(1pe13.5))
c
          endif
        enddo
      write (lun11,993)
c
      return
c
c
 15   continue
c
      write (lun11,*)'line luminosities (erg/sec/10**38) and depths'
      write (lun11,9923)
9923  format (1x,' line, wavelength, ion, ref. lum.,trn. lum.,',
     $'backward depth, forward depth')
c     step through lines
      nlpl=1
      do lnn=1,nlsvn
c
c       get line data
        ln=lnn
        ml=nplin(ln)
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        elin=abs(rdat1(np1r))
        if ((lrtyp.ne.14).and.(abs(elin).gt.0.1)
     $       .and.(abs(elin).lt.9.e+9)) then

          elin=abs(rdat1(np1r))
          ergsev=1.602197e-12
          ener=ergsev*(12398.41)/max(elin,1.e-24)
          etst=ener/ergsev
          idest1=idat1(np1i)
          idest2=idat1(np1i+1)
          aij=rdat1(np1r+2)
c
c         get ion data
          nilin=npar(ml)
          mlm=nilin-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          do ktt=1,min(8,nkdt)
            write (kinam1(ktt:ktt),'(a1)')kdat1(np1k-1+ktt)
            enddo
          do ktt=nkdt+1,9
            write (kinam1(ktt:ktt),'(a1)')kblnk
            enddo
c
          j=ln
          elmtp=elum(1,j)
          elmtpb=elum(2,j)
          write (lun11,9924)j,elin,kinam1,
     $     elmtp,elmtpb,tau0(1,j), tau0(2,j)
9924      format (1h ,i9,1pe13.5,1x,a9,1x,4(1pe13.5))
c
          endif
        enddo
      write (lun11,993)
c
      return
c
c
 18   continue
c
      lpril=0
      write (lun11,*)'line wavelengths and levels'
c     step through lines
      nlpl=1
      do lnn=1,nlsvn
c
c       get line data
        ln=lnn
        ml=nplin(ln)
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        elin=abs(rdat1(np1r))
c
c       exclude rate type 14
        if ((lrtyp.ne.14).and.(abs(elin).gt.0.1)
     $       .and.(abs(elin).lt.9.e+9)) then
c
          ergsev=1.602197e-12
          ener=ergsev*(12398.41)/max(elin,1.e-24)
          etst=ener/ergsev
          idest1=idat1(np1i)
          idest2=idat1(np1i+1)
          aij=rdat1(np1r+2)
c
c         get ion data
          nilin=npar(ml)
          mlm=nilin-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          do ktt=1,min(8,nkdt)
            write (kinam1(ktt:ktt),'(a1)')kdat1(np1k-1+ktt)
            enddo
          do ktt=nkdt+1,9
            write (kinam1(ktt:ktt),'(a1)')kblnk
            enddo
c
          if (lpril.ge.1)
     $      write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
     $          (kdat1(np1ki+mm-1),mm=1,nkdti)
c
c         now find level data
          jkk=idat1(np1i+nidt-1)
          call func2l(jkk,lpril,lun11,t,xee,xpx,
     $              idat1,rdat1,kdat1,nptrs,
     $              npar,npnxt,npfi,
     $              rniss,rlev,ilv,
     $              nlpt,iltp,nlev,klev)

          ggup=rlev(2,idest1)
          gglo=rlev(2,idest2)
          do lk=1,20
            klablo(lk)=klev(lk,idest1)
            klabup(lk)=klev(lk,idest2)
            enddo
          flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo)
          ilevlo=idest1
          ilevup=idest2
c
          j=ln
          write (lun11,9929)j,elin,kinam1,
     $      (klev(mm,ilevlo),mm=1,20),(klev(mm,ilevup),mm=1,20),
     $      rlev(1,ilevlo),rlev(1,ilevup),rlev(2,ilevlo),rlev(2,ilevup),
     $      rlev(3,ilevlo),rlev(3,ilevup),
     $      ilv(1,ilevlo),ilv(1,ilevup),ilv(2,ilevlo),ilv(2,ilevup),
     $      ilv(3,ilevlo),ilv(3,ilevup)
 9929     format (1h ,i9,1pe13.5,1x,a9,1x,2(20a1,1x),6(1pe13.5),
     $          6i6)
c
          endif
        enddo
      write (lun11,993)
c
      return
c
c
 20   continue
c
      lpril=0
      write (lun11,*)'line finding list'
      do jlk=1,nlsvn
         j=jlk
         ml=nplin(j)
         mlm=ml-1
         call drd(ltyp,lrtyp,lcon,
     $     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $     nptrs,0,lun11)
         elin=abs(rdat1(np1r))
         jpnt(j)=j
         elsv(j)=abs(elin)
         enddo
c     sort
      done=.false.
      niter=0
      do while (.not.done)
        niter=niter+1
        done=.true.
c        do jjj=1,100
        do jjj=1,nlsvn-1
          j=jpnt(jjj)
          jp1=jpnt(jjj+1)
c          write (lun11,*)jjj,j,jp1,
c     $           elsv(jp1),elsv(j)
          if (elsv(jp1).lt.elsv(j)) then
            jpnt(jjj)=jp1
            jpnt(jjj+1)=j
            done=.false.
            endif
          enddo
        enddo
c
c     print out sorted list
      do jlk=1,nlsvn
        j=jpnt(jlk)
        lnn=j
c
c       get line data
        ln=lnn
        ml=nplin(ln)
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
c
c       exclude rate type 14
        if ((lrtyp.ne.14).and.(abs(elin).gt.0.1)
     $       .and.(abs(elin).lt.9.e+9)) then
c
          elin=abs(rdat1(np1r))
          ergsev=1.602197e-12
          ener=ergsev*(12398.41)/max(elin,1.e-24)
          etst=ener/ergsev
          idest1=idat1(np1i)
          idest2=idat1(np1i+1)
          aij=rdat1(np1r+2)
c
c         get ion data
          nilin=npar(ml)
          mlm=nilin-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          do ktt=1,min(8,nkdt)
            write (kinam1(ktt:ktt),'(a1)')kdat1(np1k-1+ktt)
            enddo
          do ktt=nkdt+1,9
            write (kinam1(ktt:ktt),'(a1)')kblnk
            enddo
c
          if (lpril.ge.1)
     $      write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
     $          (kdat1(np1ki+mm-1),mm=1,nkdti)
c
c         now find level data
          jkk=idat1(np1i+nidt-1)
          call func2l(jkk,lpril,lun11,t,xee,xpx,
     $              idat1,rdat1,kdat1,nptrs,
     $              npar,npnxt,npfi,
     $              rniss,rlev,ilv,
     $              nlpt,iltp,nlev,klev)

          ggup=rlev(2,idest1)
          gglo=rlev(2,idest2)
          do lk=1,20
            klablo(lk)=klev(lk,idest1)
            klabup(lk)=klev(lk,idest2)
            enddo
          flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo)
          ilevlo=idest1
          ilevup=idest2
c
          write (lun11,9929)j,elin,kinam1,
     $      (klev(mm,ilevlo),mm=1,20),(klev(mm,ilevup),mm=1,20),
     $      rlev(1,ilevlo),rlev(1,ilevup),rlev(2,ilevlo),rlev(2,ilevup),
     $      rlev(3,ilevlo),rlev(3,ilevup),
     $      ilv(1,ilevlo),ilv(1,ilevup),ilv(2,ilevlo),ilv(2,ilevup),
     $      ilv(3,ilevlo),ilv(3,ilevup)
           endif
         enddo
      write (lun11,993)
c
      return
c
 21   continue
c
      lpril=0
      write (lun11,*)' level opacities and emissivities'
      write (lun11,*)'index,energy,ion,level,index,emiss in,emiss out,th
     $reshold opacity,absorbed energy,depth in, depth out'
c
C     First look for element data (jk is element index)
      klel=11
      mlel=npfirst(klel)
      jk=0
      kk=0
      jkk=0
c
c     step through elements
      do while (mlel.ne.0)
c
c       get element data
        jk=jk+1
        mt2=mlel-1
        call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,
     $        nptrs,0,lun11)
        mllel=idat1(np1i+nidt-1)
        xeltp=rdat1(np1r)
        xeltp=abel(mllel)
        nnz=idat1(np1i)
        if (lpril.ge.1)
     $        write (lun11,*)'element:',jk,mlel,mllel,nnz,
     $                  (kdat1(np1k-1+mm),mm=1,nkdt)
c
C       ignore if the abundance is small
        if (xeltp.lt.1.e-10) then
            jkk=jkk+nnz
          else
c
c           now step thru ions (jkk is ion index)
            klion=12
            mlion=npfirst(klion)
            jkk=0
            kl=0
            do while ((mlion.ne.0).and.(kl.lt.nnz))
              jkk=jkk+1
c
C             retrieve ion name from kdati
              mlm=mlion-1
              call drd(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidti,np1i,nkdti,np1ki,mlm,
     $            nptrs,0,lun11)
              ethi=rdat1(np1r)
c
C             if not accessing the same element, skip to the next element
              mlleltp=idat1(np1i+nidti-2)
              if (mlleltp.eq.mllel) then
c
                kl=kl+1
                if (lpril.ge.1)
     $            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
     $                        (kdat1(np1ki+mm-1),mm=1,nkdti)
c
c               now find level data
                call func2l(jkk,lpril,lun11,t,xee,xpx,
     $              idat1,rdat1,kdat1,nptrs,
     $              npar,npnxt,npfi,
     $              rniss,rlev,ilv,
     $              nlpt,iltp,nlev,klev)
c
c               now step through rate type 7 data
                mltype=7
                ml=npfi(mltype,jkk)
                mllz=0
                if (ml.ne.0) mllz=npar(ml)
                mlpar=0
                if (ml.ne.0) mlpar=npar(ml)
                do while ((ml.ne.0).and.(mlpar.eq.mllz))
c
c                 get rrc data
                  kkkl=npconi2(ml)
                  if (lpril.ne.0) write (lun11,*)kkkl,ml,idest1,
     $                    elumab(1,kkkl),elumab(2,kkkl)
c
c                 test for non-zero rrc data
                  if ((kkkl.gt.0).and.(kkkl.le.ndat2)
     $                .and.((elumab(1,kkkl).gt.1.e-36)
     $                .or.(elumab(2,kkkl).gt.1.e-36))) then
c
c                   get rrc  data
                    mlm=ml-1
                    call drd(ltyp,lrtyp,lcon,
     $                nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $                nptrs,0,lun11)
                    idest1=idat1(np1i+nidt-2)
                    nlevp=nlev
                    idest2=nlevp+idat1(np1i-1+nidt-3)-1
c
c                   label for lower level
                    do lk=1,20
                      write (ktmp20(lk:lk),'(a1)')klev(lk,idest1)
                      enddo
                    klevl=ktmp20
c
c                   label for upper level
                    write (ktmp20(1:20),'(a20)')'continuum           '
                    klevu=ktmp20
c
c                   ion label
                    do lk=1,nkdti
                      write (ktmp8(lk:lk),'(a1)')kdat1(np1ki+lk-1)
                      enddo
                    do lk=nkdti+1,8
                      write (ktmp8(lk:lk),'(a1)')kblnk
                      enddo
c
                    eth=rlev(4,idest1)-rlev(1,idest1)
                    ett=eth
c
c                   get upper level data
                    if (idest2.gt.nlevp) then
                      jkk3=jkk+1
                      if (lpril.gt.1)
     $                  write (lun11,*)jkk3,ndtmp,nlevp,idest2
                      ndtmp=npfi(13,jkk3)
                      if (lpril.gt.1)
     $                  write (lun11,*)jkk3,ndtmp,nlevp,idest2
                      if (ndtmp.le.0) stop 'ndtmp error'
                      mllz=npar(ndtmp)
                      iltmp=0
                      do while ((ndtmp.ne.0).and.
     $                    (iltmp.ne.(idest2-nlevp+1)).and.
     $                    (npar(ndtmp).eq.mllz))
                        mlm=ndtmp-1
                        call drd(ltyp2,lrtyp2,lcon2,
     $                    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $                    nptrs,0,lun11)
                        iltmp=idat1(np1i2+nidt2-2)
                        if (lpril.gt.1) write (lun11,*)nidt2,iltmp,ndtmp
                        ndtmp=npnxt(ndtmp)
                        if (ndtmp.le.0) stop 'ndtmp error'
                        enddo
c                     NB fix to excited level PI and rec
                      ett=ett+rdat1(np1r2)
                      eth=ett
                      if (lpril.gt.1)
     $                  write (lun11,*) ndtmp,iltmp,idest2,ett
c                     label for lower level
                      ktmp20=kblnk20
                      do lk=1,nkdt2
                        write (ktmp20(lk:lk),'(a1)')kdat1(np1k2+lk-1)
                        enddo
                      klevu=ktmp20
                      endif
c
c                   other data
                    mmlv=npilev(idest1,jkk)
                    cabcompare=bremsint(nbinc(ethc,epi,ncn2))
                    mlcu=kkkl
                    write (lun11,969)kkkl,mmlv,ktmp8,idest1,idest2,
     $                  klevl,klevu,eth,
     $                  cemab(1,mlcu),cemab(2,mlcu),opakab(mlcu),
     $                  cabab(mlcu),tauc(1,mlcu),tauc(2,mlcu)
 969                format (1x,2i6,1x,a8,1x,2i6,1x,2(a20,1x),
     $                  8(1pe13.5),2i6)
c
c                   done with this rrc
                    endif
c
c                 end of loop over rrcs
                  ml=npnxt(ml)
                  if (ml.ne.0) mlpar=npar(ml)
                  enddo
c
c               end of test for element
                endif
c
C             Go to next ion
              mlion=npnxt(mlion)
              enddo
c
c         end of test for non-zero element abund
          endif
c
        mlel=npnxt(mlel)
C       Go to next element
        enddo
c
      write (lun11,993)
c
      return
C
c
7     continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Write level populations
c
      write (lun11,9985)
9985  format (1x,' level populations ')
      write (lun11,9986)
 9986 format (1x,' ion                      level              '
     $,' e_exc population')

C     lpril is flag for printing debug information
      lpril=0
      if (lpril.ne.0) then
        write (lun11,*)'raw data'
        do j=1,nnml
          if (xilev(j).gt.1.e-37)
     $     write (lun11,*)j,xilev(j),elumab(1,j)
          enddo
        endif
c
C     First look for element data (jk is element index)
      klel=11
      mlel=npfirst(klel)
      jk=0
c
c     step through elements
      do while (mlel.ne.0)
c
c       get element data
        jk=jk+1
        mt2=mlel-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,
     $    nptrs,0,lun11)
        mllel=idat1(np1i+nidt-1)
        nnz=idat1(np1i)
        xeltp=rdat1(np1r)
        xeltp=abel(mllel)
        if (lpril.ne.0)
     $        write (lun11,*)'element:',jk,mlel,mllel,nnz,
     $                    (kdat1(np1k-1+mm),mm=1,nkdt),xeltp
c
C       ignore if the abundance is small
        if (xeltp.lt.1.e-10) then
            jkk=jkk+nnz
          else
c
c           now step thru ions (jkk is ion index)
            klion=12
            mlion=npfirst(klion)
            jkk=0
            kl=0
            do while ((mlion.ne.0).and.(kl.lt.nnz))
c
              jkk=jkk+1
C             retrieve ion name from kdati
              mlm=mlion-1
              call drd(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidt,np1i,nkdti,np1ki,mlm,
     $            nptrs,0,lun11)
c
C             if not accessing the same element, skip to the next element
              mlleltp=idat1(np1i+nidt-2)
              if (mlleltp.eq.mllel) then
c
                kl=kl+1
                if (lpril.ne.0)
     $            write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
     $                        (kdat1(np1ki+mm-1),mm=1,nkdti)
                do ktt=1,min(8,nkdti)
                  write (kinam1(ktt:ktt),'(a1)')kdat1(np1ki-1+ktt)
                  enddo
                do ktt=nkdti+1,9
                  write (kinam1(ktt:ktt),'(a1)')kblnk
                  enddo
c
c               get level data
                call func2l(jkk,lpril,lun11,t,xee,xpx,
     $              idat1,rdat1,kdat1,nptrs,
     $              npar,npnxt,npfi,
     $              rniss,rlev,ilv,
     $              nlpt,iltp,nlev,klev)
c
c               step thru levels
                do mm2=1,nlev
c
c                 get level pointer
                  mmtmp=npilev(mm2,jkk)
                  if (mmtmp.ne.0) then
                    kkkl=mmtmp
                    mmlv=mmtmp
c
c                   test for level pop
                    if (xilev(kkkl).gt.1.d-64) then
c
c                     get data
                      eth=rlev(1,mm2)
                      dep=xilev(kkkl)/(rniss(kkkl)+1.e-34)
                      write (lun11,9296)kkkl,kinam1,
     $                   (klev(lk,mm2),lk=1,20),eth,xilev(kkkl),
     $                   rniss(kkkl),dep
 9296                 format (1x,i6,1x,a8,1x,(20a1),7(1pe13.5))
c

c                     end of test for level pop
                      endif

c                   end of test for level pointer
                    endif
c
c                 end of step thru levels
                  enddo
c
c               end of test for element
                endif
c
C             Go to next ion
              mlion=npnxt(mlion)
              enddo
c
C           end of test for abundance
            endif
c
C       Go to next element
        if (mlel.ne.0) mlel=npnxt(mlel)
        enddo
c
      write (lun11,*)'done with 7'
      write (lun11,993)
c
c
      return
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Print a short summary line of the radial calculation
C
9     continue
c
      elsum=0.
      ergsev=1.602197e-12
      do jlk=1,nlsvn
         ln=jlk
         ml=nplin(ln)
         mlm=ml-1
         call drd(ltyp,lrtyp,lcon,
     $     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $     nptrs,0,lun11)
         elin=abs(rdat1(np1r))
         if ((elin.lt.1.e+8).and.(elin.gt.1.)) then
           elsum=elsum+elum(1,ln)+elum(2,ln)
           endif
         enddo
      sumtmp1=0.
      sumtmp2=0.
      ergsev=1.602197e-12
      r19=r*1.e-19
      tmp1=zremsz(1)
      tmp2=zrems(1,1)
      do jk=2,ncn2
         tmp1o=tmp1
         tmp1=zremsz(jk)
         sumtmp1=sumtmp1+(tmp1+tmp1o)*(epi(jk)-epi(jk-1))*ergsev/2.
         tmp2o=tmp2
         tmp2=zrems(1,jk)
         sumtmp2=sumtmp2+(tmp2+tmp2o)*(epi(jk)-epi(jk-1))*ergsev/2.
         enddo
c      terr=(sumtmp1-sumtmp2-elsum)/(sumtmp1+1.e-24)
      terr=(sumtmp1-sumtmp2)/(sumtmp1+1.e-24)
      uu1=enlum/(12.56*xpx*r19*r19)/3.e+10
      alguu1=log10(max(1.e-24,uu1))
      skse=xlum/(xpx*r19*r19)
      zeta=log10(max(1.e-24,skse))
      ecc=2.998e+10
      ekt=t*(0.861707)*ergsev
      sksec=skse/12.56/((1.+xee)*ekt*ecc)
      zetac=log10(max(1.e-24,sksec))
      enn0=xpx
      nlyc=nbinc(13.7d0,epi,ncn2)
      nry=nlyc+1
      egam=zremsz(nry)/(2.*12.56*enn0*ecc*r19*r19+1.e-24)
      nry=nbinc(13.6d0,epi,ncn2)+1
c      write (lun11,9969)r,rdel,zeta,xee,xpx,t,hmctot,
c     $ dpthc(1,nry),dpthc(2,nry),ntotit,lnerrd
c9969  format (1x,9(1pe10.3),2i3)
      tmp1=log10(r)
      tmp2=log10(max(1.e-36,min(99.,rdel/r)))
      tmp2c=log10(max(xcol,1.e-10))
      tmp3=log10(xpx)
      tmp4=log10(t)+4.
      tmp5=log10(max(dpthc(1,nry),1.e-10))
      tmp6=log10(max(dpthc(2,nry),1.e-10))
      tmp7=min(99.99,max(-99.99,hmctot*100.))
      tmp8=min(99.99,max(-99.99,terr*100.))
      write (tmpst,9889)tmp1,tmp2,tmp2c,zeta,xee,tmp3,tmp4,tmp7,
     $ tmp8,tmp5,tmp6,ntotit
      write (lun11,9889)tmp1,tmp2,tmp2c,zeta,xee,tmp3,tmp4,tmp7,
     $ tmp8,tmp5,tmp6,ntotit
 9889  format (1x,11(1x,f6.2),2i3)
      call xwrite(tmpst,10)
c
      return
c
c
 16   continue
c
c     times
!      write (lun11,*)'times:',tread,tloop,tfunc,trates1,thcor,trates2,    !jg
!     $          theat
      ttot=0.
      do ll=1,ntyp
!        ttmpi=tucalc(ll)/max(1,ncall(ll))   !jg
!        ttot=ttot+tucalc(ll)   !jg
!        write (lun11,9892)ll,ncall(ll),tucalc(ll),ttmpi   !jg
! 9892   format (1x,2i8,2(1pe11.3))
        enddo
      write (lun11,*)'total ucalc=',ttot
      write (lun11,993)
c
      return
c
c
 17   continue
c
c     column headings
      klabs(1)='log(r)'
      klabs(2)='delr/r'
      klabs(3)='log(N)'
      klabs(4)='log(xi)'
      klabs(5)=' x_e  '
      klabs(6)='log(n)'
      klabs(7)='log(t)'
      klabs(8)='h-c(%)'
      klabs(9)='h-c(%)'
      klabs(10)='log(tau)'
c      klabs(10)='ntotit'
      write (lun11,9979)(klabs(mm),mm=1,10)
      write (tmpst,9979)(klabs(mm),mm=1,10)
9979  format (2x,3(1x,a6),1x,a7,a6,4(1x,a6),(1x,a9))
      call xwrite(tmpst,10)
      klabs(1)='      '
      klabs(2)='      '
      klabs(3)='      '
      klabs(4)='      '
      klabs(5)='      '
      klabs(6)='      '
      klabs(7)='      '
      klabs(8)='      '
      klabs(9)='      '
      klabs(10)='fwd   '
      klabs(11)='rev   '
      write (lun11,9989)(klabs(mm),mm=1,11)
      write (tmpst,9989)(klabs(mm),mm=1,11)
      call xwrite(tmpst,10)
9989  format (3x,11a7)
c
      return
c
c
 22   continue
c
c     ionization parameter etc.
      rdum=delr
      delr=rdum
      ergsev=1.602197e-12
      r19=r*(1.e-19)
c      write (lun11,*)enlum,xpx,r,xlum
      uu1=enlum/(12.56*xpx*r19*r19)/3.e+10
c      write (lun11,*)uu1
      enlumx=0.
      nb1=nbinc(100d0,epi,ncn2)
      nb10=nbinc(10000.d0,epi,ncn2)
c      write (lun11,*)'nb1=',nb1,nb10
      do kl=nb1,nb10
c        write (lun11,*)kl,epi(kl),zremsz(kl),enlumx
        enlumx=enlumx+(zremsz(kl)/epi(kl)+zremsz(kl-1)/epi(kl-1))
     $                *(epi(kl)-epi(kl-1))/2.
        enddo
      uux=enlumx/(12.56*xpx*r19*r19)/3.e+10
      alguux=log10(max(1.e-24,uux))
      alguu1=log10(max(1.e-24,uu1))
      skse=xlum/(xpx*r19*r19)
      zeta=log10(max(1.e-24,skse))
      ecc=2.998e+10
      ekt=t*(0.861707)*ergsev
c      sksec=skse/(12.56*((1.+xee)*ekt+pradl/(1.e-24+xpx))*ecc)
      sksec=skse/12.56/((1.+xee)*ekt*ecc)
      zetac=log10(max(1.e-24,sksec))
      enn0=xpx
      nlyc=nbinc(13.7d0,epi,ncn2)
      nry=nlyc+1
      egam=zremsz(nry)/(2.*12.56*enn0*ecc*r19*r19+1.e-24)
      nry=nbinc(13.6d0,epi,ncn2)+1
9968  format (1x,' log(Xi)=',1pe11.3, ' log(u1)=',1pe11.3,
     $ ' log(ux)=',1pe11.3,' gamma=',1pe11.3, ' rdel=',1pe11.3)
 9965 format (1x,' r=',1pe11.3,' t=',1pe11.3,' log(xi)=',1pe11.3,
     $ ' n_e=',1pe11.3,' n_p=',1pe11.3)
9966  format (1x,'httot=',1pe11.3,' cltot=',1pe11.3,
     $      'taulc=',1pe11.3,'taulcb=',1pe11.3)
      write(lun11,9965)r,t,zeta,xnx,xpx
      write(lun11,9966)httot,cltot,dpthc(1,nry),dpthc(2,nry)
      write(lun11,9968)zetac,alguu1,alguux,egam,rdel
      write (lun11,993)
c
      return
c
 25   continue
c
c      write (lun11,*)'outputting to the common block',nlsvn
c      do mm=1,ncn2
c        epi2(mm)=epi(mm)
c        do ll=1,3
c          zrems2(ll,mm)=zrems(ll,mm)
c          enddo
c        enddo
c      lpril=1
c      nilino=0
c      jkktmp=0
c      do j=1,nlsvn
c          kk=j
c          ln=nplin(j)
c          ml=ln
c          if (ml.ne.0) then
cc            write (lun11,*)'   ',j,ml
c            mlm=ml-1
c            call drd(ltyp,lrtyp,lcon,
c     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
c     $        idat1,rdat1,kdat1,nptrs,0,lun11)
c            elin=rdat1(np1r)
c            llo=idat1(np1i)
c            lup=idat1(np1i+1)
c            elin=rdat1(np1r)
c            aij=rdat1(np1r+2)
c            nilin=npar(ml)
c            if ((nilin.gt.0).and.(nilin.lt.ndat2)) then
c                if (nilin.ne.nilino) jkktmp=jkktmp+1
c                nilino=nilin
c                mlm=nilin-1
c                call drd(ltyp,lrtyp,lcon,
c     $            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
c     $            idat1,rdat1,kdat1,nptrs,0,lun11)
c                do mm=1,nkdt
c                  kdtmp(mm)=kdat1(np1k-1+mm)
c                  enddo
c                do mm=nkdt+1,9
c                  kdtmp(mm)=kblnk
c                  enddo
c                nilin=idat1(np1i+2)
cc                write (lun11,*)ml,nilin,npar(ml)
c                newout=newout+1
c                newout=min(newout,nnnl)
c                lnewo(newout)=j
c                ml=npfi(13,jkktmp)
c                mllz=npar(ml)
c                lupfnd=0
c                llofnd=0
c                mlpar=npar(ml)
c                do while ((ml.ne.0).and.(mlpar.eq.mllz)
c     $            .and.((llofnd.ne.1).or.(lupfnd.ne.1)))
cc                    write (lun11,*)ml,nptrs(2,ml),mltype,jkk
c                  mlm=ml-1
c                  call drd(ltyp,lrtyp,lcon,
c     $              nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
c     $              idat1,rdat1,kdat1,nptrs,0,lun11)
c                  nlevmx=nlevmx+1
c                  nlev=idat1(np1i+nidt-2)
cc                  write (lun11,*)ml,nlev,llo,lup,(kdat1(np1k-1+mm),mm=1,nkdt)
c                  if (nlev.eq.llo) then
c                    do mm=1,20
c                      if (mm.le.nkdt) then
c                          klevl(mm)=kdat1(np1k-1+mm)
c                        else
c                          klevl(mm)=kblnk
c                        endif
c                      enddo
cc                   write (lun11,*)kk,ktmp2
c                    llofnd=1
c                    gglo=rdat1(np1r+1)
c                    endif
c                  if (nlev.eq.lup) then
c                    do mm=1,20
c                      if (mm.le.nkdt) then
c                          klevu(mm)=kdat1(np1k-1+mm)
c                        else
c                          klevu(mm)=kblnk
c                        endif
c                      enddo
c                    lupfnd=1
c                    ggup=rdat1(np1r+1)
c                    endif
c                  ml=npnxt(ml)
c                  if (ml.ne.0) mlpar=npar(ml)
c                  enddo
c                if ((llofnd.eq.1).and.(lupfnd.eq.1)) then
c                  flinewo(newout)=(1.e-16)*aij*ggup*elin*elin
c     $                             /((0.667274)*gglo)
c                  aijewo(newout)=aij
c                  ggloewo(newout)=gglo
c                  ggupewo(newout)=ggup
c                  do mm=1,8
c                    kdewo(mm,newout)=kdtmp(mm)
c                    enddo
c                  do mm=1,20
c                    kdewol(mm,newout)=klevl(mm)
c                    enddo
c                  do mm=1,20
c                    kdewou(mm,newout)=klevu(mm)
c                    enddo
c                  elewo(newout)=elin
c                  tau0ewo(newout)=tau0(1,j)
c                  elout(1,newout)=elum(1,j)
c                  elout(2,newout)=elum(2,j)
cc                  write (lun11,*)kk,ln,j,(kdtmp(mm),mm=1,8),elin,
cc     $             tau0(1,j),elum(1,j),elum(2,j),newout
cc9955             format (1x,2i8,1x,8a1,3(1pe11.3))
c                  endif
c              endif
c            endif
c          enddo
c      call commonprint(lun11)
c
      return
c
 26   continue
c
      return
c
c     ferland print
      lpril=0
c     print 500 strongest emission lines
      write (lun11,*)'log(emission line fluxes (erg/sec/cm^2))'
      kltmpo=0
      nlplmx=500
      eliml=0.1
      elimh=1.0e10
c     find the strongest lines.
      do  lm=1,nlplmx
        kltmp(lm)=0
        enddo
      nlpl=1
      do lnn=1,nlsvn
        ln=lnn
        ml=nplin(ln)
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        elin=abs(rdat1(np1r))
        if (lrtyp.ne.14) then
          nilin=npar(ml)
          mlm=nilin-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          nilin2=idat1(np1i+nidt-1)
          elmmtpp=(elum(2,ln)+elum(1,ln))/2.
          if (lpril.ne.0)
     $       write (lun11,*)lnn,elin,nilin,elmmtpp,ln,ml
          if ((ln.gt.0).and.(ln.lt.nnnl)
     $       .and.(elin.ge.eliml).and.(elin.le.elimh)
     $       .and.(elin.le.8.9e+6)
     $       .and.(elmmtpp.gt.1.e-36)
     $       .and.(nilin2.gt.0).and.(nilin2.le.nni))
     $        then
            lmm=0
            elcomp=1.e+10
            do while ((lmm.lt.nlpl).and.(elmmtpp.lt.elcomp))
              lmm=lmm+1
              kl2=kltmp(lmm)
              elcomp=0.
              if (kl2.gt.0)
     $          elcomp=(elum(2,kl2)+elum(1,kl2))/2.
              enddo
            if (lpril.ne.0)
     $       write (lun11,8516)ln,elin,elmmtpp,lmm,nlpl,kl2,elcomp
c 8516       format (1h ,i4,2e12.4,3i4,e12.4)
            kltmpo=ln
            do  k=lmm,min(nlplmx,nlpl)
              if ((lpril.ne.0).and.(kltmp(k).ne.0))
     $          write (lun11,*)'in 557 loop',k,kltmp(k),kltmpo
              kltmpn=kltmp(k)
              kltmp(k)=kltmpo
              kltmpo=kltmpn
              enddo
           nlpl=min(nlplmx,nlpl+1)
           if (lpril.ne.0)
     $       write (lun11,*)'done with 557 loop',lm
            endif
          endif
        enddo
       if (nlpl.gt.0) kltmp(nlpl)=kltmpo
c      nlpl=nlpl-1
      write (lun11,9599)
 9599 format (1x,'index, ion, wavelength, transmitted, reflected,total')
      r19=r*1.e-19
      do  kk=1,nlpl
        if (lpril.ne.0)
     $    write (lun11,*)'kk=',kk
        ln=kltmp(kk)
        if (ln.ne.0) then
          ml=nplin(ln)
          if (ml.ne.0) then
            if (lpril.ne.0)
     $      write (lun11,*)'   ',ln,ml
            mlm=ml-1
            call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $        nptrs,0,lun11)
            elin=abs(rdat1(np1r))
            nilin=npar(ml)
            mlm=nilin-1
            call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $        nptrs,0,lun11)
            do mm=1,nkdt
              kdtmp(mm)=kdat1(np1k-1+mm)
              enddo
            do mm=nkdt+1,9
              kdtmp(mm)=kblnk
              enddo
             flux1=elum(1,ln)/12.56/r19/r19
             flux2=elum(2,ln)/12.56/r19/r19
c            nilin=idat1(np1i+2)
            if (lpril.ne.0)
     $      write (lun11,*)ml,nilin,npar(ml)
            write (lun11,9956)kk,ln,(kdtmp(mm),mm=1,9),elin,
     $      log10(flux1),log10(flux2),log10(flux1+flux2)
 9956       format (1x,2i8,1x,9a1,4(1pe13.5))
            endif
          endif
        enddo
      write (lun11,993)
      return
c

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Write  ion column densities
c     requires that zrtmp be filled by calling pprint(12)
C
 27   continue
c
      write (lun11,*)'ion column densities'
      write (lun11,9447)
9447   format (1x,'index, ion, column density')
c
      do lk=1,nni
        xcoltmp(lk)=0.
        enddo
c
c     step thru ions
      klion=12
      mlion=npfirst(klion)
      lk=0
      do while (mlion.ne.0)
c
c        get ion data
         lk=lk+1
         ltyp=klion
         mlm=mlion-1
         call drd(ltyp,lrtyp,lcon,
     $     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $     nptrs,0,lun11)
         do mm=1,nkdt
           kdtmp(mm)=kdat1(np1k-1+mm)
           enddo
         do mm=nkdt+1,9
           kdtmp(mm)=kblnk
           enddo
c
c        get element data
         nell=npar(mlion)
         mlm=nell-1
         call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
c         write (lun11,*)mlion,lk,np1i,nidt,np1i+nidt-1,
c     $      idat1(np1i+nidt-1),ababs(idat1(np1i+nidt-1)),mlm
         if ((idat1(np1i+nidt-1).gt.0)
     $     .and.(idat1(np1i+nidt-1).le.nl)) then
           abundel=ababs(idat1(np1i+nidt-1))
           xeltp=abundel
c
           do jkl=2,numrec
c             write (lun11,*)jkl,lk,zrtmp(2,jkl),zrtmp(8+lk,jkl),
c     $                             zrtmp(5,jkl),xeltp,xcoltmp(lk)
             xcoltmp(lk)=xcoltmp(lk)
     $         +(zrtmp(8+lk,jkl)*zrtmp(5,jkl)
     $             +zrtmp(8+lk,jkl-1)*zrtmp(5,jkl-1))
     $         *(zrtmp(2,jkl)-zrtmp(2,jkl-1))*xeltp/2.
             enddo
c
c          print out
           if (xcoltmp(lk).gt.1.e-15)
     $      write (lun11,9446)lk,(kdtmp(mm),mm=1,9),
     $      xcoltmp(lk)
9446       format (1x,i4,1x,9a1,1pe16.8)
c
           endif
c
         mlion=npnxt(mlion)
         enddo
c
c
      return
      end
      subroutine printerror(lun11,status)
      implicit none

c     print out the fitsio error messages to the user
c     author:  T. Bridgman
c
      integer status, lun11
      character errtext*30,errmessage*80
c
c     check if status is ok (no error); if so, simply return
      if (status .le. 0)return
c
c     get the text string which describes the error
      call ftgerr(status,errtext)
      write (lun11,*)'fitsio error status =',status,': ',errtext

c     read and print out all the error messages on the fitsio stack
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          write (lun11,*)errmessage
          call ftgmsg(errmessage)
      end do
c
      end
      subroutine rdflo(flo,qry,ios,lpri,ind)

c       reads integer from file param corresponding to character string
c       'char'. if no integer there, prompts user for parameter ...
c      author:  T. Kallman
c
       implicit none

       real*8 sum, float, flo, sum2, flo2

       integer lun11, lpri, ind, kk, lfnde, lfndd
       integer lnon, lfndm, lfnd2, lfnd, lfnd2o, ll
       integer lfndb, lneg, ind2, idec, iexp, itmp, ios

       character(1) qtmp
       character(72) qry
       character(1) chtst(16)
c
       data chtst/'0','1','2','3','4','5','6','7','8','9','+','-',
     $            ' ','.','e','e'/
c
       lun11=6
       if (lpri.gt.2) write (lun11,*) 'in rdflo',ind
       if (lpri.gt.2) write (lun11,*)qry
c
c      scan for e and dot
c       kk=ind
       kk=0
       lfnde=0
       lfndd=0
       lnon=0
       lfndm=0
       lfnd2=0
102    kk=kk+1
       if (kk.ge.72) return
       read (qry(kk:kk),'(a1)') qtmp
       if (qtmp.eq.'.') lfndd=kk
       if ((qtmp.eq.'-').and.(lfndd.eq.0)) lfndm=kk
       if ((qtmp.eq.'e').or.(qtmp.eq.'e')) lfnde=kk
       lfnd=0
       lfnd2o=lfnd2
       lfnd2=0
       do 3011 ll=1,16
         if (qtmp.eq.chtst(ll)) lfnd=1
         if ((qtmp.eq.chtst(ll)).and.(ll.ne.13)) lfnd2=1
 3011    continue
       if ((lfnd2.eq.1).and.(lfnd2o.eq.0))
     $      lfndb=kk
       if (lfnd.eq.0) lnon=1
       if (kk.lt.ind) go to 102
       if (lpri.gt.2) write (lun11,9905)lfnde,lfndd
9905   format (1x,' e and d switches',2i4)
       if (lfndd.eq.0) go to 10
       if (lnon.ne.0) go to 1011
       lneg=1
       if (lfndm.ne.0) lneg=-1
       read (qry(lfndb:lfndb),'(a1)') qtmp
       if ((qtmp.eq.'+').or.(qtmp.eq.'-'))lfndb=lfndb+1
c
c      scan for mantissa
       ind2=ind
       kk=lfndb
       idec=lfndd
       if (lpri.gt.2) write (lun11,*)'scanning mantissa'
       if (lfnde.ne.0) ind2=lfnde-1
       if (lpri.gt.2) write (lun11,*)'ind,ind2,lfnde,lfndb:',
     $    ind,ind2,lfnde,lfndb
       sum=0.
       iexp=1-(kk-idec+1)
       kk=kk-1
104       kk=kk+1
          if (kk.eq.idec) go to 301
          iexp=iexp-1
          read (qry(kk:kk),'(i1)') itmp
          sum=sum+float(itmp)*10.**iexp
301       continue
          if (lpri.gt.2) write (lun11,9907)kk,itmp,iexp,sum
9907      format (1x,' kk,itmp,iexp,sum ',3i8,e12.4)
          if (kk.lt.ind2) go to 104
       flo=sum*float(lneg)
       ios=0
       if (lfnde.eq.0) return
c
c      scan for exponent
       if (lpri.gt.2) write (lun11,*)'scanning exponent'
       kk=lfnde+1
       lneg=1
       sum2=0.
       read (qry(kk:kk),'(a1)') qtmp
       if (qtmp.ne.'-') go to 1104
          kk=kk+1
          lneg=-1
 1104     continue
       if (qtmp.eq.'+') kk=kk+1
       iexp=1-(kk-ind)
       kk=kk-1
 105      kk=kk+1
          iexp=iexp-1
          read (qry(kk:kk),'(i1)') itmp
          sum2=sum2+float(itmp)*10.**iexp
          if (lpri.gt.2) write (lun11,9907)kk,itmp,iexp,sum2
          if (kk.lt.ind) go to 105
       flo2=sum2*float(lneg)
       flo=flo*10.**flo2
       if (lpri.gt.2) write (lun11,*)'returning:',flo2,flo
       ios=0

       ios=0
       return
c
 10    continue
       ios=-1
       return
c
 1011  continue
       ios=999
c
       return
       end
       subroutine rdint(lo,qry,ios,lpri,ind)

c       reads integer from file param corresponding to character string
c       'char'. if no integer there, prompts user for parameter ...
c      author:  T. Kallman
c
       implicit none

       integer lun11, lpri, ind, kk
       integer lnon, lfndm, lfnd2, lfnd, lfnd2o, ll
       integer lfndb, lneg, ind2, idec, iexp, itmp, ios
       integer isum, lo

       character(1) qtmp
       character(72) qry
       character(1) chtst(15)
c
       data chtst/'0','1','2','3','4','5','6','7','8','9','+','-',
     $           ' ','.','e'/
c
       lun11=6
       if (lpri.gt.2) write (lun11,*) 'in rdint',ind
       if (lpri.gt.2) write (lun11,*)qry
c
c      scan for e and dot
       kk=0
       lfndm=0
       lnon=0
       lfnd2=0
102    kk=kk+1
       if (kk.ge.72) return
       read (qry(kk:kk),'(a1)') qtmp
       if (qtmp.eq.'-') lfndm=kk
       lfnd=0
       lfnd2o=lfnd2
       lfnd2=0
       do 3011 ll=1,13
         if (qtmp.eq.chtst(ll)) lfnd=1
         if ((qtmp.eq.chtst(ll)).and.(ll.le.10)) lfnd2=1
 3011    continue

       lfndb=0
       if ((lfnd2.eq.1).and.(lfnd2o.eq.0))
     $      lfndb=kk
       if (lfnd.eq.0) lnon=1
       if (kk.lt.ind) go to 102
       if (lnon.ne.0) go to 1011
       lneg=1
       if (lfndm.ne.0) lneg=-1
c
c      scan for mantissa
       kk=lfndb
       ind2=ind
       idec=ind+1
       if (lpri.gt.2) write (lun11,*)'scanning mantissa'
       isum=0
       iexp=-1-(kk-idec-1)
       kk=kk-1
104       kk=kk+1
          iexp=iexp-1
          read (qry(kk:kk),'(i1)') itmp
          isum=isum+itmp*10**iexp
          if (lpri.gt.2) write (lun11,9907)kk,itmp,iexp,isum
9907      format (1x,' kk,itmp,iexp,sum ',4i8)
          if (kk.lt.ind2) go to 104
       lo=isum*lneg
       ios=0
       return
c
c 10    continue
c       ios=-1
c       return
c
 1011  continue
       ios=999
c
       return
       end
      subroutine readtbl(nptrs,np1r,np1i,np1k,np2,npdat2,
     &           rdat1,idat1,kdat1,nidat1,filename,credate,lpri,lun11)
c
c     reads in atomic data
c       written by Ke Zhang, Nov. 9, 2001
c
      implicit none
c
      integer nidatt
      parameter (nidatt=35000000)
c
      integer nidat1,npdat2,lpri,lun11,np1r,np1i,np1k,np2
      integer idat1(nidat1)
      real rdat14(nidatt)
      real*8 rdat1(nidat1)
      character(1) kdat1(nidat1)
      integer nptrs(10,npdat2),ntptr(nidatt),mm

      real nulle
      integer nullj
      character filename*256,nullstr*1
      character credate*63, comment*50
      logical anynull

      integer status,unit,readwrite,blocksize,hdutype
      integer row,col,j,i,lenact

      status=0
      nullstr=' '
      nullj=0
      nulle=0.

      if (lpri.ne.0) write (lun11,*)'in readtbl'
c
c get an unused logical unit number and open the fits file
c      call ftgiou(unit,status)
      call getlun(unit)
      readwrite=0
c     print *, 'in readtbl:'
c     print *, '   ', unit
      call ftopen(unit,filename,readwrite,blocksize,status)
c
c     print *, '   ', status
c

c Read the primary header & get the file creation date
      call ftmahd(unit,1,hdutype,status)
      call ftgkys(unit,'DATE',credate,comment,status)
      if(status .gt.0) call printerror(lun11,status)
c     write(lun11,*)'Atomic Data Version: ',credate(1:lenact(credate))

      col=1
      row=1

c move to the next extension : POINTERS
      call ftmrhd(unit,1,hdutype,status)

c read LENGTH keywords: total records #
      call ftgkyj(unit,'LENGTH',np2,comment,status)

c read POINTERS data
      call ftgcvj(unit,col,row,1,10*np2,nullj,
     &               ntptr,anynull,status)
      do i=1,np2
        do j=1,10
          nptrs(j,i)=ntptr(10*(i-1)+j)
        enddo
      enddo


c move to the next extension : REALS
      call ftmrhd(unit,1,hdutype,status)

c read LENGTH keywords: total reals #
      call ftgkyj(unit,'LENGTH',np1r,comment,status)

c read REAL data
      call ftgcve(unit,col,row,1,np1r,nulle,
     &               rdat14,anynull,status)
      do mm=1,nidat1
        rdat1(mm)=rdat14(mm)
        enddo


c move to the next extension : INTEGERS
      call ftmrhd(unit,1,hdutype,status)

c read LENGTH keywords: total integers #
      call ftgkyj(unit,'LENGTH',np1i,comment,status)

c read INTEGER data
      call ftgcvj(unit,col,row,1,np1i,nullj,
     &               idat1,anynull,status)


c move to the next extension : CHARS
      call ftmrhd(unit,1,hdutype,status)

c read LENGTH keywords: total chars #
      call ftgkyj(unit,'LENGTH',np1k,comment,status)

c read CHAR data
      call ftgcvb(unit,col,row,1,np1k,nullstr,
     &             kdat1,anynull,status)


c close the file and free the unit number
      call ftclos(unit, status)
c      call ftfiou(unit, status)
      close(unit)

c check for any error, and if so print out error messages
      if (status .gt. 0) call printerror(lun11,status)

      return
      end
      subroutine remtms(ct)
      implicit none
c
      real*8 ct
      real a(2)
c     real etime
c
c      ct=etime(a)
       ct=0.
c      write (6,*)'in remtms:',ct,a
c
      return
      end
      subroutine rread1(trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,numrec,npass,
     $ nlimd,rmax,xpxcol,xi,zeta,lfix,
     $ lun11,abel,cfrac,emult,taumax,xeemin,spectype,specfile,specunit,
     $ kmodelname,nloopctl,critf,vturbi,eptmp,zrtmp,numcon2,ncn2,radexp)
c
c     this routine handles reading of the input data.
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
      real*8 eptmp(ncn),zrtmp(ncn),abel(nl),abel2(30)
      character(8) stringst,kblnk8
      character(80) specfile,spectype,stringsl,kblnk80,stringst2
      character(30) kmodelname
      integer nloopctl,specunit,ierr,ll,lcdd2,lun13,nenergy,ncn2
      integer lwri,lpri,lcdd,numrec,npass,nlimd,lfix,lun11,numcon2,mm
      real*8 trad,xlum,r,t,xpx,p,rmax,xpxcol,xi,zeta,cfrac,emult,taumax,
     $     xeemin,critf,vturbi,ccc,xlum2,xpcol,r19,radexp
c
      data kblnk8/'        '/
      data kblnk80/
     $'
     $               '/
c
            ierr=0
            call uclgsr8('cfrac',cfrac,ierr)
c
c           temperature
            call uclgsr8('temperature',t,ierr)
c
c           pressure/density switch
            call uclgsi('lcpres',lcdd2,ierr)
            lcdd=1-lcdd2
c
c           pressure
            call uclgsr8('pressure',p,ierr)
c
c           density
            call uclgsr8('density',xpx,ierr)
c
c           spectrum
            specfile=kblnk80
            spectype=kblnk80
            stringst2=kblnk80
            specunit=0
            stringst=kblnk8
            call uclgst('spectrum',spectype,ierr)
            xlum=1.
            read (spectype(1:8),'(a8)')stringst
            if (stringst.eq.'file    ') then
              specfile=kblnk80
              call uclgst('spectrum_file',specfile,ierr)
              read (specfile(1:80),'(a80)')stringst2
              call getlun(lun13)
              open (unit=lun13,file=stringst2)
              call uclgsi('spectun',specunit,ierr)
              read (lun13,*)nenergy
              numcon2 = nenergy
              do ll=1,nenergy
                read (lun13,*)eptmp(ll),zrtmp(ll)
                if (specunit.eq.1) zrtmp(ll)=zrtmp(ll)*eptmp(ll)
                if (specunit.eq.2) zrtmp(ll)=10.**zrtmp(ll)
                enddo
              endif
c
c           trad
            call uclgsr8('trad',trad,ierr)
c
c           luminosity
            call uclgsr8('rlrad38',xlum,ierr)
c
c           column density
            call uclgsr8('column',xpcol,ierr)
c
c           ionization parameter
            call uclgsr8('rlogxi',zeta,ierr)
c
c           number of steps
            call uclgsi('nsteps',numrec,ierr)
c
c           number of iterations
            call uclgsi('niter',nlimd,ierr)
c
c           write switch
            call uclgsi('lwrite',lwri,ierr)
c
c           print switch
            call uclgsi('lprint',lpri,ierr)
c
c           step size choice
            call uclgsi('lstep',lfix,ierr)
c
c           abundances
            call uclgsr8('habund',abel2(1),ierr)
            call uclgsr8('heabund',abel2(2),ierr)
            call uclgsr8('liabund',abel2(3),ierr)
            call uclgsr8('beabund',abel2(4),ierr)
            call uclgsr8('babund',abel2(5),ierr)
            call uclgsr8('cabund',abel2(6),ierr)
            call uclgsr8('nabund',abel2(7),ierr)
            call uclgsr8('oabund',abel2(8),ierr)
            call uclgsr8('fabund',abel2(9),ierr)
            call uclgsr8('neabund',abel2(10),ierr)
            call uclgsr8('naabund',abel2(11),ierr)
            call uclgsr8('mgabund',abel2(12),ierr)
            call uclgsr8('alabund',abel2(13),ierr)
            call uclgsr8('siabund',abel2(14),ierr)
            call uclgsr8('pabund',abel2(15),ierr)
            call uclgsr8('sabund',abel2(16),ierr)
            call uclgsr8('clabund',abel2(17),ierr)
            call uclgsr8('arabund',abel2(18),ierr)
            call uclgsr8('kabund',abel2(19),ierr)
            call uclgsr8('caabund',abel2(20),ierr)
            call uclgsr8('scabund',abel2(21),ierr)
            call uclgsr8('tiabund',abel2(22),ierr)
            call uclgsr8('vabund',abel2(23),ierr)
            call uclgsr8('crabund',abel2(24),ierr)
            call uclgsr8('mnabund',abel2(25),ierr)
            call uclgsr8('feabund',abel2(26),ierr)
            call uclgsr8('coabund',abel2(27),ierr)
            call uclgsr8('niabund',abel2(28),ierr)
            call uclgsr8('cuabund',abel2(29),ierr)
            call uclgsr8('znabund',abel2(30),ierr)
            do mm=1,nl
              abel(mm)=abel2(mm)
              enddo
c
c
            call uclgsi('npass',npass,ierr)
C           Test if npass is even.  If it is, change it to next lowest odd
c            if(mod(npass,2) .eq. 0) then
c              write(lun11,*)'rread1: npass should always be odd.'
c              write(lun11,*)'rread1: resetting to ',npass-1
c              npass=npass-1
c            endif

            stringsl=kblnk80
            call uclgst('modelname',stringsl,ierr)
            read (stringsl(1:30),'(a30)')kmodelname
c
c
c           step parameters
            call uclgsr8('emult',emult,ierr)
            call uclgsr8('taumax',taumax,ierr)
c
c           min xee
            call uclgsr8('xeemin',xeemin,ierr)
c
c           critf
            call uclgsr8('critf',critf,ierr)
c
c           vturbi
            call uclgsr8('vturbi',vturbi,ierr)
c
c           ncn2
            call uclgsi('ncn2',ncn2,ierr)
            ncn2=max(999,min(999999,ncn2))
            if (ierr.ne.0) ncn2=9999
c
c           radexp
            call uclgsr8('radexp',radexp,ierr)
            if (ierr.ne.0) radexp=0.
c
            call uclgsi('loopcontrol',nloopctl,ierr)
c
c
            ccc = 3.e+10
            xlum2=xlum
            xpxcol=xpcol
            xi=10.**zeta
            if (lcdd.ne.1) then
               xpx = p/1.38e-12/max(t,1.e-24)
               r19 = sqrt(xlum2/12.56/ccc/max(1.e-24,p*xi))
            else
               r19 = sqrt(xlum2/max(1.e-24,xpx*xi))
            endif
            rmax = xpxcol/(max(xpx,1.e-36))
            r = r19*(1.e+19)
c
      return
      end
      subroutine rstepr(unit,hdunum,radin,radout,rdel,t,prs,
     $             xcol,xee,xpx,xi,
     $             idat1,rdat1,kdat1,nptrs,npnxt,npfi,
     $             npfirst,npar,npilev,
     $             xilev,bilev,rniss,nloopctl,
     $             lun11,status)
C
C     Reads a FITS extension binary table containing
C     nrhs columns and at most nrhdimj rows
c     author: T. Kallman
C
C     Parameters:
C        unit    integer            File unit number
C        hdunum  integer            Number of last HDU written
C        radin   real               inner radius of shell
C        radout  real               outer radius of shell
C                                   nb but now it is delr in the call
C        t    real               temperature of shell
C        prs    real               pressure in shell
C        nrhdimj  integer            Maximum number of rows
C        idat1   integer(nidat1)    Needed by the atomic database
C        rdat1   real(nidat1)       Needed by the atomic database
C        kdat1   char*nidat1        Needed by the atomic database
C        nptrs                      Needed by the atomic database
C        npnxt                      Needed by the atomic database
C        npfi                       Needed by the atomic database
C        npfirst                    Needed by the atomic database
C        npcon                      Needed by the atomic database
C        npconi                     Needed by the atomic database
C        npcon2                     Needed by the atomic database
C        xilev   real(nrhdimj)       Fractional level population array
C        cemab   real(2,nrhdimj)     Recombination emission
C        opakab  real(nrhdimj)       Opacity
C        tauc    real(2,nrhdimj)     Optical depth
C        poptol  real               Tolerance for population level
C        nloopctl integer           Loop control variable
C        nzone   integer            Pass number through iteration process
C        status  integer            Returned status code
C
      implicit none
      include './PARAM'
      integer mllz

      integer nptmpdim
      parameter (nptmpdim=400000)
c
C     Allocation for passed parameters
      real*8 xilev(nnml),bilev(nnml),rniss(nnml)
      real*8 rdat1(nrdat1)
      real rtmp
      real*8 radin, radout, rdel,t, prs, xcol,xee,xpx,xi
      integer unit,hdunum, nrows, status, nloopctl
      integer idat1(nidat1),nptrs(nptt,ndat2)

c     pointers to master data
      integer npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni),npar(ndat2)
      integer npilev(nd,nni)

C     Internal work areas
      integer ntptr
      integer tfields,varidat
      character(16) ttype(5),tform(5),tunit(5)
      integer colnum,frow,felem,hdutype,ll, ltyp
      integer lrtyp, lcon, nrdt, nidt, mm, lun11, lpril,lpri
      integer jkk,nidti,nlines
      real eliml,elimh,elin,elmmtpp,elcomp,eth,xeltp
      character(33) extname
      character(20) kcom
      character(1) kdat1(nkdat1)

C     Database manipulation quantities
      integer nelems,nullj,nkeys,irow2,nspace
      real anynull,nulle
      character(1) kblnk,kdtmp(200),nullstr
      logical done
      integer nhdu

      data kblnk/' '/
c
      data tform/'1J','1E','1E','1E','1E'/

      data ttype/'index','energy','opacity','fwd dpth',
     $ 'bck dpth'/

      data tunit/' ','ev','/cm',' ',' '/

      varidat=0
c
      lpri=0
      lpril=lpri
c
      if (lpri.ne.0)
     $ write (lun11,*)'in rstepr ',hdunum
c

      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu
      call ftmahd(unit,1,hdutype,status)
      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu
c
C     Move to the appropriate HDU (hdunum) in the file
      mm=hdunum
c
      if (lpri.ne.0)
     $ write(lun11,*)'rstepr2: Moving to extension',mm
      call ftmahd(unit,mm,hdutype,status)
      if (lpri.ne.0)
     $ write (lun11,*)unit,mm,hdutype,status
      if (status .gt. 0)call printerror(lun11,status)
      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu

C     Determine the number of keywords in the header
      nkeys=0
      call ftghsp(unit,nkeys,nspace,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftghsp:',unit,nkeys,nspace,status
c
c
C     Read each 80-character keyword record, and print it out
      call ftgkyj(unit,'NAXIS2',nrows,kcom,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkyj:',nrows,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'RINNER',rtmp,kcom,status)
      radin=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',radin,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'ROUTER',rtmp,kcom,status)
      radout=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',radout,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'RDEL',rtmp,kcom,status)
      rdel=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',rdel,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'TEMPERAT',rtmp,kcom,status)
      t=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',t,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'PRESSURE',rtmp,kcom,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',prs,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'COLUMN',rtmp,kcom,status)
      xcol=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye, xcol=',xcol,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'XEE',rtmp,kcom,status)
      xee=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xee,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'DENSITY',rtmp,kcom,status)
      xpx=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xpx,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'LOGXI',rtmp,kcom,status)
      xi=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xi,kcom,status
      if (status .gt. 0)call printerror(lun11,status)


      felem=1
      nelems=1
      nullstr=' '
      nullj=0
      nulle=0.
      do irow2=1,nrows
        if (lpri.ne.0)
     $   write (lun11,*)'row=',irow2
        colnum=1
        call ftgcvj(unit,colnum,irow2,felem,nelems,nullstr,
     $                  ntptr,anynull,status)
        colnum=7
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
c       note that here we switch the inward and outward
        xilev(ntptr)=rtmp
        colnum=8
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        rniss(ntptr)=rtmp
        enddo
C

c----------------------------------------------------------------
C     Compute checksums
      call ftpcks(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

      return
      end
      subroutine rstepr2(unit,hdunum,radin,radout,rdel,t,prs,
     $             xcol,xee,xpx,xi,
     $             idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $             nplin,nlsvn,rcem,oplin,tau0,nloopctl,
     $             lun11,status)
C
C     Reads a FITS extension binary table containing
C     nrhs columns and at most nrhdimj rows
c     author: T. Kallman
C
C     Parameters:
C        unit    integer            File unit number
C        hdunum  integer            Number of last HDU written
C        radin   real               inner radius of shell
C        radout  real               outer radius of shell
C                                   nb but now it is delr in the call
C        t    real               temperature of shell
C        prs    real               pressure in shell
C        nrhdimj  integer            Maximum number of rows
C        idat1   integer(nidat1)    Needed by the atomic database
C        rdat1   real(nidat1)       Needed by the atomic database
C        kdat1   char*nidat1        Needed by the atomic database
C        nptrs                      Needed by the atomic database
C        npnxt                      Needed by the atomic database
C        npfi                       Needed by the atomic database
C        npfirst                    Needed by the atomic database
C        npcon                      Needed by the atomic database
C        npconi                     Needed by the atomic database
C        npcon2                     Needed by the atomic database
C        xilev   real(nrhdimj)       Fractional level population array
C        cemab   real(2,nrhdimj)     Recombination emission
C        opakab  real(nrhdimj)       Opacity
C        tauc    real(2,nrhdimj)     Optical depth
C        poptol  real               Tolerance for population level
C        nloopctl integer           Loop control variable
C        nzone   integer            Pass number through iteration process
C        status  integer            Returned status code
C
      implicit none
      include './PARAM'
      integer mllz

      integer nptmpdim
      parameter (nptmpdim=400000)
c
C     Allocation for passed parameters
      real*8 tau0(2,nnnl), rcem(2,nnnl)
      real*8 rdat1(nrdat1)
      real rtmp
      real*8 radin, radout,rdel, t, prs, xcol,xee,xpx,xi
      integer unit,hdunum, nrows, status, nloopctl
      integer idat1(nidat1),nptrs(nptt,ndat2)
c     line opacities
      real*8 oplin(nnnl)

c     pointers to master data
      integer npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni),npar(ndat2)
      integer nplin(nnnl),nlsvn

C     Internal work areas
      integer ntptr
      integer tfields,varidat
      character(16) ttype(5),tform(5),tunit(5)
      integer colnum,frow,felem,hdutype,ll, ltyp
      integer lrtyp, lcon, nrdt, nidt, mm, lun11, lpril,lpri
      integer jkk,nidti,nlines
      real eliml,elimh,elin,elmmtpp,elcomp,eth,xeltp
      character(33) extname
      character(20) kcom
      character(1) kdat1(nkdat1)
      integer nhdu

C     Database manipulation quantities
      integer nelems,nullj,nkeys,irow2,nspace
      real anynull,nulle
      character(1) kblnk,kdtmp(200),nullstr
      logical done

      data kblnk/' '/
c
      data tform/'1J','1E','1E','1E','1E'/

      data ttype/'index','energy','opacity','fwd dpth',
     $ 'bck dpth'/

      data tunit/' ','ev','/cm',' ',' '/

      varidat=0
c
      lpri=0
      lpril=lpri
c
      if (lpri.ne.0)
     $ write (lun11,*)'in rstepr2 ',hdunum

      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu
      call ftmahd(unit,1,hdutype,status)
      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu
c
C     Move to the appropriate HDU (hdunum) in the file
      mm=hdunum
c
      if (lpri.ne.0)
     $ write(lun11,*)'rstepr2: Moving to extension',mm
      call ftmahd(unit,mm,hdutype,status)
      if (lpri.ne.0)
     $ write (lun11,*)unit,mm,hdutype,status
      if (status .gt. 0)call printerror(lun11,status)
      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu

C     Determine the number of keywords in the header
      nkeys=0
      call ftghsp(unit,nkeys,nspace,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftghsp:',unit,nkeys,nspace,status
c
c
c
C     Read each 80-character keyword record, and print it out
      call ftgkyj(unit,'NAXIS2',nrows,kcom,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkyj:',nrows,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'RINNER',rtmp,kcom,status)
      radin=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',radin,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'ROUTER',rtmp,kcom,status)
      radout=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',radout,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'RDEL',rtmp,kcom,status)
      rdel=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',rdel,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'TEMPERAT',rtmp,kcom,status)
      t=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',t,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'PRESSURE',rtmp,kcom,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',prs,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'COLUMN',rtmp,kcom,status)
      xcol=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye, xcol=',xcol,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'XEE',rtmp,kcom,status)
      xee=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xee,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'DENSITY',rtmp,kcom,status)
      xpx=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xpx,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'LOGXI',rtmp,kcom,status)
      xi=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xi,kcom,status
      if (status .gt. 0)call printerror(lun11,status)


      felem=1
      nelems=1
      nullstr=' '
      nullj=0
      nulle=0.
      do irow2=1,nrows
      if (lpri.ne.0)
     $   write (lun11,*)'row=',irow2
        colnum=1
        call ftgcvj(unit,colnum,irow2,felem,nelems,nullstr,
     $                  ntptr,anynull,status)
        colnum=6
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        rcem(1,ntptr)=rtmp
        colnum=7
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        rcem(2,ntptr)=rtmp
        colnum=8
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        oplin(ntptr)=rtmp
        colnum=9
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        tau0(1,ntptr)=rtmp
        colnum=10
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        tau0(2,ntptr)=rtmp
        enddo
C

c----------------------------------------------------------------
C     Compute checksums
      call ftpcks(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

      return
      end
      subroutine rstepr3(unit,hdunum,radin,radout,rdel,t,prs,
     $             xcol,xee,xpx,xi,
     $             idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $             npfirst,npilev,npconi2,ncsvn,
     $             rniss,cemab,cabab,opakab,tauc,nloopctl,
     $             lun11,status)
C
C     Reads a FITS extension binary table containing
C     nrhs columns and at most nrhdimj rows
c     author: T. Kallman
C
C     Parameters:
C        unit    integer            File unit number
C        hdunum  integer            Number of last HDU written
C        radin   real*8               inner radius of shell
C        radout  real*8               outer radius of shell
C                                   nb but now it is delr in the call
C        t    real*8               temperature of shell
C        prs    real*8               pressure in shell
C        nrhdimj  integer            Maximum number of rows
C        idat1   integer(nidat1)    Needed by the atomic database
C        rdat1   real(nidat1)       Needed by the atomic database
C        kdat1   char*nidat1        Needed by the atomic database
C        nptrs                      Needed by the atomic database
C        npnxt                      Needed by the atomic database
C        npfi                       Needed by the atomic database
C        npfirst                    Needed by the atomic database
C        npcon                      Needed by the atomic database
C        npconi                     Needed by the atomic database
C        npcon2                     Needed by the atomic database
C        xilev   real(nrhdimj)       Fractional level population array
C        cemab   real(2,nrhdimj)     Recombination emission
C        opakab  real(nrhdimj)       Opacity
C        tauc    real(2,nrhdimj)     Optical depth
C        poptol  real*8               Tolerance for population level
C        nloopctl integer           Loop control variable
C        nzone   integer            Pass number through iteration process
C        status  integer            Returned status code
C
      implicit none
      include './PARAM'
      integer mllz

      integer nptmpdim
      parameter (nptmpdim=400000)
c
C     Allocation for passed parameters
      real*8 rdat1(nrdat1)
      real rtmp
      real*8 radin, radout, rdel,t, prs, xcol,xee,xpx,xi
      integer unit,hdunum, nrows, status, nloopctl
      integer idat1(nidat1),nptrs(nptt,ndat2)
c     line opacities
      real*8 rniss(nnml)
      real*8 tauc(2,nnml)
      real*8 cemab(2,nnml),opakab(nnml),cabab(nnml)

c     pointers to master data
      integer npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni),npar(ndat2)
      integer npilev(nd,nni)
      integer nplin(nnnl)
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer ncsvn

C     Internal work areas
      integer ntptr
      real rwrk1(nptmpdim)
      integer tfields,varidat
      character(16) ttype(5),tform(5),tunit(5)
      integer colnum,frow,felem,hdutype,ll, ltyp
      integer lrtyp, lcon, nrdt, nidt, mm, lun11, lpril,lpri
      integer jkk,nidti,nlines
      real*8 eliml,elimh,elin,elmmtpp,elcomp,eth,xeltp
      character(33) extname
      character(20) kcom
      character(1) kdat1(nkdat1)

C     Database manipulation quantities
      integer nelems,nullj,nkeys,irow2,nspace,nhdu
      real anynull,nulle
      character(1) kblnk,kdtmp(200),nullstr
      logical done

      data kblnk/' '/
c
      data tform/'1J','1E','1E','1E','1E'/

      data ttype/'index','energy','opacity','fwd dpth',
     $ 'bck dpth'/

      data tunit/' ','ev','/cm',' ',' '/

      varidat=0
c
      lpri=0
      lpril=lpri
c
      if (lpri.ne.0)
     $ write (lun11,*)'in rstepr3 ',hdunum
c
      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu
      call ftmahd(unit,1,hdutype,status)
      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu
c
C     Move to the appropriate HDU (hdunum) in the file
      mm=hdunum
c
      if (lpri.ne.0)
     $ write(lun11,*)'rstepr3: Moving to extension',mm
      call ftmahd(unit,mm,hdutype,status)
      if (lpri.ne.0)
     $ write (lun11,*)unit,mm,hdutype,status
      if (status .gt. 0)call printerror(lun11,status)
      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu

C     Determine the number of keywords in the header
      nkeys=0
      call ftghsp(unit,nkeys,nspace,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftghsp:',unit,nkeys,nspace,status
c
c
c
C     Read each 80-character keyword record, and print it out
      call ftgkyj(unit,'NAXIS2',nrows,kcom,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkyj:',nrows,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'RINNER',rtmp,kcom,status)
      radin=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',radin,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'ROUTER',rtmp,kcom,status)
      radout=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',radout,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'RDEL',rtmp,kcom,status)
      rdel=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',rdel,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'TEMPERAT',rtmp,kcom,status)
      t=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',t,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'PRESSURE',rtmp,kcom,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',prs,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'COLUMN',rtmp,kcom,status)
      xcol=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye, xcol=',xcol,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'XEE',rtmp,kcom,status)
      xee=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xee,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'DENSITY',rtmp,kcom,status)
      xpx=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xpx,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'LOGXI',rtmp,kcom,status)
      xi=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xi,kcom,status
      if (status .gt. 0)call printerror(lun11,status)


      felem=1
      nelems=1
      nullstr=' '
      nullj=0
      nulle=0.
      do irow2=1,nrows
      if (lpri.ne.0)
     $   write (lun11,*)'row=',irow2
        colnum=1
        call ftgcvj(unit,colnum,irow2,felem,nelems,nullstr,
     $                  ntptr,anynull,status)
        colnum=7
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        cemab(1,ntptr)=rtmp
        colnum=8
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        cemab(2,ntptr)=rtmp
        colnum=9
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        cabab(ntptr)=rtmp
        colnum=10
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        opakab(ntptr)=rtmp
        colnum=11
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        tauc(1,ntptr)=rtmp
        colnum=12
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        tauc(2,ntptr)=rtmp
        enddo
C

c----------------------------------------------------------------
C     Compute checksums
      call ftpcks(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

      return
      end
      subroutine rstepr4(unit,hdunum,radin,radout,rdel,t,prs,
     $             xcol,xee,xpx,xi,
     $             idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $             epi,ncn2,dpthc,opakc,rccemis,nloopctl,
     $             lun11,status)
C
C     Reads a FITS extension binary table containing
C     nrhs columns and at most nrhdimj rows
c     author: T. Kallman
C
C     Parameters:
C        unit    integer            File unit number
C        hdunum  integer            Number of last HDU written
C        radin   real*8               inner radius of shell
C        radout  real*8               outer radius of shell
C                                   nb but now it is delr in the call
C        t    real*8               temperature of shell
C        prs    real*8               pressure in shell
C        nrhdimj  integer            Maximum number of rows
C        idat1   integer(nidat1)    Needed by the atomic database
C        rdat1   real(nidat1)       Needed by the atomic database
C        kdat1   char*nidat1        Needed by the atomic database
C        nptrs                      Needed by the atomic database
C        npnxt                      Needed by the atomic database
C        npfi                       Needed by the atomic database
C        npfirst                    Needed by the atomic database
C        npcon                      Needed by the atomic database
C        npconi                     Needed by the atomic database
C        npcon2                     Needed by the atomic database
C        xilev   real(nrhdimj)       Fractional level population array
C        cemab   real(2,nrhdimj)     Recombination emission
C        opakab  real(nrhdimj)       Opacity
C        tauc    real(2,nrhdimj)     Optical depth
C        poptol  real*8               Tolerance for population level
C        nloopctl integer           Loop control variable
C        nzone   integer            Pass number through iteration process
C        status  integer            Returned status code
C
      implicit none
      include './PARAM'
      integer mllz

      integer nptmpdim
      parameter (nptmpdim=400000)
c
C     Allocation for passed parameters
      real*8 rdat1(nrdat1)
      real rtmp
      real*8 radin, radout,rdel, t, prs, xcol,xee,xpx,xi
      integer unit,hdunum, nrows, status, nloopctl
      integer idat1(nidat1),nptrs(nptt,ndat2)
c     energy bins
      real*8 epi(ncn)
c     continuum opacities
      real*8 opakc(ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn)
      integer ncn2

c     pointers to master data
      integer npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni),npar(ndat2)

C     Internal work areas
      real ntptr(nptmpdim)
      real rwrk1(nptmpdim)
      integer tfields,varidat
      character(16) ttype(5),tform(5),tunit(5)
      integer colnum,frow,felem,hdutype,ll, ltyp
      integer lrtyp, lcon, nrdt, nidt, mm, lun11, lpril,lpri
      integer jkk,nidti,nlines
      real*8 eliml,elimh,elin,elmmtpp,elcomp,eth,xeltp
      character(33) extname
      character(20) kcom
      character(1) kdat1(nkdat1)

C     Database manipulation quantities
      integer nelems,nullj,nkeys,irow2,nspace,nhdu
      real anynull,nulle
      character(1) kblnk,kdtmp(200),nullstr
      logical done

      data kblnk/' '/
c
      data tform/'1J','1E','1E','1E','1E'/

      data ttype/'index','energy','opacity','fwd dpth',
     $ 'bck dpth'/

      data tunit/' ','ev','/cm',' ',' '/

      varidat=0
c
      lpri=0
      lpril=lpri
c
      if (lpri.ne.0)
     $ write (lun11,*)'in rstepr4 ',hdunum
c
      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu
      call ftmahd(unit,1,hdutype,status)
      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu
c
C     Move to the appropriate HDU (hdunum) in the file
      mm=hdunum
c
      if (lpri.ne.0)
     $ write(lun11,*)'rstepr4: Moving to extension',mm
      call ftmahd(unit,mm,hdutype,status)
      if (lpri.ne.0)
     $ write (lun11,*)unit,mm,hdutype,status
      if (status .gt. 0)call printerror(lun11,status)
      call FTGHDN(unit, nhdu)
      if (lpri.ne.0)
     $ write (lun11,*)'current hdu ',nhdu

C     Determine the number of keywords in the header
      nkeys=0
      call ftghsp(unit,nkeys,nspace,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftghsp:',unit,nkeys,nspace,status
c
c
c
C     Read each 80-character keyword record, and print it out
      call ftgkyj(unit,'NAXIS2',nrows,kcom,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkyj:',nrows,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'RINNER',rtmp,kcom,status)
      radin=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',radin,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'ROUTER',rtmp,kcom,status)
      radout=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',radout,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'RDEL',rtmp,kcom,status)
      rdel=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',rdel,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'TEMPERAT',rtmp,kcom,status)
      t=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',t,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'PRESSURE',rtmp,kcom,status)
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',prs,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'COLUMN',rtmp,kcom,status)
      xcol=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye, xcol=',xcol,kcom,status
      if (status .gt. 0)call printerror(lun11,status)

      call ftgkye(unit,'XEE',rtmp,kcom,status)
      xee=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xee,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'DENSITY',rtmp,kcom,status)
      xpx=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xpx,kcom,status
      if (status .gt. 0)call printerror(lun11,status)
c
      call ftgkye(unit,'LOGXI',rtmp,kcom,status)
      xi=rtmp
      if (lpri.ne.0)
     $ write (lun11,*)'after ftgkye',xi,kcom,status
      if (status .gt. 0)call printerror(lun11,status)


      felem=1
      nelems=1
      nullstr=' '
      nullj=0
      nulle=0.
      do irow2=1,nrows
      if (lpri.ne.0)
     $   write (lun11,*)'row=',irow2
        colnum=3
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        opakc(irow2)=rtmp
        colnum=4
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        dpthc(1,irow2)=rtmp
        colnum=5
        call ftgcve(unit,colnum,irow2,felem,nelems,nullstr,
     $                  rtmp,anynull,status)
        dpthc(2,irow2)=rtmp
        enddo
C

c----------------------------------------------------------------
C     Compute checksums
      call ftpcks(unit,status)
      if (status .gt. 0)call printerror(lun11,status)

      return
      end
      subroutine savd(jkstep,ldir,
     $       lpri,iunit,iunit2,iunit3,iunit4,
     $       idat1,rdat1,kdat1,nptrs,npnxt,npfi,
     $       npfirst,npar,npilev,npconi2,ncsvn,
     $       t,p,r,rdel,delr,xcol,xee,xpx,zeta,
     $       xilev,bilev,rniss,abel,
     $       nplin,nlsvn,rcem,oplin,tau0,
     $       cemab,cabab,opakab,tauc,
     $       epi,ncn2,dpthc,opakc,rccemis,nloopctl,
     $       lunlog,status)
c
c     this routine  saves only depths for iterative calculation
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
C     Allocation for passed parameters
      real*8 rdat1(nrdat1)
      real*8 r,delr,rdel, t, p, xcol,xee,xpx,zeta
      integer unit,hdunum, nrows, status, nloopctl
      integer iunit,iunit2,iunit3,iunit4
      integer idat1(nidat1),nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
      integer nlsvn,ncsvn
c     energy bins
      real*8 epi(ncn)
c     continuum opacities
      real*8 opakc(ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn)
      integer ncn2
c     line opacities
      real*8 oplin(nnnl)
      real*8 abel(nl)
      real*8 tauc(2,nnml)
      real*8 cemab(2,nnml),opakab(nnml),cabab(nnml)
c     pointers to master data
      integer npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni),npar(ndat2)
      integer npilev(nd,nni)
      integer nplin(nnnl)
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      real*8 xilev(nnml),bilev(nnml),rniss(nnml)
      real*8 tau0(2,nnnl), rcem(2,nnnl)

c     continuum optical depths
      integer ldir,lpri,lun,lunlog,jkstep
      integer lind1,lind2,kl,ll,nlyc,nry,nbinc
c
c      write (lunlog,*)'in savd',jkstep,iunit,iunit2,iunit3,iunit4
      call fstepr(iunit,jkstep,r,delr,rdel,t,p,
     $          xcol,xee,xpx,zeta,
     $          idat1,rdat1,kdat1,nptrs,npnxt,npfi,
     $          npfirst,npar,npilev,
     $          xilev,bilev,rniss,nloopctl,
     $          lunlog,status)
      call fstepr2(iunit2,jkstep,r,delr,rdel,t,p,
     $          xcol,xee,xpx,zeta,
     $          idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $          nplin,nlsvn,rcem,oplin,tau0,nloopctl,
     $          lunlog,status)
      call fstepr3(iunit3,jkstep,r,delr,rdel,t,p,abel,
     $          xcol,xee,xpx,zeta,
     $          idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $          npfirst,npilev,npconi2,ncsvn,
     $          rniss,cemab,cabab,opakab,tauc,nloopctl,
     $          lunlog,status)
      call fstepr4(iunit4,jkstep,r,delr,rdel,t,p,
     $          xcol,xee,xpx,zeta,
     $          idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $          epi,ncn2,dpthc,opakc,rccemis,nloopctl,
     $          lunlog,status)
c      write (lunlog,*)' r=',r,xpx,t,delr,xcol,xee,xpx,zeta
      if (status .gt. 0)call printerror(lunlog,status)
c      write(lun,REC=jkstep)t,r,rdel,delr,xcol,xee,xpx,tau0d,dpthcd,taucd
      nlyc=nbinc(13.7d0,epi,ncn2)
      nry=nlyc+1
      if (lpri.ne.0) write (lunlog,*)'in savd',rdel,t,tauc(1,25),
     $                ldir,dpthc(1,nry),dpthc(2,nry)
c
      return
      end
      subroutine setptrs(lun11,lpri,
     $ idat1,rdat1,kdat1,nptrs,np2,
     $ npnxt,npfi,npar,npfirst,nplin,
     $ nplini,npcon,npconi,npilev,npilevi,
     $ npconi2,nlevs,nlsvn,ncsvn,abcosmic,abel)
c
c     this program set the pointers of the database
c       Written by Ke Zhang, Oct.8 2001
c
c     data structures are:
c      data: the database arrays (integer, real, character)
c       idat1(nidat1)
c       rdat1(nrdat1),
c       kdat1(nkdat1)
c     descriptions of database entries, and pointers
c       nptrs(nptt,ndat2)
c         nptrs(2,nx)=data type
c         nptrs(3,nx)=rate type
c         nptrs(4,nx)=continuation flag
c                       (n=number of continuations to come)
c         nptrs(5,nx)=number of reals
c         nptrs(6,nx)=number of integers
c         nptrs(7,nx)=number of characters
c         nptrs(8,nx)=pointer to reals
c         nptrs(9,nx)=pointer to integers
c         nptrs(10,nx)=pointer to characters
c
c       pointers:
c       next record:
c         npnxt(ndat2)
c       parent record (=ion header or element header)
c         npar(ndat2)
c       first record of a given rate type
c         npfirst(ntyp)
c       first record of rate type ntyp for ion nni
c         npfi(ntyp,nni)
c       pointer for line data from array containing luminosities
c         nplin(nnnl)
c       (inverse) pointer for line data to array containing luminosities
c          from database array
c         nplini(ndat2)
c       pointer for continuum data (pi xsection) from array containing luminosities
c         npcon(nnml)
c       pointer to abundance array to first level of element nni
c         npconi2(ndat2)
c       (inverse) pointer for continuum data (pi xsection) from array containing
c           luminosities
c         npconi(ndat2)
c
c
      implicit none
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
      real*8 rdat1(nrdat1)
      real*8 abel(nl),abcosmic(30)
c
      integer idat1(nidat1)
      integer nptrs(nptt,ndat2)
      integer npnxt(ndat2),npar(ndat2)
      integer npfirst(ntyp)
      integer npnxt2(ndat2)
      integer npfi(ntyp,nni)
      integer nplin(nnnl),nplini(ndat2),npcon(nnml)
      integer npilev(nd,nni),npilevi(nnml)
      integer npconi2(ndat2)
      integer npconi(ndat2)
      integer nlevs(nni)
      integer melpt(nl)
      integer mlold(ntyp)
      integer indx, iion, ilev, icon, iline, i, j
      integer iel2, iel, lpri, lun11, np2, lrtp
      integer iilev, mltmpn, mlfnd, nclev, mltst
      integer mltmp, npartmpn, nlsvn, ncsvn
      integer lsrt, mml, niter, melptmp, npfirst2
      integer mllo, mlloo, itst, ltyp, lrtyp2, lcon
      integer nrdt, nidt, nkdt, mll,mlm
      integer itmp,mm,np1i,np1r,np1k,ntptmp
c
      character(1) kdat1(nkdat1)

c            pointer structure
c     type    desc         nr  ni  nk      daught  par
c     1       rr, a&p      2   1   0               14
c     2       hcx          4   1   0               14
c     3       ai           2   1   0               14
c     4       line dat 1   2   3   0        5      14
c     5       line dat 2   4   3   0                4
c     6       lev dat  1   4   3   0               14
c     7       dr a&p       5   1   0               14
c     8       dr a&r       0   0   0               14
c     9       hecx         4   1   0               14
c     10      lev dat 2    0   2  30                6
c     11      2 ph         2   2   0               14
c     12      pixc, bpl    5   2   0               14
c     13      el           2   2  30       14       0
c     14      ion          1   2   8       all     13
c     15      pixc bkh 1   5   1   0       20      14
c     16      pixc bkh     0   0   0               14
c     17      cx: cota     4   3   0               14
c     18      rr: cota     3   1   0               14
c     19      pixc hullac  0   0   0               14
c     20      pixc bkh 2   5   1   0       21      15
c     21      pixc bkh 3   4   4  11               20
c     22      dr stroey    5   1   0               14
c     23      pixc clark   5   2   0       24      14
c     24      pixc clark 2 4   4   0               23
c     25      ci r&s       0   0   0               14
c     26      ci cota      2   2   0               14
c
        indx=1
c
c the main data index
      go to 9009
      if (lpri.ne.0) then
c       first an experimental print
        write (lun11,*)'np2=',np2
        do itmp=1,np2
          CALL DRD(ltyp,lrtyp2,lcon,nrdt,np1r,nidt,np1i,nkdt,np1k,
     &          itmp-1,Nptrs,0,Lun11)
          write (lun11,*)'itmp=',itmp
c          write (lun11,*)'nkdt=',nkdt,(kdat1(np1k-1+mm),mm=1,nkdt)
          call dprinto(ltyp,lrtyp2,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,rdat1,idat1,kdat1,lun11)
          enddo
        endif
 9009   continue

c the ion index
      iion=1

c the level index
      ilev=1

c the continum index
      icon=1

c the line index
      iline=1

c initialize pointers

      do i=1,ntyp
        npfirst(i)=0
      enddo

      do i=1,nni
        do j=1,ntyp
          npfi(j,i)=0
        enddo
        nlevs(i)=0
      enddo

      do i=1,ndat2
        npnxt(i)=0
      enddo
c
      do i=1,nl
        abcosmic(i)=0.
        enddo

      mlold(11)=0
      iel2=0
      do while ((iel2.le.nl).and.(indx.lt.np2))
        iel2=iel2+1
        iel=iel2
c        if (abel(iel).lt.1.e-15) then
c
          if (lpri.ne.0)
     $     write (lun11,*)'iel=',iel2,iel,abel(iel)
c  pass by elements that has neglectable abundance
c          indx=indx+1
c          do while((nptrs(3,indx).ne.11).and.(indx.lt.np2))
c            indx=indx+1
c          enddo
c
c        else

c  register element record
c  npfirst,npnxt,npar,mlold

          if (npfirst(11).eq.0) then
            npfirst(11)=indx
          else
            npnxt(mlold(11))=indx
          endif
c
          CALL DRD(ltyp,lrtyp2,lcon,nrdt,np1r,nidt,np1i,nkdt,np1k,
     &          indx-1,Nptrs,0,Lun11)
          iel=idat1(np1i)
          abcosmic(iel)=rdat1(np1r)

          if (lpri.ne.0)
     $     write (lun11,*)'registering element:',iel,abel(iel),indx,
     $                         mlold(11),abcosmic(iel)
          mlold(11)=indx
          npar(indx)=0
          indx=indx+1

c  go through ions

          ltyp=nptrs(2,indx)
          lrtp=nptrs(3,indx)
          if (lpri.ne.0)
     $     write (lun11,*)'lrtp=',lrtp,indx
          do while(lrtp.eq.12)

c
          if (lpri.ne.0)
     $     write(lun11,*) iel,iion

c  register ion record
c  npfirst,npnxt,npar,mlold

            if (npfirst(12).eq.0) then
              npfirst(12)=indx
            else
              npnxt(mlold(12))=indx
            endif
            if (lpri.ne.0)
     $       write (lun11,*)'npfirst(12)=',npfirst(12),indx
            mlold(12)=indx
            npar(indx)=mlold(11)
            indx=indx+1

c  level records, rate type 13
c  npfirst,npnxt,npar,mlold,npfi,npilev,npilevi

            if (nptrs(3,indx).eq.13) then
              npfi(13,iion)=indx
              iilev=1
              if (npfirst(13).eq.0) then
                npfirst(13)=indx
              else
                npnxt(mlold(13))=indx
              endif
              if (lpri.ne.0)
     $         write (lun11,*)'filling npilev:'
              do while(nptrs(3,indx).eq.13)
                npar(indx)=mlold(12)
                npnxt(indx)=indx+1
                npilev(iilev,iion)=ilev
                npilevi(ilev)=iilev
            if(npilev(iilev,iion).eq.0)print *, 'AJA **** ', iilev,iion
                if (lpri.ne.0)
     $           write (lun11,*)ilev,iilev,indx,iion
                iilev=iilev+1
                ilev=ilev+1
                indx=indx+1
              enddo
              mlold(13)=indx-1
              npnxt(indx-1)=0
            endif


            do i=1,2
              if (i.eq.1) then
                lrtp=7
              else
                lrtp=1
              endif
              if (nptrs(3,indx).eq.lrtp) then
                npfi(lrtp,iion)=indx
                if (npfirst(lrtp).eq.0) then
                  npfirst(lrtp)=indx
                else
                  npnxt(mlold(lrtp))=indx
                endif
                if (lpri.ne.0)
     $           write (lun11,*)'npconi loop',indx
                do while(nptrs(3,indx).eq.lrtp)
c                 npcon points from the array of continuum emissivities
c                    to the photoionization data
c                 npconi points from the levels to the arrays of
c                    array of continuum emissivities
c                 npconi2 points from the photoionization data
c                    to the array of continuum emissivities
c                    (inverse if npcon)
c                 icon is the index of the continuum emissivity array
c                    element
c                 indx is the index of the photoionization data
                  npar(indx)=mlold(12)
                  npnxt(indx)=indx+1
                  npcon(icon)=indx
                  if (lpri.ne.0)
     $             write (lun11,*)'index into continuum  array:',
     $                icon
                  if (lpri.ne.0)
     $             write (lun11,*)'index of photoionization element:',
     $                indx
c                 now search for the level that goes with this
c                    photoionization data
                  mltmpn=npfi(13,iion)
                  mlfnd=0
                  nclev=idat1(nptrs(6,indx)+nptrs(9,indx)-2)
                  if (nclev.gt.nlevs(iion)) nlevs(iion)=nclev
                  mltst=nclev
                  if (lpri.ne.0)
     $             write (lun11,*)'searching for level:'
                  mltmp=mltmpn
                  if (mltmpn.ne.0) then
                      npartmpn=npar(mltmpn)
                    else
                      npartmpn=0
                    endif
                  do while ((mlfnd.ne.mltst).and.(mltmpn.ne.0)
     $             .and.(indx.ne.0).and.(npartmpn.eq.npar(indx)))
                    mltmp=mltmpn
                    mlfnd=idat1(nptrs(6,mltmp)+nptrs(9,mltmp)-2)
                    mltmpn=npnxt(mltmp)
                    if (mltmpn.ne.0) then
                        npartmpn=npar(mltmpn)
                      else
                        npartmpn=0
                      endif
                    if (lpri.ne.0)
     $              write (lun11,*)mltmp,mlfnd,mltmpn,npartmpn,
     $                               npar(indx),nclev
                    enddo
                  npconi2(indx)=icon
c                  npconi(icon)=npfi(13,iion)-1+nclev
                  if (mltmp.ne.0) then
                    npconi(mltmp)=icon
                    endif
                  if (lpri.ne.0)
     $             write (lun11,*)indx,npar(indx),icon,nclev,
     $                nptrs(3,indx),lrtp,mltmp
                  indx=indx+1
                  icon=icon+1
                enddo
                mlold(lrtp)=indx-1
                npnxt(indx-1)=0
              endif
           enddo
c
c  lines data and lines pointers, rate type 4, 9 & 14
c  npfirst,npnxt,npar,mold,npfi,nplin,nplini

            if (lpri.ne.0)
     $       write (lun11,*)'nplin,nplini,:'
            do i=1,3
              if (i.eq.1) then
                lrtp=4
              elseif (i.eq.2) then
                lrtp=9
              else
                lrtp=14
              endif
              if (lpri.ne.0) write (lun11,*)' indx=',indx,lrtp,iion
              if (nptrs(3,indx).eq.lrtp) then
                npfi(lrtp,iion)=indx
                if (npfirst(lrtp).eq.0) then
                  npfirst(lrtp)=indx
                else
                  npnxt(mlold(lrtp))=indx
                endif
                do while(nptrs(3,indx).eq.lrtp)
                  npar(indx)=mlold(12)
                  npnxt(indx)=indx+1
                  nplin(iline)=indx
                  nplini(indx)=iline
                  if (lpri.ne.0)
     $             write (lun11,*)indx,iline
                  indx=indx+1
                  iline=iline+1
                enddo
                mlold(lrtp)=indx-1
                npnxt(indx-1)=0
              endif
            enddo

c  pointers for rate types 6,8,3,5,40
c  npfirst,npnxt,npar,mold,npfi

            do i=1,5
              if (i.eq.1) then
                lrtp=6
              elseif (i.eq.2) then
                lrtp=8
              elseif (i.eq.3) then
                lrtp=3
              elseif (i.eq.4) then
                lrtp=5
              else
                lrtp=40
              endif
              if (nptrs(3,indx).eq.lrtp) then
                npfi(lrtp,iion)=indx
                if (npfirst(lrtp).eq.0) then
                  npfirst(lrtp)=indx
                else
                  npnxt(mlold(lrtp))=indx
                endif
                do while(nptrs(3,indx).eq.lrtp)
                  npar(indx)=mlold(12)
                  npnxt(indx)=indx+1
                  indx=indx+1
                enddo
                mlold(lrtp)=indx-1
                npnxt(indx-1)=0
              endif
            enddo

c  pointers for other rate types
c  npfirst,npnxt,npar,mold,npfi

            lrtp=nptrs(3,indx)
            do while((lrtp.ne.12).and.(lrtp.ne.11).and.(lrtp.ne.0))
              npar(indx)=mlold(12)
              if (npfirst(lrtp).eq.0) then
                npfirst(lrtp)=indx
              else
                npnxt(mlold(lrtp))=indx
              endif
              mlold(lrtp)=indx
              if (npfi(lrtp,iion).eq.0) npfi(lrtp,iion)=indx
c              write (lun11,*)iion,lrtp,indx,npfi(lrtp,iion)
              indx=indx+1
              lrtp=nptrs(3,indx)
            enddo

c  ionization data and continum pointers, rate type 7 & 1
c  npfirst,npnxt,npar,mlold,npfi,npcon,npconi,npconi2,nlevs


            iion=iion+1

          enddo
c        endif
      enddo

      nlsvn=iline-1
      ncsvn=icon-1
c     write (lun11,*)'nlsvn=',nlsvn,ncsvn
c
      go to 9000
c
c     sort the element abundances
      lsrt=0
      do mml=1,nl
        melpt(mml)=mml
      enddo
      niter=0
      do while (lsrt.eq.0)
        lsrt=1
        niter=niter+1
        do mml=1,nl-1
          if (abel(melpt(mml)).lt.abel(melpt(mml+1))) then
            melptmp=melpt(mml)
            melpt(mml)=melpt(mml+1)
            melpt(mml+1)=melptmp
            lsrt=0
          endif
        enddo
      enddo
c
c
c     now redo the element pointers
c     zero the new next pointers
      do mml=1,np2
        npnxt2(mml)=0
        enddo
      npfirst2=0
      mllo=0
c     step thru elements
      do mml=1,nl
        mlloo=mllo
        mll=npfirst(11)
        itst=0
        do while ((mll.ne.0).and.(itst.ne.melpt(mml)))
          mlm=mll-1
          call drd(ltyp,lrtyp2,lcon,nrdt,np1r,nidt,np1i,nkdt,np1k,
     &      mlm,nptrs,0,lun11)
          itst=idat1(np1i-1+nidt)
          mllo=mll
          mll=npnxt(mll)
          enddo
        if (mllo.ne.0) then
          if (npfirst2.eq.0) then
              npfirst2=mllo
            else
              npnxt2(mlloo)=mllo
            endif
          endif
        enddo
      npnxt2(mlloo)=0
      npfirst(11)=npfirst2
      do mml=1,np2
        if ((npnxt2(mml).ne.0).or.(mml.eq.mlloo)) then
          npnxt(mml)=npnxt2(mml)
          endif
        enddo
c
c
 9000  continue
c
       return
c
c     now print stuff sorted
      ntptmp=11
        mll=npfirst(ntptmp)
        write (lun11,*)'ntptmp=',ntptmp
        do while (mll.ne.0)
          CALL DRD(ltyp,lrtyp2,lcon,nrdt,np1r,nidt,np1i,nkdt,np1k,
     &          mll-1,Nptrs,0,Lun11)
          write (lun11,*)'mll=',mll
          call dprinto(ltyp,lrtyp2,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,rdat1,idat1,kdat1,lun11)
          mll=npnxt(mll)
          enddo
c

c
      return
      end
      function splinem(p1,p2,p3,p4,p5,x)
      implicit none
c
      real*8 p1, p2, p3, p4, p5, x
      real*8 s, s2, s3, s4, x0
      real*8 t0, t1, t2, t3, splinem

c
c 5-point spline interpolation of y(x), for x in the range (0,1)
c knot values p1=y(0), p2=y(1/4), p3=y(1/2), p4=y(3/4), p5=y(1)
c     author:  M. Bautista
c
       s=1./30.
       s2=32.*s*(19.*p1-43.*p2+30.*p3-7.*p4+p5)
       s3=160.*s*(-p1+7.*p2-12.*p3+7.*p4-p5)
       s4=32.*s*(p1-7.*p2+30.*p3-43.*p4+19.*p5)
       if (x.gt.0.25) goto 1
       x0=x-0.125
       t3=0.0
       t2=0.5*s2
       t1=4.*(p2-p1)
       t0=0.5*(p1+p2)-0.015625*t2
       goto 4
 1     if (x.gt.0.5) goto 2
       x0=x-0.375
       t3=20.*s*(s3-s2)
       t2=0.25*(s2+s3)
       t1=4.*(p3-p2)-0.015625*t3
       t0=0.5*(p2+p3)-0.015625*t2
       goto 4
 2     if (x.gt.0.75) goto 3
       x0=x-0.625
       t3=20.*s*(s4-s3)
c      nb this was an error
c       t2=0.25*(s3-s4)
       t2=0.25*(s3+s4)
       t1=4.*(p4-p3)-0.015625*t3
       t0=0.5*(p3+p4)-0.015625*t2
       goto 4
 3     x0=x-0.875
       t3=0.0
       t2=0.5*s4
       t1=4.*(p5-p4)
       t0=0.5*(p4+p5)-0.015625*t2
 4     splinem=t0+x0*(t1+x0*(t2+x0*t3))
       return
       end
      subroutine starf(tp,xlum,epi,ncn2,zremsz,lpri,lun11)
c
c     this subroutine generates the initial spectrum.
c      optically thin bremsstrahlung spectrum
c     brems stores the flux to be used
c     author:  T. Kallman (from xstar1)
c
      implicit none
c
      include './PARAM'
c
      real*8 epi(ncn),zremsz(ncn)
      real*8 zremsi(ncn)
      real*8 ergsev,del,q,xkt,sum,tempp,sum2,xlum,tp,const
      integer numcon,ncn2,lpri,lprisv,i,lun11
c
      data ergsev/1.602197e-12/
c
      numcon=ncn2
      del=1.
      q=7.49e+08*del*xlum/(1.e-37+tp)
      xkt=1.16e-03/(1.e-37+tp)
      sum=0.
      lprisv=lpri
c      lpri=2
      if (lpri.gt.1) write (lun11,*)'in starf',tp,xlum,q,xkt
      do i=1,numcon
         tempp=epi(i)*xkt
         zremsi(i)=0.
c         zremsi(i)=(3.1415e+22)*epi(i)**3/exp(tempp)
         if (tempp.lt.1.e-3) then
             zremsi(i)=epi(i)**3/tempp
           else
             if (tempp.lt.150.)
     $         zremsi(i)=epi(i)**3/(exp(tempp)-1.)
             if (tempp.gt.150.) zremsi(i)=xlum/epi(i)/ergsev/1.e+37
           endif
         zremsi(i)=(3.1415e+22)*zremsi(i)
         if (lpri.gt.1) write (lun11,*)i,epi(i),zremsi(i)
         if (.not.((epi(i).lt.13.6).or.(epi(i).gt.1.36e+4)
     $        .or.(i.le.1)))
     $    sum=sum+(zremsi(i)+zremsi(i-1))*(epi(i)-epi(i-1))/2.
         enddo
c
      const=xlum/sum/ergsev
      sum2=0.
      do  i=1,numcon
         zremsz(i)=zremsz(i)+zremsi(i)*const
         if (i.gt.1)
     $    sum2=sum2+(zremsz(i)+zremsz(i-1))*(epi(i)-epi(i-1))/2.
         if (lpri.gt.1)
     $        write (lun11,*)i,epi(i),zremsi(i),const,zremsz(i)
         enddo
      sum2=sum2*ergsev
c      write (lun11,*)'normalization:',sum2
      lpri=lprisv
c
      return
      end
      subroutine step(ectt,emult,epi,ncn2,opakc,rccemis,fline,
     $  zrems,lpri,delr,dpthc,r,
     $  xpxcol,xcol,xpx,taumax,numrec0,lun11)
c
c
c
c     this routine claculates step sizes
c
      implicit none
c
      include './PARAM'
c
      real*8 opakc(ncn),epi(ncn),dpthc(2,ncn)
      real*8 rccemis(2,ncn)
      real*8 zrems(4,ncn)
      real*8 fline(2,nnnl)
      integer lpri,lprisv,kl,klmn,ncn2,numrec0,lun11
      real*8 taumax,emult,xpx,xcol,xpxcol,ectt,delr,r,optp2,dell,tst,
     $     delrmn,fpr2,r19,rmax
c
c
      lprisv=lpri
c      lpri=1
      if (lpri.ge.1) write (lun11,*)'in step',taumax
c
c
      rmax=xpxcol/xpx
c      delr=rmax/float(max(1,numrec0))
      delr=min(rmax,r/numrec0)
      r19=r*(1.e-19)
      fpr2=12.56*r19*r19
      klmn = 1
      do  kl = 1,ncn2
         optp2 = max(opakc(kl),1.e-24)
         dell=max(optp2*zrems(1,kl),
     $    (rccemis(1,kl)+rccemis(2,kl))*fpr2)
         tst = emult*zrems(1,kl)/(abs(dell)+1.e-24)
         tst = emult/optp2
         if ((epi(kl).gt.ectt).and.(dpthc(1,kl).le.taumax)
     $     .and.(zrems(1,kl).gt.1.e-12)) then
            if ( tst.lt.delr ) klmn = kl
            delr = min(delr,tst)
            endif
         if (lpri.ne.0) write (lun11,*)kl,epi(kl),opakc(kl),zrems(1,kl),
     $          rccemis(1,kl),rccemis(2,kl),fline(1,kl),dell,tst,
     $          dpthc(1,kl),delr
         enddo
      delrmn=(xpxcol-xcol)/xpx
      delr=min(delr,delrmn)
      if ( lpri.ge.1 ) write (lun11,*)'in step',emult,
     &   delr,epi(klmn),ectt,rmax,xpxcol,xpx,xcol,delrmn
      lpri=lprisv
c
c
      return
      end
      subroutine stpcut(ldirt,lpri,lun11,vturbi,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,ncsvn,
     $      epi,ncn2,opakc,oplin,opakab,delr,t,
     $      dpthc,tau0,tauc,eliml,elimh)
c
c     this routine updates.  calculates depths, etc.
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
      integer nbtpp
      parameter (nbtpp=10000)
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     master data
      integer idat1(nidat1),
     $      nptrs(nptt,ndat2)
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
      integer nplin(nnnl),nplini(nnnl)
      real*8 epi(ncn),opakc(ncn),dpthc(2,ncn)
      real*8 tauc(2,nnml)
      real*8 opakab(nnml)
c     line opacities
      real*8 oplin(nnnl)
c     line optical depths
      real*8 tau0(2,nnnl)
      real*8 optmp(ncn)
      integer ldon(2)
c
      real*8 vturbi,delr,t,eliml,elimh,
     $  dpcrit,bbb,optpp,delea,aatmp,elin,etmp,vth,
     $  vturb,deleturb,deleth,dele,aasmall,
     $  deleused,deleepi,delet,deletpp,e00,dpthmx,
     $  dpthtmp,e0,tst,opsum,optmpo,profile,optp2,
     $  tmpopmx,tmpopo,etptst
      integer ldirt,lpri,lun11,
     $  ncn2,llk,lnn,ln,ml,lup,nilin,nelin,
     $  nbtmp,iion,nitmp,ndtmp,mllz,
     $  iltmp,nlsvn,nlsv,i,ij,lind,
     $  lcon,ldir,mlm,ml1m,ltyp,lrtyp,ml1,ml2,ml1min,
     $  ml1max,mlc,mloff,mlmin,mlmax,ncut,nidt,
     $  nkdt,nrdt,np2,ncsvn
      integer nbinc,mlpar
      real*8 voigte
c     arrays containing printout info
      integer ilsv(nnnl)
      real*8 ewsv(nnnl)
      integer np1i,np1r,np1k
      real*8 oplin2(nnnl)
c
c     temporary grid for use in calculating profile
      real*8 etpp(nbtpp),optpp2(nbtpp)
c
c      real*8  tmpew,tmpewo,tmpop,tmpe,sum,sume
      real*8 tmpew,tmpewo,tmpop,tmpe,sum,sume,rnormchk
c
c     Not used
      integer javi
c
      javi=np2
      javi=npfirst(1)
      javi=nplini(1)
      javi=npcon(1)
      javi=npconi(1)
      javi=npilev(1,1)
      javi=npilevi(1)
      javi=npconi2(1)
      npconi2(1)=javi
c
c     if (lpri.ne.0)
c    $ write (lun11,*)'in stpcut:',delr,nlsvn
c
      do llk=1,ncn2
        optmp(llk)=0.
        enddo
c
c
c     calculate continuum depths
      if (lpri.ne.0) write (lun11,*)'calculating depths in stpcuta'
      lind=1
      if (ldirt.gt.0) lind=2
      dpthmx=0.
      optpp=0.
      do  i = 1,ncn2
c         opakc(i)=opakc(i)+optmp(i)
         optpp=min(optmp(i),1.e+3/(1.e-24+delr))
         dpthtmp=(opakc(i)+optpp)*delr
         if (lpri.ne.0) write (lun11,*)i,epi(i),opakc(i),optmp(i),
     $        dpthc(lind,i),dpthtmp
         dpthc(lind,i) = dpthc(lind,i) + dpthtmp
         if (dpthtmp.gt.dpthmx) then
           dpthmx=dpthtmp
           endif
         enddo
c
c     calculate line depths
      do  i = 1,nlsvn
         tau0(lind,i) = tau0(lind,i) + oplin(i)*delr
         if (lpri.ne.0)
     $    write (lun11,*)i,lind,oplin(i),delr,tau0(lind,i)
         enddo
c
c     calculate level depths
c      write (lun11,*)'in stpcut:',delr,lind
      do i = 1,ncsvn
         tauc(lind,i) = tauc(lind,i) + opakab(i)*delr
c         write (lun11,*)i,opakab(i),tauc(lind,i)
         enddo
c
c
      return
      end
      subroutine szcoll(ni,nj,tt,rate,ic)
c
c     calculates electron impact excitation rates from semiempirical
c     formula (eq.35) from smpson & zhang (1988, apj 335, 516)
c       real*8  abethe(11), hbethe(11), rbethe(11)
c       real*8  fvg1(5),fvg2(5),fvg3(5)
c     author:  M. Bautista
c
       implicit none
c
       real*8 abethe(11), hbethe(11), rbethe(11)
       real*8 fvg1(5),fvg2(5),fvg3(5)
       integer ni,nj,ic,i
       real*8 tt,rate,eion,const,rn2,g1,g2,g3,xx,gaunt,an,hn,rrn,
     $      ann,dnn,cnn,yy,e2,e3,eint1,fnn
c
       data(abethe(i),i=1,11)/ 1.30, 0.59, 0.38, 0.286, 0.229, 0.192,
     1       0.164, 0.141, 0.121, 0.105, 0.100 /
       data(hbethe(i),i=1,11)/ 1.48, 3.64, 5.93, 8.32, 10.75, 12.90,
     1       15.05, 17.20, 19.35, 21.50, 2.15 /
       data(rbethe(i),i=1,11)/ 1.83, 1.60, 1.53, 1.495, 1.475, 1.46,
     1       1.45, 1.45, 1.46, 1.47, 1.48 /
       data(fvg1(i),i=1,5)/ 1.133, 1.0785, 0.9935, 0.2328, -0.1296/
       data(fvg2(i),i=1,5)/ -0.4059, -0.2319, 0.6282, -0.5598,0.5299/
       data(fvg3(i),i=1,5)/ 0.07014, 0.02947, 0.3887, -1.181, 1.47/
c
       eion=1.578203e+5
       const=8.63e-6
c
       rn2=(float(ni)/float(nj))**2
c  computes fvalue as in johnson 1972)
       g1=0.
       g2=0.
       g3=0.
       if (ni.eq.1 ) then
        g1=fvg1(1)
        g2=fvg2(1)
        g3=fvg3(1)
       endif
       if (ni.eq.2 ) then
        g1=fvg1(2)
        g2=fvg2(2)
        g3=fvg3(2)
       endif
       if (ni.ge.3) then
        g1=fvg1(3)+fvg1(4)/ni+fvg1(5)/ni/ni
        g2=(fvg2(3)+fvg2(4)/ni+fvg2(5)/ni/ni)/ni*(-1.)
        g3=(fvg3(3)+fvg3(4)/ni+fvg3(5)/ni/ni)/ni/ni
       endif
       xx=1.-rn2
       gaunt=g1+g2/xx+g3/xx/xx
       fnn=1.9603*gaunt/(xx**3)*ni/(nj**3)
c
       if (ni.lt.11) then
         an=abethe(ni)
         hn=hbethe(ni)
         rrn=rbethe(ni)
       else
         an=abethe(11)/float(ni)
         hn=hbethe(11)*float(ni)
         rrn=rbethe(11)
       endif
       ann=fnn*4.*(ni**4)/(1.-rn2)
       dnn=ann*hn*((1.-rn2)**rrn - an*rn2)
c       write (lun11,*)ann,hn,rn2,rrn,an,fnn,dnn
       cnn=1.12*ni*ann*(1.-rn2)
       if ((nj-ni).eq.1) cnn=cnn*exp(-0.006*((ni-1)**6)/ic)
       yy=eion*ic*ic*(1./float(ni*ni)-1./float(nj*nj))/tt
       call eint(yy,eint1,e2,e3)
       rate=const/sqrt(tt)/ni/ni/ic/ic*(dnn*exp(-yy)+(ann+
     1   yy*(cnn-dnn))*eint1)
c       write (lun11,*)'in szcoll:',const,ni,nj,tt,ic,dnn,yy,
c     $     ann,cnn,dnn,eint1
c
       return
       end
      subroutine szirc(nn,T,rz,rno,cii,lpri,lun11)
c
c     calculates electron impact ionizition rates from semiempirical
c     formula (eq.35) from Smpson & Zhang (1988, ApJ 335, 516)
c     author:  M. Bautista
c
       implicit none
c
       real*8 abethe(11), hbethe(11), rbethe(11)
       integer nn,i
       real*8 t,rz,rno,cii,boltz,eion,const,rc,an,hn,rrn,tt,rn,yy,
     $      e1,e2,e3
       integer lpri,lun11
c
       DATA(abethe(i),i=1,11)/ 1.134, 0.603, 0.412, 0.313, 0.252,
     1       0.211, 0.181, 0.159, 0.142, 0.128, 1.307 /
       DATA(hbethe(i),i=1,11)/ 1.48, 3.64, 5.93, 8.32, 10.75, 12.90,
     1       15.05, 17.20, 19.35, 21.50, 2.15 /
       DATA(rbethe(i),i=1,11)/ 2.20, 1.90, 1.73, 1.65, 1.60, 1.56,
     1       1.54, 1.52, 1.52, 1.52, 1.52 /
c
       Boltz=1.38066e-16
       Eion=2.179874e-11
       const=4.6513e-3
C
       rc=float(int(rno))
       if (nn.lt.11) then
         an=abethe(nn)
         hn=hbethe(nn)
         rrn=rbethe(nn)
       else
         an=abethe(11)/float(nn)
         hn=hbethe(11)*float(nn)
         rrn=rbethe(11)
       endif
       tt= T*Boltz
       rn=float(nn)
c      yy=rz*rz/(rn*rn)*Eion/tt
       yy=rz*rz*Eion/tt*(1./rn/rn-1./rc/rc-.25*(1./(rc-1.)**2-
     c    1./rc/rc))
       call eint(yy,e1,e2,e3)
       cii=const*sqrt(tt)*(rn**5)/(rz**4)*an*yy* (
     1   e1/rn-(exp(-yy)-yy*e3)/(3.*rn)+(yy*e2-2.*yy*e1+exp(-yy))*
     2   3.*hn/rn/(3.-rrn)+(e1-e2)*3.36*yy)
      if (lpri.ne.0) write (lun11,*)'in szirc',nn,t,rz,rno,an,hn,
     $  rrn,tt,rn,yy,e1,e2,e3,cii
c
c       write (lun11,*)YY,1./yy,cii
C      give result
c
      end
      subroutine szirco(nn,t,rz,cii)
c
c     calculates electron impact ionizition rates from semiempirical
c     formula (eq.35) from smpson & zhang (1988, apj 335, 516)
c     author:  M. Bautista
c
c
       implicit none
c
       real*8 abethe(11), hbethe(11), rbethe(11)
       integer nn,i
       real*8 t,rz,cii,boltz,eion,const,an,hn,rrn,tt,rn,yy,
     $      e1,e2,e3,term1,term2,term3
c
       data(abethe(i),i=1,11)/ 1.134, 0.603, 0.412, 0.313, 0.252,
     1       0.211, 0.181, 0.159, 0.142, 0.128, 1.307 /
       data(hbethe(i),i=1,11)/ 1.48, 3.64, 5.93, 8.32, 10.75, 12.90,
     1       15.05, 17.20, 19.35, 21.50, 2.15 /
       data(rbethe(i),i=1,11)/ 2.20, 1.90, 1.73, 1.65, 1.60, 1.56,
     1       1.54, 1.52, 1.52, 1.52, 1.52 /
c
       boltz=1.38066e-16
       eion=2.179874e-11
       const=4.6513e-3
c
       if (nn.lt.11) then
         an=abethe(nn)
         hn=hbethe(nn)
         rrn=rbethe(nn)
       else
         an=abethe(11)/float(nn)
         hn=hbethe(11)*float(nn)
         rrn=rbethe(11)
       endif
       tt= t*boltz
       rn=float(nn)
c       rz=float(nz)
       yy=rz*rz/(rn*rn)*eion/tt
       call eint(yy,e1,e2,e3)
       term1=e1/rn-(exp(-yy)-yy*e3)/(3.*rn)
       term2=(yy*e2-2.*yy*e1+exp(-yy))*3.*hn/rn/(3.-rrn)
       term3=(e1-e2)*3.36*yy
       cii=const*sqrt(tt)*(rn**5)/(rz**4)*an*yy* (
     1   term1+term2+term3)
c       write (lun11,*)'in szirc:',nn,t,an,hn,rrn,rn,yy,e1,e2,e3,term1,
c     $    term2,term3,cii
c      give result
c
      end
      subroutine trnfrc(lpri,lun11,ldir,
     $      r,xpxcol,xpx,
     $      epi,ncn2,zremsz,dpthc,opakc,
     $      zrems,bremsa,bremsint,flinel)
c
c     this routine calculates continuum transfer
c     author:  T. Kallman (from xstar1)
c
      implicit none
c
      include './PARAM'
c
      real*8 epi(ncn),zremsz(ncn),dpthc(2,ncn),bremsa(ncn)
      real*8 bremsint(ncn),opakc(ncn),
     $          zrems(4,ncn),flinel(ncn)
      integer lpri,lun11,ldir,ncn2,jkp,jk,ncnm
      real*8 r,xpxcol,xpx,r19,fpr2,sum,rmax,sumtmp
c
      r19=r/(1.e+19)
      fpr2=(12.56)*r19*r19
      ncnm=ncn2-1
      bremsa(ncn2)=0.
      bremsa(ncnm)=0.
      bremsint(ncn2)=0.
      bremsint(ncnm)=0.
      sum=0.
      rmax=xpxcol/xpx
      if (lpri.ne.0) write (lun11,*)'in trnfrc:',rmax,xpxcol,xpx
      do 1 jkp=1,ncnm
         jk=ncnm+1-jkp
c
c        for outward only
c
         if (ldir.lt.0) then
             bremsa(jk)=zrems(1,jk)/fpr2
c     $                      +flinel(jk)
           else
             bremsa(jk)=zremsz(jk)*exp(-dpthc(1,jk))/fpr2
           endif
         sumtmp=(bremsa(jk)+bremsa(jk+1))*(epi(jk+1)-epi(jk))/2.
         bremsint(jk)=bremsint(jk+1)+sumtmp*(1.602197e-12)
         if (lpri.ne.0) write (lun11,*)jk,epi(jk),dpthc(1,jk),
     $       zremsz(jk),opakc(jk),
     $       zremsz(jk)*exp(-dpthc(1,jk))/fpr2,bremsa(jk)
     $       ,bremsint(jk)
1        continue
c
      return
      end
      subroutine trnfrn(lpri,lun11,
     $       nlsvn,ncsvn,ncn2,
     $       zrems,zremso,elumab,elumabo,elum,elumo)
c
c     this routine updates escaping continua and lines
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'

      real*8 elumab(2,nnml),elumabo(2,nnml)
      real*8 zrems(4,ncn),zremso(4,ncn)
      real*8 elum(3,nnnl),elumo(3,nnnl)
      integer numcon,ncn2,kl,ll,jkk,jk,lpri,lun11,
     $     ncsvn,nlsvn
c
c     transfer continuum
      if (lpri.ge.1) write (lun11,*)'in trnfrn'
      numcon=ncn2
      do kl=1,numcon
        do ll=1,4
          zremso(ll,kl) = zrems(ll,kl)
          enddo
        enddo
c
c     transfer lines
      do jkk=1,nlsvn
        jk=jkk
        do ll=1,2
          elumo(ll,jk)=elum(ll,jk)
          if (lpri.ge.1) write (lun11,*)jk,elum(1,jk),elumo(1,jk)
          enddo
        enddo
c
c     transfer RRCs
      do jkk=1,ncsvn
        jk=jkk
        do ll=1,2
          elumabo(ll,jk)=elumab(ll,jk)
          if (lpri.ge.1) write (lun11,*)jk,elumab(1,jk),elumabo(1,jk)
          enddo
        enddo
c
c
      return
      end
      subroutine ucalc(ndesc,nrdesc,ml,lcon,jkion,vturbi,
     $   nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,
     $   ans3,ans4,idest1,idest2,idest3,idest4,
     $   abund1,abund2,ptmp1,ptmp2,xpx,opakab,
     $   opakc,opakscatt,rccemis,fline,lpriu,kdesc2,
     $   rr,t,trad,tsq,xee,xh1,xh0,
     $   epi,ncn2,bremsa,bremsint,
     $   rniss,rlev,ilev,nlpt,iltp,nlev,klev,lfast,lun11,
     $          idat1,rdat1,kdat1,nptrs,np2,
     $   npar,npnxt,npfi,npfirst,
     $   nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $   npconi2,ncsvn,rates,vsav,idrates)
c
c     this routine calculates rates for all atomic processes
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
c
      integer nptmpdim
      parameter (nptmpdim=ncn)
c
      character(1) klev(100,nd)
      character(49) kdesc(ntyp),kdesc2
      character(29) krdesc(ntyp)
c     master data
      integer idat1(nidat1),nptrs(nptt,ndat2)
      integer np2
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     the saved rates
      real*8 rates(4,ndat2),vsav(4,ndat2)
      integer idrates(2,ndat2)
      real*8 epi(ncn)
      real*8 rniss(nd)
      real*8 bremsa(ncn),bremsint(ncn)
      real*8 rccemis(2,ncn),opakc(ncn),opakscatt(ncn)
      real*8 fline(2,nnnl)
      integer ilev(10,nd),nlpt(nd),iltp(nd)
      real*8 rlev(10,nd)
      real*8 aa(11),aaa(11,10),bbb(10),sg(ncn)
      real*8 rstorey(5),dcfe(8),defe(8),alhe(2),alh(2)
      real*8 etmpp(nptmpdim),stmpp(nptmpdim),ttmp(400),xsec(100)
      real*8 stmpe(nptmpdim)
      real*8  zc,eion,far,gam,scal,etmp8(ncn),stmp8(ncn)
      real*8 scal2
      real*8 a, aa1,aarec, aax, abund1, abund2, adi, aij,
     &     airt, al, algt, alm, alp, alph
      real*8 alpha, alppp, ans1, ans1o, ans2, ans2d, ans3,
     &     ans3d, ans4, ans4s, ansar2, ansar2o, ap, arad, atan, atmp, b,
     &     bb
      real*8 bbb2, bbrec, bbx, bdi, beta, bethe, c, beth,
     &     cai, ccrec, ccx, ch, ch2, ch3, chi, chir, chitmp,
     &     cii
      real*8 cij, cijpp, cion, citmp1, citmp2, cji, clu, cn, cno, crate,
     &     crec, crit53, csum, csum2, cul, d, ddd, ddx, delea
      real*8 del1, del2, dele, delev, delt, den, dirt, dirtemp,
     &     e, e0, e1, eai, ecm, ediff, ee1exp, ee1expo, eelo, eeup
      real*8 eex, eexc, efnd, eij, eijry, ekt, elammu, elin, elo, em1,
     &     em2ph, emax, enelec, ener, enn, ep, epii, erel, ergsev,eijkev
      real*8 eta, eth, etkh, etmp, ett, ett2, ettry, eup, exp10,
     &     expo, exptmp, f1, f2, fchi, ff, ff2, fh2lke, fi, flin
      real*8 float, fudge, gamma, gflin, ggl, gglo, ggu, ggup, hecxrt,
     &     hij, opakab, texp, flinabs, opakb1,
     $     p1, p2, p3, p4, p5, phi, phi1, phi2
      real*8 pi, pp, ppp, psi, ptmp1, ptmp2, q2, qq, r19, rate,
     &     rcemsum, rctmp1, rctmp2, rec,
     &     rinf,rcem1,rcem2,rnist
      real*8 rm, rr, rrrt, rs, s0, scale, sd, se, sg0,
     &     sgth, sigma, sigvtherm, sqrt, sscal, sth, sum
      real*8 swrat, t, t0, t1, t3s2, t6, tbig, temp, tfnd,
     &     time1, time2, tk, tm, tmr, tq, trad
      real*8 tsq, tst, ttz, tz,
     &     upsil, vth, vtherm, vturb, vturbi, wav, xee, xh0, xh1, xhe1
      real*8 xkt, xkto, xnx, xpx, xx, y, ya, ypow, yw, ywsq, yy, z1,
     &     zap, zeff, zz, zzz, y0,y1,yyqq
      real*8 dc,dt4,t2,term1,term2,term3,optst,opcrit,hcxrt
      real*8 er,ee1,ee2,ee3,er1,co,z2s,qij,sig,cr,crp,cr1,
     $     tmin,tmax,alphamilne,amilnerr
      real*8 tstr(100),cstr(100)
      real*8 tt,e3,e2,rho,ee
c      real*8 min       !jg
      integer nspline
      integer i57, ic, idest1, idest2, idest3, idest4,
     &     ierr, ik, il, iltmp, int, iq, ist,
     &     itmp, iz
      integer jj, jkion, jkk, jkk2, jkk3, jkkl, jlo, kdim, kl, l2, lcon,
     &     lcon2, lf, lfast, lfastl, lfasto, lff,  lforce, li, li1,lii
      integer lk, ll, lm, lorb, lpri, lprib, lpric, lpril, lprim,
     &     lprisv, lprit, lpriu, lrcalc, lrtyp, lrtyp2, lskp, ltyp,
     &     ltyp2, lun11, luse8, nkdti
      integer lz, m, ml, ml2, ml3, mlion, mllz, mlp,
     &     mm, mm5, mml, n, na, nb1, nbinc, nbmx,mlm
      integer ncn2, ncsvn, ndesc, ndtmp, nelin, nf, ni,
     &     nidt, nidt2, nidti, nilin, nind, nistage, njj, nptmp
      integer nkdt, nkdt2, nlev, nlevp, nll, nlsvn, nmin, nmx,
     &     nn,  nnz, nphint, npr, nprn, ndtmpo
      integer nq, nrdt, nrdti, nrdesc, nrdt2, nsh, nskp, ntcs,
     &     ntmp, ntmp2, nu, numcon2, nzel, nterm
      integer lunsv,lfnd
      integer np1r,np1i,np1k,np1r2,np1i2,np1k2
      integer lctype,ncase,npts
c
c     Not used
      real*8 javir
      integer javi
c      character(80) javik
c
      save aa,bb,ddd,ett,ggup,gglo,hij,opcrit,
     $         swrat,elin,pi,c,ergsev,etmp8,luse8
c
c      data opcrit/1.e-39/
      data opcrit/1.e-26/
      data ergsev/1.602197e-12/
      data pi/3.1415927/,c/2.997925e10/,luse8/0/
      data krdesc(1)/'ground state ionization      '/
      data krdesc(2)/'level ionization/recombinatio'/
      data krdesc(3)/'bound-bound collision        '/
      data krdesc(4)/'bound-bound radiative        '/
      data krdesc(5)/'bound-free collision (level) '/
      data krdesc(6)/'total recombination          '/
      data krdesc(8)/'total recombination          '/
      data krdesc(7)/'bound-free radiative (level) '/
      data krdesc(9)/'2 photon decay               '/
      data krdesc(11)/'element data                 '/
      data krdesc(12)/'ion data                     '/
      data krdesc(13)/'level data                   '/
      data krdesc(23)/'collisional superlevel->spect'/
      data krdesc(14)/'radiative superlevel->spect  '/
      data krdesc(15)/'CI total rate                '/
      data krdesc(40)/'CI from superlevels          '/
      data krdesc(41)/'non-radiative auger transtion'/
      data krdesc(42)/'Inner shell photoabsorption  '/
      data kdesc(1)/'radiative recombination:  aldrovandi and pequign '/
      data kdesc(2)/'charge exch. h0: Kingdon and Ferland             '/
      data kdesc(3)/'autoionization: hamilton, sarazin chevalier      '/
      data kdesc(4)/'line data radiative: mendosa; raymond and smith  '/
      data kdesc(5)/'2 photon transition collisional                  '/
      data kdesc(6)/'level data                                       '/
      data kdesc(7)/'dielectronic recombination: aldrovandi and pequi '/
      data kdesc(8)/'dielectronic recombination: arnaud and raymond   '/
      data kdesc(9)/'charge exch. H0 Kingdon and Ferland              '/
      data kdesc(10)/'charge exchange H+ Kingdon and Ferland          '/
      data kdesc(11)/'2 photon radiative                              '/
      data kdesc(12)/'photoionization, excited levels: hydrogenic     '/
      data kdesc(13)/'element data:                                   '/
      data kdesc(14)/'ion data:                                       '/
      data kdesc(15)/'photoionization: barfield koontz and huebner    '/
      data kdesc(16)/'arnaud and raymond ci                           '/
      data kdesc(17)/'collisional excitation hydrogenic: cota         '/
      data kdesc(18)/'radiative recombination hydrogenic: cota        '/
      data kdesc(19)/'photoionization: hullac                         '/
      data kdesc(20)/'charge exchange H+ Kingdon and Ferland          '/
      data kdesc(21)/'pixc bkh continued 3                            '/
      data kdesc(22)/'dielectronic recombination: storey              '/
      data kdesc(23)/'photoionization, excited levels: clark          '/
      data kdesc(24)/'pi xc clark continued                           '/
      data kdesc(25)/'collisional ionization: raymond and smith       '/
      data kdesc(26)/'collisional ionization hydrogenic: cota         '/
      data kdesc(27)/'photoionization: hydrogenic                     '/
      data kdesc(28)/'line data collisional: mendosa; raymond and smi '/
      data kdesc(29)/'collisional ionization data: scaled hydrogenic  '/
      data kdesc(30)/'radiative recombination hydrogenic: gould and t '/
      data kdesc(31)/'line data no levels                             '/
      data kdesc(32)/'collisional ionization: cota                    '/
      data kdesc(33)/'line data collisional: hullac                   '/
      data kdesc(34)/'line data radiative: mendosa; raymond and smitha'/
      data kdesc(35)/'photoionization: table (from bkh)               '/
      data kdesc(36)/'photoionization, excited levels:hydrogenic(no l)'/
      data kdesc(37)/'iron 3pq dr data from badnell                   '/
      data kdesc(38)/'total rr  from badnell amdpp.phys.strath.ac.uk  '/
      data kdesc(39)/'total dr  from badnell amdpp.phys.strath.ac.uk  '/
      data kdesc(40)/'                                                '/
      data kdesc(41)/'                                                '/
      data kdesc(42)/'                                                '/
      data kdesc(43)/'total photoionization cross sections tabulated  '/
      data kdesc(44)/'                                                '/
      data kdesc(45)/'                                                '/
      data kdesc(46)/'                                                '/
      data kdesc(47)/'                                                '/
      data kdesc(48)/'                                                '/
      data kdesc(49)/'op pi xsections for inner shells                '/
      data kdesc(50)/'op line rad. rates                              '/
      data kdesc(51)/'op and chianti line coll rates                  '/
      data kdesc(52)/'same as 59 but rate type 7                      '/
      data kdesc(53)/'op pi xsections                                 '/
      data kdesc(54)/'h-like cij, bautista (hlike ion)                '/
      data kdesc(55)/'hydrogenic pi xsections, bautista format        '/
      data kdesc(56)/'tabulated collision strength, bautista          '/
      data kdesc(57)/'effective charge to be used in coll. ion.       '/
      data kdesc(58)/'hlike rec rates, bautista                       '/
      data kdesc(59)/'verner pi xc                                    '/
      data kdesc(60)/'calloway h-like coll. strength                  '/
      data kdesc(62)/'calloway h-like coll. strength                  '/
      data kdesc(61)/'h-like cij, bautista (non-hlike ion)            '/
      data kdesc(63)/'h-like cij, bautista (hlike ion)                '/
      data kdesc(64)/'hydrogenic pi xsections, bautista format        '/
      data kdesc(65)/'effective charge to be used in coll. ion.       '/
      data kdesc(66)/'Like type 69 but, data in fine structure.       '/
      data kdesc(67)/'Effective collision strengths from Keenan et al.'/
      data kdesc(68)/'coll. strength He-like ions by Zhang & Sampason '/
      data kdesc(69)/'Kato & Nakazaki (1996) fit to Helike coll. strgt'/
      data kdesc(70)/'Coefficients for phot x-section of suplevels    '/
      data kdesc(71)/'Transition rates from superlevel to spect. lvls '/
      data kdesc(72)/'Autoinization rates (in s^-1) for satellite lvls'/
      data kdesc(73)/'Fit to coll. strengths satellite lvls Helike ion'/
      data kdesc(74)/'Delta functions to add to phot. x-sections  DR  '/
      data kdesc(75)/'autoionization data for Fe XXiV satellites      '/
      data kdesc(76)/'2 photon decay                                  '/
      data kdesc(77)/'coll rates from 71                              '/
      data kdesc(78)/'Auger level data                                '/
      data kdesc(79)/'fluorescence line data                          '/
      data kdesc(80)/' Collisional ionization rates gnd of Fe and Ni  '/
      data kdesc(81)/' Bhatia Fe XIX collision strengths              '/
      data kdesc(82)/' Fe UTA rad rates                               '/
      data kdesc(83)/' Fe UTA level data                              '/
      data kdesc(84)/' Iron K Pi xsections, spectator Auger binned    '/
      data kdesc(85)/' Iron K Pi xsections, spectator Auger summed    '/
      data kdesc(86)/' Iron K Auger data from Patrick                 '/
      data kdesc(88)/' Iron inner shell resonance excitation (Patrick)'/
      data kdesc(91)/' aped line wavelengths same as 50               '/
      data kdesc(92)/' aped collision strengths                       '/
      data kdesc(93)/' OP PI xsections?                               '/
      data kdesc(94)/' OP PI xsections?                               '/
      data kdesc(95)/' Bryans CI rates                                '/

      javir=trad
c      trad=javir
      javi=nlpt(1)
      javi=iltp(1)
c      javik=klev(1,1)
      javi=nlsvn
c      nlsvn=javi
      javi=npcon(1)
      javi=npconi(1)
      javi=npilev(1,1)
      javi=npilevi(1)
      javi=npconi2(1)
      javi=ncsvn
      javi=idrates(1,1)
c      javik=krdesc(1)
c
      call remtms(time1)
c
      xnx=xpx*xee
c
      lpri=lpriu
      if (lpri.gt.1)
     $  write (lun11,*)'in ucalc:',ndesc,lcon,nrdt,nidt,nkdt,
     $  ml,(rdat1(np1r+mm-1),mm=1,nrdt),(idat1(np1i+mm-1),mm=1,nidt),
     $  (kdat1(np1k+mm-1),mm=1,nkdt)
       if (lpri.gt.1) write (lun11,*)'in ucalc, inputs:',
     $   t,xee,xpx,xnx
c
      vturb=vturbi
c
      kdesc2=kdesc(ndesc)
c
      if (luse8.eq.0) then
        luse8=1
        do mm=1,ncn2
          etmp8(mm)=dble(epi(mm)/13.605692)
          enddo
        endif
      nlevp=nlev
      ans1=0.
      ans2=0.
      ans3=0.
      ans4=0.
      idest1=0
      idest2=0
      idest3=idat1(np1i+nidt-1)
      idest4=idat1(np1i+nidt-1)+1
      opakab=0.
      lforce=1
c
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
     $  17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,
     $  36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,
     $  56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,
     $  76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95),
     $  ndesc
c
c
c     rr, a&p formula
1     continue
c      write (lun11,*)'in ucalc, ndesc=1'
      arad=rdat1(np1r)
      eta=rdat1(np1r+1)
      rrrt=arad/t**eta
      ans1=rrrt*xnx
      idest1=1
      idest2=0
c      write (lun11,*)'in ucalc, ndesc=1'
      go to 9000
c
c     h charge exchange recombination
 2    continue
      aax=rdat1(np1r)
      bbx=rdat1(np1r+1)
      ccx=rdat1(np1r+2)
      ddx=rdat1(np1r+3)
      rate=aax*expo(log(t)*bbx)*max(0.,(1.+ccx*expo(ddx*t)))*(1.e-9)
      ans1=rate*xh0
c      if (lpri.ge.1) write (lun11,*)'note turning off charge exchange'
c      ans1=0.
      ans2=0.
      if (nrdesc.eq.5) then
        ans2=rate*xh0
        ans1=0.
        endif
      idest1=1
      idest2=nlevp
      if (lpri.gt.1) write (lun11,*)'type 2 data',aax,bbx,ccx,ddx,rate,
     $                               xh0,ans1,idest2
      go to 9000
c      beth=rdat1(np1r)
c      alh(1)=rdat1(np1r+1)
c      alh(2)=rdat1(np1r+2)
c      ntcs=2
c      if (t.lt.1.) ntcs=1
c      hcxrt = beth*t**alh(ntcs)
c      xh1 = xiin(1)*xpx
c      xh1=0.
c      xh2 =max(0.,(1.-xiin(1)))*xpx
c      ans1=hcxrt*xh1
c      idest1=1
c      idest2=0
c      go to 9000
c
 3    continue
c     autoionization rates
      ekt = t*(0.861707)
      cai=rdat1(np1r)
      eai=rdat1(np1r+1)
      airt = cai*expo(-eai/ekt)/tsq
      ans1=airt*xnx
      idest1=1
      idest2=1
      go to 9000
c
 4    continue
c     line rates, coll and rad
c           write (lun11,*)'level data'
c           do 1906 ll=1,nlev
c             write (lun11,*)ll,(rlev(mm,ll),mm=1,3),
c     $          (ilev(mm,ll),mm=1,3),(klev(mm,ll),mm=1,3)
c 1906        continue
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      elin=abs(rdat1(np1r))
      flin=rdat1(np1r+1)
c      if (flin.le.1.e-10) flin=1.
      eeup=rlev(1,idest1)
      eelo=rlev(1,idest2)
      if (eeup.lt.eelo) then
         itmp=idest1
         idest1=idest2
         idest2=itmp
         endif
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      ggup=rlev(2,idest1)
      gglo=rlev(2,idest2)
      a=rdat1(np1r+4)
      hij=elin*1.e-8
      elammu=elin*1.e-4
      aij=(6.67e+7)*gglo*flin/ggup/elammu/elammu
c     this is a fudge to avoid badnumerics from fine structure.
      if (flin.le.1.01e-12) aij=1.e+5
      if (elin.ge.1.e+9) aij=1.e+5
      ans1=aij*(ptmp1+ptmp2)
      ans4=aij*(ptmp1+ptmp2)*ergsev*12398.41/abs(elin)
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5)
      sigma=(0.02655)*flin*elin*(1.e-8)/vtherm
      sigvtherm=sigma
      ener=12398.41/abs(elin)
      nb1=nbinc(ener,epi,ncn2)
      ans2=0.
c      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10
      if (elin.gt.0.99e+9) then
         ans2=0.
         sigvtherm=0.
         endif
      ans1=ans1+ans2*ggup/(1.e-36+gglo)
      opakab=sigvtherm
      ans3=ans2*ener*ergsev
      ans4=ans1*ener*ergsev
c      write (lun11,*)'ltyp=4',idest1,idest2,elin,flin,ggup,gglo
      go to 9000
c
 5    continue
c     2photon rates, col
      idest1=idat1(np1i+1)
      idest2=idat1(np1i)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      ans1=0.
      ans2=0.
      ggup=rlev(2,idat1(np1i+1))
      gglo=rlev(2,idat1(np1i))
      hij=elin*1.e-8
      ekt=t*(0.861707)
      eex=abs(rlev(1,idat1(np1i))-rlev(1,idat1(np1i+1)))
      ans2=(8.629e-8)*rdat1(np1r+4)*t**rdat1(np1r+5)/ggup
      exptmp=expo(-eex/ekt)
      ans1=(8.629e-8)*rdat1(np1r+4)*t**rdat1(np1r+5)*exptmp/gglo
c      write (lun11,*)'ltyp=5',idest1,idest2,elin,flin,ggup,gglo,
c     $       ans2,eex,ekt,exptmp,ans1
      go to 9000
c      write (lun11,*)'in ucalc, ltyp=17',idat1(np1i),idat1(np1i+1),ggup,
c     $          gglo,eex,ans2,ans1
c
c     level quantities, partition function
 6    continue
      go to 9000
c
c     dr, a&p formula
 7    continue
      adi=rdat1(np1r)
      bdi=rdat1(np1r+1)
      t0=rdat1(np1r+2)
      t1=rdat1(np1r+3)
      ap=1.
      dirt=adi*ap*(1.e-06)*expo(-t0/t)
     $  *(1.+bdi*expo(-t1/t))/(t*sqrt(t))
      ans1=dirt*xnx
c      if (lpri.ne.0) write (lun11,*)'type 7 data:',
c     $  adi,bdi,t0,t1,dirt,xnx,ans1
      idest1=1
      idest2=0
      go to 9000
c
c     dr, arnaud and raymond
 8    continue
      dirt=0.
      ekt=0.861707*t
      t3s2=t**(-1.5)
      tmr = 1.e-6*t3s2
      do 820 n = 1,4
        dcfe(n)=rdat1(np1r+n-1)
        defe(n)=rdat1(np1r-1+n+4)
        dirt = dirt + dcfe(n)*expo(-defe(n)/ekt)
 820     continue
      dirt = dirt*tmr
      ans1=dirt*xnx
      idest1=1
      idest2=0
      go to 9000
c
c     he charge exchange
 9    continue
      aax=rdat1(np1r)
      bbx=rdat1(np1r+1)
      ccx=rdat1(np1r+2)
      ddx=rdat1(np1r-1+4)
      texp=min(t,1000.)**bbx
      rate=aax*texp*(1.+ccx*expo(ddx*t))*(1.e-9)
      ans2=rate*xh0
      ans1=0.
      idest1=1
      idest2=nlevp
      if (nidt.gt.1) then
        idest1=idat1(np1i)
        idest2=nlevp+idat1(np1i+1)-1
        ans2=ans2/6.
        endif
c      if (jkion.eq.18) idest1=3
      go to 9000
      bethe=rdat1(np1r)
      alhe(1)=rdat1(np1r+1)
      alhe(2)=rdat1(np1r+2)
      ntcs=2
      if (t.lt.1.) ntcs=1
      hecxrt = bethe*t**alhe(ntcs)
c      xhe1 = xiin(2)*xpx
      xhe1=0.
      ans1=hecxrt*xhe1
      idest1=1
      idest2=0
      go to 9000
c
c
 10   continue
c     charge transfer ionzation as used in func2, for level rates
      aax=rdat1(np1r)
      bbx=rdat1(np1r+1)
      ccx=rdat1(np1r+2)
      ddx=rdat1(np1r-1+4)
      eex=rdat1(np1r-1+5)
      rate=aax*t**bbx*(1.+ccx*expo(ddx*t))
     $             *expo(-eex/t)*(1.e-9)
      ans1=rate*xh1
      ans2=0.
      idest1=1
      idest2=nlevp
      go to 9000
c
 11   continue
      idest2=idat1(np1i)
      idest1=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      dele=abs(rlev(1,idest2)-rlev(1,idest1))
      ggl=rdat1(np1r+2)
      ggu=rdat1(np1r-1+4)
      ans1=(6.669e+15)*rdat1(np1r+1)*ggl/(ggu*rdat1(np1r)*rdat1(np1r))
      ans4=ans1*dele*ergsev
      go to 9000
c
 12   continue
      go to 36
c
 13   continue
      go to 9000
c
 14   continue
      go to 9000
c
 15   continue
      lprisv=lpri
      if (lpri.gt.1) write (lun11,*)'ltyp=15',ml,npar(ml)
c      if (lpri.gt.1) write (lun11,*)(rdat1(np1r-1+jj),jj=1,nrdt)
c      if (lpri.ne.0) write (lun11,*)(idat1(np1i-1+jj),jj=1,nidt)
c      if (lpri.ne.0) write (lun11,*)(kdat1(np1k-1+jj),jj=1,nkdt)
      if (ml.le.0) go to 9000
      nilin=npar(ml)
      if (nilin.le.0) go to 9000
      ntmp=nrdt/2
      do ml2=1,ntmp+1
        etmpp(ml2)=rdat1(np1r-1+2*ml2-1)
        stmpp(ml2)=rdat1(np1r-1+2*ml2)*1.e-18
        enddo
      mlm=nilin-1
      call drd(ltyp2,lrtyp2,lcon2,
     $  nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $  nptrs,0,lun11)
      idest1=idat1(np1i+nidt-2)
      idest2=idat1(np1i+nidt-3)-idat1(np1i+nidt-1)
      if (lpri.gt.1)
     $ write (lun11,*)ml,nilin,rdat1(np1r),idest1
      ett=rdat1(np1r2)
      if (lpri.gt.1)
     $ write (lun11,*)'ett=',ett
      nb1=nbinc(ett,epi,ncn2)
      gglo=rlev(2,1)
      ggup=rlev(2,nlevp)
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      d=rdat1(np1r+1)
      do mm=1,11
        aa(mm)=rdat1(np1r-1+3+mm)
        enddo
c      aa(1)=min(max(aa(1),-6.),6.)
      ekt=t*(0.861707)
      if (lpri.gt.1)
     $  write (lun11,*)'in ucalc, ind=15:',lcon,
     $               nrdt,nidt,nkdt
      if (lpri.gt.1)
     $ write (lun11,891)(rdat1(np1r-1+mm),mm=1,nrdt)
 891  format (1x,10(1pe10.3))
      if (lpri.gt.1)
     $ write (lun11,892)(idat1(np1i-1+mm),mm=1,nidt)
 892  format (1x,10(i6))
      if (lpri.gt.1)
     $ write (lun11,893)(kdat1(np1k-1+mm),mm=1,nkdt)
 893  format (1x,100a1)
      na=idat1(np1i-1+nidt-5)
      nsh=idat1(np1i-1+nidt-4)
      do lk=1,na
        ll=idat1(np1i-1+9*lk-5)
        if (lpri.gt.1)
     $   write (lun11,*)'ll=',ll,lk
        lz=15*(lk-1)
        ett=rdat1(np1r-1+1+lz)
        ddd=rdat1(np1r-1+2+lz)
        bb=rdat1(np1r-1+3+lz)
        aa(1)=rdat1(np1r-1+4+lz)
        aa(2)=rdat1(np1r-1+5+lz)
        if (lpri.gt.1)
     $   write (lun11,*)'ltest=2',ett,ddd,bb,aa(1)
        do 1011 mml=1,5
          aa(2+mml)=rdat1(np1r-1+mml+5+lz)
 1011     continue
        do 1012 mml=1,4
          aa(7+mml)=rdat1(np1r-1+mml+10+lz)
 1012     continue
        bbb(lk)=bb
        do 1013 mml=1,11
          aaa(mml,lk)=aa(mml)
          if (lpri.gt.1)
     $     write (lun11,*)mml,aa(mml)
 1013     continue
        if (lpri.gt.1)
     $   write (lun11,*)'ltest=0',ll,na,nsh,aa(7)
        enddo
      lprib=0
      if (lpri.gt.1) lprib=lpri
      if (lpri.gt.1)
     $ write (lun11,*)'calling bkhsgo:',ett,t,
     $ ddd,(bbb(mm),mm=1,3),na,
     $ aaa(1,1),aaa(7,1)
      lfastl=1
      call bkhsgo(sg,ett,ddd,bbb,na,
     $         aaa,epi,ncn2,t,lprib,lfastl,lun11)
      lprib=0
c      if (lpri.gt.1) lprib=lpri
      gglo=rlev(2,1)
      ggup=rlev(2,nlevp)
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      if (lpri.ge.1) then
        npr=nbinc(ett,epi,ncn2)
        write (lun11,*)'bkh threshold xsection:',
     $         npr,ett,sg(npr)
        endif
      lpri=lprisv
      go to 9000
c
 16   continue
      ekt=t*(0.861707)
      njj=int(nrdt/5)
      lpriu=0
      if (lpriu.ne.0)
     $ write (lun11,*)'ltyp=16:',idat1(np1i+nidt-1)
      csum=0.
      csum2=0.
      do mm=1,njj
        mm5=5*(mm-1)
        eth=rdat1(np1r-1+mm5+1)
        a=rdat1(np1r-1+mm5+2)
        b=rdat1(np1r-1+mm5+3)
        c=rdat1(np1r-1+mm5+4)
        d=rdat1(np1r-1+mm5+5)
        xx=eth/ekt
        if (lpriu.ne.0)
     $   write (lun11,*)'xx=',xx,eth,ekt,a,b,c,d,np1r,mm5
        em1=ee1expo(xx)
        f1=em1/xx
        if (lpriu.ne.0)
     $   write (lun11,*)'before ff2:',f1,em1,xx
        f2=ff2(xx,lpriu,lun11)
        if (lpriu.ne.0)
     $   write (lun11,*)xx,a,b,c,d,em1,f1,f2
        fi=a*(1.-xx*f1)+b*(1.+xx-xx*(2.+xx)*f1)
     $   +c*f1+d*xx*f2
        fi=max(fi,0.)
        csum=csum+fi*expo(-xx)/xx
        csum2=csum2+fi/xx
        if (lpriu.ne.0)
     $   write (lun11,*)mm,mm5,a,b,c,d,xx,f1,fi,csum
        enddo
      citmp1=csum*(6.69e-7)/ekt**(1.5)
      ans1=citmp1*xnx
      citmp2=csum2*(6.69e-7)/ekt**(1.5)
      idest1=1
      idest2=1
      ggup=rlev(2,nlevp)
      gglo=rlev(2,1)
c     note that rinf has exponential removed
      rinf=(2.08e-22)*gglo/ggup/t/tsq
      ans2=citmp2*xnx
      ans2=ans2*rinf*xnx
      if (nrdesc.eq.5) then
          idest2=nlevp
        else
          idest2=1
        endif
      if (lpriu.ne.0)
     $   write (lun11,*)csum,citmp1,citmp2,ans1,
     $     ggup,gglo,rinf,ans2,idest1,idest2
      go to 9000
c
 17   continue
c     line rates, col
      ans1=0.
      ans2=0.
      hij=elin*1.e-8
c      write (lun11,*)'ltyp=4',idest1,idest2,elin,flin,ggup,gglo
      ekt=t*(0.861707)
      idest1=idat1(np1i+1)
      idest2=idat1(np1i)
      eeup=rlev(1,idest2)
      eelo=rlev(1,idest1)
      if (eelo.gt.eeup) then
        idest1=idat1(np1i)
        idest2=idat1(np1i+1)
        eeup=rlev(1,idest2)
        eelo=rlev(1,idest1)
        endif
      eex=eeup-eelo
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      ans2=(8.629e-8)*rdat1(np1r)*t**rdat1(np1r+1)/ggup
      ans1=0.
      exptmp=expo(-eex/ekt)
      exptmp=1.
      if (ekt.gt.eex/20.)
     $ ans1=ans2*ggup*exptmp/gglo
      if (lpri.gt.1)
     $ write (lun11,*)'in ucalc, ltyp=17',idat1(np1i),
     $   idat1(np1i+1),ggup,
     $   gglo,eex,ans2,ans1
      go to 9000
c
 18   continue
      aarec=rdat1(np1r)
      bbrec=rdat1(np1r+1)
      ccrec=rdat1(np1r+2)
      ttz=rdat1(np1r-1+4)
      algt=log10(t/(1.e-32+ttz))+4.
      algt=max(algt,3.5)
      algt=min(algt,7.5)
      idest1=idat1(np1i)
      ans1=exp10(aarec+bbrec*(algt-ccrec)**2)/t/1.e+4
      ans1=ans1*xnx
c      ans1=0.
      ans2=0.
      idest2=0
      go to 9000
c
 19   continue
      etkh=rdat1(np1r-1+5)
      enelec=1.
      eth=etkh
      nb1=nbinc(eth,epi,ncn2)
      idest1=idat1(np1i)
      idest2=nlevp
      ggup=rlev(2,nlevp)
      gglo=rlev(2,idest1)
      swrat=gglo/ggup
      ekt=t*(0.861707)
      lm=nb1
      do while (lm.le.nphint)
         bbb2=epi(lm)/max(etkh,1.e-30)
         etmp=log(bbb2)
         alppp=rdat1(np1r)+etmp*(rdat1(np1r+1)+etmp*
     $         (rdat1(np1r+2)+etmp*rdat1(np1r-1+4)))
         ppp=expo(alppp)
         sg(lm)=(1.e-18)*enelec*ppp*(13.606)/etkh
         call enxt(eth,nb1,lpri,epi,ncn2,t,lfastl,lun11,
     $                  lm,nskp,nphint,lrcalc)
        lm=lm+nskp
        enddo
      call phintfo(sg,eth,ans1,ans2,ans3,ans4,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      go to 9000
c
 20   continue
c     charge transfer ionzation as used in func1, for total rate
      aax=rdat1(np1r)
      bbx=rdat1(np1r+1)
      ccx=rdat1(np1r+2)
      ddx=rdat1(np1r-1+4)
      eex=rdat1(np1r-1+5)
      rate=aax*t**bbx*(1.+ccx*expo(ddx*t))
     $             *expo(-eex/t)*(1.e-9)
      ans1=rate*xh1
      ans2=0.
      idest1=1
      idest2=nlevp
      idest3=idat1(np1i+nidt-1)
      idest4=idest3+1
      go to 9000
c
 21   continue
c     h charge exchange dalgarno and butler
      beth=rdat1(np1r)
      alh(1)=rdat1(np1r+1)
      alh(2)=rdat1(np1r+2)
      ntcs=2
      if (t.lt.1.) ntcs=1
      hcxrt = beth*t**alh(ntcs)
      ans1=hcxrt*xh1
      idest1=1
      idest2=0
      go to 9000
c
c     dr storey
 22   continue
      ans1=0.
      ans2=0.
      idest1=1
      idest2=0
      if (t.gt.6.) go to 9000
      do 221 kl=1,5
        rstorey(kl)=rdat1(np1r-1+kl)
 221    continue
c      if (rstorey(5).lt.0.) go to 9000
      t3s2=t**(-1.5)
      dirtemp=
     $   (1.e-12)*(rstorey(1)/t+rstorey(2)
     $   +t*(rstorey(3)+t*rstorey(4)))*t3s2
     $   *expo(-rstorey(5)/t)
      dirtemp=max(dirtemp,0.)
      if (lpri.gt.1) write (lun11,*)'in ucalc, ltyp=22:',
     $   ndesc,lcon,nrdt,nidt,nkdt,
     $  ml,(rdat1(np1r-1+mm),mm=1,nrdt),(idat1(np1i-1+mm),mm=1,nidt),
     $  (kdat1(np1k-1+mm),mm=1,nkdt),dirtemp,xnx
      ans1=dirtemp*xnx
      idest1=1
      idest2=0
      go to 9000
c
 23   continue
      lfastl=1
      lprisv=lpri
c      lpri=2
      if (lpri.gt.1)
     $ write (lun11,*)'in ucalc, 23:',rdat1(np1r),rdat1(np1r+1),
     $  rdat1(np1r+2),rdat1(np1r-1+4),
     $  rdat1(np1r-1+5),swrat,idat1(np1i),idat1(np1i+1),idat1(np1i+2),
     $  idat1(np1i-1+4),idat1(np1i-1+5)
      ett=rlev(1,idat1(np1i))
      eth=ett
      nb1=nbinc(eth,epi,ncn2)
      gglo=rlev(2,idat1(np1i))
      ggup=rlev(2,nlevp)
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      ekt=t*(0.861707)
      jkk2=idat1(np1i+nidt-1)
      nilin=npar(ml)
      if (nilin.le.0) go to 9000
      jkk3=0
      jkk=0
      ndtmp=npfirst(12)
      do while ((jkk.ne.jkk2).and.(ndtmp.ne.0))
        jkk3=jkk3+1
        mlm=ndtmp-1
        call drd(ltyp2,lrtyp2,lcon2,
     $    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $    nptrs,0,lun11)
        jkk=idat1(np1i2+nidt2-1)
        ndtmp=npnxt(ndtmp)
        enddo
      if (ndtmp.le.0) go to 9000
      zzz=float(idat1(np1i2))
      ndtmp=npfi(13,jkk3)
      mllz=npar(ndtmp)
      if (lpri.gt.1) write (lun11,*)jkk,jkk2,jkk3,zzz,ndtmp
      iltmp=1
      mlm=ndtmp-1
      call drd(ltyp2,lrtyp2,lcon2,
     $    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $    nptrs,0,lun11)
      ndtmpo=ndtmp
      ndtmp=npnxt(ndtmp)
      do while ((ndtmp.ne.0).and.(iltmp.ne.idat1(np1i))
     $      .and.(npar(ndtmp).eq.mllz))
        mlm=ndtmp-1
        call drd(ltyp2,lrtyp2,lcon2,
     $    nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $    nptrs,0,lun11)
        iltmp=idat1(np1i2+nidt2-2)
        ndtmpo=ndtmp
        ndtmp=npnxt(ndtmp)
        enddo
      ndtmp=ndtmpo
      nprn=idat1(np1i)
      enn=float(nprn)
      if ((enn.le.1.e-24).or.(zzz.le.1.e-24).or.(ett.le.1.e-6))
     $  go to 9000
      sg0=6.3e-18*enn/zzz/zzz
      if (lpri.gt.1)
     $ write (lun11,*)'ind=23:',ml,nilin,zzz,jkk,nprn,sg0,ett
      ll=nb1
      do while (ll.le.nphint)
        epii=epi(ll)
        sg(ll)=sg0*(epii/ett)**(-3)
        call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,
     $                  ll,nskp,nphint,lrcalc)
        ll=ll+nskp
        enddo
      lprib=0
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      lpri=lprisv
c      ans2=ans2*xnx
      idest1=idat1(np1i+nidt-2)
      idest2=nlevp
      go to 9000
c
 24   continue
      go to 9000
c
 25   continue
      if (nrdesc.eq.5) then
          idest2=nlevp
        else
          idest2=1
        endif
      e=rdat1(np1r)
      a=rdat1(np1r+1)
      b=rdat1(np1r+2)
      c=rdat1(np1r-1+4)
      d=rdat1(np1r-1+5)
      cion = 0.
      chir = (t*1.e+4)/(11590.*e)
      citmp1=cion
      ans1=citmp1*xnx
      ans2=0.
      idest1=1
c      idest2=1
      if ( chir.le..0115 ) go to 9000
      chi = max(chir,0.1)
      ch2 = chi*chi
      ch3 = ch2*chi
      alpha = (.001193+.9764*chi+.6604*ch2+.02590*ch3)
     &        /(1.0+1.488*chi+.2972*ch2+.004925*ch3)
      beta = (-.0005725+.01345*chi+.8691*ch2+.03404*ch3)
     &       /(1.0+2.197*chi+.2457*ch2+.002503*ch3)
      ch = 1./chi
      fchi = 0.3*ch*(a+b*(1.+ch)+(c-(a+b*(2.+ch))*ch)*alpha+d*beta*ch)
      chitmp=expo(-1./chir)
      cion = 2.2e-6*sqrt(chir)*fchi/(e*sqrt(e))
      citmp1=cion
      ans1=citmp1*xnx
      ggup=rlev(2,nlevp)
      gglo=rlev(2,nidt-1)
c     note that rinf has exponential removed
      rinf=(2.08e-22)*gglo/ggup/t/tsq
      ans2=ans1*rinf*xnx
      ans1=ans1*chitmp
c      idest1=idat1(np1i+nidt-2)
      idest4=idat1(np1i+nidt-1)+1
      idest3=idat1(np1i+nidt-1)
      go to 9000
c
 26   continue
      go to 9000
c      ekt=t*(0.861707)
c      idest1=idat1(np1i)
c      gglo=rlev(2,idest1)
c      edelt=abs(rlev(1,idest1)-rlev(1,nlev))
c      exptmp=expo(-edelt/ekt)
c      ans1=(4.1416e-9)*rdat1(np1r)*t**rdat1(np1r+1)/gglo
c      ggup=rlev(2,nlev)
c      rinf=(2.08e-22)*gglo/ggup/t/tsq
c      ans2=ans1*rinf
c      ans1=ans1*exptmp
c      write (lun11,*)'ltyp=26',idest1,gglo,ggup,
c     $   edelt,rdat1(np1r),rdat1(np1r+1),ans1
c      idest2=nlev
c      go to 9000
c
 27   continue
      lprisv=lpri
c      ett=rdat1(np1r+1)
      lfastl=1
      idest1=1
      if (nrdesc.eq.1) then
          idest2=0
        else
          idest2=nlevp
        endif
      ett=abs(rlev(1,nlev)-rlev(1,1))
c      if (lpri.gt.1)
c      write (lun11,*)'in ucalc, ind=27:',rlev(1,nlev),
c     $     rlev(1,1),nlev,ett
      if (ett.le.1.e-5) go to 9000
      eth=ett
      nb1=nbinc(eth,epi,ncn2)
      gglo=rlev(2,idat1(np1i))
      swrat=gglo
      ekt=t*(0.861707)
      ll=nb1
      do while (ll.le.nphint)
        epii=epi(ll)
        e=epii
        eth=ett
        zap = e/eth - 1.
        y = e/eth
        yy=sqrt(zap)
        yy=max(yy,1.e-04)
        fh2lke=((6.3e-18)/rdat1(np1r)/rdat1(np1r))
     $   *y**(-4)*expo(4.-4.*atan(yy)/yy)
     $   /(1.-expo(-6.2832/yy))
c        fh2lke=((6.3e-18)/rdat1(np1r)/rdat1(np1r))*y**(-3)
        sg(ll)=fh2lke
        if (lpri.ge.2) write (lun11,*)ll,epii,zap,y,yy,fh2lke
        call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,
     $                  ll,nskp,nphint,lrcalc)
        ll=ll+nskp
        enddo
      lprib=0
      if (lpri.gt.1) lprib=1
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      lpri=lprisv
      go to 9000
c
 28   continue
c     line rates, col
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if (rlev(1,idat1(np1i+1)).lt.rlev(1,idat1(np1i))) then
        idest2=idat1(np1i)
        idest1=idat1(np1i+1)
        endif
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      nind=5
      if (nrdt.ge.12) then
        do ll=1,4
          ttmp(ll)=rdat1(np1r-1+nrdt-4+ll)
          enddo
        jlo=0
        call hunt3(ttmp,4,t,jlo,0,lun11)
        nind=nrdt-8+jlo
        endif
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      elin=abs(rdat1(np1r))
      hij=elin*1.e-8
      if (elin.le.1.e-24) go to 9000
c      nind=nrdt-2
      cijpp=rdat1(np1r-1+nind)
      ekt=0.861707*t
      delt=12398.41/elin/ekt
      cji=(8.626e-8)*cijpp/tsq/ggup
      cij=0.
      exptmp=expo(-delt)
      cij=cji*ggup*exptmp/tsq/gglo
      ans1=cij*xnx
      ans2=cji*xnx
      if (lpri.gt.1) then
        write (lun11,*)'ltyp=28',idest1,idest2,elin,flin,ggup,gglo
        write (lun11,*)'       ',nrdt,(rdat1(np1r-1+mm),mm=1,8),nind,jlo
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp
        endif
      elin=0.
      go  to 9000
c
 29   continue
      go to 9000
c      anstmp=rdat1(np1r+1)*(8.626e-8)/tsq
c      ans2=anstmp*(2.08e-22)*(rdat1(np1r+2)/rdat1(np1r-1+4))/t/tsq
c      ans1=0.
c      delt=rdat1(np1r)/t
c      if (delt.lt.50.) then
c         exptmp=1.
c         exptmp=expo(-delt)
c         ans1=anstmp*exptmp
c         endif
c      idest1=idat1(np1i+1)
c      idest2=nlev
c      write (lun11,*)'ltyp=29',ans1,ans2,(rdat1(np1r-1+ii),ii=1,4),anstmp,xnx
c      go to 9000
c
 30   continue
c      write (lun11,*)'ltyp=30',idat1(np1i)
        nmx=idat1(np1i)
        t6=t/100.
        zeff=float(nmx)
        beta=zeff*zeff/(6.34*t6)
        yy=beta
        vth=(3.10782e+7)*sqrt(t)
c       fudge factor makes the 2 expressions join smoothly
        ypow=min(1.,(0.06376)/yy/yy)
        fudge=0.9*(1.-ypow)+(1./1.5)*ypow
        phi1=(1.735+log(yy)+1./6./yy)*fudge/2.
        phi2=yy*(-1.202*log(yy)-0.298)
        phi=phi1
        if (yy.lt.0.2525) phi=phi2
        rrrt=2.*(2.105e-22)*vth*yy*phi
        ans1=rrrt*xnx
        ans2=0.
        idest1=1
        idest2=0
      go to 9000
c
 31   continue
c     line rates, coll and rad
c           write (lun11,*)'level data'
c           do 1906 ll=1,nlev
c             write (lun11,*)ll,(rlev(mm,ll),mm=1,3),
c     $          (ilev(mm,ll),mm=1,3),(klev(mm,ll),mm=1,3)
c 1906        continue
      idest1=idat1(np1i+1)
      idest2=idat1(np1i)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      elin=abs(rdat1(np1r))
      flin=rdat1(np1r+1)
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
c      a=rdat1(np1r-1+5)
      ans1=0.
      ans2=0.
      if (ml.le.0) go to 9000
      nilin=npar(ml)
      if (nilin.le.0) go to 9000
      nelin=npar(nilin)
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000
      mlm=nelin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      a=rdat1(np1r2+1)
      hij=elin*1.e-8
      aij=(0.02655)*flin*8.*pi/hij/hij*gglo/(1.e-24+ggup)
      ans1=aij*(ptmp1+ptmp2)
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5)
      sigma=(0.02655)*flin*elin*(1.e-8)/vtherm
      sigvtherm=sigma
      ener=12398.41/abs(elin)
      nb1=nbinc(ener,epi,ncn2)
      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10
      if (elin.gt.0.99e+9) then
         ans2=0.
         sigvtherm=0.
         endif
      ans1=ans1+ans2*ggup/(1.e-36+gglo)
c     notice that opakab does not have abundance in
      opakab=sigvtherm
      ans3=ans2*ener*ergsev
      ans4=ans1*ener*ergsev
      delea=0.
      lfasto=4
      if (lfasto.ge.4) ans2=0.
c      if (opakab.gt.1.e-34)
c     $  call linopac(lpri,lun11,opakab,rcem1,rcem2,elin,vturb,t,a,
c     $               delea,epi,ncn2,opakc,opakscatt,rccemis,fline,
c     $               lfasto)
      ans4=ans1*ener*ergsev
c      write (lun11,*)'ltyp=31',idest1,idest2,elin,flin,ggup,gglo
      go to 9000
c
 32   continue
      idest1=idat1(np1i)
      gglo=rdat1(np1r-1+4)
      ans1=0.
      ans2=0.
      go to 9000
c      if (gglo.lt.1.e-24) go to 9000
c      ekt=t*(0.861707)
c      edelt=rdat1(np1r+2)
c      ans1=(4.1416e-9)*rdat1(np1r)*t**rdat1(np1r+1)*expo(-edelt/ekt)
c     $        /gglo
c      write (lun11,*)'ltyp=26',idest1,gglo,edelt,rdat1(np1r),rdat1(np1r+1),ans1
c      idest2=nlev
c      go to 9000
c
 33   continue
c     line rates, col
      idest1=idat1(np1i+1)
      idest2=idat1(np1i)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      ggup=rlev(2,idat1(np1i+1))
      gglo=rlev(2,idat1(np1i))
      elin=abs(rdat1(np1r))
      hij=elin*1.e-8
      if (elin.le.1.e-24) go to 9000
      nind=4
      cijpp=rdat1(np1r-1+nind)
      ekt=0.861707*t
      delt=12398.41/elin/ekt
      exptmp=expo(-delt)
      cij=(8.626e-8)*cijpp*exptmp/tsq/gglo
      cji=(8.626e-8)*cijpp/tsq/ggup
      ans1=cij*xnx
      ans2=cji*xnx
      go to 9000
c
 34   continue
c     line rates, coll and rad
c           write (lun11,*)'level data'
c           do 1906 ll=1,nlev
c             write (lun11,*)ll,(rlev(mm,ll),mm=1,3),
c     $          (ilev(mm,ll),mm=1,3),(klev(mm,ll),mm=1,3)
c 1906        continue
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      elin=abs(rdat1(np1r))
      aij=rdat1(np1r+1)
      eeup=rlev(1,idest1)
      eelo=rlev(1,idest2)
      if (eeup.lt.eelo) then
         itmp=idest1
         idest1=idest2
         idest2=itmp
         endif
      ggup=rlev(2,idest1)
      gglo=rlev(2,idest2)
c      ggup=rdat1(np1r-1+4)
c      gglo=rdat1(np1r+2)
      a=rdat1(np1r-1+5)
      hij=elin*1.e-8
      elammu=elin*1.e-4
c      flin=aij*hij*hij*ggup/((0.02655)*8.*pi*gglo)
      flin=aij*hij*hij*ggup/((0.667274)*gglo)
      ans1=aij*(ptmp1+ptmp2)
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5)
      sigma=(0.02655)*flin*elin*(1.e-8)/vtherm
      sigvtherm=sigma
      ener=12398.41/abs(elin)
      nb1=nbinc(ener,epi,ncn2)
      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10
      if (elin.gt.0.99e+9) then
         ans2=0.
         sigvtherm=0.
         endif
      ans1=ans1+ans2*ggup/(1.e-36+gglo)
c     notice that opakab does not have abundance in
      opakab=sigvtherm
      delea=0.
      lfasto=4
      if (lfasto.ge.4) ans2=0.
      ans3=ans2*ener*ergsev
      ans4=ans1*ener*ergsev
c      if (opakab.gt.1.e-34)
c     $  call linopac(lpri,lun11,opakab,rcem1,rcem2,elin,vturb,t,a,
c     $               delea,epi,ncn2,opakc,opakscatt,rccemis,fline,
c     $               lfasto)
c      if (lpri.ne.0)
c     $ write (lun11,*)'ltyp=34',idest1,idest2,elin,flin,ggup,gglo,
c     $                         a,aij,hij,pi
      go to 9000
c
 35   continue
      lprisv=lpri
      if (lpri.gt.1)
     $  write (lun11,*)'in ucalc, ind=15:',lcon
      ett=rdat1(np1r)
      if (ett.le.(1.e-24)) go to 9000
      ntmp=(nrdt-1)/2
      do ml2=1,ntmp
        etmpp(ml2)=rdat1(np1r-1+1+2*ml2)
        stmpp(ml2)=rdat1(np1r-1+2*ml2)
        enddo
      nb1=nbinc(ett,epi,ncn2)
      gglo=rlev(2,1)
      ggup=rlev(2,nlevp)
      if (ggup.le.1.e-24) then
         write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      numcon2=max(2,ncn2/50)
      nphint=ncn2-numcon2
      idest1=idat1(np1i-1+6)
      idest2=idat1(np1i-1+5)-idat1(np1i-1+7)
      ekt=t*(0.861707)
      jlo=0
      ll=nb1
      lfastl=1
      do while (ll.le.nphint)
          epii=epi(ll)
          efnd=(epii-ett)/13.605692
          call hunt3(etmpp,ntmp,efnd,jlo,0,lun11)
          ml2=jlo
          mlp=ml2+1
          del1=(efnd-etmpp(ml2))/(etmpp(mlp)-etmpp(ml2))
          del2=(efnd-etmpp(mlp))/(etmpp(mlp)-etmpp(ml2))
          sg(ll)=-stmpp(ml2)*del2+stmpp(mlp)*del1
c          if (lpri.gt.1)
c     $    write (lun11,*)ll,epii,sg(ll),ml2,stmpp(ml2),stmpp(mlp),
c     $              del1,del2
          call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,
     $                  ll,nskp,nphint,lrcalc)
          ll=ll+nskp
          enddo
      lprib=0
      if (lpri.gt.1) lprib=lpri
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      if (lpri.gt.1) then
        npr=nbinc(ett,epi,ncn2)+2
        write (lun11,*)'bkh threshold xsection:',
     $         npr,ett,sg(npr)
        endif
      lpri=lprisv
      go to 9000

c
 36   continue
c      photoionization, excited levels:hydrogenic(no l)
      lprisv=lpri
c      lpri=2
      if (lpri.gt.1)
     $  write (lun11,*)'in ucalc, ind=36:',(idat1(np1i-1+mm),mm=1,5)
      idest1=idat1(np1i+nidt-2)
      ett=rlev(1,nlevp)-rlev(1,idest1)
      idest2=nlevp
      if (ett.le.1.e-5) go to 9000
      eth=ett
      nb1=nbinc(eth,epi,ncn2)
      gglo=rlev(2,idest1)
      ggup=rlev(2,nlevp)
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      ekt=t*(0.861707)
      if (ml.le.0) go to 9000
      nilin=npar(ml)
      if (nilin.le.0) go to 9000
      nelin=npar(nilin)
      if (lpri.gt.1)
     $ write (lun11,*)'in ucalc, ind=36:',
     $   ml,nilin,nelin
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000
      mlm=nilin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      nistage=idat1(np1i2)
      mlm=nelin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      nzel=idat1(np1i2)
      nq=idat1(np1i)
      nq=min(10,nq)
      zz=float(nzel-nistage+1)
      sgth=(6.3e-18)*nq*nq/zz/zz
      if (lpri.gt.1) write (lun11,*)nb1,nq,nzel,nistage,zz,
     $                              ett,sgth,idest1,gglo,ggup
      ll=nb1
      lfastl=1
      do while (ll.le.nphint)
        epii=epi(ll)
        sg(ll)=sgth*(epii/ett)**(-3)
        call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,
     $                  ll,nskp,nphint,lrcalc)
        ll=ll+nskp
        enddo
      lprib=0
      if (lpri.gt.1) lprib=lpri
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      lpri=lprisv
      go to 9000
c
 37   continue
c     total dr for fe 3pq ions from badnell 2006 Ap. J. Lett 651 L73
      dirt=0.
      ekt=0.861707*t
      t3s2=t**(-1.5)
      tmr = 1.e-6*t3s2
      nterm=idat1(np1i)
      do  n = 1,nterm
        dcfe(n)=rdat1(np1r-1+n)
        defe(n)=rdat1(np1r-1+n+4)
        dirt = dirt + dcfe(n)*expo(-defe(n)/ekt)
        enddo
      dirt = dirt*tmr
      ans1=dirt*xnx
      idest1=1
      idest2=0
      go to 9000
c
 38   continue
c     total rr  from badnell http://amdpp.phys.strath.ac.uk/tamoc/DATA/DR/
      a=rdat1(np1r)
      b=rdat1(np1r+1)
      t0=rdat1(np1r+2)/1.e+4
      t1=rdat1(np1r-1+4)/1.e+4
      if (nrdt.gt.4) then
        c=rdat1(np1r-1+5)
        t2=rdat1(np1r-1+6)/1.e+4
        b=b+c*exp(-t2/t)
        endif
      term1=(T/T0)**(0.5)
      term2=(1.+(T/T0)**(0.5))**(1.-b)
      term3=(1.+(T/T1)**(0.5))**(1.+b)
      rrrt=a/(1.e-34+term1*term2*term3)
      ans1=rrrt*xnx
      if (lpri.gt.1) write (lun11,*)a,b,c,t0,t1,t2,
     $         term1,term2,term3,rrrt,ans1
      idest1=1
      idest2=0
      go to 9000
c
 39   continue
c     total dr  from badnell http://amdpp.phys.strath.ac.uk/tamoc/DATA/DR/
      dirt=0.
      ekt=0.861707*t
      t3s2=t**(-1.5)
      tmr = 1.e-6*t3s2
      nterm=nrdt/2
      do  n = 1,nterm
        dc=rdat1(np1r-1+n)
        dt4=rdat1(np1r-1+n+nterm)/1.e+4
        dirt = dirt + dc*exp(-dt4/t)
        if (lpri.gt.1) write (lun11,*)n,dc,dirt
        enddo
      dirt = dirt*tmr
      ans1=dirt*xnx
      if (lpri.gt.1) write (lun11,*)nterm,dirt,ans1
      idest1=1
      idest2=0
      go to 9000
c
 40   continue
      go to 9000
c
 41   continue
      go to 9000
c
 42   continue
      go to 9000
c
 43   continue
c     total photoionization cross sections tabulated in
c     format like 53 (not used)
      go to 9000
c
 44   continue
      go to 9000
c
 45   continue
      go to 9000
c
 46   continue
      go to 9000
c
 47   continue
      go to 9000
c
 48   continue
      go to 9000
c
 49   continue
 499  continue
c     op pi xsections
c     old version
      lprisv=lpri
c      if (lpri.ge.1) lpri=2
c     these are the initial and final levels and indeces
c     notice that these are relative to the current ion
c     (not relative to the element as a whole)
      idest1=idat1(np1i+nidt-2)
      idest4=idat1(np1i+nidt-3)
      idest2=nlevp+max(0,idat1(np1i-1+nidt-3))-1
      if (lpri.gt.1) write (lun11,*)'idest1=',idest1,idest2
      if ((idest1.ge.nlevp).or.(idest1.le.0)) go to 9000
      if (ml.le.0) go to 9000
      eth=rlev(4,idest1)-rlev(1,idest1)
      ett=eth
      nilin=npar(ml)
      if (lpri.gt.1) write (lun11,*)'nilin=',nilin,ml
      if (nilin.le.0) go to 9000
      ntmp=nrdt/2
      do ml2=1,ntmp+1
        etmpp(ml2)=rdat1(np1r-1+2*ml2-1)
        stmpp(ml2)=rdat1(np1r-1+2*ml2)*1.e-18
        stmpp(ml2)=max(stmpp(ml2),0.)
        enddo
c      ett=ett+max(0.,13.605692*etmpp(1))
      optst=abund1*stmpp(1)
c      if ((optst.lt.opcrit).and.(lfast.eq.2)) go to 9000
      if (lpri.gt.1) write (lun11,*)'ett=',ett,etmpp(1)
      if (ett.le.0.) go to 9000
      ntmp2=nptmpdim
      call phextrap(etmpp,stmpp,ntmp,ntmp2,ett,ncn2,lpri,lun11)
      nb1=nbinc(ett,epi,ncn2)
      tst=abs(bremsint(nb1)/max(1.e-24,vsav(1,ml))-1.)
      xkt=ett/(0.861707*t)
      r19=rr/1.e+19
c      if ((tst.le.0.01).and.(lforce.ne.1)) then
      if (lforce.ne.1) then
        xkto=vsav(4,ml)
        tq=ee1exp(xkt)/max(1.e-24,ee1exp(xkto))
        ans1=rates(1,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans2=rates(2,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        ans3=rates(3,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans4=rates(4,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        if (lpri.gt.1) write (lun11,*)'type 49 scaling:',
     $    ml,vsav(2,ml),r19,tq,xnx,vsav(3,ml),xkt,xkto,tst,
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml),ett,t
        go to 9000
        endif
      vsav(2,ml)=r19
      vsav(1,ml)=bremsint(nb1)
      vsav(3,ml)=xnx
      vsav(4,ml)=xkt
      mlm=nilin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdti,np1r2,nidti,np1i2,nkdti,np1k2,mlm,
     $  nptrs,0,lun11)
      emax=etmpp(ntmp)*13.6+eth
      gglo=rlev(2,idest1)
      ggup=rlev(2,nlevp)
      idest3=idat1(np1i-1+nidti)
      idest4=idest3+1
      if (idest2.gt.nlevp) then
        jkk3=jkion+1
        if (lpri.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        ndtmp=npfi(13,jkk3)
        if (lpri.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        mllz=npar(ndtmp)
        iltmp=0
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))
     $      .and.(npar(ndtmp).eq.mllz))
           mlm=ndtmp-1
           call drd(ltyp2,lrtyp2,lcon2,
     $       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $       nptrs,0,lun11)
           iltmp=idat1(np1i2+nidt2-2)
           if (lpri.gt.1) then
             write (lun11,*)nidt2,iltmp,ndtmp
             write (lun11,*)np1r2,np1i2,np1k2,mlm
             call dprinto(ltyp2,lrtyp2,lcon2,
     $          nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,
     $          rdat1,idat1,kdat1,lun11)
             endif
           ndtmp=npnxt(ndtmp)
           enddo
         ggup=rdat1(np1r2+1)
         if (lpri.gt.1)
     $    write (lun11,*) ndtmp,iltmp,idest2,ggup
         endif
      if (lpri.gt.1) write (lun11,*)'before phint53'
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      if (lpri.gt.1) then
        write (lun11,*)'type 49 data:',idat1(np1i),
     $    idat1(np1i+nidt-1),t,xnx,
     $    eth,gglo,ggup,swrat
        call dprinto(ndesc,nrdesc,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,rdat1,idat1,kdat1,lun11)
        endif
      lprib=0
      if (lpri.gt.1) lprib=lpri
      rnist=rniss(idest1)
     $  *exp(-(ett+max(0.,13.605692*etmpp(1)))/(0.861707)/t)
     $  /rniss(nlevp)
      if (lpri.gt.1)
     $  write (lun11,*)'ett=',ett,etmpp(1),
     $  ett+max(0.,13.605692*etmpp(1)),
     $  rniss(idest1)/rniss(nlevp),rnist
      call phint53(stmpp,etmpp,ntmp,ett,ans1,ans2,ans3,ans4,
     $  abund1,abund2,ptmp1,ptmp2,xpx,opakab,rnist,
     $  opakc,rccemis,lprib,epi,ncn2,bremsa,t,swrat,xnx,
     $  lfast,lun11)
      if (lpri.gt.1) then
        npr=nb1
        write (lun11,*)'bautista threshold xsection:',
     $         npr,ett,eth,rdat1(np1r),sg(npr),ans2,swrat
        endif
      rates(1,ml)=ans1
      rates(2,ml)=ans2
      rates(3,ml)=ans3
      rates(4,ml)=ans4
      if (lpri.gt.1) write (lun11,*)'rates:',
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml)
      lpri=lprisv
      go to 9000
c
 50   continue
c     op line rad. rates
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
c     nb check this out:  no bound-bound decays from continuum
      if ((idest1.le.0).or.(idest1.ge.nlev)
     $  .or.(idest2.le.0).or.(idest2.ge.nlev))
     $      go to 9000
      aij=rdat1(np1r+2)
c      aij=min(aij,1.e+10)
      eeup=rlev(1,idest1)
      eelo=rlev(1,idest2)
      if (eeup.lt.eelo) then
         itmp=idest1
         idest1=idest2
         idest2=itmp
         endif
      elin=abs(rdat1(np1r))
      if (elin.le.1.e-34) go to 9000
      ggup=rlev(2,idest1)
      gglo=rlev(2,idest2)
      if (ml.le.0) go to 9000
      nilin=npar(ml)
      if (nilin.le.0) go to 9000
      nelin=npar(nilin)
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000
      flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo)
c
      mlm=nelin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      a=rdat1(np1r2+1)
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5)
      ener=12398.41/elin
      dele=ener*vtherm/3.e+10
      elammu=elin*1.e-4
      ans1=aij*(ptmp1+ptmp2)
      sigma=(0.02655)*flin*elin*(1.e-8)/vtherm
      sigvtherm=sigma
      jkkl=nplini(ml)
      if (jkkl.le.0) go to 9000
      ml3=nplin(jkkl)
      if (ml3.le.0) go to 9000
      mlm=ml3-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $  nptrs,0,lun11)
      elin=abs(rdat1(np1r))
      ener=12398.41/abs(elin)
      nb1=nbinc(ener,epi,ncn2)
      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10
     $      *flinabs(ptmp1)
c
c     turning off rex
      ans2=0.
c
      if (elin.gt.0.99e+9) then
         ans2=0.
         sigvtherm=0.
         endif
c      ans1=ans1+ans2*ggup/(1.e-36+gglo)
c     note that now opakab does not have abundance in
      opakab=sigvtherm
      lfasto=2
c      lfasto=4
      delea=0.
      lfnd=0
      lpriu=0
c      if (lpri.ge.1) lpriu=3
      call deleafnd(jkion,idest1,ml,
     $   nrdt,np1r,nidt,np1i,nkdt,np1k,
     $   idat1,rdat1,kdat1,nptrs,np2,
     $   npar,npnxt,npfi,npfirst,
     $   nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $   npconi2,ncsvn,delea,lfnd,lpriu,lun11)
c
      if (lfnd.eq.0) delea=rdat1(np1r+2)*(4.136e-15)
      ans4=ans1*ener*ergsev
      ans3=ans2*ener*ergsev
      rcem1=abund2*ans4*ptmp1/(1.e-34+ptmp1+ptmp2)
      rcem2=abund2*ans4*ptmp2/(1.e-34+ptmp1+ptmp2)
      opakb1=sigvtherm*abund1
c     this test should prevent calculation when called from func2
c     since abund1 will be zero
c      lpriu=lpri
      lpriu=0
c     if ((nrdesc.ne.9).and.(lfasto.le.4).and.(opakb1.gt.1.e-34))
c    $ call linopac(lpriu,lun11,opakb1,ans2,sigvtherm,vtherm,bremsa,
c    $               rcem1,rcem2,elin,vturb,t,a,delea,epi,ncn2,
c    $               opakc,opakscatt,rccemis,fline,lfasto)
      if (lpri.gt.1)
     $  write (lun11,*)'in ucalc, ind=50:',
     $  ml,nilin,nelin,elin,flin,rdat1(np1r+2),gglo,ggup,a,vtherm,vturb,
     $  ans1,ans2,idest1,idest2,idest3,idest4,nlev,sigvtherm,
     $  bremsa(nb1),nb1,abund1,abund2,delea,lfnd
      if (nrdesc.ne.9) go to 9000
c
c       special for 2 photon
        ansar2=0.
        em2ph=aij
        lskp=1
        emax=ener
        nbmx=nbinc(emax,epi,ncn2)
        if (lpri.gt.1)
     $  write (lun11,*)'in ucalc, ind=50:',
     $  ml,nilin,nelin,elin,flin,rdat1(np1r),gglo,ggup,a,vtherm,ans2
     $   ,nbmx
        rcemsum=0.
        lfastl=0
        ll=2
        do while (ll.le.nbmx)
          ansar2o=ansar2
          ansar2=epi(ll)*epi(ll)*max(0.,(epi(nbmx)-epi(ll)))
          rcemsum=rcemsum+(ansar2+ansar2o)
     $                   *(epi(ll)-epi(ll-lskp))/2.
          call enxt(epi(1),nb1,0,epi,ncn2,t,lfastl,lun11,
     $                  ll,lskp,nphint,lrcalc)
          ll=ll+lskp
          enddo
        rctmp1=0.
        rctmp2=0.
        ll=2
        do while (ll.le.nbmx)
          ansar2=epi(ll)*epi(ll)*max(0.,(epi(nbmx)-epi(ll)))
          ansar2=ansar2*em2ph*emax/(1.e-24+rcemsum)
          rctmp1=abund2*ansar2*ptmp1/12.56
          rctmp2=abund2*ansar2*ptmp2/12.56
          rccemis(1,ll)=rccemis(1,ll)+rctmp1
          rccemis(2,ll)=rccemis(2,ll)+rctmp2
          call enxt(epi(1),nb1,0,epi,ncn2,t,lfastl,lun11,
     $                  ll,nskp,nphint,lrcalc)
          ll=ll+nskp
          enddo
      go to 9000
c
 51   continue
c     line rates, col, burgess and tully from manuel
      idest1=idat1(np1i+2)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      eeup=rlev(1,idest1)
      eelo=rlev(1,idest2)
      if (eeup.lt.eelo) then
         itmp=idest1
         idest1=idest2
         idest2=itmp
         eeup=rlev(1,idest1)
         eelo=rlev(1,idest2)
         endif
      ggup=rlev(2,idest1)
      gglo=rlev(2,idest2)
      eijry=rdat1(np1r)
      eij=eijry*13.605692
      elin=12398.41/eij
      hij=elin*1.e-8
c      if (lpri.ne.0)
c     $ write (lun11,*)'type 51 data:',elin
      if (elin.le.1.e-24) go to 9000
      ekt=0.861707*t
      delt=12398.41/elin/ekt
      if (lpri.gt.1)
     $ write (lun11,*)elin,ekt,delt
c      if (delt.gt.50.) go to 9000
      c=rdat1(np1r+1)
      p1=rdat1(np1r+2)
      p2=rdat1(np1r-1+4)
      p3=rdat1(np1r-1+5)
      p4=rdat1(np1r-1+6)
      p5=rdat1(np1r-1+7)
      tk=t*1.e+4
c      tk=max(tk,(1.e+4)*12398.54/elin/(0.861707)/50.)
      tk=max(tk,2.8777e+6/elin)
      ik=idat1(np1i)
      cijpp=upsil(ik,eijry,c,p1,p2,p3,p4,p5,tk)
      ekt=0.861707*t
      delt=12398.41/elin/ekt
      cji=(8.626e-8)*cijpp/tsq
     $      /ggup
      exptmp=expo(-delt)
      cij=cji*ggup*exptmp/gglo
      if (lpri.gt.1)
     $ write (lun11,*)'ltyp=51',c,p1,p2,p3,p4,p5,ik,
     $      eij,idest1,idest2,cij,cji,xnx,cijpp
      ans1=cij*xnx
      ans2=cji*xnx
      elin=0.
      go to 9000
c
 52   continue
c     same as 59 but rate type 7
      go to 59
c
 53   continue
 533   continue
c     op pi xsections
       lprisv=lpri
c      if (lpri.ge.1) lpri=2
c     these are the initial and final levels and indeces
c     notice that these are relative to the current ion
c     (not relative to the element as a whole)
      idest1=idat1(np1i+nidt-2)
      idest2=nlevp+idat1(np1i-1+nidt-3)-1
      if (lpri.gt.1) write (lun11,*)'idest1=',idest1,idest2,nlevp,ml
      if ((idest1.ge.nlevp).or.(idest1.le.0)) go to 9000
      if (ml.le.0) go to 9000
      eth=rlev(4,idest1)-rlev(1,idest1)
      eexc=rlev(1,idest1)
      ett=eth
      nilin=npar(ml)
      if (lpri.gt.1) write (lun11,*)'nilin=',nilin,ml
      if (nilin.le.0) go to 9000
      ntmp=nrdt/2
c      ett=ett+max(0.,13.605692*etmpp(1))
      if (lpri.gt.1) write (lun11,*)'ett=',ett,etmpp(1)
      if (ett.le.0.) go to 9000
      nb1=nbinc(ett,epi,ncn2)
      tst=abs(bremsint(nb1)/max(1.e-24,vsav(1,ml))-1.)
      xkt=ett/(0.861707*t)
      r19=rr/1.e+19
c      if ((tst.le.0.01).and.(lforce.ne.1)) then
      if (lforce.ne.1) then
        xkto=vsav(4,ml)
        tq=ee1exp(xkt)/max(1.e-24,ee1exp(xkto))
        ans1=rates(1,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans2=rates(2,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        ans3=rates(3,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans4=rates(4,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        if (lpri.gt.1) write (lun11,*)'type 53 scaling:',
     $    ml,vsav(2,ml),r19,tq,xnx,vsav(3,ml),xkt,xkto,tst,
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml),ett,t
        go to 9000
        endif
      vsav(2,ml)=r19
      vsav(1,ml)=bremsint(nb1)
      vsav(3,ml)=xnx
      vsav(4,ml)=xkt
      mlm=nilin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      emax=etmpp(ntmp)*13.6+eth
      gglo=rlev(2,idest1)
      ggup=rlev(2,nlevp)
      idest3=idat1(np1i+nidt-1)
      idest4=idest3+1
      if (idest2.gt.nlevp) then
        jkk3=jkion+1
        if (lpri.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        ndtmp=npfi(13,jkk3)
        if (lpri.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        if (ndtmp.le.0) go to 9000
        mllz=npar(ndtmp)
        iltmp=0
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))
     $      .and.(npar(ndtmp).eq.mllz))
           mlm=ndtmp-1
           call drd(ltyp2,lrtyp2,lcon2,
     $       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $       nptrs,0,lun11)
           iltmp=idat1(np1i2+nidt2-2)
           if (lpri.gt.1) write (lun11,*)nidt2,iltmp,ndtmp
           ndtmp=npnxt(ndtmp)
           if (ndtmp.le.0) go to 9000
           enddo
c        NB fix to excited level PI and rec
         ett=ett+rdat1(np1r2)
         eth=ett
         ggup=rdat1(np1r2+1)
         if (lpri.gt.1)
     $    write (lun11,*) ndtmp,iltmp,idest2,ggup,ett
         endif
      sscal=1.
      do ml2=1,ntmp
        etmpp(ml2)=rdat1(np1r-1+2*ml2-1)
        stmpp(ml2)=rdat1(np1r-1+2*ml2)*1.e-18*sscal
        stmpp(ml2)=max(stmpp(ml2),0.)
        stmpe(ml2)=stmpp(ml2)*(etmpp(ml2)*13.605692+ett)
        if (lpri.gt.1) write (lun11,9819)ml2,etmpp(ml2),stmpp(ml2)
 9819   format (1x,i6,2(1pe11.3))
        enddo
      optst=abund1*stmpp(1)
c      if ((optst.lt.opcrit).and.(lfast.eq.2)) go to 9000
      ntmp2=nptmpdim
c     nb includes extrapolation
c     this is dangerous.  It does the right thing for ground-ground,
c       but some cross sections should not be extrapolated.
c      call phextrap(etmpp,stmpp,ntmp,ntmp2,ett,ncn2,lpri,lun11)
      if (lpri.gt.1) write (lun11,*)'before phint53',eexc,eth,lfast
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      if (lpri.gt.1) then
        write (lun11,*)'type 53 data:',idat1(np1i),
     $    idat1(np1i+nidt-1),t,xnx,
     $    eth,gglo,ggup,swrat
        call dprinto(ndesc,nrdesc,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,rdat1,idat1,kdat1,lun11)
        endif
      lprib=0
      if (lpri.gt.1) lprib=lpri
      rnist=rniss(idest1)
     $  *exp(-(ett+max(0.,13.605692*etmpp(1)))/(0.861707)/t)
     $  /rniss(nlevp)
      if (lpri.gt.1)
     $  write (lun11,*)'ett=',ett,etmpp(1),
     $  ett+max(0.,13.605692*etmpp(1)),
     $  rniss(idest1)/rniss(nlevp),rnist
      call phint53(stmpp,etmpp,ntmp,ett,ans1,ans2,ans3,ans4,
     $  abund1,abund2,ptmp1,ptmp2,xpx,opakab,rnist,
     $  opakc,rccemis,lprib,epi,ncn2,bremsa,t,swrat,xnx,
     $  lfast,lun11)
      if (lpri.gt.1) then
        npr=nb1
        write (lun11,*)'bautista threshold xsection:',
     $         npr,ett,eth,rdat1(np1r),sg(npr),ans2,swrat
        endif
      if (lpri.gt.1) then
        temp=t*1.e+4
        do ml2=1,ntmp
          etmpp(ml2)=rdat1(np1r-1+2*ml2-1)
          stmpp(ml2)=rdat1(np1r-1+2*ml2)
          enddo
        lprim=0
        call milne(temp,ntmp,etmpp,stmpp,ett/13.6,alphamilne,
     $     lun11,lprim)
        alphamilne=alphamilne*xnx
        amilnerr=(log10(alphamilne/max(1.e-34,ans2)))
        if ((abs(amilnerr).gt.0.05)
     $    .and.((alphamilne.gt.1.e-28).or.(ans2.gt.1.e-28))
     $    .and.(lfast.gt.1))
     $     write (lun11,*)'milne error',alphamilne,ans2,amilnerr
        endif
      rates(1,ml)=ans1
      rates(2,ml)=ans2
c      rates(2,ml)=alphamilne
      rates(3,ml)=ans3
      rates(4,ml)=ans4
      if (lpri.gt.1) write (lun11,*)'rates:',
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml)
      lpri=lprisv
      go to 9000
c
 54   continue
c     h-like cij, bautista (hlike ion)
      idest1=idat1(np1i-1+nidt-3)
      idest2=idat1(np1i+nidt-3)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      ans3=0.
      ans4=0.
      lprisv=lpri
c      if (lpri.ge.1) lpri=2
      if (lpri.gt.1) write (lun11,*)'type 54 data:',
     $  idat1(np1i-1+nidt-3),idat1(np1i+nidt-3)
      if (rlev(1,idest2).lt.rlev(1,idest1)) then
        itmp=idest2
        idest2=idest1
        idest1=itmp
        endif
      eeup=rlev(1,idest2)
      eelo=rlev(1,idest1)
      elin=12398.41/abs(eeup-eelo+1.e-24)
      hij=elin*1.e-8
      ekt=0.861707*t
      delt=12398.41/elin/ekt
c      if (delt.gt.50.) go to 9000
      ni=ilev(1,idest2)
      li=ilev(3,idest2)
      nf=ilev(1,idest1)
      lf=ilev(3,idest1)
      if (lpri.gt.1) write (lun11,*)
     $  eeup,eelo,elin,ni,li,nf,lf
      if (ni.eq.nf) go to 9000
      if (ni.lt.nf) then
        ntmp=ni
        ni=nf
        nf=ntmp
        endif
      iq=idat1(np1i+nidt-2)
      if (lpri.gt.1)
     $ write (lun11,*)'before anl1:',ni,nf,li,lf,iq,idest1,idest2,
     $  eelo,eeup,idat1(np1i-1+nidt-3),idat1(np1i+nidt-3)
      call anl1(ni,nf,lf,iq,alm,alp,lpri,lun11)
      ans1=alp
      if (li.lt.lf) ans1=alm
      lpri=lprisv
      go to 9000
c
 55   continue
c      hydrogenic pi xsections, bautista format
      lprisv=lpri
      if (lpri.gt.1)
     $  write (lun11,*)'in ucalc, ind=55:',(idat1(np1i-1+mm),mm=1,5)
      idest1=idat1(np1i+nidt-2)
      ett=rlev(1,nlevp)-rlev(1,idest1)
      idest2=nlevp
      if (ett.le.1.e-5) go to 9000
      eth=ett
      nb1=nbinc(eth,epi,ncn2)
      gglo=rlev(2,idest1)
      ggup=rlev(2,nlevp)
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      ekt=t*(0.861707)
      if (ml.le.0) go to 9000
      nilin=npar(ml)
      if (nilin.le.0) go to 9000
      nelin=npar(nilin)
      if (lpri.gt.1)
     $ write (lun11,*)'in ucalc, ind=55:',
     $   ml,nilin,nelin
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000
      mlm=nilin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      nistage=idat1(np1i2)
      mlm=nelin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      nzel=idat1(np1i2)
      zz=float(nzel-nistage)
      sgth=(6.3e-18)/zz/zz
      ll=nb1
      lfastl=1
      do while (ll.le.nphint)
        epii=epi(ll)
        sg(ll)=sgth*(epii/ett)**(-3)
        call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,
     $                  ll,nskp,nphint,lrcalc)
        ll=ll+nskp
        enddo
      lprib=0
      if (lpri.gt.1) lprib=lpri
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      lpri=lprisv
      go to 9000
c
 56   continue
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      if (rlev(1,idat1(np1i+1)).lt.rlev(1,idat1(np1i))) then
        idest2=idat1(np1i)
        idest1=idat1(np1i+1)
        endif
      lprisv=lpri
c      if (lpri.ge.1) lpri=2
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      eeup=rlev(1,idest2)
      eelo=rlev(1,idest1)
      dele=abs(eeup-eelo)
      if (dele.le.1.e-16) go to 9000
      ntmp=nrdt/2
      do kl=1,ntmp
        ttmp(kl)=rdat1(np1r-1+kl)
        enddo
      tfnd=log10(t*1.e+4)
      jlo=0
      call hunt3(ttmp,ntmp,tfnd,jlo,0,lun11)
      jlo=min(jlo,ntmp-1)
      nind=ntmp+jlo
      if (lpri.gt.1) write (lun11,*)'type 56:',
     $  idest1,idest2,ggup,gglo,dele,jlo,nind,
     $  rdat1(np1r-1+nind),rdat1(np1r-1+nind+1),
     $  tfnd,ttmp(jlo+1),ttmp(jlo)
      cijpp=(rdat1(np1r-1+nind+1)-max(1.e-36,rdat1(np1r-1+nind)))
     $   *(tfnd-ttmp(jlo))/(ttmp(jlo+1)-ttmp(jlo)+1.e-24)
     $     +max(1.e-36,rdat1(np1r-1+nind))
      cijpp=max(0.,cijpp)
      ekt=0.861707*t
      delt=dele/ekt
      cij=0.
      exptmp=expo(-delt)
      if (lpri.gt.1) write (lun11,*)'type 56:',
     $  idest1,idest2,ggup,gglo,dele
      cij=(8.626e-8)*cijpp*exptmp/tsq/gglo
      cji=(8.626e-8)*cijpp/tsq/ggup
      ans1=cij*xnx
      ans2=cji*xnx
      if (lpri.gt.1) write (lun11,*)'type 56 data:',
     $  idest1,idest2,dele,cijpp,delt,exptmp,cij,cji,
     $  nind,gglo,ggup,tfnd,ntmp,jlo,ntmp,ttmp(jlo)
      lpri=lprisv
      go to 9000
c
 57   continue
c     same as  65 (?)
c     effective charge to be used in coll. ion.
      lprisv=lpri
      lpri=0
c      if (lprisv.ge.1) lpri=2
      tz=t*1.e+4
      idest1=idat1(np1i+nidt-2)
      idest2=nlevp
      if (lpri.gt.1)
     $ write (lun11,*)'in ucalc at 57:',idest1,idat1(np1i),rdat1(np1r)
      if ((idat1(np1i).le.0).or.(idest1.le.1).or.(idest1.gt.nlevp))
     $        go to 9000
      i57=idat1(np1i)
      eth=max(0.,rlev(1,nlevp)-rlev(1,idest1))
      ekt=0.861707*t
c      tz=max(tz,(1.e+4)*eth/(0.861707)/50.)
c      tz=max(tz,2.320975e+02*eth)
      e1=rlev(1,idest1)
      ep=rlev(4,idest1)
      if (ep.le.0.) go to 9000
      call calt57(tz,xnx,e1,ep,i57,cion,crec,lun11,lpri)
      if (lpri.gt.1)
     $ write (lun11,*)'ltype=57:',cion,crec,gglo,ggup,nlevp,idest1,rinf,
     $  eth,ekt,ans1,ans2
c
c     trying a fudge to test cloudy's ci
c      if (lprisv.ge.1) write (lun11,*)'fudging ci for test'
c      if (i57.eq.2) cion=cion*2.
c      if (i57.eq.2) crec=crec*2.
c      if (i57.eq.3) cion=cion*5.
c      if (i57.eq.3) crec=crec*5.
c      if (i57.eq.4) cion=cion*15.
c      if (i57.eq.4) crec=crec*15.
c
      ans1=cion*xnx
      ggup=rlev(2,nlevp)
      gglo=rlev(2,idest1)
c     note that rinf has exponential removed
      rinf=gglo/(1.e-36+ggup)
      ans2=crec*rinf*xnx*xnx
c     set to zero for ground state because we have more accurate rates
c     for these levels: types 95 or 25
      if (idest1.eq.1) then
        ans1=0.
        ans2=0.
        endif
      go to 9000
c
 58   continue
c      bautista cascade rates. defunct.
      go to 9000
c
 59   continue
      lprisv=lpri
      lpril=lpri
c      if (lpri.ge.1) lpril=2
      if (lpril.gt.1) write (lun11,*)'ltyp=59',ml,npar(ml)
      if (lpril.gt.1) write (lun11,*)(rdat1(np1r-1+jj),jj=1,nrdt)
      if (lpril.gt.1) write (lun11,*)(idat1(np1i-1+jj),jj=1,nidt),nidt
      if (lpril.gt.1) write (lun11,*)(kdat1(np1k-1+jj),jj=1,nkdt)
      if (ml.le.0) go to 9000
c
c
c     experiment with only vfky
c      if (nrdt.le.6) go to 9000
c
      lfastl=1
      nilin=npar(ml)
      idest3=idat1(np1i+nidt-1)
      idest4=idat1(np1i+nidt-3)
c     why was this statement here?
      if (idest4.gt.idest3+1) go to 9000
      idest1=idat1(np1i+nidt-2)
      idest2=nlevp+idat1(np1i-1+nidt-3)-1
      idest2=max(idest2,1)
c      nb must uncomment these if func2a is called
c      if (nrdesc.eq.7) then
c        idest2=nlevp
c        endif
      if ((nilin.le.0).or.(nilin.gt.np2)) go to 9000
      mlm=nilin-1
      call drd(ltyp2,lrtyp2,lcon2,
     $  nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $  nptrs,0,lun11)
      if (lpril.gt.1)
     $ write (lun11,*)ml,nilin,rdat1(np1r),idest1,rdat1(np1r2),nlevp
      ett=rdat1(np1r2)
      if ((idest1.gt.nlevp).or.(idest1.le.0)) go to 9000
      if (ml.le.0) go to 9000
      if (ett.le.0.) go to 9000
      nb1=nbinc(ett,epi,ncn2)
      numcon2=max(2,ncn2/50)
      nphint=ncn2-numcon2
      nphint=max(nphint,nb1+1)
      if (nb1.ge.nphint-1) go to 9000
      ett=rdat1(np1r)
      nb1=nbinc(ett,epi,ncn2)
      if (nb1.ge.(ncn2-1)) go to 9000
      tst=abs(bremsint(nb1)/max(1.e-24,vsav(1,ml))-1.)
      if (lpril.gt.1)
     $ write (lun11,*)ett,nb1,bremsint(nb1),ml,vsav(1,ml),tst,
     $  lforce
      xkt=ett/(0.861707*t)
      r19=rr/1.e+19
      if (lforce.ne.1) then
        xkto=vsav(4,ml)
        tq=ee1exp(xkt)/max(1.e-24,ee1exp(xkto))
        ans1=rates(1,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans2=rates(2,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        ans3=rates(3,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans4=rates(4,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
c        write (lun11,*)'in ucalc2, scaling used'
        go to 9000
        endif
      vsav(2,ml)=r19
      vsav(1,ml)=bremsint(nb1)
      vsav(3,ml)=xnx
      vsav(4,ml)=xkt
      gglo=rlev(2,1)
      ggup=rlev(2,nlevp)
      if (idest2.gt.nlevp) then
        jkk3=jkion+1
        if (lpril.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        ndtmp=npfi(13,jkk3)
        if (lpril.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        if (ndtmp.le.0) go to 9000
        mllz=npar(ndtmp)
        iltmp=0
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))
     $      .and.(npar(ndtmp).eq.mllz))
           mlm=ndtmp-1
           call drd(ltyp2,lrtyp2,lcon2,
     $       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $       nptrs,0,lun11)
           iltmp=idat1(np1i2+nidt2-2)
           if (lpril.gt.1) write (lun11,*)nidt2,iltmp,ndtmp
           ndtmp=npnxt(ndtmp)
           if (ndtmp.le.0) go to 9000
           enddo
c        NB fix to excited level PI and rec
         ett=ett+rdat1(np1r2)
         eth=ett
         ggup=rdat1(np1r2+1)
         if (lpril.gt.1)
     $    write (lun11,*) ndtmp,iltmp,idest2,ggup,ett
         endif
      if (lpril.gt.1) write (lun11,*)nlevp,ggup
      if (ggup.le.1.e-24) then
        if (lpril.gt.1) write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      ett=rdat1(np1r)
      nb1=nbinc(ett,epi,ncn2)
      if (lpril.gt.1)
     $ write (lun11,*)'ett=',ett,nb1,nphint,swrat,gglo,ggup
      if (nb1.ge.(nphint-1)) go to 9000
      if ((bremsint(nb1).lt.1.e-20).and.(lpril.gt.1))
     $    write (lun11,*)'skipping 59',
     $         nb1,bremsint(nb1)
      if (bremsint(nb1).lt.1.e-20) go to 9000
      if (nrdt.eq.9) then
          ett=rdat1(np1r)
          emax=rdat1(np1r+1)
          e0=rdat1(np1r+2)
          s0=rdat1(np1r-1+4)
          ya=rdat1(np1r-1+5)
          pp=rdat1(np1r-1+6)
          yw=rdat1(np1r-1+7)
          y0=rdat1(np1r-1+8)
          y1=rdat1(np1r-1+9)
          l2=0
        else
          e0=rdat1(np1r+1)
          s0=rdat1(np1r+2)
          ya=rdat1(np1r-1+4)
          pp=rdat1(np1r-1+5)
          yw=rdat1(np1r-1+6)
          y0=0.
          y1=0.
          l2=idat1(np1i+2)
        endif
      ywsq=yw*yw
      qq=5.5+l2-pp/2.
      if (lpril.gt.1) write (lun11,*)'qq=',
     $   l2,qq,ya,ywsq,pp,yw,s0
      ll=nb1
      do while (ll.le.nphint)
        epii=epi(ll)
        xx=epii/e0-y0
        if (nrdt.eq.9) then
            yy=sqrt(xx*xx+y1*y1)
          else
            yy=xx
          endif
        yyqq=qq*log(max(1.e-34,yy))
        yyqq=exp(-min(60.,max(-60.,yyqq)))
        term1=((xx-1.)*(xx-1.)+ywsq)
        term2=yyqq
        term3=(1.+sqrt(yy/ya))**(-pp)
        ff=term1*term2*term3
        sg(ll)=s0*ff*(1.e-18)
        if (lpril.gt.1) write (lun11,*)ll,epii,sg(ll),
     $    yy,yyqq,xx,term1,term2,term3,qq,ff
        call enxt(ett,nb1,0,epi,ncn2,t,lfastl,lun11,
     $                  ll,nskp,nphint,lrcalc)
        ll=ll+nskp
        enddo
      ekt=t*(0.861707)
      lprib=0
      if (lpril.gt.1) lprib=lpril
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      if (lpril.gt.1) then
        npr=nb1
        write (lun11,*)'verner threshold xsection:',
     $         npr,ett,sg(npr),opakab
        endif
c     nb this turns off all recombination into excited levels
c     for type 59...
c      if ((nrdesc.eq.1).or.(idest1.gt.1)) then
c        ans4=0.
c        ans2=0.
c        endif
      rates(1,ml)=ans1
      rates(2,ml)=ans2
      rates(3,ml)=ans3
      rates(4,ml)=ans4
      lpri=lprisv
      go to 9000
c
 60   continue
c      go to 9000
c     calloway h-like coll. strength
      lpril=0
c      if (lpri.ge.1) lpril=2
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      if (rlev(1,idat1(np1i+1)).lt.rlev(1,idat1(np1i))) then
        idest2=idat1(np1i)
        idest1=idat1(np1i+1)
        endif
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      dele=abs(rlev(1,idest2)-rlev(1,idest1))
      if (dele.le.1.e-24) go to 9000
      ekt=0.861707*t
      delt=dele/ekt
      temp=t*1.e+4
      temp=max(temp,0.02*dele*1.e+4/(0.861707))
      call calt60_62(temp,nrdt,ndesc,np1r,np1i,rdat1,idat1,cijpp)
c      cijpp=cijpp/2./2.
      cji=(8.626e-8)*cijpp/tsq/(1.e-16+ggup)
      exptmp=expo(-delt)
      cij=cji*ggup*exptmp/(1.e-16+gglo)
      ans1=cij*xnx
      ans2=cji*xnx
      if (lpril.gt.1) then
        write (lun11,*)'ltyp=60',idest1,idest2,temp,flin,ggup,gglo
        write (lun11,*)'       ',nrdt,(rdat1(np1r-1+mm),mm=1,8),jlo
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp,dele,delt
        endif
      go to 9000
c
 61   continue
      go to 9000
c
 62   continue
      go to 60
c
 63   continue
c      if (lpri.ne.0) write (lun11,*) 'type 63 data not implemented'
c      go to 9000
      lpril=0
c      if (lpri.ge.1) lpril=2
      idest1=idat1(np1i-1+nidt-3)
      idest2=idat1(np1i+nidt-3)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      eeup=rlev(1,idest2)
      eelo=rlev(1,idest1)
      elin=12398.41/abs(eeup-eelo+1.e-24)
      hij=elin*1.e-8
      ekt=0.861707*t
      delt=12398.41/elin/ekt
      if (lpril.ne.0) write (lun11,*)'delt=',delt
      if (delt.gt.50.) go to 9000
      ni=ilev(1,idest1)
      li=ilev(3,idest1)
      nf=ilev(1,idest2)
      lf=ilev(3,idest2)
      sum=0.
      iq=idat1(np1i+nidt-2)
      if (lpril.ne.0)
     $ write (lun11,*)'ltyp=63',idest1,idest2,ni,li,nf,lf
      if (nf.eq.ni) then
        if (lpril.ne.0)
     $     write (lun11,*)'nf=ni',lf,li
        if (abs(lf-li).eq.1) then
          lff=min(lf,li)
          lii=max(lf,li)
          li1=max(1,lff)     ! mab
          do  nn=li1,ni-1
            if (lpril.ne.0)
     $        write (lun11,*)'before anl1'
           if (lii.ge.1) then
            if (lpril.ne.0)
     $        write (lun11,*)'li=1',ni,nn,lii-1,iq
            call anl1(ni,nn,lii-1,iq,alm,alp,lpril,lun11)
c              write (lun11,*)'li=1',ni,nn,lii-1,iq,alp
            sum=sum+alp
           endif
           if (nn.gt.lii+1) then
            if (lpril.ne.0)
     $        write (lun11,*)'nn=li+1',ni,nn,lii+1,iq
            call anl1(ni,nn,lii+1,iq,alm,alp,lpril,lun11)
            sum=sum+alm
           endif
          enddo
          if (lpril.ne.0)
     $     write (lun11,*)'after anl1',sum
          ecm=abs(rlev(1,idest1)-rlev(1,idest2))*8059.9
          ecm=0.
          nnz=idat1(np1i-1+4)
          tbig=t*1.e+4
          z1=1.
          rm=1800.
          il=0
          psi=0.75/nnz/nnz*lii/(2*lii+1)*ni*ni*(ni*ni-lii*lii)
          if (lpril.ne.0)
     $     write (lun11,*)'before amcrs',ecm,ni,lii,sum
          call amcrs(ni,lii,tbig,nnz,z1,rm,xnx,sum,ecm,psi,il,cn,
     $        lpril,lun11)
          cno=cn
          iz=idat1(np1i-1+4)
          if (lf.lt.li) then
            ans1=cn
            ans2=cn*rlev(2,idest1)/rlev(2,idest2)
          else
            ans2=cn
            ans1=cn*rlev(2,idest2)/rlev(2,idest1)
          endif
          if (lpril.ne.0)
     $     write (lun11,*)'after amcrs',cn,iz,cno,ans1,ans2
        endif
      else
        if (lpril.ne.0)
     $     write (lun11,*)'nf.ne.ni'
        aa1=0.
        if (abs(lf-li).eq.1) then
          sum=0.
          nu=max(ni,nf)
          nll=min(ni,nf)
          do lff=0,nll-1
            call anl1(nu,nll,lff,iq,alm,alp,lpri,lun11)
             sum=sum+alp*(2*lff+3)
             if (lff.gt.0)  then
              sum=sum+alm*(2*lff-1)
             endif
             if (lff.eq.lf .and. li.gt.lf) aa1=alp
             if (lff.eq.lf .and. li.lt.lf) aa1=alm
             if (lpril.ne.0) write (lun11,*)'after anl1',
     $           lff,li,lf,sum,alp,alm,aa1
          enddo
          if (lpril.ne.0)
     $     write (lun11,*)'after anl1',sum,alp,alm,aa1
          nnz=idat1(np1i-1+4)
          tbig=t*1.e+4
          call erc(nll,nu,tbig,nnz,se,sd,sum,lun11,lpril)
c ***** check if ans1 and ans2 are correct or inverted
          ans1=se*(2*lf+1)*aa1/sum
          ans2=sd*(2*li+1)*aa1/sum
          if ((nf.gt.ni).or.(lf.gt.li)) then
           atmp=ans1
           ans1=ans2
           ans2=atmp
          endif
          if (lpril.ne.0)
     $     write (lun11,*)'after erc',se,sd,ans1,ans2
        endif
      endif
c
      ans1=ans1*xnx
      ans2=ans2*xnx
      go to 9000
c
 64   continue
c     hydrogenic pi xsections, bautista format
      lprisv=lpri
      idest1=idat1(np1i+nidt-2)
      ett=abs(rlev(1,nlevp)-rlev(1,idest1))
      if (lpri.gt.1)
     $  write (lun11,*)'in ucalc, ind=64:',(rdat1(np1r-1+mm),mm=1,5)
      if (ett.le.1.e-5) go to 9000
      zzz=float(idat1(np1i+2))
      enn=float(idat1(np1i))
      eth=ett
      nb1=nbinc(eth,epi,ncn2)
      gglo=rlev(2,idest1)
      swrat=gglo
      idest2=nlevp
      ekt=t*(0.861707)
      ll=nb1
      lorb=idat1(np1i+1)
      ic=idat1(np1i+2)
      nq=idat1(np1i)
      mm=0
      lfastl=1
      do while (ll.le.nphint)
        mm=mm+1
        epii=epi(ll)
        e=epii
        eth=ett
        erel=max(0.,(e-eth)/13.605692)
        call hphotx(erel,ic,nq,xsec,lun11,lpri)
        sg(ll)=xsec(lorb+1)*(1.e-18)
        stmpp(mm)=xsec(lorb+1)
        etmpp(mm)=erel
        call enxt(ett,nb1,lpri,epi,ncn2,t,lfastl,lun11,
     $                  ll,nskp,nphint,lrcalc)
        ll=ll+nskp
        enddo
      lprib=0
      if (lpri.gt.1) lprib=lpri
      call phintfo(sg,ett,ans1,ans2,ans3,ans4,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lprib,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      lprim=0
      ntmp=ll-nb1
      temp=t*1.e+4
      ntmp=mm
      call milne(temp,ntmp,etmpp,stmpp,eth/13.6,ans2,lun11,lprim)
      ans2=ans2*swrat
      lpri=lprisv
      go to 9000
c
c
 65   continue
c     effective charge to be used in coll. ion.
      tz=t*1.e+4
      idest1=idat1(np1i+nidt-2)
      idest2=nlevp
      ggup=rlev(2,nlevp)
      gglo=rlev(2,1)
      eth=max(0.,rlev(1,nlevp)-rlev(1,idest1))
      ekt=0.861707*t
c      if (eth/ekt.gt.50.) go to 9000
      call szirco(idat1(np1i),tz,rdat1(np1r),cii)
      ans1=cii*xnx
c     note that rinf has exponential removed
      rinf=(2.08e-22)*gglo/ggup/t/tsq
      ans2=ans1*rinf*expo(eth/ekt)
      go to 9000
c
 66   continue
c     Like type 69 but, data in fines tructure
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      if (rlev(1,idat1(np1i+1)).lt.rlev(1,idat1(np1i))) then
        idest2=idat1(np1i)
        idest1=idat1(np1i+1)
        endif
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      elin=rdat1(np1r)
      if (elin.le.1.e-24) go to 9000
      elin=12398.41/elin
      ekt=0.861707*t
      delt=12398.41/elin/ekt
c      if (delt.gt.50.) go to 9000
      hij=elin*1.e-8
      temp=t*1.e+4
c      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)
      temp=max(temp,2.8777e+6/elin)
      call calt66(temp,np1r,rdat1,gamma)
      cijpp=gamma
      cji=(8.626e-8)*cijpp/tsq/ggup
        exptmp=expo(-delt)
        cij=cji*ggup*exptmp/gglo
      ans1=cij*xnx
      ans2=cji*xnx
      if (lpri.gt.1) then
        write (lun11,*)'ltyp=66',idest1,idest2,elin,flin,ggup,gglo
        write (lun11,*)'       ',nrdt,(rdat1(np1r-1+mm),mm=1,8),nind,jlo
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp
        endif
      elin=0.
      go to 9000
c
 67   continue
c     Effective collision strengths from Keenan et al.
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      if (rlev(1,idat1(np1i+1)).lt.rlev(1,idat1(np1i))) then
        idest2=idat1(np1i)
        idest1=idat1(np1i+1)
        endif
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      elin=abs(rdat1(np1r))
      hij=elin*1.e-8
      if (elin.le.1.e-24) go to 9000
      ekt=0.861707*t
      delt=12398.41/elin/ekt
      temp=t*1.e+4
c      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)
      temp=max(temp,2.8777e+6/elin)
      call calt67(temp,np1r,rdat1,gamma)
      cijpp=gamma
      cijpp=max(0.,cijpp)
      cji=(8.626e-8)*cijpp/tsq/ggup
        exptmp=expo(-delt)
        cij=cji*ggup*exptmp/gglo
      ans1=cij*xnx
      ans2=cji*xnx
      if (lpri.gt.1) then
        write (lun11,*)'ltyp=69',idest1,idest2,elin,flin,ggup,gglo
        write (lun11,*)'       ',nrdt,(rdat1(np1r-1+mm),mm=1,8),nind,jlo
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp
        endif
      elin=0.
      go to 9000
c
 68   continue
c     coll. strength He-like ions by Zhang & Sampason
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      if (rlev(1,idat1(np1i+1)).lt.rlev(1,idat1(np1i))) then
        idest2=idat1(np1i)
        idest1=idat1(np1i+1)
        endif
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      eeup=rlev(1,idest2)
      eelo=rlev(1,idest1)
      elin=12398.41/abs(eeup-eelo+1.e-24)
      hij=elin*1.e-8
      if (elin.le.1.e-24) go to 9000
      ekt=0.861707*t
      delt=12398.41/elin/ekt
      temp=t*1.e+4
c      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)
      temp=max(temp,2.8777e+6/elin)
      if (lpri.gt.1) then
        write (lun11,*)'ltyp=68',idest1,idest2,elin,flin,ggup,gglo
        write (lun11,*)'       ',nrdt,(rdat1(np1r-1+mm),mm=1,8),nind,jlo
        endif
      call calt68(temp,np1r,np1i,rdat1,idat1,gamma)
      cijpp=gamma
      cijpp=max(cijpp,0.)
      cji=(8.626e-8)*cijpp/tsq/ggup
      ekt=0.861707*t
      delt=12398.41/elin/ekt
        exptmp=expo(-delt)
        cij=cji*ggup*exptmp/gglo
      ans1=cij*xnx
      ans2=cji*xnx
      if (lpri.gt.1) then
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp
        endif
      elin=0.
      go to 9000
c
 69   continue
c     Kato & Nakazaki (1996) fit to Helike coll. strgt
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      if (rlev(1,idat1(np1i+1)).lt.rlev(1,idat1(np1i))) then
        idest2=idat1(np1i)
        idest1=idat1(np1i+1)
        endif
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      eeup=rlev(1,idest2)
      eelo=rlev(1,idest1)
      elin=12398.41/abs(eeup-eelo+1.e-24)
      hij=elin*1.e-8
      if (elin.le.1.e-24) go to 9000
      ekt=0.861707*t
      delt=12398.41/elin/ekt
      m=nrdt
      temp=t*1.e+4
c      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)
      temp=max(temp,2.8777e+6/elin)
      if (lpri.gt.1) then
        write (lun11,*)'ltyp=69',idest1,idest2,elin,flin,ggup,gglo
        write (lun11,*)'       ',nrdt,(rdat1(np1r-1+mm),mm=1,8),nind,jlo
        endif
      call calt69(temp,m,np1r,rdat1,gamma)
      cijpp=gamma
      cijpp=max(cijpp,0.)
      cji=(8.626e-8)*cijpp/tsq/ggup
      ekt=0.861707*t
      delt=12398.41/elin/ekt
        exptmp=expo(-delt)
        cij=cji*ggup*exptmp/gglo
      ans1=cij*xnx
      ans2=cji*xnx
      if (lpri.gt.1) then
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp
        endif
      elin=0.
      go to 9000
c
 70   continue
c     Coefficients for phot x-section of suplevels
c      lfastl=lfast
      lfastl=3
      temp=t*1.e+4
      ans3=0.
      ans4=0.
      den=xpx
      m=1000
      lpric=0
c      if (lpri.ge.1) lpric=2
      mlion=npar(ml)
      idest1=idat1(np1i+nidt-2)
      idest1=min(idest1,nlev-1)
      idest2=nlev+idat1(np1i-1+nidt-3)-1
      idest2=max(idest2,nlev)
      ggup=rlev(2,nlevp)
      ett=abs(rlev(1,idest1)-rlev(1,nlevp))
      if (lpric.ge.1)
     $ write (lun11,*)'rlev:',idest1,nlevp,rlev(1,idest1),rlev(1,nlevp)
      if (idest2.gt.nlevp) then
        jkk3=jkion+1
        if (lpric.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        ndtmp=npfi(13,jkk3)
        if (lpric.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        mllz=npar(ndtmp)
        iltmp=0
        nptmp=mllz
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))
     $      .and.(nptmp.eq.mllz))
           mlm=ndtmp-1
           call drd(ltyp2,lrtyp2,lcon2,
     $       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $       nptrs,0,lun11)
           iltmp=idat1(np1i2+nidt2-2)
           if (lpric.gt.1) write (lun11,*)nidt2,iltmp,ndtmp
           ndtmp=npnxt(ndtmp)
           nptmp=0
           if (ndtmp.ne.0) nptmp=npar(ndtmp)
           enddo
         ggup=rdat1(np1r2+1)
         ett=abs(rlev(1,idest1)+rdat1(np1r2))
         endif
       if (lpric.ge.1)
     $    write (lun11,*) ndtmp,iltmp,idest2,ggup,ett
      xkt=ett/(0.861707*t)
      nb1=nbinc(ett,epi,ncn2)
      if (lforce.ne.1) then
        xkto=vsav(4,ml)
        tq=ee1exp(xkt)/max(1.e-24,ee1exp(xkto))
        ans1=rates(1,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans2=rates(2,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        ans3=rates(3,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans4=rates(4,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        if (lpric.gt.1) write (lun11,*)'type 53 scaling:',
     $    ml,vsav(2,ml),r19,tq,xnx,vsav(3,ml),xkt,xkto,tst,
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml),ett,t
        go to 9000
        endif
      vsav(2,ml)=r19
      vsav(1,ml)=bremsint(nb1)
      vsav(3,ml)=xnx
      vsav(4,ml)=xkt
      mlm=mlion-1
      call drd(ltyp2,lrtyp2,lcon2,
     $  nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $  nptrs,0,lun11)
      ist=idat1(np1i2)
      ic=ist
      eth=ett
      gglo=rlev(2,idest1)
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      if (lpric.ne.0) then
        write (lun11,*)'type 70 data:',idat1(np1i),idat1(np1i+nidt-1),t,
     $           xnx,eth,gglo,ggup,swrat
        call dprinto(ndesc,nrdesc,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,rdat1,idat1,kdat1,lun11)
        endif
      ettry=ett/13.6
      call calt70(temp,den,ettry,ic,m,np1r,np1i,rdat1,idat1,
     1             ntmp,etmpp,stmpp,rec,al,lun11,lpric)
      if (lpric.ne.0) write (lun11,*)'after  calt70:',rec,stmpp(1)
      crit53=0.01
      do mm=1,ntmp
        stmpp(mm)=stmpp(mm)*1.e-18
        stmpp(mm)=max(stmpp(mm),0.)
        enddo
      call phint53hunt(stmpp,etmpp,ntmp,ett,ans1,ans2d,ans3d,ans4s,
     $ lpric,epi,ncn2,bremsa,t,swrat,xnx,crit53,lfastl,lun11)
      if (ans2d.le.1.e-36) then
        ans1=0.
        ans2=0.
        go to 9000
        endif
      scale=rec*xnx/ans2d
      ans1=ans1*scale
c     does the swrat not belong?
c      ans2=rec*xnx*swrat
      ans2=rec*xnx
c      ans2=ans2d
      tm=t*1.e4
      q2=2.07e-16*xnx*(tm**(-1.5))
      rs=q2/swrat
      ans1o=ans1
c      ans1=min(ans1,ans2/rs)
      if (lpric.ge.2)
     $ write (lun11,*)'type 70 limit:',ans2,rs,swrat,
     $   xnx,tm,q2,ans1o,ans1,scale,rec
c
c     testing superlevel phot.
c      ans1=0.
c
      rates(1,ml)=ans1
      rates(2,ml)=ans2
      rates(3,ml)=ans3
      rates(4,ml)=ans4
      if (lpric.ge.1) write (lun11,*)'rates:',
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml)
      go to 9000
c
 71   continue
c     Transition rates from superlevel to spect. lvls
      temp=t*1.e+4
      lpril=0
      den=xpx
      m=1000
      if (lpril.ne.0)
     $  write (lun11,*)'before calt71:',rdat1(np1r),
     $    rdat1(np1r+1),rdat1(np1r+2)
      call calt71(temp,den,ic,m,np1r,np1i,rdat1,idat1,
     $            wav,aij,lun11,lpril)
      idest1=idat1(np1i-1+nidt-3)
      idest2=idat1(np1i+nidt-3)
      if ((idest1.le.0).or.(idest1.gt.nlev).or.
     $   (idest2.le.0).or.(idest2.gt.nlev)) go to 9000
      if (lpril.ne.0)
     $ write (lun11,*)idest1,idest2,aij,wav,ml
      ans1=aij
c
c
      ans2=0.
      if (ml.le.0) go to 9000
      nilin=npar(ml)
      if (nilin.le.0) go to 9000
      nelin=npar(nilin)
      elin=wav
      ggup=rlev(2,idest1)
      gglo=rlev(2,idest2)
      flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo)
      mlm=nelin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      a=rdat1(np1r2+1)
      elammu=elin*1.e-4
      ans1=aij*(ptmp1+ptmp2)
c     special fudge for ca i and ca ii
      if ((idat1(np1i-1+6).eq.96).or.(idat1(np1i-1+6).eq.97))
     $ ans1=min(ans1,1.e+10)
c
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5)
      sigma=(0.02655)*flin*elin*(1.e-8)/vtherm
      sigvtherm=sigma
      ener=12398.41/abs(elin)
      nb1=nbinc(ener,epi,ncn2)
      ans2=0.
c      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10
      if (elin.gt.0.99e+9) then
         ans2=0.
         sigvtherm=0.
         endif
      ans1=ans1+ans2*ggup/(1.e-36+gglo)
      opakab=sigvtherm
      ans3=ans2*ener*ergsev
      ans4=ans1*ener*ergsev
      if (elin.gt.0.1) then
        dele=12398.41/(elin+1.e-24)
        ans4=ans1*dele*(1.602197e-12)
        endif
c      ans4=0.
      if (lpril.ne.0)
     $ write (lun11,*)' ',vtherm,ans2,ans4,flin
      go to 9000
c
 72   continue
c     Autoinization rates (in s^-1) for satellite lvls
      idest1=idat1(np1i-1+nidt-3)
      idest2=idat1(np1i+nidt-3)
      temp=t*1.e+4
      call calt72(temp,np1r,rdat1,nrdt,rate)
      ans1=rate*xnx
      ans2=0.
      ggup=rlev(2,nlevp)
      gglo=rlev(2,1)
c     note that rinf has exponential removed
      rinf=(2.08e-22)*gglo/ggup/t/tsq
      dele=rdat1(np1r+1)
      ans2=rate*xnx*rinf*xnx*expo(dele/temp)
      go to 9000
c
 75   continue
c     Autoinization rates (in s^-1) for satellite lvls
c        now including final ion stage
      idest3=idat1(np1i+nidt-1)
      idest4=idat1(np1i+nidt-3)
      idest2=idat1(np1i+nidt-2)+nlev-1
      idest1=idat1(np1i-1+nidt-3)
      idest1=max(idest1,1)
      idest2=max(idest2,1)
      temp=t*1.e+4
      call calt72(temp,np1r,rdat1,nrdt,rate)
      ans1=rate*xnx
      ans2=0.
      go to 9000
c
 73   continue
c     Fit to coll. strengths satellite lvls Helike ion
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      if (rlev(1,idat1(np1i+1)).lt.rlev(1,idat1(np1i))) then
        idest2=idat1(np1i)
        idest1=idat1(np1i+1)
        endif
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      elin=abs(rdat1(np1r))
      hij=elin*1.e-8
      ekt=0.861707*t
      delt=12398.41/elin/ekt
      if (elin.le.1.e-24) go to 9000
      m=1000
      temp=t*1.e+4
c      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)
      temp=max(temp,2.8777e+6/elin)
      crate=0.
      call calt73(temp,np1r,np1i,rdat1,idat1,crate)
c      write (lun11,*)'type 73 calc:',
c     $  (rdat1(np1r-1+lk),lk=1,7),(idat1(np1i-1+lk),lk=1,4),crate,
c     $  gglo,ggup
      cijpp=crate/gglo
      cijpp=max(cijpp,0.)
      cji=(8.626e-8)*cijpp/tsq/ggup
        exptmp=expo(-delt)
       cij=cji*ggup*exptmp/gglo
      ans1=cij*xnx
      ans2=cji*xnx
      if (lpri.gt.1) then
        write (lun11,*)'ltyp=69',idest1,idest2,elin,flin,ggup,gglo
        write (lun11,*)'       ',nrdt,(rdat1(np1r-1+mm),mm=1,8),nind,jlo
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp
        endif
      elin=0.
      go to 9000
c
 74   continue
c     Delta functions to add to phot. x-sections  DR
      temp=t*1.e+4
      den=xpx
      m=1000
      rec=0.
      lprisv=lpri
c      if (lpri.ge.1) lpri=2
      if (lpri.gt.1) write (lun11,*)'type 74 data:',den,temp,
     $ (rdat1(np1r-1+mm),mm=1,nrdt),(idat1(np1i-1+mm),mm=1,nidt)
      call calt74(temp,ncn2,epi,bremsa,nrdt,np1r,rdat1,rate,
     $       alpha)
      idest1=idat1(np1i+nidt-2)
      idest2=nlevp
      idest3=idat1(np1i+nidt-1)
      idest4=idest3+1
      gglo=rlev(2,idest1)
      ggup=rlev(2,idest2)
      if (lpri.gt.1) write (lun11,*)'returning from calt74:',
     $  rate,alpha,idest1,idest2,gglo,ggup
      ans1=rate
      alpha=alpha*gglo/ggup
      ans2=alpha
      lpri=lprisv
      go to 9000
c
 81   continue
c     bhatia Fe XIX
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      if (rlev(1,idat1(np1i+1)).lt.rlev(1,idat1(np1i))) then
        idest2=idat1(np1i)
        idest1=idat1(np1i+1)
        endif
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      eeup=rlev(1,idest2)
      eelo=rlev(1,idest1)
      elin=12398.41/abs(eeup-eelo+1.e-24)
      hij=elin*1.e-8
      if (elin.le.1.e-24) go to 9000
      ekt=0.861707*t
      delt=12398.41/elin/ekt
      m=nrdt
      temp=t*1.e+4
c      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)
      temp=max(temp,2.8777e+6/elin)
      if (lpri.gt.1) then
        write (lun11,*)'ltyp=75',idest1,idest2,elin,flin,ggup,gglo
        write (lun11,*)'       ',nrdt,(rdat1(np1r-1+mm),mm=1,8),nind,jlo
        endif
      cijpp=rdat1(np1r)
      cijpp=max(cijpp,0.)
      cji=(8.626e-8)*cijpp/tsq/ggup
      ekt=0.861707*t
      delt=12398.41/elin/ekt
        exptmp=expo(-delt)
        cij=cji*ggup*exptmp/gglo
      ans1=cij*xnx
      ans2=cji*xnx
      if (lpri.gt.1) then
        write (lun11,*)'       ',cij,cji,xnx,cijpp,exptmp
        endif
      elin=0.
      go to 9000
c
 76   continue
c     2 photon decay (just  like 50)
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      lpril=lpri
      aij=rdat1(np1r)
      eeup=rlev(1,idest1)
      eelo=rlev(1,idest2)
      if (eeup.lt.eelo) then
         itmp=idest1
         idest1=idest2
         idest2=itmp
         endif
      elin=12398.41/abs(eeup-eelo)
      ggup=rlev(2,idest1)
      gglo=rlev(2,idest2)
      if (ml.le.0) go to 9000
      nilin=npar(ml)
      if (nilin.le.0) go to 9000
      nelin=npar(nilin)
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000
      flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo)
      mlm=nelin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      a=rdat1(np1r2+1)
      elammu=elin*1.e-4
c      if (flin.le.1.e-10) flin=1.
      ans1=aij
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5)
      ans2=0.
      ans4=aij*ergsev*12398.41/abs(elin)
      ansar2=0.
      em2ph=aij
      lskp=1
      emax=12398.41/elin
      nbmx=nbinc(emax,epi,ncn2)
      if (lpri.gt.1)
     $  write (lun11,*)'in ucalc, ind=76:',
     $  ml,nilin,nelin,elin,flin,rdat1(np1r),gglo,ggup,a,vtherm,ans2
     $   ,nbmx
        rcemsum=0.
        lfastl=0
        ll=1+lskp
        do while (ll.le.nbmx)
          ansar2o=ansar2
          ansar2=epi(ll)*epi(ll)*max(0.,(epi(nbmx)-epi(ll)))
          rcemsum=rcemsum+(ansar2+ansar2o)
     $                   *(epi(ll)-epi(ll-lskp))/2.
          call enxt(epi(1),nb1,lpril,epi,ncn2,t,lfastl,lun11,
     $                  ll,lskp,nphint,lrcalc)
          ll=ll+lskp
          enddo
c        rcemsum=(emax**3)/12.
        rctmp1=0.
        rctmp2=0.
        ll=2
        do while (ll.le.nbmx)
          ansar2=epi(ll)*epi(ll)*max(0.,(epi(nbmx)-epi(ll)))
          ansar2=ansar2*em2ph*emax/(1.e-24+rcemsum)
          rctmp1=abund2*ansar2*ptmp1/12.56
          rctmp2=abund2*ansar2*ptmp2/12.56
          rccemis(1,ll)=rccemis(1,ll)+rctmp1
          rccemis(2,ll)=rccemis(2,ll)+rctmp2
          call enxt(epi(1),nb1,lpril,epi,ncn2,t,lfastl,lun11,
     $                  ll,nskp,nphint,lrcalc)
          ll=ll+nskp
          enddo
        if (lpri.gt.1)
     $  write (lun11,*)'in ucalc, ind=76:',
     $  ml,nilin,nelin,elin,flin,rdat1(np1r+2),gglo,ggup,a,vtherm,ans2
        ans4=aij*ergsev*12398.41/abs(elin)
        go to 9000
c
 77   continue
c     coll rates from 71
c     Transition rates from superlevel to spect. lvls
c      go to 9000
      den=xpx
      m=1000
      clu=0.
      cul=0.
      idest1=idat1(np1i-1+nidt-3)
      idest2=idat1(np1i+nidt-3)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      eup=rlev(1,idest2)
      elo=rlev(1,idest1)
      wav=12398.41/(eup-elo+1.e-24)
      ekt=0.861707*t
      delt=wav/ekt
      lprit=0
c      if (lpri.ne.0) lprit=1
      temp=t*1.e+4
c      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)
      temp=max(temp,2.8777e+6/wav)
      call calt77(lprit,lun11,temp,den,m,np1r,np1i,rdat1,idat1,cul,clu)
      ans1=clu
      ans2=cul
      go to 9000
c
 78   continue

      go to 9000
c
 79   continue
c     fluorescence lines
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      if ((idest1.le.0).or.(idest1.gt.nlev)
     $  .or.(idest2.le.0).or.(idest2.gt.nlev))
     $      go to 9000
      elin=abs(rdat1(np1r))
      flin=rdat1(np1r+1)
c      if (flin.le.1.e-10) flin=1.
      eeup=rlev(1,idest1)
      eelo=rlev(1,idest2)
      if (eeup.lt.eelo) then
         itmp=idest1
         idest1=idest2
         idest2=itmp
         endif
      ggup=rlev(2,idest1)
      gglo=rlev(2,idest2)
      a=rdat1(np1r-1+5)
      hij=elin*1.e-8
      elammu=elin*1.e-4
      aij=(6.67e+7)*gglo*flin/ggup/elammu/elammu
c     this is a fudge to avoid badnumerics from fine structure.
      if (flin.le.1.01e-12) aij=1.e+5
      if (elin.ge.1.e+9) aij=1.e+5
      ans1=aij*(ptmp1+ptmp2)
      ans4=aij*(ptmp1+ptmp2)*ergsev*12398.41/abs(elin)
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5)
      sigma=(0.02655)*flin*elin*(1.e-8)/vtherm
      sigvtherm=sigma
c     notice that opakab does not have abundance in
      opakab=sigvtherm
c      ans2=(0.02655)*flin*elin*(1.e-8)/vtherm
      ans2=0.
      go to 9000
c
 80   continue
c Collisional ionization rates gnd of Fe and Ni
      go to 9000
c
 82   continue
c     Fe UTA rad rates
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
c     nb check this out:  no bound-bound decays from continuum
      if ((idest1.le.0).or.(idest1.ge.nlev)
     $  .or.(idest2.le.0).or.(idest2.ge.nlev))
     $      go to 9000
      gflin=rdat1(np1r+2)
      aij=rdat1(np1r-1+4)
      eeup=rlev(1,idest1)
      eelo=rlev(1,idest2)
      if (eeup.lt.eelo) then
         itmp=idest1
         idest1=idest2
         idest2=itmp
         endif
      elin=abs(rdat1(np1r))
      ggup=rlev(2,idest1)
      gglo=rlev(2,idest2)
      if (ml.le.0) go to 9000
      nilin=npar(ml)
      if (nilin.le.0) go to 9000
      nelin=npar(nilin)
      if ((nilin.le.0).or.(nelin.le.0)) go to 9000
c      flin=(1.e-16)*aij*ggup*elin*elin/((0.667274)*gglo)
      flin=gflin
      mlm=nelin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      a=rdat1(np1r2+1)
      vtherm=((vturb*1.e+5)**2+(1.29e+6/sqrt(a/t))**2)**(0.5)
      ener=12398.41/elin
      dele=ener*vtherm/3.e+10
      delev=vtherm/(elin*(1.e-8))
      delea=rdat1(np1r-1+6)*(4.14e-15)
      elammu=elin*1.e-4
      ans1=aij*(ptmp1+ptmp2)
      sigma=(0.02655)*flin/delev
c      sigvtherm=(0.02655)*flin*elin*(1.e-8)/3.e+10
      sigvtherm=sigma
      jkkl=nplini(ml)
      if (jkkl.le.0) go to 9000
      ml3=nplin(jkkl)
      if (ml3.le.0) go to 9000
      mlm=ml3-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $  nptrs,0,lun11)
      elin=abs(rdat1(np1r))
      ener=12398.41/abs(elin)
      nb1=nbinc(ener,epi,ncn2)
      ans2=sigvtherm*bremsa(nb1)*vtherm/3.e+10
c
c     turning off rex
      ans2=0.
c
c      ans4=ans1*ener*ergsev
c     notice that opakab does not have abundance in
      opakab=sigvtherm
      lfasto=2
      ans3=ans2*ener*ergsev
c     this is a cheat.  there is still an error in the 82/83 data that
c       makes some fluorescence emission
c     this test should prevent calculation when called from func2
c     since abund1 will be zero
      opakb1=sigvtherm*abund1
c      lpriu=lpri
      lpriu=0
      rcem1=0.
      rcem2=0.
c     if (opakb1.gt.1.e-34)
c    $ call linopac(lpriu,lun11,opakb1,ans2,sigvtherm,vtherm,bremsa,
c    $               rcem1,rcem2,elin,vturb,t,a,delea,epi,ncn2,
c    $               opakc,opakscatt,rccemis,fline,lfasto)
      if (lpri.gt.1)
     $  write (lun11,*)'in ucalc, ind=82:',
     $  ml,nilin,nelin,elin,flin,rdat1(np1r+2),gglo,ggup,a,vtherm,ans2,
     $  idest1,idest2,idest3,idest4,nlev,sigvtherm,bremsa(nb1),nb1
      go to 9000
c
c
 83   continue
c     Fe UTA level data
      go to 9000
c
 84   continue
      lprisv=lpri
      lpril=0
      go to 9000
c      if (lpri.ge.1) lpril=2
      if (lpril.gt.1) write (lun11,*)'ltyp=84',ml,npar(ml)
      if (lpril.gt.1) write (lun11,*)(rdat1(np1r-1+jj),jj=1,nrdt)
      if (lpril.gt.1) write (lun11,*)(idat1(np1i-1+jj),jj=1,nidt)
      if (lpril.gt.1) write (lun11,*)(kdat1(np1k-1+jj),jj=1,nkdt)
      if (ml.le.0) go to 9000
      lfastl=lfast
      nilin=npar(ml)
      idest3=idat1(np1i+nidt-1)
      idest4=idest3+1
      idest1=idat1(np1i+nidt-2)
      idest2=1
      ntmp=nrdt/2-1
      ett2=rdat1(np1r)
      ett=rdat1(np1r+2)*13.605692
      ediff=rdat1(np1r-1+2*ntmp)*13.605692-ett2
      scal2=rdat1(np1r+1)
      do ml2=1,ntmp
        etmpp(ml2)=rdat1(np1r-1+2*ml2+1)-rdat1(np1r+2)
        stmpp(ml2)=rdat1(np1r-1+2*ml2+2)*1.e-18*scal2
        stmpp(ml2)=max(stmpp(ml2),0.)
        if (lpril.gt.1) write (lun11,*)ml2,etmpp(ml2),stmpp(ml2)
        enddo
      ett=ett2-(rdat1(np1r-1+2*ntmp+1)-rdat1(np1r+2))*13.6
      ntmp2=nptmpdim
      call phextrap(etmpp,stmpp,ntmp,ntmp2,ett,ncn2,lpri,lun11)
      nb1=nbinc(ett,epi,ncn2)
      numcon2=max(2,ncn2/50)
c         numcon2=200
      nphint=ncn2-numcon2
      nphint=max(nphint,nb1+1)
      if (lpril.gt.1)
     $ write (lun11,*)'ltyp=84:',ett,ett2,ediff,ntmp,
     $  etmpp(1),stmpp(1),etmpp(ntmp+1),stmpp(ntmp+1),nb1,nphint
      if (nb1.ge.nphint-1) go to 9000
      if (lpril.gt.1)
     $ write (lun11,*)ett,nb1,bremsint(nb1),ml,vsav(1,ml),tst,
     $  lforce
       lprib=0
       lprib=lpril
c       call phint5384(stmpp,etmpp,ntmp,ett,ans1,ans2,ans3,ans4,
c     $   abund1,abund2,ptmp1,ptmp2,xpx,opakab,delr,
c     $   opakc,rccemis,lprib,epi,ncn2,bremsa,t,trad,swrat,xnx,crit53,
c     $    lfast,lun11)
      ans1=0.
      ans3=0.
      ans4=0.
      ans2=0.
      rates(1,ml)=ans1
      rates(2,ml)=ans2
      rates(3,ml)=ans3
      rates(4,ml)=ans4
      lpri=lprisv
      go to 9000

c
 85   continue
      lprisv=lpri
      lpril=0
c      if (lpri.ge.1) lpril=2
      if (lpril.gt.1) write (lun11,*)'ltyp=85',ml,npar(ml)
      if (lpril.gt.1) write (lun11,*)(rdat1(np1r-1+jj),jj=1,nrdt)
      if (lpril.gt.1) write (lun11,*)(idat1(np1i-1+jj),jj=1,nidt)
      if (lpril.gt.1) write (lun11,*)(kdat1(np1k-1+jj),jj=1,nkdt)
      if (ml.le.0) go to 9000
      lfastl=1
      nilin=npar(ml)
      idest3=idat1(np1i+nidt-1)
      idest4=idest3+1
      idest1=idat1(np1i+nidt-2)
      idest2=1
      ett2=rdat1(np1r+1)*13.605692
      nmin=idat1(np1i)
      jkk=idest3
      zc=dfloat(jkk-114)
      eion=dble(rdat1(np1r+1))
      kdim=ncn2
      far=dble(rdat1(np1r+2))
      gam=dble(rdat1(np1r-1+4))
      scal=dble(rdat1(np1r-1+5))
      call pexs(nmin,kdim,zc,eion,far,gam,scal,
     +                etmp8,stmp8,ierr,lpril,lun11)
      do mm=1,ncn2
        stmpp(mm)=sngl(stmp8(mm))*1.e-18
        enddo
      call phintfo(stmpp,ett2*0.8,ans1,ans2,ans3,ans4,
     $ abund1,abund2,xpx,opakab,
     $ opakc,lpril,epi,ncn2,bremsa,t,swrat,xnx,lfastl,lun11)
      opakab=0.
      ans4=0.
      ans2=0.
      rates(1,ml)=ans1
      rates(2,ml)=ans2
      rates(3,ml)=ans3
      rates(4,ml)=ans4
      lpri=lprisv
      go to 9000

 86   continue
c     iron auger data
c     this statement causes pileup of populations in some superlevels.
c      if (idat1(np1i-1+nidat-1).ne.idat1(np1i-1+nidat)+1) go to 9000
      ans1=rdat1(np1r+1)
      ans2=0.
      idest1=idat1(np1i-1+nidt-3)
c      idest2=nlevp
      idest2=nlevp+idat1(np1i-1+nidt-4)-1
      idest3=idat1(np1i+nidt-1)
      idest4=idest3+1
      go to 9000
c
 87   continue
      go to 9000
c
 88   continue
c     op inner shell photoexcitation
      lprisv=lpri
      idest1=idat1(np1i+nidt-2)
c      idest2=idat1(np1i+nidt-3)
      idest2=nlevp
c      if (lpri.ge.1) lpri=2
      lunsv=lun11
      if (lpri.gt.1) write (lun11,*)'ltyp=88,idest1=',idest1,idest2
      if ((idest1.ge.nlevp).or.(idest1.le.0)) go to 9000
      if (ml.le.0) go to 9000
      eth=rlev(4,idest1)-rlev(1,idest1)
      ett=eth
      nilin=npar(ml)
      if (lpri.gt.1) write (lun11,*)'nilin=',nilin,ml
      if (nilin.le.0) go to 9000
      ntmp=nrdt/2
      do ml2=1,ntmp
        etmpp(ml2)=rdat1(np1r-1+2*ml2-1)
        stmpp(ml2)=rdat1(np1r-1+2*ml2)*1.e-18
        stmpp(ml2)=max(stmpp(ml2),0.)
        enddo
      ntmp2=nptmpdim
      call phextrap(etmpp,stmpp,ntmp,ntmp2,ett,ncn2,lpri,lun11)
c      ett=ett+max(0.,13.605692*etmpp(1))
      if (lpri.gt.1) write (lun11,*)'ett=',ett,etmpp(1)
      if (ett.le.0.) go to 9000
      nb1=nbinc(ett,epi,ncn2)
      tst=abs(bremsint(nb1)/max(1.e-24,vsav(1,ml))-1.)
      xkt=ett/(0.861707*t)
      r19=rr/1.e+19
      if (lforce.ne.1) then
        xkto=vsav(4,ml)
        tq=ee1exp(xkt)/max(1.e-24,ee1exp(xkto))
        ans1=rates(1,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans2=rates(2,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        ans3=rates(3,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans4=rates(4,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        if (lpri.gt.1) write (lun11,*)'type 49 scaling:',
     $    ml,vsav(2,ml),r19,tq,xnx,vsav(3,ml),xkt,xkto,tst,
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml),ett,t
        go to 9000
        endif
      vsav(2,ml)=r19
      vsav(1,ml)=bremsint(nb1)
      vsav(3,ml)=xnx
      vsav(4,ml)=xkt
      mlm=nilin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdti,np1r2,nidti,np1i2,nkdti,np1k2,mlm,
     $  nptrs,0,lun11)
      emax=etmpp(ntmp)*13.6+eth
      gglo=rlev(2,idest1)
      ggup=rlev(2,idest2)
      idest3=idat1(np1i-1+nidti)
      idest4=idest3+1
      if (lpri.gt.1) write (lun11,*)'before phint53',gglo,ggup
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      if (lpri.gt.1) then
        write (lun11,*)'type 88 data:',idat1(np1i),idat1(np1i+nidt-1),
     $           t,xnx,eth,gglo,ggup,swrat
        call dprinto(ndesc,nrdesc,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,rdat1,idat1,kdat1,lun11)
        endif
c      if ((lpri.ge.1).and.(idest1.le.4).and.(jkion.eq.29)
c     $    .and.(abund1.gt.1.e-34))  then
c        lun99=99
c        write (lun99,*)'type 88 data:', idest1, idest2,
c     $           eth,gglo,ggup,swrat, abund1,abund2
c        call dprinto(ndesc,nrdesc,lcon,
c     $          nrdt,rdat,nidt,idat,nkdt,kdat,lun99)
c        do mm=1,ncn2
c          opaksv(mm)=opakc(mm)
c          opakc(mm)=0.
c          enddo
c        endif
      lprib=lpri
      lprib=0
      if (lpri.gt.1) lprib=lpri
      rnist=rniss(idest1)*exp(-ett/(0.861707)/t)/rniss(nlevp)
      call phint53(stmpp,etmpp,ntmp,ett,ans1,ans2,ans3,ans4,
     $  abund1,abund2,ptmp1,ptmp2,xpx,opakab,rnist,
     $  opakc,rccemis,lprib,epi,ncn2,bremsa,t,swrat,xnx,
     $  lfast,lun11)
c      if ((lpri.ge.1).and.(idest1.le.4).and.(jkion.eq.29)
c     $    .and.(abund1.gt.1.e-34))  then
c        nhit=0
c        do mm=1,ncn2
c          if ((opakc(mm).gt.1.e-34).and.(epi(mm).gt.500.)
c     $        .and.(epi(mm).lt.800.)) then
c            write (lun99,919)mm,epi(mm),opakc(mm),
c     $        opakc(mm)/max(1.e-34,abund1)/xpx,opaksv(mm)+opakc(mm)
c            nhit=1
c            endif
c          opakc(mm)=opaksv(mm)+opakc(mm)
c          enddo
c        if (nhit.eq.0) write (lun99,*)'no cross section'
c        endif
      if (lpri.gt.1) then
        npr=nb1
        write (lun11,*)'bautista threshold xsection:',
     $         npr,ett,eth,rdat1(np1r),sg(npr),ans2,swrat
        endif
      rates(1,ml)=ans1
      ans2=0.
      ans4=0.
      rates(2,ml)=ans2
      rates(3,ml)=ans3
      rates(4,ml)=ans4
      if (lpri.gt.1) write (lun11,*)'rates:',
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml)
      lun11=lunsv
      lpri=lprisv
      go to 9000
c
 89   continue
      go to 9000
c
 90   continue
      go to 9000
c
 91   continue
c     a values from atomdb.  same as 50.
      go to 50
      go to 9000
c
 92   continue
c     collision strengths from atomdb
c
      lpril=0
c      if (lpri.ge.1) lpril=2
      if (lpril.gt.1) write (lun11,*)'ltyp=92',ml,npar(ml)
      if (lpril.gt.1) write (lun11,*)(rdat1(np1r-1+jj),jj=1,nrdt)
      if (lpril.gt.1) write (lun11,*)(idat1(np1i-1+jj),jj=1,nidt)
      if (lpril.gt.1) write (lun11,*)(kdat1(np1k-1+jj),jj=1,nkdt)
c
c
c     general stuff
      lctype=idat1(np1i+3-1)
      idest1=idat1(np1i)
      idest2=idat1(np1i+1)
      tmin=rdat1(np1r)
      tmax=rdat1(np1r+1)
      do mml=1,20
        tstr(mml)=rdat1(np1r+1+mml)
        cstr(mml)=rdat1(np1r+21+mml)
        enddo
      ggup=rlev(2,idest2)
      gglo=rlev(2,idest1)
      eeup=rlev(1,idest2)
      eelo=rlev(1,idest1)
      elin=12398.41/abs(eeup-eelo+1.e-24)
      hij=elin*1.e-8
      if (elin.le.1.e-24) go to 9000
      ekt=0.861707*t
      temp=t*1.e+4
      eij=abs(eeup-eelo)
      eijkev=eij/1.e+3
      tk=t*1.e+4
      mlm=npar(ml)-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      nistage=idat1(np1i2)
      mlm=npar(mlm)-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      nzel=idat1(np1i2)
      if (lpril.gt.1)
     $ write (lun11,*)'before calc_maxwell_rates',lctype,tmin,tmax,
     $   eijkev,tk,zzz,gglo,ggup
      call calc_maxwell_rates(lun11,lpril,lctype,tmin,tmax,
     $  Tstr,cstr, eijkev,  tk, nzel,  gglo,  ggup,  cij, cji)
      ans1=cij*xnx
      ans2=cji*xnx
      if (lpril.gt.1) then
        write (lun11,*)'type 92 data',lctype,cij,cji,xnx
        endif
       go to 9000
c
c      old code for 92 not used
c      delt=12398.41/elin/ekt
c      temp=max(temp,(1.e+4)*12398.54/elin/(0.861707)/50.)
c      eijry=eij/13.605692
c      tsq=sqrt(t)
c      cijpp=0.
c      temp=max(temp,2.8777e+6/elin)
c      if (lctype.eq.11) then
cc       chianti type 1 (dere et al. 1997)
c        cijpp=upsil(1,eijry,cijpp,cstr(1),cstr(2),
c     $      cstr(3),cstr(4),cstr(5),tk)
c        endif
c      if (lctype.eq.12) then
cc       chianti type 2 (dere et al. 1997)
c        cijpp=upsil(2,eijry,cijpp,cstr(1),cstr(2),
c     $      cstr(3),cstr(4),cstr(5),tk)
c        endif
c      if (lctype.eq.13) then
cc       chianti type 3 (dere et al. 1997)
c        cijpp=upsil(3,eijry,cijpp,cstr(1),cstr(2),
c     $      cstr(3),cstr(4),cstr(5),tk)
c        endif
c      if (lctype.eq.14) then
cc       chianti type 4 (dere et al. 1997)
c        cijpp=upsil(4,eijry,cijpp,cstr(1),cstr(2),
c     $      cstr(3),cstr(4),cstr(5),tk)
c        endif
c      if (lctype.eq.31) then
cc       sampson goett and clark 1983 type 1
c        if (nrdt.lt.7) go to 9000
c        y=eij/ekt
c        aa=rdat1(np1r+2)
c        co=rdat1(np1r+3)
c        cr=rdat1(np1r+4)
c        crp=rdat1(np1r+5)
c        rr=rdat1(np1r+6)
c        sig=rdat1(np1r+6)
c        z2s=rdat1(np1r+7)
c        zeff=float(idat1(np1i+2))-sig
c        if (y.gt.40.)  go to 9000
c        call expint(y,em1)
c        e1=em1/y*exp(-y)
c        if (y*a+y.le.80) then
c            call eint(y*a+y,ee1,ee2,ee3)
c          else
c            ee1=0.
c            ee2=0.
c            ee3=0.
c          endif
c        er=0.
c        er1=0.
c        if (rr.eq.1.) then
c          er=ee1
c          er1=ee2
c          endif
c        if (rr.eq.2.) then
c          er=ee2
c          er1=ee3
c          endif
c        if (y*a+y.le.40) then
c            qij=co*exp(-y)+1.55*z2s*e1+y*exp(y*a)*(cr*er/(a+1.)**(rr-1.)
c     #      +cr1*er1/(a+1.)**rr)
c          else
c            qij=co*exp(-y)+1.55*z2s*e1
c          endif
c        cijpp=qij*exp(y)/zeff/zeff
c        endif
c      if (lctype.eq.32) then
cc       sampson goett and clark 1983 type 2
c        endif
c      if (lctype.eq.33) then
cc       sampson goett and clark 1983 type 3
c        endif
c      if (lctype.eq.41) then
cc       kato and nakazaki 1989 type 1
c        call calt66(temp,np1r+2,rdat1,gamma)
c        cijpp=gamma
c        endif
c      if (lctype.eq.42) then
cc       kato and nakazaki 1989 type 3
c        go to 9000
c        endif
c      if (lctype.gt.100) then
c        ncase=int(lctype/50)
cc       ltype=100 -->ncase=2
cc       ltype=150 -->ncase=3
cc       ltype=200 -->ncase=4
cc       ltype=250 -->ncase=5
cc       ltype=300 -->ncase=6
cc       ltype=350 -->ncase=7
cc       ltype=400 -->ncase=8
cc       ltype=450 -->ncase=9
cc       ltype=500 -->ncase=10
cc       ltype=550 -->ncase=11
cc       ltype=600 -->ncase=12
cc       ltype=650 -->ncase=13
cc       ltype=700 -->ncase=14
cc       ltype=750 -->ncase=15
cc       ltype=800 -->ncase=16
cc       ltype=850 -->ncase=17
cc       ltype=900 or greater -->ncase=18
cc       don't do qs
c        if ((ncase.ge.6).and.(ncase.le.9)) stop 'ncase=6-9'
cc        if ((ncase.ge.6).and.(ncase.le.9)) go to 9000
c        if (ncase.ge.14) stop 'ncase=14'
cc        if (ncase.ge.14) go to 9000
cc
c        npts=lctype-50*ncase
c        if ((tk.le.tmin).or.(tk.ge.tmax)) go to 9000
c        mm=1
c        do while ((tk.lt.tstr(mm)).and.(mm.lt.npts))
c          mm=mm+1
c          enddo
c        mm=max(mm-1,1)
c        cijpp=cstr(mm)+(cstr(mm+1)-cstr(mm))*(tk-tstr(mm))
c     $                 /(tstr(mm+1)-tstr(mm)+1.e-38)
c        endif
c
c      cji=(8.626e-8)*cijpp/tsq/ggup
c      ekt=0.861707*t
c      delt=12398.41/elin/ekt
c      exptmp=expo(-delt)
c      cij=cji*ggup*exptmp/gglo
c
 93   continue
      go to 9000
c     op pi xsections
      lprisv=lpri
      if (nrdt.gt.3) go to 533
      idest1=idat1(np1i+nidt-2)
      idest2=nlevp+idat1(np1i-1+nidt-3)-1
      if (lpri.gt.1) write (lun11,*)'idest1=',idest1,idest2,nlevp,ml
      if ((idest1.ge.nlevp).or.(idest1.le.0)) go to 9000
      if (ml.le.0) go to 9000
      eth=rlev(4,idest1)-rlev(1,idest1)
      eexc=rlev(1,idest1)
      ett=eth
      nilin=npar(ml)
      if (lpri.gt.1) write (lun11,*)'nilin=',nilin,ml
      if (nilin.le.0) go to 9000
      ntmp=nrdt/2
c      ett=ett+max(0.,13.605692*etmpp(1))
      if (lpri.gt.1) write (lun11,*)'ett=',ett,etmpp(1)
      if (ett.le.0.) go to 9000
      nb1=nbinc(ett,epi,ncn2)
      tst=abs(bremsint(nb1)/max(1.e-24,vsav(1,ml))-1.)
      xkt=ett/(0.861707*t)
      r19=rr/1.e+19
c      if ((tst.le.0.01).and.(lforce.ne.1)) then
      if (lforce.ne.1) then
        xkto=vsav(4,ml)
        tq=ee1exp(xkt)/max(1.e-24,ee1exp(xkto))
        ans1=rates(1,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans2=rates(2,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        ans3=rates(3,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans4=rates(4,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        if (lpri.gt.1) write (lun11,*)'type 93 scaling:',
     $    ml,vsav(2,ml),r19,tq,xnx,vsav(3,ml),xkt,xkto,tst,
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml),ett,t
        go to 9000
        endif
      vsav(2,ml)=r19
      vsav(1,ml)=bremsint(nb1)
      vsav(3,ml)=xnx
      vsav(4,ml)=xkt
      mlm=nilin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      gglo=rlev(2,idest1)
      ggup=rlev(2,nlevp)
      idest3=idat1(np1i+nidt-1)
      idest4=idest3+1
      if (idest2.gt.nlevp) then
        jkk3=jkion+1
        if (lpri.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        ndtmp=npfi(13,jkk3)
        if (lpri.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        if (ndtmp.le.0) go to 9000
        mllz=npar(ndtmp)
        iltmp=0
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))
     $      .and.(npar(ndtmp).eq.mllz))
           mlm=ndtmp-1
           call drd(ltyp2,lrtyp2,lcon2,
     $       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $       nptrs,0,lun11)
           iltmp=idat1(np1i2+nidt2-2)
           if (lpri.gt.1) write (lun11,*)nidt2,iltmp,ndtmp
           ndtmp=npnxt(ndtmp)
           if (ndtmp.le.0) go to 9000
           enddo
c        NB fix to excited level PI and rec
         ett=ett+rdat1(np1r2)
         eth=ett
         ggup=rdat1(np1r2+1)
         if (lpri.gt.1)
     $    write (lun11,*) ndtmp,iltmp,idest2,ggup,ett
         endif
      if (lpri.gt.1) write (lun11,*)'before phint53',eexc,eth,lfast
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      if (lpri.gt.1) then
        write (lun11,*)'type 93 data:',idat1(np1i),idat1(np1i+nidt-1),
     $           t,xnx,eth,gglo,ggup,swrat
        call dprinto(ndesc,nrdesc,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,rdat1,idat1,kdat1,lun11)
        endif
      lprib=0
      if (lpri.gt.1) lprib=lpri
      sth=1.e-18*rdat1(np1r+1)
      alph=rdat1(np1r+2)
      e1=rdat1(np1r)
      lfastl=1
      call phint53pl(sth,e1,alph,ett,ans1,ans2,ans3,ans4,
     $  abund1,abund2,ptmp1,ptmp2,xpx,opakab,
     $  opakc,rccemis,lprib,epi,ncn2,bremsa,t,swrat,xnx,
     $  lfastl,lun11)
      if (lpri.gt.1) then
        npr=nb1
        write (lun11,*)'bautista threshold xsection:',
     $         npr,ett,eth,rdat1(np1r),sg(npr),ans2,swrat
        endif
      rates(1,ml)=ans1
      rates(2,ml)=ans2
      rates(3,ml)=ans3
      rates(4,ml)=ans4
      if (lpri.gt.1) write (lun11,*)'rates:',
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml)
      lpri=lprisv
      go to 9000
c
 94   continue
      go to 9000
c     op pi xsections
c     old version
      lprisv=lpri
c      if (lpri.ge.1) lpri=2
      if (nrdt.gt.3) go to 499
      idest1=idat1(np1i+nidt-2)
      idest4=idat1(np1i+nidt-3)
      idest2=nlevp+idat1(np1i+nidt-4)-1
      if (lpri.gt.1) write (lun11,*)'idest1=',idest1,idest2
      if ((idest1.ge.nlevp).or.(idest1.le.0)) go to 9000
      if (ml.le.0) go to 9000
      eth=rlev(4,idest1)-rlev(1,idest1)
      ett=eth
      nilin=npar(ml)
      if (lpri.gt.1) write (lun11,*)'nilin=',nilin,ml
      if (nilin.le.0) go to 9000
      if (lpri.gt.1) write (lun11,*)'ett=',ett,etmpp(1)
      if (ett.le.0.) go to 9000
      nb1=nbinc(ett,epi,ncn2)
      tst=abs(bremsint(nb1)/max(1.e-24,vsav(1,ml))-1.)
      xkt=ett/(0.861707*t)
      r19=rr/1.e+19
c      if ((tst.le.0.01).and.(lforce.ne.1)) then
      if (lforce.ne.1) then
        xkto=vsav(4,ml)
        tq=ee1exp(xkt)/max(1.e-24,ee1exp(xkto))
        ans1=rates(1,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans2=rates(2,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        ans3=rates(3,ml)*vsav(2,ml)*vsav(2,ml)/r19/r19
        ans4=rates(4,ml)*tq*xnx/max(1.e-24,vsav(3,ml))
        if (lpri.gt.1) write (lun11,*)'type 94 scaling:',
     $    ml,vsav(2,ml),r19,tq,xnx,vsav(3,ml),xkt,xkto,tst,
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml),ett,t
        go to 9000
        endif
      vsav(2,ml)=r19
      vsav(1,ml)=bremsint(nb1)
      vsav(3,ml)=xnx
      vsav(4,ml)=xkt
      mlm=nilin-1
      call drd(ltyp,lrtyp,lcon,
     $  nrdt,np1r2,nidt,np1i2,nkdt,np1k2,mlm,
     $  nptrs,0,lun11)
      ntmp=nrdt/2
      gglo=rlev(2,idest1)
      ggup=rlev(2,nlevp)
      idest3=idat1(np1i+nidt-1)
      idest4=idest3+1
      if (idest2.gt.nlevp) then
        jkk3=jkion+1
        if (lpri.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        ndtmp=npfi(13,jkk3)
        if (lpri.gt.1)
     $    write (lun11,*)jkk3,ndtmp,nlevp,idest2
        mllz=npar(ndtmp)
        nptmp=mllz
        do while ((ndtmp.ne.0).and.(iltmp.ne.(idest2-nlevp+1))
     $      .and.(nptmp.eq.mllz))
           mlm=ndtmp-1
           call drd(ltyp2,lrtyp2,lcon2,
     $       nrdt2,np1r2,nidt2,np1i2,nkdt2,np1k2,mlm,
     $       nptrs,0,lun11)
           iltmp=idat1(np1i2+nidt2-2)
           if (lpri.gt.1) write (lun11,*)nidt2,iltmp,ndtmp
           ndtmp=npnxt(ndtmp)
           nptmp=0
           if (ndtmp.ne.0) nptmp=npar(ndtmp)
           enddo
         ggup=rdat1(np1r2+1)
         if (lpri.gt.1)
     $    write (lun11,*) ndtmp,iltmp,idest2,ggup
         endif
      if (lpri.gt.1) write (lun11,*)'before phint53'
      if (ggup.le.1.e-24) then
        write (lun11,*) 'ggup error'
        return
        endif
      swrat=gglo/ggup
      if (lpri.gt.1) then
        write (lun11,*)'type 94 data:',idat1(np1i),idat1(np1i+nidt-1),
     $           t,xnx,eth,gglo,ggup,swrat
        call dprinto(ndesc,nrdesc,lcon,
     $          nrdt,np1r,nidt,np1i,nkdt,np1k,rdat1,idat1,kdat1,lun11)
        endif
      lprib=0
      if (lpri.gt.1) lprib=lpri
      sth=1.e-18*rdat1(np1r+1)
      alph=rdat1(np1r+2)
      e1=rdat1(np1r)
      lfastl=1
      call phint53pl(sth,e1,alph,ett,ans1,ans2,ans3,ans4,
     $  abund1,abund2,ptmp1,ptmp2,xpx,opakab,
     $  opakc,rccemis,lprib,epi,ncn2,bremsa,t,swrat,xnx,
     $  lfastl,lun11)
      if (lpri.gt.1) then
        npr=nb1
        write (lun11,*)'bautista threshold xsection:',
     $         npr,ett,eth,rdat1(np1r),sg(npr),ans2,swrat
        endif
      rates(1,ml)=ans1
      rates(2,ml)=ans2
      rates(3,ml)=ans3
      rates(4,ml)=ans4
      if (lpri.gt.1) write (lun11,*)'rates:',
     $    rates(1,ml),rates(2,ml),rates(3,ml),rates(4,ml)
      lpri=lprisv
      go to 9000
c
 95   continue
c     bryans ci rates
      if (nrdesc.eq.5) then
          idest2=nlevp
        else
          idest2=1
        endif
      ee=rdat1(np1r)
      tmin=rdat1(np1r+1)
      nspline=(nrdt-2)/2
      ekt=0.861707*t
      tt=ekt/ee
c     constant is ln(2)
      xx=1.-(0.693147)/log(tt+2.)
      mm=1
      do while ((mm.lt.nspline).and.(xx.gt.rdat1(np1r+1+mm)))
        if (lpri.gt.1) write (lun11,*)mm,xx,rdat1(np1r+1+mm)
        mm=mm+1
        enddo
c      this illustrates the storage scheme for the splines
c      do mm=1,nspline
c        tspline(mm)=rdat1(np1r+1+mm)
c        vspline(mm)=rdat1(np1r+1+nspline+mm)
c        enddo
c     linear interpolation
      rho=(rdat1(np1r+1+nspline+mm-1)
     $  +(xx-rdat1(np1r+1+mm-1))*
     $    (rdat1(np1r+1+nspline+mm)-rdat1(np1r+1+nspline+mm-1))
     $    /(rdat1(np1r+1+mm)-rdat1(np1r+1+mm-1)))
c      dere equation 7
      call eint(1./tt,e1,e2,e3)
      citmp1=1.e-6*e1*rho/sqrt(tt*ee**3)
      ans1=citmp1*xnx
      ans2=0.
      idest1=1
c      idest2=1
      ggup=rlev(2,nlevp)
      gglo=rlev(2,nidt-1)
c     note that rinf has exponential removed
      rinf=(2.08e-22)*gglo/ggup/t/tsq
      ans2=ans1*rinf*xnx
c      idest1=idat1(np1i+nidt-2)
      idest4=idat1(np1i+nidt-1)+1
      idest3=idat1(np1i+nidt-1)
      if (lpri.gt.1)
     $ write (lun11,*)'ltype=95:',ee,tt,nspline,xx,mm,rho,e1,citmp1,
     $  ekt,ans1,ans2
      go to 9000


 9000 continue
c
      call remtms(time2)
c      write (lun11,*)'ndesc=',ndesc
!jg      tucalc(ndesc)=tucalc(ndesc)+abs(time2-time1)
!jg      ncall(ndesc)=ncall(ndesc)+1
c
      if (lpri.gt.1)
     $ write (lun11,9931)krdesc(nrdesc),kdesc(ndesc),ndesc,ans1,ans2,
     $     ans3,ans4,idest1,idest2,rdat1(np1r)
 9931 format (1x,'in ucalc :',a28,a56,i4,4x,4(1pe10.2),2i4,3(1pe10.2))
c
      return
      end
      subroutine uclgsi(kdum,iresult,ierr)
      implicit none
      integer iresult, ierr
      character*(*) kdum
c
      ierr=0
c
      read (5,*)iresult
!      write (6,*)'in uclgsi, iresult=',iresult
      return
      end
      subroutine uclgsr(kdum,result,ierr)
      implicit none
      integer ierr
      real result
      character*(*) kdum
c
      ierr=0
c
      read (5,*)result
!      write (6,*)'in uclgsr, result=',result
      return
      end
      subroutine uclgst(kdum,kresult,ierr)
      implicit none
      integer ierr, ll
      character*(*) kdum,kresult
      character*80 kres2
      character*1 ktmp
c
      ierr=0
c
      read (5,'(a80)')kres2
      do ll=1,80
        read (kres2(ll:ll),'(a1)')ktmp
        write(kresult(ll:ll),'(a1)')ktmp
        enddo
      return
      end
      subroutine uclgsr8(kdum,result,ierr)
      implicit none
      integer ierr
      real*8 result
      real result4
      character*(*) kdum
c
      ierr=0
c
      call uclgsr(kdum,result4,ierr)
      result=result4
c
      return
      end
      subroutine unsavd(jkstep,ldir,
     $       lpri,iunit,iunit2,iunit3,iunit4,
     $       idat1,rdat1,kdat1,nptrs,npnxt,npfi,
     $       npfirst,npar,npilev,npconi2,ncsvn,
     $       t,p,r,rdel,delr,xcol,xee,xpx,zeta,
     $       xilev,bilev,rniss,
     $       nplin,nlsvn,rcem,oplin,tau0,
     $       cemab,cabab,opakab,tauc,
     $       epi,ncn2,dpthc,opakc,rccemis,nloopctl,
     $       lunlog,status)
c
c     this routine  saves only depths for iterative calculation
c     author:  T. Kallman
c
      implicit none
c
      include './PARAM'
c
C     Allocation for passed parameters
      real*8 rdat1(nrdat1)
      real*8 r,delr,rdel, t, p, xcol,xee,xpx,zeta
      integer unit,hdunum, nrows, status, nloopctl
      integer iunit,iunit2,iunit3,iunit4
      integer idat1(nidat1),nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
      integer nlsvn,ncsvn
c     energy bins
      real*8 epi(ncn)
c     continuum opacities
      real*8 opakc(ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn)
      integer ncn2
c     line opacities
      real*8 oplin(nnnl)
      real*8 abel(nl)
      real*8 tauc(2,nnml)
      real*8 cemab(2,nnml),opakab(nnml),cabab(nnml)
c     pointers to master data
      integer npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni),npar(ndat2)
      integer npilev(nd,nni)
      integer nplin(nnnl)
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      real*8 xilev(nnml),bilev(nnml),rniss(nnml)
      real*8 tau0(2,nnnl), rcem(2,nnnl)

c     continuum optical depths
      integer ldir,lpri,lun,lunlog,jkstep
      real*8 tau0d(2,nnnl),dpthcd(2,ncn),taucd(2,nnml)
      integer lind1,lind2,kl,ll,nlyc,nry,nbinc
c
      r=0.
      delr=0.
      t=0.
      p=0.
      if (status .gt. 0)call printerror(lunlog,status)
      call rstepr(iunit,jkstep,r,delr,rdel,t,p,
     $          xcol,xee,xpx,zeta,
     $          idat1,rdat1,kdat1,nptrs,npnxt,npfi,
     $          npfirst,npar,npilev,
     $          xilev,bilev,rniss,nloopctl,
     $          lunlog,status)
      call rstepr2(iunit2,jkstep,r,delr,rdel,t,p,
     $          xcol,xee,xpx,zeta,
     $          idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $          nplin,nlsvn,rcem,oplin,tau0d,nloopctl,
     $          lunlog,status)
      call rstepr3(iunit3,jkstep,r,delr,rdel,t,p,
     $          xcol,xee,xpx,zeta,
     $          idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $          npfirst,npilev,npconi2,ncsvn,
     $          rniss,cemab,cabab,opakab,taucd,nloopctl,
     $          lunlog,status)
 9000 continue
      call rstepr4(iunit4,jkstep,r,delr,rdel,t,p,
     $          xcol,xee,xpx,zeta,
     $          idat1,rdat1,kdat1,nptrs,npnxt,npfi,npar,
     $          epi,ncn2,dpthcd,opakc,rccemis,nloopctl,
     $          lunlog,status)
c      write (lunlog,*)' r=',r,xpx,t,delr,xcol,xee,xpx,zeta
      if (status .gt. 0)call printerror(lunlog,status)
c      read(lun,REC=jkstep)t,r,rdel,delr,xcol,xee,xpx,tau0d,dpthcd,taucd
      lind1=1
      lind2=2
      if (ldir.gt.0) lind2=1
      if (ldir.lt.0) lind1=2
      do ll=lind1,lind2
        do kl=1,nnnl
          tau0(ll,kl)=tau0d(ll,kl)
          enddo
        do kl=1,ncn2
          dpthc(ll,kl)=dpthcd(ll,kl)
          enddo
        do kl=1,nnml
          tauc(ll,kl)=taucd(ll,kl)
          enddo
        enddo
      nlyc=nbinc(13.7d0,epi,ncn2)
      nry=nlyc+1
      if (lpri.ne.0) write (lunlog,*)'in unsavd',rdel,t,tauc(1,25),
     $                ldir,dpthc(1,nry),dpthc(2,nry)
c
      return
      end
      function upsil(k,eij,c,p1,p2,p3,p4,p5,t)
      implicit none
      real*8 e, eij, c, p1, p2, p3, p4, p5, t
      real*8 y, splinem, upsil, x
      integer k
c
c     this routine calculates upsilons for Burgess and Tully
c     author:  M. Bautista
c
c     t = electron temperature in Kelvin
c     p# = spline knot values
c     c = abscissa scale parameter
c     k = transition type
c     eij = transition energy (Ryd)
c
       e=abs(t/(1.57888e5*eij))                  !<<<<<< CORRECTED LINE
       if ((k.eq.1).or.(k.eq.4.)) x=log((e+c)/c)/log(e+c)
       if ((k.eq.2).or.(k.eq.3)) x=e/(e+c)
       y=splinem(p1,p2,p3,p4,p5,x)
       if (k.eq.1) y=y*log(e+2.71828)
       if (k.eq.3) y=y/(e+1)
       if (k.eq.4) y=y*log(e+c)
       upsil=y
c
      return
      end
      subroutine velimp(n,l,temp,ic,z1,rm,ne,sum,cn)
      implicit none
c
c     impact parameter collision rate calculated following the method of
c     pengelly & seaton (1964) but using the lowest cross-section at every
c     velocity.
c     note that cn is the rate for nl -> nl-1 and hence l > 0 *
c     cne(l+1)=cn
c     cen(l)=cn*(2.*l+1)/(2.*l-1)
c     author:  M. Bautista
c
      real*8 ne, temp, z1, rm, sum, cn
      real*8 pi, pa, pd, alfa, b, dnl, bb
      real*8 va, vd, ava, vb, avb, avd, eb
      real*8 ea, ed, xa, xb, xd, expo, den
      real*8 ca, cad, cd
      integer n, l, ic
c
      cn=0.
      if((l.eq.0).or.(sum.eq.0.)) go to 50
      den=l*(n*n-l*l)+(l+1)*(n*n-(l+1)*(l+1))
      dnl=6.*z1/ic*z1/ic*n*n*(n*n-l*l-l-1)
      pi=2.*acos(0.)
      pa=0.72/sum
      pd=6.90*sqrt(temp/ne)
      alfa=3.297e-12*rm/temp
      b=1.157*sqrt(dnl)
      bb=b*b
c
      va=pd/pa
      vd=b/pd
      vb=sqrt(va*vd)
c
           ava=alfa*va*va
           avb=alfa*vb*vb
           avd=alfa*vd*vd
      ea=0.
      ed=0.
      xa=expo(-ava)
      xb=expo(-avb)
      xd=expo(-avd)
           if(ava.lt.50.) call expint(ava,ea)
c           call expint(ava,ea)
      ea=ea/ava*xa
           call expint(avb,eb)
      eb=eb/avb*xb
           if(avd.lt.50.) call expint(avd,ed)
c           call expint(avd,ed)
      ed=ed/avd*xd
c
      if(va.gt.vd) then
      if(avb.gt.1.e-3) then
      cn=sqrt(pi*alfa)*(pa*pa*(2./alfa/alfa-xb*(vb**4+2.*vb*vb/alfa+
     #2./alfa/alfa))+bb*xb+2.*bb*eb-bb*ea)
           else
      cn=sqrt(pi*alfa)*bb*(1.+avb*(1./3.-avb/4.)+2.*eb-ea)
                endif
c
           else
      if(ava.gt.1.e-3) then
      ca=sqrt(pi*alfa)*pa*pa*(2./alfa/alfa-xa*(va**4+2.*va*va/alfa+
     #2./alfa/alfa))
           else
      ca=sqrt(pi*alfa)*pd*pd*va**4*alfa*(1/3.-ava/4.+ava*ava/10.)
                 endif
c
      cad=sqrt(pi*alfa)*pd*pd/alfa*(xa*(1.+ava)-xd*(1.+avd))
      cd=sqrt(pi*alfa)*bb*(xd+ed)
      cn=ca+cad+cd
                 endif
c
      cn=cn*l*(n*n-l*l)/den
c
 50      return
c
      end
      function voigte(vs,a)
c      real*8 function voigtevoigt(vs,a)
      implicit none
c
c     computes a voigt function  h = h(a,v)
c     a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
c     traving (landolt-b\rnstein, p. 449).
c     author:  tlusty
c
c     the integral of this function over vs (-infinity -> infinity) is sqrt(pi)
c
      real*8 vs, a, un, two, voigte, sqp, sq2
      real*8 v, u, v2, ex, quo, h1, h, pqs, h1p
      real*8 h2p, h3p, h4p, psi, a2, u2
      real*8 ak(19),a1(5)
      integer k, m, i
      parameter (un=1., two=2.)
c
      data ak      /-1.12470432, -0.15516677,  3.28867591, -2.34357915,
     ,  0.42139162, -4.48480194,  9.39456063, -6.61487486,  1.98919585,
     , -0.22041650, 0.554153432, 0.278711796,-0.188325687, 0.042991293,
     ,-0.003278278, 0.979895023,-0.962846325, 0.532770573,-0.122727278/
      data sqp/1.772453851/,sq2/1.414213562/
c
      v = abs(vs)
      u = a + v
      v2 = v*v
      if (a.eq.0.0) go to 140
      if (a.gt.0.2) go to 120
      if (v.ge.5.0) go to 121
c
      ex=0.
      if(v2.lt.100.)ex = exp(-v2)
      k = 1
c
  100 quo = un
      if (v.lt.2.4) go to 101
      quo = un/(v2 - 1.5)
      m = 11
      go to 102
c
  101 m = 6
      if (v.lt.1.3) m = 1
  102 do 103 i=1,5
         a1(i) = ak(m)
         m = m + 1
  103 continue
      h1 = quo*(a1(1) + v*(a1(2) + v*(a1(3) + v*(a1(4) + v*a1(5)))))
      if (k.gt.1) go to 110
c
c a le 0.2  and v lt 5.
c
      h = h1*a + ex*(un + a*a*(un - two*v2))
      voigte=h
      return
c
  110 pqs = two/sqp
      h1p = h1 + pqs*ex
      h2p = pqs*h1p - two*v2*ex
      h3p = (pqs*(un - ex*(un - two*v2)) - two*v2*h1p)/3. + pqs*h2p
      h4p = (two*v2*v2*ex - pqs*h1p)/3. + pqs*h3p
      psi = ak(16) + a*(ak(17) + a*(ak(18) + a*ak(19)))
c
c 0.2 lt a le 1.4  and  a + v le 3.2
c
      h = psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p))))
      voigte=h
      return
c
  120 if (a.gt.1.4.or.u.gt.3.2) go to 130
      ex=0.
      if(v2.lt.100.)ex = exp(-v2)
      k = 2
      go to 100
c
c a le 0.2  and  v ge 5.
c
  121 h = a*(15. + 6.*v2 + 4.*v2*v2)/(4.*v2*v2*v2*sqp)
      voigte=h
      return
c
  130 a2 = a*a
      u = sq2*(a2 + v2)
      u2 = un/(u*u)
c
c a gt 1.4  or  a + v gt 3.2
c
      h = sq2/sqp*a/u*(1. + u2*(3.*v2 - a2) +
     ,        u2*u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2))
      voigte=h
      return
c
c a eq 0.
c
  140 h=0.
      if(v2.lt.100.)h=exp(-v2)
      voigte=h
      return
      end
      real*8 function voigtedoppler(vs,a)
c
c     fake version does doppler
c     the integral of this function over vs (-infinity -> infinity) is sqrt(pi)
c
      real*8 a,vs
c
      voigtedoppler=exp(-vs*vs)
c
      return
      end
c      real*8 function voigte(vs,a)
      real*8 function voigtelorentzian(vs,a)
c
c     fake version does lorentzian
c
c     the integral of this function over vs (-infinity -> infinity) is 1
c
      real*8 a,vs
c
      voigtelorentzian=a/(vs*vs+a*a)/3.14
C
      Return
      end
      subroutine writespectra(lun11,lpri,nparms,
     $       parname,partype,parval,parcomm,atcredate,
     $       t,vturbi,epi,ncn2,dpthc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,elum,zrems,zremsz,kmodelname,nloopctl)
c
C     Write extension containing the spectra for
C     this particular model.
C
C     Modifications:
C       04/01/1999,WTB: Disabled appending loop control value to
C               extension name due to changes in xstar2xspec design
c       051/17/2003 TK added auger damping
c     author:  T. Bridgman
C
c
      implicit none
      include './PARAM'
c
c     passed parameters
      character(30) kmodelname
      integer nparms, nloopctl, lun11
      character(20) parname(55)
      character(10) partype(55)
      real*8 parval(55)
      character(30) parcomm(55)
c     master data
      integer idat1(nidat1),nptrs(nptt,ndat2)
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     line luminosities
      real*8 elum(3,nnnl)
c     energy bins
      real*8 epi(ncn)
c     continuum lum
      real*8 zrems(4,ncn),zremsz(ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
      real*8 zrtmp(6,ncn)
      real rtmp(ncn)
      character(16) knam,klabs(6),kunits(6),kform(6),kblnk16
      character(30) extname
      integer unit,istatus
      integer nlsvn, ll
      integer np2, tbcol(6), nrows, rowlen, kk
      integer frow, felem, colnum, tfields, status, verbose,mm
      real*8 eliml, elimh
      real*8 vturbi
      integer ilsv(nnnl),nlsv
      real*8 ewsv(nnnl),elsv(nnnl)
c     the atomic data creation date
      character(63) atcredate
c
c jg
      real*8 t,xlum

      integer lpri, ncn2
c
c     Not used
      integer javi
c
      data kblnk16/'                '/
c
      javi=np2
      np2=javi
      javi=npfirst(1)
      javi=nplini(1)
      javi=npcon(1)
      javi=npconi(1)
      javi=npilev(1,1)
      javi=npilevi(1)
      javi=npconi2(1)

c
c
      verbose=lpri
      eliml=0.1
      elimh=1.e+5
      elimh=min(elimh,8.9e+4)
c
c     open and prepare the fits file for spectral data
      if(verbose.gt.0) write (lun11,*)'writespectra: opening header',
     $  kmodelname
      knam='xout_spect1.fits'
      call fheader(unit,knam,atcredate,kmodelname,istatus)
      if(istatus.gt.0) call printerror(lun11,istatus)
c
c
c     write extension of parameter values
      if(verbose.gt.0)
     $     write (lun11,*)'writespectra: write parameter list'
      call fparmlist(unit,1,kmodelname,nparms,parname,partype,parval,
     $               parcomm,nloopctl,istatus,lun11)
      if(istatus.gt.0) call printerror(lun11,istatus)
      if(verbose.gt.0)
     $  write (lun11,*)'writespectra: building data tables'


      xlum=parval(10)
      call binemis(lun11,lpri,xlum,
     $       t,vturbi,epi,ncn2,dpthc,
     $       idat1,rdat1,kdat1,nptrs,
     $       npar,npnxt,npfi,
     $       nplin,nlsvn,
     $       eliml,elimh,elum,zrems,zremsz,ilsv,
     $       zrtmp,ewsv,elsv,nlsv)

c
c
c     write the spectral data to the extension
      do mm=1,6
        kunits(mm)=kblnk16
        klabs(mm)=kblnk16
        kform(mm)=kblnk16
        enddo
      klabs(1)='energy          '
      kform(1)='E13.5'
      kunits(1)='eV'
      klabs(2)='incident        '
      kform(2)='E13.5'
      kunits(2)='erg/s/erg'
      klabs(3)='transmitted     '
      kform(3)='E13.5'
      kunits(3)='erg/s/erg'
      klabs(4)='emit_inward     '
      kform(4)='E13.5'
      kunits(4)='erg/s/erg'
      klabs(5)='emit_outward    '
      kform(5)='E13.5'
      kunits(5)='erg/s/erg'
      klabs(6)='scattered       '
      kform(6)='E13.5'
      kunits(6)='erg/s/erg'
c     build extension name
      extname='XSTAR_SPECTRA'
C      if(nloopctl.gt.0) then
C          write(ktmp2,'(i4.4)')nloopctl
C          extname='xstar_spectra_' // ktmp2
C          endif
      if(verbose.gt.0)
     $   write (lun11,*)'writespectra: writing spectral data'

c     append a new empty extension onto the end of the primary array
      status=0
      call ftcrhd(unit,status)
      if(verbose.gt.0)
     $    write (lun11,*)'writespectra: writing header table'

      tfields=6
      nrows=ncn2
      rowlen=0
      do mm=1,6
      tbcol(mm)=0
      enddo

c     write the required header parameters for the ascii table
      status=0
      call ftphtb(unit,rowlen,nrows,tfields,klabs,tbcol,kform,kunits,
     &            extname,status)
      if (status .gt. 0)call printerror(lun11,status)
      status=0
c
c     map each column to a 1-d array before writing to the file
      do kk=1,tfields
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra: building column ',kk
        frow=1
        felem=1
        colnum=kk
        do ll=1,nrows
          rtmp(ll)=zrtmp(kk,ll)
          enddo
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status)
        if (status .gt. 0)call printerror(lun11,status)
        enddo

c     compute checksums
      if(verbose.gt.0) write (lun11,*)'writespectra: writing checksum'
      status=0
      call ftpcks(unit,status)
c     check for any error, and if so print out error messages
      if (status .gt. 0)call printerror(lun11,status)

      if(verbose.gt.0) write (lun11,*)'writespectra: closing file'
      call fitsclose(lun11,unit,istatus)
c
c
      return
      end
      subroutine writespectra2(lun11,lpri,nparms,parname,partype,parval,
     $       parcomm,atcredate,epi,ncn2,dpthc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,elum,tau0,kmodelname,nloopctl)

C
C     Write extension containing the spectra for
C     this particular model.
C
C     Modifications:
C       04/01/1999,WTB: Disabled appending loop control value to
C               extension name due to changes in xstar2xspec design
C
c     author:  T. Bridgman
c
      implicit none
c
      integer ncn2
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     passed parameters
      character(30) kmodelname
      integer nparms, nloopctl, lun11
      character(20) parname(55)
      character(10) partype(55)
      real*8 parval(55)
      character(30) parcomm(55)
c     master data
      integer idat1(nidat1),nptrs(nptt,ndat2)
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     line luminosities
      real*8 elum(3,nnnl)
c     energy bins
      real*8 epi(ncn)
c     the atomic data creation date
      character(63) atcredate
c     continuum optical depths
      real*8 dpthc(2,ncn)
c     line optical depths
      real*8 tau0(2,nnnl)
      real rtmp(ncn)
      character(16) knam,klabs(9),kunits(9),kform(9),kblnk16
      character(30) extname
      integer unit,istatus, nilin, nkdt,nidt,lcon,lrtyp,ltyp,ml,status
      integer nlsvn, ln, ll, lnn, nrdt,mllz
      integer tbcol(9), nrows, rowlen
      integer np2, kk
      integer frow, felem, colnum, tfields, verbose,mm
      real*8 eliml, elimh, elmmtpp,elin
C     Internal work areas
      integer ntptr(nnnl)
      character(10) kion(nnnl)
      character(20) klevl(nnnl),klevu(nnnl)
      integer lpri,lpril
      integer jkk, nlev
      integer nlplmx,nilin2,nlpl,lmm,kltmpn,kltmpo,
     $         llo,lup,llofnd,lupfnd,
     $         k,kl2,lm,kk2,mlpar,mlm,np1i,np1k,np1r
      real*8 elcomp
      character(20) ktmp2,kblnk20
c     Database manipulation quantities
      character(1) kblnk,kdtmp(200)
      integer kltmp(1000)
      real elsv(1000)
      logical done
c
c     Not used
      real*8 javir
      integer javi
c
      data kblnk/' '/
      data kblnk16/'                '/
      data kblnk20/'                    '/
c
      if (nlsvn.le.5) return

      javir=epi(1)
      epi(1)=javir
      javir=dpthc(1,1)
      javi=ncn2
      np2=javi
      javi=npfirst(1)
      javi=nplini(1)
      javi=npcon(1)
      javi=npconi(1)
      javi=npilev(1,1)
      javi=npilevi(1)
      javi=npconi2(1)
c

      verbose=lpri
      eliml=0.1
      elimh=1.0e10
c
c     open and prepare the fits file for spectral data
      if(verbose.gt.0) write (lun11,*)'writespectra2: opening header',
     $  kmodelname
      knam='xout_lines1.fits'
      call fheader(unit,knam,atcredate,kmodelname,istatus)
      if(istatus.gt.0) call printerror(lun11,istatus)

c     write extension of parameter values
      if(verbose.gt.0)
     $ write (lun11,*)'writespectra2: write parameter list'
      call fparmlist(unit,1,kmodelname,nparms,parname,partype,parval,
     $               parcomm,nloopctl,istatus,lun11)
      if(istatus.gt.0) call printerror(lun11,istatus)
      if(verbose.gt.0)
     $  write (lun11,*)'writespectra2: building data tables'
c
c     build spectra data tables
      if (verbose.gt.0) write (lun11,*)' '
      kltmpo=0
      lpril=0
      if (verbose.gt.0)
     $  write (lun11,*)'emission line luminosities (erg/sec/10**38))'
      nlplmx=1000
      eliml=0.1
      elimh=1.0e10
c     find the strongest lines.
      do  lm=1,nlplmx
       kltmp(lm)=0
       enddo
      nlpl=1
      do lnn=1,nlsvn
        ln=lnn
        ml=nplin(ln)
        mlm=ml-1
        call drd(ltyp,lrtyp,lcon,
     $    nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $    nptrs,0,lun11)
        elin=abs(rdat1(np1r))
        if (lrtyp.ne.14) then
          nilin=npar(ml)
          mlm=nilin-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $      nptrs,0,lun11)
          nilin2=idat1(np1i-1+nidt)
          elmmtpp=(elum(2,ln)+elum(1,ln))/2.
          if (verbose.gt.0)
     $         write (lun11,*)lnn,elin,nilin,elmmtpp,ln,ml
          if ((ln.gt.0).and.(ln.lt.nnnl)
     $         .and.(elin.ge.eliml).and.(elin.le.elimh)
     $         .and.(elin.le.8.9e+6)
     $         .and.(elmmtpp.gt.1.e-36)
     $         .and.(nilin2.gt.0).and.(nilin2.le.nni))
     $           then
c
            lmm=0
            elcomp=1.e+10
            do while ((lmm.lt.nlpl).and.(elmmtpp.lt.elcomp))
              lmm=lmm+1
              kl2=kltmp(lmm)
              elcomp=0.
              if (kl2.gt.0)
     $          elcomp=(elum(2,kl2)+elum(1,kl2))/2.
              enddo
c
            if (verbose.gt.0)
     $       write (lun11,8516)ln,elin,elmmtpp
8516        format (1h ,i4,2e12.4)
            kltmpo=ln
            do  k=lmm,min(nlplmx,nlpl)
              if ((lpril.ne.0).and.(kltmp(k).ne.0))
     $         write (lun11,*)'in 557 loop',k,kltmp(k),kltmpo
              kltmpn=kltmp(k)
              kltmp(k)=kltmpo
              kltmpo=kltmpn
              enddo
             nlpl=min(nlplmx,nlpl+1)
            if (verbose.gt.0)
     $       write (lun11,*)'done with 557 loop',lm
            endif
          endif
        enddo
      if (nlpl.gt.0) kltmp(nlpl)=kltmpo
c      nlpl=nlpl-1
      if (verbose.gt.0)
     $    write (lun11,959)
959   format (1x,'index, ion, wavelength, transmitted, reflected')
      kk2=0
      do  kk=1,nlpl
        ln=kltmp(kk)
        if (ln.ne.0) then
          ml=nplin(ln)
          klevl(kk)=kblnk20
          klevu(kk)=kblnk20
          if (ml.ne.0) then
            if (verbose.gt.0)
     $        write (lun11,*)'   ',ln,ml
            mlm=ml-1
            call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $        nptrs,0,lun11)
            llo=idat1(np1i)
            lup=idat1(np1i+1)
            elsv(kk)=abs(rdat1(np1r))
            nilin=npar(ml)
            mlm=nilin-1
            call drd(ltyp,lrtyp,lcon,
     $        nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $        nptrs,0,lun11)
            do mm=1,nkdt
              kdtmp(mm)=kdat1(np1k-1+mm)
              enddo
            do mm=nkdt+1,10
              kdtmp(mm)=kblnk
              enddo
c            nilin=idat(nidt)
            write(kion(kk),'(10a1)')(kdtmp(mm),mm=1,10)
            done=.false.
            jkk=1
            do while (.not.done)
              ml=npfi(13,jkk)
              if (ml.ne.0) then
                if ((npar(ml).eq.nilin).or.(jkk.gt.nni))
     $            done=.true.
                endif
              jkk=jkk+1
              enddo
            if (jkk.gt.nni) ml=0
            if (ml.ne.0) then
              mllz=npar(ml)
              mlpar=npar(ml)
              lupfnd=0
              llofnd=0
              do while ((ml.ne.0).and.(mlpar.eq.mllz)
     $           .and.((llofnd.ne.1).or.(lupfnd.ne.1)))
                mlm=ml-1
                call drd(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $            nptrs,0,lun11)
                nlev=idat1(np1i+nidt-2)
                if (lpri.ne.0)
     $            call dprinto(ltyp,lrtyp,lcon,
     $            nrdt,np1r,nidt,np1i,nkdt,np1k,rdat1,idat1,kdat1,lun11)
                if (lpri.ne.0)
     $            write (lun11,*)nlev,llo,lup,llofnd,lupfnd
                if (nlev.eq.llo) then
                  do mm=1,20
                    if (mm.le.nkdt) then
                        write(ktmp2(mm:mm),'(a1)')kdat1(np1k-1+mm)
                      else
                        write(ktmp2(mm:mm),'(a1)')kblnk
                      endif
                    enddo
                  klevl(kk)=ktmp2
                  llofnd=1
                  endif
                if (nlev.eq.lup) then
                  do mm=1,20
                    if (mm.le.nkdt) then
                        write(ktmp2(mm:mm),'(a1)')kdat1(np1k+mm-1)
                      else
                        write(ktmp2(mm:mm),'(a1)')kblnk
                      endif
                    enddo
                  klevu(kk)=ktmp2
                  lupfnd=1
                  endif
                ml=npnxt(ml)
                mlpar=0
                if (ml.ne.0) mlpar=npar(ml)
                enddo
              kk2=kk2+1
              ntptr(kk2)=ln
              if (verbose.gt.0) then
                write (lun11,*)ml,nilin,npar(ml)
                write (lun11,9955)kk,ln,(kdtmp(mm),mm=1,9),elsv(kk),
     $               elum(1,ln),elum(2,ln)
                write (lun11,*)klevu(kk)
                write (lun11,*)klevl(kk)
                endif
 9955           format (1x,2i8,1x,9a1,3(1pe11.3))
              endif
            endif
          endif
        enddo
c      if (nlpl.le.0) return
c
      nlpl=kk2
      nlpl=max(nlpl,1)
c
      if (verbose.gt.0) then
        do kk=1,nlpl
          write (lun11,*)kk,ntptr(kk)
          enddo
        endif
c

c     write the spectral data to the extension
      do mm=1,9
        kunits(mm)=kblnk16
        klabs(mm)=kblnk16
        kform(mm)=kblnk16
        enddo
      klabs(1)='index           '
      kform(1)='I6'
      kunits(1)='  '
      klabs(2)='ion             '
      kform(2)='A9'
      kunits(2)=' '
      klabs(3)='lower_level     '
      kform(3)='A20'
      kunits(3)='  '
      klabs(4)='upper_level     '
      kform(4)='A20'
      kunits(4)='  '
      klabs(5)='wavelength      '
      kform(5)='F10.2'
      kunits(5)='A'
      klabs(6)='emit_inward     '
      kform(6)='E11.3'
      kunits(6)='erg/s/10**38'
      klabs(7)='emit_outward    '
      kform(7)='E11.3'
      kunits(7)='erg/s/10**38'
      klabs(8)='depth_inward    '
      kform(8)='E11.3'
      kunits(8)='  '
      klabs(9)='depth_outward   '
      kform(9)='E11.3'
      kunits(9)='  '
c     build extension name
      extname='XSTAR_LINES'
C      if(nloopctl.gt.0) then
C          write(ktmp2,'(i4.4)')nloopctl
C          extname='xstar_spectra_' // ktmp2
C          endif
      if(verbose.gt.0)
     $   write (lun11,*)'writespectra2: writing spectral data'

c     append a new empty extension onto the end of the primary array
      status=0
      call ftcrhd(unit,status)
      if(verbose.gt.0)
     $    write (lun11,*)'writespectra2: writing header table'

      tfields=9
      nrows=nlpl
      rowlen=0
      do mm=1,9
        tbcol(mm)=0
        enddo


c     write the required header parameters for the ascii table
      status=0
      call ftphtb(unit,rowlen,nrows,tfields,klabs,tbcol,kform,kunits,
     &            extname,status)
      if (status .gt. 0)call printerror(lun11,status)
      status=0
c
c     map each column to a 1-d array before writing to the file
      kk=1
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk,nrows
        status=0
        call ftpclj(unit,colnum,frow,felem,nrows,ntptr,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=2
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk,nrows
        status=0
        call ftpcls(unit,colnum,frow,felem,nrows,kion,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=3
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcls(unit,colnum,frow,felem,nrows,klevl,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=4
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcls(unit,colnum,frow,felem,nrows,klevu,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=5
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,elsv,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=6
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        do ll=1,nrows
          rtmp(ll)=elum(1,ntptr(ll))
          enddo
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=7
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        do ll=1,nrows
          rtmp(ll)=elum(2,ntptr(ll))
          enddo
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=8
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        do ll=1,nrows
          rtmp(ll)=tau0(1,ntptr(ll))
          enddo
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=9
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        do ll=1,nrows
          rtmp(ll)=tau0(2,ntptr(ll))
          enddo
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status)
        if (status .gt. 0)call printerror(lun11,status)

c     compute checksums
      if(verbose.gt.0) write (lun11,*)'writespectra2: writing checksum'
      status=0
      call ftpcks(unit,status)
c     check for any error, and if so print out error messages
      if (status .gt. 0)call printerror(lun11,status)

      if(verbose.gt.0) write (lun11,*)'writespectra2: closing file'
      call fitsclose(lun11,unit,istatus)
c
c
      return
      end
      subroutine writespectra3(lun11,lpri,nparms,parname,partype,parval,
     $       parcomm,atcredate,epi,ncn2,dpthc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,elum,zrems,zremsz,kmodelname,nloopctl)

C
C     Write extension containing the spectra for
C     this particular model.
C
C     Modifications:
C       04/01/1999,WTB: Disabled appending loop control value to
C               extension name due to changes in xstar2xspec design
C
c     author:  T. Bridgman
c
      implicit none

      integer ncn2

      include './PARAM'
c
c
c     passed parameters
      character(30) kmodelname
      integer nparms, nloopctl, lun11
      character(20) parname(55)
      character(10) partype(55)
      real*8 parval(55)
      character(30) parcomm(55)
c     master data
      integer idat1(nidat1),nptrs(nptt,ndat2)
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     line luminosities
      real*8 elum(3,nnnl)
c     energy bins
      real*8 epi(ncn)
c     the atomic data creation date
      character(63) atcredate
c     continuum lum
      real*8 zrems(4,ncn),zremsz(ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
      real*8 zrtmp(5,ncn)
      real rtmp(ncn)
      character(16) knam,klabs(5),kunits(5),kform(5),kblnk16
c      character(30) ktmp2
      character(30) extname
      integer unit,istatus, kl,lpri
      integer nlsvn, ll, numcon
      integer np2, tbcol(5), nrows, rowlen, kk
      integer frow, felem, colnum, tfields, status, verbose,mm
c
c     Not used
      real*8 javir
      integer javi
c      character(80) javik
c
      data kblnk16/'                '/

      javir=epi(1)
c      epi(1)=javir
      javir=dpthc(1,1)
      javi=ncn2
      javi=np2
c      np2=javi
      javi=npfirst(1)
      javi=nplini(1)
      javi=npcon(1)
      javi=npconi(1)
      javi=npilev(1,1)
      javi=npilevi(1)
      javi=npconi2(1)

      javi=idat1(1)
      javir=rdat1(1)
c      javik=kdat1(1)
      javi=nptrs(1,1)
      javi=npar(1)
      javi=npnxt(1)
      javi=npfi(1,1)
      javi=nplin(1)
      javi=nlsvn
      javir=elum(1,1)

c

      verbose=lpri
c
c     open and prepare the fits file for spectral data
      if(verbose.gt.0) write (lun11,*)'writespectra3: opening header',
     $  kmodelname
      knam='xout_cont1.fits'
      call fheader(unit,knam,atcredate,kmodelname,istatus)
      if(istatus.gt.0) call printerror(lun11,istatus)

c     write extension of parameter values
      if(verbose.gt.0)
     $  write (lun11,*)'writespectra: write parameter list'
      call fparmlist(unit,1,kmodelname,nparms,parname,partype,parval,
     $               parcomm,nloopctl,istatus,lun11)
      if(istatus.gt.0) call printerror(lun11,istatus)
      if(verbose.gt.0)
     $  write (lun11,*)'writespectra: building data tables'

c     build spectra data tables
      numcon=ncn2
      do ll=1,ncn2
        zrtmp(4,ll)=0.
        zrtmp(5,ll)=0.
        enddo
      do kl=1,numcon
         zrtmp(4,kl)=zrtmp(4,kl)+zrems(2,kl)
         zrtmp(5,kl)=zrtmp(5,kl)+zrems(3,kl)
         zrtmp(3,kl)=zremsz(kl)*exp(-dpthc(1,kl))
c         write (lun11,968)kl,epi(kl),zremsz(kl),
c     $          zrtmp1(kl),zrtmp2(kl)
         zrtmp(2,kl)=zremsz(kl)
         zrtmp(1,kl)=epi(kl)
         enddo

c     write the spectral data to the extension
      do mm=1,5
        kunits(mm)=kblnk16
        klabs(mm)=kblnk16
        kform(mm)=kblnk16
        enddo
      klabs(1)='energy          '
      kform(1)='E11.3'
      kunits(1)='eV'
      klabs(2)='incident        '
      kform(2)='E11.3'
      kunits(2)='erg/s/erg'
      klabs(3)='transmitted     '
      kform(3)='E11.3'
      kunits(3)='erg/s/erg'
      klabs(4)='emit_inward     '
      kform(4)='E11.3'
      kunits(4)='erg/s/erg'
      klabs(5)='emit_outward    '
      kform(5)='E11.3'
      kunits(5)='erg/s/erg'
c     build extension name
      extname='XSTAR_SPECTRA'
C      if(nloopctl.gt.0) then
C          write(ktmp2,'(i4.4)')nloopctl
C          extname='xstar_spectra_' // ktmp2
C          endif
      if(verbose.gt.0)
     $  write (lun11,*)'writespectra: writing spectral data'
c      call writespectra(unit,ktmp1,zrtmp,5,999,ncn,
c     $                 klabs,kform,kunits)

c     append a new empty extension onto the end of the primary array
      status=0
      call ftcrhd(unit,status)
      if(verbose.gt.0)
     $   write (lun11,*)'writespectra: writing header table'

      tfields=5
      nrows=ncn2
      rowlen=0
      do mm=1,5
      tbcol(mm)=0
      enddo

c     write the required header parameters for the ascii table
      status=0
      call ftphtb(unit,rowlen,nrows,tfields,klabs,tbcol,kform,kunits,
     &            extname,status)
      if (status .gt. 0)call printerror(lun11,status)
      status=0
c
c     map each column to a 1-d array before writing to the file
      do kk=1,tfields
        if(verbose.gt.0)
     $    write (lun11,*)'writespectra: building column ',kk
        frow=1
        felem=1
        colnum=kk
        do ll=1,nrows
          rtmp(ll)=zrtmp(kk,ll)
          enddo
        if(verbose.gt.0)
     $    write (lun11,*)'writespectra: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rtmp,status)
        if (status .gt. 0)call printerror(lun11,status)
        enddo

c     compute checksums
      if(verbose.gt.0) write (lun11,*)'writespectra: writing checksum'
      call ftpcks(unit,status)
c     check for any error, and if so print out error messages
      if (status .gt. 0)call printerror(lun11,status)

      if(verbose.gt.0) write (lun11,*)'writespectra: closing file'
      call fitsclose(lun11,unit,istatus)
c
      return
      end
      subroutine writespectra4(lun11,lpri,nparms,parname,partype,parval,
     $       parcomm,atcredate,epi,ncn2,dpthc,abel,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,elumab,tauc,kmodelname,nloopctl)
C
C     Write extension containing the spectra for
C     this particular model.
C
C     Modifications:
C       04/01/1999,WTB: Disabled appending loop control value to
C               extension name due to changes in xstar2xspec design
C
c     author:  T. Bridgman
c
      implicit none
c
      integer ncn2

      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     passed parameters
      character(30) kmodelname
      integer nparms, nloopctl, lun11
      character(20) parname(55)
      character(10) partype(55)
      real*8 parval(55)
      character(30) parcomm(55)
c     master data
      integer idat1(nidat1),nptrs(nptt,ndat2)
      real*8 rdat1(nrdat1)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
      real*8 elumab(2,nnml),tauc(2,nnml)
c     energy bins
      real*8 epi(ncn)
c     the atomic data creation date
      character(63) atcredate
c     continuum optical depths
      real*8 dpthc(2,ncn)
      real*8 abel(nl)
      real*8 rlev(10,nd)
      real rsv1(nnml),rsv2(nnml),rsv3(nnml),rsv4(nnml),rsv5(nnml)
      integer ntptr(nnml)
      character(8) kdtmpi(nnml),kdtmp8
      character(20) kdtmpl(nnml),kdtmp20
      character(1) klev(100,nd),kblnk
      character(16) knam,klabs(8),kunits(8),kform(8),kblnk16
      character(30) extname
      integer unit,istatus,kl,nkdt,nidt,lcon,lrtyp,ltyp,ml,nkdti
      integer nlsvn,nrdt
      integer np2, tbcol(8), nrows, rowlen, kk
      integer frow, felem, colnum, tfields, status, verbose,mm
      real*8 eth,xeltp
      integer lpril,lpri,klel,mlel,jk,mt2,mllel,nnz,jkk,klion,mlion,
     $        mlleltp,nlevmx,mltype,mllz,nlev,lk,kkkl,idest1,
     $        kksv,mlpar,mlm,np1k,np1ki,np1i,np1r
c
c     Not used
      real*8 javir
      integer javi
c
      data kblnk/' '/
      data kblnk16/'                '/

      javi=lpri
      javir=epi(1)
c      epi(1)=javir
      javir=dpthc(1,1)
      javi=ncn2
      javi=np2
c      np2=javi
      javi=nlsvn
      javi=nplini(1)
      javi=nplin(1)
      javi=npcon(1)
      javi=npconi(1)
      javi=npilev(1,1)
      javi=npilevi(1)

c

      lpril=lpri
      verbose=lpri
c     open and prepare the fits file for spectral data
      if(verbose.gt.0) write (lun11,*)'writespectra4: opening header',
     $  kmodelname
      knam='xout_rrc1.fits'
      call fheader(unit,knam,atcredate,kmodelname,istatus)
      if(istatus.gt.0) call printerror(lun11,istatus)

c     write extension of parameter values
      if(verbose.gt.0)
     $     write (lun11,*)'writespectra4: write parameter list'
      call fparmlist(unit,1,kmodelname,nparms,parname,partype,parval,
     $               parcomm,nloopctl,istatus,lun11)
      if(istatus.gt.0) call printerror(lun11,istatus)
      if(verbose.gt.0)
     $  write (lun11,*)'writespectra4: building data tables'

c     build spectra data tables
      lpril=verbose
c     print 500 strongest recombination continua
c      write (lun11,*)'recombination continuum luminosities',
c     $  '(erg/sec/10**38))'
c      write (lun11,*)'ion, level, energy (eV), RRC luminosity '
C     lpril is flag for printing debug information
C      initialize line counter
      kksv=0
      jkk=0
C      First look for element data (jk is element index)
        klel=11
        mlel=npfirst(klel)
        jk=0
        do while (mlel.ne.0)
          jk=jk+1
          mt2=mlel-1
          call drd(ltyp,lrtyp,lcon,
     $      nrdt,np1r,nidt,np1i,nkdt,np1k,mt2,
     $      nptrs,0,lun11)
          if (nidt.gt.0) then
            mllel=idat1(np1i-1+nidt)
            xeltp=rdat1(np1r)
            xeltp=abel(mllel)
            nnz=idat1(np1i)
            if (lpril.ne.0)
     $        write (lun11,*)'element:',jk,mlel,mllel,nnz,
     $                    (kdat1(np1k-1+mm),mm=1,nkdt)
C           ignore if the abundance is small
            if (xeltp.lt.1.e-10) then
                jkk=jkk+nnz
              else
c               now step thru ions (jkk is ion index)
                klion=12
                mlion=npfirst(klion)
                jkk=0
                kl=0
                do while ((mlion.ne.0).and.(kl.lt.nnz))
                  jkk=jkk+1
C                 retrieve ion name from kdati
                  mlm=mlion-1
                  call drd(ltyp,lrtyp,lcon,
     $              nrdt,np1r,nidt,np1i,nkdti,np1ki,mlm,
     $              nptrs,0,lun11)
C                 if not accessing the same element, skip to the next element
                  mlleltp=idat1(np1i+nidt-2)
                  if (mlleltp.eq.mllel) then
                    kl=kl+1
                    if (lpril.ne.0)
     $                write (lun11,*)'  ion:',kl,jkk,mlion,mlleltp,
     $                            (kdat1(np1k-1+mm),mm=1,nkdti)
c                   now find level data
c                   step thru types
                    nlevmx=0
                    mltype=13
                    ml=npfi(mltype,jkk)
                    mllz=npar(ml)
                    mlpar=npar(ml)
c                   step thru records of this type
                    do while ((ml.ne.0).and.(mlpar.eq.mllz))
                      mlm=ml-1
                      call drd(ltyp,lrtyp,lcon,
     $                  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $                  nptrs,0,lun11)
                      nlev=idat1(np1i+nidt-2)
                      nlevmx=max(nlevmx,nlev)
                      if ((nlev.gt.0).and.(nlev.le.nd)) then
                        if (lpril.ne.0)
     $                    write (lun11,*)'level quantities:',
     $                    ml,nlev,ltyp,lrtyp,rdat1(np1r),rdat1(np1r+1)
                        do  lk=1,nrdt
                          rlev(lk,nlev)=rdat1(np1r-1+lk)
                          enddo
                        do lk=1,nkdt
                          klev(lk,nlev)=kdat1(np1k-1+lk)
                          enddo
                        do lk=nkdt+1,20
                          klev(lk,nlev)=kblnk
                          enddo
                        endif
                      ml=npnxt(ml)
                      if (ml.ne.0) mlpar=npar(ml)
                      enddo
                    nlev=nlevmx
                    mltype=7
                    ml=npfi(mltype,jkk)
                    mllz=npar(ml)
                    mlpar=npar(ml)
                    do while ((ml.ne.0).and.(mlpar.eq.mllz))
c                     step thru records of this type
                      mlm=ml-1
                      call drd(ltyp,lrtyp,lcon,
     $                  nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $                  nptrs,0,lun11)
                      kkkl=npconi2(ml)
                      idest1=idat1(np1i+nidt-2)
                      if ((kkkl.gt.0).and.(kkkl.le.ndat2)
     $                  .and.((elumab(1,kkkl).gt.1.e-36)
     $                  .or.(elumab(2,kkkl).gt.1.e-36))) then
                        kksv=kksv+1
                        eth=rlev(4,idest1)-rlev(1,idest1)
                        ntptr(kksv)=kkkl
                        rsv1(kksv)=eth
                        rsv2(kksv)=elumab(1,kkkl)
                        rsv3(kksv)=elumab(2,kkkl)
                        rsv4(kksv)=tauc(1,kkkl)
                        rsv5(kksv)=tauc(2,kkkl)
                        do mm=1,nkdti
                          write (kdtmp8(mm:mm),'(a1)')kdat1(np1ki-1+mm)
                          enddo
                        do mm=nkdti+1,8
                          write (kdtmp8(mm:mm),'(a1)')kblnk
                          enddo
                        kdtmpi(kksv)=kdtmp8
                        do mm=1,20
                          write (kdtmp20(mm:mm),'(a1)')klev(mm,idest1)
                          enddo
                        kdtmpl(kksv)=kdtmp20
                        if (lpril.ne.0)
     $                   write (lun11,*)jkk,idest1,
     $                   eth,elumab(1,kkkl),elumab(2,kkkl)
                        if (lpril.ne.0)
     $                   write (lun11,9293)kdtmpi(kksv),
     $                          (klev(lk,idest1),lk=1,20),eth,
     $                          elumab(1,kkkl)
 9293                   format(1x,20a1,20a1,2(1pe11.3))
                        endif
                      ml=npnxt(ml)
                      if (ml.ne.0) mlpar=npar(ml)
                      enddo
                    endif
C                 Go to next ion
                  mlion=npnxt(mlion)
                  enddo
              endif
            endif
          if (mlel.ne.0) mlel=npnxt(mlel)
C         Go to next element
          enddo

c     write the spectral data to the extension
      do mm=1,8
        kunits(mm)=kblnk16
        klabs(mm)=kblnk16
        kform(mm)=kblnk16
        enddo
      klabs(1)='index          '
      kform(1)='I6'
      kunits(1)='  '
      klabs(2)='ion            '
      kform(2)='A9'
      kunits(2)='  '
      klabs(3)='level          '
      kform(3)='A20'
      kunits(3)='  '
      klabs(4)='energy         '
      kform(4)='E11.3'
      kunits(4)='eV'
      klabs(5)='emit_outward    '
      kform(5)='E11.3'
      kunits(5)='erg/s'
      klabs(6)='emit_inward     '
      kform(6)='E11.3'
      kunits(6)='erg/s'
      klabs(7)='depth_outward   '
      kform(7)='E11.3'
      kunits(7)='  '
      klabs(8)='depth_inward    '
      kform(8)='E11.3'
      kunits(8)='  '
c     build extension name
      extname='XSTAR_SPECTRA'
C      if(nloopctl.gt.0) then
C          write(ktmp2,'(i4.4)')nloopctl
C          extname='xstar_spectra_' // ktmp2
C          endif
      if(verbose.gt.0)
     $   write (lun11,*)'writespectra: writing spectral data'

c     append a new empty extension onto the end of the primary array
      status=0
      call ftcrhd(unit,status)
      if(verbose.gt.0)
     $    write (lun11,*)'writespectra: writing header table'

      tfields=8
      nrows=kksv
      if (lpril.ne.0) write (6,*)'nrows=',nrows
      rowlen=0
      do mm=1,8
      tbcol(mm)=0
      enddo

c     write the required header parameters for the ascii table
      status=0
      call ftphtb(unit,rowlen,nrows,tfields,klabs,tbcol,kform,kunits,
     &            extname,status)
      if (status .gt. 0)call printerror(lun11,status)
      status=0
c
c     map each column to a 1-d array before writing to the file
      kk=1
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpclj(unit,colnum,frow,felem,nrows,ntptr,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=2
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcls(unit,colnum,frow,felem,nrows,kdtmpi,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=3
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcls(unit,colnum,frow,felem,nrows,kdtmpl,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=4
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rsv1,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=5
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rsv2,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=6
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rsv3,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=7
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rsv4,status)
        if (status .gt. 0)call printerror(lun11,status)
      kk=8
        if(verbose.gt.0)
     $      write (lun11,*)'writespectra2: building column ',kk
        frow=1
        felem=1
        colnum=kk
        if(verbose.gt.0)
     $     write (lun11,*)'writespectra2: writing column ',kk
        status=0
        call ftpcle(unit,colnum,frow,felem,nrows,rsv5,status)
        if (status .gt. 0)call printerror(lun11,status)

c     compute checksums
      if(verbose.gt.0) write (lun11,*)'writespectra: writing checksum'
      status=0
      call ftpcks(unit,status)
c     check for any error, and if so print out error messages
      if (status .gt. 0)call printerror(lun11,status)

      if(verbose.gt.0) write (lun11,*)'writespectra: closing file'
      call fitsclose(lun11,unit,istatus)
c
c
      return
      end
      subroutine xstarcalc(lpri2,lnerrd,nlimdt,
     $       lpri,lprid,lun11,tinf,vturbi,critf,
     $       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       decomp,ecomp,sxcomp,
     $       tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,nlevs,ncsvn,rates,vsav,idrates,
     $       ntotit,
     $       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,
     $         cllines,clcont,htcomp,clcomp,clbrems,
     $       xilev,bilev,rniss,nmat,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,fline,flinel, t_min, mhd_heat, clsup)

      implicit none
c
c
      include './PARAM'
c
c
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     line luminosities
      real*8 elum(3,nnnl),elumo(3,nnnl)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
      real*8 fline(2,nnnl),flinel(ncn)
c     line optical depths
      real*8 tau0(2,nnnl)
c     energy bins
      real*8 epi(ncn)
c      continuum lum
      real*8 zrems(4,ncn),zremso(4,ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn),brcems(ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     level populations
      real*8 xilev(nnml),bilev(nnml),rniss(nnml)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
      real*8 elumab(2,nnml),elumabo(2,nnml)
      real*8 tauc(2,nnml)
c     ion abundances
      real*8 xii(nni)
c     heating and cooling
      real*8 htt(nni),cll(nni)
      real*8 rrrt(nni),pirt(nni)
      integer nlevs(nni)
c     element abundances
      real*8 ababs(nl)
c     the saved rates
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
c     compton heating data
      real*8 decomp(ncomp,ncomp),ecomp(ncomp),sxcomp(ncomp)
c
c     state variables
      real*8 p,r,t,xpx,delr
c     heating-cooling variables
      real*8 httot,cltot,htcomp,clcomp,clbrems,elcter,cllines,
     $     clcont,hmctot
      real*8 trad
      real*8 cfrac,critf,vturbi,xee,tinf
      integer lcdd,ncn2
c     variables associated with thermal equilibrium solution
      integer nmat,ntotit,lnerrd
c     switches
      integer lprid,lpri,nlimdt
c     strings for atomic data read
      integer nlsvn,ncsvn,lun11,np2,lprisv,lpri2
c
      real*8 t_min
      real*8 mhd_heat
      real*8 clsup
c
      lprisv=lpri
      lpri=0
      if (nlimdt.ne.0) then
        call dsec(lnerrd,nlimdt,
     $       lpri,lprid,lun11,tinf,vturbi,critf,
     $       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       decomp,ecomp,sxcomp,
     $       tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,nlevs,ncsvn,rates,vsav,idrates,
     $       ntotit,
     $       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,
     $         cllines,clcont,htcomp,clcomp,clbrems,
     $       xilev,bilev,rniss,nmat,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,fline,flinel, t_min, mhd_heat)
        endif
c       do ll=1,nnml
c         xilev(ll)=0.
c         enddo
c
       if (lpri2.eq.1) lpri=1
       call func(lpri,lun11,vturbi,critf,
     $       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       decomp,ecomp,sxcomp,
     $       tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,nlevs,ncsvn,rates,vsav,idrates,
     $       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,
     $         cllines,clcont,htcomp,clcomp,clbrems,
     $       xilev,bilev,rniss,nmat,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,fline,flinel, mhd_heat)
       call funcsyn(lpri,lun11,vturbi,critf,
     $       t,trad,r,delr,xee,xpx,ababs,cfrac,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,
     $       zrems,zremso,elumab,elumabo,elum,elumo,
     $       decomp,ecomp,sxcomp,
     $       tau0,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,nlevs,ncsvn,rates,vsav,idrates,
     $       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,
     $         cllines,clcont,htcomp,clcomp,clbrems,
     $       xilev,bilev,rniss,nmat,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,fline,flinel, clsup)
c
c
       lpri=lprisv
c
c
      return
      end
      subroutine xstarsetup(lnerrd,nlimd,
     $       lpri,lprid,lunlog,tinf,critf,
     $       t,tp,r,delr,xee,xpx,ababs,abel,cfrac,xlum,p,lcdd,
     $       epi,ncn2,bremsa,bremsint,atcredate,
     $       decomp,ecomp,sxcomp,
     $       zrems,zremsz,
     $       tau0,dpthc,tauc,
     $       idat1,rdat1,kdat1,nptrs,np2,
     $       npar,npnxt,npfi,npfirst,nlevs,
     $       nplin,nplini,nlsvn,npcon,npconi,npilev,npilevi,
     $       npconi2,ncsvn,rates,vsav,idrates,
     $       ntotit,
     $       xii,rrrt,pirt,htt,cll,httot,cltot,hmctot,elcter,
     $       xilev,bilev,rniss,nmat,elum,
     $       rcem,oplin,rccemis,brcems,opakc,opakscatt,cemab,
     $       cabab,opakab,nlin,elin)
c
c      this routine does many of the setup chores: read in atomic
c        data, set up pointers, zeroing variables.
c      NB: no input parameters are affected
c
      implicit none
c
      include './PARAM'
c
      integer nptmpdim
      parameter (nptmpdim=400000)
c
c     global xstar data
c     master data
      integer idat1(nidat1)
      real*8 rdat1(nrdat1)
      integer nptrs(nptt,ndat2)
      character(1) kdat1(nkdat1)
c     pointers to master data
      integer npar(ndat2),npnxt(ndat2),
     $      npfirst(ntyp)
      integer npfi(ntyp,nni)
c     pointers to line data
      integer nplin(nnnl),nplini(ndat2)
c     pointers to line data
      integer npcon(nnml),npconi2(ndat2),npconi(ndat2)
      integer npilev(nd,nni),npilevi(nnml)
c     line luminosities
      real*8 elum(3,nnnl),elumo(3,nnnl)
c     line emissivities
      real*8 rcem(2,nnnl)
c     line opacities
      real*8 oplin(nnnl)
      real*8 fline(2,nnnl),flinel(ncn)
c     line optical depths
      real*8 tau0(2,nnnl)
c     energy bins
      real*8 epi(ncn)
c      continuum lum
      real*8 zrems(4,ncn),zremso(4,ncn),
     $          zremsz(ncn)
c     continuum optical depths
      real*8 dpthc(2,ncn)
c     continuum flux
      real*8 bremsa(ncn),bremsint(ncn)
c     continuum emissivities
      real*8 rccemis(2,ncn),brcems(ncn)
c     continuum opacities
      real*8 opakc(ncn),opakscatt(ncn)
c     level populations
      real*8 xilev(nnml),bilev(nnml),rniss(nnml)
      real*8 cemab(2,nnml),cabab(nnml),opakab(nnml)
      real*8 elumab(2,nnml),elumabo(2,nnml)
      real*8 tauc(2,nnml)
c     ion abundances
      real*8 xii(nni)
c     heating and cooling
      real*8 htt(nni),cll(nni)
      real*8 rrrt(nni),pirt(nni)
      integer nlevs(nni)
c     element abundances
      real*8 abel(nl),abcosmic(30),ababs(nl)
c     the saved rates
      real*8 rates(4,ndat2)
      integer idrates(2,ndat2)
      real*8 vsav(4,ndat2)
c     compton heating data
      real*8 decomp(ncomp,ncomp),ecomp(ncomp),sxcomp(ncomp)
      integer nlin(nnnl)
      real*8 elin(nnnl)
c     the atomic data creation date
      character(63) atcredate
c
c     local variables
c     state variables
      real*8 p,r,t,xpx,delr,tp
c     heating-cooling variables
      real*8 httot,cltot,htcomp,clcomp,clbrems,elcter,cllines,
     $     clcont,hmctot
c     input parameters
      real*8 xlum,xpxcol
      real*8 cfrac,critf,xee
      integer lcdd,ncn2
c     variables associated with thermal equilibrium solution
      integer nmat,ntotit,lnerrd
c     switches
      integer lprid,lpri
      integer  nlimd,lunlog,lun11,lun25
c     temporary for xwrite
      character(133) tmpst
c     strings for atomic data read
      character(256) datafil3,datafil4,datafile
      logical ex3,ex4
c     temporary integers
      integer ll,mm,ldir,jk,mlm
      integer nlsvn,ncsvn, lenact
c     times
      real*8 tinf,t1s
      integer np1r,np1i,np1k,np2
      integer jlk,j,ml,ltyp,lrtyp,lcon,
     $        nrdt,nidt,nkdt
C     storing info for parameters
      character(20) parname(55)
c
c     Not used
      integer javi
      real*8 javir
c
c     these are the anders and grevesse abundances from the xspec
c       manual
c      data abcosmic/'1.00d+00 ','9.77d-02 ','1.45d-11 ','1.41d-11 ',
c     $ '3.98d-10 ','3.63d-04 ','1.12d-04 ','8.51d-04 ','3.63d-08 ',
c     $ '1.23d-04 ','2.14d-06 ','3.80d-05 ','2.95d-06 ','3.55d-05 ',
c     $ '2.82d-07 ','1.62d-05 ','1.88d-07 ','3.63d-06 ','1.32d-07 ',
c     $ '2.29d-06 ','1.26d-09 ','9.77d-08 ','1.00d-08 ','4.84d-07 ',
c     $ '2.45d-07 ','4.68d-05 ','8.60d-08 ','1.78d-06 ','1.62d-08 ',
c     $ '3.98d-08 '/
c     old values
      data abcosmic/1.00E+000, 1.00D-001,
     $ 1.d-10, 1.d-10, 1.d-10, 3.70D-004, 1.10D-004,
     $ 6.80D-004, 3.98D-008, 2.80D-005, 1.78D-006, 3.50D-005,
     $ 2.45D-006, 3.50D-005, 3.31D-007, 1.60D-005, 3.98D-007,
     $ 4.50D-006, 8.91D-008, 2.10D-006, 1.66D-009, 1.35D-007,
     $ 2.51D-008, 7.08D-007, 2.51D-007, 2.50d-005, 1.26D-007,
     $ 2.00D-006, 3.16D-008, 1.58D-008/
c
c
C     Parameter Names
C
      data parname/'cfrac','temperature',
     $   'lcpres','pressure','density','spectrum',
     $   'spectrum_file','spectun','trad',
     $   'rlrad38','column','rlogxi',
     $   'nsteps','niter','lwrite',
     $   'lprint','lstep',
     $   'habund','heabund',
     $   'liabund','beabund','babund','cabund',
     $   'nabund','oabund','fabund','neabund',
     $   'naabund','mgabund','alabund',
     $   'siabund','pabund','sabund',
     $   'clabund','arabund','kabund',
     $   'caabund','scabund','tiabund',
     $   'vabund','crabund','mnabund ',
     $   'feabund','coabund','niabund',
     $   'cuabund','znabund','emult','taumax','xeemin',
     $   'critf','vturbi','npass','modelname',
     $   'loopcontrol'/

c
      javir=t
      javir=tp
      javi=nlimd
      javi=lcdd
      lcdd=javi
      javir=delr
      javir=cfrac
      javir=xlum
      javir=p
      p=javir
      javir=critf
c
      call remtms(t1s)
c
c     opening message
c     write (lunlog,*)'xstar version 2.2.1bn26'
c
c     Test if atomic database files are available.  Abort if not.
      call getenv('LHEA_DATA', datafile)
c     datafil4 = datafile(1:lenact(datafile))//'/atdb.fits'
      datafil4 = './ptxlib/xstarcomm/atdb.fits'
c     datafil3 = datafile(1:lenact(datafile))//'/coheat.dat'
      datafil3 = './ptxlib/xstarcomm/coheat.dat'
c     inquire(file=datafil3,exist=ex3)
c     inquire(file=datafil4,exist=ex4)
c     if (.not.(ex3 .and. ex4 )) then
c        write(tmpst,*)'xstar: One or more of the Atomic Database files'
c        'Etmpst
c        call xwrite(tmpst,10)
c        write(tmpst,*)'xstar: are missing.'
c        write(lunlog,*)tmpst
c        call xwrite(tmpst,10)
c        write(tmpst,*)'xstar: ',datafil4(1:lenact(datafil4))
c        write(lunlog,*)tmpst
c        call xwrite(tmpst,10)
c        write(tmpst,*)'xstar: ',datafil3(1:lenact(datafil3))
c        write(lunlog,*)tmpst
c        call xwrite(tmpst,10)
c        write(tmpst,*)'Program aborting...'
c        write(lunlog,*)tmpst
c        call xwrite(tmpst,10)
c        close(lunlog)
c        return
c     endif
c
c
c
!      tread=0.
!      trates1=0.
!      thcor=0.
!      trates2=0.
!      theat=0.
!      do kl=1,ntyp
!         tucalc(kl)=0.
!         ncall(kl)=0
!         enddo
c
c
c     read in
c     write (lunlog,*)'Loading Atomic Database...'
c     write (tmpst,*)'Loading Atomic Database...'
      call xwrite(tmpst,10)
c
      call readtbl(nptrs,np1r,np1i,np1k,np2,ndat2,
     &      rdat1,idat1,kdat1,nidat1,datafil4,atcredate,lpri,lunlog)
c
c     print *, lun25
      call getlun(lun25)
      lun25 = 25
c     print *, lun25
c     open(unit=lun25,file=datafil3)
c     open(unit=lun25,file='/lustre/scratch4/turquoise/.mdt3/kinch/h3d_10_a09_10/run_300/ptxlib/xstarcomm/coheat.dat')
c     print *, lun25
c     rewind(lun25)
c     print *, lun25
c     read (lun25,901)
 901  format (1x)
c     do mm=1,ncomp
c       do ll=1,ncomp
c         read (lun25,902)sxcomp(mm),
c    $            ecomp(ll),decomp(mm,ll)
 902      format (9x,e12.4,12x,2(e12.4))
c         enddo
c       enddo
c     close(lun25)
c
c
c     Initialize the database
c     write (lunlog,*)'initializng database...'
c     write (tmpst,*)'initializng database...'
      call xwrite(tmpst,10)
      lpri=0
      call setptrs(lunlog,lpri,
     $ idat1,rdat1,kdat1,nptrs,np2,
     $ npnxt,npfi,npar,npfirst,nplin,
     $ nplini,npcon,npconi,npilev,npilevi,
     $ npconi2,nlevs,nlsvn,ncsvn,abcosmic,abel)
      lpri=0
c
c
c
c     read in parameter values
c     write(lunlog,*)'Atomic Abundances'
c     write(lunlog,*)'Element      Solar    Hydrogen'
      do ll=1,nl
        ababs(ll)=abel(ll)*abcosmic(ll)
c       write(lunlog,9990)parname(17+ll),abel(ll),ababs(ll)
        enddo
c     write(lunlog,*)' '
 9990 format(A10,2(E12.4))
c
c
c     set up and initialize
      tinf=0.31
      call init(lunlog,bremsa,bremsint,tau0,dpthc,tauc,
     $   xii,rrrt,pirt,htt,cll,httot,cltot,
     $   cllines,clcont,htcomp,clcomp,clbrems,
     $   xilev,rcem,oplin,rccemis,brcems,opakc,opakscatt,
     $   cemab,cabab,opakab,elumab,elumabo,elum,elumo,
     $   zrems,zremso,rates,vsav,idrates,fline,flinel)

c
      do jk=1,ncn
        zrems(1,jk)=0.
        zrems(2,jk)=0.
        dpthc(1,jk)=0.
        dpthc(2,jk)=0.
        enddo
      do jk=1,nnnl
        elum(1,jk)=0.
        elum(2,jk)=0.
        tau0(1,jk)=0.
        tau0(2,jk)=0.
        enddo
      do jk=1,nnml
        bilev(jk)=0.
        rniss(jk)=0.
        tauc(1,jk)=0.
        tauc(2,jk)=0.
        enddo
c
      lnerrd=0
      ntotit=0
      lprid=0
      xee=1.21
      elcter=0.
      hmctot=0.
c
          ldir=1
          lun11=lunlog
          call trnfrc(lpri,lun11,ldir,
     $      r,xpxcol,xpx,
     $      epi,ncn2,zremsz,dpthc,opakc,
     $      zrems,bremsa,bremsint,flinel)
c
c
c
      nmat=nds
c
      do jlk=1,nlsvn
         j=jlk
         ml=nplin(j)
         mlm=ml-1
         call drd(ltyp,lrtyp,lcon,
     $     nrdt,np1r,nidt,np1i,nkdt,np1k,mlm,
     $     nptrs,0,lun11)
         elin(jlk)=0.
         nlin(jlk)=0
         if (lrtyp.ne.14) then
           elin(jlk)=rdat1(np1r)
           nlin(jlk)=idat1(np1i+nidt-1)
           endif
         enddo
c
      return
      end
      subroutine xwrite(string,ll)
      implicit none
      integer ll
      character*120 string
c
c     Not used
      integer javi
      javi=ll
c      ll=javi
c
c     write (6,*)string
      return
      end
