	basic_machine=i860-stratus
		os=-sysv4
		;;
	strongarm-* | thumb-*)
		basic_machine=arm-`echo "$basic_machine" | sed 's/^[^-]*-//'`
		;;
	sun2)
		basic_machine=m68000-sun
		;;
	sun2os3)
		basic_machine=m68000-sun
		os=-sunos3
		;;
	sun2os4)
		basic_machine=m68000-sun
		os=-sunos4
		;;
	sun3os3)
		basic_machine=m68k-sun
		os=-sunos3
		;;
	sun3os4)
		basic_machine=m68k-sun
		os=-sunos4
		;;
	sun4os3)
		basic_machine=sparc-sun
		os=-sunos3
		;;
	sun4os4)
		basic_machine=sparc-sun
		os=-sunos4
		;;
	sun4sol2)
		basic_machine=sparc-sun
		os=-solaris2
		;;
	sun3 | sun3-*)
		basic_machine=m68k-sun
		;;
	sun4)
		basic_machine=sparc-sun
		;;
	sun386 | sun386i | roadrunner)
		basic_machine=i386-sun
		;;
	sv1)
		basic_machine=sv1-cray
		os=-unicos
		;;
	symmetry)
		basic_machine=i386-sequent
		os=-dynix
		;;
	t3e)
		basic_machine=alphaev5-cray
		os=-unicos
		;;
	t90)
		basic_machine=t90-cray
		os=-unicos
		;;
	tile*)
		basic_machine=$basic_machine-unknown
		os=-linux-gnu
		;;
	tx39)
		basic_machine=mipstx39-unknown
		;;
	tx39el)
		basic_machine=mipstx39el-unknown
		;;
	toad1)
		basic_machine=pdp10-xkl
		os=-tops20
		;;
	tower | tower-32)
		basic_machine=m68k-ncr
		;;
	tpf)
		basic_machine=s390x-ibm
		os=-tpf
		;;
	udi29k)
		basic_machine=a29k-amd
		os=-udi
		;;
	ultra3)
		basic_machine=a29k-nyu
		os=-sym1
		;;
	v810 | necv810)
		basic_machine=v810-nec
		os=-none
		;;
	vaxv)
		basic_machine=vax-dec
		os=-sysv
		;;
	vms)
		basic_machine=vax-dec
		os=-vms
		;;
	vpp*|vx|vx-*)
		basic_machine=f301-fujitsu
		;;
	vxworks960)
		basic_machine=i960-wrs
		os=-vxworks
		;;
	vxworks68)
		basic_machine=m68k-wrs
		os=-vxworks
		;;
	vxworks29k)
		basic_machine=a29k-wrs
		os=-vxworks
		;;
	w65*)
		basic_machine=w65-wdc
		os=-none
		;;
	w89k-*)
		basic_machine=hppa1.1-winbond
		os=-proelf
		;;
	x64)
		basic_machine=x86_64-pc
		;;
	xbox)
		basic_machine=i686-pc
		os=-mingw32
		;;
	xps | xps100)
		basic_machine=xps100-honeywell
		;;
	xscale-* | xscalee[bl]-*)
		basic_machine=`echo "$basic_machine" | sed 's/^xscale/arm/'`
		;;
	ymp)
		basic_machine=ymp-cray
		os=-unicos
		;;
	none)
		basic_machine=none-none
		os=-none
		;;

# Here we handle the default manufacturer of certain CPU types.  It is in
# some cases the only manufacturer, in others, it is the most popular.
	w89k)
		basic_machine=hppa1.1-winbond
		;;
	op50n)
		basic_machine=hppa1.1-oki
		;;
	op60c)
		basic_machine=hppa1.1-oki
		;;
	romp)
		basic_machine=romp-ibm
		;;
	mmix)
		basic_machine=mmix-knuth
		