HD_COMPONENT_NAME	= nicer

HD_COMPONENT_VERS	=

HD_SUBDIRS		= nimergetime nicer-multi-maketime

HD_TEST_SUBDIRS		= $(HD_SUBDIRS) .

include ${HD_STD_MAKEFILE}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 infile,s,a,"file.mkf",,,"Name of FITS mkf file or (@files.lis) "
outfile,s,a,"maketime.gti",,,"Name of output FITS file"
expr,s,a,"",,,"Selection Expression"
name,s,h,"NAME",,,"Column containing HK parameter names"
value,s,h,"VALUE",,,"Column containing HK parameter values"
time,s,h,"TIME",,,"Column containing HK parameter times"
start,s,h,"START",,,"Column containing GTI start times"
stop,s,h,"STOP",,,"Column containing GTI stop times"
compact,b,h,NO,,,"Flag, yes if HK format is compact"
copykw,b,h,yes,,,"Copy all other keywords?"
histkw,b,h,yes,,,"Print history keyword flag"
prefr,r,h,0.5,-1.0,1.0,"Pre-Time Interval factor [0,1]"
postfr,r,h,0.5,-1.0,1.0,"Post-Time Interval factor [0,1]"
premax,r,h,-1.0,-1.0,,"Maximum time amount to extend pre-time interval (or -1) "
postmax,r,h,-1.0,-1.0,,"Maximum time amount to extend post-time interval (or -1) "
mingti,r,h,0,0,,"Minimum GTI entry size in TIME units (or 0 to accept all GTIs) "
emptygti,s,h,"APPLY",,,"Whether to IGNORE empty GTIs or APPLY an empty value "
extname,s,h,"STDGTI",,,"Extension name of output GTI extension "
clobber,b,h,no,,,"Delete outfile if it exists? (override with !filename)"
chatter, i, h, 2, 0, 5, "Verbosity level "
mode,s,h,"ql",,,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     