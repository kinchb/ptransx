that is error
  prone (no change on VMS; editable pseudo-headers are included for
  informational purposes, but the real headers are those of the
  account running Lynx - FM).
* Don't hardcode the Unix mv command; use MV_PATH define in userdefs.h.
* Moved the Unix COPY_COMAND define to userdefs.h as COPY_PATH.
* Completed the Unix ZIP support.
02-10-95
* Added code for converting the headers of downloaded binaries on VMS
  to indicate FIXED 512 record format.  See the documentation in
  userdefs.h, lynx.cfg and FIXED512.COM for more information. - FM
02-08-95
* Fixed glitch in HTInit.c which mis-casted the default MIME type for
  files with a .sh extension, causing the EXEC_SCRIPT function to be
  disfunctional. - FM
* I think I finally have Lynx showing "(p# of N)" properly in the title
  lines (we'll see if perfection has really been achieved 8-). - FM
* Extended directory browsing on VMS to the -homepage specification,
  if included on the command line, e.g., lynx -homepage=sys$login
  will start up Lynx with the default startfile, but the 'm'ain menu
  command will yield listings of the home directory. - FM
* Added -fileversions switch on VMS for including all versions of files
  in directory browser listings (otherwise, only the highest version is
  listed, with no version numbers displayed). - FM
* Modified LYEdit.c to use fopen(filename,"a") instead of access(filename,2)
  to verify write access for editing files (checking "append" access appears
  to be a more reliable way to do it across platforms/flavors). - FM
02-07-95
* Implemented directory browsing for VMS. - FM
* Fixed printer command handling to actually prompt for and use a second
  filename argument if two "%s" strings are in the command map. - FM
* Added anti-Unix-shell-spoofing code for all of the filename argument
  handling. - FM
* Fixed glitch in getenv(DISPLAY) calls on VMS.
02-04-95
* Fixed a glitch in BASE support.  Should now work properly when both
  standard and Lynx-specific hrefs are in the document. - FM
02-03-95
* Added support for posting to newsgroups from Lynx on VMS via the
  ANU-NEWS software. - FM
* Numerous href parsing enhancements. - FM
* Enhancements and bug fixes of news displays. - FM
* Enabled display of Linknames with angle brackets (as for news URLs)
  in the showinfo page. - FM
* Ensured that a file URL which is really an ftp URL will yield an
  "FTP is Disabled" statusline message when it's disabled. - FM
02-01-95
* Defined out the LYK_VERSION code and made Control-V a dead key again,
  since including the Lynx name and version strings in the showinfo ('=')
  and 'o'ptions displays is adequate, and people on different platforms
  or flavors associate Control-V with other functions. - FM
01-31-95
* Beefed up the BASE support. - FM
* Made all WWWLib cover pages and menus HTTP/1.0 com