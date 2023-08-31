, and added equivalent code for wais
  (port 210). - FM
* Updated the INSTALLATION file (I hope people read it 8-). - FM
* Hadn't included last night's typo fix in the archived fileset. - FM
01-27-95
* Added Danny Mayer's (mayer@ljo.dec.com) mods for proper handling of
  no_proxy directives for news (NNTPSERVER) hosts. - FM
* Worked in Ari Luotonen's (luotonen@dxcern.cern.ch) mods for sending
  a "Pragma: no-cache" header for use by a proxy server in conjunction
  with the RELOAD command. - FM
* Fixed to handle file URL's appropriately in conjunction with proxying
  of ftp URL's.  For file URL's on the local host, direct access is
  attempted, with no ftp attempt if that fails, whether or not proxying
  is in effect, and whether or not no_proxy directives have been set.
  Both ftp URL's, and file URL's on a remote host, are sought via ftp
  without a direct access attempt.  File URL's for remote hosts are
  converted to ftp URL's before submission to a proxy server, so no
  special procedure need be implemented to induce the proxy server to
  act on them (You shouldn't continue using file URL's when you intend
  ftp, but Lynx will handle such URL's properly when encountered in
  old documents that use file when ftp is intended.) - FM
* Fixed typo in yesterday's "psychotherapy". - FM
01-26-95
* Escaping of ISINDEX queries needed more tweaks to work properly with
  high value IsoLatin1 characters. - FM
* Applied psychotherapy to the schizophrenic behavior that was exhibited
  when download, upload, print, history or showinfo commands were used
  during displays of each other's menus or temporary files. - FM
* Enabled downloading of links from the history page. - FM
01-25-95
* Added LYK_VERSION command (^V) for showing version of Lynx. - FM
* Include file fixes for TCPWARE. - FM
* Fixed bug in showinfo() which caused core dumps when invoked while
  positioned on any form field which has a "linkname" but a NULL
  "filename". - FM
* Fixed glitch in 'p'rint menu. - FM
01-24-95
* Enabled SOCKSification for any Unix flavor via the SITE-LYDEFS,
  SITE-DEFS and SOCKSLIB definitions in the top