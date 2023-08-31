				03-Jun-1994

	The bug associated with using non-blocking socket handling with
TCPWare's UCX emulation has been traced by Bernie Votz of Process Software
to their TCPDRIVER.EXE, and has been fixed.  The kit is on FTP.PROCESS.COM
and is called DRIVERS_V405B.INC (it is only valid for TCPware for OpenVMS
Version 4.0-5, which is the current version at the time of this announcement).
The fix will be included in the next TCPWare upgrade.  If you have problems,
contact:

	Bernie Volz
	Process Software Corporation
	VOLZ@PROCESS.COM

	If you have an old version of TCPWare, and no current Process
Software support, modify [.WWW.Library.VMS]libmake.com to include NO_IOCTL in
the $ cc := cc/define=(...) list.  You will not be able to interrupt stalled
connects (i.e., you'll have to wait for those to time out), but you still will
we able to interrupt long or stalled document transfers via the 'z' command.

	To build for TCPWare, enter  @BUILD TCPWARE   or just   @BUILD
and then answer the prompt with the appropriate number for that TCPIP 
package. 

				Fote

=========================================================================
 Foteos Macrides           Worcester Foundation for Experimental Biology
 MACRIDES@SCI.WFEB.EDU     222 Maple Avenue, Shrewsbury, MA 01545
=========================================================================
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             