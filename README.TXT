Assuming that you've setup your CUDA environment all ready...

	You may need to adjust:
		NVIDIA_SDK -> The location of the nvidia sdk.
		LIBS -> The location of cudart.lib and the cufft.lib
		CUFLAGS -> check the following flags
			-ccbin --> set to the location of your compiler.
			-IC:/CUDA/v3.2/include 			<-- set to the location of cuda header files
			-IC:/CUDA/v3.2/include/cudart	<-- set to the location of cudart header files
	

#########################################################################
# Installation (Linux, CUDA >= 3.0 (might work on older releases, too)) #
#########################################################################

Rename makefile.linux to Makefile

Adjust the path to your CUDA installation as desribed above in "Makefile" and run 'make'.



#########################################################################
# Installation  Windows
#########################################################################

A GNU compatible version of make is also required as
the Makefile is not compatible with nmake. GNU Make for Win32 can be dowloaded
from http://gnuwin32.sourceforge.net/packages/make.htm.


Rename makefile.win64 or makefile.win32 to Makefile 

Adjust the path to your CUDA installation as desribed above in "Makefile" and run 'make'.

Run the Visual Studio command prompt for 32 bit builds or the  x64 Command Prompt for 64 bit builds and change into the 
subdirectory where the source files are and run make.

Note: You can not compile using a 64 bit cuda installation on a 32 bit platform.