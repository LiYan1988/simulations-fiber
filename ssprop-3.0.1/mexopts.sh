#
# gccopts.sh	Shell script for configuring MEX-file creation script,
#               mex.  These options were tested with gcc 2.95.2.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  No spaces are allowed
#               around the '=' in the variable assignment.
#
#               Note: only the gcc side of this script was tested.
#               The FORTRAN variables are lifted directly from
#               mexopts.sh; use that file for compiling FORTRAN
#               MEX-files.
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_MEX_OPT: Template Options file for building gcc MEX-files
#
# Copyright 1984-2000 The MathWorks, Inc.
# $Revision: 1.43 $  $Date: 2002/06/10 18:56:23 $
#----------------------------------------------------------------------------
#
    TMW_ROOT="$MATLAB"
    MFLAGS=''
    if [ "$ENTRYPOINT" = "mexLibrary" ]; then
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmatlb -lmat -lmwservices -lut -lm"
    else  
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat -lm"
    fi
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The cmex script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
#
# Determine the location of the GCC libraries
#
	    GCC_LIBDIR=`gcc -v 2>&1 | awk '/.*Reading specs.*/ {print substr($4,0,length($4)-6)}'`
            ;;
        alpha)
#----------------------------------------------------------------------------
            CC='gcc'
            CFLAGS='-mieee -pthread'
            CLIBS="$MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CXX='gcc'
            CXXFLAGS='-mieee -pthread'
            CXXLIBS="$MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-pthread -shared -Wl,-expect_unresolved,'*',-hidden,-exported_symbol,$ENTRYPOINT,-exported_symbol,mexVersion,-exported_symbol,'__*'"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        hpux)
#----------------------------------------------------------------------------
            CC='gcc'
            CFLAGS='-fPIC -D_POSIX_C_SOURCE=199506L -D_HPUX_SOURCE'
            CLIBS="$MLIBS -lm -L$GCC_LIBDIR -lgcc"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CXX='gcc'
            CXXFLAGS='-fPIC -D_POSIX_C_SOURCE=199506L -D_HPUX_SOURCE'
            CXXLIBS="$MLIBS -lm -L$GCC_LIBDIR -lgcc"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            LDCXX='ld'
            LDCXXFLAGS="-b +e $ENTRYPOINT +e mexVersion"
            LDCXXOPTIMFLAGS=''
            LDCXXDEBUGFLAGS=''
#
            LD='ld'
            LDFLAGS="-b +e $ENTRYPOINT +e mexVersion"
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        hp700)
#----------------------------------------------------------------------------
            CC='gcc'
#           Remove -mpa-risc-1-0 from CFLAGS if you wish to optimize
#           for target machine
            CFLAGS='-fPIC -D_HPUX_SOURCE -mpa-risc-1-0'
            CLIBS="$MLIBS -L$GCC_LIBDIR -lgcc"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CXX='gcc'
#           Remove -mpa-risc-1-0 from CFLAGS if you wish to optimize
#           for target machine
            CXXFLAGS='-fPIC -D_HPUX_SOURCE -mpa-risc-1-0'
            CXXLIBS="$MLIBS -L$GCC_LIBDIR -lgcc"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            LDCXX="$COMPILER"
            LDCXXFLAGS="-b -Wl,+e,$ENTRYPOINT,+e,mexVersion,+e,_shlInit,+e,errno"
            LDCXXOPTIMFLAGS=''
            LDCXXDEBUGFLAGS=''
#
            LD='ld'
            LDFLAGS="-b +e $ENTRYPOINT +e mexVersion +e errno"
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        ibm_rs)
#----------------------------------------------------------------------------
            CC='gcc'
            CFLAGS='-D_THREAD_SAFE -D_ALL_SOURCE'
            CLIBS="$MLIBS -lm "
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CXX='gcc'
            CXXFLAGS='-D_THREAD_SAFE -D_ALL_SOURCE'
            CXXLIBS="$MLIBS -lm "
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            LDCXX="$COMPILER"
            LDCXXFLAGS="-shared -Wl,-bE:$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDCXXOPTIMFLAGS='-O -Wl,-s'
            LDCXXDEBUGFLAGS=''
#
            LD="$COMPILER"
            LDFLAGS="-shared -Wl,-bE:$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O -Wl,-s'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        glnx86)
#----------------------------------------------------------------------------
            RPATH="-Wl,--rpath-link,$TMW_ROOT/extern/lib/$Arch,--rpath-link,$TMW_ROOT/bin/$Arch"
            CC='gcc'
            CFLAGS='-fPIC -ansi -D_GNU_SOURCE -pthread'
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O3 -malign-double -fno-exceptions -march=pentium4 -mfpmath=sse -ffast-math'
            CDEBUGFLAGS='-g'
#
            CXX='gcc'
            CXXFLAGS='-fPIC -ansi -D_GNU_SOURCE -pthread'
            CXXLIBS="$RPATH $MLIBS -lm"
            CXXOPTIMFLAGS='-O3 -malign-double -fno-exceptions -march=pentium4 -mfpmath=sse -ffast-math'
            CXXDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O3'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sgi)
#----------------------------------------------------------------------------
            CC='gcc'
            CFLAGS='-D_POSIX_C_SOURCE=199506L -D__EXTENSIONS__ -D_XOPEN_SOURCE -D_XOPEN_SOURCE_EXTENDED'
            CLIBS="-dont_warn_unused $MLIBS -lm -L$GCC_LIBDIR -lgcc"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CXX='gcc'
            CXXFLAGS='-D_POSIX_C_SOURCE=199506L -D__EXTENSIONS__ -D_XOPEN_SOURCE -D_XOPEN_SOURCE_EXTENDED'
            CXXLIBS="-dont_warn_unused $MLIBS -lm -L$GCC_LIBDIR -lgcc"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            LD='ld'
            LDFLAGS="-n32 -shared -exported_symbol $ENTRYPOINT -exported_symbol mexVersion"
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sol2)
#----------------------------------------------------------------------------
            CC='gcc'
            CFLAGS='-fPIC'
            CLIBS="$MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CXX='gcc'
            CXXFLAGS='-fPIC'
            CXXLIBS="$MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-shared -Wl,-M,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        mac)
#----------------------------------------------------------------------------
            CC='cc'
            CFLAGS='-fno-common -traditional-cpp'

            CLIBS="$MLIBS"

            COPTIMFLAGS='-O3 -DNDEBUG'
            CDEBUGFLAGS='-g'

            if [ -f /usr/bin/g++2 ]; then 
                CXX=g++2
            else
                CXX=c++
            fi
            CXXFLAGS='-fno-common -traditional-cpp'
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O3 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            LD="$CC"
            LDFLAGS="-bundle -Wl,-flat_namespace -undefined suppress"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS='nmedit -s $TMW_ROOT/extern/lib/$Arch/$MAPFILE $mex_file'
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
