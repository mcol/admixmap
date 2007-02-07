#!/bin/bash
# 
# Automatic builds for different architectures

REVISION="-r 1899"

function walton_serial_gcc {
ssh -C walton.ichec.ie "(cd genepi; svn up $REVISION; ./bootstrap \
&& export CXXFLAGS=\"-march=k8 -msse2 -mfpmath=sse -O3 -fomit-frame-pointer -pipe\" \
&& ./configure && make clean && make)" \
&& scp walton.ichec.ie:genepi/src/admixmap/hapmixmap hapmixmap-dynamic-serial-gcc-x86_64 \
|| echo "Can't build serial version on Walton."
}

function walton_parallel {
ssh -C walton.ichec.ie "(cd genepi; svn up $REVISION; \
./bootstrap \
&& . /etc/profile \
&& . /usr/share/modules/init/sh \
&& module load pathscale \
&& export MPICH=\"/usr/local/mpich/path\" \
&& export CC=\"\$MPICH/bin/mpiCC\" \
&& export CXX=\"\$CC\" \
&& export CXXFLAGS=\"-I\$MPICH/include -L\$MPICH/lib64 \
-L/opt/packages/sprng-2.0/lib \
-mcpu=auto -O3 -OPT:Ofast -fno-math-errno -ffast-math -fexceptions \
\" \
&& ./configure --enable-parallel \
&& make clean && make)" \
&& scp walton.ichec.ie:genepi/src/admixmap/hapmixmap hapmixmap-dynamic-parallel-pathCC-x86_64 \
|| echo "Can't build parallel version with Pathscale on Walton."
}

function hamilton_serial {
ssh -C hamilton.ichec.ie "(cd genepi; svn up $REVISION; ./bootstrap \
&& export CXXFLAGS=\"-mcpu=itanium2\" \
&& ./configure --enable-static-binary && make clean && make)" \
&& scp hamilton.ichec.ie:genepi/src/admixmap/hapmixmap hapmixmap-static-serial-gcc-ia64 \
|| echo "Can't build serial version with gcc on Hamilton."
}

function hamilton_parallel {
ssh -C hamilton.ichec.ie "(cd genepi; svn up $REVISION; ./bootstrap \
&& . /etc/profile \
&& export PATH=\"$PATH:/opt/intel/icc_80/bin\" \
&& export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:/opt/intel/icc_80/lib\" \
&& CC=\"/usr/bin/mpiCC\" \
CXX=\"/usr/bin/mpiCC\" \
CXXFLAGS=\"-L/opt/packages/sprng-2.0/lib -I/usr/local/include \
-mcpu=itanium2 \
\" ./configure --enable-parallel && make clean && make)" \
&& scp hamilton.ichec.ie:genepi/src/admixmap/hapmixmap hapmixmap-dynamic-parallel-icpc-ia64 \
|| echo "Can't build parallel version with ICPC on Hamilton."
}

walton_serial_gcc
walton_parallel
hamilton_serial
hamilton_parallel
