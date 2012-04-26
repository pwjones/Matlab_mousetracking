MATLAB="/Applications/MATLAB_R2011b.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/pwjones/.matlab/R2011b"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for convertMovieToBinaryArray" > convertMovieToBinaryArray_mex.mki
echo "CC=$CC" >> convertMovieToBinaryArray_mex.mki
echo "CFLAGS=$CFLAGS" >> convertMovieToBinaryArray_mex.mki
echo "CLIBS=$CLIBS" >> convertMovieToBinaryArray_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> convertMovieToBinaryArray_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> convertMovieToBinaryArray_mex.mki
echo "CXX=$CXX" >> convertMovieToBinaryArray_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> convertMovieToBinaryArray_mex.mki
echo "CXXLIBS=$CXXLIBS" >> convertMovieToBinaryArray_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> convertMovieToBinaryArray_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> convertMovieToBinaryArray_mex.mki
echo "LD=$LD" >> convertMovieToBinaryArray_mex.mki
echo "LDFLAGS=$LDFLAGS" >> convertMovieToBinaryArray_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> convertMovieToBinaryArray_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> convertMovieToBinaryArray_mex.mki
echo "Arch=$Arch" >> convertMovieToBinaryArray_mex.mki
echo OMPFLAGS= >> convertMovieToBinaryArray_mex.mki
echo OMPLINKFLAGS= >> convertMovieToBinaryArray_mex.mki
echo "EMC_COMPILER=unix" >> convertMovieToBinaryArray_mex.mki
echo "EMC_CONFIG=optim" >> convertMovieToBinaryArray_mex.mki
"/Applications/MATLAB_R2011b.app/bin/maci64/gmake" -B -f convertMovieToBinaryArray_mex.mk
