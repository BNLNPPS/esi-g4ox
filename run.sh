git clone https://github.com/BNLNPPS/esi-g4ox
cd esi-g4ox/
git checkout GenSteps
cd ..
cmake -S esi-g4ox -B build
cmake --build build
export OPTICKS_EVENT_MODE=Minimal
export OPTICKS_MAX_PHOTON=M60
QCurandState_SPEC=60:0:0 /usr/local/opticks/lib/QCurandStateTest
./build/src/simg4ox -g  esi-g4ox/geom/pfrich_min_added_parameters.gdml -m esi-g4ox/run.mac
