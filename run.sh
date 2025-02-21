git clone https://github.com/BNLNPPS/esi-g4ox
esi-g4ox/
git checkout GenStepsMT
cd ..
cmake -S esi-g4ox -B build
cmake --build build
QCurandState_SPEC=3000000:0:0 /usr/local/opticks/lib/QCurandStateTest
./build/src/simg4ox -g  esi-g4ox/geom/pfrich_min_added_parameters.gdml -m esi-g4ox/run.mac
