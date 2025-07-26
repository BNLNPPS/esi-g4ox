git clone https://github.com/BNLNPPS/esi-g4ox
cd esi-g4ox/
git checkout PerformanceAnalysis
cd ..
QCurandState_SPEC=3000000:0:0 /usr/local/opticks/lib/QCurandStateTest
export OPTICKS_EVENT_MODE=Minimal
export OPTICKS_MAX_PHOTON=M60
export OPTICKS_MAX_SLOT=M60
python scripts/cone_opticks_intersection_modify.py
cmake --build /usr/local/eic-opticks/build --parallel --target install
cmake -S esi-g4ox -B build
cmake --build build
./build/src/simg4ox -g  esi-g4ox/geom/pfrich_min_FINAL.gdml -m esi-g4ox/run.mac
