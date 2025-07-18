git clone https://github.com/BNLNPPS/esi-g4ox
cd esi-g4ox/
git checkout PrincetonEnd
poetry install
cd ..
python modify_zsphere_intersect.py
cmake -S esi-g4ox -B build
cmake --build build
export OPTICKS_EVENT_MODE=Minimal
export OPTICKS_MAX_PHOTON=M60
export OPTICKS_MAX_SLOT=M60
QCurandState_SPEC=60:0:0 /usr/local/opticks/lib/QCurandStateTest
cmake --build /usr/local/eic-opticks/build --parallel --target install
cmake --build build
./build/src/simtox
./build/src/simg4ox -g esi-g4ox/geom/princeton_pmt_spherical.gdml -m esi-g4ox/vis.mac -i
