### Quick Start

To get the source code for `esi-g4ox`, download or clone the repository from [GitHub](https://github.com/BNLNPPS/esi-g4ox). If you have the code in your `$HOME` directory and prefer not to install all external dependencies manually, you can use [`esi-shell`](https://github.com/BNLNPPS/esi-shell) for development. The following command mounts your `$HOME` inside the container:

```bash
esi-shell -- -v $HOME -e HOME=$HOME -w $HOME -u $(id -u ${USER}):$(id -g ${USER})
```

Once inside the container, run the following commands:

```bash
cd $HOME
cmake -S esi-g4ox -B build
cmake --build build
```

Similarly, to start a Geant4-based simulation of optical photons, run:

```bash
./build/src/simg4ox
```

Before runnning `simg4ox` it is recommended to set the following environment variables:

```
export TMP=/tmp/myname
export GEOM=mygeom
export OPTICKS_EVENT_MODE=DebugLite
QCurandState_SPEC=3:0:0 /usr/local/opticks/lib/QCurandStateTest
opticks_running_mode=SRM_TORCH
export storch_FillGenstep_pos=0,-270,0
export storch_FillGenstep_type=rectangle
export storch_FillGenstep_zenith=-20,20
export storch_FillGenstep_azimuth=-20,20
export OPTICKS_RUNNING_MODE=${OPTICKS_RUNNING_MODE:-$opticks_running_mode}
export OPTICKS_NUM_PHOTON=100
unset export SGenerate__GeneratePhotons_RNG_PRECOOKED=NO
```


## Visualization

Plot any volume serialized 

```
scripts/plot-csg.py ../out/csg/CSGFoundry/SSim/scene/meshmerge/0/
```
