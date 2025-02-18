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
./build/src/simg4ox -g  esi-g4ox/geom/pfrich_min_added_parameters.gdml -m esi-g4ox/run.mac
```

Before runnning `simg4ox` it is recommended to set the following environment variables:

```
export TMP=/tmp/myname
export GEOM=mygeom
export OPTICKS_EVENT_MODE=DebugLite
QCurandState_SPEC=3:0:0 /usr/local/opticks/lib/QCurandStateTest
```

If the number of max photons simulated in Opticks needs to be increased use the following commands:

```
export OPTICKS_MAX_PHOTON=M60
QCurandState_SPEC=60:0:0 /usr/local/opticks/lib/QCurandStateTest
```
