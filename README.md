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

To run a test that generates a NumPy file (`out/photons.npy`) with the simulated photons, execute:

```bash
./build/src/simtox
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
```


## Visualization

Plot any volume serialized 

```
scripts/plot-csg.py ../out/csg/CSGFoundry/SSim/scene/meshmerge/0/
```

## GPU code profiling

nsys software can be used to profile GPU code. By adding 

```
nsys profile
```
An examples as follows:

```
nsys profile ./build/src/simg4ox -g esi-g4ox/geom/opticks_raindrop.gdml -m esi-g4ox/run.mac
```

The resulting file has nsys-rep extension. This can be opened either by nsys-ui or can be exported with "nsys stats" command into an sqlite database.

In order to achieve a more detailed profiling one can use the nsight compute framework. The user needs to download nsight compute, accept the license agreement. 

If one is using docker the argument needs to be added.

An example of running the nsight compute profiling command:

```
/usr/local/NVIDIA-Nsight-Compute/ncu --set full --target-processes all -o report --import-source on ./build/src/simg4ox -g  esi-g4ox/geom/pfrich_min_added_parameters.gdml -m esi-g4ox/run.mac
```

Afterwards the resulting profile file can opened in the following GUI:

```
/usr/local/NVIDIA-Nsight-Compute/ncu-ui
```
