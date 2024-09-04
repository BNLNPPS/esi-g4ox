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
