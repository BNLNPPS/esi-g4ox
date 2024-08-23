#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <curand_kernel.h>

#include "SysRap/NP.hh"
#include "SysRap/SEvent.hh"
#include "SysRap/sphoton.h"
#include "SysRap/srng.h"
#include "SysRap/storch.h"

using namespace std;


int main(int argc, char **argv)
{
  unsigned n_torches = 1; // create one torch
  unsigned n_photons = 100;
  bool dump = true;

  NP* gs = NP::Make<float>(n_torches, 6, 4);

  storch* torch = reinterpret_cast<storch*>(gs->bytes());

  for(unsigned torch_id = 0; torch_id < gs->shape[0]; torch_id++)
  {
    cout << "torch_id: " << torch_id << endl;
    // Set all torch parameters
    storch::FillGenstep(torch[torch_id], torch_id, n_photons, dump);
  }

  gs->dump();
  cout << torch->desc() << endl;

  NP* ph = NP::Make<float>(n_photons, 4, 4);

  const quad6* qtorch = reinterpret_cast<quad6*>(gs->bytes());
  sphoton* photons = reinterpret_cast<sphoton*>(ph->bytes());
  int unused = -1;

  srng rng(12345);

  for(unsigned photon_id = 0; photon_id < n_photons; photon_id++)
  {
    storch::generate(photons[photon_id], rng, *qtorch, unused, unused);
  }

  ph->set_meta({"my photon meta info 1", "additional meta info"});
  ph->dump();
  ph->save("out/photons.npy");

  return EXIT_SUCCESS;
}
