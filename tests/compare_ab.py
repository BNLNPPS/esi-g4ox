import numpy as np
import opticks.sysrap.sevt as s

a = s.SEvt.Load("/tmp/GEOM/fakegeom/simg4ox/ALL0/A000/")
b = s.SEvt.Load("/tmp/GEOM/fakegeom/simg4ox/ALL0/B000/")

#np.allclose(a.f.record[:,1:], b.f.record[:,0:-1], rtol=0, atol=1e-5)

assert len(a.f.record) == len(b.f.record)

diff = [i for i, (a, b) in enumerate(zip(a.f.record[:, 1:], b.f.record[:, 0:-1])) if not np.allclose(a, b, rtol=0, atol=1e-5)]
print(diff)

assert diff == [14, 21, 22, 36, 40, 54, 64, 81]

