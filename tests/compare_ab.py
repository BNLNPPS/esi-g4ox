import numpy as np
import opticks.sysrap.sevt as s

a = s.SEvt.Load("/tmp/GEOM/fakegeom/simg4ox/ALL0_none/A000/")
b = s.SEvt.Load("/tmp/GEOM/fakegeom/simg4ox/ALL0_none/B000/f000")

assert a.f.record.shape == b.f.record.shape

diff = [i for i, (a, b) in enumerate(zip(a.f.record[:, 1:], b.f.record[:, 0:-1])) if not np.allclose(a, b, rtol=0, atol=1e-5)]
print(diff)

assert diff == [14, 22, 32, 34, 40, 81, 85]
