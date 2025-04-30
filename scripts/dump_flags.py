from opticks.ana.base import PhotonCodeFlags

def ffs(x):
    """return 1-based index of least-significant bit in x (0 if x==0)"""
    return (x & -x).bit_length()

pcf = PhotonCodeFlags()
print(pcf)

