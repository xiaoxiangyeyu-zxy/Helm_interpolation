from scipy.io import FortranFile

abar = 4.0
zbar = 2.0

fi = FortranFile('helm_table.bdat', 'r')
f = fi.read_reals(dtype=float)
fd = fi.read_reals(dtype=float)
ft = fi.read_reals(dtype=float)
fdd = fi.read_reals(dtype=float)
ftt = fi.read_reals(dtype=float)
fdt = fi.read_reals(dtype=float)
fddt = fi.read_reals(dtype=float)
fdtt = fi.read_reals(dtype=float)
fddtt = fi.read_reals(dtype=float)
dpdf = fi.read_reals(dtype=float)
dpdfd = fi.read_reals(dtype=float)
dpdft = fi.read_reals(dtype=float)
dpdfdt = fi.read_reals(dtype=float)
ef = fi.read_reals(dtype=float)
efd = fi.read_reals(dtype=float)
eft = fi.read_reals(dtype=float)
efdt = fi.read_reals(dtype=float)
xf = fi.read_reals(dtype=float)
xfd = fi.read_reals(dtype=float)
xft = fi.read_reals(dtype=float)
xfdt = fi.read_reals(dtype=float)
# print(f)
# print(fd)
# print(ft)
# print(len(xfdt))


