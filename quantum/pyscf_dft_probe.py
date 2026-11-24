import time, sys, resource, numpy as np
from pyscf import lib
from pyscf.pbc import gto, dft
from pyscf.pbc.scf.addons import smearing_
lib.num_threads(2)
k = int(sys.argv[1]); basis = sys.argv[2]
t=time.time()
a=4.05/2
cell = gto.Cell(); cell.atom='Al 0 0 0'; cell.a=[[0,a,a],[a,0,a],[a,a,0]]; cell.unit='A'
cell.basis=basis; cell.pseudo='gth-pbe'; cell.verbose=0; cell.build()
kpts = cell.make_kpts([k,k,k])
mf = smearing_(dft.KRKS(cell,kpts).set(xc='pbe', conv_tol=1e-7), sigma=0.01, method='fermi')
print(f"k={k} basis={basis} nao={cell.nao} mesh={cell.mesh} nk={len(kpts)}", flush=True)
e = mf.kernel(); t1=time.time()-t
print(f"  E={e:.6f} Ha conv={mf.converged} scf={t1:.1f}s maxrss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.2f}GB", flush=True)
g = mf.nuc_grad_method().kernel(); t2=time.time()-t-t1
print(f"  grad={t2:.1f}s  |g|={np.abs(g).max():.2e}", flush=True)
from pyscf.pbc.grad import krks_stress
s = krks_stress.kernel(mf.nuc_grad_method()); t3=time.time()-t-t1-t2
print(f"  stress={t3:.1f}s  diag={np.diag(s)}", flush=True)
