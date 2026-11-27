"""Finite-difference check of the pyscf worker backend: is the energy consistent with its forces and stress?

Runs quantum/worker.py pyscf on fcc Al (1-atom primitive cell and 4-atom cubic cell, gth-szv, 2x2x2 k-points) and
compares the trace of the analytic virial with a finite strain, and one analytic force with a finite displacement.
Takes about 8 minutes on 2 cores. Output: results/pyscf_bridge_fd_log.txt.
"""
import json, os, subprocess, sys, time
WORKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "worker.py")
w = subprocess.Popen([sys.executable, WORKER, "pyscf", "--kpts", "2,2,2", "--basis", "gth-szv", "--threads", "2"],
                     stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
def rq(o):
    w.stdin.write(json.dumps(o) + "\n"); w.stdin.flush(); return json.loads(w.stdout.readline())
print(rq({"op": "hello", "protocol": 1, "elements": ["Al"]}))
a = 4.05 / 2
def ev(cell, pos, i=[0]):
    i[0] += 1; t = time.time()
    r = rq({"op": "evaluate", "id": i[0], "cell": cell, "positions": pos, "species": [0]*len(pos), "symbols": ["Al"]*len(pos)})
    r["t"] = time.time() - t; return r
def strain(c, e):  # r -> (I+e) r, cell rows
    F = [[(1.0 if i == j else 0.0) + e[i][j] for j in range(3)] for i in range(3)]
    return [[sum(F[i][k] * row[k] for k in range(3)) for i in range(3)] for row in c]
cell = [[0, a, a], [a, 0, a], [a, a, 0]]
r0 = ev(cell, [[0, 0, 0]])
V = abs(cell[0][1]*(cell[1][2]*cell[2][0]-cell[1][0]*cell[2][2]) + cell[0][2]*(cell[1][0]*cell[2][1]-cell[1][1]*cell[2][0]))
print("E0 = %.6f eV  stress diag %s  t=%.1fs" % (r0["energy"], [r0["stress"][k][k] for k in range(3)], r0["t"]))
h = 2e-3
ep = ev(strain(cell, [[h,0,0],[0,h,0],[0,0,h]]), [[0,0,0]]); em = ev(strain(cell, [[-h,0,0],[0,-h,0],[0,0,-h]]), [[0,0,0]])
fd = (ep["energy"] - em["energy"]) / (2*h)   # = tr W
an = V * sum(r0["stress"][k][k] for k in range(3))
print("trace W: FD %.5f eV, analytic %.5f eV" % (fd, an))
# 2-atom cell with a displaced atom: force vs FD
c2 = [[4.05,0,0],[0,4.05,0],[0,0,4.05]]
p = [[0,0,0],[0,a,a],[a,0,a],[a,a+0.05,0]]
r = ev(c2, p); d = 1e-3
pp = [list(x) for x in p]; pp[3][1] += d; pm = [list(x) for x in p]; pm[3][1] -= d
fdf = -(ev(c2, pp)["energy"] - ev(c2, pm)["energy"]) / (2*d)
print("4-atom fcc, atom 3 displaced 0.05 A: F_y FD %.5f, analytic %.5f eV/A  (t=%.1fs)" % (fdf, r["forces"][3][1], r["t"]))
w.stdin.close(); w.wait()
