# Cu-Au forge.materials.v0 fixture

Produced at commit b4d2b53 (branch `scs-bot/dft-bridge`) with

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
./build/atomforge search Cu Au --quick --mc-sweeps 20 --md-steps 200 --out /tmp/fx --materials
```

(23 s on one core). Only the sidecar and the `report.json` it references are kept; the 46 structure files are not.
`produced_by.commit` and `produced_by.timestamp` differ from run to run, so readers' tests should not pin them.
For qpuforge's `forge.materials.v0` reader tests.
