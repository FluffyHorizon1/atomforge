# Contributing to atomforge

Thank you for helping with atomforge. This repository is **header-only C++17** (the implementation lives under `include/atomforge/`, with a small CLI in `src/`).

## Licence

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

atomforge is proprietary software, licensed under the [Secured Cyber Solutions Proprietary Software Licence v1.0](LICENSE). It is **not open source**: no use, copying, modification or distribution is permitted without written permission from SCS. Third-party components are listed in [NOTICE](NOTICE).

By contributing, you confirm that you have the right to submit the change and that SCS may use it under that licence.

## Branch and pull-request workflow

- Do your work on a **new branch**. Open a pull request for review.
- **Never push to `main` or to any `release/*` branch.** Do not merge, close, tag, or change repository settings unless the owner has approved that specific action.
- Prefer a short, descriptive branch name (`fix/…`, `docs/…`, `chore/…`, `feat/…`, `ci/…`, `test/…`).
- Keep pull requests focused. Describe what changed and how you tested it.

## Build and test

Current `main` builds the CLI only. There is no CTest suite on `main` yet (that work lives on unmerged branches such as `chore/atomforge-tests-suite` and `ci/actions-ctest`).

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/atomforge validate          # 17 elements vs experiment
./build/atomforge relax L1_2 Cu Au  # one structure + elastic constants
./build/atomforge search Cu Au      # full binary search -> results/CuAu/
```

Please run `./build/atomforge validate` before opening a pull request. If your branch also has the CTest suite:

```bash
ctest --test-dir build --output-on-failure
```

## Code style

- C++17, no compiler extensions required of contributors. Keep the library **header-only**.
- Prefer small, focused headers under `include/atomforge/`. Match the existing 4-space indentation, naming, and include style.
- Do not add a build-time dependency without a clear reason in the pull request.
- Leave committed `results/` fixtures in place unless a change is specifically about them.

Licensing enquiries: https://securedcybersolutions.co.uk
