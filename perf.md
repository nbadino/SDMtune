## Performance optimizations (nbadino fork)

This fork adds parallel training and several micro-optimizations to SDMtune.
All changes are **opt-in**: `parallel = FALSE` preserves the original serial
behavior.

### Core parallelization (`50e94ff8`)

| Function | What's parallelized | How to enable |
|----------|---------------------|---------------|
| `gridSearch()` | All hyperparameter combinations | `parallel = TRUE` |
| `train()` | K CV folds | `parallel = TRUE` |
| `randomSearch()` / `optimizeModel()` | Initial population + GA children | `parallel = TRUE` |
| `varImp()` | CV fold permutations | `parallel = TRUE` |
| `doJk()` | Jackknife variable loop | `parallel = TRUE` |

On Linux/macOS, forked processes via `parallel::mclapply`. On Windows, falls
back to `makePSOCKcluster`. Each worker sets `OMP_NUM_THREADS=1` to prevent
OpenMP oversubscription with `glmnet`.

### Micro-optimizations

| Optimization | Commit | Effect |
|--------------|--------|--------|
| Remove `Sys.sleep(0.1)` from chart updates | `50e94ff8` | Saves N_models × 0.1s |
| Mutate SWD columns in-place in permutation importance | `f2310947` | Avoids copying SWD per permutation |

### Estimated speedup (8 cores, Maxnet)

| Operation | Before | After |
|-----------|--------|-------|
| `gridSearch` (9 combos × 5 fold CV) | ~35 min | ~3-5 min |
| `train` CV 5-fold | 5× T | ~1.3× T |
| `optimizeModel` (pop=15, gen=2) | ~30 min | ~5-8 min |
| `doJk` (10 vars + with_only) | ~10 min | ~2-3 min |

### Companion optimization in maxnet

This fork works best with [nbadino/maxnet](https://github.com/nbadino/maxnet)
which auto-imputes NA values and adds range-zero jitter, eliminating the need
for monkey-patches that previously added O(rows²) overhead per training call.
