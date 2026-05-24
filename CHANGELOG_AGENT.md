# CHANGELOG_AGENT.md

## Session Log

### Session 3 — Engine Layer Bug Fixes + Validation (2026-05-23)

**Date**: 2026-05-23
**Task**: Fix critical bugs in engine.c + add tests, complete Phase 6.1

**Bug Fixes**:
- **Segfault in `engine_random_init` / `engine_file_init`**: `char[MAX_MODES][32]` 2D stack array was cast to `const char *const *` and passed to `build_initialization_result`. The function called `strlen(amplitude_keys[i])` on each entry, reading the byte data of a 2D char array as if it were an array of pointers → garbage address → SEGV in `strlen`. Fixed by constructing a proper pointer array (`raw_keys[MAX_MODES][32]` + `amplitude_keys[MAX_MODES]` where `amplitude_keys[i] = raw_keys[i]`).
- **Inverted return value check**: `random_initialization()`, `manual_initialization()`, and `file_initialization()` all return N > 0 on success (number of matched coeffs/keys) and 0 on failure (no matches). The engine functions checked `if (ret != 0) return -1`, treating success as failure. Fixed to check `if (ret == 0) return -1`.
- **Typo `ctx.basis` vs `ctx->basis`** in `engine_random_init` (would have caused compilation error).

**What Was Tested**: Full build + `make test` — all 46 tests passing across 6 test suites

**Test Results**: 46/46 tests passing (domain: 3, symmetry_ops: 16, basis: 2, initializers: 9, field: 9, engine: 13)

**Remaining Issues**:
- Pm-3m space group yields 0 modes at N=5 due to `family_planes_info` hard limit at max_index=50 (pre-existing)
- 2 compiler warnings remain: unused `nx` parameter in `fftw3d_at`, const qualifier discard in `engine_manual_init`

**Next Recommended Task**: Phase 6.2 — `Analytic.py` translation (`calculate_square_norm()`)

**Git Commit**: `2235cdf fix: engine.c pointer array bug + hash-based family key dedup + full pipeline test`
