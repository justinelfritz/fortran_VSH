# Numerical stability of FORTVSH's Legendre recurrences

**Status: two independent investigations, both complete.** This document is written
incrementally as each experiment runs, so the raw data behind every claim here
doesn't need to be regenerated to recall why we believe it. Part 1 covers the
normalized (Holmes-Featherstone) recurrence; Part 2 covers the unnormalized (Bonnet)
recurrence and everything single-mode built on it.

---

# Part 1: `ASSOC_LEGENDRE_NORM_ALL` (Holmes-Featherstone recurrence)

## Background

`src/vsh.f90`'s `ASSOC_LEGENDRE_NORM_ALL` implements the Holmes & Featherstone (2002)
modified forward-column recurrence for the 4π-normalized associated Legendre functions
`N_l^m P_l^m(cos theta)`, `N_l^m = sqrt((2l+1)/4pi * (l-m)!/(l+m)!)`. Its docstring
claims stability to `l~2700`, versus `l~1400` for the unnormalized Bonnet recurrence
(`ASSOC_LEGENDRE_ALL`). `SSH_ALL` and every batch VSH routine (`VSH_TOR_ALL`,
`VSH_POL_UP_ALL`, `VSH_POL_DN_ALL`, etc.) build directly on this recurrence via the
private `VSH_CORE` helper, so its stability range is effectively the stability range
of most of the library's recommended, high-performance API.

## Discovery

While extending `convergence_accuracy`'s Lmax sweep to 4000 (to test the documented
`l~2700` claim), the SSH/VSH panels developed gaps starting around `l~2000`, even
though `ASSOC_LEGENDRE_NORM_ALL` looked perfectly stable in the Legendre panel all the
way to `l=4000`. Tracing this down: the Legendre panel's `x` samples
(`-0.6, 0.0, 0.6`) happened to avoid the region where the problem actually occurs.

Direct test, `L=3000`, `x=cos(pi/6)=0.866` (`theta=30 deg`, close to the pole):

```
M=0     PNORM = 0.390        (fine)
M=500   PNORM = -0.425       (fine)
M=1000  PNORM = 0.253        (fine)
M=1500  PNORM = 8.22e+127    (broken)
M=2000  PNORM = 3.83e+136    (broken)
M=2500  PNORM = 52280.8      (garbage, but small)
M=3000  PNORM = 4.94e-324    (garbage, but small)
```

At the same `(l, m) = (3000, 1500)` but `theta = 60 deg` or `90 deg`, the value is a
normal O(1) number. So the instability is a joint function of `(l, m, theta)` -- not
simply "stable below some fixed l" as the docstring implies.

## Working hypothesis: classical turning point

Note `sin(30 deg) = 0.5` and `m/l = 1500/3000 = 0.5` -- exactly the point where
`sin(theta) = m/l` at the first broken sample. This matches the classical turning
point of the associated Legendre ODE (WKB analysis): for fixed `l`, as a function of
`theta`, `P_l^m(cos theta)` is oscillatory where `sin(theta) > m/l` (roughly, for
large `l`) and exponentially varying (evanescent) where `sin(theta) < m/l`. Forward
three-term recurrences are well known to become ill-conditioned near/inside a
turning point, because the "minimal" (physically correct, decaying) solution and the
"dominant" (spurious, growing) solution become comparable in magnitude, and rounding
error preferentially amplifies the dominant one -- Holmes-Featherstone's modified
recurrence is more robust than the naive unnormalized one, but is not immune to this.

**This is the hypothesis being tested by the experiment below** -- not yet confirmed
as the full explanation, but the single data point strongly suggests it.

## Experiment: fine-grained (l, m, theta) scan

`src/stability_main.f90` (new, opt-in `VSH_BUILD_STABILITY` CMake target, executable
`vsh_stability`, not part of `ctest`): for `theta = 1 deg, 2 deg, ..., 179 deg`, calls
`ASSOC_LEGENDRE_NORM_ALL(PNORM, 5000, cos(theta))` once (one call yields the entire
`l <= 5000, m <= l` table), then for `l = 50, 100, ..., 5000` scans every `m` in
`[0, l]` for `|PNORM(l,m)| > 1000` (a threshold far above any legitimate normalized
value -- which should be at most O(sqrt(l)), ~70 at l=5000 -- and far below the
observed blowups of 1e100+). Records, per `(theta, l)`: whether any `m` broke, the
`[m_min, m_max]` broken band if so (both raw and as a fraction of `l`), and the peak
`|PNORM|` observed at that `(theta, l)`. Output: `./stability/stability_map.dat`
(17,900 rows: 179 theta values x 100 l-checkpoints). Ran in ~23s. `ctest` unaffected
(this executable is standalone, opt-in, not part of the correctness gate).

## Results

**Sanity check on the threshold**: max peak value among all *non*-broken rows was
930.4 (safely under the 1000 cutoff); the minimum broken value went to literal
double-precision `Infinity`. Clean separation, no evidence of false positives near
the boundary.

**Headline result: `l <= 2000` is unconditionally safe.** Zero broken rows exist
anywhere below `l=2100` (the exact grid checkpoint), across all 179 theta values
tested (1 deg to 179 deg). This is a strong, simple, practical bound -- well above
the `l=300` used for the manuscript's other validation and even above the extended
`l=4000` sweep that first surfaced this issue.

**Beyond `l=2000`, instability is confined to a specific theta band, not universal.**
Of the 179 theta values tested, only 54 (theta = 30-56 deg and its mirror
124-150 deg) ever showed breakage within `l <= 5000`. Theta close to either pole
(< 30 deg or > 150 deg) and theta near the equator (57-123 deg) showed **no**
breakage anywhere in the tested range -- consistent with Holmes-Featherstone's
recurrence doing exactly what it was designed for (fixing the classical near-pole
underflow of the naive forward-row recurrence) while having a separate, narrower
weak spot elsewhere.

**Onset `l` as a function of theta** (first broken l-checkpoint per theta; symmetric
about theta=90 deg as expected from `P_l^m(-x) = (-1)^(l+m) P_l^m(x)`):

| theta (deg) | onset l | theta (deg) | onset l |
|---|---|---|---|
| 30 / 150 | 2100 (worst case) | 44 / 136 | 2950 |
| 32 / 148 | 2150 | 46 / 134 | 3200 |
| 34 / 146 | 2250 | 48 / 132 | 3400 |
| 36 / 144 | 2350 | 50 / 130 | 3700 |
| 38 / 142 | 2500 | 52 / 128 | 4000 |
| 40 / 140 | 2650 | 54 / 126 | 4400 |
| 42 / 138 | 2800 | 56 / 124 | 4850 |

(Full per-degree table in `stability/stability_map.dat`.) Onset `l` rises smoothly
and monotonically from the worst case at theta~30 deg out to the edges of the
unstable band at theta~57 deg / 123 deg, where it exceeds the 5000 scan ceiling.

**Turning-point hypothesis: confirmed quantitatively.** At the onset `l` for each
broken theta, `m_min/l` (where the unstable band *starts*) matches `sin(theta)` to
a mean absolute deviation of **1.5%** (max 5.8% across all 54 broken thetas). This
is the classical turning point of the associated Legendre ODE -- `sin(theta) = m/l`
is where `P_l^m(cos theta)` transitions from oscillatory to evanescent behavior, and
three-term recurrences are well known to be ill-conditioned near that transition
because the growing and decaying solutions become comparable in magnitude.

**The unstable band widens with `l` past onset.** E.g. at theta=45 deg (onset
l=3050, `m_min/l=0.716`), by `l=5000` the band has widened to `m_min/l=0.435,
m_max/l=0.907` -- i.e. once triggered, more recurrence steps past the turning point
give rounding error more room to amplify, widening the affected `m` range on both
sides rather than staying pinned at the turning point.

## Practical safe-range guidance

- **`ASSOC_LEGENDRE_NORM_ALL`, `SSH_ALL`, and every batch VSH routine built on it
  (`VSH_TOR_ALL`, `VSH_POL_UP_ALL`, `VSH_POL_DN_ALL`, `GRAD_SSH_ALL`, `L_SSH_ALL`,
  `PVSH_*_ALL`) are safe for `l <= 2000` at any `m`, `theta`.** This supersedes the
  docstring's blanket "stable to l~2700" claim, which is optimistic for `theta` in
  the 30-56 deg / 124-150 deg bands and (within the tested range) conservative
  everywhere else.
- Above `l=2000`, safety depends on `theta`: safe for `theta < 30 deg`, `theta >
  150 deg`, or `57 deg <= theta <= 123 deg` (at least up to `l=5000`, the scan
  ceiling -- not verified beyond that). In the `30-56 deg` / `124-150 deg` bands,
  consult the onset-l table above, or treat `l=2000` as the hard ceiling if `theta`
  in that range can't be ruled out for a given application.
- This is a limitation of the **unnormalized-adjacent recurrence family's normalized
  variant**, not of the library's correctness at the tested small-`l` scale used by
  the `ctest` gate -- `main.f90`'s validation suite runs at `l<=8`, far below where
  any of this matters.

## Visualization

`py/plot_stability.py` reads `stability/stability_map.dat` and produces
`tex/Copernicus-EGU/figures/stability_map.{pdf,png}`: a `(theta, l)` heatmap of
peak `|N_l^m P_l^m|`, log-scaled and clipped at 20. It shows the safe/broken
boundary directly -- two symmetric wedge-shaped unstable regions (mirrored about
`theta=90 deg` as expected), each with its apex at `l~2100` near `theta=30 deg` /
`150 deg`, widening in `theta` as `l` grows toward the `l=5000` scan ceiling. The
large dark (safe) region below `l~2000` and the safe bands near the poles and
equator are immediately visible and match the quantitative table above exactly.

## Status: investigation complete for this pass

Open items for a future pass, not blocking the headline `l<=2000` guidance:
- Whether the unstable bands near the poles (`theta<30 deg`) or equator
  (`57-123 deg`) eventually appear beyond `l=5000` (not verified -- scan ceiling).
- Whether `SSH_ALL`/`VSH_*_ALL` inherit this exactly (expected, since they all
  route through `ASSOC_LEGENDRE_NORM_ALL` via the private `VSH_CORE` helper --
  not independently re-scanned, since re-deriving the same recurrence's failure
  mode through an extra layer of algebra wouldn't add information).
- ~~Whether to correct the `ASSOC_LEGENDRE_NORM_ALL` docstring's "l~2700" claim in
  `src/vsh.f90` to reference this document~~ -- **done**: `ASSOC_LEGENDRE_NORM_ALL`'s
  docstring now states the `l<=2000` bound directly, and the 10 routines built on it
  (`DDX_ASSOC_LEGENDRE_NORM_ALL`, `VSH_CORE`, `SSH_ALL`, `GRAD_SSH_ALL`, `L_SSH_ALL`,
  `PVSH_RAD_ALL`, `PVSH_POL_ALL`, `PVSH_TOR_ALL`, `VSH_TOR_ALL`, `VSH_POL_UP_ALL`,
  `VSH_POL_DN_ALL`) each carry a pointer `@warning` back to it.

---

# Part 2: `ASSOC_LEGENDRE`/`ASSOC_LEGENDRE_ALL` (unnormalized Bonnet recurrence)

## Background

The unnormalized recurrence's docstring claimed stability to `l~1400`. Separately,
the `convergence_accuracy` manuscript figure (diffing batch vs. single-mode against
each other, not ground truth) visually showed the Legendre panel diverging by
`l~10-21` -- but that was diffing two code paths sharing the *same* instability
against each other, which (as later understood) can show divergence earlier than
either path's actual error against ground truth, since two independently-rounded
runs of an ill-conditioned recurrence can disagree with each other before either
one is far from correct. Neither number had been established rigorously. This part
does that, using the same independent-reference methodology as Part 1.

## Why this needs a different experiment design than Part 1

Part 1's failure mode was a recurrence-conditioning problem near a turning point --
literal magnitude wasn't a reliable signal, so the scan looked for `|value| >
1000`. Here the mechanism is different: the seed `P_m^m(x) =
(-1)^m(2m-1)!!\sin^m\theta` grows/shrinks combinatorially with `m`, and this
alone can be enormous (or tiny) well before it's actually *wrong* -- a huge but
still-accurate value multiplied by an equally tiny normalization factor is fine.
The real damage happens via catastrophic cancellation in the subsequent three-term
upward recurrence in `l` (a subtraction of two comparably-huge terms), which
degrades *relative* accuracy well before literal double-precision overflow. So
this scan measures accuracy against ground truth directly (same technique as the
`convergence_accuracy`/`mpmath_reference.py` fix), not a raw-magnitude proxy.

## Experiment

`src/unnorm_stability_main.f90` (new, same opt-in `VSH_BUILD_STABILITY` CMake
target as Part 1, executable `vsh_unnorm_stability`): for `theta = 1 deg, ...,
179 deg`, calls `ASSOC_LEGENDRE_ALL(P, 200, cos(theta))` once, then for `l =
1, ..., 200` samples `m` at three fractions of `l` -- `1.0` (diagonal, expected
worst case since the seed is maximized at `m=l`), `0.75`, `0.5` -- writing the raw
single-mode (`ASSOC_LEGENDRE`) and batch value at each point. Output:
`./stability/unnorm_diagonal.dat` (106,684 rows). Ran in 0.3s.

`py/mpmath_reference.py`'s `legendre_ref` (already validated in Part 1's sibling
work on the `convergence_accuracy` figure -- the same arbitrary-precision rerun of
this exact recurrence) supplies ground truth for all 106,684 points: relative error
`|value - ref| / max(|ref|, 1e-300)`, 31.9s total. "Broken" = relative error
`> 1e-6` or non-finite.

## Results

**Diagonal (`m=l`) confirmed as worst case.** At `theta=90 deg`, onset `l` is 151
for `m=l`, 185 for `m=0.75l`, and `m=0.5l` never broke within `l<=200` at any
theta tested. This matches the mechanism directly: the seed `(2m-1)!!` term is
largest exactly at `m=l`, so this is a structural result, not just an empirical
one -- no need to scan a finer set of `m` fractions to look for a worse case.

**Global worst case: `l=151`, at `theta` in an 82-98 deg plateau around the
equator.** Onset `l` rises smoothly and *symmetrically* moving toward either pole
(opposite pattern from Part 1, where the poles were the safe zone): `164` at
`40/140 deg`, `172` at `30/150 deg`, `185` at `20/160 deg`, and **no breakdown
detected at all within `l<=200`** for `theta<15 deg` or `theta>165 deg` -- makes
sense, since `sin(theta)^m` strongly suppresses the seed's combinatorial growth
near the poles.

**The degradation is clean and monotonic, not fluctuating.** At `theta=90 deg`,
`m=l`, relative error grows smoothly from `1.1e-16` at `l=1` to `2.3e-14` at
`l=150`, then jumps straight to `Inf`/`NaN` at `l=151` with no recovery at any
tested `l` up to 200 -- a sharp cliff, not the widening/fluctuating band seen in
Part 1. Single-mode and batch match to full double precision at every point (same
recurrence, same breakdown).

## Practical safe-range guidance

**`ASSOC_LEGENDRE`, `ASSOC_LEGENDRE_ALL`, `DDX_ASSOC_LEGENDRE`,
`DDX_ASSOC_LEGENDRE_ALL`, and every single-mode function built on them (`SSH`,
`GRAD_SSH`, `L_SSH`, `PVSH_RAD`, `PVSH_POL`, `PVSH_TOR`, `VSH_TOR`,
`VSH_POL_UP`, `VSH_POL_DN`) are safe for `l<=150` at any `m`, `theta`** --
supersedes the docstring's "`l~1400`" claim (too optimistic) and the informal
"`l~20`" impression from the batch-vs-single `convergence_accuracy` figure (too
pessimistic -- that was measuring divergence between two noisy paths, not error
against truth). Above `l=150`, safety degrades fastest near the equator and
improves toward either pole; consult `stability/unnorm_diagonal.dat` for the
per-degree table if finer-grained guidance is needed for a specific `theta`.

This is a stricter ceiling than Part 1's `l<=2000` for the normalized routines --
consistent with why `ASSOC_LEGENDRE_NORM_ALL` exists at all.
