# Fuzz findings — collection viewport harness

**Task**: collection-arch-hardening, step 6. Producer: `test/flutter_cljd/internal/collection/fuzz_test.cljd`.
Consumer: step 7 (turn NEW signatures into red tests), steps 8–13 (the fixes).

## Campaign

| | |
|---|---|
| Episodes | 120 (12 `deftest` batches x 10 seeds) |
| Profiles | `:full` (design §6.2 weights), `:no-layout`, `:no-jump` — 40 seeds each |
| Episode | 24 fully parameterized ops, 600 items, viewport 400–900 x 600 |
| Failing episodes | 62 (51 `:check`, 9 `:o6`, 2 `:teardown`) |
| Raw oracle signatures | 15 → **10 distinct classes** (2 known, 8 NEW) |
| Runtime | ~59 s of `flutter test` (81 s wall incl. ClojureDart compile) |
| Shrinking | prefix cut + greedy removal, ≤24 replays per failure; typical result 1–6 ops |

Reproduce: `clojure -M:cljd test -- -t fuzz --timeout 20x`.
Every vector below is paste-ready for `harness/replay!` (`:o6` ones need the
anchor guard, so use `fuzz-test/run-ops!`, which is what found them).

**Standalone verification.** Each representative vector was re-run as the FIRST
episode of its own `deftest` (clean tester, no preceding rig). 12/12 reproduced.
Two earlier candidates did NOT and were discarded as cross-episode leakage:
a masonry `[[:pump 1]]` "cold-start tripwire" and a wrap `:o4` at seed 23. Both
were `flutter_test`'s deferred exceptions from the *previous* episode's trailing
pumps; `quiesce!` now drains them into a `:teardown` failure of the episode that
actually caused them. **Any future signature must be verified standalone before
being treated as an engine bug.**

---

## Known (do not re-analyze)

### known-A — far jump then layout morph loses the segment set-point
- Oracle: `:o1` `assert-materialization-bounded!`, as a `:teardown` failure (the
  tripwire fires on the settling pass *after* the last op).
- Seeds: 40, 19, +2. Shrunk (seed 40, `:masonry` `:cross0 280`, `:animate? true`):
  ```clojure
  [[:jump 16572.216153144836] [:layout :list]]
  ```
- Excerpt: `flow layout [:list 6.0] laid out 186 children in ONE pass — budget 92.88 for the 1267.66px band`.
- Same chain as `morph_red_test/known-red-far-morph-wrap-to-list-animated`, and
  the fuzzer reduced it to **two ops** — a far jump immediately followed by a
  layout swap, exactly the adjacency step 5 predicted.

### known-B — insert/remove above the window is not compensated
- Oracle: `:o6` `key-moved` (index held, key shifted by the count delta) or
  `:o6` `intra-drift` (key held, content slid by the inserted item's height).
- Seeds: 16, 205, 10, 17, +. Shrunk (seed 10, `:list`, `:approx 80.0`):
  ```clojure
  [[:drag -251.21960163116455] [:jump 17847.88191318512] [:drag -420.46699047088623]
   [:fling 3703.7041187286377] [:drag 352.2779846191406] [:fling 4489.896488189697]
   [:drag -253.9003086090088] [:drag -672.1532535552979] [:pump 3] [:pump 2]
   [:insert 2 900010]]
  ```
- Excerpts: seed 16 `before {:key 383 :index 384} after {:key 382 :index 384}`;
  seed 10 `key 205 held, intra 86.006 -> -2.994`.
- Same defect as `harness_test/known-red-insert-remove-above-window-anchor`.

---

## NEW

### NEW-1 — a layout morph loses the anchor with no far jump involved
- Oracle: `:o6` `key-moved`. Seeds 26, 207 (both `:animate? true`).
- Shrunk (seed 26, `:grid` -> `:list`, `:cross0 640`, `:approx 80.0`) — **3 ops**:
  ```clojure
  [[:settle] [:drag 576.4958953857422] [:layout :list]]
  ```
  `before {:key 6 :intra 144.496 :index 6}` -> `after {:key 7 :intra 29.491 :index 7}`
- Seed 207 (`:list` -> `:masonry` at ~6 viewports) slides 11 items: key 43 -> 54.
- Cause guess: the segment set-point is sampled from `restingTop`, which a plain
  drag leaves pointing at the leading-overscan child rather than the viewport-top
  item — the same `anchor-before` (pre-`walk!` cache) vs `capture-resting-top!`
  (post-layout children) split step 5 named, here without needing `restingTop` to
  be nil. **Not covered by the step-5 reds**, which all require a ~50-viewport jump.

### NEW-2 — wrap runs overlap after a backward drag
- Oracle: `:o4` `overlap` + `wrap-run-pitch`. Seeds 222, 22, 11.
- Shrunk (seed 222, `:wrap`, `:cross0 900`, `:approx 120.0`) — **3 ops**:
  ```clojure
  [[:drag 4606.196022033691] [:fling -2428.6102771759033] [:drag -522.4222183227539]]
  ```
  `overlap pairs [[303 304] [303 305] [319 320] [319 321]]`, `wrap-run-pitch {:from 303 :to 304}`
- Seed 11 reaches the same state in two ops from a jump: `[[:jump 2553.9976358413696] [:settle] [:drag -333.8862705230713]]`.
- Cause guess: backward seam refinement re-lays a partial run using a cross
  cursor that was not reset for the run it re-enters, so cells stack on the same
  cross offset. Physically overlapping cells at rest — a paint-visible defect.

### NEW-3 — wrap run advance is short by one run's max main
- Oracle: `:o4` `wrap-run-advance`. Seeds 240, 223.
- Shrunk (seed 240, `:list` -> `:wrap`, `:cross0 640`, `:approx nil`) — **3 ops**:
  ```clojure
  [[:drag 3357.876420021057] [:layout :wrap] [:settle]]
  ```
  `{:want 1960.84375 :from 1910.84375 :to 1954.84375}` — the next run starts 6 px
  (exactly one `run-spacing`) too early, i.e. the run's tallest cell was not used.
- Cause guess: after a morph the run-height accumulator is seeded from the frozen
  segment snapshot's per-cell mains rather than the re-measured ones, so a run
  whose tallest cell entered during the morph advances by the wrong max.

### NEW-4 — `type 'Null' is not a subtype of type 'RenderBox' in type cast`
- Oracle: `:o1` (hard cast failure inside layout). Seeds 5, 218.
- Shrunk (seed 5, `:grid` -> `:wrap`, `:cross0 900`, `:approx nil`) — 8 ops:
  ```clojure
  [[:fling 4392.159795761108] [:drag 805.7088088989258] [:drag 788.5633134841919]
   [:fling 7113.964653015137] [:remove 284] [:fling-catch -4022.814464569092 78]
   [:drag -474.3441581726074] [:layout :wrap]]
  ```
- Both instances end on a `[:layout ...]` op that follows a `[:remove ...]`.
- Cause guess: a layout pass dereferences a child the previous pass garbage-
  collected (or one the removal took out of the shadow index space) and casts the
  nil straight to `RenderBox` — a missing nil guard on a `childAfter`/`childBefore`
  walk in the morph path. Highest-severity NEW: it is an outright crash.

### NEW-5 — Flutter's own `childCount >= leadingGarbage + trailingGarbage` assert
- Oracle: `:o1`, `sliver_multi_box_adaptor.dart:594`. Seeds 6, 138 (both `:masonry`).
- Shrunk (seed 6, `:masonry`, `:cross0 900`, `:approx 120.0`) — 6 ops:
  ```clojure
  [[:drag 3345.6783771514893] [:cross 400.0] [:settle] [:cross 640.0]
   [:jump 26061.590909957886] [:to-top]]
  ```
- Cause guess: the engine asks `collectGarbage` to drop more leading+trailing
  children than are attached, i.e. the garbage counts are computed against a
  window from a *different* cross-extent generation — the `to-top` after a
  cross change re-anchors without recomputing the attached range.

### NEW-6 — a jump past the end parks `pixels` outside the advertised range
- Oracle: `:o5` `pixels-out-of-range`, 14 seeds (the single largest `:check` class).
- Shrunk (seed 121, `:wrap`, `:cross0 280`, `:approx 80.0`) — **1 op**:
  ```clojure
  [[:jump 27949.243783950806]]
  ```
  `{:pixels 18594.584 :max 18574.528 :min 0.0}` — still 20 px past `max` after the
  fuzzer's bounded 720 ms boundary wait (45 frames), so this is not the spring
  transient a forced `jumpTo` legitimately produces.
- Cause guess: `maxScrollExtent` keeps shrinking by a few px per frame as the
  estimate for the tail is replaced by measured cells, so the clamping ballistic
  never reaches a fixed point — the position chases a moving boundary.
- NOTE for step 7: assert this only AFTER a bounded wait; the first ~100 ms of a
  past-the-end `jumpTo` is legitimately out of range.

### NEW-7 — a count change plus a far jump walks the entire dataset in one pass
- Oracle: `:o9` `unbounded` on `cache-n`. Seeds 32, 122, +2.
- Shrunk (seed 32, `:list` -> `:wrap`, `:cross0 900`, `:approx 80.0`) — **3 ops**:
  ```clojure
  [[:remove 331] [:jump 16256.178617477417] [:layout :wrap]]
  ```
  `cache-n 601` (the whole dataset) with `pass-materialized 601`, `cache-first 0`,
  `anchored0 true`, `checkpoints (600)`.
- Cause guess: `update-render!` drops the checkpoints on a count change (known-B's
  root) but leaves `anchoredTo0`, so the next far jump has no checkpoint to seed
  from and walks from index 0 — O(n) per frame, the exact invariant the memory
  note "cache is accelerator only" says must never be violated. Same root as
  known-B, different symptom; worth one red test of its own.

### NEW-8 — `committed` outgrows the live window (suspect: oracle calibration)
- Oracle: `:o9` `unbounded` on `committed-n`, 17 seeds — the largest bucket.
- Shrunk (seed 134, `:wrap`, `:cross0 640`, `:approx 80.0`) — **4 ops**:
  ```clojure
  [[:remove 489] [:drag 3047.6452589035034] [:insert 583 900004] [:to-top]]
  ```
  `committed-n 103` vs limit 72; seed 135 `172` vs `160`.
- Every instance is mid-segment (`:tweening? true`, regime `:landing`/`:flying`)
  and the overshoot is small (1.1–1.5x). O9's bound is `8 * attached + 64`, which
  collapses when the window is nearly empty during a landing.
- **Classify as PLAUSIBLE, not confirmed**: decide in step 7 whether the segment's
  committed set is legitimately proportional to the morph span rather than to the
  attached window, and if so re-calibrate O9 instead of writing a red test.

---

## Harness / fuzzer defects found and fixed

1. `make-rig` reused the previous rig's element subtree (`KeyedSubtree` key was
   always `0`), so a second rig under one tester inherited the first's
   `SliverCollectionState` **and its render object's cache**. Step 4 saw this as
   "cross-rig pollution" and worked around it by splitting deftests; the fuzzer
   cannot. Fixed: `rig-gen`, a monotone counter, seeds every rig's `:cold-gen`.
   (`harness.cljd`, the only harness change in this step.)
2. Deferred exceptions raised by an episode's trailing pumps were handed to the
   NEXT episode's first `check!`. Fixed in the fuzzer: `quiesce!` drains them and
   reports a `:teardown` failure against the episode that caused them. This is
   what produced the two false candidates listed above.
3. The shrinker's failure signature was the oracle id alone, so it happily
   converged onto a *different* bug with the same id (every assert is `:o1`).
   Fixed: the signature carries the violation's `:why` / exception head.

## Oracle notes for step 7

- `check!` never runs O6 — it has no `before`. The fuzzer supplies it around
  layout swaps and strictly-above-window inserts/removes; a red test must do the
  same (capture `top-anchor` at rest, apply, `settle!`, `o6-anchor-preserved`).
- O5's `pixels-out-of-range` clause needs a bounded wait after any jump op
  (see NEW-6), otherwise it fires on the legitimate boundary spring.
