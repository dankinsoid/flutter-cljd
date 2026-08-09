# Fuzz findings — collection viewport harness

**Task**: collection-arch-hardening, step 6. Producer: `test/flutter_cljd/internal/collection/fuzz_test.cljd`.
Consumer: step 7 (turn NEW signatures into red tests), steps 8–13 (the fixes).

**Step-7 triage (2026-08-07)**: every NEW signature now has a deterministic
deftest in `test/flutter_cljd/internal/collection/fuzz_red_test.cljd`, which
replays the shrunk vector through `fuzz-test/run-ops!` and regenerates the rig
header from the seed. 6 stayed red (reproduced 3× each, `+2 -6` every run);
NEW-6 and NEW-8 were reclassified as fuzzer/oracle miscalibration, recalibrated,
and kept as GREEN tests guarding the recalibration. Default suite: 297 pass.

**Root-cause consolidation for step 8/9** — the 8 NEW signatures are 6 roots:
- **NEW-5 and NEW-7 are ONE root** (`keyed-tween-layout`'s nil→0 window fallback).
- **NEW-1 is NOT known-A's root**: known-A/NEW-7 are about the set-point being
  absent or the band being stale; NEW-1's set-point is present and correct, and
  its *model* (rigid translation) is what loses the anchor.
- **NEW-4 lives in the uncommitted `seed-cache!` WIP, not at HEAD** — step 5's
  "the WIP is inert" holds only for the far-morph scenario; here it crashes.
- NEW-2 and NEW-3 are wrap-tiling defects independent of the rebase work.

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
- Same defect as `harness_test/insert-remove-above-window-anchor`.
- Status: **FIXED (step 9b stage 6)** — a count change drops the geometry cache
  and the next `seed-cache!` re-anchors by `viewAnchor :key` on the surviving
  reconciled child; `segment-start!` resolves its set-point slot the same way.
  Green as `insert-remove-above-window-anchor` +
  `fuzz-red-test/fuzz-insert-above-window-holds-anchor` (this vector).

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
- Status: **FIXED (step 9b stage 4)** — green as `fuzz-morph-holds-anchor-without-jump`.
  The set-point lerps the anchor's EXTENT as well as its offset and preserves the
  consumed fraction (`tw/point-desired`); O6 checks that fraction, since a
  resizing anchor cannot hold the key and the absolute intra offset at once.
- Cause guess (step 6) was WRONG: `restingTop` is correct here.
- **CONFIRMED**: `segment-start!` builds `segAnchor` as `{:from old-offset :to
  new-offset :screen (from − scrollOffset)}` and the segment shift is the rigid
  `to − from` (render.cljd L2857-2881). That pins the anchor cell's TOP EDGE to a
  fixed screen position, not the viewport-top item. Here the anchor is item 6 at
  from 432 → to 516 with screen −144.5, and `pixels` follows exactly
  (576.50 → 660.49) — but item 6's own main extent shrank 208 (grid row) → 109
  (list cell), so its 144.5 px of already-consumed height no longer exists and the
  viewport top lands inside item 7. Any morph that resizes a partially
  scrolled-past anchor (intra > 0) reproduces this. The set-point is PRESENT and
  correct: this is a defect of the set-point MODEL, a different root from
  known-A/NEW-7, and step 8 must decide what the segment preserves (the anchor's
  consumed fraction, or a re-picked top item) rather than its top edge.

### NEW-2 — wrap runs overlap after a backward drag
- Oracle: `:o4` `overlap` + `wrap-run-pitch`. Seeds 222, 22, 11.
- Shrunk (seed 222, `:wrap`, `:cross0 900`, `:approx 120.0`) — **3 ops**:
  ```clojure
  [[:drag 4606.196022033691] [:fling -2428.6102771759033] [:drag -522.4222183227539]]
  ```
  `overlap pairs [[303 304] [303 305] [319 320] [319 321]]`, `wrap-run-pitch {:from 303 :to 304}`
- Seed 11 reaches the same state in two ops from a jump: `[[:jump 2553.9976358413696] [:settle] [:drag -333.8862705230713]]`.
- Status: **FIXED (step 9c stage 8)** — green as
  `fuzz-wrap-runs-stay-disjoint-after-back-drag`. `refine-seam` reports the cross
  mismatch alongside the main delta, and a cross-stale seam re-anchors on the
  layout's new `:run-start` (a FRESH run at the head) instead of on where `:place`
  puts the head. See the step-9c amendment below.
- **CONFIRMED**: observed children at rest — index 303 `{:offset 3444.45 :cross 0
  :cross-size 128}` and index 304 `{:offset 3444.45 :cross 0 :cross-size 94}`,
  with 305 at cross 100 (= 94 + spacing 6). Both 303 and 304 are in `checkpoints`,
  i.e. both were recorded as renewal points (`:renewal-point?` = cross 0 = run
  start). Root: `refine-seam-delta` (render.cljd L443-452) reconciles the backward
  re-flow against the cache head on the **MAIN axis only** — it compares
  `(:offset gm)` to `old-head` and ignores the `:cross` the same `:place` call
  returns. A re-flow whose prefix ends with an OPEN run (`{:x 128 :y 3444 :h 96}`)
  places the head at the same `y`, so `d` = 0, the seam is declared converged and
  `stitch-prefix` keeps the stale head verbatim at cross 0. Wrap's flow state
  carries a cross cursor that the seam protocol has no way to reconcile.
- **Step-9c amendment**: in this vector the dominant producer is not the stitch
  but the FROZEN (velocity-suppressed) re-anchor. The drag is fast, so `emit?` is
  false and the branch re-seeds the prefix at `:anchor(lo, first-offset − d)` with
  `d` measured from where `:place` puts the HEAD. When the prefix's last run is
  OPEN, that alignment drops the prefix's run start exactly onto the retained
  head's run — traced at cf=304: `d` = −143.54, and the re-anchored prefix put 303
  at the head's own 3572.45. Aligning on the next RUN start instead leaves a
  `:run-spacing` gap above the head, needs no correction and no in-window churn,
  and keeps the leading margin approximate as designed. The pure-stitch shape
  (`d` = 0 with a stale cross) is real too and takes the same branch.

### NEW-3 — wrap run advance is short by one run's max main
- Oracle: `:o4` `wrap-run-advance`. Seeds 240, 223.
- Shrunk (seed 240, `:list` -> `:wrap`, `:cross0 640`, `:approx nil`) — **3 ops**:
  ```clojure
  [[:drag 3357.876420021057] [:layout :wrap] [:settle]]
  ```
  `{:want 1960.84375 :from 1910.84375 :to 1954.84375}` — the next run starts 6 px
  too early: `want − from` = 63 = 57 + `run-spacing` 6, `to − from` = 57, so the
  advance is exactly the cell's own extent and the run gap is missing.
- Status: **FIXED (step 9c stage 9)** — green as
  `fuzz-wrap-run-advance-keeps-the-run-gap`. `leading-step-entry` backs off by the
  run PITCH (`:run-start`, `:flow-end` only as the fallback).
- **CONFIRMED**: the failing boundary sits in the above-window leading margin
  (`cacheFirst` 38, first checkpoint 40) and the run is single-celled, i.e. the
  `estimate-leading-step!` → `leading-step-entry` path (render.cljd L1784-1852).
  Its advance is `adv = max(ext, flow-end(env, st') − head)`; for wrap, `:anchor`
  seeds `{:x nil :y head :h 0}` and `:place` puts `i` INTO that same empty run, so
  `st'` keeps `:y = head`, `flow-end` returns `head`, the second term is 0 and
  `adv` collapses to `ext`. The cell is placed at `head − ext` instead of
  `head − ext − run-spacing`. Not the morph's snapshot — the morph only matters
  because it re-seeds the cache and forces the leading estimate to run.

### NEW-4 — `type 'Null' is not a subtype of type 'RenderBox' in type cast`
- Oracle: `:o1` (hard cast failure inside layout). Seeds 5, 218.
- Shrunk (seed 5, `:grid` -> `:wrap`, `:cross0 900`, `:approx nil`) — 8 ops:
  ```clojure
  [[:fling 4392.159795761108] [:drag 805.7088088989258] [:drag 788.5633134841919]
   [:fling 7113.964653015137] [:remove 284] [:fling-catch -4022.814464569092 78]
   [:drag -474.3441581726074] [:layout :wrap]]
  ```
- Both instances end on a `[:layout ...]` op that follows a `[:remove ...]`.
- Status: **RED-TESTED** — `known-red-fuzz-null-render-box-cast-in-morph-after-remove`.
- **CONFIRMED, and it is the UNCOMMITTED WIP, not HEAD.** Stack:
  `seed-cache! (render.dart:4221) ← flow-layout! ← segment-start! ← performLayout`.
  render.dart:4221 is `(rs.firstChild as RenderBox)` = render.cljd **L1430**,
  `^r/RenderBox fc2 (.-firstChild rs)` — a line that exists only in the working
  tree's `seed-cache!` WIP. The two lines above it
  (`(when (pos? lead) (.collectGarbage rs lead 0))`, L1428) destroy every attached
  child whose index is below `restingTop.idx`; when the whole attached window is
  above the resting anchor — which a backward `fling-catch` + drag produces —
  `lead` equals the child count, `firstChild` becomes nil, and the non-nullable
  cast on the next line throws before the `(if fc2 …)` / `(and fc2 …)` guards the
  WIP itself wrote can run. Step 5 recorded the WIP as "inert"; that holds only
  for the far-morph scenario. Fix shape: `^r/RenderBox?` plus re-checking
  emptiness after the GC — but the WIP is superseded by step 9 anyway.

### NEW-5 — Flutter's own `childCount >= leadingGarbage + trailingGarbage` assert
- Oracle: `:o1`, `sliver_multi_box_adaptor.dart:594`. Seeds 6, 138 (both `:masonry`).
- Shrunk (seed 6, `:masonry`, `:cross0 900`, `:approx 120.0`) — 6 ops:
  ```clojure
  [[:drag 3345.6783771514893] [:cross 400.0] [:settle] [:cross 640.0]
   [:jump 26061.590909957886] [:to-top]]
  ```
- Status: **FIXED (step 9b stage 5)** — green as
  `fuzz-garbage-counts-stay-within-child-count`. **Same root as NEW-7.**
  `keyed-tween-layout`'s window queries clamp a nil answer to the source's own
  edge; `indexed-layout!` refuses an inverted window (debug assert / release
  clamp).
- **CONFIRMED**: stack is `indexed_layout! → collectGarbage` — mid-segment the
  engine runs `indexed-layout!` over `segTween` (render.cljd L259), and its
  garbage counts come from `first-index` / `target-last` (L2504-2522).
  `keyed-tween-layout` resolves both through the frozen snapshot with a nil
  fallback of **0** (tween.cljd L352-357: `(int (or (to-first …) 0))`). After the
  `to-top` the viewport sits entirely BELOW the frozen window: `frozen-frame-source`
  answers `first-index(0)` = the snapshot base (73) but `last-index(850)` = nil →
  0. The window is INVERTED (73 > 0), so `calculateLeadingGarbage(firstIndex: 73)`
  counts indices 0..72 and `calculateTrailingGarbage(lastIndex: 0)` counts 1..98 —
  the same children twice, 171 > childCount 99. Neither the base window nor
  `widen-window-indices` can invert on its own; only the nil→0 fallback can.
  The fix must distinguish "no opinion" from "index 0", and `indexed-layout!`
  should refuse an inverted window regardless.

### NEW-6 — a jump past the end parks `pixels` outside the advertised range
- Oracle: `:o5` `pixels-out-of-range`, 14 seeds (the single largest `:check` class).
- Shrunk (seed 121, `:wrap`, `:cross0 280`, `:approx 80.0`) — **1 op**:
  ```clojure
  [[:jump 27949.243783950806]]
  ```
  `{:pixels 18594.584 :max 18574.528 :min 0.0}` — still 20 px past `max` after the
  fuzzer's bounded 720 ms boundary wait (45 frames), so this is not the spring
  transient a forced `jumpTo` legitimately produces.
- Status: **RECLASSIFIED — fuzzer miscalibration, not an engine leak.** Green test
  `fuzz-red-test/fuzz-jump-past-end-converges-into-range`; the fix is
  `fuzz-test/boundary-settle!`'s bound, 45 → 150 frames.
- **Evidence** (per-frame trace of the 1-op repro): `pixels` DOES converge, at
  frame 68, exactly onto `max`. The first half of the flight is the estimate
  collapsing — `maxScrollExtent` 23760 → 18335 over 18 frames as measured tail
  cells replace the `:approximate-item-size 80` estimate, then back UP to 18573 at
  frame 26 — and every dimension change restarts the boundary spring
  (`IdleScrollActivity`/`BallisticScrollActivity.applyNewDimensions` →
  `goBallistic`). Once the extent stops moving the spring converges
  asymptotically over ~40 more frames. At the old 45-frame bound the sample was
  18592.95, which is precisely the number the campaign reported. The extent's
  ~29 % transient error at the far end of an estimated wrap is by design
  (estimates outside the window), so nothing here is an engine defect.

### NEW-7 — a count change plus a far jump walks the entire dataset in one pass
- Oracle: `:o9` `unbounded` on `cache-n`. Seeds 32, 122, +2.
- Shrunk (seed 32, `:list` -> `:wrap`, `:cross0 900`, `:approx 80.0`) — **3 ops**:
  ```clojure
  [[:remove 331] [:jump 16256.178617477417] [:layout :wrap]]
  ```
  `cache-n 601` (the whole dataset) with `pass-materialized 601`, `cache-first 0`,
  `anchored0 true`, `checkpoints (600)`.
- Status: **FIXED (step 9b stage 5)** — green as
  `fuzz-count-change-then-far-jump-stays-window-bounded`. **Same root as NEW-5**,
  and NOT the same root as known-B. Beyond the window clamp, a segment now carries
  a validity `:domain`; a window that leaves it settles the segment in place so
  the resting far-jump inverse seed takes over (O(window)).
- Cause guess (step 6) was WRONG: `update-render!` DOES clear `anchoredTo0` on a
  count change (render.cljd L971-975).
- **CONFIRMED**: the `[:remove …]` is load-bearing because `:animate? true` makes
  it OPEN A SEGMENT, not because of its checkpoint invalidation. Engine state
  traced op by op: after the remove and after the far jump it is still
  `{:cache-first 0 :cache-n 11 :anchored0 true :checkpoints (0..9) :tweening? true}`
  — the jump moved `pixels` by 16 256 px and the band did not follow. Mid-segment
  the window comes from `segTween`, i.e. `keyed-tween-layout`'s nil→0 fallback
  (tween.cljd L352): with the viewport ABOVE the frozen snapshot `to-first`
  returns nil → `first-index` 0, so the band stays pinned at the top for the whole
  segment. `seed-cache!` cannot repair it either — its far-jump re-seed requires
  `(nil? fc)` (L1360) and a non-empty cache short-circuits to nil (L1374). The
  `[:layout :wrap]` then clears the cache, `seed-cache!` falls through to the
  `(zero? idx)` from-0 seed, and the capture pass walks the 16 000 px gap:
  `cache-n 600`, `pass-materialized 600`, `anchored0 true`.
  (Letting the segment settle first repairs it: the next resting pass re-seeds to
  `cache-first 178`. That is why only the jump→morph adjacency reproduces.)

### NEW-8 — `committed` outgrows the live window (suspect: oracle calibration)
- Oracle: `:o9` `unbounded` on `committed-n`, 17 seeds — the largest bucket.
- Shrunk (seed 134, `:wrap`, `:cross0 640`, `:approx 80.0`) — **4 ops**:
  ```clojure
  [[:remove 489] [:drag 3047.6452589035034] [:insert 583 900004] [:to-top]]
  ```
  `committed-n 103` vs limit 72; seed 135 `172` vs `160`.
- Status: **RECLASSIFIED — O9 miscalibration.** `harness/o9-bounded-state` now
  bounds `committed-n` by `8 * max(attached, cache-n) + 64` while a segment runs;
  `cache-n` and `checkpoints` keep the attached-window bound. Green test
  `fuzz-red-test/fuzz-committed-stays-band-proportional-mid-segment`.
- **Evidence**: `commit-pass!` upserts every attached child's rect and prunes to
  the cache band ONLY when `commit-prune?` allows it — never mid-segment, so that
  the segment's `from` capture and its entering/exiting/leaving participants
  survive (render.cljd L838-875). Mid-segment retention is therefore proportional
  to the BAND times the segment's pass count, not to the attached window, which a
  landing collapses: the repro shows `committed-n 103` with **one** attached child
  and `cache-n 55`. The oracle was measuring the wrong denominator; genuine
  unbounded growth (the whole dataset) still fails, and NEW-7's `cache-n 600` is
  unaffected.

---

---

# Campaign 2 (step 10, 2026-08-09) — post-rebase re-campaign

Run after rebase stages 0-9 (9a/9b/9c) and the TEMP-instrumentation revert.

| | |
|---|---|
| Episodes | **220** (22 `deftest` batches x 10 seeds) |
| Profiles | `:full` 1-60, `:no-layout` 101-160, `:no-jump` 201-260, **`:churn` 301-340** |
| Failing episodes | 55 raw → **43** after two fuzzer recalibrations (campaign 1: 62/120) |
| Distinct classes | **7** (5 red-tested as NEW-9..NEW-13, 2 singletons recorded only) |
| Runtime | 74 s of `flutter test`, ~4 min wall incl. the ClojureDart recompile |

**Campaign-1 signatures that did NOT return**: known-A, known-B's plain form,
NEW-1..NEW-5, NEW-7. Their ratchet tests stay green and their campaign-1 seeds
now run to completion. The `:no-jump` `:o4 wrap-run-pitch` shape 9c flagged as a
watch item did **not** reappear — the frozen seam's loose leading margin is
inside the oracle's tolerance.

**`:churn` profile** (new): `:layout` starved and `:animate?` forced on, so a
count change is the only segment opener and the next one lands mid-segment. Its
`:remove-anchor` op deletes the item the key re-anchor looks for. 6 of its 40
episodes failed, all reproducing classes the other profiles also produce — no
churn-only signature, which answers rebase-design open Q3 in the negative: the
count-change fallback when the anchor's child does not survive did not produce a
class of its own.

## NEW-9 — the capture pass produces a rebase with no absorption channel

- Oracle: `:o1`, the stage-3 absorption assert
  (`capture produced rebase N with no absorption channel (segAnchor nil)`).
- **17 seeds, all four profiles** — the single largest class, and the only one
  whose tripwire did not exist before step 9a.
- Shrunk to **3 ops** by six independent seeds, one shape:
  ```clojure
  [[:remove 222] [:drag 1211.795711517334] [:layout :list]]     ; seed 256
  [[:cross 900.0] [:drag 897.3258590698242] [:layout :wrap]]    ; seed 250
  [[:insert 467 900005] [:fling 6780.926990509033] [:layout :list]] ; seed 228
  ```
  A cache invalidator (count change or cross change) → a motion op → a morph.
- Red: `known-red-fuzz-capture-rebase-has-no-absorption-channel` (seed 256;
  reproduced standalone 3/3 across three seeds).
- **Classification: STOP-worthy.** The assert is 9a stage-2/3's own exactly-once
  contract, so the invariant step 9 set out to establish does not hold. The
  underlying condition (a produced rebase that `segAnchor` cannot carry) is
  step 5's link (c) in a new dress — 9a made `viewAnchor` non-nil while children
  exist, but a non-nil `viewAnchor` does not guarantee a RESOLVABLE `segAnchor`:
  traced state shows `viewAnchor {:idx 242}` against `cacheFirst 250`, i.e. the
  anchor index is outside the capture pass's own re-seeded cache, so the frozen
  `to` lookup answers nil. In release the assert is compiled out and the rebase
  is silently discarded — the far-morph symptom class.
- Not fixed here (step 10 is verification only). It blocks step 11: the
  capture-mode design returns the rebase as a VALUE, which presumes a channel.

## NEW-10 — a morph holds the fraction exactly and applies it to the next item

- Oracle: `:o6` `key-moved`. Seeds 241, 207, 52, 235.
- Shrunk (seed 241, `:no-jump`, `:list` -> `:grid`) — 5 ops:
  ```clojure
  [[:drag 770.4910182952881] [:fling 1578.5015106201172]
   [:drag 739.6548843383789] [:pump 1] [:layout :grid]]
  ```
  `before {:key 22 :frac 0.6774226815552667}` -> `after {:key 21 :frac 0.6774226815552665}`
- The consumed fraction is preserved to 1e-15 — 9b's stage-4 set-point works —
  but it is applied to the anchor's NEIGHBOUR, so the slot resolution is off by
  one. Distinct from NEW-1, whose fraction was not preserved at all.
- Red: `known-red-fuzz-morph-slides-anchor-one-item`.

## NEW-11 — known-B's shape returns when a far jump precedes the insert

- Oracle: `:o6` `key-moved`. Seed 16 (known-B's own campaign-1 seed).
- Shrunk — the `[:jump 17467.40162372589] [:insert 74 900008]` adjacency is
  load-bearing; the same insert without the jump is green (9b's ratchet).
  `before {:key 383 :index 384}` -> `after {:key 382 :index 384}`, intra
  identical: the re-anchor held the INDEX, not the key.
- Red: `known-red-fuzz-insert-above-window-after-far-jump`.

## NEW-12 — a rotation advertises a scrollExtent below the content it laid

- Oracle: `:o5` `extent-below-content`. Seeds 4, 17, 23, 38, 140, 322.
- Shrunk (seed 4, `:full`) — 4 ops:
  ```clojure
  [[:cross 400.0] [:jump 17421.205043792725] [:drag 3253.394079208374] [:cross 640.0]]
  ```
  `scroll-extent 27192.0` vs `last-end 32392.0`; seed 38 is worse (19288 vs 46048).
- Same oracle as step 5's `far-morph-capture-extent-truncated`, which stage 3
  closed for the flow capture path. The cross-change path still under-advertises.
- Red: `known-red-fuzz-rotation-extent-covers-laid-content`.

## NEW-13 — a morph after a far jump still caches most of the dataset

- Oracle: `:o9` `unbounded` on `cache-n`. Seeds 37, 222, 220, 2, 29.
- Shrunk (seed 37, `:full`) — 5 ops:
  ```clojure
  [[:remove 45] [:layout :list] [:jump 20821.444988250732] [:layout :wrap] [:remove 131]]
  ```
  `cache-n 370-403` against a limit of 72, with `pass-materialized` equal to it —
  genuine per-pass work, not retention.
- NEW-7's shape minus the count-change adjacency its diagnosis relied on, so the
  stage-5 `:domain` settle does not cover it.
- Red: `known-red-fuzz-morph-after-far-jump-stays-window-bounded`.

## Recorded, not red-tested (single seed, long vector)

- `:o4 hole` at seed 110 (`:no-layout`, 9 ops, ends on a drag after a far jump).
- `:o6 intra-drift` at seed 55 (`:full`, 7 ops, ends on an above-window insert).
- Flutter's `sliver.dart` assert at seed 251 (`:no-jump`, masonry, 8 drags) —
  NEW-5's neighbourhood, but a different assert line.

## Open calibration question — O9's mid-segment `cache-n` denominator

4 of the 9 `:o9` failures are `committed-n` 3-11 % over its bound and 4 more are
`cache-n` over a limit of 72, i.e. `8 * attached + 64` with **one** attached
child: mid-segment a landing collapses the attached window, which is exactly the
denominator problem step 7 fixed for `committed-n` and not for `cache-n`.
Loosening it was rejected here — every candidate denominator
(`max(attached, cache-n)`, `max(attached, pass-materialized)`, `committed-n`)
also masks seed 37's genuine `pass-materialized 403`. The strengthening move is
to give O9 a separate per-pass WORK probe on `pass-materialized` and let the
retention bounds relax; that is oracle design, deferred out of step 10.

## Fuzzer recalibrations (step 10)

6. `boundary-settle!` ran only after the jump ops, so a drag or fling released
   past a boundary — and any op that MOVES `maxScrollExtent` under a resting
   offset — was sampled mid-spring. 10 false `:o5 pixels-out-of-range` seeds;
   two were traced converging exactly onto the boundary (seed 123 −525 → 0.0,
   seed 23 8740 → 8578 = max). The bounded wait now runs after every op.
7. The O6 anchor guard fired on count changes taken with `pixels` resting ON
   `maxScrollExtent`, where a REMOVAL shrinks the range and drags the viewport
   with it — the top item legitimately changes (seed 105: max 18334.84 → 18033.84,
   `pixels` clamped with it). 3 false `:o6` seeds. `:remove` and `:layout` are
   now unguarded at the end of the range; `:insert` stays guarded, since growth
   never forces a clamp. Coverage note: a bottom-pinned anchor is now unguarded
   for those two ops, which is why seed 146's genuine engine miss (insert above
   the window at `pixels == max`, index held) is recorded here rather than caught
   by the guard.

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
4. (step 7) `boundary-settle!`'s 45-frame wait was shorter than a past-the-end
   jump's boundary spring on an estimated extent — 14 false `:o5` seeds. Now 150.
5. (step 7) O9 bounded `committed-n` by the ATTACHED window, which a landing
   collapses to one child, while the engine scopes `committed` to the cache band
   between segments — 17 false `:o9` seeds. Mid-segment the bound is now the band.

## Oracle notes for step 7

- `check!` never runs O6 — it has no `before`. The fuzzer supplies it around
  layout swaps and strictly-above-window inserts/removes; a red test must do the
  same (capture `top-anchor` at rest, apply, `settle!`, `o6-anchor-preserved`).
- O5's `pixels-out-of-range` clause needs a bounded wait after any jump op
  (see NEW-6), otherwise it fires on the legitimate boundary spring.
