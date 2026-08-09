# Task: Collection engine architectural hardening

**Slug**: collection-arch-hardening
**Started**: 2026-08-07
**Status**: in-progress

## Goal
Eliminate the class of cross-subsystem device-round bugs in the collection
virtualization engine by (1) building a RenderViewport simulation harness +
invariant fuzzer so inter-frame bugs become CI failures, (2) making the
coordinate rebase a first-class exactly-once structure, (3) giving the segment
capture pass an explicit corrections-off mode, and (4) migrating the engine to
anchor-primary scroll truth (anchor index + intra-item offset as source of
truth, absolute offset derived). The active far-scroll wrap→list morph bug is
reproduced in the harness first (no device run) and closed by the rebase work.

## Raw request
> и так я снова хочу сделать ревью нашего collection. у него уже ~10 device-round
> фиксов - вот я думаю может у нас где то всё таки архитектурные проблема рас
> всплывают разные баги. [...review agreed...] /orchestrate
> Scope: всё, включая anchor-primary. Активный баг: воспроизвести в harness.

## Code map
(from this session's review — verified 2026-08-07)
- Engine: `src/flutter_cljd/internal/collection/render.cljd` (3004 lines) —
  `CollectionRenderSliver` deftype ~40 mutable fields (L49-143); flow driver
  `flow-layout!` (L2305), indexed driver (L2487), `flying-flow!` (L2203),
  `segment-start!` (L2733), `seed-cache!` (L1333), `backfill-leading!` (L1900).
- Rebase mechanisms to unify: `reanchorShift`/`crossLayoutReanchor` (L121-128,
  consumed L1380-1447, emitted L2401-2405), seam refinement (L1993-2024),
  exact-landing one-shots (L509-531, L1930-1964), top-underflow (L481-507).
- Layout model: `internal/collection/layout.cljd`; tween/frame sources:
  `internal/collection/tween.cljd`; host: `internal/sliver_collection.cljd`;
  box host: `internal/collection/box.cljd`.
- Debug tripwires (fuzz oracles): `assert-window-canonical!` (L2075),
  `assert-materialization-bounded!` (L2133).
- Active bug chain (hypothesis): capture pass ends on reanchorShift
  correction-only geometry → `live-only-flow-window` reads scrollExtent=0
  (L1133-1135) → snapshot :extent truncated; the correction itself is discarded
  by `from-relay!` overwrite; if seg-anchor is nil (L2857-2868) the rebase is
  lost entirely. TEMP instrumentation at HEAD (fb24f35) — revert when closed.
- Build/verify: `bin/check` (compile example + dart analyze); tests run via
  `clojure -M:cljd test` (bare `flutter test` runs STALE code). Widget gotcha:
  `(= w w)` false for const-foldable widgets under --track-widget-creation.
- Existing tests: 183+ pure/detached-render-object tests, NO viewport harness.

## Decisions log
- 2026-08-07: Review conclusions accepted by user: (a) window-bounded concept is
  sound; (b) real defect = rebase spread over 4 mechanisms with one-shot flags,
  no exactly-once owner; (c) capture pass runs full armed machinery, gating is
  scattered; (d) offset-primary should become anchor-primary; (e) edge-slide is
  NOT overengineering (pays for the window-bound invariant) — only its
  classification complexity is; velocity sensor has worst complexity/value ratio
  but suppression gate is a required consumer. (by: user)
- 2026-08-07: Scope = full (harness + rebase + capture-mode + anchor-primary),
  anchor-primary last, gated on harness equivalence. (by: user)
- 2026-08-07: Active wrap→list bug: reproduce in harness first, no device run;
  rebase refactor closes it; then revert TEMP commit fb24f35. (by: user)
- 2026-08-07: Harness design (.claude/orchestrate/harness-design.md, 501fe36):
  widget-level attach through PUBLIC host via testWidgets; rig-atom + re-pumpWidget
  (no stateful shim); introspection = allRenderObjects + Dart field reads, zero
  render.cljd changes; RenderProxySliver probe for corrections; oracles O1-O9;
  LCG-seeded fuzzer, 14 ops, paste-ready replay vector. (by: harness-design agent)
- 2026-08-07: Fuzz dataset starts at 2200 items, drop to ~600 if baseline >60s;
  no :fast alias for now; no debug counters in render.cljd without re-opening
  the design. (by: coordinator)
- 2026-08-07: rebase-design open Q4 (:expect-ws slack) resolved per its default:
  one viewport (remainingPaintExtent) for all causes. (by: step-9a)
- 2026-08-07: rebase-design open Q5 (O9 corrected? input): commit-prune? keeps
  reading `(some? scrollOffsetCorrection)` off the FINAL geometry — with the
  no-emit capture that is exactly "the emit arm ran"; the mid-segment committed
  bound re-verifies in step 10's re-campaign as planned. (by: step-9a)
- 2026-08-07: pendingRebase arms only for :cross-layout (pre-emit, by
  update-render! — the stale-offsets fact originates at the mutation) and
  :landing (post-emit, :seed :init); seam/underflow corrected re-runs need no
  seed contract (cache prefix/dead-reckon serve them), so they arm nothing —
  minimal faithful mapping of the three deleted one-shots. (by: step-9a)
- 2026-08-07: segment-tail residual is a rebase producer the design did not
  enumerate: the host clock can complete before a t=1 pass runs, leaving
  desired(1) − segPrevDesired un-emitted; consume-seg-tail! publishes it via
  the accumulator on the first resting pass (segAnchor finally has an explicit
  end-of-life). Required for the far-morph-animated flip (O6 eps 0.5px vs
  ~5.8px tail on an 89k-px rebase). (by: step-9a)
- 2026-08-07: segAnchor flow `to` reads the capture-placed cache entry; the
  indexed target keeps the exact to-src frame (never nil in range); the
  exiting-anchor re-pick ships with segAnchor v2 in 9b's stage 4 per the
  design's stage map. (by: step-9a)
- 2026-08-08: rebase-design open Q2 (fraction set-point endpoint source) resolved
  per its default: `tw/point-desired` lerps the FROZEN endpoint extents; the
  segment tween is never read back, and the exiting-anchor re-pick covers the
  leave-slide case the alternative was meant to. (by: step-9b)
- 2026-08-08: rebase-design open Q1 (domain exit vs host clock) resolved per its
  default: the engine settles silently and the host clock plays out; no callback,
  no coupling. Step 10's re-campaign is the revisit trigger. (by: step-9b)
- 2026-08-08: segAnchor v2 keeps a `:screen` residual
  (`(from-off + frac·from-ext) − scrollOffset`). The design's identity
  "desired(0) == the segment-start scrollOffset" holds only when `frac` is
  unclamped AND the anchor was sampled at this pass's offset; the residual makes
  it exact for a clamped frac and for a drag landing between the sample and the
  segment start, preserving v1's drag composition. (by: step-9b)
- 2026-08-08: the viewAnchor epilogue (and the exactly-once rebase assert) gate on
  `segTween` nil, not `tweenAnim` nil — a domain-settled segment (§3.2) runs the
  resting drivers while the host clock plays out, and without resampling there the
  next segment starts from a pre-jump anchor whose capture rebase has no
  absorption channel (stage 3's assert fired on exactly that). (by: step-9b)
- 2026-08-08: §3.3's key re-anchor also applies to `segment-start!`'s set-point
  slot, which the design did not enumerate: a count change that OPENS a segment
  reaches segment-start! with a viewAnchor whose `:idx` is one behind, so the
  set-point aimed the viewport at the neighbouring item's slot and re-created the
  one-item slide the cache re-anchor had just prevented. (by: step-9b)
- 2026-08-08: O6's intra check is now the consumed FRACTION, not the absolute
  intra offset. Under a morph that resizes the anchor the two are incompatible,
  and holding the pixel offset IS NEW-1; with an unchanged extent they coincide,
  so no other scenario's strength changes. (by: step-9b)
- 2026-08-08: `settle-segment!` omits the design's `clips := {}` — both drivers
  clear `clips` at pass entry, so the write would be dead. (by: step-9b)

- 2026-08-09: the wrap seam's cross mismatch is NOT repairable by dropping the
  stale entries: `backfill-leading!` runs AFTER `walk!`, and a cross-only seam has
  no correction to force a re-run, so the pass would end with the window unplaced.
  The repair is the re-anchor ALIGNMENT (`:run-start`), which keeps the retained
  head's run and leaves the leading margin approximate. (by: step-9c)
- 2026-08-09: rebase-design open Q3 (count-change fallback when the anchor's child
  does NOT survive) resolved by measurement, not argument: the `:churn` profile —
  layouts starved, `:animate?` forced, `:remove-anchor` deleting the item the key
  re-anchor looks for — produced no churn-only signature in 40 episodes. Not a
  distinct defect class. (by: step-10)
- 2026-08-09: O9's mid-segment `cache-n` bound is left miscalibrated on purpose.
  Every candidate denominator that admits a landing's collapsed attached window
  also admits NEW-13's `pass-materialized 403`; the honest fix is a separate
  per-pass WORK probe, deferred as oracle design. (by: step-10)
- 2026-08-09: layouts gain an optional `:run-start (fn [env state] double)` — the
  next independent run's main offset. `:flow-end` answers where the next CHILD may
  land, which inside an open run is the run's own start, so it cannot express the
  inter-run gap that BACKWARD reasoning needs (leading dead-reckon, seam
  alignment). Defined for wrap; list's `:flow-end` already is its run start and
  masonry keeps the fallback. (by: step-9c)

- 2026-08-09: NEW-9's durable rule (step 10b): (a) viewAnchor samples the
  DISPLAYED frame after EVERY pass, mid-segment included — the "pre-segment
  truth" gate is wrong the moment the user scrolls during a live segment, and
  the displayed frame is the same frame `committed` records, so the anchor and
  the next segment's `from` capture stay coherent; (b) the :absorb arm is
  TOTAL: when no anchor resolves, a produced capture rebase is absorbed as a
  pure shift over the first placed slot (desired telescopes to exactly the
  rebase). The emit-from-capture alternative was rejected — a correction would
  rebase scrollOffset while the tween's `from` side stays in the old frame;
  forcing the capture window to cover a stale anchor slot would re-open the
  O(gap) walk. (by: step-10b)
- 2026-08-09: O6's subject is the BEFORE key's own child, never the reported
  top child: in a multi-column target the viewport top is a SET of row-mates
  sharing one span, and `top-anchor` reports its lowest index. Both campaign-2
  :o6 key-moved classes (NEW-10, NEW-11) were the engine holding the anchor's
  fraction point exactly while O6 compared representative keys. (by: step-10b)
- 2026-08-09: O5's extent-below-content bounds its content scan to children
  starting within pixels + 2 viewports: a far-moving segment pre-materializes
  its END window at target offsets past the lerping extent by design (NEW-12's
  cross-change segment was healthy — the set-point followed the anchor row into
  the re-packed grid). At rest the horizon covers the whole retained band, so
  the step-5 truncation class keeps its oracle. (by: step-10b)
- 2026-08-09: the measured aggregates are ENGINE-level truth, not a flow-driver
  cache: indexed resting passes feed them too, and the gate is `segTween` nil
  (9b's lesson) — a domain-settled pass under a playing clock must feed the EMA
  or the next capture's `:approx-offset` has no size basis. A basis-less
  estimate (`:avg-main` absent) is rejected by the cross-layout re-anchor:
  spacing-only arithmetic is not an estimate (NEW-13 link 2). (by: step-10b)
- 2026-08-09: the far-window seed is decided by GEOMETRY alone, never by what
  happens to be attached: pre-gc! can only drop a prefix, so cacheless
  leftovers (kept leave-slides past the band) survive it and must not veto the
  inverse seed — cache is an accelerator, never the algorithm (NEW-13 link 3).
  (by: step-10b)
- 2026-08-09: masonry defines `:run-start` = max(heights) — per the contract's
  own words ("what :anchor would be seeded with there"), the level-column
  anchor's next run line is the furthest column end. Fixes the leading-step
  run-gap drop and gives refine-seam's cross-stale arm the alignment repair 9c
  left degraded for masonry. (by: step-10b)

## Open questions
- ~~Step 11 is blocked on step 10's NEW-9~~ **RESOLVED (step 10b)**: NEW-9 is
  closed; the absorption channel exists whenever a capture produces a rebase
  (anchor resolution first, pure-shift fallback second), so capture-mode's
  "rebase returned as a value" premise holds. Campaign-3 shows 0 recurrence.
- Residual campaign-3 classes (fuzz-findings.md §Campaign 3) are the next
  triage set: the o6 above-window count-change family (7 seeds, the largest
  un-triaged engine family), two o5 remove-op extent seeds, and the
  deliberately-open O9 denominator. None block step 11. (raised by: step-10b)

## Checklist
- [x] 1. Explore test infra + example app usage of collection — agent: Explore, model: sonnet
- [x] 2. Design viewport harness API — agent: Plan, model: opus
- [x] 3. Implement harness core — agent: general-purpose, model: opus
- [x] 4. Green baseline scenario tests (scroll/jump/rotate/morph basics) — agent: general-purpose, model: sonnet
- [x] 5. Red repro: far-scroll wrap→list morph + correction-only capture extent — agent: general-purpose, model: opus
- [x] 6. Invariant fuzzer over random op sequences — agent: general-purpose, model: opus
- [x] 7. Triage fuzz findings → additional red tests — agent: general-purpose, model: sonnet
- [x] 8. Design pending-rebase structure (exactly-once protocol) — agent: general-purpose, model: fable
- [x] 9a. Rebase stages 0-3: WIP preserve+revert, viewAnchor, transport, no-emit capture — agent: general-purpose, model: fable
- [x] 9b. Rebase stages 4-6: fraction set-point, window sanity/domain, count-change — agent: general-purpose, model: opus
- [x] 9c. Wrap kernel fixes NEW-2/NEW-3 (design stages 8-9) — agent: general-purpose, model: opus
- [x] 10. Fuzz re-campaign (design stage 7); verify reds green; revert TEMP fb24f35 — agent: general-purpose, model: sonnet
- [x] 10b. Close campaign-2 defects (NEW-9..NEW-13 + masonry kernel) — agent: general-purpose, model: fable
- [ ] 11. Capture-mode: design + implement (corrections statically off, rebase returned as value) — agent: general-purpose, model: opus
- [ ] 12. Design anchor-primary migration (truth = anchor idx+intra offset; equivalence criteria; kill-list of subsystems) — agent: Plan, model: fable
- [ ] 13. Implement anchor-primary in stages behind harness equivalence — agent: general-purpose, model: fable
- [ ] 14. Cleanup: remove dead subsystems, update docs/CollectionRectAnimator.md + memory — agent: general-purpose, model: sonnet

## Step results

### 10b. Campaign-2 defect closure
- Status: done (commits ae01352, a2365e1, dcd3896, 0e22b24, 93ccf4a — one per
  defect; NEW-10 and NEW-11 share one commit, one shared root: the O6
  representative-key gap).
- Files changed: render.cljd (viewAnchor every-pass sampling, total :absorb
  fallback, EMA gate + indexed feed, basis-less estimate guard, geometric far
  seed), layout.cljd (masonry `:run-start`), harness.cljd (O6 before-key check,
  O5 horizon), fuzz_red_test.cljd (5 flips + docstring), render_test.cljd
  (masonry flip, refine-seam masonry alignment pins), fuzz-findings.md
  (per-class dispositions + Campaign 3), this plan.
- Suite: default `-x "known-red || fuzz"` 312 → **318** green after the final
  commit (each intermediate commit green: 313/315/316/317/318); O7 in-suite
  each run; `-t known-red` 6 → **0** (empty — every campaign-2 red closed);
  every flip verified 3x standalone; `bin/check` clean.
- **NEW-9 (engine, FIXED)** — the exactly-once breach was anchor STALENESS, not
  resolution: the cache invalidator opens a segment, the motion op drags the
  viewport mid-segment where the §2.1 gate froze viewAnchor, and the morph's
  capture consumed an anchor ~1200px behind the window (traced: va {:idx 0}
  against capture 20..32, rebase 840.5 discarded). Rule in the Decisions log:
  displayed-frame sampling every pass + a total :absorb arm (pure-shift
  fallback). Fixing staleness at the source also repairs the OTHER viewAnchor
  consumer (seed-cache!'s cross-layout re-anchor aims at the anchor slot).
- **NEW-13 (engine, FIXED, root ≠ NEW-9)** — a three-link starvation chain
  (EMA gate on tweenAnim; indexed passes never feeding aggregates → wrap
  estimate degenerating to 30px for index 234, trusted by the re-anchor;
  attachment vetoing the far-window inverse seed). Each link is its own
  standing-principle violation; details in fuzz-findings.md.
- **NEW-10 + NEW-11 (oracle, RECLASSIFIED)** — traced healthy: the set-point
  resolved the anchor's own slot (by key after the count change) and held its
  fraction point to 1e-15; the "moved"/"index held" keys were lower-index
  row-mates of the anchor in the re-packed grid row. O6 now checks the
  before-key's own child (strict subset of the old failure set — no masking).
- **NEW-12 (oracle, RECLASSIFIED)** — the cross-change segment was following
  the anchor row correctly; O5 counted the segment's designed end-window
  arrivals as under-advertised content. Horizon-bounded now; the guard test
  settles to verify the honest resting total.
- **Masonry kernel (FIXED)** — `:run-start` = max(heights); the refine-seam
  test's nil-alignment pin became value assertions (uneven prefix aligns its
  tallest column above the head's line).
- **Campaign-3 recurrence run**: full 220 episodes. Fixed classes: 0
  recurrence (NEW-9 0/17; 241/207, 16, 4/140, 37/29 all clean). 16 failing
  batches / 20 distinct seeds remain (campaign-2: 43 episodes), all
  pre-existing shapes: the deliberately-open O9 denominator (window-scale
  cache-n/committed-n mid-segment), the o6 above-window count-change family
  (seed-55/146 shape, 7 seeds — largest un-triaged engine family), two o5
  remove-op extent seeds, boundary-spring o5 singles, seed-251's sliver.dart
  assert. Recorded in fuzz-findings.md §Campaign 3.
- Next-step impact: **step 11 is unblocked** — the rebase-as-value premise
  holds (channel guaranteed). The o6 count-change family is the natural next
  triage before/with step 12's anchor-primary design, whose key re-anchor is
  exactly the machinery those seeds stress.

### 10. Fuzz re-campaign + verification + TEMP cleanup
- Status: done (commits f609581, 1586454, + this step's test/doc commits)
- Files changed: render.cljd (TEMP revert only — no engine fix); fuzz_test.cljd
  (`:churn` profile, `:remove-anchor` op kind, 10 new batch deftests, two oracle
  recalibrations), harness.cljd (`remove-anchor!` + its `apply-op!` arm),
  fuzz_red_test.cljd (5 new reds + `:churn` in `profile-of`), render_test.cljd
  (masonry leading-gap red), fuzz-findings.md (Campaign 2 section).
- Suite: default `-x "known-red || fuzz"` **312** green (unchanged — every new
  test is tagged). `-t known-red` **0 → 6**, all six added by this step.
  `bin/check` clean. Campaign runtime 74 s of `flutter test`.
- **Campaign 2**: 220 episodes (22 batches x 10 seeds) over four profiles —
  `:full` 1-60, `:no-layout` 101-160, `:no-jump` 201-260 (campaign-1 ranges kept
  as ratchets, plus 20 fresh seeds each) and the new `:churn` 301-340.
  55 failing episodes raw, **43** after two fuzzer recalibrations, in **7**
  classes. Campaign 1 was 62/120. known-A, known-B's plain form and NEW-1..NEW-5,
  NEW-7 did NOT return; 9c's `:no-jump` `:o4 wrap-run-pitch` watch item did not
  reappear either.
- **STOP-worthy — NEW-9, `capture produced rebase N with no absorption channel
  (segAnchor nil)`.** 17 seeds, ALL FOUR profiles, reduced by six independent
  seeds to the same three ops: a cache invalidator (count or cross change) → a
  motion op → a morph. This is stage 2/3's OWN exactly-once contract failing, so
  the invariant step 9 was built to establish does not hold. `viewAnchor` being
  non-nil (9a stage 1) does not make `segAnchor` RESOLVABLE: traced state has
  `viewAnchor {:idx 242}` against `cacheFirst 250`, so the frozen `to` lookup
  answers nil and the rebase is discarded — step 5's link (c) in a new dress. In
  release the assert compiles out and the content simply does not follow.
  **Blocks step 11**: capture-mode returns the rebase as a VALUE, which presumes
  the channel this class proves absent. Not fixed here (step 10 is verification).
- Other engine classes, red-tested, each verified standalone first:
  NEW-10 (`:o6`, 4 seeds) a morph preserves the consumed fraction to 1e-15 and
  applies it to the anchor's NEIGHBOUR — 9b's set-point is right, its slot
  resolution is off by one. NEW-11 (`:o6`, seed 16) known-B's shape returns when
  a far jump precedes the above-window insert; the re-anchor holds the index, not
  the key. NEW-12 (`:o5`, 6 seeds) a cross change advertises a scrollExtent below
  the content of the same frame — stage 3 closed this for the flow capture path
  only. NEW-13 (`:o9`, 5 seeds) a morph after a far jump caches ~370-403 of 600
  items with `pass-materialized` equal to it — NEW-7 minus the count-change
  adjacency its diagnosis relied on, so the stage-5 `:domain` settle misses it.
  Recorded without reds (single seed, long vector): `:o4 hole` seed 110,
  `:o6 intra-drift` seed 55, a Flutter `sliver.dart` assert at seed 251.
- Oracle miscalibrations fixed (both fuzzer-side, both traced before changing):
  (1) `boundary-settle!` ran only after jump ops, so a drag/fling released past a
  boundary was sampled mid-spring — 10 false `:o5`; it now runs after every op.
  (2) the O6 anchor guard fired on removals taken with `pixels` resting ON
  `maxScrollExtent`, where the shrinking range legitimately drags the viewport —
  3 false `:o6`; `:remove`/`:layout` are unguarded there, `:insert` still guarded.
  **Left alone deliberately**: O9's mid-segment `cache-n` denominator is
  `8*attached+64` with attached collapsed to 1 by a landing — the same problem
  step 7 fixed for `committed-n` — but every candidate loosening also masks
  NEW-13's genuine `pass-materialized 403`. The fix is a separate per-pass WORK
  probe, which is oracle design; recorded in fuzz-findings.md.
- **rebase-design open Q3 answered (negative)**: the `:churn` profile — layout
  starved, `:animate?` forced, `:remove-anchor` deleting the very item the key
  re-anchor looks for — produced 6 failures in 40 episodes and NOT ONE
  churn-only signature. The count-change fallback when the anchor's child does
  not survive is not a distinct defect class; mid-segment count changes surface
  NEW-9, the same class every other profile does.
- 9c's latent masonry gap is now pinned: `known-red-leading-step-masonry-backs-off
  -by-the-run-gap` (pure kernel test) — masonry has no `:run-start` and its
  `:flow-end` answers the head, so `leading-step-entry` returns offset 950.0 where
  944.0 is correct, dropping `:main-axis-spacing`. No harness oracle covers it:
  the estimated cell only TOUCHES the head, so O4 sees neither overlap nor hole.
- Next-step impact: step 11 must not start on NEW-9's premise. Either NEW-9 is
  closed first (a segAnchor resolution that cannot answer nil while a rebase is
  outstanding, or an absorption arm that does not need one), or step 11's
  rebase-as-value contract is re-designed around a channel that may be missing.
  NEW-10/11/13 all point at the same seam — slot/key resolution against a cache
  the same pass re-seeded — and are worth fixing as one group rather than singly.

### 9c. Wrap kernel fixes NEW-2/NEW-3
- Status: done (commits f3d8d54, 7899ef3 — one per fix, bisectable)
- Files changed: render.cljd (`refine-seam`, `reflow-from-checkpoint!`,
  `leading-step-entry`), layout.cljd (`:run-start` hook + contract docstring);
  render_test.cljd (+2 pure tests, `refine-seam-delta` call sites), fuzz_red_test.cljd
  (both reds un-tagged + ns docstring). TEMP fb24f35 instrumentation kept (step 10
  reverts).
- Suite: default `-x "known-red || fuzz"` 308 → **310** (stage 8) → **312**
  (stage 9), green after each; `-t known-red` 2 → 1 → **0** (no tests match — the
  known-red set is empty for the first time); both flips verified 3x standalone;
  `bin/check` clean at the end.
- Stage 8 (NEW-2, seam cross): `refine-seam` (was `refine-seam-delta`) returns
  `{:delta :cross-stale? :align-delta}`. **The design's fix shape was half the
  story**: the fuzz vector's overlap is produced by the FROZEN re-anchor, not by
  the stitch. `emit?` is false under the fast drag, so the seam re-anchors the
  prefix at `:anchor(lo, first-offset − d)` with `d` = where `:place` puts the
  head — and when the prefix's last run is OPEN that lands the prefix's run start
  ON the retained head's run (traced cf=304: d=−143.54, 303 re-anchored onto the
  head's 3572.45). So the branch order is: converged → stitch; resting with a real
  main delta → adopt/emit (its re-run re-walks the window canonically, cross
  included); cross-stale → re-anchor by `:align-delta` (the next RUN start), which
  keeps the head's run line free with no correction; else today's main re-anchor.
  Dropping the stale entries instead is NOT viable: `backfill-leading!` runs after
  `walk!`, and with no correction to force a re-run the pass would end with the
  window unplaced (verified by tracing the repro).
- Stage 9 (NEW-3, leading run advance): `leading-step-entry`'s back-off reads
  `(:run-start layout)` and keeps `:flow-end` as the fallback, so the cell lands one
  run PITCH above the head instead of flush against it.
- Not done (out of scope, latent): masonry's leading step drops
  `:main-axis-spacing` the same way — it has no run structure, hence no
  `:run-start`, and no oracle covers it. `refine-seam`'s cross check does fire for
  masonry, where it degrades to today's behavior (no `:run-start` ⇒ the stale seam
  is kept) rather than to a wrong alignment.
- Next-step impact: step 10's re-campaign should re-run the `:no-jump` profile
  (both vectors came from it) and watch `:o4` `wrap-run-pitch` — the frozen
  alignment deliberately leaves the leading margin one run gap "loose", which is a
  legal approximate margin but a new shape for the oracle. `-t known-red` is empty,
  so any new red is a genuine regression, not a leftover.

### 9b. Rebase stages 4-6
- Status: done (commits 2ae7abe, 04334dc, 4b0d614 — one per stage, bisectable)
- Files changed: render.cljd, tween.cljd (src); harness.cljd (top-anchor
  `:extent`/`:frac`, O6 fraction rule), tween_test.cljd, fuzz_red_test.cljd,
  harness_test.cljd (flip re-tags + one new red-turned-green), docs
  /CollectionRectAnimator.md (§7b set-point delta wording). TEMP fb24f35
  instrumentation kept (step 10 reverts).
- Suite: default `-x "known-red || fuzz"` 302 → **308** green after each stage
  (+1 pure `point-desired` test, +4 un-tagged reds, +1 new known-B fuzz vector);
  O7 (harness-smoke cold-vs-warm) green each stage; `-t known-red` 6 → **2**
  (`+0 -4`: NEW-1, NEW-5, NEW-7, known-B — exactly the stage 4-6 targets; NEW-2
  and NEW-3 remain, they are 9c's); `bin/check` clean at the end.
- Stage 4 (fraction set-point): segAnchor v2
  `{:from-off :from-ext :to-off :to-ext :frac :screen}`; `tw/point-desired` is the
  ONE formula and `point-correction-delta`, `consume-seg-tail!`'s residual and the
  leave-slide `shift` all route through it, so the 9a lockstep hazard is
  structural now, not a convention. `to-ext` comes from the same place as
  `to-off` (capture-placed cache entry for flow, `to-src` frame for indexed).
  Exiting-anchor re-pick (`repick-view-anchor`) takes the first surviving attached
  slot at/after the anchor, framed by its pre-segment committed rect. Flipped:
  NEW-1 (`fuzz-morph-holds-anchor-without-jump`, verified 3x).
  **O6 had to change with it**: `baseline-layout-morph-shallow-wrap-grid` and
  `-grid-masonry` went red on `intra-drift` because the anchor RESIZES there
  (wrap 109 → grid 128) and the fraction was preserved exactly — key and absolute
  intra offset cannot both hold. O6 compares `frac-before × extent-after`.
- Stage 5 (window sanity + domain): `keyed-tween-layout`'s window queries clamp a
  nil answer to the source's own edge (`first-index` → its last index,
  `last-index` → its base) instead of fabricating 0; `indexed-layout!` asserts
  non-inversion in debug and clamps in release. `segment-start!` stamps `:domain`
  = the union of the t=0 and t=1 window spans; the SEGMENT CONTINUE branch runs
  `segment-domain-exited?` (the `window-far-from-span?` kernel factored out of
  `window-far-from-band?`) and, on exit, `settle-segment!` drops
  segTween/segAnchor/leaving, sets segGen := curGen and dispatches to the resting
  drivers. Flipped: NEW-5, NEW-7 (both verified 3x).
  En route: NEW-7 first moved from `cache-n 601` to a stage-3 absorption-assert
  failure — the settled passes never resampled viewAnchor because the epilogue
  gated on `tweenAnim`. Gate is `segTween` now (Decisions log).
- Stage 6 (count change): `update-render!` drops cache + baseState and arms
  `{:cause :count-change :seed :anchor}`; `seed-cache!` resolves the anchor via
  `attached-index-of-key` (O(window)), GCs above it (NEW-4's guards intact) and
  anchors at {new-index, viewAnchor `:off`} with no rebase — the offsets stay in
  frame, only the estimate above the anchor is stale and the standard producers
  repair it. `segment-start!` resolves its set-point slot the same way (see the
  Decisions log — the design named only seed-cache!). Flipped: known-B
  (`insert-remove-above-window-anchor`, key-moved + intra-drift, verified 3x) and
  a new `fuzz-insert-above-window-holds-anchor` ratchet over known-B's own
  campaign vector (seed 10).
- Next-step impact: 9c (NEW-2/NEW-3) is untouched by all of this — both are seam
  kernels. Step 10's re-campaign should exercise (a) the domain-exit path under
  the `:no-layout` heavy-churn profile (design open Q3's count-change fallback is
  still unverified), and (b) mid-segment count changes, the one path where the
  cache drop + key re-anchor runs against tween-frame offsets.

### 9a. Rebase stages 0-3
- Status: done (commits 07c9ea3, 640daa3, 6832572, a3c2276, e540dab — one per stage,
  bisectable; WIP preserved as .claude/orchestrate/wip-seed-cache-reanchor.patch)
- Files changed: render.cljd; harness.cljd (engine-map keys), render_test.cljd
  (reanchor-band 0-sentinel), morph_red_test.cljd + fuzz_red_test.cljd (flip
  re-tags). TEMP fb24f35 instrumentation kept (step 10 reverts).
- Suite: default `-x "known-red || fuzz"` 297 → **302** green after each stage;
  O7 (harness-smoke cold-vs-warm) green each stage; `-t known-red` 11 → **6**
  (`+0 -6`: NEW-1/2/3/5/7 + known-B — exactly the stages 4-6/8-9 targets);
  bin/check clean at stages 1-3.
- Stage 0 (WIP revert): flipped NEW-4 → `fuzz-null-render-box-cast-in-morph-after-remove`
  (kept green as a guard against an unguarded re-anchor GC).
- Stage 1 (viewAnchor): field replaces restingTop + both capture paths
  (`anchor-before` resting write, `capture-resting-top!`); sampled in the
  performLayout epilogue from ATTACHED children in the corrected frame
  (`so' = scrollOffset + emitted correction`), every resting pass incl.
  re-seeding ones; non-nil-while-sized-children debug assert. Cross-layout
  re-anchor anchors on viewAnchor with the NEW-4 guards (skip GC that would
  empty the window; nullable firstChild re-read). Flipped:
  `view-anchor-survives-far-jump` (was resting-top-lost),
  `far-morph-wrap-to-list-not-animated`, and — **unexpected early flip,
  investigated** — `far-morph-capture-extent-truncated`: its oracle
  (mid-segment scrollExtent >= laid content) is a proxy; the restored set-point
  keeps the laid window inside the still-truncated snapshot extent. Un-tagged as
  `far-morph-capture-extent-covers-laid-content`; the extent root was fixed
  properly in stage 3.
- Stage 2 (passRebase + pendingRebase): all producers publish via `rebase+!`;
  ONE emission site in flow-layout!; `reanchor-band` reads the accumulator
  (0 = none). pendingRebase replaces crossLayoutReanchor/landingEmitted/
  pendingLandingWs; seed-cache! clears it unconditionally on entry; the
  consumed :landing marker feeds :violated detection and is wiped by the next
  seed — no disarm dance. Exactly-once epilogue assert (resting passes):
  emitted correction == produced rebase. No flips, by design.
- Stage 3 (no-emit capture): flow-layout! grew a two-valued mode axis
  (:resting/:capture — step 11 extends it); :capture never writes correction
  geometry and returns {:rebase :extent}. `flow-total-extent` factored out of
  finish-flow-geometry! supplies the frozen snapshot :extent from the TARGET
  layout (live-only-flow-window no longer reads .-geometry back). segAnchor
  flow `to` = the capture walk's placed cache entry (post-rebase frame);
  snapshot lookup only a fallback. Absorption asserts: capture geometry carries
  no correction; produced rebase with segAnchor nil fails loudly.
  **New defect found + fixed en route**: the segment clock can complete before
  a t=1 layout pass runs, so the set-point's un-emitted tail (~5.8px here)
  snapped content at segment end — `consume-seg-tail!` retires the lingering
  segAnchor on the first resting pass and publishes the residual through the
  accumulator (indexed driver honors the accumulator at rest via its correction
  exit). Flipped: `far-morph-wrap-to-list-animated` (verified stable 3x).
- Field-count delta: −5 (restingTop, crossLayoutReanchor, reanchorShift,
  landingEmitted, pendingLandingWs) +3 (viewAnchor, passRebase, pendingRebase)
  = net −2, as designed.
- Deferred to 9b per the design's own stage map: exiting-anchor re-pick (§2.5,
  stage 4), fraction set-point, tween window sanity/domain, count-change
  re-anchor.
- Next-step impact: 9b's stage 4 replaces segAnchor v1 {:from :to :screen} with
  the fraction record; `consume-seg-tail!`'s residual formula
  (`(to − screen) − segPrevDesired`) must be updated in lockstep with the new
  `desired` formula or the tail emission silently regresses.

### 1. Explore test infra
- Status: done
- Key findings:
  - Run: `clojure -M:cljd test [ns...]` (recompiles, then dart test; `--`/`++` pass dart-test flags). Never bare `flutter test`.
  - No live-tree harness exists; render_test.cljd (~1300 lines) tests pure kernels only. box_test.cljd detaches RenderCollectionBox via manual `.add` — does NOT transfer to the sliver (needs real RenderSliverBoxChildManager).
  - KEY: `cljd.test` deftest supports `:runner (ft/testWidgets [tester])` (clojuredart doc TESTING.md) — flutter_test is already a dev-dep; async bodies native. Lowest-friction route to a real RenderViewport (pumpWidget, tester.fling, jumpTo).
  - Manual repro rig for the active bug: `example/src/repl.cljd` `ws2-preview` (~L1139-1241): 2200 items, Wrap → "Jump ~50vp" → List ⇒ "top element #1700, nothing above, clamps at 0".
  - Tests are plain clojure.test style; compiled artifacts in test/cljd-out (generated).
- Next-step impact: design fork resolved toward testWidgets/widget-level harness; fake child manager route documented as fallback.

### 2. Design viewport harness API
- Status: done
- Files changed: .claude/orchestrate/harness-design.md (commit 501fe36)
- Key findings:
  - Tree: Directionality → MediaQuery → Align → SizedBox(cross×main) → KeyedSubtree(cold-gen) → custom-scroll(:clamping) → sliver-collection; Align mandatory (pumpWidget tight constraints); rotation = SizedBox width change.
  - Items {:id n}, sizes from id (not index), text-free. Every op = rig-atom write + re-pumpWidget so didUpdateWidget fires like in the app.
  - CLI ns narrowing narrows COMPILATION only — use -t/-x/-N to narrow the run; helper ns must carry a smoke deftest; no deftest timeout ⇒ split fuzz tests + `--timeout 4x`.
  - O8 (correction convergence) depends on RenderProxySliver probe; degrades to O1+O3 if unworkable.
- Next-step impact: steps 3-6 implement strictly from harness-design.md.

### 3. Implement harness core
- Status: done
- Files changed: test/flutter_cljd/internal/collection/harness.cljd (new),
  .claude/orchestrate/harness-design.md (run-command correction). Zero src/ changes.
- Suite: 283 pass (baseline 281 + 2 smoke tests). Harness run ≈6 s.
- Deviations from the design (all verified against Flutter 3.38):
  - Run flags: `flutter test` rejects `-N`; use `--plain-name <substring>`.
    `-t`/`-x`/`--timeout` unchanged. Design §7 amended.
  - `scroll-by!` uses a raw gesture (down / one moveBy / up), not `dragFrom`.
    The scrollable owns the ONLY arena member, so it is accepted at pointer-down
    with a zero pending offset and no touch slop is lost — the drag lands exactly.
    `dragFrom` would silently eat 20 px.
  - Masonry args are `:cross-axis-spacing`/`:main-axis-spacing`; the design's
    `:spacing`/`:run-spacing` are ignored by `masonry-layout` (would have meant 0 gaps).
  - O4 `:grid`/`:masonry` check equal columns derived from the OBSERVED cell size
    instead of re-deriving the grid solver's track math (same strength, no coupling).
  - O4 wrap run-advance skips run 0 (leading GC can cut it, so its max main is unreliable).
  - `check!` is async (O3 pumps a frame); `check-report!` also returns `:ran`/`:skipped`.
  - O7 cold-restarts the SAME rig (one tester ⇒ one root); documented as destructive.
  - §1.3's "no-op rebuild passes the same vector instance" is unreachable: the host
    does `(vec source)`, which always allocates, so `data-changed` is true on every
    rebuild — in the app too. Harmless (`pump!` never rebuilds), but the corollary buys nothing.
- Confirmed working: probe sliver logs corrections (O8 non-vacuous), all four
  layouts pass O4 at rest, `advance-segment!` walks a live wrap→list morph, O7
  warm-vs-cold agrees, `replay!` round-trips an op vector.
- Next-step impact: the animated smoke test already trips the TEMP fb24f35
  instrumentation (`reanchorShift=1028.0 ... snap-extent` on a wrap→list morph),
  so step 5's red has the right code path under the harness.

### 4. Green baseline scenario tests
- Status: done
- Files changed: test/flutter_cljd/internal/collection/harness_test.cljd (new,
  11 deftests). Zero src/ changes (the pre-existing uncommitted WIP diff on
  render.cljd from an earlier session was left untouched and unstaged).
- Suite: 294 pass with `-- -x known-red` (283 baseline + 11). Full run
  (compile + `flutter test`) ≈ 28-30s, inside the <60s budget.
- Coverage: cold start × 4 layouts, small scrolls both directions, deep fling,
  far jumpTo both directions + back to 0, exact jump-to-top landing, cross
  (rotation) change at rest ×2 and at depth, insert/remove/swap/rotate at rest,
  layout morph at shallow depth × all 4 combos (list→wrap, wrap→grid,
  masonry→list, grid→masonry).
- Known-red (2 deftests, `:tags [:known-red]`, excluded via `-x known-red`,
  reproduced 3x in isolation each — real engine gaps, not harness defects):
  - `known-red-insert-remove-above-window-anchor` — `update-render!`
    (render.cljd L963-975) resets checkpoints on an item-count change but
    never clears the per-index geometry cache; inserting/removing above an
    already-materialized window leaves the cache answering for the OLD item
    at that index, so the viewport-top key silently shifts by the count delta
    (O6 fails: key 354→353 with index and intra-offset unchanged).
  - Investigated but NOT kept red: a `:grid`→`:masonry` shallow morph inside a
    shared 4-combo deftest failed once (O6), but reproduces standalone 3/3 as
    GREEN — root cause was cross-rig pollution (undisposed controllers/tickers
    across sequential `make-rig` calls sharing one `testWidgets` tester within
    one deftest, not the engine). Fixed by giving every independent scenario
    its own deftest (also applied to the other 3 morph combos and to
    data-mutation, so a red assertion never masks untested siblings —
    `is` throws in this cljd port, aborting the rest of the deftest body).
- Harness changes: none needed; zero defects found in harness.cljd itself.
- Next-step impact: the insert/remove-above-window cache-staleness gap is a
  distinct architectural finding (checkpoints/cache invalidation split) worth
  folding into step 9's rebase-unification scope, alongside the wrap→list
  far-morph bug step 5 targets.

### 5. Red repro: far-scroll wrap→list morph
- Status: done
- Files changed: test/flutter_cljd/internal/collection/morph_red_test.cljd (new,
  4 deftests, all `:tags [:widget :known-red]`). Zero src/ changes; zero harness
  changes (no harness defect surfaced). The pre-existing uncommitted WIP on
  `seed-cache!` (restingTop re-anchor) was left in place — **the repro includes
  that WIP** and, per the diagnosis below, the WIP is inert in this scenario.
- Suite: `-- -x known-red` 294 pass (unchanged); `-- -t known-red` 5 fail
  (the 4 new + step 4's insert/remove anchor red). Reproduced 2× identically.
- Reds added:
  - `known-red-far-morph-wrap-to-list-animated` — the full recipe. Fails on the
    settled oracle battery (O1: engine tripwire) and on O6 (anchor 1346 → 342).
  - `known-red-far-morph-capture-extent-truncated` — mid-segment `scrollExtent`
    falls below the content the same frame laid out.
  - `known-red-resting-top-lost-after-far-jump` — root cause, engine-state level:
    a drag leaves `restingTop` set, a far jump leaves it nil.
  - `known-red-far-morph-wrap-to-list-not-animated` — plain layout swap at depth;
    O6 only (anchor 1346 → 1340).
- Reference numbers (deterministic, 2200 items, cross 400, wrap offset 30115,
  wrap anchor {:idx 1346 :intra 31}):
  - capture pass: `reanchorShift=88778.079`, `cacheFirst=1337`, `snap-base=1337`,
    `snap-frames=13`, `snap-extent=119744.079`; true resting list total 196259.419.
  - `resting=nil`, `seg-anchor=nil`; `pixels` stays 30115.0 for the whole segment.
  - mid-segment `scrollExtent` 76616 / 101795 / 113259 / 117568 / 119046 against
    `content-end` 118663 — 4 of 5 frames advertise less than they laid;
    `maxScrollExtent` bottoms at 76016 (correct value 195659).
  - first resting pass after the segment trips
    `assert-materialization-bounded!`: "flow layout [:list 6.0] laid out 999
    children in ONE pass — budget 92.98 for the 1265.4px band".
- Confirmed diagnosis (hypothesis links (a)/(b)/(c) from the code map):
  - **(c) CONFIRMED — dominant link, but for a different reason than hypothesized.**
    `seg-anchor` is nil not because the slot's target frame is off the snapshot
    window, but because `restingTop` is *itself* nil. `flow-layout!` samples it
    from `anchor-before`, which reads `.-cache` **before this pass's `walk!`
    refills it**; every re-seeding pass therefore records nil, and the last pass
    of a far jump always re-seeds (`inverse-seed!` and `backfill-leading!`'s
    top-underflow branch both do `(.-cache! rs [])`). Nothing re-lays while idle,
    so nil survives until the next real drag. Empirically: after `scroll-by!`
    restingTop = {:idx 1352 …}; after `jumpTo` + settle + extra pumps it stays
    nil. That is exactly why the device recipe needs *jump then morph*.
    Consequence: no set-point ⇒ no scroll correction for the whole segment ⇒
    cells settle ~88 000 px below the viewport, which is the "top element #1700,
    nothing above it" symptom.
    Note the contrast: the indexed driver captures the anchor from the ATTACHED
    CHILDREN post-layout (`capture-resting-top!`); the flow driver captures it
    from the pre-walk cache. Same concept, two implementations, one of them lossy.
  - **(b) CONFIRMED.** The capture pass ends on the reanchorShift correction-only
    branch (`flow-layout!` L2401-2405) and `from-relay!` overwrites that geometry
    — the first segment frame shows `offset` unchanged at 30115 with
    `:reanchor-shift 88778.079` still sitting in the field. The rebase survives
    only via the seg-anchor set-point, which (c) has already nulled.
  - **(a) CONFIRMED with a caveat.** `live-only-flow-window`'s `shadow-ext` reads
    `(.-scrollExtent g)` off that correction-only geometry ⇒ 0.0 ⇒ `:extent`
    degenerates to `frames-main-extent` = the window end (119744) instead of the
    list total (196259). Caveat for step 9/11: fixing the zero read is NOT
    sufficient — the honest `shadow-ext` there is the OLD wrap total (49837),
    still below the window end. The snapshot has no access to the TARGET
    layout's `:max-extent`; that must be supplied explicitly.
  - **New finding — the resting path is not healthy either.** Without `:animate`
    the morph re-anchors on `firstChild` (index 1336 = the leading overscan head)
    rather than the viewport-top item (1346), so content slides by the overscan:
    anchor 1346 → 1340, intra 31 → 46. No clamp, no extent collapse, no tripwire,
    O2/O3/O4/O5 all clean — a small, isolated defect, not the catastrophic one.
    So the fix is NOT scoped to the segment protocol alone: `seed-cache!`'s
    anchor choice is independently wrong. The uncommitted WIP targets precisely
    this, and is inert here because the `restingTop` it reads is nil per (c).
- Next-step impact (step 5): step 8's pending-rebase design must own three things the
  current code splits: (1) a single anchor-capture point that is valid after ANY
  pass (including re-seeding ones), (2) a rebase that survives `from-relay!`
  instead of racing it, (3) a target-layout total extent carried into the frozen
  segment snapshot. Step 6's fuzzer should include `[:jump …] [:layout …]` with
  no intervening drag — the ordering that makes the anchor nil.

### 6. Invariant fuzzer over random op sequences
- Status: done
- Files changed: test/flutter_cljd/internal/collection/fuzz_test.cljd (new,
  13 deftests), .claude/orchestrate/fuzz-findings.md (new),
  test/flutter_cljd/internal/collection/harness.cljd (one fix, below).
  Zero src/ changes; the uncommitted WIP on render.cljd left untouched/unstaged.
- Suite: `-- -x "known-red || fuzz"` 295 pass (unchanged). NOTE: two `-x` flags
  do NOT compose — the last wins; pass one boolean selector.
- Campaign: 120 episodes (12 batches x 10 seeds, 3 weight profiles), 24 ops each,
  600 items. 62 failing episodes, 15 raw oracle signatures, **10 distinct classes:
  2 known + 8 NEW**. ~59 s of `flutter test` time. Shrinking (prefix cut + greedy
  removal, <=24 replays) reduced most failures to 1–6 ops.
- Design deviations (recorded in the ns docstring): three weight profiles
  (`:full` = design §6.2 verbatim, plus `:no-layout` / `:no-jump` so the two known
  bugs stop masking the rest of the op space); a batch keeps going after a failing
  episode and reports every seed at once (`is` throws here); prefix binary search
  dropped — `replay!` already reports the failing op index; `check!` cannot run O6,
  so the fuzzer captures the anchor itself around layout swaps and strictly
  above-window inserts/removes.
- Headline NEW signatures (full list + shrunk vectors: fuzz-findings.md):
  - NEW-4 `type 'Null' is not a subtype of type 'RenderBox'` — an outright crash
    in a morph pass that follows a remove. Highest severity.
  - NEW-5 Flutter's `childCount >= leadingGarbage + trailingGarbage` assert after
    cross change + `to-top` on masonry.
  - NEW-1 a morph loses the anchor with NO far jump (3 ops: settle, drag, layout)
    — the step-5 reds all need a ~50-viewport jump, so this is uncovered.
  - NEW-2/NEW-3 wrap tiling: runs physically overlap after a backward drag; run
    advance short by one run's max main after a morph.
  - NEW-7 count change + far jump walks the WHOLE dataset in one pass
    (`cache-n 601`) — same root as known-B, violates the O(window) invariant.
  - NEW-6 a jump past the end parks `pixels` ~20 px outside `maxScrollExtent` and
    never converges (14 seeds, 1-op repro).
  - NEW-8 (17 seeds) `committed-n` above O9's window-proportional bound —
    flagged PLAUSIBLE; may be O9 miscalibration for a landing window.
  - known-A reduced to TWO ops: `[[:jump 16572.2] [:layout :list]]`.
- Harness fix (genuine defect): `make-rig` reused the previous rig's element
  subtree — its `KeyedSubtree` key was always 0, so a second rig under one tester
  inherited the first's `SliverCollectionState` and render-object cache. Step 4
  worked around this by splitting deftests; a fuzzer cannot. Now seeded from a
  monotone `rig-gen`. Two fuzzer-side attribution fixes are listed in
  fuzz-findings.md §"Harness / fuzzer defects".
- Method note for step 7: every reported signature was re-verified as the FIRST
  episode of its own deftest (12/12 reproduced). Two candidates that failed that
  check were deferred exceptions leaking from the previous episode — always
  verify standalone before calling something an engine bug.
- Next-step impact: step 7 has 8 NEW signatures with paste-ready vectors; NEW-4
  and NEW-5 are crash-class and should lead. NEW-1 widens step 8/9's scope: the
  anchor-capture defect is not gated on `restingTop` being nil.

### 7. Triage fuzz findings → red tests
- Status: done
- Files changed: test/flutter_cljd/internal/collection/fuzz_red_test.cljd (new,
  8 deftests), harness.cljd (O9 recalibration), fuzz_test.cljd
  (`boundary-settle!` 45 → 150 frames), .claude/orchestrate/fuzz-findings.md.
  Zero src/ changes; the uncommitted WIP on render.cljd left untouched/unstaged.
- Suite: `-- -x "known-red || fuzz"` 297 pass (295 + the 2 reclassified greens);
  `-- -t known-red` 11 fail (4 morph + 1 harness_test + 6 new);
  `--plain-name fuzz-red-test` is `+2 -6`, identical across 3 runs.
- Method: each red replays its shrunk vector through `fuzz-test/run-ops!` (the
  driver that found it) with the rig header REGENERATED from the seed, so
  :animate?/:approx/:layout0/:cross0 cannot drift from the generator.
- Dispositions:
  - NEW-1 red `known-red-fuzz-morph-loses-anchor-without-jump`
  - NEW-2 red `known-red-fuzz-wrap-runs-overlap-after-back-drag`
  - NEW-3 red `known-red-fuzz-wrap-run-advance-short-after-morph`
  - NEW-4 red `known-red-fuzz-null-render-box-cast-in-morph-after-remove`
  - NEW-5 red `known-red-fuzz-garbage-counts-exceed-child-count`
  - NEW-6 reclassified (fuzzer wait too short) → green `fuzz-jump-past-end-converges-into-range`
  - NEW-7 red `known-red-fuzz-count-change-then-far-jump-walks-dataset`
  - NEW-8 reclassified (O9 denominator) → green `fuzz-committed-stays-band-proportional-mid-segment`
- Confirmed diagnoses (full text + evidence in fuzz-findings.md):
  - **NEW-5 and NEW-7 are ONE root**: `keyed-tween-layout`'s
    `:first-index`/`:last-index` map "the frozen snapshot has no opinion at this
    offset" to index 0 (tween.cljd L352-357). Viewport BELOW the snapshot ⇒
    last-index 0 with first-index inside it ⇒ an inverted window whose two
    garbage walks count the same children twice (NEW-5). Viewport ABOVE it ⇒
    first-index 0 ⇒ the band stays pinned at the top for the whole segment while
    `pixels` moves 16 000 px, and the next layout swap re-flows the whole gap in
    one pass (NEW-7). NEW-7 is NOT known-B's root — `update-render!` does clear
    `anchoredTo0` on a count change; the `[:remove]` matters only because it
    opens a segment.
  - **NEW-1 is NOT known-A's root**: `segAnchor` is present and correct, but the
    segment shift is the rigid `to − from`, which pins the anchor cell's TOP EDGE
    to a screen position. A morph that resizes a partially scrolled-past anchor
    (grid 208 → list 109) therefore slides the viewport top into the next item.
    Step 8 must decide what the set-point preserves, not only where it comes from.
  - **NEW-4 is a defect of the uncommitted `seed-cache!` WIP, not of HEAD**: its
    leading `collectGarbage` can destroy every attached child, and the next line
    casts `.-firstChild` to a non-nullable `RenderBox` (render.cljd L1428-1430)
    before its own nil guards run. Step 5's "the WIP is inert" holds only for the
    far-morph scenario.
  - NEW-2: `refine-seam-delta` reconciles the backward re-flow on the MAIN axis
    only; a re-flow ending with an OPEN run yields d = 0, so `stitch-prefix`
    keeps the stale head at cross 0 and two adjacent indices share a run line.
  - NEW-3: `leading-step-entry`'s `adv = max(ext, flow-end(st') − head)` collapses
    to `ext` for wrap (`:anchor` seeds an empty run at `head`, `:place` puts i
    into it, so `flow-end` returns `head`) — the run gap is dropped.
- Next-step impact: step 8's pending-rebase design gains two consumers beyond
  step 5's three — (4) window queries must distinguish "no opinion" from index 0
  and the window must be refused when inverted, and (5) a far jump taken WHILE a
  segment runs must re-anchor the band instead of deferring to the frozen
  snapshot. NEW-1 adds a design question to the set-point itself. NEW-2/NEW-3 are
  wrap-tiling defects orthogonal to the rebase — schedule them separately.

### 8. Design pending-rebase structure
- Status: done
- Files changed: .claude/orchestrate/rebase-design.md (commit 1f09e71)
- Key findings: viewAnchor {idx key off extent frac} sampled post-layout every
  resting pass (never nil while children exist); passRebase per-pass accumulator,
  one publication site, three arms (:emit/:absorb/:value) + exactly-once debug
  assert; pendingRebase record replaces landingEmitted/pendingLandingWs/
  crossLayoutReanchor; capture extent via pure flow-total-extent (never reads
  .-geometry back); set-point preserves intra-anchor FRACTION; tween windows
  clamp to snapshot edge + segment :domain with engine-side settle; count change
  = drop cache + re-anchor by viewAnchor key; NEW-2/NEW-3 freestanding.
  Step 9 = 10 ordered stages, each default-green + O7-gated; fields 5->3.
- Next-step impact: step 9 split into 9a (stages 0-3, fable), 9b (stages 4-6),
  9c (stages 8-9); fuzz re-campaign (stage 7) folded into step 10. WIP diff
  preserved as a patch file before stage-0 revert.
