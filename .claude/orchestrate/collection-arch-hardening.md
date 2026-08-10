# Task: Collection engine architectural hardening

**Slug**: collection-arch-hardening
**Started**: 2026-08-07
**Status**: done (2026-08-10)

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

- 2026-08-09: the pass mode is a FOUR-value axis, not the design's three: the
  `tweenAnim`-vs-`segTween` split the engine had been re-deriving at every gate
  IS the `:settled` mode (§3.2's domain-settled pass — resting drivers running
  while the host clock plays out). Naming it retires the bug class 9b stage 5
  hit; the alternative (three modes plus a residual boolean) would have kept
  exactly the re-derivation step 11 exists to remove. (by: step-11)
- 2026-08-09: `:settled` deliberately keeps today's answers rather than the
  "natural" ones: refinement / origin corrections / overscan / commit prune stay
  OFF there, even though its emission arm is live and no set-point competes with
  it. Turning them on is a behavior change, not a consolidation, and step 11 is
  a consolidation; the table now makes the choice visible in one place instead
  of implying it from `tweenAnim`. Candidate for step 12/13 to revisit — a
  jump-to-top taken during a segment defers its exact landing until the clock
  ends. (by: step-11)
- 2026-08-10: anchor-primary Q1/Q3/Q6 resolved per their leanings, unexercised
  here (Q1's leadExact bool and Q3's 1e-3 hysteresis land with stage 3; Q6's
  far-jump frac 0 is untouched by stage 2). Q5 resolved AGAINST its leaning by
  measurement: the Delta-epilogue does not close the o6 family, and the family's
  root is the segment set-point, not the leading side. (by: step-13a stage 2)
- 2026-08-10: leading-extent measurements SUPERSEDE, they never accumulate. Once
  no step moves the window, every step compares the same rigid frame against its
  own reference, so each report is the WHOLE displacement restated one run closer
  to the origin — a k-seam chain reports bias, 2*bias, 3*bias. Summing them (the
  old telescoping sum, correct only because each shift realigned the frame)
  translates the window by a multiple of its own error. (by: step-13a stage 2)
- 2026-08-10: translate-window! drops the cache + checkpoints instead of
  rewriting them: a flow state is layout-opaque and carries absolute offsets, so
  only a new layout hook could translate one, and the design forbids extending
  the contract. Cost is nil — walk! lays every child every pass anyway; the cache
  only saves its :place arithmetic. The band's FRONTIER is kept and baseState
  re-derived at its translated offset, which is exact because the backfill leaves
  cacheFirst on a renewal point (re-seeding at the first attached child instead is
  mid-run and loses a run). (by: step-13a stage 2)
- 2026-08-10: reanchor-band and anchor-before/anchor-delta SURVIVE stage 2
  (the design listed them as stage-2 deaths). Both serve producers that lay the
  window in the CORRECTED frame by construction — the cross-layout re-anchor and
  the seg-tail before the walk, an in-window resize above the anchor during it —
  which is the opposite class from the epilogue's translate. They die with
  stage 3's anchor-seeded walk, where the anchor's in-pass delta is 0 by
  construction. (by: step-13a stage 2)
- 2026-08-10: 9a's "segAnchor flow `to` reads the capture-placed cache entry"
  SURVIVES, re-founded: the cache entry is not a seeding artifact but the frame
  the capture pass actually placed, and the tween's own `to-src` is now built to
  reproduce it. The stage-2 diagnosis inverted this — it read a constant
  `point-desired` as the set-point going silent, when the silence was right and
  the tween's target was the moving side. A `to` slot and a `to-src` frame that
  can disagree in-window is the bug; agreeing by construction is the fix, and
  which of the two the set-point reads then stops mattering. (by: step-13a stage 2b)
- 2026-08-10: a frozen segment must never lay a child it has no frame for. An
  attached child above the `to-src` window gets `zero-frame` — offset 0 — which
  strands firstChild at the scroll top and costs the next capture an O(gap) walk
  (`:o9` seed 340). The window sources are already nil-clamped (§3.2 layer 1); the
  per-child frame is not. (by: step-13a stage 2b)

- 2026-08-10: the per-child frame clamp is the SAME rule as the window clamp, not
  a new one: a query the source cannot answer resolves to the source's own edge.
  `parked-frame` probes those edges at offsets a bounded source must cover (its
  start and its own reported extent) rather than ±infinity — the ±inf form the
  window queries use would reach a layout's index math (`(.floor infinity)`,
  `getMaxChildIndexForScrollOffset(infinity)`) as soon as an INDEXED target could
  answer nil. Parking is expressed as a substitution of `to` before the
  enter/exit/leave cond, so all four branches inherit it instead of each growing
  a nil case. (by: step-13a stage 2c)
- 2026-08-10: a parked cell is collapsed on the LAYOUT's own axis, not the main
  axis: only the collapse axis has a `:full-*` slot, so zeroing anything else
  would lay the child's content out at extent 0 (§8a's own rule). For a
  `:cross` layout the parked cell keeps the edge row's main span — it is a
  row-mate, which is what sharing that slot means. (by: step-13a stage 2c)
- 2026-08-10: `:o6` 34 / 55 / 205 are NOT a boundary-clamp class — 34 and 205
  drift identically with the viewport nowhere near an edge. The class is the
  count-change re-anchor (anchor pinned, no rebase) against the leading
  re-measure the Δ-epilogue then translates the window by: at seed 55 the
  translation is 2991.8 while the anchor's own displacement is 2889.8. A sliver
  cannot clamp its own emission — `SliverConstraints` has `precedingScrollExtent`
  and nothing about what follows, so `maxScrollExtent` is not knowable here and a
  self-computed clamp would fire early in a composed viewport. Deferred to
  stage 3, where the anchor-seeded walk makes the anchor's in-pass delta 0 by
  construction. (by: step-13a stage 2c)

- 2026-08-10: a flow pass's frame is a run-chain from wherever its walk STARTS,
  so it is not reproducible from a different start. `translate-window!` empties
  the band, which moves the next pass's walk start to the frontier and re-packs
  the runs below it; `stitch-prefix` glues a re-anchored prefix onto a retained
  suffix that matches only in the head's OFFSET, so the chain the next pass
  re-derives is a different one. Nothing in the pass can prevent that — the
  anchor is what holds ACROSS it, which is why anchor-follow had to stop reading
  the emptied cache. (by: step-13b stage 3)
- 2026-08-10: re-deriving the window inside `translate-window!` (so the next pass
  would reproduce it) is stage 2b's rejected move in another costume: a forward
  re-flow seeded at the frontier levels a masonry's columns and fabricates a wrap
  run break. MEASURED: seed 34 and seed 52 both drift worse with it. The rigid
  translation stays; the anchor absorbs the re-derivation. (by: step-13b stage 3)
- 2026-08-10: the anchor reference must be an IDENTITY, not "the first cached
  entry overlapping the top". A bsearch reference re-picks the run's lowest
  row-mate every pass, so a re-pack that splits the run holds the run and drops
  the child O6 tracks. viewAnchor is that identity, re-resolved by key across a
  count change. (by: step-13b stage 3)
- 2026-08-10: a SYNTHESIZED checkpoint's seam is not a measurement — it sits at
  `:approx-offset`, which is a function of the measured EMA that this very window
  feeds. Honouring it under a latched anchor makes the frame chase itself one
  `1 − ema-alpha` = 0.75 step per pass, and the viewport's 10-cycle correction
  budget runs out first (seed 52 turned :o6 into :o1 that way). This is design
  §2.1's latching invariant, discovered from the other end: "aggregate drift never
  moves a latched lead". A pass that REASSIGNED the anchor still honours it —
  there the estimate is the truth it just chose. (by: step-13b stage 3)
- 2026-08-10: seeding EVERY `:anchor`-plan re-seed at viewAnchor's slot (design
  §2.2 step 2, the anchor-seeded walk) was MEASURED AND REJECTED: it breaks NEW-2
  (wrap runs stop tiling after a backward drag) because the anchor is generally
  mid-run and `:anchor` opens a fresh one there. The slot resolution stays scoped
  to the causes that renumber or re-space it. (by: step-13b stage 3)
- 2026-08-10: anchor-primary Q1 resolved AGAINST its leaning: no `leadEstimate`
  field. The latch it describes already exists as `viewAnchor :off` (sampled every
  pass) + `anchoredTo0`, and a second copy of the anchor's offset would be state
  to keep in sync, not single ownership — against the design's own "every stage
  should delete more state than it adds". Q3 (1e-3 hysteresis) confirmed: no new
  tunable was needed. (by: step-13b stage 3)
- 2026-08-10: an O(gap) walk is not always a seeding bug — `inverse-seed!`'s
  `off-start >= ws/2` guard is one, applied where it does not belong. The guard
  rejects a seed the walk would have to travel to; at the LAST index there is
  nothing to travel, so a window past the content end is O(1) from there whatever
  the estimate says. Its `:degenerate` fallback ("keep the band") contradicts the
  very plan that called it: `:far-inverse` was chosen BECAUSE the band is far.
  (by: step 13c)
- 2026-08-10: the leading edge of a pass is what it LAID, not what its index math
  says. `first-index`'s own frame is the leading edge only for a monotone source;
  a masonry column and a morph target are not, and the leading walk may stop short
  of `first-index` entirely — `calculateCacheOffset` then gets `from > to`. Same
  rule on the other end: the advertised extent is `max(total, laid trail)`, which
  is §2.4's honesty rule the flow driver already had and the segment drivers did
  not. (by: step 13c)
- 2026-08-10: the leading side of a correction IS clampable, unlike the trailing
  side (stage 2c). `SliverConstraints` carries `precedingScrollExtent`, so
  "content space is pinned at the leading edge" is a bound the sliver can compute:
  below it the viewport parks `pixels` out of range and — the sliver's own
  scrollOffset already clamped at 0 — no further pass runs to walk it back. The
  accumulator is trimmed before the emit arm reads it, so the exactly-once
  identity survives; the trimmed remainder is design §3's boundary rule (at the
  scroll top the content lands and the viewport shows index 0). (by: step 13c)
- 2026-08-10: a set-point built on an ESTIMATE is worse than none, but a set-point
  built on a frame the capture actually placed and then discarded is free. The
  frozen `to` snapshot starts at the capture window, so a morph that packs the
  anchor above it lost the set-point entirely; the capture walk had placed the
  slot, so it now returns that frame with its rebase and extent. When the walk did
  not cover the slot the value stays nil. (by: step 13c)
- 2026-08-10: closing the anchor's TARGET frame does not close the anchor's
  IDENTITY. Seed 235 keeps drifting because the capture materializes the window at
  the CURRENT offset: after the set-point moves the viewport there is no child at
  the new top, and the anchor re-samples onto the window's first child. The END
  window is only knowable after the set-point the capture produces — a chicken-
  and-egg the §7b/§7c machinery would have to break, out of this plan's scope.
  (by: step 13c)

- 2026-08-09: the INDEXED capture pass no longer feeds the measured EMA. It fed
  it iff no PREVIOUS segment tween happened to still be around — an accident of
  the `segTween` gate, and the opposite of the flow capture. A capture
  materializes a widened, transient window at frames that may never display;
  10b's "aggregates are ENGINE truth" is about DISPLAYED resting truth, which
  the surrounding resting passes supply for the same cells. (by: step-11)

## Open questions
- ~~Step 11 is blocked on step 10's NEW-9~~ **RESOLVED (step 10b)**: NEW-9 is
  closed; the absorption channel exists whenever a capture produces a rebase
  (anchor resolution first, pure-shift fallback second), so capture-mode's
  "rebase returned as a value" premise holds. Campaign-3 shows 0 recurrence.
- ~~Residual campaign-3 classes~~ **PARTLY RESOLVED (step-13a stage 2b)**: the O9
  denominator closed in stage 0; the o6 count-change family is fixed at its real
  root (the segment's `to` FRAMES, not the `to` slot — see stage 2b). Fuzz is 11
  seeds, a strict subset of stage 2's 14, so stage 3 proceeds on a better tree.
  **NEW-14 CLOSED (step-13b stage 3)**: `:o6` 34/55/205 and 52@16 plus `:o5` 7
  are green, the red test is un-tagged. Fuzz is 6 seeds. **Step 13c**: `:o1`
  56/251 and `:o5` 23/133 green; fuzz is **2 seeds**, `:o6` 235 + 337, both
  carried as documented reds (NEW-16/NEW-17 below).
  (raised by: step-13a stage 2, updated by: stage 2b, stage 2c, step-13b stage 3,
  step 13c)
- ~~**NEW-15** (raised by step-13b stage 3, seed 56)~~ **RESOLVED (step 13c)**,
  and NOT where stage 3 placed it: `align-start-fallback` never runs. Traced
  per-pass — the jump past the content end takes `seed-plan :far-inverse`,
  `inverse-seed!` saturates the inversion at the LAST index, and its
  `off-start >= ws/2` guard then rejects that seed and returns `:degenerate`,
  which KEEPS the far band. `align-start!` attaches at the band's end and the
  walk covers the whole gap (558 children, 557 budget). The guard now exempts the
  last index — nothing follows it, so the seed is O(1) however far the estimate
  undershot.
- **NEW-16 (raised by step 13c, seed 235, `:no-jump` op 12)**: a morph into a
  DENSER layout. The capture window is materialized at the CURRENT offset, so
  once the set-point moves the viewport to the anchor's new (much smaller) target
  offset there is no attached child at the new viewport top; `sample-view-anchor!`
  re-anchors onto the window's first child (25px below the top), and the first
  resting pass after the segment holds THAT one — +115px against the anchor the
  user was looking at. Closing it needs the capture to materialize the END window,
  which is only knowable after the set-point the capture itself produces. Red
  test `known-red-morph-loses-the-anchor-above-the-capture-window`.
- **NEW-17 (raised by step 13c, seed 337, `:churn` op 8, untraced)**: an insert
  above the window with `pixels` resting ON `maxScrollExtent` drifts the masonry
  anchor's consumed fraction by ~6px (intra 38.02 -> 32.00). Same shape as NEW-14,
  which step 13b closed for 34/55/205; this residual survives it. Red test
  `known-red-insert-above-window-at-max-drifts-the-fraction`.
- ~~Frameless attached children lay at offset 0~~ **RESOLVED (stage 2c)**:
  `parked-frame` clamps the per-child frame to the source's own edge. Seed 340 is
  clean and the `to-src` set-point flip is now inert on the fuzz tree (identical
  10 seeds with and without it) — stage 3 may take it on its own merits.

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
- [x] 11. Capture-mode: design + implement (corrections statically off, rebase returned as value) — agent: general-purpose, model: opus
- [x] 12. Design anchor-primary migration (truth = anchor idx+intra offset; equivalence criteria; kill-list of subsystems) — agent: Plan, model: fable
- [x] 13. Implement anchor-primary in stages behind harness equivalence — agent: general-purpose, model: fable
- [x] 14. Cleanup: remove dead subsystems, update docs/CollectionRectAnimator.md + memory — agent: general-purpose, model: sonnet

## Outcome

**Shipped.** The engine's coordinate rebase, its pass-mode gating and its scroll
truth are each owned by one structure instead of being emergent properties of
several mechanisms agreeing:

- **One pass mode, one capability table** (`pass-mode` / `pass-caps` /
  `pass-allows?`): `:resting | :settled | :capture | :segment`, decided once per
  `performLayout`. No engine code gates on `tweenAnim`/`segTween` any more.
  `:settled` is the mode the engine used to re-derive at every gate and get wrong.
- **One rebase transport**: `passRebase` accumulator (one write path, `rebase+!`),
  one emission site (`emit-correction!`, capability-asserted), three arms — emit /
  value / absorb — and a debug exactly-once assert. A capture returns
  `{:rebase :extent :anchor-frame}` as plain values; its emit arm is structurally
  unreachable. The absorption arm is total, so a produced rebase can never vanish.
- **One leading-side producer**: the Δ-epilogue. Measurements supersede, the
  window translates rigidly once, and the ANCHOR's own displacement is what gets
  emitted, under the latching invariant (aggregate drift never moves a latched
  anchor).
- **One scroll truth**: `viewAnchor {:idx :key :off :extent :frac}`, sampled from
  displayed frames every pass, re-resolved by KEY across a count change. The
  segment set-point preserves its consumed FRACTION, not its top edge.
- **One boundary rule**: an answer a source does not have resolves to the source's
  own edge, never to 0 — window queries, per-child frames (`parked-frame`), and
  the correction floor at the content start.

**Numbers**

| | before | after |
|---|---|---|
| default suite (`-x "known-red || fuzz"`) | 183 | **326** green |
| viewport-level tests | 0 | harness + 5 fuzz profiles, 230 episodes |
| oracles | — | **O1–O12** |
| fuzz failing seeds | 62/120 episodes (campaign 1) | **2** |
| `-t known-red` | — | **2** (both documented below) |
| render-object fields | 44 | **42** |

**Deleted subsystems**: `restingTop` + `capture-resting-top!`, `crossLayoutReanchor`,
`reanchorShift`, `landingEmitted`, `pendingLandingWs`, `pendingRebase` (→ a bare
`reseedCause` keyword), `landing-reseed-decision`, the landing cache-drop/`:init`
dance, `frontier-after-prune`, `origin-refine-emit?` (→ one `:lead-emit?`
capability), the emit-vs-frozen branch pair, the `wsx` shifted-frame loop,
`shift-attached-from!`-as-producer, `fromRects`/`fromExtent`.

**Docs**: `docs/CollectionRectAnimator.md` §9–§9f (the scroll model) and
`docs/CollectionTesting.md` (harness, oracles, fuzz, red-test workflow).

**Carried deferrals** — both red-tested, both `:o6`, both about the anchor's
IDENTITY rather than its target frame:

- **NEW-16** `known-red-morph-loses-the-anchor-above-the-capture-window`
  (seed 235). A morph into a denser layout: the capture materializes the window at
  the CURRENT offset, so after the set-point moves the viewport there is no child
  at the new top and the anchor re-samples onto the window's first child. Closing
  it needs the capture to materialize the END window, which is knowable only after
  the set-point the capture itself produces — the §7b/§7c chicken-and-egg. Pick it
  up at `segment-start!`'s capture-window choice, not at the set-point.
- **NEW-17** `known-red-insert-above-window-at-max-drifts-the-fraction`
  (seed 337, untraced). An insert above the window with `pixels` resting ON
  `maxScrollExtent` drifts the masonry anchor's consumed fraction by ~6px. Same
  shape as NEW-14 (count-change re-anchor vs leading re-measure), which step 13b
  closed for seeds 34/55/205. Pick it up by tracing the epilogue per pass at the
  bottom boundary, where the RANGE clamp is designed to win over the anchor.

## Step results

### 14. Documentation, verification, close-out
- Status: done
- Files changed: docs/CollectionRectAnimator.md (header, §4 note, §5 state, §9
  rewritten into §9–§9f, §14), docs/CollectionTesting.md (new),
  layout.cljd (contract docstring: `:collapse-axis` also picks the parked-cell
  axis; `:approx-offset` with no aggregate basis is treated as absent), this plan.
  Zero engine behavior changes.
- Verification, all on this commit's tree:
  - `-x "known-red || fuzz"` → **326 green**, "All tests passed!"
  - `-t known-red` → **0 passed, 2 failed** — NEW-16 and NEW-17, nothing else.
  - `-t fuzz --timeout 20x` (all five profiles, 230 episodes) → 21 batches green,
    **2 failing seeds**: `:no-jump` 235 op 12 `:o6`, `:churn` 337 op 8 `:o6`.
    Identical to step 13c's recorded set.
  - `bin/check` → no Dart errors.
  - `git grep` over `src/`: no TEMP / debug instrumentation anywhere; fb24f35's
    prints were reverted in f609581 and did not creep back.
- Half-done review of the plan's ~55 commits: none. Every stage that claimed a
  fix has its commit, its suite number and its fuzz delta recorded; the two
  survivors are tagged reds, not silent skips.
- The layout contract needed no new hooks: `:run-start` was documented when it
  landed (9c), `:renewal-index` / `:approx-offset` / `:collapse-axis` kept their
  meaning. Only the two facts above were unstated.

### 11. Capture-mode boundary (pass mode + capability table)
- Status: done (commit 97b1815 — one commit; the mode, the capability table and
  the emission assert are one indivisible consolidation, and every intermediate
  split would have left the engine deriving segment-ness two ways at once).
- Files changed: render.cljd (`pass-mode`, `pass-caps`, `pass-allows?`,
  `emit-correction!`, `passMode` field, 17 call sites; `origin-refine-emit?`
  deleted, `refine-emit?`/`commit-prune?` take the capability),
  render_test.cljd (+2 pure deftests, polarity flips), harness.cljd
  (`:pass-mode` in the engine map), this plan.
- Suite: default `-x "known-red || fuzz"` 318 → **320** green (+2 pure tests);
  O7 in-suite; `-t known-red` **still empty**; `bin/check` clean.
- **The mode**: `:resting | :settled | :capture | :segment`, decided ONCE per
  `performLayout` from (tween-anim?, curGen vs segGen, segTween?, domain-exited?)
  — the dispatch cond, lifted into a pure kernel — and stored on the render
  object. `pass-caps` gives each mode 12 named capabilities (`:emit-rebase?`
  `:refine?` `:origin-refine?` `:overscan?` `:prune-commit?` `:measure-feed?`
  `:tripwire?` `:canonical-assert?` `:epilogue-asserts?` `:union-window?`
  `:set-point?` `:seg-tail?`); `pass-allows?` debug-asserts that a capability
  HAS a row, so a new feature cannot inherit an accidental OFF.
- **17 sites consolidated**: performLayout's overscan gate / dispatch / epilogue
  asserts / prune gate; flow-layout!'s seg-tail, refine + origin-refine, union
  window, EMA, tripwire, value-vs-emit arm; finish-flow-geometry!'s canonical
  assert; seg-scroll-correction; indexed-layout!'s seg-tail, union window, EMA,
  correction exit; segment-start!'s two capture calls. No engine code reads
  `tweenAnim`/`segTween` as a GATE any more — the six remaining reads are state
  ownership (listener swap, tween value, settle/clear).
- **Boundary asserts**: `emit-correction!` is the single correction-writing path
  and debug-asserts `:emit-rebase?`, so the emit arm is structurally unreachable
  in `:capture` — the complement of 9a's absorption asserts, which fire from the
  other side (a produced rebase with no channel). `flow-layout!` loses its `mode`
  argument (the pass carries it) and keeps `{:rebase :extent}` as the design's
  `:value` arm, documented at the definition as step 13's input.
- Latent divergences found (2), both logged in the Decisions log: (1) the
  indexed capture's EMA feed depended on whether a previous segment tween was
  still around — FIXED (capture never feeds, matching the flow capture);
  (2) `:settled` defers refinement / origin corrections / overscan / commit
  prune although its emission arm is live and no set-point competes — PRESERVED
  as-is (a behavior change, not a consolidation) and handed to step 12/13.
- **Fuzz: byte-identical to the pre-change commit over all four profiles**, 220
  episodes (`:full` a-f, `:churn` a-d, `:no-layout` a-f, `:no-jump` a-f): the
  same 20 failing seeds, at the same op indices, with the same oracles
  (`:full` 2/17/23/34/38/52/55, `:churn` 322/337, `:no-layout`
  111/122/133/135/137/151, `:no-jump` 205/220/222/235/251). Campaign-3's
  residual set, unchanged — the consolidation neither fixed nor broke anything.
- Next-step impact: step 12's anchor-primary design gets the mode axis as its
  staging surface — an anchor-primary driver is a new capability column, not a
  new set of `tweenAnim` checks — and `{:rebase :extent}` as the already-pure
  capture interface. The `:settled` divergence is the first thing to decide
  there, since anchor-primary makes "the offset model is frozen" mean something
  different.

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

### 13a. Anchor-primary stages 0-2
- Status: stages 0-2 done; stage 2's unmet gate closed by stage 2b (below)
- Stage 0 (20867d1): oracles O10 (anchor rest stability: (key,frac) bit-stable
  across idle pumps), O11 (truth equation off + frac*extent == scrollOffset),
  O12 (extent quiescence); O9 split into an unconditional per-pass WORK probe
  (materialization-budget arithmetic, x2 mid-segment) + relaxed mid-segment
  retention denominators — this closed step-10's deliberately-open :o9 family
  (7 seeds) as calibrated green. BouncingScrollPhysics rig option + :bounce
  fuzz batch (seeds 401-410) green. All three new oracles hold on the current
  engine — no known-red additions. Suite 320 green.
- Stage 1 (50b169e): seed-cache! re-expressed as resolve/reassign — the anchor
  slot resolves FIRST (resolve-anchor-slot: by key on count change, by index
  cross-layout), classification is one pure kernel (seed-plan: :keep/:far-top/
  :far-inverse/:landing-init/:cold-deep/:at-0/:anchor/:restart-0), each old cond
  arm becomes an effect of that plan. align-start!/walk! mechanically unchanged
  (at rest the band IS the anchor's frame). Suite 321 green; fuzz byte-identical
  to the stage-0 baseline over all five profiles.
- Fuzz baseline carried into stage 2 (13 seeds): :o6 111/151/17/205/235/34/52/55,
  :o5 23/133/322/337, :o1 251; bounce batch green.
- Stage 2 (bdf7228): the Delta-epilogue. backfill-leading! measures and returns
  ONE leading-extent delta; reflow-from-checkpoint! keeps only the ex-frozen
  mechanics (prefix re-anchored to abut the retained head) and reports the shift
  it used; translate-window! moves the whole laid window rigidly once, at pass
  end, and the accumulator emits it. Dead: the emit-vs-frozen branch pair, the
  shifted-frame `wsx` loop, shift-attached-from!-as-producer, the landing
  cache-drop/`:init` dance, landing-reseed-decision, seed-plan's :landing-init,
  pendingRebase (-> reseedCause, a bare keyword). New pure kernel: lead-delta.
  Suite 321 green; -t known-red empty; bin/check clean.
- Stage 2 fuzz: **gate NOT met**. 13 -> 14 seeds. Green: :o6 111. NEW: :o6 123
  (:no-layout, insert above window), :o5 7 (:full, extent-below-content
  mid-segment). Surviving :o6 17/34/52/55/151/205/235; :o5 23/133/322/337;
  :o1 251. The set is not a strict subset, so Q5's fallback does not carry it.
- Stage 2 triage (traced, seeds 52 + 123): the o6 count-change family is NOT the
  leading-repair class the design assumed — it is the SEGMENT set-point.
  segment-start! reads the anchor's `to` slot from the CAPTURE-placed cache
  entry, and a count-change capture is seeded by the key re-anchor at the
  anchor's OLD offset, so `to-off == from-off`: point-desired is constant and
  the set-point emits nothing all segment while the tween moves the content to
  its live-re-flow target. Seed 52 (masonry, `[:remove 90]`): segAnchor
  {:from-off 4044.5 :to-off 4044.5}, zero corrections over the segment, anchor
  ends 280px high. Identical shape at seed 123. Stage 2 changes nothing here,
  which is why the family neither closes nor shrinks; the membership churn (111
  out, 123 in) is a chaotic re-roll of the same latent defect.

### 13a stage 2b. The segment `to` frame (o6 count-change family)
- Status: done (one commit). Fuzz **14 -> 11 seeds, a strict subset**; suite 321
  green; `-t known-red` empty; `bin/check` clean.
- Traced (seeds 52 + 51, per-pass): stage 2's diagnosis named the right symptom
  and the wrong side. `to-off == from-off` is CORRECT for seed 52 — `[:remove 90]`
  is far above the window, so nothing in it may move. The drift came from the
  OTHER endpoint: `live-only-flow-window` re-flowed the whole window from its head
  with `flow-seed-state`, i.e. `:anchor` — which for masonry LEVELS the columns and
  for wrap fabricates a run break. The capture walk had placed that window
  anchor-seeded + backfilled, a placement no single forward re-flow reproduces, so
  the tween's `to` frames sat 126px (seed 52) / 165px (seed 51) off the placement
  the pass had just made, and the segment slid the content there while the
  set-point — correctly — emitted nothing.
- Fix: the capture walk's own placement IS the target frame. Survivors above the
  first in-window dying cell keep their cache entries verbatim; only the suffix
  below it re-flows, seeded with `state-before` at that cell — the gap closes,
  everything else stays exactly where the capture put it. O(window), no new state.
- Seeds: `:o6` 17 / 123 / 151 green; 52's `[:remove 90]` green (it now fails at
  op 16, `[:layout :wrap]` — a morph, not a count change). `:o6` 34 / 55 / 205
  survive: traced (55) to a DIFFERENT mechanism — the far-jump estimate refinement
  emits a +2991px correction that lands the offset past `maxScrollExtent`, and the
  boundary clamp eats 63px. `:o6` 235 + 52@16 are layout-morph key-moved.
- `:o5` 7 (`:full` op 19, extent-below-content mid-segment, introduced by stage 2)
  SURVIVES and is still untraced — out of stage-2b budget, and unrelated to the
  segment frames by its oracle. Carried as a fuzz seed, no red test: reproducing
  it needs the trace stage 2 never took either.
- MEASURED AND REJECTED: reading `to` from `to-src` (the task's proposed fix, and
  the flip of 9a's decision). With the frames fixed the two sources agree
  everywhere except below an in-window dying cell, so it closes NOTHING extra on
  today's seeds, and it adds seed 340 (`:churn`, `:o9` work: 163 materializations
  vs a 158 budget). That seed exposes a pre-existing latent, not the flip itself: a
  child attached ABOVE the frozen `to-src` window has neither a committed rect nor
  a `to` frame, so `keyed-tween-layout` lays it at offset **0** (`zero-frame`); the
  next capture then finds firstChild at 0, inverse-seeds by estimate, and walks 163
  children in one pass. Fix that first, then the flip is free.

### 13a stage 2c. Frameless-child-at-0 + the `:o6` boundary class retriaged
- Status: done (2 commits — one per defect). Fuzz **11 -> 10 seeds**, a strict
  subset (`:o5` 322 green); seed 340 clean. Default suite 321 -> **323** green
  (+2 pure tween tests), O7 in-suite; `-t known-red` **1** — this step's
  documented deferral; `bin/check` clean.
- **Blocker FIXED (1d5fd70)**: `keyed-tween-layout` gave a cell with no committed
  history AND no target frame `zero-frame` — offset 0. `parked-frame` resolves the
  source's nearest EDGE frame instead and parks the cell there, collapsed on the
  layout's own axis, carrying that frame's extent as the stable full size the
  child is laid out at. Substituted into `to` before the enter/exit/leave cond, so
  every branch inherits it. **The principle**: rebase stage 5's window clamp — an
  answer the source does not have resolves to the source's own edge, never to 0 —
  now covers all three source queries, not two.
- Traced (seed 340, `:churn`, per pass): capture at cf=145 ends with the window
  [145..157] and child **144** attached with no cache entry; from-relay! lays it
  at 0.0, it becomes firstChild@0, and the next capture inverse-seeds from there
  and walks **163** children (budget ~78). After the fix that pass materializes
  **14** and no `zero-frame` call is reached at all in the episode.
- The flip was applied locally to reproduce (stage 2b's rejected `to-src`
  set-point read) and reverted, per the task. With the fix in place the campaign
  is the SAME 10 seeds with and without it — the flip is now inert, so stage 3 can
  decide it on its own merits rather than on seed 340.
- **`:o6` 34/55/205 RETRIAGED as NEW-14, deferred (b6239e5)**: not the
  boundary-clamp class stage 2b named. 34 (offset 2745 of 7967) and 205 (11399 of
  26042) drift identically nowhere near an edge; all three fail on
  `[:insert <above the window>]`. The count-change re-anchor pins the anchor at
  its old offset (no rebase, by design), `backfill-leading!` measures the
  re-flowed leading extent and the Δ-epilogue translates the whole window by it —
  2991.8px at seed 55, against the anchor's own 2889.8px displacement over the
  same pass. Seed 55's overshoot additionally lands past `maxScrollExtent` and the
  spring takes another 39px, which is what made the clamp look causal.
  Clamping the emission is **not available to a sliver**: `SliverConstraints`
  carries `precedingScrollExtent` and nothing about the slivers that follow, so
  `maxScrollExtent` is not knowable and a self-computed clamp would fire early in
  a composed viewport — and it would leave the 102px untouched. Red-tested at
  seed 34 (`known-red-insert-above-window-drifts-by-the-leading-measure`),
  reproduced 3x, deferred to stage 3's anchor-seeded walk.
- Carried untouched, as instructed: `:o6` 235 + 52@16 (layout-morph key-moved),
  `:o5` 7/23/133/337, `:o1` 251.
- Next-step impact: **stage 3 is unblocked** — the anchor can be widened without
  the segment stranding children at the origin behind it. Stage 3 inherits NEW-14
  as its first correctness target: it is the leading-estimate-vs-anchor
  disagreement the anchor-seeded walk exists to remove, and its red test is the
  acceptance criterion.

### 13b stage 3. NEW-14, the latching invariant, and `:lead-emit?`
- Status: done (2 commits — 4f9f13c the anchor/latching fix, 54edbf4 the
  capability collapse). Default suite 323 -> **326** green, O7 in-suite;
  `-t known-red` **empty** (NEW-14 un-tagged); `bin/check` clean.
- Fuzz **10 -> 6 seeds**. Green: `:o6` 34 / 52 / 55 / 205, `:o5` 7. Carried:
  `:o6` 235, `:o5` 23 / 133, `:o1` 251. Flipped: `:o5` 337 op 6 -> `:o6` 337
  op 8 (same seed, later op). **NEW: `:o1` 56** — the one gate deviation, logged
  as NEW-15 in Open questions; it is a jump-past-the-end `:restart-0` walk the
  work probe catches, on a path stage 3 does not touch. The gate asked for a
  strict subset; the set is a subset plus that one seed, and no oracle or probe
  was weakened to get there.
- **NEW-14, traced (seed 34, per pass)**: the count-change re-anchor pins the
  anchor at its old offset (by design), the backfill measures the leading extent
  against a checkpoint, and the Δ-epilogue translates the window by it — that
  part is CORRECT and holds the anchor. The drift arrives on the NEXT pass: the
  translation empties the band, so the pass re-derives the window from the
  frontier, the run-chain re-packs, index 201 moves 115px against the anchor, and
  `anchor-before` — which could only read the cache the translation had just
  dropped — returned nil. Nothing held the anchor across the re-derivation.
- Fix, three parts, all in the epilogue:
  1. `anchor-before` falls back to the engine's own `viewAnchor` when the seed
     emptied the band, re-resolved by KEY across a count change and lifted into
     the frame the seed already rebased to. Gated on `anchor-preserved?` — the
     new pure kernel over `seed-plan` — so a far jump / cold start still
     teleports the anchor instead of fighting the teleport.
  2. `anchor-delta` reads the LAID CHILD, not the cache: the translation drops it.
  3. The epilogue is ONE producer — translate first, then emit the anchor's own
     displacement over the whole pass. A pass with no anchor to hold falls back
     to the raw measure. `seed-cache!`/`normalize-cache!` now return their plan.
- **The latching invariant, arrived at empirically**: part 1 alone closed NEW-14
  but turned seed 52 into an `:o1` (10 layout cycles exhausted) and added seed 33
  — the anchor-follow was emitting once per round of an estimate cascade that
  already emitted once per round. The cascade's engine is `synth-renewal-
  checkpoint`: its `:offset` is `:approx-offset`, a function of the measured EMA
  that the window itself feeds, so honouring its seam makes the frame chase its
  own aggregate at 0.75/pass. A synthesized checkpoint's seam now reports 0.0
  under a latched anchor (its re-flow still runs — real runs above the window is
  what it is for) and is honoured only when the pass reassigned the anchor.
- MEASURED AND REJECTED, both logged in the Decisions log: re-deriving the window
  inside `translate-window!` (stage 2b's levelled-runs defect again — seeds 34
  and 52 both drift worse), and seeding every `:anchor` re-seed at viewAnchor's
  slot per design §2.2 step 2 (breaks NEW-2's wrap tiling: the anchor is mid-run
  and `:anchor` opens a fresh run there).
- **Capability collapse + `:settled` ON**: `:refine?` + `:origin-refine?` ->
  `:lead-emit?`; the ballistic sub-gate stays in the epilogue (`refine-emit?`
  filters the seam, the origin residuals ride the capability unfiltered, #B).
  `:settled` gets `:lead-emit?` true — step 11's deferred divergence, dissolved
  as design §3 owns it. `:overscan?`/`:prune-commit?` stay OFF there. Pinned by
  `settled-lands-exactness-during-the-playing-clock` (kernel-level: flipping the
  row turns it red) plus a harness scenario for the mid-clock jump-to-top. The
  change is fuzz-neutral — identical 6 seeds before and after.
- **Field delta: 0** (`anchoredTo0` kept, no `leadEstimate` added — Q1 resolved
  against its leaning, see the Decisions log). Deleted: `frontier-after-prune`
  (the identity kernel) and its write site, one pass-caps column. Added:
  `anchor-preserved?` (pure, no state).
- Not done, out of stage-3 scope as scoped here: stages 4 (the `:o5` remove-extent
  family) and 5 (cleanup handoff).

### 13c stage 4-5. The residual fuzz set + cleanup
- Status: done (6 commits). Fuzz **6 -> 2 seeds**, a strict subset. Default suite
  **326** green, O7 in-suite; `-t known-red` **2** — this step's two documented
  deferrals; `bin/check` clean.
- Green: `:o1` 56 (NEW-15), `:o1` 251, `:o5` 23, `:o5` 133. Carried: `:o6` 235
  (NEW-16), `:o6` 337 (NEW-17). No new seeds on any of the five profiles.
- **`:o1` 56 / NEW-15 (9a347d0)**: NOT `align-start-fallback` — see the Open
  questions entry. `inverse-seed!`'s guard rejected the saturated last-index seed
  and kept the far band. The last index is now exempt.
- **`:o1` 251 (b72ea77)**: `indexed-layout!` took the leading edge from
  `first-index`'s frame while the leading walk had stopped at index 85 of a
  non-monotone masonry morph — `calculateCacheOffset` asserted `from <= to`
  (leading-off 4845 vs trail-end 3844). It is now the minimum laid offset.
- **`:o5` 23 (b72ea77)**: the segment advertised the lerping total (9345) while
  laying edge-sliding enter cells out to 10472 — 41 cells the capture
  pre-materialized below the window, entering from the cache-window edge. Both
  segment writers now advertise `max(total, laid trail)`.
- **`:o5` 133 (d76573b)**: a landing at the scroll top emitted the anchor's full
  -267 displacement; `pixels` parked at -267 for all 150 settle frames with no
  activity to spring it back (the sliver sees scrollOffset 0, so no further pass
  runs). `correction-floor` + `clamp-rebase!`.
- **`:o6` 235 (66db72a, PARTIAL)**: traced per-pass. Two distinct defects in one
  seed. (1) `segment-start!` read the anchor's `to` slot from the frozen snapshot,
  which starts at the capture window; a morph into a denser layout packs the
  anchor above it, the lookup answered nil and the set-point went silent — the
  viewport-top key moved 22 items. FIXED: the capture returns the frame its own
  walk placed (`anchor-slot-frame`, read before `trim-front!`, riding the
  epilogue's rigid translation). The set-point now converges exactly on the
  anchor's target (traced: so 489.59 -> 102.84). (2) NEW-16, deferred: the
  segment's window never covers that new top, so the anchor re-samples onto the
  window's first child and the next resting pass holds it (+115px).
- **`:o6` 337 (NEW-17)**: not traced — out of budget after 235. Recorded with its
  measured signature and red-tested.
- **Cleanup (955b927)**: `fromRects` + `fromExtent` are written at the end of
  `segment-start!` and never read — the `from` rects reach the tween as a local
  and the frozen extent as its `:max-extent` lerp input. **Field delta -2**
  (44 -> 42). **Kernel delta +2**: `correction-floor` (pure) and
  `anchor-slot-frame`; nothing else in render.cljd is unreferenced — every
  `defn` has a live call site, `pass-caps` has no unused column, and
  `landing-decision`/`underflow-decision` are still `backfill-leading!`'s two
  origin measurements. `anchoredTo0` stays (stage 3's Q1 resolution).
- Not done: step 14's doc pass (docs/CollectionRectAnimator.md §9/§2.4, memory).
