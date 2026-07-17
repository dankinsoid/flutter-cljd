# WS2 — Consolidated design: approximate-layout protocol & phased plan

Final design gate before WS2 implementation. Consolidates the three settled pieces
(2.1 overscan — `ws2-1-overscan-design.md`; 2.3 backward re-flow + refinement —
`ws2-3-offset-refinement-research.md`; the protocol-hook decision recorded in
`collection-scroll-anim-followups.md` §2.2) and adds the missing piece: the exact
shape of the approximate-layout protocol (2.2) and the unified architecture that
ties them together.

Nothing here re-litigates a recorded decision; where this doc goes *beyond* the
earlier ones (a memory bound on the checkpoint cache, cold-start index seeding by
inverting `:approx-offset`, the monotonicity constraint that inversion needs, the
temporal separation of segment vs. refinement corrections) it is flagged in §6.

File:line pointers re-located against the current tree (WS1.1/1.2/1.3 landed —
78e1690 spring retarget, fe412a1 host dedup, 862a5ca). `render.cljd` unless noted.

---

## 0. The defect, restated as one root cause

There is no true from-index-0 flow geometry until you physically scroll to the
start. Every window above the first-visited region renders at an **estimated** main
offset, and the leading margin degenerates because `backfill-leading!` (703-741) is
list-shaped: one child per main slot, `:cross` hardcoded `0.0` (728), never runs the
layout's `:place`/state. Two symptoms, both flow-only (indexed targets are exact):

- **#5** — a wrap window above the anchor collapses to a single left column; correct
  wrap only snaps in when the backfill finally reaches index 0 and re-flows from the
  top (`finish-flow-geometry!` 765).
- **#6** — `segAnchor :to` (segment-start! ~1257) is the survivor's offset in a
  `to-src` anchored at the survivor's **current (approximate)** offset. The extent
  shrinkage lives *above* that window and is invisible to the capture, so `shift =
  to − from ≈ 0`; scrollOffset barely moves during the segment and the true shorter
  extent only appears at the post-settle top re-flow → `ScrollPosition` ballistic-
  clamps the residual (the "offset settles after the animation" the user sees).

Fix shape (settled): give every flow layout an O(1) offset estimate and exact
backward re-flow via renewal-point checkpoints; reconcile absolute offset by
`scrollOffsetCorrection` (pixel-stable), scoped to a velocity-widened overscan
window.

---

## 1. Protocol spec (WS2.2)

### 1.1 Parts & boundaries

The flow layout stays **pure geometry** (`layout.cljd`); the engine (`render.cljd`)
owns materialisation, caching, correction, and now the checkpoint store and the
measured aggregates. WS2.2 adds exactly **two optional layout hooks** plus one
**engine-owned option** (`approximate-item-size`) that feeds a pair of measured
aggregates into `env`. No new engine→layout coupling beyond that; overscan and
refinement *consume* these hooks and add no further layout surface.

Two aggregates the engine maintains and puts on `env` (built in `make-env`,
`layout.cljd:72`, refreshed each pass from the live window):

- `:avg-main`  — EMA of measured child main-extent over the live window.
- `:avg-cross` — EMA of measured child cross-extent over the live window.

Both are **seeded** by `approximate-item-size` before any cell is measured, and
**measured reality overrides** the seed as soon as the live window has cells (the
seed is a cold prior, nothing more). These are the only inputs an `:approx-offset`
implementation needs beyond the layout's own baked opts (`:spacing`, `:run-spacing`,
`:cross-axis-count`, …).

### 1.2 The two hooks — signatures

```clojure
;; OPTIONAL. Called by walk! at each placed index. Returns true iff the flow state
;; BEFORE i (== state after i-1) is EXACTLY reconstructible from (i, offset_i) via
;; :anchor — i.e. i begins a fresh, independent run. The engine retains such indices
;; as checkpoints (see §2.1). `frame` is i's just-placed geometry {:offset :cross}.
:renewal-point? (fn [env ^int i frame] -> bool)

;; OPTIONAL. O(1) main-axis offset estimate for an ARBITRARY index, with NO walk.
;; Must be MONOTONE NON-DECREASING in i (the engine binary-searches it to invert
;; offset->index for cold-start / far-jump seeding — see §2.2, §4.1). Reads the
;; engine's measured aggregates off env (:avg-main/:avg-cross) plus baked opts.
:approx-offset  (fn [env ^int i] -> double)
```

**nil semantics (fallback behaviour):**

- `:renewal-point?` absent → the engine has no cheap exact boundaries, so it retains
  a **periodic full-`:state` checkpoint** at an index stride (the masonry path).
  Retention is still exact (it stores the real `:state`), just not aligned to
  natural boundaries. A layout that supplies the hook gets minimal, well-distributed
  checkpoints; one that doesn't still gets correct backward re-flow, at a coarser
  stride.
- `:approx-offset` absent → the engine falls back to a **layout-agnostic linear
  density**: `off_ref + (i − i_ref)·d`, where `(i_ref, off_ref)` is the nearest
  known frontier and `d` is the per-index main advance measured over the live window
  (or, before any window, `:avg-main`/`:estimate`). Good enough for the indicator;
  worse than a structural estimate for multi-column layouts, hence the hook.

Neither hook changes an existing layout: both are optional with derivable defaults.
The checkpoint `:state` payload is always the exact flow state `walk!` already
computes (658-660) — storing it uniformly is the simple design; the wrap
"state is `:anchor`-derivable from offset alone" observation is an optional payload
shrink we need not take.

### 1.3 `:anchor` — the two-mode contract (doc comment as it will become)

```
:anchor  (fn [env index offset] state) — reconstruct the flow state as if children
         0..index-1 ended at `offset`. TWO MODES, by whether `index` is a renewal
         point (see :renewal-point?):

           - EXACT at a renewal point. `index` begins an independent run, so the
             true state-before-`index` is a function of `offset` alone (wrap: an
             empty run at y=offset; list: the offset itself). A forward walk seeded
             here reproduces true geometry — this is what makes backward re-flow
             exact (engine seeds :anchor at a cached run boundary, walks forward).

           - DEAD-RECKONING elsewhere. At a non-renewal index the real state depends
             on items < index that :anchor cannot see, so it fabricates a PLAUSIBLE
             state (wrap: a spurious run break; masonry: level columns). Good enough
             to resume a walk after cache loss; NOT exact. Anchoring mid-run at
             one-child-per-slot is precisely the old single-left-column bug.

         Optional but strongly recommended; the engine dead-reckons from the first
         attached child after a cache loss instead of measuring from index 0.
```

### 1.4 Per-layout implementations (sketch)

| layout | `:renewal-point?` | `:approx-offset` | `:anchor` mode | backward re-flow |
|---|---|---|---|---|
| **list** (variable) | `(fn [_ _ _] true)` — every index | `(* i (+ avg-main spacing))` | exact everywhere (state = scalar offset) | trivially exact per item |
| **list** (fixed `:item-extent`) | n/a — routes to the **indexed** engine (`fixed-list-layout`, layout.cljd:115), O(1) exact | via `:geometry` | n/a | none needed |
| **wrap** | `(fn [_ _ frame] (= 0.0 (:cross frame)))` — run starts | `run-idx·(+ avg-main run-spacing)`, `run-idx = ⌊i / items-per-run⌋`, `items-per-run = max(1, ⌊cross-extent/(avg-cross+spacing)⌋)` | exact at run starts, dead-reckon mid-run | exact once run boundaries cached |
| **masonry** | absent → engine periodic full-`:state` stride | `⌊i/n⌋·(+ avg-main main-spacing)` | approximate everywhere (assumes level columns) | exact from a retained `:state` checkpoint |
| **grid** (`:grid`, delegate) | **indexed — no flow protocol.** `:geometry` is O(1) exact everywhere | subsumed by `:geometry` (exact) | n/a | none — why indexed targets never showed #5/#6 |

**Grid, explicitly:** nothing changes. It has `:geometry`/`:first-index`/
`:last-index`/`:max-extent`, so every index is exact math. The engine's approximate
path (cold-start seed, flying-regime placement) just calls `:geometry` for indexed
layouts instead of `:approx-offset`. Grid gains no hooks.

**`approximate-item-size` flow:** an **engine/collection option** (threaded like
`:spacing` through `sliver-collection*` args, not a per-layout knob), accepted as a
main-extent double or a Size. It seeds `:avg-main`/`:avg-cross` before measurement;
once the live window has real cells, the measured EMA overrides it. It reaches
`:approx-offset` only via those env aggregates, so it improves cold-window and
flying-regime placement uniformly (default-density path included) without any
layout needing to know about it.

---

## 2. Unified architecture

### 2.1 Run-boundary (checkpoint) cache

- **Owner:** `CollectionRenderSliver`, new mutable field `checkpoints` — a sparse,
  index-sorted map `{^int index → {:offset ^double :state <flow-state>}}`. `:state`
  is the flow state *before* that index (the seed for a forward re-walk from it).
- **Populated:** in `walk!` (623-677), when a child is (re)placed. Record index `i`
  when `(:renewal-point? layout)` returns true, or — when the hook is absent — when
  `i` hits the periodic stride. The `:state` to store is `state-before-i` (= the
  `state` local passed into `place-fn`, i.e. `:state` of `cache[rel-1]` or
  `baseState`). Source data already exists; this is the sparse GC-surviving subset
  of walk!'s per-index `:state`.
- **Data shape rationale:** `{:offset :state}` is all a backward mini-walk needs —
  seed `:anchor`/`:state` at the checkpoint, re-measure real children forward. Wrap's
  `:state` is `{:x nil :y offset :h 0}` (tiny); masonry's is the column vector;
  list's is a scalar.
- **Invalidation triggers:**
  1. **Item resize / measured-size mismatch at index j** — `walk!` already truncates
     the live cache at the mismatch `rel` (649); extend that to drop **checkpoints
     with index ≥ j** (their `:offset` prefix-sums shifted).
  2. **Cross-extent change** — the `crossCached` guard (158-161) already drops the
     live cache; drop **all** wrap/masonry checkpoints too (packing depends on
     cross-extent). List checkpoints are cross-independent but dropping all is
     simplest and cheap.
  3. **Item-count / layout `:id` change** — `update-render!` (348-350) drops the
     cache on layout change; drop checkpoints in lockstep. On a count change with a
     diff available (the host supplies inserts/removes), drop **checkpoints ≥ the
     lowest touched index**; without index info, drop all.
- **Memory bound (this doc's addition — earlier docs said only "sparse subset past
  GC"):** checkpoints are pruned on the **same geometric rule as the committed cache**
  — a checkpoint whose index left the overscan window is dropped, **plus** the single
  nearest checkpoint strictly above the leading margin is retained as its re-flow
  seed. Never pruned by count or age. Bound: `O(overscan_extent / min_run_height)` —
  a band around the viewport, not the full `[0, current]` prefix. This is sound
  because near-0 exact landing only runs when the first materialised index is already
  small (a short walk to 0, whose ultimate seed is `:init` at index 0, needs no
  stored checkpoint); when scrolled deep we never try to reach 0 exactly, we use
  `:approx-offset` for the indicator.

### 2.2 `backfill-leading!` rebuilt on the checkpoint cache

Replace the list-shaped estimator (703-741) with checkpoint-seeded forward
mini-walks (reusing `tw/flow-frames` / the `live-only-flow-window` pattern):

```
;; materialise the leading overscan margin [ws_o, cache[0].offset) as REAL runs
while (pos? cacheFirst) and cache[0].offset > ws_o:
  target = cacheFirst - 1
  (k, seed-state, off_k) = nearest-checkpoint(<= target)   ; else (0, init, 0.0)
  ; 1. insert+measure the run's children as leading children (backward, cf-1..k)
  ; 2. place FORWARD k..cf-1 via flow-frames seeded {(:anchor|:state) at k, off_k}
  ;    -> real :cross, real per-run :offset (kills #5's single-left-column)
  ; 3. write offsets/cross; prepend entries; cacheFirst := k; baseState := seed-state
```

- **Cold-start / far-jump seed (this doc's addition):** when the cache is empty and
  `scrollOffset` is large with no anchor child, the engine **inverts** the monotone
  `:approx-offset` by binary search to pick the start index, attaches it at
  `:approx-offset(i)`, walks the visible window, and marks `anchoredTo0 = false`.
  This is what stops a cold deep start from walking O(n) from index 0. Indexed
  layouts use `:first-index` for this (already exact).
- **Absolute offset stays reconciled, not renewed:** `off_k` from a checkpoint may be
  stale (mutation above); the re-flowed run lands at the stale absolute offset and
  the accumulated delta is emitted as `scrollOffsetCorrection` (existing backfill
  correction path 832-834 + `anchor-delta` 835-836). Structure is exact; position is
  reconciled — the split the WS2.3 research proved.

### 2.3 Refinement loop × overscan regimes × committed-cache scoping — the interlock

Per-pass flow in **resting / landing** regime (flow driver, `flow-layout!` 798-838):

1. **velocity sensor** updates `(v, dir)` once per frame, correction-compensated
   (WS2.1 §2.1).
2. **overscan resolver** → window `[ws_o, we_o] = [ws − behind, we + ahead]` (WS2.1
   §2.2). Used for materialisation/GC only; geometry advertised to the viewport
   stays on the true `[ws, we]`.
3. **fidelity gate** → `:resting`/`:flying`/`:landing` (WS2.1 §2.3).
4. `pre-gc!`/`walk!`/`widen-window!`/`post-gc!` run over `[ws_o, we_o]`. `walk!`
   records checkpoints.
5. **refinement** = `backfill-leading!` (§2.2) over the leading overscan margin;
   each real run replaces its estimated extent; the accumulated delta → one
   `scrollOffsetCorrection`. **Near the start** (first materialised index small):
   keep re-flowing up to index 0, land 0 at offset 0, residual = final correction,
   set `anchoredTo0 = true`. **Suppressed** while `dir == up` and `|v|` above the
   resting threshold (WS2.3 tier 3 / TanStack's rule) — defer to the next low-
   velocity pass so a same-frame offset jump doesn't race an upward drag/fling.
6. `commit-pass!` (232-258) upserts committed and **prunes geometrically** to the
   overscan window (WS2.1 §2.4), exempting active-segment participants
   (`entering ∪ exiting ∪ keys of leaving`) while `tweenAnim` is non-nil. The
   `committed-max-age = 600` age prune (204-207) is removed. Checkpoints prune on the
   same rule (§2.1).

In **flying** regime: overscan disabled; place the visible window via `:approx-offset`
(indexed: `:geometry`); no backfill, no checkpoint retention beyond stride;
refinement corrections suppressed. Geometry still on the true window. As velocity
decays → landing regime re-enables real materialisation + refinement so it lands
exactly (WS2.1 §4).

### 2.4 Estimate-error accounting — the known-exact frontier

- Field `anchoredTo0 ^bool` on the sliver: **true** once a pass placed index 0 at
  offset 0 with an unbroken cache down to the live window; **false** on cold deep
  start, far jump, or any invalidation above the frontier.
- While `anchoredTo0`: reported `scrollExtent`/offsets are **truth** — the scroll
  indicator is exact.
- While `!anchoredTo0`: offsets are checkpoint-relative and `scrollExtent` comes from
  `:estimate` / `:approx-offset` extrapolation — honest-but-approximate (the
  indicator may drift slightly toward truth as runs materialise; monotone-ish beats
  frozen, matching RV's `computeVerticalScrollOffset`).
- Refinement's contract is exactly "drive `anchoredTo0` toward true": every dirty
  resting pass lowers `cacheFirst` via real re-flow; reaching 0 flips the flag and
  emits the residual. The frontier is `cacheFirst` + `anchoredTo0`; the "error" is
  the gap between the checkpoint-relative offset and truth, bounded by the
  correction and reconciled incrementally.

### 2.5 Overscan open questions — resolved or carried

1. **Advertised `cacheExtent` vs private overscan children** — **carry to device
   check (overscan phase).** Start by keeping `SliverGeometry.cacheExtent` on the
   true window and holding overscan children privately (engine owns
   `collectGarbage`). Fallback if keep-alive/semantics misbehave: advertise
   cacheExtent over the overscan region (framework-consistent, costs sibling cache
   budget).
2. **Velocity from `scrollOffset` delta vs host-pushed** — **decided: render-object-
   local, correction-compensated**, device-verify the compensation doesn't leak
   WS1.3/WS2.3 corrections into velocity. Fallback: host `NotificationListener`.
3. **Thresholds (`v_slow`, `v_fast`, `lookahead-ms`, `ahead-max`, EMA α, hysteresis)**
   — **carry to the calibration/verify pass**; ship as named tunable consts.
4. **Landing detection precision** — **resolved.** Landing engages on **deceleration
   alone**, independent of *where* it lands: real materialisation + refinement always
   re-anchor on the actual landed geometry. Exact-at-0 is not a separate path — it is
   the special case where the landed window contains index 0. No "near element 0"
   trigger needed.
5. **Reorder / box host** — **resolved: option accepted-and-ignored.** Overscan is a
   no-op for the box host (no scroll; `remainingCacheExtent = main-max`) and for the
   list-only `SliverReorderableList` reorder path. Device-check that a present-but-
   ignored `:overscan` doesn't break reorder.

---

## 3. Phased implementation plan

Each phase is independently landable and separately verifiable. `bin/check`
(`bin/check --clean` if a stale mixin surfaces) gates every phase; device checks use
the example app.

**Phase A — Protocol hooks + `:anchor` two-mode doc + measured aggregates.**
Add `:renewal-point?` / `:approx-offset` to list/wrap/masonry; rewrite the `:anchor`
docstring (§1.3); thread `approximate-item-size` through `sliver-collection*`; have
`make-env` carry `:avg-main`/`:avg-cross` (seeded by the option, EMA from the live
window). No engine behaviour change beyond the default-density fallback being
available.
*Verify:* unit tests per layout — `:renewal-point?` truth pattern (wrap only at
`cross==0`), `:approx-offset` monotone and within tolerance of a real walk, seed→
measured override. `bin/check`. No visible device change.

**Phase B — Checkpoint store + real backward re-flow + cold-start seed. [FIXES #5]**
Add `checkpoints` field; record in `walk!`; rewrite `backfill-leading!` as
checkpoint-seeded forward mini-walks (reuse `flow-frames`); invert monotone
`:approx-offset` for cold-start/far-jump index seeding; checkpoint invalidation
(resize / cross / count-or-id). Scope to the plain cache window `[ws, we]` for now
(Phase E widens it to overscan).
*Verify:* unit — backfill produces correct multi-column `:cross` for a wrap window
above the anchor; checkpoint drop on resize above. Device — scroll a wrap collection
deep, scroll up: leading rows are correct wrap immediately (no single left column);
cold-start at a large initial offset attaches near the target without an O(n) hitch.

**Phase C — Progressive offset refinement + exact landing + suppression. [FIXES #6]**
Emit the accumulated re-flow delta as `scrollOffsetCorrection` every dirty resting/
landing pass; `anchoredTo0` frontier; exact landing at index 0; suppress above-anchor
corrections during backward drag/ballistic-up. Depends on B (real leading structure);
composes with the shipped WS1.3 relative segment correction (§4 scenario 6).
*Verify:* unit — telescoping deltas sum to the true residual; suppression gate honours
`dir==up`. Device — list→wrap morph while scrolled deep: scrollOffset moves *during*
the segment, no post-settle snap; `animateTo` a far target lands exactly on item 0;
no visible pixel jerk during any correction.

**Phase D — Committed/checkpoint geometric scoping (drop the age prune).**
Prune committed geometrically to the current window (segment-participant exemption);
remove `committed-max-age`; apply the same rule to checkpoints. Independent of the
velocity work — uses the plain window until E.
*Verify:* unit — `prune-committed` drops by window not age; participants exempt
mid-segment. Device — re-insert a long-gone key: enter-wipe (not a stale-`from`
glide); committed map stays bounded under a long feed.

**Phase E — Overscan: velocity sensor + regime gate + directional widening + flying
approximate path.** WS2.1 in full: `update-velocity!` + fields; `resolve-overscan`/
`overscan-window`/`fidelity-regime` (pure); widen `[ws,we]→[ws_o,we_o]` for
materialisation/GC in both drivers (geometry stays true); flying-regime approximate
placement via `:approx-offset`; re-scope B/C/D windows from cache to overscan; the
`:overscan` config surface. Depends on A (`:approx-offset`).
*Verify:* unit — `overscan-window` velocity law (rest→symmetric base, motion→capped
ahead), regime hysteresis. Device — fast fling doesn't degenerate; `animateTo` far
then immediate reverse drag stays smooth; threshold calibration (open Q3); resolve
open Q1/Q2/Q5 on device.

**WS3 attaches after E:**
- **WS3.1** (double before→after layout only when a segment will actually run) —
  guard the two-pass compute behind "a segment is pending". Independent of window
  mechanics; lands cleanly once the segment/compute structure is stable post-E.
- **WS3.2** (window-scoped diffing, gated on animations + key-fn) — diff only the
  current window rather than the whole item set. **Depends on the materialised-window
  boundary that E defines**, so it naturally follows E (and reuses the overscan window
  as the diff scope).

---

## 4. Scenario audit (end-to-end through the design)

**1. Cold start at large initial offset.** Cache empty, no children. Engine inverts
monotone `:approx-offset` (binary search) → start index near the offset, attaches at
`:approx-offset(i)`, walks the visible window, `anchoredTo0 = false`. Indicator uses
`:estimate`/`:approx-offset` (roughly right). No O(n) walk from 0. As the user scrolls
up, `walk!` records checkpoints and refinement lowers the frontier. **Holds** (Phase
B adds the inverse-seed; earlier docs assumed `:approx-offset` handled this but didn't
specify the inversion — see §6).

**2. list→wrap morph while scrolled deep (#5/#6 path).** Layout `:id` changes → cache
+ checkpoints dropped. The **visible** window flows correctly forward from the anchored
top (`walk!` + `:anchor(wrap top)` dead-reckon = a valid run start) → multi-column,
not single-column: #5's visible degeneration is gone. The **leading margin** above the
top has no wrap checkpoints yet, so it is transiently approximate — but it is off-
screen and only matters on subsequent upward scroll, where `walk!` rebuilds wrap
checkpoints; acceptable. For #6: the extent shrinkage lives above the window; Phase C
refinement re-flows the leading runs (from checkpoints, or from 0 when near start),
making the survivor top's wrap offset reflect the real shortening, so `segAnchor`'s
`shift = to − from` is the true (negative) delta → scrollOffset moves during the
segment via the WS1.3 relative correction → no post-settle snap. **Holds** (#5 = B,
#6 = B structure + C offset reconciliation feeding segAnchor).

**3. Item mutation far above the window.** Viewing index 10000, insert at index 5.
Diff gives index 5; checkpoints ≥ 5 dropped — but none exist near 5 (only in the
overscan band ~10000), so nothing local drops. `anchoredTo0` is already false (deep),
so offsets are estimates; `:approx-offset`/`:estimate` pick up `itemCount+1` and the
indicator adjusts slightly. Visible content is anchored on a visible child (`anchor-
delta`) → no pixel jump. Scrolling to top later re-flows from 0 including the new item
→ exact. **Holds.**

**4. Fast fling to top.** Flying regime: overscan off, approximate placement,
refinement suppressed (also `dir==up`). `scrollOffsetCorrection` is layout-phase and
never re-targets the active `Simulation` → fling stays smooth. As `|v|` decays →
landing regime: real materialisation + refinement re-flow toward 0, land index 0 at
offset 0, `anchoredTo0 := true`. No mid-fling offset yank. **Holds.**

**5. `animateTo` far target then immediate reverse drag.** `animateTo` far → flying
(approx, overscan off). The reverse drag cancels the `ScrollPosition` activity; the
velocity sensor (EMA + hysteresis) tracks the new direction — the engine reacts to
*motion*, never asks "programmatic or user". A slow reverse drag → resting/landing →
real materialisation at the new position via checkpoints/`:approx-offset`. Each frame's
window is re-derived (no stale window survives the flying stretch). The interrupt frame
needs no special handling. **Holds** (this is exactly why the WS2.1 gate keys off
velocity, not source).

**6. Segment playing during refinement (WS1.3 + refinement deltas composing).** During
a morph the driver is `indexed-layout!(segTween)`, which emits **only** the WS1.3
relative `seg-scroll-correction`. Flow-path refinement runs on resting/scrolling passes
*outside* the segment (the segment-start capture's internal `flow-layout!` correction
is overwritten by `from-relay!`, so it isn't emitted). So the two corrections are
**temporally separated and mutually exclusive per pass** — never double-counted. They
share the incremental/relative accumulator model (`segPrevDesired` /
`point-correction-delta` and `anchor-delta`), so switching between them across passes
is safe. The one accepted seam: at segment end, `segPrevDesired` may have drifted (edge
clamp); the first resting pass's `anchor-delta` reconciles it in one pass (WS2.3 §3.1).
**Holds**, with that documented one-pass reconciliation. Design clarification captured:
"refinement runs every dirty pass within the overscan window" means every *resting/
scrolling* pass — the segment freezes the current best offset into `segAnchor` rather
than refining live during segment-continue passes.

No scenario breaks the design; the audit produced three clarifications (cold-start
inversion, memory bound, segment/refinement temporal separation) now folded in above.

---

## 5. Integration points (by fn/field name, re-located)

| concern | fn / field | file:line | change |
|---|---|---|---|
| protocol hooks | `list-layout`/`wrap-layout`/`masonry-layout` | layout.cljd:140/158/187 | add `:renewal-point?`, `:approx-offset`; rewrite `:anchor` doc |
| aggregates | `make-env` | layout.cljd:72 | add `:avg-main`/`:avg-cross` (seed + EMA) |
| checkpoint store | new `checkpoints` field + record in `walk!` | render.cljd:47-86 / 623-677 | sparse `{idx→{:offset :state}}` |
| backward re-flow | `backfill-leading!` | render.cljd:703-741 | checkpoint-seeded forward mini-walks (reuse `flow-frames`) |
| cold-start seed | `seed-cache!` / `align-start!` | render.cljd:534-621 | invert `:approx-offset` for the start index |
| refinement / frontier | `flow-layout!` + new `anchoredTo0` | render.cljd:798-838 | emit accumulated delta; land at 0; suppression gate |
| committed prune | `prune-committed` / `commit-pass!` | render.cljd:223-258 | geometric prune; drop `committed-max-age` (204) |
| velocity + overscan | new `update-velocity!`/`resolve-overscan`/`overscan-window`/`fidelity-regime` | render.cljd `performLayout` 156 + both drivers | WS2.1 §2/§7 |
| config | `sliver-collection*` + `SliverCollection`/`CollectionSliverWidget` | sliver_collection.cljd:357 / 66 / render.cljd:1324 | `:overscan`, `:approximate-item-size` alongside `:animate` |

---

## 6. Where this doc goes beyond the earlier three (flag)

- **Memory bound on the checkpoint cache.** WS2.3 said "keep a sparse subset past
  GC" without a bound. Resolved here: checkpoints are geometrically scoped to the
  overscan window (+ one seed above the leading margin), the *same* rule as the
  committed cache. Bounded `O(overscan/min_run_height)`; sound because near-0 exact
  landing only runs when already near 0.
- **Cold-start / far-jump index seeding by inverting `:approx-offset`.** The recorded
  decision framed `:approx-offset` as a per-index estimate for the indicator and
  flying-regime placement; it did not specify that the engine must *invert* it to pick
  a start index and avoid the O(n) from-0 walk on a cold deep start. Added here.
- **Monotonicity constraint on `:approx-offset`.** Required for the above inversion
  (binary search). A strengthening of the recorded "O(1) estimate", not a conflict.
- **Segment vs. refinement corrections are temporally separated, not co-emitted.** The
  overscan doc's phrasing "refinement runs on every dirty pass within the overscan
  window" is sharpened: refinement runs on resting/scrolling passes; a playing segment
  uses the WS1.3 correction with the offset frozen at capture. No pass emits two
  corrections. This clarifies rather than contradicts.
```
