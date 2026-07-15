# WS2.3 — Progressive offset refinement & backward wrap re-flow (research)

Read-only analysis feeding WS2.3 (and the WS1.3 scroll-model choice it must be resolved with).
Scope set by the plan in `collection-scroll-anim-followups.md` §§1.3, 2.1–2.3 (+ the WS2 header
defect paragraph). Two decisions already recorded there are treated as fixed constraints:

- **WS2.1** — the `committed` cache is scoped to the overscan window, pruned geometrically
  ("left the window"), never by age; active-segment participants are never pruned mid-segment.
- **WS1.2** — `.from 0.0` on segment retrigger is REQUIRED (accept the ~2-vsync freeze).

File:line pointers are as of this reading (`src/flutter_cljd/internal/collection/…`); re-verify
before editing.

---

## 0. The defect in one sentence

There is no true from-index-0 flow geometry until you physically scroll to the start, so every
window above the first-visited region is rendered at an **estimated** main offset. Two visible
symptoms fall out (both flow-only; indexed targets are exact — `tween.cljd:130-141`):

- **#5/#6 leading degeneration.** `backfill-leading!` (`render.cljd:678-716`) is list-shaped: one
  child per main slot, `noff = head − ext`, `:cross` hardcoded `0.0` (lines 698–703). It never
  calls the layout's `:place`/state, so a wrap window above the anchor collapses to a single left
  column. Correct wrap only appears when the backfill finally reaches index 0, drops the cache,
  and the next pass re-flows from the top (`finish-flow-geometry!` `render.cljd:740-771`).
- **#6 offset settles after the segment.** `seg-anchor :to` is the survivor's frame offset in a
  `to-src` that, for a flow target, is a frozen live-only window anchored at the survivor's
  **current (approximate)** offset (`render.cljd:1139-1181`). So `shift = to − from` is
  approximate; scrollOffset barely moves during the segment and the true (shorter) extent only
  appears at the post-settle top re-flow → `ScrollPosition` ballistic-clamps the residual.

---

## 1. OPEN QUESTION — backward wrap re-flow without walking from 0

### Answer: **YES for run STRUCTURE; the absolute offset stays a prefix-sum from 0.**

A run boundary **is** a renewal point. Split the claim into the two things a placement produces:

**(a) Run structure (which items share a row, their cross positions, the run's main height) —
exactly recoverable from a run boundary, independent of everything above.**

Wrap's `:place` state is `{:x cross-end-of-current-run, :y run-start-offset, :h run-height-so-far}`
(`layout.cljd:156-170`). For a run-**start** index `i` (the first item of its run), the true
state-before-`i` is exactly `{:x nil :y y_i :h 0.0}` — which is precisely what wrap's `:anchor`
returns (`layout.cljd:172`). Seed `:place` with `{:x nil …}` and item `i` takes the `nil` branch:
`cross = 0`, `offset = y_i`; then `i+1, i+2, …` pack forward using only their own measured sizes
and `:cross-extent`. **No term reads any item < i.** The intra-run dependency (`x`, `h`) is fully
reset at each run boundary — that reset *is* the definition of a run. So:

> Given that `i` starts a run and the main offset `y_i` where that run begins, the layout of
> `[i, i+1, …]` is completely independent of items `< i`.

The `:anchor` is therefore **exact at a run boundary and only there**. Anchoring mid-run (the
current backfill's implicit assumption, one-item-per-slot) is wrong precisely because `x=nil`
forces a spurious run break — which is the visible "single left column" bug.

**(b) Absolute main offset — still global.** `y_i = Σ(heights of all runs 0..i-1)`. That prefix
sum is inherently from-0. A run-boundary cache makes each upward step *exact and O(run)* instead
of re-derivable only from 0, but the **value** of `y_i` is only as trustworthy as the pass that
computed it: mutate any item above `i` and the stored `y_i` is stale. This is why refinement +
scroll-offset **correction** (not pixel motion) is unavoidable — see §3. Structure is renewable;
absolute position is reconciled.

### Backward re-flow is a sequence of forward mini-walks, one per run

You cannot place `cf-1` in isolation: it may sit mid-run, whose start `k < cf-1` you must know to
lay the whole run `[k…cf-1]`. With a run-boundary cache you look up `k` and `y_k`, seed
`:anchor(env k y_k)`, and walk **forward** `k → cf-1` (re-measuring the real children). Then
recurse to the boundary above `k`. Total work = O(items above the window) = exactly what you were
going to materialize anyway. **Without** the cache, finding `k` has no lower bound short of index
0 — that is what forces run-boundary caching.

### Renewal-point property across the other layouts

| layout | state | renewal points | backward-reflow exactness |
|---|---|---|---|
| **list** (`layout.cljd:139-143`) | running `offset` (scalar) | **every index** | trivially exact per item (what backfill does today) |
| **wrap** | `{:x :y :h}`, resets per run | **run starts** (detectable: placed `:cross == 0`) | exact once run boundaries are cached |
| **masonry** (`layout.cljd:196-204`) | vector of per-column ends | **index 0 only** — the column vector never resets at any natural index | NOT exact from a lone `(index, offset)`; `:anchor`'s `(vec (repeat n offset))` is a genuine approximation (assumes columns level) |
| **grid** (indexed) | none — O(1) math | every index | exact everywhere; no backfill needed (why indexed targets don't show the bug) |

**Minimal extra state to make a non-renewable layout renewable:** store the layout's own `:state`
snapshot at the checkpoint. The engine **already computes and caches `:state` per index** in
`walk!` (`render.cljd:633-635`, entry `{… :state st}`), where `st` is the state *after* placing
`i` = the exact `:anchor`-equivalent for `i+1`. The only thing missing is **retention past GC**.
So:

- **wrap:** the checkpoint payload is trivial — `(index, y_i)` suffices, `:anchor` reconstructs the
  rest. Detect a boundary by `:cross == 0`.
- **masonry:** the checkpoint must carry the full column vector (the `:state`), because there is no
  cheaper renewal. Boundary = arbitrary stride (e.g. every N items), not natural.
- **list:** no checkpoint needed; any offset is a renewal point.

This is a strict generalization of `:anchor`: today it means "reconstruct a *plausible* state from
an offset (dead-reckon)"; run-boundary caching upgrades selected indices to "reconstruct the
*exact* state" by retaining the real snapshot.

---

## 2. Cross-framework pattern — estimated vs. measured reconciliation

Every mature virtualizer converges on the **same** shape:

> **Anchor on a visible item. Absorb estimation error by correcting the SCROLL OFFSET, not the
> pixels. Refine opportunistically, and suppress the correction in the one case where it would be
> visible (scrolling up into freshly measured content).**

### Flutter slivers
`RenderSliverList` dead-reckons offsets of unmaterialized children; when the running estimate is
inconsistent (the zeroth child would land at a non-zero offset, or it runs out of children before
the scroll offset) it returns `SliverGeometry(scrollOffsetCorrection: −firstChildScrollOffset)`
(or `−scrollOffset`). The correction is consumed by the viewport **within the same layout pass**:
`performLayout` returns early, the parent re-runs layout with the adjusted offset, and — crucially
— **no scroll listener fires**. `ScrollPosition.correctBy`/`correctPixels` change `pixels` without
`notifyListeners`, so content is pixel-stable while the programmatic value jumps. `applyContentDimensions`
clamps once real min/max extents are known. `CustomScrollView.center`/`anchor` pick which sliver
holds offset 0 — the mechanism behind bidirectional (chat-style) lists. Corrects **during the
layout pass, every pass an estimate is dirty**, at edge-approach (zeroth-child consistency).
[SliverGeometry](https://api.flutter.dev/flutter/rendering/SliverGeometry-class.html) ·
[RenderSliverList.performLayout](https://api.flutter.dev/flutter/rendering/RenderSliverList/performLayout.html) ·
[RenderSliverList](https://api.flutter.dev/flutter/rendering/RenderSliverList-class.html)

### RecyclerView / LinearLayoutManager
Layout is re-anchored each pass on `AnchorInfo` — a **currently visible child** and its coordinate
— never on an absolute pixel. `computeVerticalScrollOffset` (via `ScrollbarHelper`) **estimates**
the bar from the first visible child's index, its pixel position, and average item height — honest
enough for the scrollbar, not a source of truth for layout. `scrollToPosition` lands the target
just barely visible; `scrollToPositionWithOffset` pins it to an exact edge offset. Reconciliation
happens **at every layout pass**, anchored on the visible child, so off-screen size error never
moves on-screen content. [LinearLayoutManager source](https://androidx.tech/artifacts/recyclerview/recyclerview/1.1.0-source/androidx/recyclerview/widget/LinearLayoutManager.java.html) ·
[computeVerticalScrollOffset accuracy](https://medium.com/@rituel521/improving-accuracy-of-computeverticalscrolloffset-for-linearlayoutmanager-38699a9d03b)

### UICollectionView (self-sizing)
`estimatedItemSize` lets layout run before cells are measured. A measured cell returns
`preferredLayoutAttributes`; `shouldInvalidateLayout(forPreferredLayoutAttributes:withOriginal:)`
compares sizes and, if different, `invalidationContext(forPreferredLayoutAttributes:…)` reports the
invalid index paths **and a `contentSize` adjustment**. Pre-iOS-13 a measurement *above* the
viewport shifted `contentOffset` → the notorious "jump on scroll up". iOS 13+ `selfSizingInvalidation`
(`.enabledIncludingConstraints`) makes UIKit **compensate content offset automatically** so
already-visible cells stay put. `targetContentOffset(forProposedContentOffset:)` lets the layout
snap the *rest* position honestly. Corrects **on measurement**, with automatic offset compensation
to hide it. [selfSizingInvalidation](https://developer.apple.com/documentation/uikit/uicollectionview/4001100-selfsizinginvalidation) ·
[self-sizing + invalidationContext walkthrough](https://medium.com/@mfc83/uicollectionviewlayout-with-self-sizing-cells-part-3-79b1c3f0dc8e)

### Web — TanStack Virtual / react-virtuoso / CSS
`measureElement` compares measured vs `estimateSize`; for an item **above** the current offset it
accumulates the `delta` into `scrollAdjustments` and calls `_scrollToOffset(offset + adjustments)`
— it moves the scroll container by the delta to keep visible items visually fixed (programmatic
offset correction, not browser anchoring). It **skips this during backward/smooth scroll** (only
measures a buffer around the target) so a resize far from the anchor can't yank the animation.
Dynamic examples set CSS `overflow-anchor: none` to disable the browser's own scroll-anchoring so
the two mechanisms don't fight. react-virtuoso is the same idea: a size cache + a "sticky"
anchor item whose position is held while offsets around it are re-estimated. Corrects **during
scroll for items above, suppressed while scrolling backward / mid-smooth-scroll**.
[TanStack scroll management](https://deepwiki.com/TanStack/virtual/2.5-scroll-management) ·
[TanStack Virtualizer API](https://tanstack.com/virtual/latest/docs/api/virtualizer) ·
[reversed dynamic list discussion](https://github.com/TanStack/virtual/discussions/195)

### Distilled common pattern
1. **Anchor = a visible item**, tracked by index + its current offset (never an absolute pixel).
2. **Correct the offset, not the pixels.** Estimation error above the anchor is absorbed by moving
   the *programmatic* scroll value (Flutter `scrollOffsetCorrection`/`correctBy`, RV re-anchor,
   UIKit content-offset compensation, TanStack `scrollAdjustments`). Pixels of visible content do
   not move.
3. **Refine opportunistically, incrementally** — each pass/measurement applies its own delta; no
   single big snap.
4. **Suppress the correction where it would be seen** — i.e. while scrolling *up into* freshly
   measured content (TanStack's rule; iOS's pre-13 bug is exactly the un-suppressed case). Where
   each corrects: Flutter/RV **every layout pass** (edge-approach + consistency); UIKit **on
   measurement**; web **during scroll for above-anchor items, gated off on backward/smooth scroll**.

---

## 3. Recommended model for this engine

### 3.1 WS1.3 verdict — **(a) relative correction. Not (b) hands-off.**

Hands-off (b) directly contradicts the WS2 note: when the survivor set shrinks, scrollOffset *must*
advance during the segment to keep the anchored top under the viewport, otherwise the extent only
settles post-segment (= bug #6). And absolute pinning (current `set-point-correction`,
`tween.cljd:103-123`) reverts drag: `desired = lerp(from,to,t) − screen`, `delta = desired − scroll-off`
— a user drag of `+D` becomes `−D` extra correction next pass, so it snaps back ("scroll
unavailable"). Relative correction is what every framework above actually does (§2 pt 2/3) and what
this engine's own resting/flow path already does via `anchor-delta` (`render.cljd:728-738`: correct
by *how far the anchored child moved this pass*).

**Make the segment path incremental like `anchor-delta`.** Keep the `seg-anchor {:from :to :screen}`
capture; change the correction math from absolute to per-pass delta:

```
prev-desired  seeded = scroll-off at segment start   ; first pass delta = 0, no jump
desired(t)    = lerp(from, to, t) − screen
correction    = desired(t_now) − prev-desired        ; independent of scroll-off ⇒ drag survives
prev-desired := desired(t_now)
```

Now `scrollOffset = user-drag + Σ(animation deltas)`; the animation contributes exactly `(to−from)`
over the segment and **composes** with drag instead of overwriting it. This is a ~5-line change to
`seg-scroll-correction` (`render.cljd:838-848`) plus one mutable field (`prevDesired`) on the
sliver; `set-point-correction` is **reused** (same anchor, incremental accumulator) not replaced.
`anchor-delta` and the segment path then share one model (relative), which is the right invariant.

Pixel stability holds because the correction is emitted as `SliverGeometry.scrollOffsetCorrection`
inside the layout pass (`render.cljd:808-812`) — consumed same-frame by the viewport with no
listener notification, exactly the Flutter mechanic in §2. The programmatic value may jump; pixels
do not.

### 3.2 When refinement runs

Three tiers, cheapest first:

- **Every pass while an estimate is dirty (bounded).** A pass whose leading region is estimated
  runs `backfill-leading!` in its new re-flow form (§3.3) for the overscan window only (WS2.1
  scope). Each materialized real run replaces its estimated extent; the accumulated delta is
  emitted as one `scrollOffsetCorrection` (as today, `render.cljd:807-808`). Bounded because the
  overscan window is bounded — this is not a from-0 walk.
- **On approach to the start window.** When the first cached index is small and offset is near 0,
  keep re-flowing upward until index 0 lands at offset 0; the residual becomes the final
  correction so **scroll stops exactly on element 0** (the WS2.3 acceptance criterion). This is the
  existing "reached 0 away from offset 0 → drop cache, return correction" branch
  (`render.cljd:711-716`), now fed by real run geometry instead of one-per-slot estimates.
- **Suppress during active backward user scroll / ballistic-up**, per TanStack's rule (§2 pt 4):
  applying an above-anchor correction while the user is dragging *up into* that region is the one
  case a correction is visible. Defer it to the next resting pass. (See §3.4 edge cases.)

### 3.3 `backfill-leading!` — real backward re-flow via run cache

Replace the list-shaped estimator (`render.cljd:678-716`) with run-cache-seeded forward mini-walks:

- **Cache payload** (retained past GC, keyed by absolute index): `{:index :offset :state}` —
  wrap stores `(index, y_i)` (state is `:anchor`-derivable); masonry stores the full `:state`
  vector; list needs none. Populated during `walk!` whenever a boundary is crossed
  (wrap: placed `:cross == 0`; masonry: index stride). Source data already exists — `walk!` caches
  `:state` per index (`render.cljd:633-635`); this is the sparse, GC-surviving subset.
- **Backfill step:** to materialize above `cache[0]` (index `cf`), find the nearest cached boundary
  `k ≤ cf-1`, seed `:anchor(env k y_k)` (exact — `k` is a renewal point), walk **forward**
  `k…cf-1` with real `:place`/measure, prepend the produced entries (real `:cross`, real per-run
  `:offset`). Recurse to the boundary above `k`. This yields correct multi-column wrap rows above
  the anchor — kills #5's single-left-column and #6's post-settle snap (the true shorter extent is
  known *before* the segment ends).
- **Invalidation:** a checkpoint is trusted only if a forward re-flow from it re-measures children
  to the sizes that produced its `:offset`. `walk!` already truncates the cache on an extent
  mismatch (`render.cljd:616-624`); checkpoints inherit the same rule — a size change at/above `k`
  invalidates `k` and every checkpoint below it (their `:offset` prefix-sums shifted). Cross-extent
  change invalidates **all** wrap/masonry checkpoints (packing depends on it) — piggyback on the
  existing `crossCached` guard (`render.cljd:49`, drops the cache on cross change). Font/text-scale
  change is just a measured-size change → caught by the re-measure/mismatch path. Item mutation
  (insert/remove/resize at index `j`) invalidates checkpoints with `:index ≥ j`.

### 3.4 Edge cases

- **Ballistic scroll in flight during refinement.** `scrollOffsetCorrection` is layout-phase and
  does not perturb the active `Simulation` velocity — the correction shifts content-space under a
  velocity-driven pixel trajectory, so the fling continues smoothly (this is standard sliver
  behaviour). Do **not** re-target the ballistic simulation. But *suppress* above-anchor
  corrections while the fling is heading **up** into the estimated region (§3.2 tier 3) — otherwise
  the same-frame offset jump races the simulation and can look like a stutter. Defer to the first
  low-velocity/resting pass.
- **Overscroll / bounce.** During `BouncingScrollPhysics` overscroll, `scrollOffset` is out of
  `[min,max]` and pixels are rubber-banded. Emit corrections but expect `applyContentDimensions`
  to re-clamp when the true extent lands; never emit a correction that would fight the bounce
  spring at the edge. Gate refinement of the *start* window until overscroll settles, so a
  bounce-at-top doesn't get an offset yank mid-spring.
- **Scroll-indicator honesty.** The scrollbar reads `scrollExtent` from `finish-flow-geometry!`
  (`render.cljd:761-768`), today derived from `:estimate`/linear extrapolation. As real runs
  materialize, the total shrinks/grows; feed the refined extent so the thumb size/position
  converge to truth (matches RV's `computeVerticalScrollOffset` estimate → exact as items are
  measured). Accept minor thumb drift during refinement — honest and monotone-ish beats frozen.

### 3.5 Flag — impact on the WS2.2 approximate-layout protocol shape

The renewal-point analysis says the **flow layout protocol should gain two optional hooks**;
WS2.2 should design them in, not bolt them on later:

1. **`:renewal-point?`** — `(fn [env i frame state] → bool)` telling the engine which walked
   indices are exact renewal points worth retaining as checkpoints. Default derivable: list → every
   index; wrap → `(= 0.0 (:cross frame))`; masonry → false (engine falls back to a periodic
   full-`:state` checkpoint). Without this the engine can't cheaply know where the exact backward
   seeds are, and backward re-flow degrades to from-0.
2. **`:approx-offset`** — `(fn [env i] → double)`, an O(1) main-offset estimate for an arbitrary
   index, so a cold window (large initial offset, or list→wrap remap) renders at a plausible offset
   *before* any real reflow, and scroll-to-offset lands near the target (WS2.2's stated bar).
   Today only `:estimate` (total extent) and `:anchor` (state-from-offset) exist; neither gives a
   per-index offset guess. Wrap can approximate from average measured run height; masonry from
   average column height. Should be optional (falls back to linear from `:estimate`).

Neither breaks existing layouts (both optional with derivable defaults), but both are protocol
surface — cheaper to shape in WS2.2 than to retrofit. The `:anchor` contract should also be
**documented as two-mode**: exact at renewal points, dead-reckoning elsewhere — the backward
re-flow relies on the exactness half, which today's docstring (`layout.cljd:43-45`) only frames as
approximate.

---

## 4. TL;DR for the coordinator

- **Open question: YES** — backward wrap re-flow is exact without a from-0 walk *for run structure*,
  because a run boundary is a true renewal point (`:place` state resets per run; `:anchor` at a
  run-start is exact). It requires **run-boundary caching** to locate the seed; the absolute offset
  remains a prefix-sum reconciled by scroll-offset correction. Masonry has no natural renewal point
  (cache full `:state` at a stride); list is renewable everywhere; grid is exact (indexed).
- **Cross-framework pattern:** anchor on a visible item, correct the **offset not the pixels**,
  refine incrementally, suppress the correction while scrolling backward into fresh measurements.
  (Flutter `scrollOffsetCorrection`/`correctBy`, RV `AnchorInfo`, UIKit iOS-13 `selfSizingInvalidation`,
  TanStack `scrollAdjustments`.)
- **WS1.3 verdict: relative correction** — make `seg-scroll-correction` emit the per-pass animation
  delta (reuse `set-point-correction`'s anchor, add a `prevDesired` accumulator), unifying it with
  `anchor-delta`. Fixes drag-revert *and* moves the offset during the segment (WS2 requirement);
  pixel-stable via layout-phase `scrollOffsetCorrection`.
- **Flag for WS2.2:** add optional `:renewal-point?` and `:approx-offset` hooks to the flow
  protocol, and re-document `:anchor` as exact-at-renewal-points / dead-reckon-elsewhere.
