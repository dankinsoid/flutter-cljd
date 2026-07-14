# Collection engine — scroll & animation follow-ups

Follow-up plan after the v2 rect-animator (see `docs/CollectionRectAnimator.md`). Groups six
device-verified observations plus the user's architectural direction into three workstreams.
Investigation done read-only; file:line pointers are current as of this writing — re-verify before editing.

Legend: **[bug]** observed defect · **[design]** architectural direction · **[opt]** optimization ·
`file:line` from the read-only analysis.

---

## Workstream 1 — Interruption smoothness (retrigger of an in-flight segment)
Self-contained, low risk, single shared segment tween is kept (no per-item controllers).
This is what the user reported first: insert-during-insert / remove-during-insert breaks the
cell's own collapse/expand while neighbour movement stays smooth.

### 1.1 Capture completeness for the wipe (collapse/expand) — **[bug] PRIMARY (#2)**
Root cause: on a new `segment-start!` each cell's `from` is read from `committed`, but
`commit-pass!`/`child-rect` (render.cljd:204-219, 179-188) write the **full** laid extent
(`.-size child`), not the lerped **visible window** — that window lives only in `clips`
(render.cljd:961-975) and never reaches `committed`. Worse, `keyed-tween-layout`
(tween.cljd:328-341) ignores the captured `from` for wipes: entering restarts the window from
0, exiting from full. → a cell 60% open re-collapses to 0 and re-opens on the next insert;
offset stays smooth only because offset-from is captured live.
- Fix: write the **visible** extent into `committed` (read from `clips` by index, else laid size)
  plus a separate `:full-main`/`:full-cross` — box already does this (box.cljd:370-371).
- Then `keyed-tween-layout`: entering → `from* = (or from (zero-axis to))`; exiting → `from`
  with the visible extent, `full` = its `:full-*`. Wipe continues from the current window.
- Verify: insert-during-insert, remove-during-insert, insert-during-remove all continue the
  wipe from where it was (no snap).

### 1.2 One-frame freeze at retrigger — **[bug/question] (#3)**
`start-segment!` calls `.forward .from 0.0` (sliver_collection.cljd:356): the retrigger frame
renders t=0 (zero progress vs. previous frame); first Ticker tick lands next vsync → one frozen
frame. No velocity/time handoff.
- User question: is `.from 0.0` actually necessary? **Investigate first.** Options:
  (a) accept (platform-normal — UICollectionView fixes on retarget too);
  (b) carry remaining time / don't reset progress;
  (c) velocity handoff via spring-simulation instead of a fixed-duration restart (costlier,
      touches the segment controller).
- Decide after 1.1 lands — the snap fix may make the freeze far less noticeable.

### 1.3 Scroll behaviour during & after a segment — **[bug] (#4, + user pt 5)**
During a segment `seg-scroll-correction` (render.cljd:807-818) emits an **absolute**
`set-point-correction` (tween.cljd:103-123) that re-pins scrollOffset to the anchor lerp every
pass → any user drag is reverted next pass ("scroll unavailable"). Also wrong: scroll behaviour
right after a segment.
- Decide the model (open):
  (a) **relative** correction — apply only the animation delta, add the user drag on top so drag
      composes with the animation; or
  (b) **don't touch offset during the segment at all** (let cells animate, leave scroll free).
- Note: (b) interacts with WS2 — if offset must move to reveal the target window, it can't be
  fully hands-off. Resolve alongside 2.3.

---

## Workstream 2 — Approximate-any-layout virtualization + progressive offset refinement
The big architectural piece. Reframes bugs #5/#6 (list→wrap while scrolled: off-window cells
degenerate to a single left column; correct wrap only snaps in at the top; scroll offset settles
in a second pass after the segment) as one root defect: **there is no true from-index-0 flow
geometry until you scroll to the start.** Needs a design pass before coding — has an open
question (below). Independent of the animator.

Concrete current defect (both #5 and #6): `backfill-leading!` (render.cljd:648-686) is
list-shaped — one child per main slot, `cross` hardcoded `0.0` (render.cljd:669,673-674), never
runs the layout's `:place`/state. Wrap's `:anchor` (layout.cljd:172) seeds only a forward walk;
everything above the anchor is estimated, not re-flowed. #6: `segAnchor :to` (render.cljd:1146)
is anchored at the survivor's **current** offset (approximate), so scrollOffset barely moves
during the segment; the true (shorter) extent appears only at the post-settle top re-flow
(`finish-flow-geometry!` render.cljd:710-741) → `ScrollPosition` ballistic-clamps the residual =
the "offset settles after the animation" the user sees. Indexed targets are exact
(tween.cljd:130-141) and don't show this — flow-only.

### 2.1 Overscan window — layout slightly beyond the viewport — **[design] (user pt 1)**
Materialize/measure a margin of cells past the visible edges so entering cells and off-window
neighbours are real, not estimated. How wide is hard — may scale with scroll velocity.
- Research: RecyclerView `prefetch`/extra-layout-space, iOS `UICollectionView` prefetching,
  react-virtuoso / TanStack Virtual overscan, Flutter `cacheExtent`.

### 2.2 Approximate-layout protocol — **[design] (user pt 3)**
Every layout must estimate offset/extent for an arbitrary index window **without** a full
from-0 walk, so any window renders at an approximate offset. Optional user-supplied
`approximate-item-size` for better estimates.
- Scope of the estimate: it only needs to be good enough that (a) the scroll indicator is
  roughly right and (b) scroll-to-offset lands near the target.
- **Out of scope (noted):** scroll-to-**key** — nice but needs nontrivial changes.

### 2.3 Progressive offset refinement during scroll — **[design/bug] (user pt 3, 5; fixes #5, #6)**
When only an approximate window is known (large initial offset, or list→wrap remap), refine the
true offset either continuously while scrolling or when reaching the start window, so scroll
stops **exactly** on the first element. Constraint: **no visual jerk**; the programmatic scroll
value may change discontinuously (compensation is allowed to be non-smooth as long as pixels
aren't).
- Requires the leading backfill to reconstruct real rows via the layout's `:place`/state
  (a genuine backward wrap re-flow), likely caching run boundaries.
- **OPEN QUESTION:** is a correct backward wrap re-flow expressible without walking from index 0,
  given `:place` is inherently forward? May force run-boundary caching. Resolve in design.
- May touch scroll physics / the `ScrollPosition`. Research: how RecyclerView, UICollectionView,
  virtualized web lists reconcile estimated vs. measured offsets without visible jumps.

---

## Workstream 3 — Compute & diff gating — **[opt]**

### 3.1 Double (before→after) layout only when a segment will actually run — **(user pt 2)**
Guard the two-pass layout compute behind "an animation is pending", so the static case pays for
one pass.

### 3.2 Window-scoped diffing, gated on animations + key-fn — **(user pt 6)**
Run the diff **only** when animations are enabled AND a key-fn is present, and ideally only over
the current window rather than the whole item set — feasible via the before→after compute.

---

## Suggested ordering & dependencies
1. **WS1.1** first — isolated, highest-impact, fixes the primary reported bug.
2. **WS1.2 / WS1.3** right after — but 1.3's model choice interacts with WS2.3, so pick the
   scroll model with WS2 in mind.
3. **WS2** is a design-first block: settle the approximate-layout protocol (2.2), overscan (2.1),
   and progressive-offset (2.3) — including the open backward-reflow question — **before** coding.
   Bugs #5/#6 fall out of this, not out of a point patch to `backfill-leading!`.
4. **WS3** last: 3.1 depends on the before→after path; 3.2 depends on windowing from WS2.

## Cross-cutting research to feed the design agent
Overscan sizing vs. scroll speed; estimated-vs-measured offset reconciliation without visual
jumps; progressive offset refinement — across RecyclerView, UICollectionView, Flutter slivers
(`cacheExtent`), and virtualized web lists (react-virtuoso, TanStack Virtual).
