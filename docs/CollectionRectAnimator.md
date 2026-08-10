# Collection Rect Animator (v2)

> Design note for the shipped engine (`internal/collection/render.cljd`,
> `tween.cljd`). Supersedes the v1 two-path scheme (`MoveRender` position-glide +
> `tween-layout` size-morph, committed at 636cf22).
>
> §1–§8a are the rect-animator model: what is animated and why. §9 is the scroll
> model the sliver host needs to make it hold — pass modes, the rebase transport,
> the anchor. A `>` note under a heading marks a section a later one refines;
> where they disagree, the later one is the engine.

## 1. What changes and why

v1 had **two** independent animators that between them covered "things move" and
"things resize":

- **MoveRender** (`animation.cljd`) — captured each attached child's *global*
  position before a data change and glided position-only from the delta. Wrapped
  every item in a `MoveLayer`. Fired on data changes.
- **tween-layout** (`tween.cljd`) — on a resolved-layout `:id` change (list↔grid,
  grid column count) the engine drove an indexed tween that lerped two
  *index-addressed* frame sources, animating position **and** size via real
  relayout at interpolated tight constraints. Fired on `:id` changes.

Plus a third path — **insert/remove** via `SizeTransition`+`FadeTransition` on a
per-item controller (diff-classified).

Three code paths, three trigger conditions, three notions of "from". They compose
badly: a simultaneous data + layout change had to *snap* because neither path
could own both. Move was global-space (so it accidentally fought scroll in edge
cases and could not animate size). Size-morph was index-addressed (so a shuffle,
where a key jumps index, was invisible to it).

**v2 collapses all three into one keyed, engine-owned per-cell rect animator.**
The engine always animates each visible cell's content-space RECT (offset + size)
from its last-committed rect toward its current target, whatever the cause — data
shift, layout swap, grid column change, one sibling resizing. There is no
`:id`-change detection. `:move` is the single public knob; `:layout` is removed.
The diff shrinks to enter/exit classification only.

## 2. Why this is *one* mechanism, not a merge of three

The unifying observation: **position-move, size-morph, and layout-swap are all the
same event — a cell's target rect differs from where it currently sits.** Once the
engine owns a persistent `{key → last-committed-rect}` and lays every cell at
tight constraints derived from a *lerped* rect, all three fall out for free:

- **Move only** — from/target differ in offset, agree in size. Lerp of equal sizes
  is bit-identical (`lerp(a,a,t)=a`), so the child gets *identical* tight
  constraints every tick. A cell laid at tight constraints is its own relayout
  boundary, and `RenderObject.layout` early-returns when constraints are unchanged
  and the child isn't dirty ⇒ `performLayout` never runs; only `layoutOffset`/paint
  update. This is exactly as cheap as the old `MoveRender`, with no proxy widget.
- **Size morph** — from/target differ in size ⇒ tight constraints change each tick
  ⇒ the cell reflows for real, *per cell*, automatically. Text rewraps. This is the
  old size-morph, now per-cell instead of whole-layout.
- **Layout swap / column change** — just a target-rect change for every visible
  key. No special trigger.

So "just move vs. move+relayout" is **free, not a branch we write** (constraint 1).
The only requirement: lerp of equal endpoints must be bit-identical so stable-size
cells get constraint *identity*, not constraint *jitter*. Component-wise
`lerp(a,b,t)=a+(b-a)t` already satisfies `lerp(a,a,t)=a`.

## 3. Content-space, not screen-space (constraint 2)

Rects are stored and animated in **content space** — the offset along the scroll
axis, the same coordinate `layoutOffset` lives in — never global/screen coords.

Reason: global coordinates move under ordinary scrolling. If `from` were a captured
screen position (as `MoveRender` did), every scroll would look like a move to
animate. Content-space offsets are invariant under scroll, so a plain scroll pass
(constraints change, no widget update) lerps `from→target` with `from==target` ⇒
no animation. This is also why the persistent state is keyed content-space rects,
not captured global positions.

## 4. Resting vs. animating: a mode switch, and why it must be one

> **Refined by §9a.** The switch is a FOUR-value axis — `:resting`, `:settled`,
> `:capture`, `:segment` — decided once per pass and read through a capability
> table. The two-value framing below is *why* it exists, not its shape.

A flow layout (list/wrap/masonry) *self-sizes* children — it measures each at
loose constraints to learn its natural size, and that measurement is how the engine
detects data changes (self-healing cache: re-measure, compare extent). Tight
constraints *impose* a size, so under tight constraints a child can no longer tell
us it wants a different size. These two needs are irreconcilable in one pass —
hence a mode switch (already the shape of v1, kept):

- **Resting** (no animation in flight): the existing flow/indexed drivers run
  unchanged. Flow measures naturally, self-heals, corrects scroll. The engine
  records each visible cell's committed rect `{key → rect}` at the end of the pass.
- **Animating**: tight lerped constraints, **no natural re-measure**. The target
  rects were captured at animation start (one natural measure pass to seed them)
  and the tween interpolates committed→target with a single global clock. At `t=1`
  we settle back to resting, which re-measures naturally so the next data change is
  seen again.

The natural-size measurement therefore happens **once per segment** (at start /
retarget), not per tick — this is what lets a stable-size cell reach
zero-relayout per tick during the glide.

### 4.1 The target is captured by the ENGINE, from the NEW tree (critical)

The target must be measured **after** the widget update, over the *new* children —
and only the engine's first layout pass has that tree. At `didUpdateWidget` time the
render tree still holds the **old** children; the new/inserted cells aren't built or
laid out yet. So the host **cannot** build the target (a v1 mistake we don't inherit:
v1 got away with a host-built target only because its morph fired on a pure layout
swap with *unchanged data* — same cells, so the old tree measured the right sizes,
and it explicitly *snapped* on any simultaneous data change).

For v2, where a data change (flow insert/remove/shuffle) *is* animated, a
host-built-from-old-tree target mis-keys: the inserted cell isn't measurable, and the
survivors' new positions aren't derivable from the old sequence. Therefore:

**The engine captures the target on the first pass of a segment, from the new tree,
keyed.** On the gen-change pass the new children are attached; the engine runs its
normal driver's measurement (flow: dry-project the new window under the target
layout; indexed: exact math) to produce the target frames, freezes them, then lays
this pass at `t≈0` (= `from`, no jump). Because the target is captured from the new
tree in the new index space, it aligns byte-for-byte with what the engine queries,
and keying both sides fixes shuffle alignment for free.

## 5. State on the render object

```
committed : {key → rect}   ;; rect last laid this cell got; updated every pass
segTween  : layout         ;; the built keyed lerp; owns the frozen `from` map,
                           ;; the frozen `to-src` and the segment's :domain
segAnchor : set-point      ;; §9d; nil outside a segment
viewAnchor: {:idx :key :off :extent :frac}  ;; §9c — the viewport-top truth
curGen / segGen : int      ;; host bumps curGen on each (re)start; segGen is the
                           ;; gen the engine last BUILT segTween for
```

`rect = {:offset :cross :main-extent :cross-extent}` in content space.

The frozen `from` map is *not* a render-object field: it is built at segment
start and reaches the tween as a closed-over local, so it cannot drift from the
tween that reads it. `target` likewise never outlives the pass that queries it.

**Keying.** The render layer is index-addressed, but animation identity is by key.
The host passes `key-of : (fn [^int index] → key)` (closure over data + `:key-fn`).
`:key-fn` is **optional**: with it, a cell's rect history follows its data identity
across shuffles/inserts; without it, `key-of` is the index — a first-class mode,
not a degraded one (a positional collection animates position/size correctly, it
just can't tell a shuffle from a set of coincident moves). The engine maps `index →
key-of → key` per child. During removes the index space is shadow-data (live +
dying), so `key-of` reads it, exactly like the reorderable `key-of` already does.

**Persistence & GC** (constraint 3). `committed`/`from` survive child GC: a cell
scrolled off and back keeps its rect history under its key. They are pruned when a
key leaves the *data* (gone after an enter/exit settle) or has been untouched for
N passes (LRU-ish bound so an infinite feed doesn't leak). Never-seen or
long-pruned keys have no `from` ⇒ they appear at target (no glide), which is the
correct behavior for a cell that shuffles in from far off-screen.

## 6. Segment lifecycle (jump-free retarget — constraint 5)

The host owns the `AnimationController` (needs a vsync) + curve and the *trigger*;
the render object owns the target capture and the clock's *effect*. On any
`didUpdateWidget` where a cell could have moved (data changed, layout changed, count
changed), the host bumps `gen`, calls `controller.forward(from: 0)`, and hands the
engine the plain target layout + the animation + `gen` + the `entering`/`exiting`
key sets + `key-of`. It builds **no** target (see §4.1). Then:

1. First layout pass of the new segment sees `gen` changed ⇒ **freeze**
   `from := committed` (the currently displayed rects — which, mid-previous-glide,
   are the interpolated rects; jump-free retarget) and **capture the target** from
   the new tree (§4.1): `to-src` = a frame source over the target layout (indexed
   math, or a dry-projection of the new flow window), plus `from-extent` = the
   scroll extent at segment start. The engine builds the keyed lerp
   (`keyed-tween-layout`) internally and uses it as its layout for the segment.
2. Every pass: `displayed = lerp(from[key], to-src(index), t)` keyed via
   `key-of(index)`, lay at tight constraints from `displayed` size, set offset;
   `committed[key] := displayed`.
3. At `t=1`: `finish` — host drops animating mode, engine settles to resting.

No `capture-positions!` step is needed (v1 needed it because `MoveRender` read
screen coords lazily on paint); the engine already holds content-space `committed`
continuously, so freezing `from := committed` is a map copy.

## 7. Windowing during a segment (union — constraint 4 & 6)

Which indices to materialize while animating = **union** of:

- **target window** — `first/last-index` of the `to-src` over `[ws, we]`, and
- **attached children still overlapping `[ws, we]`** — any currently-attached cell
  whose laid (lerped) content-space rect still overlaps the window is kept, even if
  its index falls outside the target window. (Implemented driver-side as
  `widen-window!`, since the keyed `from`/`committed` maps have no cheap
  index-inverse; the attached set is the practical, bounded stand-in.)

The union makes a cell leaving the window under the new layout, or entering it,
materialize for the whole glide instead of popping at an edge. A cell whose target
is visible but whose `from` was far off-screen (a long shuffle from a pruned key)
has no `from` ⇒ appears at target; we do **not** drag the entire list across the
viewport. That bounds shuffle cost to cells with live history.

Because a simultaneous data + layout change is *also* just a per-key rect change,
it no longer needs to snap (constraint 6): both hosts animate it cleanly.

> **Superseded by §7a for indexed layouts.** The `widen-window!`/attached-children
> union above was the first cut. Runtime testing (2026-07-13) showed it can't hold:
> it leaves the *from*-visible cells unmaterialized during a large non-monotonic move
> (shuffle) — viewport-top goes empty → overscroll bounce — and it interacts fatally
> with the reserve-slot exit of §8 to make indexed remove animate *nothing*. §7a is
> the canonical model; it subsumes both windowing and enter/exit for indexed hosts.

## 7a. Segment mechanics, canonical: pure-new-data target + overlay exits

Two runtime bugs turned out to share one root cause:

- **Indexed remove animates nothing.** The host reinserts the dying cell at its
  *old index* in shadow-data (to keep it in the tree for the collapse). That leaves
  the index space unchanged, so `indexed-frame-source` yields `to-frame(i) ==
  committed(i)` for every live cell ⇒ **stay-glide is dead** (nothing below slides
  up), and the collapse shrinks inside a *reserved slot no neighbor vacates* ⇒ reads
  as "vanish + snap." A flow layout reflows the collapsed cell out and glides;
  indexed never does, because the reinsert freezes the indices.
- **Shuffle bounces.** The segment materializes only the *target* window; the
  *from*-visible cells (old on-screen keys, scattered to new indices by the shuffle)
  are never laid at their committed positions ⇒ viewport-top empty ⇒ bounce.

Both are the same defect: **cells that are leaving the on-screen window still occupy
the indexed target's index space.** The fix removes them from it.

### Principle

> During a segment the indexed target must present **zero reserved extent** for a
> leaving cell, so survivors get their new (shifted) target frames and **stay-glide
> revives on its own** — remove, shuffle, insert alike. But a leaving cell leaves for
> two different reasons, and each keeps that principle a **different** way (below).

This is why §7's brute-force `[min..max]` span is unnecessary: a leaving cell never
needs its scattered real index materialized. It is either a zero-extent shadow slot
(data-remove) or an off-index overlay (reindex); arriving cells are the ordinary
**contiguous** window near `scrollOffset`.

### Two exit/enter semantics — different reasons, different mechanisms

A cell crosses the window edge for two different reasons; conflating them was the v1
mistake. They use **different** mechanisms, each validated by a real framework:

1. **Data collapse** (key removed from / added to the *data*) → **deferred-delete +
   remapped target, in-tree — NOT an overlay.** Keep the dying cell as an ordinary
   shadow child (it is built by the user builder, hit-tested, GC'd, contiguous — no
   adaptor fight). The fix is in the *target source*, not the tree: wrap `to-src` so a
   live shadow index `i` answers `to-frame(data-index i)` (pure-new-data math ⇒
   survivors glide) and a dying slot answers a **synthesized gap frame** (offset = the
   next survivor's new frame, cross from committed, `main-extent → 0` via the existing
   `:exiting` path). Net: the dying slot presents zero reserved extent to survivors
   while the cell itself shrinks in place. This is exactly Flutter `AnimatedList`'s
   `_itemIndexToIndex` remap; the un-remapped version is `AnimatedGrid`'s
   neighbors-don't-glide bug — i.e. our 2b symptom. The `collapse-wrap` of §8 owns the
   shrink; **enter** is the same, `0 → final`.
2. **Viewport slide** (key exists both before and after; only *reindexing* pushed it
   off-screen — shuffle, or a layout change moved its slot past the edge) → a
   **full-size slide** from the committed rect past the near window edge, **no size
   change** — the viewport clips it. This is `UICollectionView`'s
   `finalLayoutAttributesForDisappearingItem`. Deferred-delete cannot apply: the key
   is live at a real off-window index. The mechanism shipped IN-TREE, not as an
   overlay: the child list tolerates index gaps, so the cell stays materialized at
   its scattered index and the drivers never bridge the gap (§7b).

### Classifier (per cell, per segment)

`from = committed[key]`, `target = to-src(index)`, both content-space; window `[ws,we]`.

| old→ / new↓ | key ∉ new (removed) | key ∈ new, target ∉ window | key ∈ new, target ∈ window | key ∈ new, from ∉ window |
|---|---|---|---|---|
| **key ∈ old** | **exit-collapse** (→0, in place) | **leave-slide** (overlay→edge) | **stay-glide** (lerp from→target) | **arrive-slide** (offscreen from→target) |
| **key ∉ old** | — | — | **enter-collapse** (0→final) | enter-collapse (0→final) |

- **exit-collapse** — deferred-delete + remap (semantics §1): in-tree shadow child,
  zero-extent gap frame; no overlay, no adaptor fight.
- **leave-slide** — in-tree kept child at its real (scattered) index, full-size
  slide to the near window edge, viewport clips. No overlay needed (see §7b).
- **arrive-slide** — **zero new machinery**: an ordinary indexed child of the
  contiguous target window, laid at `t≈0` on its off-screen `committed[key]`, lerping
  to its on-screen target, clipped by the viewport until it crosses the edge.
- **stay-glide** and **enter-collapse** are unchanged from §6/§8.

### Shared core (both hosts, no duplication)

Pure functions in `tween.cljd` express all of the above once; both hosts consume
them (§10). (The five-class matrix above is a design vocabulary — in code the
classification is carried implicitly by the `entering`/`exiting` key sets plus the
`augment-from-edges`/`augment-to-edges` builders; a standalone `classify-cell` fn
existed as a reference implementation and was removed as dead code, 2026-07-15.)

- `(shadow-entries n key-of exiting committed)` → the shadow index space as
  live/dying entries, fed to the remap below (each host supplies its own n/key-of).
- `(shadow-frame-source to-src shadow-entries)` → a frame source with the collapse
  remap baked in (same shape as `frozen-frame-source`, so each host's driver consumes
  it unchanged).

  **The collapse remap wraps BOTH layout branches, not just indexed.** A layout
  branch's only job is to emit `to-src` = a frame source over PURE-NEW (live-only)
  data — survivors gap-closed, indexed by pure-live index. The universal
  `shadow-frame-source` then remaps shadow indices onto it (dying slots → zero-extent
  gap, live shadow index → pure-live index). The branches differ only in how the
  live-only source is produced: **indexed** = exact O(1) math over live indices
  (`indexed-frame-source`); **flow** = a live-only *re-flow* (`live-only-flow-window`
  in `render.cljd`) that, after the flow measure pass, drops the dying cells, reindexes
  the survivors to contiguous pure-live indices, and re-places them via `flow-frames`
  (a real re-flow, correct for wrap/masonry where offset-subtraction would be wrong),
  frozen as a `frozen-frame-source`. Indexed is the O(1) special case of flow. A flow
  branch that instead froze the raw SHADOW snapshot (dying measured at full size)
  gives survivors `to == from` ⇒ no glide — the bug this fixes.
- `(edge-slide-block entries ws we)` → the edge-clamped slide targets: each
  edge-group translated as a block that keeps its relative layout, abutting the near
  window edge (superseded the per-cell `slide-out-frame`, removed 2026-07-15).

Host-supplied: **sliver** — `[ws,we]` from constraints, the in-tree kept-child
management of §7b, scroll correction (exists), GC exemption for sliding keys. **box** —
infinite window (so the classifier degenerates: no slide classes, collapse only),
`key-of` over child widgets, own-size lerp (exists). Box therefore **gains enter/exit
for free** from this core, but wiring its host side is a separate step (2d); 2c ships
the core + sliver.

### What this removes / revises

- 2b's indexed `rebuild-shadow!` **drop-dying** branch is **reverted**: the dying cell
  returns to shadow-data (in-tree), and the collapse comes from the `shadow-frame-source`
  remap instead. (2b proved survivors glide under pure-new-data indices; the remap
  achieves that *while* keeping the cell in-tree to shrink.)
- `build-shadow-data`'s old reinsert kept the dying cell at a *full-extent* old slot;
  the remap makes that slot **zero-extent** to survivors. No reserved gap.
- §7's `widen-window!` / attached-children union (indexed hosts) — leaving cells are
  either zero-extent shadow slots or off-index overlays, arriving cells the contiguous
  window, so no scattered union.

### Span cap (unbounded shuffle)

A full-list shuffle makes *every* on-screen key a leave-slide; the kept set is
bounded by the viewport (≈ one screenful of leaving + one of arriving) and no
bridge cells are built to reach scattered indices (§7b), so virtualization holds. If the from-visible **and** target-visible sets together
exceed a threshold (e.g. a viewport can't host them), the segment **snaps** that
change rather than animating — logged, not silent.

## 7b. Leave-slide mechanism: in-tree at the real index (overlay sketch retired)

The keepAlive-overlay sketch proved unnecessary: the adaptor's child list requires
only **ascending** indices, not contiguous ones, so a leaving cell stays materialized
at its real (scattered) index while the drivers simply never bridge the gap. Shipped
(2026-07-14), universal for flow and indexed targets:

- **`to` = edge slide** (`augment-to-edges`, tween.cljd — the mirror of §7c's
  `augment-from-edges`): at segment start, every attached cell whose committed rect
  overlaps the current window `[ws,we]` but whose index falls outside the **union of
  the current and end-of-segment target window index spans** gets
  `to` = its target frame slid past the near END-window edge (`edge-slide-block`
  over `[end-ws, end-we]`) — index above the span → `:trailing`, below →
  `:leading`. The end window is the current one shifted by the
  §9d set-point delta (`desired(1) − desired(0)`), so a cell the viewport merely
  scrolls away from keeps its real target and is NOT edge-dragged. Exiting
  (data-removed) keys are skipped — they exit-collapse in place (§7a).
  `keyed-tween-layout` takes the map as `:leaving` and lerps committed → edge at
  constant size (pure translation: no clip, no resize, and it overrides even an
  exact far indexed frame so a shuffle never drags a cell across the list).
- **Kept materialized**:
  - *flow capture pass*: `post-gc!`'s kept range is widened over the attached
    committed-overlapping cells (`widen-window!`), so the capture no longer GCs them
    at segment start; `live-only-flow-window` restricts its re-flow to the
    cache-backed (walked) window so kept cells' stale sizes can't corrupt the
    pure-live source; `from-relay!` then lays them at `t≈0` = committed (in place —
    this is what fills the viewport top on a shuffle and kills the bounce).
  - *indexed passes*: GC bounds were already widened over attached overlapping
    cells; the forward loop now creates children ONLY inside the true target window
    `[first-index0, target-last0]` and JUMPS across index gaps to (re)lay kept
    out-of-window children — no bridge cells are ever built to reach a scattered
    index, so virtualization holds (the kept set is bounded by one screenful).
- **GC'd once gone**: `committed` is re-stamped with the lerped rect every pass, so
  when the glide crosses the window edge the cell drops out of the widened bounds
  and ordinary GC collects it; returning to resting clears `:leaving`.

Known edges (accepted, self-healing at settle): scrolling mid-segment can bring a
leaving key's real index back into the window — the window path then lays it at its
frozen edge lerp until the segment settles (one snap at rest). A flow end-window
that falls outside the frozen capture snapshot (large §9d shift on a mid-scroll flow
switch) falls back to current-window classification — boundary cells may edge-slide
and re-enter at settle. A leading-scattered kept cell in the flow capture can still
anchor the walk at a stale offset (pre-existing behavior, unchanged by this).

## 7c. `from` = committed ∪ edge slide (glide cells entering the viewport)

`from` was sourced **only** from `committed` — the rects of cells that were on-screen
last pass. A cell that was *below* the old viewport and lands *inside* the new one
(e.g. grid2→grid4 packs more per row, pulling formerly-below cells up into view) has
no `committed` ⇒ `keyed-tween` returns its target directly ⇒ it **pops in without a
glide**. Wrong: it should slide in from the edge it came in through.

Fix — **universal**, layout-agnostic: an off-screen cell entering the viewport only
needs the SIDE it came from, not its exact old position. So at segment start
`from[key]` falls back to the cell's own **target rect slid to the near viewport
edge** (an arrive-slide) when `committed` is absent. The edge is chosen by the cell's
index relative to the on-screen (committed) index range: an index **above** the
committed max came from **below** (bottom edge); **below** the committed min came from
**above** (top edge); interleaved (or nothing on-screen) picks the nearer edge. This
is `augment-from-edges` in `tween.cljd`, wired unconditionally in `segment-start!`
(both the flow and indexed branches).

Because the `from` is derived from the cell's *own target frame* — which every layout
kind already produces — this works **identically for flow and indexed** (list, wrap,
grid, masonry). There is no old-layout math and no `prevLayout` field: the previous
resolved layout is neither retained nor consulted.

`augment-to-edges` (§7b) is the exact mirror on the LEAVE side: a cell bound
off-screen only needs the side it leaves through, so its `to` = its own **committed**
frame slid past that edge. Enter and leave share `edge-slide-block` and the same
index-vs-span edge choice; the enter side keys off the attached `on-screen` set, the
leave side off committed-overlap + index-outside-span.

Tradeoff vs. the earlier indexed-only approach (which sourced `from` from the exact
old-layout frame): the glide now starts at the viewport edge, not at the cell's exact
prior position — a shorter path, but smooth, and universal across every layout kind
and across simultaneous data + layout changes (no old data-index to disambiguate).

## 8. Enter / exit (constraint 4)

> **Exit is revised by §7a; the collapse *mechanism* is revised by §8a.** "Exit"
> splits into *data collapse* (key left the data → shrink toward its gap point) and
> *viewport slide* (key reindexed off-screen → full-size slide). The widget-side
> `ClipRect > OverflowBox` wrap below was the first collapse driver; §8a moves that
> clip into the engine (the widget can't know an *insert*'s target size at build).

The diff still runs — but only to classify **enter** (key present in new data,
absent in old) and **exit** (key absent in new data, present in old). Note the
distinction the diff must preserve: a key absent from the *window* but present in
the *data* is a scroll-off, **not** an exit.

- **Exit**: the engine drives the cell's rect main-extent `final → 0` (and holds
  its cross); when it reaches ~0 the dying entry is dropped (shadow-data, as today).
- **Enter**: the engine drives main-extent `0 → final`.

Cells rarely lay out correctly at 3px or any non-final size, so the **content must
not reflow through intermediate sizes**. The engine owns the *outer* rect; the
content is wrapped:

```
ClipRect > OverflowBox(final main-extent, alignment=start) > (custom :builder or content)
```

`OverflowBox` feeds the child its **final** constraints regardless of the shrinking
outer box, so the content keeps constant constraints (no per-tick cascade relayout)
while the engine's clip reveals/hides it. Per-section policy `:clip` (default) or
`:scale`; custom `:insert`/`:remove` `:builder` (opacity/transform) is preserved,
with *size* owned by the engine (the builder no longer drives `SizeTransition`).

## 8a. The engine owns the enter/exit visual (paint-clip)

The widget-side `ClipRect > OverflowBox` of §8 has a fatal gap: **it must know the
cell's stable full size to pin the content, and for an *insert* that size is the
target — which the host does not have at build time** (the target is captured by the
engine during layout, *after* build; §4.1, §timing). A `remove` works only because
its stable size is `committed`, which the host *can* snapshot. So collapse ownership
belongs where the size lives: the engine.

**Mechanism.** During a segment the engine lays each animating cell at its **full
stable main-extent** — `committed` for exit, `target` for enter, `lerp(from,to)` for a
stay whose size changes — never the shrinking lerped extent. So the content lays out
once at a real size and never reflows to its natural (e.g. text) height. The cell's
**scroll contribution** (what closes the gap for neighbors) and its **painted window**
are then decoupled at PAINT: the engine overrides `paint` and, for each animating
cell, `pushClipRect` to the cell's **lerped visible window** and paints the full-size
child clipped to it (translated so the reveal anchors at the correct side).

- **exit-collapse** — full size `committed`, clip window shrinks `committed → 0` toward
  the **gap point** (where survivors converge: main-axis for a list, the neighbour
  slot for a grid), anchored at the gap side (**start**, not center).
- **enter-collapse** — full size `target`, clip window grows `0 → target` from the same
  side. No build-time size needed — the engine holds `target`.
- The widget-side `collapse-wrap` (ClipRect>OverflowBox) is **retired**; `:insert`/
  `:remove` `:builder` (opacity/transform) still wraps the content, size owned by the
  engine.

Leave-slides (§7b) don't need this clip — a slide never resizes, so their
`full-main-extent == main-extent` and they paint unclipped; the viewport itself
clips past the edge.

## 9. Scroll correction: the rebase model

A **rebase** is a coordinate-frame event: the engine decides that the content the
viewport is anchored to now lives at a different content-space offset (new-frame
= old-frame + Δ), and `scrollOffset` must move by the same Δ so the visible pixels
do not. Content above the viewport collapsing (a remove, list→grid packing the
region before the anchor), a leading estimate turning out wrong, a cache re-seed
after a far jump, a segment's set-point — all of them are that one event.

Flutter offers exactly one valve for it, `scrollOffsetCorrection`, and it is
re-entrant: a correcting pass is discarded and re-run, up to ten times per frame.
Everything in §9a–§9f exists to make sure a Δ is produced by one owner, travels by
one route, and is consumed exactly once.

### 9a. The pass mode and the capability table

A layout pass is one of four modes, decided ONCE per `performLayout` (`pass-mode`,
a pure kernel over the dispatch it feeds) and stored on the render object:

| mode | when |
|---|---|
| `:resting` | no host clock running |
| `:settled` | the clock plays out but the engine holds no segment model — never built, or the window left the frozen `:domain` (§9f) |
| `:capture` | `segment-start!`'s target-materialization pass |
| `:segment` | a live pass over the frozen segment tween |

`:settled` is not a nicety: the resting drivers run under a playing clock, and
`tweenAnim` (a clock exists) and `segTween` (a segment model exists) are *not*
interchangeable. Every gate that re-derived segment-ness from those two fields was
a place the two could be confused, and each confusion shipped as a bug.

`pass-caps` maps each mode to twelve named capabilities — `:emit-rebase?`
`:lead-emit?` `:overscan?` `:prune-commit?` `:measure-feed?` `:tripwire?`
`:canonical-assert?` `:epilogue-asserts?` `:union-window?` `:set-point?`
`:seg-tail?` — and `pass-allows?` debug-asserts that the capability HAS a row, so
a newly added feature cannot inherit an accidental OFF. The table is the single
owner of segment-ness: no engine code reads `tweenAnim`/`segTween` as a gate, only
as state ownership (listener swap, tween value, settle).

### 9b. Transport: one accumulator, three arms, exactly once

- **`passRebase`** — a per-pass double, reset each `performLayout`, written only
  through `rebase+!`. Every producer (cross-layout re-anchor, seam measurement,
  top underflow, exact landing, anchor follow, segment tail) adds its Δ and
  continues working in the shifted frame.
- **Exactly one arm consumes it** at pass end, selected by the mode:
  - **emit** (`:resting` / `:settled` / `:segment`) — `emit-correction!` writes
    `scrollOffsetCorrection`. It is the single correction-writing path and
    debug-asserts `:emit-rebase?`.
  - **value** (`:capture`) — the emit arm is structurally unreachable, and the
    driver returns `{:rebase :extent :anchor-frame}` as plain values over no
    mutable state. A capture correction would be overwritten by `from-relay!`,
    which lays the same pass at t≈0; the delta would vanish silently.
  - **absorb** (`segment-start!`) — the returned rebase is folded into the
    set-point, whose incremental corrections deliver it over the segment. The
    absorb arm is TOTAL: when no anchor resolves, the rebase is absorbed as a pure
    shift over the first placed slot, so a produced rebase always has a channel.
- **Exactly-once assert** (debug, epilogue): a resting pass's emitted correction
  equals the produced rebase; a capture pass's displayed geometry carries none.

`reseedCause` (`:count-change` | `:cross-layout`) is the only cross-pass remnant:
it tells the next `seed-cache!` that the frame must be re-derived from the anchor,
and is cleared on consumption.

### 9c. The Δ-epilogue and the anchor latch

The viewport-top truth is `viewAnchor {:idx :key :off :extent :frac}`, sampled from
the DISPLAYED frames after EVERY pass — mid-segment included, because the user can
scroll during a segment and the displayed frame is the same frame `committed`
records. `:frac` is the fraction of the anchor item consumed above the viewport
top, and `:key` is what survives a count change (keyed reconciliation has already
moved the element; the anchor is re-resolved by key, never by index).

A resting flow pass therefore ends with ONE producer, not five:

1. `backfill-leading!` MEASURES the extent above the anchor and returns one
   leading delta. Measurements **supersede, they never accumulate**: no step moves
   the window, so every step compares the same rigid frame against its own
   reference and each report is the WHOLE displacement restated one run closer to
   the origin. Summing them translates the window by a multiple of its own error.
2. `translate-window!` moves the whole laid window rigidly, once, and drops the
   cache and checkpoints rather than rewriting them (a flow state is layout-opaque
   and carries absolute offsets). The band's frontier is kept and `baseState`
   re-derived at its translated offset.
3. The epilogue emits the ANCHOR's own displacement over the whole pass — not the
   raw measure. A pass with no anchor to hold falls back to the measure.

**The latching invariant**: aggregate drift never moves a latched anchor. A
synthesized checkpoint sits at `:approx-offset`, which is a function of the
measured EMA that this very window feeds — honouring its seam under a latched
anchor makes the frame chase itself one `1 − ema-alpha` step per pass until the
viewport's ten-cycle budget runs out. A pass that REASSIGNED the anchor (far jump,
cold start) still honours it: there the estimate is the truth it just chose.

A flow pass's frame is a run-chain from wherever its walk STARTS, so it is not
reproducible from a different start; re-deriving the window inside the translation
levels a masonry's columns and fabricates a wrap run break. The rigid translation
stays and the anchor absorbs the re-derivation — which is why the anchor must be
an IDENTITY (`viewAnchor`), not "the first cached entry overlapping the top": a
bsearch reference re-picks the run's lowest row-mate every pass, so a re-pack that
splits the run holds the run and drops the child.

**The leading side is clampable, the trailing side is not.** `SliverConstraints`
carries `precedingScrollExtent`, so "content space is pinned at the leading edge"
is a bound the sliver can compute (`correction-floor`): below it the viewport parks
`pixels` out of range and no further pass runs to walk it back. It says nothing
about the slivers that FOLLOW, so `maxScrollExtent` is not knowable here and a
self-computed trailing clamp would fire early in a composed viewport. The
accumulator is trimmed before the emit arm reads it, so exactly-once survives.

### 9d. The segment set-point preserves the anchor's FRACTION

During a segment the scroll offset is derived from the anchor's trajectory:

```
desired(t) = lerp(from-off, to-off, t) + frac · lerp(from-ext, to-ext, t) − screen
```

`tw/point-desired` is the one formula; `point-correction-delta`, the segment tail
and the leave-slide shift all route through it. What stays at the viewport top is
the content point at `frac` of the anchor item, **not its top edge** — an anchor
that resizes across a morph (a grid row becoming a list cell) no longer slides out
from under the viewport. With an unchanged extent this is pixel-identical to a
rigid `to − from` shift. `screen` is the residual
`(from-off + frac·from-ext) − scrollOffset`, which keeps `desired(0)` exact for a
clamped `frac` and for a drag landing between the sample and the segment start.

The corrections are the ANIMATION delta only — they never read the live
scrollOffset, so a concurrent drag composes (`scrollOffset = drag + Σcorrections`)
and the emitted deltas telescope to exactly `desired(1) − desired(0)`.

Both endpoints come from frames the pass actually placed:

- `to-off`/`to-ext` from the capture walk's own placement (`anchor-frame`), read
  before the band is trimmed. A `to` slot and a `to-src` frame that can disagree
  in-window is the bug; agreeing by construction is the fix. When the walk did not
  cover the slot the value is nil — a set-point built on an estimate is worse than
  none.
- If the anchor's key is `exiting`, the anchor re-picks the first surviving
  attached slot at/after the viewport top before freezing.

The host clock can complete before a t=1 pass runs, leaving `desired(1) −
segPrevDesired` un-emitted; `consume-seg-tail!` publishes that residual through the
accumulator on the first resting pass, which is `segAnchor`'s explicit end of life.

### 9e. An unanswerable query resolves to the source's own edge

One rule, three queries. A frozen snapshot has no opinion outside its span, and
mapping "no opinion" to index 0 or offset 0 is what fabricates inverted windows,
band-pinned segments and O(gap) walks:

- **`:first-index` / `:last-index`** clamp a nil answer to the snapshot's own edge
  (past the end → its last index; above the start → its base), so the window is
  always well-ordered and inside the snapshot. `indexed-layout!` additionally
  refuses inversion — debug assert, release clamp.
- **the per-child frame** (`parked-frame`): a cell with neither committed history
  nor a target frame is parked at the source's nearest EDGE frame, collapsed on the
  layout's own `:collapse-axis` so it reserves nothing, carrying that frame's
  extent as the stable full size it is laid out at (§8a). Offset 0 would strand it
  at the scroll top, where it becomes `firstChild` and costs the next capture an
  O(gap) walk. Parking is a substitution of `to` before the enter/exit/leave cond,
  so all four branches inherit it.
  The edges are probed at offsets a bounded source must cover — its start and its
  own reported extent — never ±infinity, which would reach a layout's index math.

### 9f. A frozen segment has a validity domain

`segment-start!` stamps the tween with `:domain` = the union of the t=0 and t=1
window spans. A frozen segment never extrapolates: when the current window lies
more than one window extent outside the domain, `settle-segment!` drops
`segTween`/`segAnchor`/`leaving`, sets `segGen := curGen` and dispatches to the
resting drivers, whose far-jump inverse seed takes over — O(window), band follows
pixels. The host clock keeps playing to its end; that pass is `:settled` (§9a) and
the engine settles silently, with no callback and no coupling back to the host.

Cache is an accelerator, never the algorithm. The far-window seed is decided by
GEOMETRY alone — leftovers that happen to be attached must not veto an inverse
seed — and the per-pass work bound is O(window + overscan) unconditionally,
enforced by `assert-materialization-bounded!` in debug.

## 10. Both hosts

- **Sliver** (`render.cljd` + `sliver_collection.cljd`) — windowed; owns the
  segment/gen, the `AnimationController`, and the correction above.
- **Box** (`box.cljd`) — **the same rect animator, not a separate one.** It shares
  the whole core — `{key → rect}` committed/from maps, `gen`-freeze, `lerp-frame`,
  tight lerped constraints per cell, the resting/animating mode switch, enter/exit
  wrapping — and only *drops* the parts that exist because the sliver scrolls:
  - no windowing / virtualization (it measures every child anyway) ⇒ no
    `first/last-index`, no union window, no GC;
  - no scrolling ⇒ no `anchor-delta` / scroll correction; both endpoints are exact
    because the box measures everything.
  What it *adds* over the sliver: its **own** size lerps with the cells
  (`lerp(committed-content-extent, target-content-extent, t)`) so the box's height
  animates as its content reflows. `key-of` is the child index (positional) unless
  a child carries a `ValueKey`. Concretely, `box.cljd` loses the v1 `tween-layout`
  map wrapper and `:id`-change gating; `RenderCollectionBox` gains the committed/
  from/gen fields and a lerp branch in `content-size`; `CollectionBoxState` bumps
  `gen` + `forward` on any `didUpdateWidget` that could move a child.

## 11. Retired

- `MoveRender`, `MoveLayer`, `MoveRenderWidget`, `capture-positions!`,
  `visual-position` — subsumed by the engine's keyed rect glide.
- The `SizeTransition`+`FadeTransition` default insert/remove *size* driver —
  size is now the engine's; the fade/transform builder hook stays.
- `:layout` as an `:animate` sub-key and all `:id`-change detection
  (`maybe-layout-tween!`, the `layout-changed` branches).

## 12. Axis change is just another rect change

An `:axis` swap (`:vertical` ↔ `:horizontal`) is **not** a special case in v2 — it
falls out of the unified model. A cell's content-space rect (offset, cross,
main-extent, cross-extent) simply lands at a different target under the new axis,
and `lerp-frame` interpolates it component-wise like any other move. (This was a
v1 non-goal; it now works for free.)

## 13. Non-goals

Coexistence with `:reorderable` (kept mutually exclusive — `SliverReorderableList`
animates its own reorders).

## 14. Verification

The properties above are held by a `RenderViewport` simulation harness (oracles
O1–O12) plus a seeded invariant fuzzer over random op sequences. How to run them,
what each oracle asserts, and how to turn a fuzz failure into a red test:
**`docs/CollectionTesting.md`**.

The invariants that map directly onto this document: O6 (the anchor's consumed
fraction survives an op that must not move content — §9d), O11 (the truth
equation `off + frac·extent == scrollOffset` — §9c), O7 (a warm morph settles to
the same shape a cold start at the same anchor produces — §4.1/§7a), O4 (tiling,
rest-only because §8a lays cells at full size and paint-clips them), O9 (the
per-pass work bound of §9f).
