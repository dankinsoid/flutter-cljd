# Collection Rect Animator (v2)

> Design note. Supersedes the v1 two-path scheme (`MoveRender` position-glide +
> `tween-layout` size-morph, committed at 636cf22). Written before coding —
> architecture and abstractions first, implementation second.

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
from      : {key → rect}   ;; frozen segment start; committed at (re)start time
target    : {index → rect} ;; this pass's goal, from the live layout (transient)
gen       : int            ;; host bumps on each (re)start; engine compares
```

`rect = {:offset :cross :main-extent :cross-extent}` in content space.

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

## 8. Enter / exit (constraint 4)

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

## 9. Scroll correction (the gotcha, tested first)

The **indexed** layout path emits geometry directly and does **not** call
`scrollOffsetCorrection`. When content above the viewport collapses — a remove, or
list→grid shrinking the region before the anchor — `scrollOffset` can jump and the
anchored child visibly shifts. The flow path already handles this (`anchor-before`
/ `anchor-delta` → correction); the indexed/animating path must grow the same:

Before committing geometry, snapshot the anchored child (first visible cell at the
current `scrollOffset`) and its offset; after computing this pass's rects, if that
key's offset moved by > ε, emit `scrollOffsetCorrection = delta` so the viewport
follows it. A failing test (`anchored child stays put while content above it
animates`) is written **before** the engine change.

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

## 14. Test plan

1. **Scroll-anchor stability** (write first, failing): anchored child fixed while
   content above animates its rect; engine emits the correction.
2. Content-space rect continuity across a scroll (no spurious animation).
3. Stable-size cell ⇒ zero child `performLayout` during a move (assert constraint
   identity / a layout-count probe).
4. Key-based match across a data shuffle (a key's `from` follows the key, not the
   index).
5. Enter/exit classification vs. scroll-off (windowed absence ≠ data absence).
6. Retarget continuity (interrupt mid-glide, re-aim, no jump — keyed).
7. Box: own-size lerp; exact endpoints.
