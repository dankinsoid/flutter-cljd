# Task: Unified keyed engine-owned per-cell rect animator (collection v2)

**Slug**: collection-rect-animator
**Started**: 2026-07-13
**Status**: Phase 1 DONE (bin/check --clean 0 errors; 186 tests green; landed).
Phase 2 IN PROGRESS — runtime bugs found on-device, §7a redesign (see below).

## Goal
Replace v1's three animation paths (MoveRender position-glide, tween-layout
size-morph, SizeTransition insert/remove) with ONE keyed, engine-owned per-cell
rect animator. The engine animates each visible cell's content-space rect
(offset+size) from its last-committed rect toward its live target, whatever the
cause. No `:id` detection; `:move` is the single public knob; `:layout` removed;
the diff shrinks to enter/exit only. Both hosts (sliver + box). Design approved.

Success = both demos morph cleanly (incl. shuffle + simultaneous data+layout),
stable-size cells do zero child relayouts during a move, anchored child stays put
while content above animates, `bin/check` clean, `clojure -M:cljd test` green.

## Raw request
> Refactor the collection animation system into ONE keyed, engine-owned per-cell
> rect animator, replacing the current two separate paths (MoveRender + tween
> size-morph). v2 on top of v1 (636cf22). Full spec: see the session prompt and
> docs/CollectionRectAnimator.md (the approved design note).

## Design authority
`docs/CollectionRectAnimator.md` is the source of truth (approved, with edits:
key-fn optional / index fallback first-class; axis-morph now supported, not a
non-goal; box shares the core minus windowing+correction). Every subagent reads it.

## Code map
- `src/flutter_cljd/internal/collection/tween.cljd` — v1 combinator + NEW pure
  `anchor-delta` (done). Home for host-agnostic keyed machinery.
- `src/flutter_cljd/internal/collection/render.cljd` — sliver engine
  (CollectionRenderSliver, flow-layout!/indexed-layout!, set-tween-anim!,
  snapshot-window, update-render!).
- `src/flutter_cljd/internal/collection/box.cljd` — box host (RenderCollectionBox,
  content-size, CollectionBox/State).
- `src/flutter_cljd/internal/sliver_collection.cljd` — host state, didUpdateWidget,
  maybe-layout-tween!, apply-diff*, build.
- `src/flutter_cljd/internal/collection/animation.cljd` — MoveRender + SizeTransition
  + parse-anim-config (to retire/trim).
- `src/flutter_cljd/internal/collection/diff.cljd` — compute-diff (enter/exit).
- `example/src/repl.cljd` — Morph (~748) + Box Morph (~791) demos.
- Tests: `test/flutter_cljd/internal/collection/{tween,box}_test.cljd`.
- Verify: `bin/check` (compile+analyze), `clojure -M:cljd test`.

## Decisions log
- 2026-07-13: mode switch kept (resting=natural measure/self-heal; animating=tight
  lerped, no re-measure; natural size measured once per segment). (by: coordinator, user-approved)
- 2026-07-13: union window = target-window ∪ attached-children-still-overlapping
  (lives partly in driver; committed is key-addressed, no cheap index inverse). (user-approved)
- 2026-07-13: shuffle cost bounded — only keys with live `from` glide; pruned/
  never-seen appear at target. (user-approved)
- 2026-07-13: key-fn OPTIONAL; index-as-key is first-class. (by: user)
- 2026-07-13: axis change is just another rect change; supported, not a non-goal. (by: user)
- 2026-07-13: box = same animator minus windowing/correction, plus own-size lerp;
  key-of = child index or ValueKey. (user-approved)
- 2026-07-13: enter/exit tick on the SAME global segment clock (an insert/remove is
  just a rect main-extent 0↔final change) — no per-cell controllers. (by: user)
- 2026-07-13: collapse wrap = clip only for now — ClipRect > OverflowBox(final,
  align-start); content laid at final size, engine clips. :scale deferred. Custom
  :insert/:remove :builder preserved (wraps content inside the OverflowBox), size
  owned by the engine (SizeTransition size driver retired). (by: user)
- 2026-07-13: enter/exit realized in the combinator via `entering`/`exiting` key
  sets → from/to becomes a 0-main-extent version of the target frame (grow/shrink),
  so it rides the same keyed-tween. (by: coordinator)
- 2026-07-13 PIVOT (design flaw found by Step-4 subagent): a host-built `to-src`
  reads the OLD render tree (new/inserted cells aren't built until the next build),
  so flow insert/remove/shuffle MIS-KEYS. FIX: the ENGINE captures `to-src` from the
  NEW tree on the segment's first pass; the host only triggers. Reuses Step-2
  committed/union/correction. Design note §4.1/§6/§7 updated. (by: user-approved)
- 2026-07-13: gen-change pass mechanism — for a FLOW target the new inserted cells
  only materialize when a real driver walks them, so the engine runs the normal flow
  driver once to create+measure them, snapshots that as `to-src` (index-frozen over
  the NEW indices), then re-lays attached children at `from` (t≈0) so there's no
  target-flash. INDEXED target: `to-src` = indexed-frame-source (pure math, no
  pre-walk). Subsequent passes run the indexed keyed-tween driver. Fallback if the
  flow lifecycle entanglement is too risky: LIVE-TARGET (measure natural + tight-lerp
  each pass) — correct but defers constraint-1 zero-relayout. (by: coordinator)

## Open questions
- [ ] Enter/exit content-collapse wrap: confirm ClipRect+OverflowBox(final,align-start)
      default + :scale option + custom :builder interaction (raise at that step).

## Checklist
- [x] 0. Design note + pure `anchor-delta` core + test — done (commits landed)
- [x] 1. Keyed geometry combinator in tween.cljd (pure `keyed-tween-layout`) +
        unit tests — done (+175 green, committed)
- [x] 2. Sliver engine: committed map + key-of plumbing + union windowing +
        anchor-delta correction (gated on tweenAnim) — done (bin/check clean, 179 green)
- [x] 3. Enter/exit primitives: keyed-tween-layout 7-arg {:entering :exiting} +
        collapse-wrap in animation.cljd — done (+182 green, bin/check clean, committed)
- [x] 4E. Engine segment (render.cljd): gen-change capture from new tree + from-relay
        display + subsequent keyed-tween passes — done (bin/check clean, +182, committed).
        HOST CONTRACT: `update-render! (rs target-layout item-count tween-anim ^int gen
        entering exiting key-of)`; widget opts `:seg-gen`, `:entering`, `:exiting`,
        `:key-of`. Host passes PLAIN target layout (never a tween) + tween-anim + gen++
        per segment + entering/exiting key-sets + key-of. curGen field = incoming gen.
- [x] 4H. Sliver host: trigger-only + diff→enter/exit + collapse-wrap + retire :id/
        MoveRender wiring — done (bin/check clean, 182 green, committed)
- [x] 5. Retire dead code: MoveRender/MoveLayer/capture-positions/visual-position from
        animation.cljd; SizeTransition default size driver — agent: general-purpose
- [x] 6. Box host: RenderCollectionBox committed/from/gen + content-size animating
        branch + own-size lerp; CollectionBoxState gen+forward; retire tween wrapper
        — agent: general-purpose
- [x] 7. parse-anim-config: drop :layout, :move single knob; widgets.cljd docstrings
        — agent: general-purpose
- [x] 8. Tests: content-space continuity, stable-size zero-relayout, key-match
        shuffle, enter/exit vs scroll-off, retarget continuity, box own-size
        — agent: general-purpose
- [x] 9. Demos: repl.cljd Morph + Box Morph exercise shuffle + simultaneous
        data+layout — agent: general-purpose
- [x] 10. Final verify: bin/check --clean + full test run — coordinator

## Phase 2 — §7a redesign (runtime fixes, 2026-07-13)

On-device testing of the Morph demo surfaced two bugs with ONE root cause, plus the
design fix. Design authority: `docs/CollectionRectAnimator.md` §7a (canonical).

**Bugs (both = leaving cells still occupy the indexed target's index space):**
- Indexed REMOVE animates nothing. `build-shadow-data` reinserts the dying cell at
  its old index → index space frozen → `indexed-frame-source` gives `to-frame(i) ==
  committed(i)` → stay-glide dead (nothing slides up); collapse shrinks in a reserved
  slot that never vacates → "vanish + snap". (diag: render.cljd:863-864 +
  animation.cljd:108-116; tween.cljd:291-292 else-branch is committed→committed.)
- SHUFFLE bounces. Segment materializes only the target window; from-visible cells
  (scattered to new indices) never laid at committed → viewport-top empty → bounce.

**Fix (§7a):** indexed `to-src` built from PURE NEW DATA (no dying reinsert). Cells
leaving the on-screen window → OVERLAYS at committed pos, outside index space →
survivors reindex → stay-glide revives. Two exit/enter semantics: data-collapse
(0↔final in place) vs viewport-slide (full-size slide, no squeeze, viewport clips).
leave-slide + exit-collapse = one overlay path, two drivers. arrive-slide = ordinary
contiguous indexed child laid at offscreen committed → target. Span-cap → snap if
from+target visible sets exceed a viewport.

**Retires:** `build-shadow-data` reinsert-at-old-index / reserved-slot exit;
`widen-window!` union (for indexed).

### Phase 2 checklist
- [ ] 2a. Design note §7a — DONE (docs/CollectionRectAnimator.md).
- [ ] 2b. Engine: indexed `to-src` from pure new data (drop dying reinsert in the
        target index space); survivors reindex so stay-glide/remove/insert glide.
        Files: render.cljd (segment-start!/indexed path), animation.cljd
        (build-shadow-data). Verify remove: below-items slide up.
- [ ] 2c. Overlay path: leaving cells (exit-collapse + leave-slide) rendered as
        overlays at committed rect, outside index space; two drivers (collapse=0↔final,
        slide=full-size→edge+clip). Host (sliver_collection.cljd dyingItems) +
        engine paint. Classifier per §7a matrix.
- [ ] 2d. Span-cap → snap when visible-from ∪ visible-target exceeds a viewport.
- [ ] 2e. Verify on-device: remove (below slides up + collapse), shuffle (no bounce,
        first shuffle animates), insert, list↔grid morph anchor stable. bin/check
        clean + tests green + tween/render unit tests for reindex-on-remove.

## Step results

### 2. Sliver engine support
- Status: done (committed)
- Files: render.cljd, render_test.cljd (new)
- Added: CollectionRenderSliver fields `keyOf`, `committed` ({key->rect+:seen}),
  `passTick`. Widget/opts key `:key-of`. Fns: key-for, child-rect, prune-committed,
  commit-pass!, `committed-rects` (host getter), widen-window-indices, widen-window!,
  indexed-anchor-cells. performLayout increments passTick + commit-pass! (always-on).
  Union windowing + anchor-delta correction gated on tweenAnim (resting unchanged).
- HOST CONTRACT for Step 3:
  * Engine expects opts `:key-of (fn [^int i] -> key)` in collection-sliver-widget.
  * Host reads `render/committed-rects` at segment start → the `from` map for
    keyed-tween-layout. Snapshot at fresh `forward(from:0)`, NOT mid-correction
    (commit-pass! runs even on correction/early-return passes; harmless but snapshot
    at segment start).
  * Host builds keyed-tween-layout(id, to-src, committed, from-extent, key-of, progress)
    and passes it as the layout + the Animation as :tween-anim, exactly like v1 passes
    tween-layout today.

### 1. keyed-tween-layout (pure combinator)
- Status: done
- Files: tween.cljd (+keyed-tween-layout, arity 6), tween_test.cljd (+5 deftests)
- Signature: `(keyed-tween-layout id to-src committed ^double from-extent key-of progress)`.
  :geometry lerps committed[(key-of i)] → to-frame(i); missing key → target (no glide);
  windowing passes through to-src; :max-extent lerps from-extent→to-extent.
- Next-step impact: Step 2 driver supplies committed map + from-extent + key-of, and
  owns the committed-side union (attached children still overlapping window).

### 0. Design note + anchor-delta
- Status: done
- Files: docs/CollectionRectAnimator.md, tween.cljd (anchor-delta), tween_test.cljd
- Key findings: no viewport/widget-tester harness in tests (183 pure unit tests) →
  scroll correction tested as pure fn. anchor-delta green.
- Commits: design note; anchor-delta+test; design revisions.
