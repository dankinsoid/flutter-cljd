# Task: Unified keyed engine-owned per-cell rect animator (collection v2)

**Slug**: collection-rect-animator
**Started**: 2026-07-13
**Status**: in-progress

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
- 2026-07-13: enter/exit content wrap (ClipRect+OverflowBox vs :scale) DEFERRED —
  discuss with user when Step reached.

## Open questions
- [ ] Enter/exit content-collapse wrap: confirm ClipRect+OverflowBox(final,align-start)
      default + :scale option + custom :builder interaction (raise at that step).

## Checklist
- [x] 0. Design note + pure `anchor-delta` core + test — done (commits landed)
- [ ] 1. Keyed geometry combinator in tween.cljd (pure `keyed-tween-layout` +
        union-window helper) + unit tests — agent: general-purpose
- [ ] 2. Sliver engine: committed/from/gen keyed state + key-of plumbing + indexed/
        animating driver applies keyed lerp + emits anchor-delta correction + union
        windowing — agent: general-purpose
- [ ] 3. Sliver host (sliver_collection.cljd): gen-bump+forward on any move-causing
        update; build key-of (data/key-fn/shadow); retire maybe-layout-tween!/:id
        detection — agent: general-purpose
- [ ] 4. Enter/exit: ClipRect+OverflowBox wrap; retire SizeTransition default size;
        keep :builder hook (USER CHECK-IN first) — agent: general-purpose
- [ ] 5. Retire MoveRender/MoveLayer/capture-positions from animation.cljd + callers
        — agent: general-purpose
- [ ] 6. Box host: RenderCollectionBox committed/from/gen + content-size animating
        branch + own-size lerp; CollectionBoxState gen+forward; retire tween wrapper
        — agent: general-purpose
- [ ] 7. parse-anim-config: drop :layout, :move single knob; widgets.cljd docstrings
        — agent: general-purpose
- [ ] 8. Tests: content-space continuity, stable-size zero-relayout, key-match
        shuffle, enter/exit vs scroll-off, retarget continuity, box own-size
        — agent: general-purpose
- [ ] 9. Demos: repl.cljd Morph + Box Morph exercise shuffle + simultaneous
        data+layout — agent: general-purpose
- [ ] 10. Final verify: bin/check --clean + full test run — coordinator

## Step results

### 0. Design note + anchor-delta
- Status: done
- Files: docs/CollectionRectAnimator.md, tween.cljd (anchor-delta), tween_test.cljd
- Key findings: no viewport/widget-tester harness in tests (183 pure unit tests) →
  scroll correction tested as pure fn. anchor-delta green.
- Commits: design note; anchor-delta+test; design revisions.
