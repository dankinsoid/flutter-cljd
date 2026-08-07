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

## Open questions
- (none)

## Checklist
- [x] 1. Explore test infra + example app usage of collection — agent: Explore, model: sonnet
- [x] 2. Design viewport harness API — agent: Plan, model: opus
- [x] 3. Implement harness core — agent: general-purpose, model: opus
- [ ] 4. Green baseline scenario tests (scroll/jump/rotate/morph basics) — agent: general-purpose, model: sonnet
- [ ] 5. Red repro: far-scroll wrap→list morph + correction-only capture extent — agent: general-purpose, model: opus
- [ ] 6. Invariant fuzzer over random op sequences — agent: general-purpose, model: opus
- [ ] 7. Triage fuzz findings → additional red tests — agent: general-purpose, model: sonnet
- [ ] 8. Design pending-rebase structure (exactly-once protocol) — agent: Plan, model: fable
- [ ] 9. Implement rebase unification (mechanism-by-mechanism, harness green between) — agent: general-purpose, model: opus
- [ ] 10. Verify step-5 reds green; revert TEMP fb24f35 — agent: general-purpose, model: sonnet
- [ ] 11. Capture-mode: design + implement (corrections statically off, rebase returned as value) — agent: general-purpose, model: opus
- [ ] 12. Design anchor-primary migration (truth = anchor idx+intra offset; equivalence criteria; kill-list of subsystems) — agent: Plan, model: fable
- [ ] 13. Implement anchor-primary in stages behind harness equivalence — agent: general-purpose, model: fable
- [ ] 14. Cleanup: remove dead subsystems, update docs/CollectionRectAnimator.md + memory — agent: general-purpose, model: sonnet

## Step results

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
