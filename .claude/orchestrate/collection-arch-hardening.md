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
- [x] 4. Green baseline scenario tests (scroll/jump/rotate/morph basics) — agent: general-purpose, model: sonnet
- [x] 5. Red repro: far-scroll wrap→list morph + correction-only capture extent — agent: general-purpose, model: opus
- [x] 6. Invariant fuzzer over random op sequences — agent: general-purpose, model: opus
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
