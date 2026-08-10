# Testing the collection engine

How the virtualization engine (`internal/collection/render.cljd`) is verified:
a `RenderViewport` simulation harness with an invariant oracle battery, a seeded
fuzzer over random op sequences, and a set of deterministic reproductions.

The design the tests hold is `docs/CollectionRectAnimator.md`.

## 1. Running

```sh
clojure -M:cljd test -- -x "known-red || fuzz"   # the default green suite
clojure -M:cljd test -- -t known-red             # the documented deferrals
clojure -M:cljd test -- -t fuzz --timeout 20x    # the campaign
clojure -M:cljd test -- --plain-name harness-smoke
bin/check                                        # compile example + dart analyze
```

Mechanics that bite:

- **Never bare `flutter test`** — it runs whatever Dart was compiled last.
  `clojure -M:cljd test` recompiles, then execs it.
- Namespace arguments after `test` narrow **compilation only**, not the run. A
  narrowed compile also leaves the other namespaces' `_test.dart` stale, so a run
  right after a rename fails on unrelated files. Compile everything before
  trusting a green run.
- Selection is `-t <tag>` / `-x <tag>`. **Two `-x` flags do not compose** — the
  last one wins. Pass one boolean selector: `-x "known-red || fuzz"`.
- By-name selection is `--plain-name <substring>` (or `--name <regexp>`).
  `:cljd/opts :kind :flutter` execs `flutter test`, which **rejects `-N`** — the
  flag `dart test` documents.
- `--timeout 20x` is not optional for the fuzz tag: one batch is 10 episodes ×
  ~25 ops × a full oracle pass, and every failing episode additionally pays for
  shrinking. `dart test` caps a test at 30 s.
- `deftest` has no per-test timeout knob; that is why the fuzz campaign is split
  across many small `deftest`s rather than one.

Tags: `:widget` on everything that mounts a rig, `:fuzz` on campaign batches,
`:known-red` on a reproduced defect that is deliberately not fixed yet.
`-t known-red` failing is the expected state; `-t known-red` passing means a red
was closed and its tag should come off.

## 2. The harness

`test/flutter_cljd/internal/collection/harness.cljd` mounts the **public** host
(`custom-scroll` + `sliver-collection`) under `testWidgets`, so the host↔engine
protocol — segment gens, the shadow/dying index space, `didUpdateWidget` diffing
— is exercised for real. Everything it does towards `render.cljd` is read-only;
no engine field is ever written.

- Every mutation writes the rig atom and re-`pumpWidget`s. There is no stateful
  shim: the root widget's type and key are unchanged, so the element tree is
  reused and `didUpdateWidget` fires exactly as in the app.
- Items are `{:id n}` and their sizes derive from the **id, never the index**, so
  a shuffle really re-shapes the layout. No text ⇒ no font metrics, no
  `--track-widget-creation` const-folding traps.
- A test-only `RenderProxySliver` probe between the viewport and the collection
  logs the correction of every layout pass. This is the only way to see an
  *intra-frame* correction: `RenderViewport` applies them via `position.correctBy`,
  which notifies nobody.
- `flutter test` runs in debug, so the engine's own asserts
  (`assert-window-canonical!`, `assert-materialization-bounded!`, the
  exactly-once rebase assert) fire and arrive through O1.

**Async discipline — the number one trap.** Every `!`-suffixed harness fn compiles
to an `async` Dart fn and MUST be awaited. Never await inside a lazy seq op
(`map`/`for`); use `doseq`/`loop`, which inline the await into the enclosing async
fn. A forgotten await produces frame-ordering nondeterminism that reads exactly
like an engine bug.

Every test must end at rest (`settle!`) or `flutter_test` fails the teardown with
"A Ticker was active when its provider was disposed".

## 3. The oracles

`(check! rig ctx)` returns a vec of `{:oracle :detail}` violations;
`check-report!` also returns `:ran` and `:skipped`, so a scenario cannot silently
check nothing. Rest-only oracles are skipped while a segment runs or the position
is still ballistic, and the skip is recorded.

O1 runs first and, if it fires, nothing else runs — the tree is mid-layout after
an engine assert.

| # | asserts |
|---|---|
| **O1** | no deferred exception. The carrier for the engine's debug asserts AND for Flutter's "exceeded its maximum number of layout cycles" — a non-convergent correction loop is caught for free. Must run after EVERY op: `flutter_test` defers exceptions, so an unchecked one surfaces in an unrelated later test. Reading it CLEARS it. |
| **O2** | at the top of the range, index 0 sits at offset 0 and nothing is above it. The direct oracle for "top element #1700, nothing above, clamps at 0". |
| **O3** | idle stability: pump one frame at rest and nothing moves (`eps-exact` 1e-9). Catches correction loops oscillating below the ten-cycle cap, and per-frame drift. |
| **O4** | tiling, REST-ONLY (§8a lays cells at full size and paint-clips them, so mid-segment rects legitimately overlap): no two rects overlap on both axes, no hole larger than the run spacing, the merged span covers the viewport; plus per-layout rules for list / wrap / grid / masonry. |
| **O5** | extent sanity: `scrollExtent` finite and ≥ the content laid within the reachable horizon (`pixels + 2 viewports` — a far-moving segment pre-materializes its END window past the lerping extent by design); `maxPaintExtent ≥ paintExtent`; `pixels` inside the range (out-of-range is legal mid-spring under bouncing physics). |
| **O6** | anchor preservation for ops that must not move content (layout swap, insert/remove strictly above the window, `:approx` toggle, segment settle): the content point at the before-key's **consumed fraction** is still at the viewport top. The invariant is the fraction, not the absolute intra offset — a morph that resizes the anchor cannot hold both, and holding the pixel offset is what pushes the top past a shrinking anchor. Checked on the before-key's OWN after-child, never on the reported top child: in a multi-column target the top is a set of row-mates sharing one span. Cross-extent changes are excluded by design. |
| **O7** | cold-vs-warm equivalence: settle, `cold-restart!` the same rig at the current offset, compare. Absolute offsets only when BOTH rigs report `anchoredTo0`; otherwise shape (pairwise offset differences, sizes, cross positions) plus an identical top index — a cold deep start is by design a checkpoint-relative estimate, but it may not land on a different item. **Destructive**: the rig runs on a fresh element tree afterwards. |
| **O8** | correction convergence (needs the probe): at most 10 correcting passes in a row, the last pass of a frame never leaves a correction pending, and no frame emits a correction while nothing drives the viewport. |
| **O9** | two disjoint bounds. **WORK**: `pass-materialized` stays under the engine's own `materialization-budget` over the band the pass was given — O(window+overscan) per pass, unconditional, no mode exempt (×2 mid-segment, where a capture legitimately lays a widened window). **RETENTION**: `cache-n` / `committed-n` / `checkpoints` stay proportional to the live window; mid-segment the denominator is the BAND, because a landing collapses the attached window to one child while the segment legitimately retains band-scale state. |
| **O10** | the anchor never teleports at rest: across an idle pump `viewAnchor`'s (key, frac) is bit-stable. Bites when something schedules relayout at rest. |
| **O11** | the truth equation: after a settled pass `off + frac·extent == scrollOffset` (rest-only; mid-segment the displayed frames tween against a moving set-point). |
| **O12** | extent quiescence: at rest with no mutations `maxScrollExtent` is constant across pumps. An extent that breathes at rest restarts boundary springs (`applyNewDimensions → goBallistic`) and betrays EMA drift leaking into the reported extent. |

Tolerances: `eps-exact` 1e-9 (idle stability — identical inputs must give
identical doubles), `eps-geom` 1e-3 (tiling arithmetic), `eps-equiv` 0.5
(cold-vs-warm, where estimates legitimately differ sub-pixel-ish).

O6 cannot run from `check!` alone — it needs a `before` captured at rest. The
fuzzer captures it itself around layout swaps and strictly-above-window
inserts/removes.

## 4. The fuzzer

`fuzz_test.cljd`. An episode = a seed → a rig header + 24 fully parameterized
ops, over 600 items. The PRNG is a hand-rolled 32-bit LCG, so a seed reproduces
an episode independently of the Dart SDK's `Random`, and **no PRNG draw ever
happens while ops are applied** — which is what makes shrinking work and makes
the printed `ops` vector paste-able into a test.

After every op: bounded pumps, then the whole oracle battery. An episode stops at
its first violation; a BATCH keeps going and reports every failing seed at once
(`is` throws in this cljd port, so one `is` at the end is the only way a deftest
can report more than one failure).

### Profiles and seed ranges

| profile | seeds | what it starves, and why |
|---|---|---|
| `:full` | 1–60 | the design's weight table verbatim |
| `:no-layout` | 101–160 | no layout swap — a morph bug that fires early masks the rest of the op space |
| `:no-jump` | 201–260 | no far jump / to-top, same reason |
| `:churn` | 301–340 | data mutation. `:layout` starved and `:animate?` forced, so a count change is the ONLY segment opener and `:settle` is rare enough that the next one lands mid-segment. Adds `:remove-anchor`, which deletes the very item the key re-anchor looks for. |
| `:bounce` | 401–410 | `:full`'s op space under `BouncingScrollPhysics` — real overscroll, so boundary springs run |

A full campaign is 230 episodes and ~60–90 s of `flutter test` time. Seed ranges
are ratchets: they are never renumbered, so a campaign is comparable against any
earlier one.

### Reading a failure

```
FUZZ FAIL seed=235 profile=:no-jump op=12/24 kind=:o6 oracles=[:o6]
header: {:animate? true :approx nil :layout0 :grid :cross0 400.0 :items0 600}
ops: [[:drag -420.38] [:cross 280.0] ... [:layout :wrap]]
engine: {...}   children: [...]
```

The `ops` vector is valid cljd data. Failures are shrunk automatically
(`shrink!`: prefix cut — free, since the failing index is known — then greedy
left-to-right removal, each candidate on a fresh rig), so what gets printed is
usually 1–6 ops.

## 5. Turning a fuzz failure into a red test

Paste the shrunk vector into `fuzz_red_test.cljd`:

```clojure
(deftest known-red-<what-must-hold>
  :tags [:widget :known-red]
  :runner (ft/testWidgets [tester])
  (let [f (await (run-seed! tester 235
                            [[:drag -420.3800868988037]
                             [:cross 280.0]
                             [:layout :wrap]]))]
    (is (nil? f) (str "<the invariant, stated positively>" (fail-str f)))))
```

- `run-seed!` **regenerates the header from the seed** (`profile-of` maps the seed
  to its profile) rather than transcribing it, so `:animate?` / `:approx` /
  `:layout0` / `:cross0` can never drift out of sync with the generator.
- It replays through `fuzz-test/run-ops!` — the same driver that found and shrank
  the vector, including the bounded post-jump boundary wait O5 needs and the
  anchor capture O6 needs. A red that replays through a hand-rolled loop instead
  reproduces something slightly different.
- `harness/replay!` is the lower-level form (`rig-from-header!` + an op vector,
  checking after each op) for a test that needs its own header.

Discipline that has paid for itself repeatedly:

1. **Verify standalone before calling it an engine bug.** Re-run the candidate as
   the FIRST episode of its own deftest. `flutter_test` defers exceptions, so a
   failure attributed to episode N is sometimes episode N−1's leak.
2. **Reproduce three times.** A red that only reproduces once is a rig-pollution
   artifact — sequential rigs under one tester share tickers and controllers.
3. When a red goes green, drop the `:known-red` tag rather than deleting the
   test: it becomes a ratchet.
4. Fixing an oracle counts as a fix only when the new failure set is a strict
   SUBSET of the old one. Anything else is masking.

## 6. What is not covered here

- Pure kernels (`render_test.cljd`, `tween_test.cljd`, `layout_test.cljd`) test
  the engine's pure functions directly, with no viewport. Fastest loop; prefer
  them whenever a defect can be stated over inputs and outputs.
- `box_test.cljd` drives `RenderCollectionBox` detached via manual `.add`. It does
  NOT transfer to the sliver, which needs a real `RenderSliverBoxChildManager`.
- `bin/check` compiles the example app and runs `dart analyze` on the generated
  Dart. `clojure -M:cljd compile` alone never runs the Dart front-end, so an
  internally-inconsistent build reports success and fails only at `flutter run`.
  What `dart analyze` still cannot see: calls to members ClojureDart could not
  resolve — emitted as `(x as dynamic).m()`, failing only at runtime. Dynamic
  dispatch is idiomatic here, so the compiler's `DYNAMIC WARNING:` lines are
  baseline noise, not a gate.
- No test asserts what the engine *paints*. O4 checks layout rects; the §8a
  paint-clip is verified by eye.
