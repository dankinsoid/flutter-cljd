# Design: RenderViewport simulation harness + invariant fuzzer

**Task**: collection-arch-hardening, step 2 (design). Consumers: steps 3–6.
**Status**: design accepted, unimplemented.

This document is the sole input for the implementation steps. Everything an
implementer needs — widget tree, API shapes, field names, run commands, traps —
is fixed here.

---

## 0. Ground facts (verified, do not re-derive)

- `clojure -M:cljd test [ns...]` compiles the listed namespaces (all, if none)
  and then execs **`flutter test`** (`:cljd/opts :kind :flutter` ⇒ `bin` =
  `flutter`). Bare `flutter test` alone runs stale compiled Dart.
- `flutter test` runs **every** `test/**/*_test.dart`. Narrowing the namespace
  list on the CLI narrows *compilation only*, not the run. Test selection is
  `-t <tag>` / `-x <tag>` / `-N <substring>` after `--` (or `++` to append to
  `:dart-test-args`).
- A namespace is emitted into `test/cljd-out/**_test.dart` only if it
  `is-test-ns`: it has `:dart.test/dir` ns-meta **or** contains at least one
  `deftest` (which generates the `main` entry point). A helper ns with neither
  is compiled into the *library* output; a helper ns with the meta but no
  `deftest` produces a `_test.dart` **without `main`** and fails the run.
  ⇒ **`harness.cljd` must contain at least one `deftest`.** A smoke test.
- `deftest` accepts only `:tags` and `:runner`; there is **no per-test timeout
  option**. `dart test`'s default 30 s per-test cap applies. Raise it globally
  with `-- --timeout 4x`, or keep each `deftest` short.
- `:runner (ft/testWidgets [tester])` works (clojuredart `doc/TESTING.md`,
  `samples/widget_tests`); `flutter_test` is already a dev-dependency; `await`
  in the body compiles (the enclosing fn becomes `async`).
- `clojure.test` and `cljd.test` are the same namespace in ClojureDart; existing
  tests require `clojure.test`. Keep that for consistency.
- `flutter test` runs in debug ⇒ `fd/kDebugMode` is true ⇒ the engine's
  `assert-window-canonical!` and `assert-materialization-bounded!` fire. The
  harness inherits both for free.
- Engine deftype fields are plain public Dart fields, readable from any ns
  (`sliver_collection.cljd` already does `(.-approxSize! ro ...)`).
  **⇒ introspection needs zero changes to `render.cljd`.**

---

## 1. Attach route

### 1.1 Decision: public host, not `collection-sliver-widget`

The harness pumps `flutter-cljd.widgets/sliver-collection` (the public host)
inside `flutter-cljd.widgets/custom-scroll`. Rationale: the bug class under
investigation lives in the **host↔engine protocol** — `seg-gen` bumps, the
shadow/dying index space, `key-of`, `didUpdateWidget` diffing. Driving
`render/collection-sliver-widget` directly would bypass `SliverCollectionState`
and therefore bypass the entire segment/morph path, which is where the active
wrap→list bug is. A raw-sliver variant may be added later for micro-tests; it is
not part of this design.

### 1.2 The tree

```clojure
(defn- ^w/Widget rig-tree [rig]
  (let [{:keys [ctl items layout cross main animate? approx cold-gen probe]} @rig]
    (w/Directionality
     .textDirection ui/TextDirection.ltr
     .child
     (w/MediaQuery
      .data (w/MediaQueryData .size (ui/Size cross main))
      .child
      (w/Align                       ; loosens the surface constraints
       .alignment w/Alignment.topLeft
       .child
       (w/SizedBox
        .width cross .height main    ; the viewport box — cross change = rotation
        .child
        (w/KeyedSubtree
         .key (w/ValueKey cold-gen)  ; bump = fresh element tree (cold start)
         .child
         (f/custom-scroll
          :controller ctl
          :physics :clamping         ; deterministic ballistics; no overscroll
          (probe-sliver
           probe
           (apply f/sliver-collection
                  items
                  :key-fn :id
                  (concat
                   (when animate? [:animate {:duration segment-ms
                                             :move {:duration segment-ms}}])
                   (when approx [:approximate-item-size approx])
                   (layout-args layout)
                   [(fn [_ctx item] (tile item layout))]))))))))))
```

Notes:

- `w/Align` is required: `pumpWidget` hands the root **tight** surface
  constraints (800×600), so a bare `SizedBox` would be ignored.
- `:physics :clamping` is pinned. Bouncing physics adds an overscroll spring
  whose settle time is long and whose negative `pixels` muddies top-edge
  oracles. One dedicated scenario test may opt into `:bouncing`.
- Cross-extent change (`:cross`) is a **rebuild with a different `SizedBox`
  width**, not `tester.view.physicalSize`. It exercises the same
  `(not= crossAxisExtent crossCached)` invalidation branch in `performLayout`
  and avoids `view` reset/teardown pitfalls.
- `cold-gen` is bumped only by `cold-restart!` (mirrors `ws2-preview`).

### 1.3 Rebuild = re-`pumpWidget`

There is **no** stateful shim. Every mutation writes the rig atom and calls
`(await (.pumpWidget tester (rig-tree rig)))`. The root widget's runtimeType and
key are unchanged, so the element tree is reused and `SliverCollectionState
.didUpdateWidget` fires exactly as in the app.

Corollary the implementer must respect: `didUpdateWidget` starts a segment when
`(not (identical? new-data old-data))`. A no-op rebuild must therefore pass
**the same vector instance**; every mutating op must produce a **new** vector.
The rig holds `:items` as an immutable vec and never mutates in place.

### 1.4 Item builder — deterministic, text-free

Items are `{:id ^int n}`. Sizes derive from **`:id`, never from the index**, so a
shuffle actually re-shapes the layout and the cold-layout equivalence oracle is
not trivially satisfied.

```clojure
(defn ^double item-main  [^int id] (+ 44.0  (* 13.0 (mod (* id 37) 7))))   ; 44..122
(defn ^double item-cross [^int id] (+ 60.0  (* 17.0 (mod (* id 23) 5))))   ; 60..128

(defn ^w/Widget tile [item layout]
  (w/SizedBox .width  (item-cross (:id item))   ; ignored under a list's tight cross
              .height (item-main  (:id item))))
```

No text ⇒ no font metrics, no `--track-widget-creation` const-folding traps, no
platform-dependent layout. Sizes are non-uniform on purpose (defeats any
fixed-extent fast path) and are exact doubles.

`layout-args`:

| key       | args                                                                 |
|-----------|----------------------------------------------------------------------|
| `:list`   | `[:spacing 6]`                                                       |
| `:wrap`   | `[:layout :wrap :spacing 6 :run-spacing 6]`                          |
| `:grid`   | `[:layout :grid :columns {:count 3 :gap 8} :rows {:size {:aspect 1} :gap 8}]` |
| `:masonry`| `[:layout :masonry :cross-axis-count 2 :spacing 6 :run-spacing 6]`   |

The rig records the spacings so the tiling oracle knows the legal gaps.

### 1.5 Correction probe (optional, zero release footprint)

A test-only proxy sliver between the viewport and the collection records, per
layout pass, the correction the engine emitted. This is the only way to observe
**intra-frame** corrections: `RenderViewport` applies them via
`position.correctBy`, which does not notify listeners, so nothing is visible
from outside the frame.

```clojure
(deftype RenderProbeSliver [^:mutable sink]        ; sink = atom of vec
  :extends r/RenderProxySliver
  (performLayout [this]
    (.performLayout ^super this)
    (let [^r/RenderSliver ch (.-child this)
          ^r/SliverGeometry g (.-geometry ch)]
      (swap! sink conj
             {:scroll-offset (.-scrollOffset ^r/SliverConstraints (.-constraints this))
              :correction (.-scrollOffsetCorrection g)
              :scroll-extent (.-scrollExtent g)
              :paint-extent (.-paintExtent g)}))))

(deftype ProbeSliver [^w/Widget child sink]
  :extends (w/SingleChildRenderObjectWidget .child child)
  (createRenderObject [_ _ctx] (RenderProbeSliver sink))
  (updateRenderObject [_ _ctx ^RenderProbeSliver ro] (.-sink! ro sink)))
```

`RenderProxySliver` is abstract but its `performLayout` is concrete. If
`^super` dispatch on it misbehaves under cljd, inline the two lines it does:
`(.layout this child (.-constraints this) .parentUsesSize true)` then
`(.-geometry! this (.-geometry child))`.

The probe is transparent (pass-through constraints and geometry) and lives only
in `test/`. If it turns out to fight the child manager or the growth direction,
**drop it** — every oracle except O8 degrades gracefully to the pump-boundary
form. Do not add a debug counter to `render.cljd` to replace it without
re-opening this design.

### 1.6 Fallback route (documented, not chosen)

A hand-rolled `RenderSliverBoxChildManager` + a directly constructed
`RenderViewport` driven by explicit `layout()` calls. Rejected: it means
re-implementing `SliverMultiBoxAdaptorElement`'s child lifecycle (create /
remove / keep-alive / `findChildIndexCallback`), and it loses the host state
machine entirely — the very protocol under test. Reach for it only if
`testWidgets` proves structurally unusable, and re-open this design first.

---

## 2. Rig lifecycle

```clojure
(defn make-rig
  "Creates and mounts a rig. Returns the rig atom. Must be awaited."
  [^ft/WidgetTester tester
   & {:keys [items layout cross main animate? approx seed]
      :or {layout :list cross 400.0 main 600.0 animate? true}}]
  ...)                              ; pumps the first frame, registers teardown
```

- `ctl` = one `w/ScrollController` per rig (`.initialScrollOffset` for cold
  starts). Registered for disposal: `(ft/addTearDown (fn [] (.dispose ctl)))`.
- `cold-restart! [rig ^double init-offset]` swaps in a fresh controller, bumps
  `:cold-gen`, pumps. The retired controller is disposed on teardown, not
  inline — its old subtree is still attached during that frame.
- Default dataset: 2200 items, `(mapv (fn [i] {:id i}) (range 2200))` — matches
  `ws2-preview` so the manual repro and the harness repro are the same shape.
- Every test **must** end at rest (`settle!`). A live ticker at teardown makes
  `flutter_test` fail with "A Ticker was active when its provider was disposed".

**Async discipline (the #1 implementation trap):** any harness fn containing
`await` compiles to an `async` Dart fn returning a `Future`. Every call site
must `await` it. A forgotten `await` yields silent frame-ordering nondeterminism
that looks like an engine bug. Convention: every side-effecting harness fn is
named with a trailing `!` and **must** be awaited. Never wrap awaited calls in
lazy seq operations (`map`, `for`); use `doseq`/`loop`, which inline the `await`
into the enclosing async fn.

---

## 3. Scenario API

All fns take the rig atom first and are awaited. Sign convention: **positive
delta = scroll forward** (content moves up, `pixels` increases).

### 3.1 Time

| fn | meaning |
|----|---------|
| `(pump! rig)` / `(pump! rig ms)` | one frame / one frame after `ms` |
| `(pump-n! rig n ms)` | `n` frames of `ms` each (default 16) |
| `(settle! rig)` | `pumpAndSettle` with an explicit **10 s** timeout cap |
| `(advance-segment! rig)` | 5 pumps of `segment-ms/5`, invoking `check!` at each — inspects mid-morph frames |

`settle!` must pass the timeout explicitly; the flutter default is 10 minutes,
which turns a hang into a CI stall. When a scenario risks a non-settling state
(the bug being hunted), use `pump-n!` and assert instead.

### 3.2 Input

| fn | implementation |
|----|----------------|
| `(scroll-by! rig dy)` | `(.dragFrom tester center (ui/Offset 0.0 (- dy)))` then `pump!` |
| `(fling! rig speed)` | `(.flingFrom tester center (ui/Offset 0.0 (- (sign speed) 300.0)) (abs speed))` then bounded pumps |
| `(grab! rig)` / `(move! rig dy)` / `(release! rig)` | raw `TestGesture` — for interrupt scenarios (catch a fling mid-ballistic) |
| `(jump-to! rig off)` | `(.jumpTo ctl off)` + `pump!` |
| `(animate-to! rig off ms curve)` | call `.animateTo`, **discard the Future**, then pump across `ms` |

`dragFrom`/`flingFrom` take a start `Offset` — no `Finder`, no `Type` literal,
which sidesteps cljd's awkward type-value interop. `center` is the rig's
viewport center.

`animate-to!` must **not** await `.animateTo`'s Future: it completes only when
the animation ends, which needs pumps, which the await would prevent. Deadlock.

### 3.3 Data mutation (each rebuilds + pumps)

| fn | note |
|----|------|
| `(insert! rig idx id)` | explicit `id` — replay determinism |
| `(remove-at! rig idx)` | |
| `(swap-items! rig i j)` | keyed shuffle without changing the count |
| `(rotate! rig i j)` | block move — the realistic reorder |
| `(set-items! rig v)` | wholesale replace |

### 3.4 Structural

| fn | note |
|----|------|
| `(set-layout! rig kw)` | `:list` ↔ `:wrap` ↔ `:grid` ↔ `:masonry` mid-scroll; drives the layout `:id` change ⇒ `crossLayoutReanchor` |
| `(set-cross! rig w)` | simulated rotation; drives the `crossCached` cache drop |
| `(set-approx! rig v)` | toggle `:approximate-item-size` |
| `(cold-restart! rig off)` | fresh element tree + controller at `off` |

---

## 4. Introspection

Read-only. No engine field is ever written from the harness.

```clojure
(defn ^render/CollectionRenderSliver sliver [rig]
  (first (filter #(instance? render/CollectionRenderSliver %)
                 (.-allRenderObjects ^ft/WidgetTester (:tester @rig)))))
```

`allRenderObjects` (a `WidgetController` getter) avoids `tester.renderObject<T>`,
whose Dart type argument is painful from cljd.

| fn | returns |
|----|---------|
| `(offset rig)` | `(.-pixels (.-position ctl))` |
| `(scroll-limits rig)` | `{:min :max :viewport}` from the `ScrollPosition` |
| `(geometry rig)` | `{:scroll-extent :paint-extent :cache-extent :max-paint-extent :layout-extent :overflow?}` from `SliverGeometry` |
| `(children rig)` | vec, attached order, of `{:index :offset :cross :main-size :cross-size :end}` — `layoutOffset` / `crossAxisOffset` off `SliverGridParentData`, sizes off `(.-size child)` (vertical axis ⇒ main=height, cross=width) |
| `(engine rig)` | the internals map below |
| `(probe-log rig)` / `(probe-reset! rig)` | the correction log |
| `(snapshot rig)` | `{:offset :scroll-extent :children}` — the canonical comparable |

```clojure
(defn engine [rig]
  (let [rs (sliver rig)]
    {:cache-first (.-cacheFirst rs)      :cache-n (count (.-cache rs))
     :anchored0 (.-anchoredTo0 rs)       :regime (.-curRegime rs)
     :reanchor-shift (.-reanchorShift rs) :cross-reanchor (.-crossLayoutReanchor rs)
     :landing-emitted (.-landingEmitted rs) :pending-landing-ws (.-pendingLandingWs rs)
     :cur-gen (.-curGen rs)              :seg-gen (.-segGen rs)
     :seg-anchor (.-segAnchor rs)        :resting-top (.-restingTop rs)
     :tweening? (some? (.-tweenAnim rs)) :vel (.-velocityEMA rs)
     :pass-materialized (.-passMaterialized rs)
     :checkpoints (sort (keys (.-checkpoints rs)))
     :committed-n (count (.-committed rs))
     :overscan [(.-overscanBehind rs) (.-overscanAhead rs)]}))
```

This map is what a failing test prints. It is also the diff surface for the
later rebase/anchor-primary refactors: a field that disappears from it is a
subsystem that was actually removed.

**Tolerances**: `eps-exact 1e-9` (idle stability — identical inputs must give
identical doubles), `eps-geom 1e-3` (tiling/spacing arithmetic),
`eps-equiv 0.5` (cold-vs-warm equivalence, where estimates legitimately differ
sub-pixel-ish).

---

## 5. Invariant oracles

`(check! rig ctx)` returns a vec of violation maps `{:oracle :detail}`; tests
assert it is empty and print `ctx` + `(engine rig)` + `(children rig)` on
failure. The fuzzer calls it after every op. Oracles that are meaningful only at
rest are skipped when `(:tweening? (engine rig))` or the position is still
ballistic; the driver records *why* it skipped so a scenario cannot silently
check nothing.

**O1 — no exception.** `(.takeException tester)` is nil. This is the carrier for
the engine's own debug asserts (`assert-window-canonical!`,
`assert-materialization-bounded!`) *and* for Flutter's
"RenderViewport exceeded its maximum number of layout cycles" — i.e. a
non-convergent correction loop is caught for free. Must be checked after **every**
op; `flutter_test` defers exceptions, so an unchecked one surfaces in an
unrelated later test.

**O2 — top edge.** When `pixels <= minScrollExtent + eps-geom`: the first
attached child has `index == 0` and `offset == 0` (±eps-geom), and no attached
child has a negative `offset`. This is the direct oracle for the active bug
("top element #1700, nothing above it, clamps at 0").

**O3 — idle stability.** Take `snapshot`, `(pump! rig)`, take `snapshot`; they
must be equal within `eps-exact` (offsets, sizes, `pixels`, `scrollExtent`).
Rest-only. Catches correction loops that oscillate below the viewport-cycle cap
and any per-frame drift.

**O4 — tiling.** Rest-only (during a segment §8a lays cells at full size and
paint-clips them, so layout rects legitimately overlap).
- *Generic, all layouts*: no two attached children's rects overlap on both axes
  by more than `eps-geom`; sorting the main-axis intervals and merging them
  leaves no hole larger than `run-spacing + eps-geom`; the merged span covers
  `[max(scrollOffset, first-child.offset), min(scrollOffset+paintExtent, last-child.end)]`.
  The coverage clause is what catches "nothing above item 1700".
- *`:list`*: `offset[i+1] == offset[i] + main[i] + spacing` for consecutive
  attached indices; every `cross == 0`.
- *`:wrap` / `:masonry`*: group children by equal `offset` (±eps-geom); within a
  run, `cross` strictly ascends with `cross[i+1] == cross[i] + cross-size[i] +
  spacing` and the run's total cross ≤ `crossAxisExtent + eps-geom`; the next
  run starts at `offset + max(main in run) + run-spacing`.
- *`:grid`*: cell rects match the grid spec exactly (indexed layout ⇒ pure math).

**O5 — extent sanity.** `scrollExtent` finite (finite item count) and
`>= last-attached-child.end`; `maxPaintExtent >= paintExtent`;
`minScrollExtent - eps <= pixels <= maxScrollExtent + eps` under clamping
physics.

**O6 — anchor preservation.** For ops that must not move content — layout swap,
cross change is *excluded*, insert/remove **above** the window, `:approx`
toggle, segment settle — the item key at the viewport top before the op is the
same key after, and its intra-item offset (`pixels - offset(top item)`) is
preserved within `eps-equiv`. This is the ws2-preview #6 scenario ("no offset
settle") and the wrap→list morph criterion, expressed as an invariant.

**O7 — cold-layout equivalence.** After a morph (or any op sequence) settles at
offset `X`: `cold-restart!` a rig with `initialScrollOffset X`, same items,
layout and cross; settle; compare.
- If **both** rigs report `anchoredTo0` true ⇒ compare absolute: same top index,
  same offsets for the common indices (±eps-equiv).
- Otherwise ⇒ compare **shape**: for the indices attached in both, the pairwise
  offset *differences*, the sizes and the cross positions must match
  (±eps-equiv), and the top index must match. Absolute agreement is not
  required — a cold deep start is by design a checkpoint-relative estimate.
  Also assert the top index is *identical*: an estimate may shift the absolute
  offset, it may not land on a different item.

This oracle is the acceptance gate for step 13 (anchor-primary): the refactor is
correct exactly when warm and cold produce the same shape at the same anchor.

**O8 — correction convergence** (needs the probe). Per frame: at most 10 layout
passes on the collection; the **last** pass of every frame emits
`scrollOffsetCorrection == nil`; while idle (no input driven by the harness, no
active segment) no frame emits a correction at all. Without the probe, O1+O3
cover the divergent subset only.

**O9 — bounded state.** `cache-n`, `committed-n` and `(count checkpoints)` stay
below `k * (count (children rig)) + 64` (k ≈ 8). Catches unbounded retention
across long episodes — the thing a single-frame assert cannot see.
`pass-materialized` is already covered by the engine's own tripwire (O1).

---

## 6. Fuzzer

### 6.1 PRNG

Hand-rolled 32-bit LCG in cljd (`next = (1664525*s + 1013904223) mod 2^32`), not
`dart:math/Random` — reproducibility must not depend on the Dart SDK version.
The seed is an explicit `int`, printed in every failure message.

### 6.2 Op alphabet

Ops are **fully parameterized values** — no PRNG draws happen during replay.
That is what makes shrinking and red-test pasting work.

| op | form | weight |
|----|------|--------|
| small scroll | `[:drag ^double dy]`, \|dy\| ∈ 0.3–1.5 viewport, either sign | 20 |
| idle | `[:pump ^int n]`, n ∈ 1–5 | 12 |
| fling | `[:fling ^double speed]`, \|speed\| ∈ 800–8000 | 10 |
| large scroll | `[:drag ^double dy]`, \|dy\| ∈ 2–8 viewports | 8 |
| insert | `[:insert ^int idx ^int id]` (idx anywhere, incl. above the window) | 8 |
| remove | `[:remove ^int idx]` | 8 |
| settle | `[:settle]` | 8 |
| layout swap | `[:layout :list\|:wrap\|:grid\|:masonry]` | 6 |
| far jump | `[:jump ^double off]`, off ∈ 0–50 viewports, incl. 0 and max | 6 |
| shuffle | `[:rotate ^int i ^int j]` | 4 |
| cross change | `[:cross ^double w]` from `#{280.0 400.0 640.0 900.0}` | 4 |
| to top | `[:to-top]` (jump 0) / `[:animate-top ^int ms]` | 4 |
| interrupt | `[:fling-catch ^double speed ^int ms]` — fling, pump `ms`, grab | 3 |

The episode header records the rig config the seed drew: `{:seed :animate? :approx :layout0 :cross0 :items0}`.

### 6.3 Episode

- Fixed length: **24 ops**.
- After each op: bounded pumps (op-specific), then `(check! rig ctx)` and
  `(.takeException tester)`.
- First violation ⇒ fail the test immediately with:

```
FUZZ FAIL seed=<n> op=<k>/24 oracle=<id>
header: {:seed 1234 :animate? true :approx nil :layout0 :wrap :cross0 400.0 :items0 2200}
ops: [[:jump 30240.0] [:layout :list] [:pump 3] ...]
engine: {...}   children: [...]
```

The `ops` vector is valid cljd data. `(replay! rig header ops)` consumes exactly
that shape, so turning a fuzz failure into a red test is copy-paste. `replay!`
is part of the harness, not the fuzzer, so red tests never depend on
`fuzz_test.cljd`.

### 6.4 Shrinking (documented; implementation optional)

`(shrink header ops pred)` where `pred` re-runs a **fresh** rig over an op list
and returns true iff it still fails the same oracle:

1. **Prefix minimization** — binary-search the shortest failing prefix.
2. **Greedy removal** — walk ops left to right, drop each, keep the drop if
   `pred` still fails. One pass is enough in practice; repeat until a fixed
   point if cheap.

Both are sound only because ops carry concrete parameters. Index-bearing ops
(`:insert` / `:remove` / `:rotate`) may become out-of-range after a removal —
`replay!` clamps indices into range and this is *defined behaviour*, not a
crash. Shrinking is a convenience; step 6 may ship without it.

---

## 7. File layout, run commands, performance

```
test/flutter_cljd/internal/collection/harness.cljd       ; rig, ops, introspection, oracles, replay!, smoke deftest
test/flutter_cljd/internal/collection/harness_test.cljd  ; scenario tests (step 4) + reds (step 5)
test/flutter_cljd/internal/collection/fuzz_test.cljd     ; fuzz episodes (step 6)
```

- `harness.cljd` ⇒ ns `flutter-cljd.internal.collection.harness`. It **must**
  contain `(deftest harness-smoke :tags [:widget] :runner (ft/testWidgets [tester]) ...)`
  — see §0. The smoke test mounts a rig, pumps, asserts children are attached
  and `check!` is empty.
- Tags: scenario tests `:tags [:widget]`; fuzz tests `:tags [:widget :fuzz]`.
- Fuzz is split across **several `deftest`s** (e.g. 6 tests × 3 seeds each)
  because `deftest` cannot set a per-test timeout and `dart test` caps at 30 s.

Run:

```sh
clojure -M:cljd test                                   # everything (compiles all)
clojure -M:cljd test -- -x fuzz                        # CI-fast: no fuzz
clojure -M:cljd test -- -t fuzz --timeout 4x           # fuzz only, relaxed cap
clojure -M:cljd test -- -N far-scroll-wrap-to-list     # one test by name
```

Optional deps.edn alias (implementer's call; combine as `-M:fast:cljd test`):

```clojure
:fast {:cljd/opts {:dart-test-args ["-x" "fuzz"]}}
```

Performance target: **baseline (`-x fuzz`) under ~60 s of `flutter test` time**,
excluding the ClojureDart compile step (which dominates the wall clock and is
unavoidable). Budget: a pump over a ~30-child window is sub-millisecond-ish;
~30 scenario tests × ~60 pumps ≈ 2 000 pumps ≈ a few seconds. The fuzz suite at
18 episodes × 24 ops × ~5 pumps ≈ 2 200 pumps lands in the same order. The real
risks to the target are `settle!` calls that burn their full timeout and the
2 200-item dataset's O(n) host-side diff on every rebuild — if the latter shows
up, drop the fuzz dataset to ~600 items and keep 2 200 only for the far-scroll
scenarios.

---

## 8. Risks and mitigations

| # | Risk | Mitigation |
|---|------|------------|
| 1 | `pumpAndSettle` hangs on a non-settling state (exactly the bug class hunted) | Always pass an explicit 10 s timeout; prefer `pump-n!` in scenarios expected to be unstable; a timeout is a legitimate *failure signal*, assert on it |
| 2 | Live ticker at teardown ⇒ "Ticker was active when its provider was disposed" | Every test ends with `settle!`; controller disposal via `ft/addTearDown` |
| 3 | Forgotten `await` on an async harness fn ⇒ nondeterministic frame order that mimics an engine bug | `!`-suffix convention; never await inside lazy seq ops; the smoke test asserts a known post-`scroll-by!` offset, which fails loudly if the await discipline breaks |
| 4 | `RenderProxySliver` subclassing misbehaves under cljd | Inline the two-line `performLayout`; or drop the probe entirely (only O8 degrades) |
| 5 | Engine debug asserts abort layout mid-frame, corrupting the tree for later ops | On the first violation the fuzzer aborts the episode and returns; never continue after O1 fires |
| 6 | `dart test`'s 30 s per-test cap kills long episodes | Short episodes, many `deftest`s, `--timeout 4x` for the fuzz tag |
| 7 | Namespace narrowing does not narrow the run — a broken unrelated test noise-floors the loop | Use `-N`/`-t`/`-x`; document in the test ns docstring |
| 8 | ~240 baseline `DYNAMIC WARNING:` compile lines drown a real error | Ignore compiler warnings entirely (per CLAUDE.md); gate on `bin/check`'s `dart analyze` error severity |
| 9 | `(= w w)` false for const-folded widgets under `--track-widget-creation` | The harness never compares widget identity; oracles compare geometry only |
| 10 | All `_test.dart` files import each other — one compile error breaks the whole suite | Expected; `bin/check` first when a compile fails |
| 11 | Cold-equivalence (O7) flakes because a cold deep start legitimately estimates a different absolute offset | O7 compares shape + top index, not absolute offsets, unless both rigs are `anchoredTo0` |
| 12 | Default surface is 800×600 and `pumpWidget` constraints are tight | `w/Align` before the `SizedBox` (§1.2) — omit it and the rotation op silently does nothing |

---

## 9. Immediate consumer — the step-5 red

The active bug as a harness scenario, to be written first and expected to fail:

```
rig: 2200 items, :layout :wrap, :cross 400.0, :animate? true
1. (jump-to! rig (* 50.0 viewport))     ; ~50 viewports deep
2. (settle! rig)                        ; record top key K and intra-offset
3. (set-layout! rig :list)              ; wrap → list morph mid-scroll
4. (advance-segment! rig)               ; oracles at each mid-morph frame
5. (settle! rig)
oracles: O6 (top key still K), O7 (equals a cold list at the same anchor),
         O4 (list tiling — no hole above the first child),
         then (jump-to! rig 0.0) + (settle! rig) ⇒ O2 (index 0 at offset 0)
```

Expected failure signature per the current hypothesis (plan §"Active bug
chain"): the capture pass ends on a correction-only geometry, so the snapshot's
`:extent` is truncated and the rebase is dropped — visible as O6 (top key jumps)
and O2 (top edge clamps with a non-zero first index).
