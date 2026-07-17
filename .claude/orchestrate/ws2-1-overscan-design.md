# WS2.1 — Overscan window (design)

Design pass for the overscan window that WS2.1 opened in `collection-scroll-anim-followups.md`
§2.1, sitting under the two already-recorded decisions (overscan-scoped committed cache; the
`:renewal-point?`/`:approx-offset` protocol hooks from §2.2) and feeding the progressive-offset
refinement in `ws2-3-offset-refinement-research.md`. Docs only — no `src/` changes here.

File:line pointers are as of this reading (`src/flutter_cljd/internal/collection/…`, `…/sliver_collection.cljd`,
`layouts.cljd`, `widgets.cljd`); WS1.1/1.3 landed after the follow-up plan was written, so names were
re-located rather than trusting the plan's line numbers. Re-verify before editing.

---

## 0. What overscan is, in one paragraph

The engine already materialises the window `[ws, we]` where `ws = scrollOffset + cacheOrigin` and
`we = ws + remainingCacheExtent` — i.e. Flutter's viewport **cacheExtent** window, handed in through
`SliverConstraints` (`flow-layout!` render.cljd:807-809, `indexed-layout!` render.cljd:921-922).
**Overscan** is an engine-internal, **directional, velocity-scaled** widening of that window —
`[ws − behind, we + ahead]` — used for **child materialisation, GC bounds, and real-geometry
capture only**, never for what the sliver advertises back to the viewport as painted/cache extent.
Its job (per §2.1's title "layout slightly beyond the viewport") is to keep entering cells and
off-window neighbours **real, not estimated**, so the rect animator's `from`/`to` and WS2.3's
run-boundary cache are seeded from measured geometry — and to read ahead of fast motion so cells are
laid out before they scroll in.

---

## 1. Relationship to Flutter's `cacheExtent` — compose, don't duplicate

| | Flutter `cacheExtent` | our overscan |
|---|---|---|
| owner | viewport (`RenderViewport.cacheExtent`, default 250px, or `:cache-extent` on the scroll view / `custom-scroll-view`) | this RenderSliver, per pass |
| symmetry | **symmetric** (both edges equally) | **asymmetric** (more ahead of motion) |
| velocity | **independent** | **scaled** by the current scroll velocity |
| purpose | keep children *warm*: painted-when-scrolled-in, semantics, keep-alive | keep off-window neighbours **real** for animation + WS2.3; read ahead of fling |
| surfaced in geometry | yes — `SliverGeometry.cacheExtent`, shared across sibling slivers | **no** — engine holds the extra children privately (it drives its own `collectGarbage`) |

The two **stack additively**: the effective materialisation window is
`[ws − behind, we + ahead]`, where `[ws, we]` *already* contains Flutter's cache region. At rest with
default config `ahead ≈ behind ≈ small base`, so we barely exceed the cache window — **no duplicated
work when idle**. We never set the viewport's `cacheExtent` from overscan: that knob is symmetric and
would also change painting, semantics, and every sibling sliver's `remainingCacheExtent`. A user who
wants a larger *symmetric warm* region keeps using `:cache-extent`; overscan is strictly the *dynamic
directional* layer on top.

Cross-framework grounding for that split (the §2 table in the WS2.3 research is the reconciliation
half; this is the read-ahead half):

- **RecyclerView** — `LinearLayoutManager.calculateExtraLayoutSpace(state, extra[])` fills a
  2-element `[start, end]` array and by default puts the extra layout space **in the scroll
  direction** (asymmetric), and a **full extra page** while a `LinearSmoothScroller` is running.
  `GapWorker` prefetches adjacent positions using the **last (dx, dy)** — directional, and bounded by
  the **frame deadline** (it only prefetches what fits before the next vsync). → asymmetric + bounded,
  no unbounded velocity scaling.
- **UICollectionView** — `isPrefetchingEnabled` (default on); UIKit picks cells to prefetch **ahead**
  of the scroll from velocity+direction. Read-ahead, asymmetric, magnitude-aware.
- **react-virtuoso** — `overscan: number | {main, reverse}` and `increaseViewportBy: {top, bottom}`,
  both in **pixels**, explicitly split **main (ahead)** vs **reverse (behind)**. → the closest prior
  art to what we want; px-based, asymmetric, but *not* velocity-scaled (fixed).
- **TanStack Virtual** — `overscan: number` = **item count**, **symmetric**, default 1. → simplest
  model; rejected as our currency (item-count is meaningless for variable-size cells and multi-column
  wrap — N items = wildly different px, and our whole engine speaks px/offsets).
- **Flutter cacheExtent** — symmetric px (or `CacheExtentStyle.viewport` multiplier), velocity-independent;
  `calculateCacheOffset`/`calculatePaintOffset` already used by `finish-flow-geometry!` (render.cljd:784-785).

Takeaway used below: **px is the internal currency**; **ahead≠behind**; **velocity scales ahead with
a cap** (nobody scales unbounded — RV bounds by frame budget, virtuoso by a fixed number).

---

## 2. Architecture — parts, responsibilities, boundaries

Five parts. The first three are new; the last two are existing engine points that overscan re-scopes.

### 2.1 Velocity/direction sensor  *(new, render-object-local)*
**Responsibility:** produce a smoothed `(velocity, direction)` estimate each frame, entirely inside
the RenderSliver, with no host plumbing.

**Why render-object-local, not host-pushed:** `SliverConstraints` carries `userScrollDirection`
(forward / reverse / **idle**) but **no velocity magnitude** (confirmed: sliver.dart fields are
`scrollOffset`, `cacheOrigin`, `remainingCacheExtent`, `viewportMainAxisExtent`, `precedingScrollExtent`,
`userScrollDirection` — no velocity). `ScrollPosition.activity.velocity` exists but `activity` is
`@protected`; reaching it means walking to the `Scrollable` and is fragile. The robust, self-contained
source is the **per-frame delta of `scrollOffset`** which the sliver sees every `performLayout`:

```
velocity ≈ (scrollOffset_now − scrollOffset_prevFrame − Σ corrections_this_frame) / dt
```

- `dt` from `SchedulerBinding.instance.currentFrameTimeStamp`; **update once per frame** (multiple
  layout passes share a frame — a `scrollOffsetCorrection` forces a re-layout in the same frame, so
  gate on the timestamp changing).
- **Correction-compensated:** the engine's own `scrollOffsetCorrection` (render.cljd:833, 836, 1047)
  moves `scrollOffset` without real motion; subtract the last emitted correction so a WS1.3/WS2.3
  offset reconciliation isn't read as a fling. This is the one subtle invariant.
- **Direction** = sign of the smoothed velocity; corroborated by `userScrollDirection` for the idle
  case (drag released, not yet ballistic).
- **Smoothing:** an EMA over velocity plus **hysteresis on the direction flip**, so a brief drag
  reversal or the natural zero-crossing at a fling's end doesn't thrash the asymmetric window every
  frame (see §5 direction reversal).

Captures user drag, ballistic fling, and programmatic `animateTo` **uniformly** — the engine reacts
to *motion*, not to its *source* (this is what makes the animateTo answer in §4 fall out for free).

New mutable fields on `CollectionRenderSliver`: `^:mutable ^double lastScrollOffset`,
`^:mutable ^double velocityEMA`, `^:mutable ^double lastCorrection`, `^:mutable frameStamp`, plus a
smoothed `^:mutable ^int scrollSign` for hysteresis.

### 2.2 Overscan resolver + window  *(new, pure)*
**Responsibility:** turn the user's `:overscan` config + the live `(velocity, direction, env)` into a
concrete `[behind ahead]` pixel pair. Two pure fns, mirroring how `parse-anim-config` (animation.cljd:28)
normalises `:animate`:

- `resolve-overscan` — parse the config form once (at widget build, like `parse-anim-config`) into a
  normalised `{:base :ahead-max :behind-max :lookahead-ms :viewport-frac? :fn}` record.
- `overscan-window` — pure `(config, velocity, direction, viewport-main-extent) → [behind ahead]`.

**Velocity law** (dimensionally clean, matches GapWorker's "what arrives before the next frames"):

```
lead   = clamp(base_ahead,  |velocity| · lookahead-ms/1000,  ahead-max)   ; px of motion to stay ahead of
trail  = base_behind                                                       ; Flutter cache already covers behind
ahead  = (direction forward?)  lead  : trail
behind = (direction forward?)  trail : lead                                ; asymmetry follows the sign
```

`lookahead-ms` ≈ 100–150ms ("stay one-or-two frames of motion ahead"); `ahead-max` ≈ 1× viewport
(cap — beyond a screen of read-ahead the cells fly past unseen, pure waste; this is the "cap" the
frameworks all impose one way or another). At rest `|velocity|=0` ⇒ `ahead=behind=base`, a small
symmetric margin so immediate neighbours are real for edge inserts.

### 2.3 Fidelity gate  *(new, threshold over the sensor)*
**Responsibility:** classify the current pass into a **regime** and pick materialisation fidelity.
This is the single mechanism that answers the animateTo question (§4) and the fling question (§5) —
both reduce to "how fast, and how close to rest."

| regime | trigger | overscan | window fill | WS2.3 refine |
|---|---|---|---|---|
| **resting / slow** | `|v|` below `v_slow` | base (symmetric) + velocity lead | **real** build+layout+place | on (every dirty pass) |
| **flying** | `|v|` above `v_fast` | **disabled** (capped to 0 ahead) | **approximate** placement via `:approx-offset` (WS2.2); only the visible window, no backfill | **suppressed** |
| **landing** | `|v|` decaying below `v_slow` while near a stable target | base + shrinking lead | **real** + WS2.3 upward re-flow to land exactly | on |

`v_slow`/`v_fast` are a hysteresis band (not one threshold) so a fling settling through the boundary
doesn't oscillate. "Near a stable target" for landing = velocity decaying **and** `scrollOffset`
delta shrinking frame-over-frame (the tail of any decelerating curve — fling or animateTo alike).

### 2.4 Committed-cache scoping  *(re-scope of an existing part — decision already recorded §2.1)*
Today `committed` (the persistent `{key → rect}` map, render.cljd:58) is pruned by **age**:
`prune-committed` (render.cljd:223-230) drops entries older than `committed-max-age` = 600 passes
(render.cljd:204-207), called from `commit-pass!` (render.cljd:258).

**Change (per the recorded decision):** prune **geometrically** — drop a key when its committed rect
left the **overscan** window `[ws − behind, we + ahead]`, not by age. The 600-pass constant and the
age branch go away.

- **Integration point:** `commit-pass!` (render.cljd:232-258) already runs once per pass after the
  driver; it gains the overscan window + the exemption set and calls a geometric `prune-committed`.
- **Mid-segment exemption (invariant):** never prune an active-segment participant —
  `entering ∪ exiting ∪ (keys of) leaving` (render.cljd:66-69) — while `tweenAnim` is non-nil. Prune
  only between segments. This is the persistent-rect analogue of what `widen-window!`
  (render.cljd:284-308) already does at the *child* level (keeping attached cells that still overlap
  under their committed rect). The two must agree: a key exempt from child GC must be exempt from
  committed prune, or its glide loses its `from`.
- **Why it's correct** (recorded): scoping to the window kills the stale-`from` class of bug — after
  WS1.1, re-inserting a long-gone key would otherwise animate from a stale offset. A key with no
  committed entry falls back to enter-wipe / `augment-from-edges` (render.cljd:1241) edge-slide, which
  is platform-standard (no glide from off-screen). Trade-off accepted: no long-range move glides —
  those drew stale geometry anyway.

### 2.5 Window application in the drivers  *(re-scope of existing points)*
The drivers compute `ws`/`we` then feed them to GC + fill. Overscan inserts **one widening step**
between "compute cache window" and "use it," and **only for materialisation/GC** — geometry reported
to the viewport stays on the true `[ws, we]`:

- `flow-layout!` (render.cljd:798-838): after `we = ws + remainingCacheExtent` (line 809), compute
  `[behind ahead]` and derive `ws_o = ws − behind`, `we_o = we + ahead`. Pass `ws_o/we_o` to
  `pre-gc!` (810), `walk!`'s fill target `we` (821), `post-gc!` (831), `backfill-leading!`'s `ws`
  (832), `align-start!` (818). **But** `finish-flow-geometry!` (765-796) keeps computing
  `paintExtent`/`cacheExtent` from constraints (Flutter cache) — overscan children are materialised,
  not advertised.
- `indexed-layout!` (render.cljd:912-1066): same shape — widen `ws/we` (921-922) before the
  index-range GC (938-942) and the `first-index/last-index` window (924-935); geometry (1050-1060)
  stays on the true window.

---

## 3. The decided model (justification)

**Units & shape — decided:**
- **Internal currency: pixels.** The engine speaks offsets/px everywhere (`[ws,we]`, cache-extent,
  `calculateCacheOffset`); item-counts (TanStack) are meaningless across variable-size / multi-column
  layouts. Config *accepts* a viewport-fraction form for ergonomics but resolves to px against
  `viewportMainAxisExtent` (available on `SliverConstraints`).
- **Separate leading/trailing — decided yes**, but expressed **relative to motion** (`ahead`/`behind`)
  and mapped to leading/trailing by the velocity sign each frame. Asymmetry *is* the feature (§1).
- **Velocity law — decided:** `lead = clamp(base, |v|·lookahead, max)`, time-based (`px/s · s = px`),
  capped at ~1 viewport. Rest ⇒ small symmetric base (neighbours real); motion ⇒ grows ahead, capped.
  Rejected: unbounded scaling (wasteful; no framework does it), fixed px (virtuoso — misses the
  "read ahead of a fling" requirement in user pt 1), item-counts (TanStack — wrong currency here).

**What overscan materialises — decided: full build + layout + place** within the overscan margin (not
measure-only). §2.1's whole point is that entering cells and neighbours are *real*: measure-only would
give geometry but no attachable/animatable child and no painted reveal for an edge insert. Concretely
three tiers of fidelity by distance from the viewport:
1. viewport + Flutter cache → painted / warm (existing behaviour, unchanged).
2. overscan margin (ours, directional) → **full real** build+layout+place, clipped-not-painted; this
   is what feeds WS2.3's run-boundary cache with **real run geometry** above/below (a wrap run entered
   here is a true multi-column row, not the single-left-column estimate of today's list-shaped
   `backfill-leading!`, render.cljd:703-741). Kills #5/#6 at the neighbour scale.
3. beyond overscan → **approximate offset only** (WS2.2 `:approx-offset`), no materialisation — used
   in the flying regime (§4).

**How it feeds WS2.3:** WS2.3 says "refinement runs on every dirty pass within the overscan window."
That window *is* this part's output. Real materialisation in tier 2 is exactly what produces the real
run boundaries WS2.3 caches; overscan sets refinement's scope. The coupling is deliberate — 2.1 owns
the *window*, 2.3 owns the *reconciliation within it*.

---

## 4. `animateTo` — the verdict (answering the user's open question)

**Question restated:** short programmatic scroll → precompute the whole path? very far → per-frame but
overscan off? or is per-frame always right?

**Verdict: per-frame layout is always the cadence. Never precompute a path. What varies with the
motion is *fidelity + overscan* (the §2.3 gate), not the cadence — and that variation is driven by
velocity, so a far `animateTo` and a fast fling are handled by the *same* mechanism with no need to
detect "programmatic vs ballistic" inside the engine.**

Reasoning per option:

- **Precompute the whole path for a short scroll — rejected.** A short `animateTo` (target within
  ~1–2 viewports) sweeps through offsets whose windows **heavily overlap**; the engine's per-index
  geometry cache (render.cljd:52) already memoises placement, so each intermediate frame is mostly
  cache hits + 0–1 new cells. Precomputing the union front-loads all that build work onto frame 1 (a
  jank spike right when the animation starts) to save work that per-frame already amortises to near
  zero. Net negative.
- **Far scroll, per-frame + overscan off — essentially right, with one addition.** A far `animateTo`
  travels continuously through *every* intermediate offset; at peak curve velocity the window can jump
  multiple viewports per frame → zero overlap, a fresh window every frame, and (worse) a huge index
  gap for backfill/walk. Materialising those fly-by frames for real is both **expensive and pointless**
  (cells flash past unseen). So in the flying regime: **disable overscan** (no point reading ahead of
  something nothing will see) **and place the visible window *approximately*** via WS2.2 `:approx-offset`
  rather than real backfill. As velocity decays near the target (landing regime), switch real
  materialisation + overscan + WS2.3 refinement back on so it **lands exactly** on the target (and, at
  offset≈0, exactly on element 0 — WS2.3's acceptance criterion). The user's "per-frame but overscan
  disabled" is right; the addition is *approximate placement*, not full real layout, for the fly-by.
- **Precompute a far path — doubly rejected.** Materialising thousands of intermediate cells up front
  is strictly worse than the per-frame-approximate fly-by.

**Why no programmatic-vs-fling detection is needed:** both produce a decelerating velocity profile;
the §2.3 gate keys off *velocity magnitude and deceleration*, which is identical for both. Collapsing
them is the simplification — the engine never asks "who started this scroll."

**Optional future (host-level, flagged): far-target warp.** UICollectionView `scrollToItem(animated:)`
and RecyclerView `smoothScrollToPosition` don't truly animate through a very far span — they **jump
near the target, then settle** (a short animated tail). That's a *ScrollController* concern, not an
engine one: the host owning the controller could rewrite a far `animateTo(target)` as
`jumpTo(near-target) + animateTo(target)`. The engine's flying/landing regimes already make the
straight-through far `animateTo` cheap and correct, so warp is a nicety, not required. Left to a
possible host-side helper; **not** in this engine.

### Trade-off table

| strategy | frame-1 cost | per-frame cost | memory | visual on land | verdict |
|---|---|---|---|---|---|
| precompute whole path (short) | **high burst** (build union) | ~0 | union of windows | exact | **rejected** — burst jank; per-frame already cheap via geometry cache |
| per-frame real, overscan on (short) | low | low (cache hits + 0–1 cells) | 1 window + overscan | exact | **chosen for short** |
| per-frame real, overscan on (far) | low | **high** (fresh non-overlapping window/frame + big backfill gap) | 1 window | exact | rejected for far — wasted real layout on unseen frames |
| per-frame **approx**, overscan **off** (far fly-by) → real+refine on landing | low | **low** (approx placement, no backfill) | 1 window | exact (refine on landing) | **chosen for far** |
| precompute whole path (far) | catastrophic | ~0 | **whole span** | exact | rejected — worst memory + burst |

---

## 5. Ballistic fling & direction reversal

- **Ahead grows with fling velocity — yes, capped.** Same `lead = clamp(base, |v|·lookahead, max)`
  law (§2.2); a fast fling widens the ahead margin so cells are laid out before they scroll in, up to
  the ~1-viewport cap. Above the `v_fast` threshold the fidelity gate flips to *flying* (approximate,
  overscan off) — so overscan doesn't grow without bound; it grows over the mid band and then *drops*
  once the fling is fast enough that read-ahead is pointless. (Net shape: overscan is a **hump** over
  velocity — zero-ish at rest, grows through the slow band, capped, then cut in the flying band.)
- **`scrollOffsetCorrection` does not perturb the fling** (WS2.3 §3.4, restated): it is layout-phase
  and doesn't re-target the active `Simulation`; the correction shifts content-space under a
  velocity-driven pixel trajectory, so the fling stays smooth. Do **not** re-aim the simulation.
- **Suppress above-anchor corrections while flinging *up* into estimated content** (WS2.3 §3.2 tier 3
  / TanStack's rule): a same-frame offset jump racing the simulation reads as a stutter. Defer to the
  first low-velocity/resting pass. Overscan cooperates: in the flying regime refinement is off anyway
  (§2.3), so the suppression is automatic there; it only matters in the slow band during an upward
  drag, where the fidelity gate keeps refinement on but WS2.3's directional suppression gates it.
- **Direction reversal — handled by EMA + hysteresis (§2.1).** On a sign flip the previously-*ahead*
  margin becomes *behind*: those cells are **already materialised**, so the instant of reversal is
  cheap (the window you hold already covers the new near side); the new ahead side then grows. The
  hazard is *thrash* — rapid tiny reversals reallocating the asymmetric window every frame. Guard with
  a smoothed velocity (EMA) and a **direction-flip hysteresis band** so a brief zero-crossing (natural
  at a fling's tail, or a jittery drag) doesn't swap the asymmetry until the reversed velocity clears a
  small threshold. `userScrollDirection` (idle/forward/reverse from constraints) corroborates the idle
  case so a released-but-not-yet-ballistic gap reads as "rest," not a spurious tiny velocity.

---

## 6. Config API sketch

Consistent with the house options style (kv pairs / options map, resolved once at build like
`parse-anim-config`). One new option `:overscan` on `sliver-collection` / `collection-view` /
`grid-view` / `collection`, threaded exactly like `:animate`:

```
sliver-collection*  (parse)                          ; sliver_collection.cljd:357-391
  → SliverCollection deftype field  :overscanConfig   ; new field, next to animConfig
  → build → collection-sliver-widget opts  {:overscan …}
  → CollectionSliverWidget field                       ; render.cljd:1324-1348
  → CollectionRenderSliver field  overscanConfig       ; via createRenderObject / update-render!
```

**Accepted forms** (resolver picks them apart, defaults chosen so zero-config is good):

```clojure
:overscan 300                 ; number → fixed px, symmetric, no velocity scaling
:overscan {:viewport 0.5}     ; fraction of viewport main extent → px at layout time
:overscan {:ahead 400 :behind 100 :max 800 :lookahead-ms 120}
                              ; full control; :max caps the velocity-scaled ahead
:overscan false               ; or 0 / :none → disable (Flutter :cache-extent still applies)
:overscan (fn [ctx] …)        ; ctx = {:velocity :direction :viewport :item-count :cache-extent}
                              ;   → returns px | {:ahead :behind} | {:viewport frac}
```

**Default (no `:overscan` given) — must be good with zero config:**
- `base` small (≈ one item or ~0.25 viewport) symmetric → immediate neighbours real for edge
  inserts, at-rest cost ≈ Flutter's cacheExtent only (no duplication).
- `lookahead-ms` ≈ 120, `ahead-max` ≈ 1.0 viewport → grows ahead of motion, capped.
- `behind` stays at base (Flutter's cache covers the just-passed region).

Rationale for not over-designing: a single `:overscan` accepting number | map | fn covers fixed,
fractional, and fully-custom without a widget zoo, and the number form is the 90% case. The fn form is
the escape hatch (velocity-aware custom laws) without committing protocol surface to it.

---

## 7. Integration points (by fn name)

| concern | fn / field | file:line | change |
|---|---|---|---|
| velocity sensor | new `update-velocity!` + fields `lastScrollOffset`/`velocityEMA`/`lastCorrection`/`frameStamp`/`scrollSign` | render.cljd `performLayout` 156-177 (top) | compute per-frame, correction-compensated |
| resolve config | new `resolve-overscan` (pure) | mirror `parse-anim-config` animation.cljd:28 | parse `:overscan` once at build |
| compute window | new `overscan-window` (pure) → `[behind ahead]` | called in both drivers | velocity law §2.2 |
| fidelity gate | new `fidelity-regime` (pure) → `:resting`/`:flying`/`:landing` | called in both drivers | §2.3 thresholds |
| flow widen | `flow-layout!` | render.cljd:798-838 | widen `ws/we` → `pre-gc!`(810) `walk!`(821) `post-gc!`(831) `backfill-leading!`(832) `align-start!`(818); geometry via `finish-flow-geometry!`(765) stays on true window |
| indexed widen | `indexed-layout!` | render.cljd:912-935 | widen before GC (938) + `first-index`/`last-index` window (924-935); geometry (1050) stays true |
| flying/approx fill | `backfill-leading!` / `walk!` | render.cljd:703-741 / 623-677 | in `:flying`, seed window via `:approx-offset` (WS2.2) instead of real backfill (this is the WS2.3 re-flow's off-switch) |
| committed prune | `commit-pass!` + `prune-committed` | render.cljd:232-258 / 223-230 | geometric prune to overscan window + segment-participant exemption; drop `committed-max-age`(204-207) |
| config plumbing | `sliver-collection*`, `SliverCollection`, `collection-sliver-widget`, `CollectionSliverWidget`, `update-render!`, `createRenderObject` | sliver_collection.cljd:357/66/…; render.cljd:333/1324/1338/1350 | add `:overscan`/`overscanConfig` alongside `:animate`/`animConfig` |
| public docs | `sliver-collection`, `collection-view`, `grid-view`, `collection` docstrings | widgets.cljd:2972, 2698, 2617, 2783 | document `:overscan` |

Note the box host (`collection`, box.cljd) is **not** virtualized and sets
`remainingCacheExtent = main-max` (box.cljd:61) with `userScrollDirection idle` (box.cljd:53) — it has
no scroll, so overscan is a **no-op there**; the option is accepted-and-ignored for the box host (keeps
the two hosts' option surface uniform).

---

## 8. Impact on WS2.2's protocol shape

Overscan **consumes** the two hooks WS2.2 already plans (`:renewal-point?`, `:approx-offset`) and adds
**no new layout-protocol surface** — overscan is an engine + config concern; layouts stay pure
geometry. Two clarifications to fold into WS2.2:

1. **`:approx-offset` is load-bearing for the *flying* regime, not only cold windows.** WS2.3 framed it
   for "large initial offset / list→wrap remap." §4 here adds: a far `animateTo`/fast fling places the
   visible window via `:approx-offset` **every fly-by frame**. So `:approx-offset` must be genuinely
   O(1) (no walk) — already its contract — and its fallback (linear from `:estimate`) must be
   acceptable for fly-by placement (it is: fly-by frames aren't scrutinised). No signature change; a
   note that this is a hot path in the flying regime.
2. **The run-boundary cache (WS2.3) is populated over the *overscan* window, not the viewport.** Tier-2
   real materialisation (§3) is what creates real run boundaries above/below the viewport; WS2.3's
   "refine on every dirty pass within the overscan window" is literally scoped by §2 here. WS2.2 should
   document that `:renewal-point?` checkpoints get retained for indices materialised in the overscan
   margin, so the sparse GC-surviving subset covers slightly beyond the screen (a hair more retention,
   still bounded by the overscan cap).

Neither is a breaking change; both are optional-with-derivable-defaults, exactly as WS2.2 already
frames them. `:anchor`'s two-mode re-documentation (exact at renewal points / dead-reckon elsewhere)
is unchanged from the WS2.3 recommendation.

---

## 9. Open questions (flagged)

1. **Advertised `cacheExtent` vs private overscan children (needs verify on device/framework).** The
   recommendation (§1, §2.5) is to keep `SliverGeometry.cacheExtent` on the *true* Flutter window and
   hold overscan children privately (the engine owns `collectGarbage`, so the framework won't collect
   them). Risk: `SliverMultiBoxAdaptorElement` keep-alive/semantics or the parent viewport's
   `remainingCacheExtent` accounting may assume advertised cacheExtent covers all live children.
   Alternative: advertise cacheExtent covering the overscan region (framework-consistent) at the cost
   of **stealing cache budget from sibling slivers** and inflating reported cache. **Resolve by
   testing** whether private extras survive a GC pass and whether semantics/keep-alive misbehave.
2. **Velocity from `scrollOffset` delta vs a host-pushed velocity — confirm accuracy under
   correction.** The correction-compensation (§2.1) is the fragile bit: if a frame emits a correction
   *and* real motion, the subtraction must be exact or velocity spikes. Needs a device check that
   WS1.3/WS2.3 corrections don't leak into the velocity estimate. Fallback if it proves noisy: a
   host-side `NotificationListener<ScrollNotification>` feeding `ScrollUpdateNotification` deltas into
   `update-render!` (more plumbing, cleaner signal).
3. **Thresholds are guesses (`v_slow`, `v_fast`, `lookahead-ms`, `ahead-max`, EMA α, hysteresis band).**
   All the numbers in §2.2/§2.3/§5 are reasoned starting points, not tuned. They need on-device
   calibration (especially the flying-regime cut-in — too low and short scrolls lose real layout, too
   high and long flings do wasteful real work). Flag for the implementation + verify pass.
4. **Landing detection precision.** "Velocity decaying **and** offset delta shrinking near a stable
   target" (§2.3) is a heuristic; a fling that stops mid-list (not near element 0) must still land on
   real geometry without a visible snap. Confirm the landing regime engages purely on deceleration,
   independent of *where* it lands, and that WS2.3's exact-at-0 branch is a special case of it, not a
   separate path.
5. **Interaction with reorder / box host.** Overscan is a no-op for the box host (§7) and reorder uses
   `SliverReorderableList` (list-only) — confirm the reorder path (sliver_collection.cljd:225-233)
   doesn't need overscan wiring or breaks if the option is present-but-ignored.
```
