# WS1.2 — One-frame freeze at retrigger: research

Read-only investigation. Answers the user's open question from
`collection-scroll-anim-followups.md` §1.2: **is `.forward .from 0.0` actually
necessary?** Context: WS1.1 (commit a38bbbd) landed — each cell's `from` is now
re-captured as its current *visible* window at `segment-start!`, so the t=0 frame
after a retrigger renders pixels identical to the previous frame. The "freeze" is
therefore a *duplicated frame*, not a visual snap.

**Verdict up front: `.from 0.0` is necessary given the current single-controller
design; dropping it produces a large visible SNAP, not a smoother start.** The
freeze is the price of the jump-free re-aim, and it is masked by the ease-in
curve. Recommend **accept (option a)** unless the device check shows a visible
hitch; the cheapest real fix if one is needed is a one-frame pre-advance (b2).

---

## 1. Pipeline trace (retrigger of an in-flight segment)

Files: `sliver_collection.cljd` (host), `collection/render.cljd` (engine),
`collection/tween.cljd` (pure tween), `collection/box.cljd` (box twin).

### The retrigger frame (call it frame N)

1. Parent rebuilds with new data → `didUpdateWidget` runs **during frame N's build
   phase** (`sliver_collection.cljd:114-139`). It classifies enter/exit
   (`apply-diff*`) and calls `start-segment!` (`:135`).
2. `start-segment!` (`sliver_collection.cljd:329-358`): bumps `tweenGen` (`:355`),
   sets `duration` on the **controller** (`:354`) and — first time only — builds a
   `CurvedAnimation` with the curve (`:352-353`), then calls
   **`(.forward ctrl .from 0.0)`** (`:357`).
3. `AnimationController.forward(from: 0.0)` does, synchronously:
   - `value = 0.0` → the setter calls `stop()` (kills the old segment's ticker),
     `_internalSetValue(0)`, then **`notifyListeners()` synchronously**.
   - The `CurvedAnimation` forwards that notification → the render object's
     listener (`render.cljd:307`, `set-tween-anim!`) calls **`markNeedsLayout`**
     while we are still inside frame N's build.
   - `_animateToInternal(1.0)` starts a **new** ticker → `ticker.start()`
     schedules a frame callback for the **next** vsync (frame N+1).
4. Frame N's layout phase runs (markNeedsLayout was honored this frame). First
   pass of the new gen → `segment-start!` (`render.cljd:1108-1221`):
   - `committed-from := committed-rects rs` (`:1128`) — the OLD committed rects,
     which after WS1.1 hold each cell's **visible window** (`committed-rects` doc,
     `render.cljd:246-252`; the visible extent is written every pass by
     `commit-pass!`).
   - builds `seg-tween` via `keyed-tween-layout`, whose `progress` thunk reads
     `(.-value tweenAnim)` = `CurvedAnimation.value` = `curve.transform(0) = 0`
     (`render.cljd:1210`, consumed at `tween.cljd:319`).
   - `keyed-tween-layout` geometry (`tween.cljd:321-362`): `from* = committed`,
     `to* = target`, `fr = lerp-frame(from*, to*, t=0) = from*` (`:356`,
     `lerp-frame` at `tween.cljd:56-67`). So every cell lays out at exactly its
     committed rect.
   - **Frame N paints the committed rects = exactly what frame N-1 painted.**
     Duplicate #1.

### The first ticker tick (frame N+1)

5. The ticker callback scheduled in step 3 fires in frame N+1's transient-callback
   phase. **`Ticker` sets `_startTime` to this first timestamp and delivers
   `onTick(Duration.zero)`** — the well-known Flutter rule that a Ticker's first
   tick is elapsed-0 (so animations don't jump on start).
6. `AnimationController._tick(0)` → `value = simulation.x(0) = lowerBound = 0` →
   `notifyListeners()` (fires unconditionally, even though value is unchanged) →
   `markNeedsLayout`. Frame N+1 lays out at t=0 again → **paints the committed
   rects again.** Duplicate #2.

### Motion begins (frame N+2)

7. Second ticker tick: `onTick(≈dt)` → `value = simulation.x(dt) > 0` → first
   frame with visible motion.

### Correction to the doc's count

The doc §1.2 says "first Ticker tick lands next vsync → **one** frozen frame."
That undercounts by one: because the Ticker's first tick is **elapsed-0**, the
controller sits at value 0 for **two** consecutive painted frames (the seg-start
layout of frame N and the elapsed-0 tick of frame N+1). Motion resumes on the
**second** ticker tick (frame N+2). So the stall is **~2 vsyncs** (≈33 ms @60 Hz,
≈16.7 ms @120 Hz), not one. This is inherent to any Ticker-driven controller and
holds identically for the box twin (`box.cljd:560-586`, same `.forward .from 0.0`
at `:586`; same elapsed-0 first tick).

### What breaks if `.from 0.0` is dropped

Drop it → `(.forward ctrl)` with the controller still at, say, 0.6:

- `forward()` does **not** reset value; it runs `_animateToInternal(1.0)` from the
  current 0.6. It also scales the remaining duration by the remaining fraction —
  `simulationDuration = duration * (1.0 − 0.6) = 0.4·duration` — so the retriggered
  segment would play in **40 % of its configured time**.
- Far worse than the shortened duration: the `progress` thunk feeds the raw
  controller value through the **same** `CurvedAnimation`, so the very first pass
  of the new segment reads `t = curve.transform(0.6) ≈ 0.68` (easeInOut). But the
  new `from` was just re-captured as the cell's **current** position and the new
  `to` is the **new** target. `lerp-frame(from, to, 0.68)` places every cell
  **68 % of the way from its current position to the new target on the first
  frame** → an immediate, large **visual SNAP**.

The freeze and the jump-free re-aim are **coupled**: capturing `from` fresh at the
current position is only consistent with `t = 0`. Resetting the controller to 0 is
what makes the fresh capture correct. Duration and curve live on the controller /
CurvedAnimation, **not per-segment** (`sliver_collection.cljd:353-354`), so there
is no per-segment state to re-base a nonzero start against. **`.from 0.0` is
necessary in this design.** The one-frame (really two-frame) freeze is a symptom
of a deeper property — the per-cell velocity resets to the curve's ease-in
velocity (~0) on every retrigger — which only velocity handoff (option c) removes.

---

## 2. Options analysis

### (a) Accept as-is — **recommended default**

- **Mechanics:** nothing changes. WS1.1 already removed the snap; the residual is
  ~2 duplicated frames at the retrigger instant, only when retriggering a segment
  that is genuinely mid-flight.
- **Perceptibility:** the segment curve is an ease-in-out (or ease-in) — velocity
  at t=0 is ~0, so the first ~15 % of a *correct* segment barely moves anyway.
  Freezing at t=0 for 2 frames overlaps the curve's already-near-stationary
  opening, so the missing motion is a fraction of a pixel-per-frame delta. At
  60 Hz that is ≈33 ms of "no additional motion" layered on a slow-start phase; at
  120 Hz ≈16.7 ms. With from-continuity in place this reads as a faint hitch, not
  a jump.
- **Platform precedent:** this is exactly what Flutter's own
  `ImplicitlyAnimatedWidget` family does. On a retarget it re-evaluates each
  tween's `begin` from the **current interpolated value** and calls
  `_controller.forward(from: 0.0)` — identical re-capture-then-reset-to-0 pattern,
  identical elapsed-0 duplicated first frame. `AnimatedContainer`, `AnimatedAlign`,
  etc. all ship with this behavior and it is considered acceptable. UIKit's
  `UICollectionView.performBatchUpdates` retarget and `.beginFromCurrentState`
  likewise **stop the running CAAnimation and start a new one from the presentation
  value** — same class of one-frame re-seat.
- **Cost:** zero.

### (b) Restart at 0 but avoid the stall

- **(b1) carry remaining time (shorter duration on retrigger):** doesn't touch the
  freeze — it still starts at t=0 — and makes retriggered segments play faster,
  which is a *worse* motion feel (a burst of speed on rapid retriggers). Reject as
  a freeze fix.
- **(b2) pre-advance the controller by one frame's dt:** instead of `.from 0.0`,
  start at `from = curveInverse?`—simpler: keep `from-capture at current position`
  but seed the controller at the value that one frame of progress would reach,
  i.e. `value₀ = dt / duration` (≈0.016/duration, a few %). Then `t =
  curve.transform(value₀) ≈ 0` (ease-in ⇒ curve of a small input is *tinier*), so
  `lerp(from, to, t)` ≈ from → the retrigger frame's jump is sub-pixel, **and** the
  elapsed-0 first ticker tick already shows `value₀` rather than 0, so motion is
  visible one frame earlier. Net effect: collapses the ~2-frame stall to ~1.
  - **Mechanics in this codebase:** replace `(.forward ctrl .from 0.0)` at
    `sliver_collection.cljd:357` (and `box.cljd:586`) with a helper that sets
    `value` to `min(1, dt/duration)` then `.forward()` (no `from`). `dt` is one
    frame — read `WidgetsBinding.instance.window` refresh rate or just hardcode
    1/60. Need to guard `value ≤ 1` for very short durations.
  - **Cost:** ~5 lines × 2 hosts + a constant. Low. Risk: a few % of initial jump
    on layouts where `to−from` is large per cell (a long move) — the jump is
    `(to−from)·curve.transform(dt/dur)`, which for a 300 ms segment at 60 Hz is
    `curve.transform(0.055)` ≈ 0.008 for easeInOut ⇒ <1 % of the move. Negligible.
- **(d) same-frame re-capture + first tick:** they *already* land in the same frame
  (frame N does both the from-capture and a value-0 layout). The stall is not a
  scheduling gap; it is the value sitting at 0 across frame N and the elapsed-0
  tick of N+1. A post-frame-callback capture would only *delay* the capture, not
  help. The single lever that puts motion on the retrigger frame without a jump is
  a nonzero-but-tiny start value — which is exactly (b2). So (d) reduces to (b2).

### (c) Velocity handoff (SpringSimulation seeded with per-cell rect velocity)

- **What velocity is available:** the design (WS1 header) is **one shared
  controller** driving a single scalar `t` via `lerp` for all cells; there are no
  per-item controllers. The controller exposes one scalar velocity `v =
  dt-value/dt`. A given cell's pixel velocity at retrigger is derivable —
  `(to_i − from_i)·curve'(t)·v` — from the OLD segment's `from_i/to_i` and `t`
  before re-capture. But it is **per-cell and different for every cell** (different
  `to−from` vectors, and the new targets differ per cell too).
- **Why a single spring can't carry it:** one shared `t`-spring produces one
  eased scalar; mapping it back through `lerp(from_i, to_i, t)` gives every cell
  the *same* normalized velocity profile. To honor each cell's actual entry
  velocity you need **per-cell simulations** (each springs from its own
  position+velocity to its own new target). That is precisely the per-item
  controller model the WS1 header rules out ("single shared segment tween is
  kept"). It is a core rewrite of the animator, not a tweak to the segment start.
- **Cost:** high (architectural), and it violates the stated WS1 constraint.
  **Reject for WS1.** It is the *correct* answer to the underlying velocity
  discontinuity (see external §Core Animation), but belongs to a future
  per-cell-simulation design, not this workstream.

---

## 3. External findings — how other platforms retarget in-flight

**Core Animation — additive animations (the canonical smooth retarget).** Since
iOS 8 all `UIView` animations are **additive by default**. A retarget adds a
**new** `CAAnimation` expressed relative to the model: `fromValue = old − new`,
`toValue = 0`, `isAdditive = true`, with the layer's model value set to the new
target. The **old animation is not stopped** — it keeps running to completion and
is summed with the new one, so the presentation value continues from its current
position *and velocity* with no reset and no duplicated frame. This is the
mechanism Apple presented in WWDC14 "Building Interruptible and Responsive
Interactions" and that the *Seamless* framework packages. Crucially, it needs
**per-property animations that can coexist and sum** — the exact opposite of this
codebase's single shared scalar `t`. It maps to option (c), not (a)/(b).
`.beginFromCurrentState`, by contrast, **does** stop the old CAAnimation and start
a new one from the presentation value — the same re-seat this codebase does, and
it is a supported, widely-used option (keyboard-driven layout changes).

**UIKit / UICollectionView.** `performBatchUpdates` retargets by tearing down and
re-adding animations; without additive layering you get the same one-frame re-seat.
The smooth path is the additive default above.

**Flutter.** `AnimationController.forward({from})` sets value then
`_animateToInternal`, which **scales the remaining duration by the remaining
fraction** (confirming the "40 %" analysis in §1). `ImplicitlyAnimatedWidget`
retargets by re-capturing each tween `begin` from the current value and calling
`forward(from: 0.0)` — **the identical re-capture-and-reset pattern this codebase
uses, with the identical elapsed-0 duplicated frame**, and it is Flutter's
accepted standard behavior. Flutter has **no** additive-animation primitive; the
smooth-retarget idiom on Flutter is to swap to a **physics simulation**
(`AnimationController.animateWith(SpringSimulation(...))` seeded with the current
velocity) — i.e. option (c). Flutter issue #142218 (jitter/dropped frames) is
about scheduler/raster jank, not retarget semantics — not relevant here.

**CSS transitions.** The one platform that retargets *smoothly by default without
additivity*. Per CSS Transitions Level 2, interrupting a transition **picks up the
current interpolated value** and transitions from there, and for a **reversal** it
applies a "reversing shortening factor" + "reversing-adjusted start value" so the
handoff keeps the eased feel and doesn't overrun. But note: CSS still **resets the
timing function to t=0** on the new (retargeted) transition from the current value
— so on a *non-reversing* retarget CSS also restarts the ease from zero velocity,
i.e. it accepts the same velocity discontinuity this codebase has; it just never
shows a *duplicated* frame because the compositor samples the new transition at the
real elapsed time from the interruption instant (no Ticker elapsed-0 first frame).

**Synthesis.** Two families exist: (1) **additive/simulation** (Core Animation
additive, Flutter spring handoff) — truly continuous position *and* velocity,
requires per-property/per-cell state → this codebase's option (c), out of WS1
scope; (2) **re-seat from current value + restart the ease at t=0** (Flutter
`ImplicitlyAnimatedWidget`, UIKit `.beginFromCurrentState`, CSS non-reversing
retarget) — accepts a velocity reset, and on Ticker-based systems accepts a
duplicated frame. This codebase is squarely in family (2), matching Flutter's own
built-ins. The duplicated frame is **normal for this family**; only family (1)
avoids it, at the cost of the shared-controller architecture WS1 preserves.

---

## 4. Recommendation

**Accept as-is (option a) — pending the WS1.1 device check.** It is zero-cost,
matches Flutter's own `ImplicitlyAnimatedWidget` behavior, and the stall lands on
the curve's near-stationary opening where it is least visible. `.from 0.0` is
**necessary** and must stay.

**What to look for during the on-device WS1.1 check**, to decide accept vs fix:

1. **Rapidly retrigger a mid-flight segment** (insert-during-insert /
   remove-during-insert, tap-spamming the mutation). Watch a single cell that is
   clearly moving.
2. **Judge the retrigger instant specifically:** with from-continuity there must be
   **no snap**. The only artifact allowed is a brief *pause* (~2 frames) before
   motion resumes.
3. **Decide:**
   - If the pause is **not perceptible** at the target refresh rate (likely at
     120 Hz, and likely even at 60 Hz given the ease-in mask) → **ship accept**.
     Close WS1.2.
   - If a hitch **is** perceptible at 60 Hz → apply the **one-frame pre-advance
     (b2)**: seed the controller at `min(1, dt/duration)` instead of 0 at
     `sliver_collection.cljd:357` and `box.cljd:586`. ~5 lines per host, halves
     the stall, sub-pixel jump. This is the cheapest fix consistent with keeping
     the single shared controller.
   - Do **not** reach for velocity handoff (c) from WS1 — it requires per-cell
     simulations and contradicts the single-shared-tween invariant. If, after the
     device check, the *velocity reset* (not the freeze) turns out to be the real
     complaint, log it as a separate future item against a per-cell-simulation
     redesign, not WS1.2.

### Sources

- [Building Interruptible and Responsive Interactions — WWDC14 (additive animations)](https://wwdcnotes.com/documentation/wwdcnotes/wwdc14-236-building-interruptible-and-responsive-interactions/)
- [UIView.animate(withDuration:…) — beginFromCurrentState — Apple Developer](https://developer.apple.com/documentation/uikit/uiview/1622451-animatewithduration)
- [Seamless — relative additive Core Animation (from = old − new, to = 0, additive)](https://github.com/KevinDoughty/Seamless)
- [Core Animation: a deep dive (additive property, render server)](https://chromium.googlesource.com/external/github.com/material-motion/motion-animator-objc/+/refs/tags/v4.0.0/README.md)
- [AnimationController.forward — Flutter API (sets value, then _animateToInternal)](https://api.flutter.dev/flutter/animation/AnimationController/forward.html)
- [AnimationController — Flutter API](https://api.flutter.dev/flutter/animation/AnimationController-class.html)
- [CSS Transitions Module Level 2 — reversing shortening factor / interrupted transitions](https://www.w3.org/TR/css-transitions-2/)
- [transition-timing-function — MDN (interrupted transition picks up current value)](https://developer.mozilla.org/en-US/docs/Web/CSS/transition-timing-function)
- [flutter/flutter#142218 — jitter/dropped frames (scheduler jank, not retarget)](https://github.com/flutter/flutter/issues/142218)
