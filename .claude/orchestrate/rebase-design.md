# Rebase design — first-class coordinate rebase with an exactly-once protocol

**Task**: collection-arch-hardening, step 8. Consumer: step 9 (implementation),
step 10 (red verification + TEMP revert), step 11 (capture-mode).
Engine: `src/flutter_cljd/internal/collection/render.cljd`; tween:
`internal/collection/tween.cljd`; host: `internal/sliver_collection.cljd`.

The uncommitted `seed-cache!` WIP in the working tree is SUPERSEDED by this
design (its intent lands as stage 1, without the NEW-4 crash). Stage 0 reverts it.

---

## 1. Problem

A *rebase* is a coordinate-frame event: the engine decides that the content the
viewport is visually anchored to now lives at a different content-space offset
(new-frame = old-frame + delta), and the scroll offset must move by the same
delta so the visible pixels do not. Today five mechanisms each implement a slice
of this with their own emission path and their own re-entrancy one-shots:

| mechanism | produce | transport | absorb | one-shot state |
|---|---|---|---|---|
| cross-layout re-anchor | `seed-cache!` | `reanchorShift` field → correction geometry | corrected re-run | `crossLayoutReanchor`, `reanchorShift` |
| seam refinement | `reflow-from-checkpoint!` | local `acc` in `backfill-leading!` → correction | corrected re-run re-walk | (shifted-frame discipline) |
| exact landing | `landing-decision` | correction | forced from-`:init` seed | `landingEmitted`, `pendingLandingWs` |
| top-underflow | `underflow-decision` | correction | corrected re-run re-seed | (drops cache in place) |
| segment set-point | `segment-start!` | `segAnchor` → incremental corrections | telescopes over segment | `segPrevDesired` |

No structure states "a rebase of Δ happened this pass and was consumed exactly
once". The far-morph bug chain is precisely a rebase produced (reanchorShift
88778) whose emission was overwritten (`from-relay!`) and whose absorption
channel (`segAnchor`) was nulled by an unrelated defect (`restingTop` nil) —
the delta silently vanished.

## 2. The model

Three first-class pieces replace the scattered state:

- **`viewAnchor`** — the single owner of "what the viewport top is looking at".
  Never nil while children exist. Replaces `restingTop` and the flow/indexed
  capture duplication (`anchor-before`-as-resting-capture vs
  `capture-resting-top!`).
- **pass rebase accumulator** — one per-pass sum of every coordinate-frame
  delta the pass produced (cross-layout re-anchor, seam deltas, underflow
  deficit, landing residual, anchor-follow). Replaces `reanchorShift` and
  `backfill-leading!`'s return-value plumbing. A double field reset each pass
  (like `passMaterialized`), written only through `rebase+!`.
- **`pendingRebase`** — the cross-pass record arming the corrected re-run's
  seed decision. Replaces `crossLayoutReanchor` + `landingEmitted` +
  `pendingLandingWs` (three one-shots → one record with an explicit lifecycle).

Field delta: `restingTop`, `crossLayoutReanchor`, `reanchorShift`,
`landingEmitted`, `pendingLandingWs` (5) are deleted; `viewAnchor`,
`passRebase`, `pendingRebase` (3) are added. `segAnchor` is upgraded in place.

### 2.1 `viewAnchor`

```clojure
{:idx    int      ;; index of the first attached child whose span reaches past scrollOffset
 :key    any      ;; key-for(idx) at sample time — survives index shifts (count changes)
 :off    double   ;; its content-space offset, in the pass's CORRECTED frame
 :extent double   ;; its measured main extent
 :frac   double}  ;; clamp((so' − off) / extent, 0, 1) — fraction consumed above the viewport top
```

**WHEN**: sampled in the `performLayout` epilogue of EVERY resting pass
(`tweenAnim` nil) — including passes that emitted a correction and passes that
re-seeded the cache. Mid-segment passes never overwrite it: `segment-start!`
consumed it, and it stays the pre-segment truth until the first post-segment
resting pass resamples.

**FROM WHAT**: the attached children's post-layout parentData (layoutOffset +
measured extent) — never the pre-walk cache. The comparison offset is
`so' = scrollOffset + (this pass's emitted correction or 0)`, so the stored
`:off`/`:frac` live in the frame the next pass will observe.

**Invariant (kills known-A's root)**: after any pass that leaves at least one
sized attached child, `viewAnchor` is non-nil. Debug assert in the epilogue.
The far jump's re-seeding last pass has attached children, so it samples fine;
the current nil comes only from reading `.-cache` after `inverse-seed!` emptied
it and before `walk!` refilled it.

`anchor-before` keeps its OTHER role — the pre-walk snapshot that
`anchor-delta` diffs against post-walk offsets to follow content reflow. That
is a same-pass mechanism, not the cross-pass anchor, and it now feeds the same
accumulator (below) instead of a separate emission branch.

### 2.2 Pass rebase accumulator

Every producer adds its delta via one helper and continues working in the
shifted frame (the discipline `backfill-leading!`'s `wsx` and `reanchor-band`
already follow — now universal):

- cross-layout re-anchor: `approx − stale` (from `seed-cache!`);
- each reconciled seam: `refine-seam-delta` result (shift-attached + continue);
- top-underflow: the manufactured `deficit`;
- exact landing: the `−z` residual;
- anchor-follow: `anchor-delta`'s post-walk drift (folded last; no overlap with
  the re-anchor because the pre-walk snapshot is taken after seeding).

The gating kernels (`refine-emit?`, `origin-refine-emit?`, `landing-decision`,
`underflow-decision`) are unchanged — they decide *whether* a producer runs;
the accumulator only unifies *transport*. `reanchor-band` generalizes to "the
band this pass materializes against is [ws, we] + passRebase-so-far".

### 2.3 Publication — exactly one consumer per pass

At pass end, if `passRebase ≠ 0`, exactly one arm runs, selected by the pass
mode (not by field-precedence `cond` chains):

1. **:emit** (resting pass) — geometry := `scrollOffsetCorrection passRebase`;
   `pendingRebase` := the armed record (§2.4). This is today's behavior for all
   four mechanisms, now through one exit.
2. **:absorb** (segment capture pass) — NO correction geometry is ever written.
   The capture driver runs in capture context and RETURNS the rebase as a
   value; `segment-start!` folds it into the set-point: the `to` side of
   `segAnchor` is evaluated in the post-rebase frame (the capture walk placed
   the target window there), so the set-point's telescoped total inherently
   carries the delta, and the incremental corrections deliver it over the
   segment. `pendingRebase` stays nil — the segment owns the delta now.
3. **:value** (step-11 capture-mode) — same as :absorb but the caller is
   whoever requested the capture; no engine state is armed.

**Exactly-once assert** (debug, `performLayout` epilogue):

```
produced := passRebase
consumed := (:emit  — geometry.scrollOffsetCorrection == produced)
          | (:absorb — capture context took the value AND segAnchor non-nil)
assert (zero? produced) or (exactly one consumed arm ran)
```

Concretely: a resting pass asserts `scrollOffsetCorrection` equals `passRebase`
whenever either is non-zero; a capture pass asserts the final displayed
geometry (`from-relay!`'s) carries NO correction and that
`(and (not (zero? produced)) (nil? segAnchor))` is impossible. The second
clause is the far-morph bug stated as an assertion: a produced rebase with no
absorption channel can no longer vanish silently — it fails loudly in debug.

### 2.4 `pendingRebase` — the corrected re-run contract

Emission and absorption-by-re-run happen in different passes; the record makes
that handoff explicit instead of encoding it in three one-shot booleans/NaNs:

```clojure
{:delta     double            ;; what was emitted
 :cause     :cross-layout | :seam | :landing | :underflow | :anchor-follow
 :seed      :init | :anchor   ;; how the corrected re-run must seed
 :basis     {:idx int :off double}  ;; for :anchor — where the frontier now lives (nil for :init)
 :expect-ws double}           ;; window start the re-run must arrive at (± slack)
```

Lifecycle: **armed** on :emit → **consumed** by the next `seed-cache!`, which
unconditionally clears the field and then:

- window within `slack` of `:expect-ws` → seed as demanded (`:init` = the
  landing contract's from-`:init` walk, exactly `landing-reseed-decision`'s
  `:init-seed`; `:anchor` = dead-reckon at `:basis`);
- window elsewhere → **abandoned**: seed normally. Because the record is gone,
  a later honest landing can emit again — no `landingEmitted` disarm dance;
  `:violated` remains detectable as "the re-run consumed a `:landing` record
  and index 0 is still off zero" (same assert text as today, but keyed off the
  consumed record instead of a surviving boolean).

Any state mutation the seed demands (dropping cache/checkpoints for `:init`)
stays at the EMIT site as today — the record does not change when cache is
dropped, only who decides the re-run's seed and how abandonment works.

### 2.5 `segAnchor` v2 — the fraction-preserving set-point (NEW-1)

**Invariant**: *the content point at fraction `frac` of the anchor item's
extent stays at the viewport's top edge throughout the segment* — the anchor
item can resize but never slides out from under the viewport top.

```clojure
{:from-off double  :from-ext double   ;; viewAnchor :off/:extent (old frame)
 :to-off   double  :to-ext   double   ;; the same index slot in the target frame
 :frac     double}                    ;; viewAnchor :frac
```

`desired(t) = lerp(from-off, to-off, t) + frac * lerp(from-ext, to-ext, t)`,
with `desired(0) = from-off + frac*from-ext = the segment-start scrollOffset`
by construction of `frac` — the incremental machinery
(`point-correction-delta`, `segPrevDesired`, drag composition, telescoping) is
unchanged except for the `desired` formula. When the anchor's extent does not
change, this reduces exactly to today's rigid `to − from` (pixel-identical),
so at-rest zero-config behavior is untouched.

The `to` side is read from the post-capture placed geometry — the flow cache
entry at the anchor index (the anchor is in-window by definition, so the
capture walk just placed it) or `geo-fn` for an indexed target — never from
the frozen snapshot lookup whose nil used to null the whole set-point. If the
anchor key is `exiting` (its item is being removed), the anchor re-picks the
first non-exiting slot at/after the viewport top before freezing.

## 3. Behavior decisions

### 3.1 Segment `to` extent is supplied, never read back (Step-5 a)

`live-only-flow-window` stops reading `(.-scrollExtent (.-geometry rs))`
entirely. A new pure helper `flow-total-extent` factors the total-extent
expression out of `finish-flow-geometry!` (walked end → `:estimate` /
`:approx-offset` extrapolation over item-count) and the capture supplies the
snapshot `:extent` explicitly from the TARGET layout:
`max(frames-main-extent, flow-total-extent(target) − dying-ext)`. This fixes
both halves of the step-5 finding: the zero read (correction-only geometry)
and the "honest old total is still wrong" caveat — the `to` endpoint's extent
is the target's estimate by construction. `fromExtent` keeps reading the last
resting geometry (honest). This also removes the capture path's dependence on
`.-geometry` side effects — a prerequisite for step 11's value-only mode.

### 3.2 Tween windowing degrades sanely; segments have a validity domain (NEW-5, NEW-7)

Two layers:

1. **No fabricated indices.** `keyed-tween-layout`'s `:first-index` /
   `:last-index` stop mapping a no-opinion (nil) answer to 0. A frozen
   snapshot with no opinion at `off` answers the CLAMPED snapshot edge: `off`
   past the end → the last snapshot index; `off` above the start → the base
   index. The returned window is therefore always well-ordered and inside the
   snapshot. `indexed-layout!` additionally refuses inversion defensively:
   debug assert `first-index <= target-last` and clamp `target-last := max`
   in release — Flutter's garbage-count assert (NEW-5) becomes unreachable
   from this layout regardless of source bugs.
2. **Domain exit ends the segment.** A frozen segment is valid only in the
   neighborhood of its frozen window; it never extrapolates. `segment-start!`
   stamps the tween with `:domain` = the snapshot's content span shifted by
   the set-point (start → end union). In the SEGMENT CONTINUE branch, if the
   current window lies more than one window extent outside the domain (the
   `window-far-from-band?` criterion), the engine settles the segment in
   place: `segTween`/`segAnchor` := nil, `segGen` := `curGen`, `leaving` :=
   {}, `clips` := {}, and the pass dispatches to the resting driver, whose
   normal far-jump inverse seed takes over — O(window), band follows pixels.
   The host clock keeps playing to its end; `finish-segment!` then finds
   nothing left to settle. (Open question 1 covers host notification.)

This closes NEW-7's whole chain (band pinned at 0 for the segment, then the
layout swap walking the 16 000 px gap) and NEW-5's inverted-window crash.

### 3.3 Count change: drop the cache, re-anchor by anchor key (known-B)

`update-render!` on a count change currently keeps the geometry cache while
every index shifted — the cache answers for the wrong items (known-B). Index
remapping by diff info is rejected: the diff runs only for keyed+animated
collections, so a general remap is unavailable. Instead:

- a count change DROPS the cache (checkpoints/anchoredTo0 already drop);
- the next `seed-cache!` resolves the anchor by KEY: walk the attached
  children (O(window)) for the child whose `key-for` equals `viewAnchor :key`
  — keyed reconciliation moved that element to its NEW index — and anchor the
  cache at `{new-index, viewAnchor :off}`. The visual anchor holds (O6: key
  and intra-offset preserved); the coordinate estimate above the anchor is now
  stale by the inserted/removed extent, which is exactly a rebase and is
  repaired by the standard seam/underflow/landing producers on subsequent
  resting passes — O(window) per pass, no O(n) walk.
- anchor key not found attached (item removed, or nothing attached) → today's
  fallback chain (firstChild stale-offset dead-reckon / inverse seed).

### 3.4 Cross-layout re-anchor anchors the viewport top, not firstChild

`seed-cache!`'s cross-layout dead-reckon re-anchors at `viewAnchor` (idx +
old off), not at `firstChild` (the overscan head) — the not-animated far-morph
defect (anchor 1346 → 1340 by the overscan margin). Attached children above
the anchor index are GC'd first with the guards the WIP lacked: nullable
re-read of `firstChild` after `collectGarbage`, and skip the GC when it would
empty the window (NEW-4's crash shape). The `approx − stale` delta goes
through `rebase+!` like every other producer.

### 3.5 NEW-2 and NEW-3 are freestanding — explicitly out of the rebase structure

Both are wrap-tiling defects in the seam KERNELS, not in rebase ownership or
transport; the unified accumulator transports whatever delta the kernel
computes. Step 9 fixes them as separate stages (8–9 below):

- **NEW-2**: `refine-seam-delta` declares convergence on the main axis alone;
  an open-run prefix matches on `:offset` while the head's `:cross` is stale.
  Fix shape: the seam is converged only when the re-flowed placement matches
  the cached head on BOTH `:offset` and `:cross` (and otherwise re-anchors /
  emits as today). Signature grows a cross argument; the stitch keeps the
  re-flowed entry.
- **NEW-3**: `leading-step-entry`'s advance collapses to `ext` for wrap
  because `:place` puts `i` into the `:anchor`-seeded empty run at `head`, so
  `flow-end` returns `head`. Fix shape: the backward advance must include the
  inter-run gap — derive it by seeding the anchor at `head`, placing `i`, and
  measuring `flow-end` AFTER also accounting for the run break (e.g. place a
  probe successor, or extend the layout contract with the run-advance the
  renewal hook already implies). Exact mechanism is stage 9's call; the design
  constraint is only that in-window geometry stays canonical (`:canon`).

## 4. Defect → design element map

| red test | design element |
|---|---|
| resting-top-lost-after-far-jump (known-A) | §2.1 viewAnchor: epilogue sampling from attached children |
| far-morph-wrap-to-list-not-animated | §3.4 cross-layout re-anchor at viewAnchor |
| far-morph-wrap-to-list-animated | §2.3 :absorb (no-emit capture) + §2.1 + §2.5 to-side from placed geometry |
| far-morph-capture-extent-truncated | §3.1 explicit target extent |
| fuzz-morph-loses-anchor-without-jump (NEW-1) | §2.5 fraction-preserving set-point invariant |
| fuzz-garbage-counts-exceed-child-count (NEW-5) | §3.2 layer 1 (clamped windows, no inversion) |
| fuzz-count-change-then-far-jump-walks-dataset (NEW-7) | §3.2 layer 2 (domain exit) |
| insert-remove-above-window-anchor (known-B) | §3.3 count-change key re-anchor |
| fuzz-null-render-box-cast-in-morph-after-remove (NEW-4) | stage 0 WIP revert; §3.4 redoes the intent guarded |
| fuzz-wrap-runs-overlap-after-back-drag (NEW-2) | §3.5 freestanding (stage 8) |
| fuzz-wrap-run-advance-short-after-morph (NEW-3) | §3.5 freestanding (stage 9) |

## 5. Step-9 stages

Each stage ends with the default suite green (`-- -x "known-red || fuzz"`),
the listed reds flipped green (and re-tagged out of `known-red`), all
previously-flipped reds still green, and O7 (cold-vs-warm equivalence) run
explicitly. Ordering puts the structure early so behavioral fixes land inside
it once, and prefers deleting one-shots over adding state.

0. **Revert the working-tree `seed-cache!` WIP** (uncommitted diff only; the
   TEMP instrumentation commit fb24f35 stays until step 10). Flips: NEW-4.
1. **viewAnchor** — add the field + epilogue sampler for both drivers; delete
   `restingTop`, `capture-resting-top!`, the `anchor-before`-based resting
   write; `segment-start!` reads viewAnchor. Rework the cross-layout
   dead-reckon per §3.4 (guarded). Flips: resting-top-lost,
   far-morph-not-animated.
2. **passRebase + pendingRebase** — behavior-preserving consolidation: all
   producers publish through `rebase+!`; one emission site; delete
   `reanchorShift`, `crossLayoutReanchor`, `landingEmitted`,
   `pendingLandingWs`; add the exactly-once epilogue assert. Flips: none by
   design (transport unification); the gate is the full green suite + O7.
3. **No-emit capture + explicit extent** — capture context on the flow driver
   (`:absorb` arm returns the value), `from-relay!` asserts no correction,
   `flow-total-extent` supplies the snapshot `:extent`, segAnchor `to` from
   placed geometry, absorption assert armed. Flips: far-morph-animated,
   capture-extent-truncated.
4. **Fraction set-point** — segAnchor v2 + `desired` formula (§2.5). Flips:
   NEW-1.
5. **Tween window sanity + domain exit** (§3.2). Flips: NEW-5, NEW-7.
6. **Count-change anchor rebase** (§3.3; viewAnchor `:key` is already stored
   by stage 1). Flips: known-B (insert-remove-above-window-anchor).
7. **Fuzz re-campaign** — rerun the step-6 campaign (all three profiles) over
   the unified engine before the freestanding fixes; new signatures triage as
   in step 7.
8. **NEW-2 seam cross reconcile** (freestanding, §3.5). Flips: NEW-2.
9. **NEW-3 leading run advance** (freestanding, §3.5). Flips: NEW-3.

Complexity budget: every added mechanism is O(1) per pass (accumulator,
record, domain check, fraction math) or O(window) (anchor sampling, key
resolve) — the materialization tripwire stays the enforcement.

## 6. Step-11 boundary (capture-mode)

The `:value` arm of §2.3 IS the interface: capture-mode = run the driver with
an explicit mode where the emit arm is statically unreachable and the pass
returns `{:rebase <record> :extent <double> :window <placed frames>}` as plain
values. The rebase record is already a plain map with no reference to mutable
fields, and §3.1 removed the capture path's `.-geometry` read-back — so step
11 changes dispatch (mode selection, corrections statically off), not the
record shape. Stage 3 should introduce the mode argument (even if two-valued
at first) so step 11 extends an existing axis instead of reworking call sites.

## 7. Open questions

1. **Domain exit vs host clock**: the engine settles silently while the host
   clock plays out; dying cells (off-screen after the far jump) snap at
   `finish-segment!`. Is a host callback to stop the clock early worth the
   coupling, or is the silent settle acceptable? (Default: silent settle;
   revisit if a fuzz signature shows visible artifacts.)
2. **Fraction set-point endpoint source**: `desired` lerps frozen endpoint
   extents itself; should it instead read the anchor's frame from `segTween`
   (identical for linear stay-cells, differs if the anchor is a leave-slide)?
   Default: frozen endpoints + the exiting-anchor re-pick rule.
3. **Count-change fallback quality**: when the anchor's child did not survive
   reconciliation, the firstChild dead-reckon may still drift O6's intra
   check; needs a dedicated fuzz profile (`:no-layout` with heavy churn) after
   stage 6.
4. **`:expect-ws` slack**: one viewport (today's `remainingPaintExtent`) kept
   for all causes, or derived per cause? Default: keep one viewport.
5. **O9 calibration**: `commit-prune?`'s `corrected?` input becomes "the emit
   arm ran" — semantically identical, but the mid-segment committed bound from
   step 7 should be re-verified in stage 7's re-campaign.
