# Anchor-primary design — the scroll truth migration, right-sized

**Task**: collection-arch-hardening, step 12. Consumer: step 13 (implementation),
step 14 (cleanup/docs). Engine: `src/flutter_cljd/internal/collection/render.cljd`;
tween: `internal/collection/tween.cljd`; layouts: `internal/collection/layout.cljd`.
Precedent: Compose LazyList (`firstVisibleItemIndex` + `firstVisibleItemScrollOffset`
as state), UIKit's content-offset preservation around `reloadData`. Constraint that
shapes everything below: Flutter's sliver protocol demands absolute `scrollOffset`
and `scrollExtent` every pass, pins content space to the leading edge (offset 0 =
content start, no negative offsets), and offers `scrollOffsetCorrection` as the
only reconciliation valve.

---

## 1. Gap analysis — how much of anchor-primary already shipped

The founding intent was "привязка к элементу (id/index + intra-offset) как истина;
абсолютный оффсет — производная для индикатора". Measured against the engine at
HEAD (post step 11), most of that intent is de-facto shipped:

**Already anchor-primary:**

- **`viewAnchor`** `{:idx :key :off :extent :frac}` is sampled from DISPLAYED
  frames after EVERY pass, mid-segment included (10b's NEW-9 rule), with fraction
  semantics (`:frac` = consumed portion of the anchor, the thing a resizing morph
  must preserve). It is exactly Compose's `(index, intraOffset)` pair, plus a key.
- **Key re-anchoring on count change** (9b stage 6): the anchor's identity is its
  KEY; `seed-cache!` and `segment-start!` both resolve it through
  `attached-index-of-key` after reconciliation moved the element. Truth follows
  data identity, not the slot.
- **The segment protocol is fully anchor-primary**: `segAnchor` v2 +
  `tw/point-desired` DERIVE the scroll offset from the anchor's fraction point for
  the whole segment; the incremental corrections are the derivative of the
  anchor's trajectory, composing with drag. Capture rebases are absorbed OVER the
  anchor (`:absorb`, total since 10b), never emitted.
- **`passRebase`/`pendingRebase`** (9a): one accumulator, one emission site, an
  exactly-once assert. Transport is already centralized.
- **`anchoredTo0`** already concedes the founding premise: the engine explicitly
  models absolute offsets as *estimates* unless a from-0-exact chain exists —
  "absolute offset is derived" is written into the state model.

**Still offset-primary — precisely:**

1. **The resting flow pass is seeded by geometry, not by the anchor.**
   `flow-layout!` derives the window from `scrollOffset`, `align-start!`
   bsearches the cache by `ws`, `walk!` places downward from cached absolute
   offsets. `viewAnchor` is sampled in the epilogue — a *report* of what the pass
   did, not its input. `seed-cache!` picks anchors from geometry
   (`window-far-from-band?`, `invert-approx-offset` at `ws`), consulting
   `viewAnchor` only in the cross-layout/count-change branches.
2. **Leading-side reconciliation is distributed and moves the window.** A change
   in the estimated extent above the anchor is discovered by five producers
   (seam refinement, underflow, exact landing, cross-layout re-anchor,
   anchor-follow), and each repairs it by SHIFTING the window content
   (`shift-attached-from!`, cache rewrite, cache drop + `:init` re-flow) mid-pass,
   threading the shifted frame through `reanchor-band` / backfill's `wsx` loop.
   The anchor stays put only because every one of those shifts is paired with a
   correction — anchor stability is an *emergent property of five mechanisms
   agreeing*, not a construction. The campaign-3 `:o6` above-window count-change
   family (7 seeds, the largest residual) is exactly this class: the count-change
   re-anchor keeps the anchor's old offset while the standard producers repair
   the stale estimate above it, and the repairs drift the anchor's fraction point.
3. **The origin machinery is three special cases of one fact.** `landing-decision`
   + `landing-reseed-decision` + the `pendingRebase` `:landing`/`:init-seed`
   dance, `underflow-decision`'s manufactured deficit, and `anchoredTo0`'s ~14
   write sites all encode "the leading estimate changed / became exact" as
   separate mechanisms with their own lifecycles.

**Right-sizing.** The plan's original framing ("migrate the engine to
anchor-primary scroll truth") is now over-sized: the truth records, the key
identity, the fraction set-point, and the absorption channel — the hard,
bug-prone 80% — shipped in steps 9–10b and survived two fuzz campaigns. What
remains is **making the resting flow pass mono-sourced on that truth**: seed the
pass from the anchor, and collapse the five leading-side repair mechanisms into
ONE derived quantity ("the anchor's leading extent estimate") whose per-pass
delta is the only correction. Call it the *leading-side consolidation*, not a
truth migration. It is smaller than the framing but NOT skippable: the o6
residual family is structural to the distributed-repair model, and patching it
pointwise would be the sixth one-shot-style fix in the leading side — the exact
class this task exists to end.

**Explicitly not worth migrating (cost > payoff):**

- **Anchor-relative content space** (rewriting cache/committed/parentData offsets
  relative to the anchor). The framework pins parentData `layoutOffset`, paint
  math and the geometry protocol to absolute content space; a shadow coordinate
  system buys no observable behavior and doubles every read/write. Rejected.
- **`flying-flow!`**. Flying deliberately abandons anchoring for approximate
  fly-by frames nobody scrutinizes; there is no visual anchor to preserve at
  fling speed. Its landing already re-derives everything. Untouched.
- **The indexed driver's interior.** Indexed layouts have exact O(1) leading
  extent — the "estimate delta" is identically zero, the epilogue is the
  identity. Only the shared epilogue/sampling paths change.

## 2. Target model

### 2.1 The truth equation

At every resting/settled pass end the engine maintains:

```
scrollOffset' (= scrollOffset + emitted correction)
  == lead(anchor.idx) + frac · extent(anchor) − screenResidual
```

where `lead(idx)` — the content extent above the anchor — is a DERIVED estimate
with a latched value and an exactness flag:

```clojure
;; replaces anchoredTo0's scattered writes + the landing/underflow machinery
leadEstimate  ;; {:off double  — the anchor's latched content offset
              ;;  :exact? bool — chained from a measured walk to index 0}
```

**Latching invariant** (the anti-oscillation rule): the anchor's `lead` changes
ONLY from (a) measurement — backfill walked/checkpoint-reflowed real runs above
the anchor, or reached index 0; (b) anchor reassignment — far jump, cold start,
cross-layout re-anchor; (c) an explicit invalidation (count/cross change).
Aggregate drift (EMA updates) NEVER moves a latched lead — it feeds only the
extent reported OUTSIDE the latched region and future reassignments. This is
today's de-facto behavior (offsets latched in the cache) stated as law; without
it the Δ-epilogue below becomes a correction storm.

### 2.2 The pass, anchor-seeded

Per resting/settled flow pass:

1. **Resolve the anchor**: `viewAnchor` by key (count change) or index; clamp to
   `item-count`. No anchor / anchor unresolvable / window far from the anchor's
   frame ⇒ **reassignment**: invert the estimate at the window (today's
   `inverse-seed!` math), `lead := {:off approx :exact? false}`, `frac := 0`.
   A far jump is not a repair — it is the one legal anchor teleport, driven by
   user intent.
2. **Lay out from the anchor outward, in the current frame**: seed
   `:anchor(idx, lead.off)`, `walk!` downward (unchanged mechanics), backfill
   upward. NO producer shifts attached content or emits mid-pass.
3. **Backfill measures instead of reconciling**: `reflow-from-checkpoint!`'s
   suppressed-branch mechanics become the ONLY mechanics — the produced prefix is
   re-anchored so its runs abut the head (`:run-start` alignment; exact, because
   a renewal-point `:anchor` state is a function of offset alone), and the
   discrepancy between the checkpoint-implied absolute offset and the
   anchor-relative placement is *accumulated into the lead delta* instead of
   shifting the window. The emit/frozen branch pair, `shift-attached-from!` as a
   seam tool, and the `wsx` telescoping loop all collapse. Reaching index 0 at
   measured offset `z` relative to the anchor sets
   `lead' := {:off measured :exact? true}` — landing is not an event, it is the
   estimate becoming exact.
4. **The Δ-epilogue — the single correction producer**:
   `Δ := lead'.off − lead.off`. If `|Δ| > 1e-3` and the mode/velocity gate
   allows: translate the window ONCE (attached offsets + cache + checkpoints,
   one pass over O(window) state), emit `Δ` through `emit-correction!`, latch
   `lead'`. Otherwise keep the old latch (estimate stays latent, exactly like
   today's frozen suppression). The gate keeps today's two levels: mid-band
   measurements (seams) defer under ballistic suppression; exactness Δ (backfill
   reached 0) ignores it — a fast up-fling must still land index 0 exactly (#B).
5. **Report**: `scrollExtent := lead.off + window + trailing estimate`
   (`flow-total-extent` survives with `lead` replacing the `anchored0` input).
   Sample `viewAnchor` (unchanged). The truth equation is a debug assert.

The sliver protocol makes step 4's translation unavoidable — content space is
pinned at the leading edge, so a lead change IS a window translation — but
anchor-primary makes it *one atomic epilogue op computed from one number*,
instead of five producers each pairing their own shift with their own correction.
The window's internal geometry is decided anchor-relative before the
translation, so the anchor structurally cannot drift against its neighbors: the
o6 count-change family's mechanism is gone by construction.

### 2.3 What dies, what survives

| subsystem | fate |
|---|---|
| `pendingRebase` — ALL three arms (`:landing`/`:init`, `:cross-layout`, `:count-change`) + `:consumed` marker + `:expect-ws` slack | **dies** — the anchor (key + latched lead) IS the cross-pass contract; no corrected re-run needs a seed instruction because the cache is correct after the translation (measured in-frame, translated in-frame) |
| `landing-decision`, `landing-reseed-decision`, the landing cache-drop + `:init` re-flow | **die** — landing = lead becoming exact in the Δ-epilogue; `:violated` becomes "lead exact but Δ ≠ 0 with no mutation" (same layout-contract assert, one site) |
| `underflow-decision` + the manufactured deficit | **dies** — a stalled backfill with `head ≤ 0` and a positive honest estimate is just `lead' > lead`, the ordinary Δ |
| `anchoredTo0` (~14 write sites) | **dies as a flag** — becomes `leadEstimate :exact?`, written by the epilogue and reassignment only; `frontier-after-prune` (the identity) goes with it |
| `backfill-leading!`'s seam machinery: the emit-vs-frozen branch pair, `shift-attached-from!` as producer tool, the `wsx` shifted-frame loop, `reanchor-band` | **die** — one measurement path, discrepancies accumulate into the lead delta |
| `anchor-before` / `anchor-delta` (anchor-follow) | **die** — the anchor is the seed; its in-pass delta is 0 by construction |
| `seed-cache!`'s cond forest | **shrinks** to: anchor resolve → reassignment (far/cold/degenerate) → count-change/cross-layout re-derivation feeding the same epilogue |
| `landing one-shot remnants` (the consumed-record `:violated` detection) | **die** with pendingRebase |
| `reflow-from-checkpoint!` | **survives, simplified** — measurement + `:run-start` re-anchor only (9c's alignment repair becomes the universal path) |
| checkpoint store (`cp-*`, `trim-front`, prune) | **survives unchanged** — it IS the lead-estimate memory; entries ride the epilogue translation |
| layout contract (`:init/:place/:anchor/:renewal-point?/:renewal-index/:run-start/:approx-offset`) | **survives unchanged** — no new hooks; `:anchor` + `:run-start` become more central |
| tween layer (`point-desired`, `keyed-tween-layout`, windows, `:domain`), segment protocol, `segment-start!`, capture `{:rebase :extent}` value, seg-tail | **survive unchanged** — already anchor-primary; the capture's `:rebase` value becomes "the Δ the capture pass computed", absorbed as today |
| committed map, clips, §7b/§7c, `widen-window!` | **survive unchanged** |
| velocity sensor, overscan, regime gate, `flying-flow!` | **survive unchanged** — suppression gates the Δ-epilogue as it gates refinement today; the epilogue must keep feeding `refineLastOffset`'s correction compensation or the sensor misreads Δ as user motion |
| `pass-mode`/`pass-caps` | **survive** — `:refine?`/`:origin-refine?` collapse into one `:lead-emit?` capability with the ballistic sub-gate carried by the epilogue; `emit-correction!` + the exactly-once assert unchanged |
| `walk!`, `align-start!`, `pre/post-gc!`, `trim-front!`, tripwires | **survive** — same mechanics, anchor-provided seed |

Field delta: −1 (`pendingRebase`) − `anchoredTo0` + `leadEstimate` = net −1,
plus the deleted kernels. Every stage should delete more state than it adds.

## 3. Boundary behaviors — specified up front

These are the places naive anchor-primary designs die; each gets an explicit
rule:

- **Top boundary / clamp at 0 / bounce.** The anchor walks to index 0 via
  backfill; lead becomes exact; after the final Δ the origin is stable and NO
  further corrections occur near the top — the boundary spring runs over a
  frozen frame. Rule: while the viewport is overscrolled past the top (pixels
  < 0; the sliver sees scrollOffset 0), the epilogue emits NOTHING except an
  exactness Δ — an estimate-shuffling correction mid-bounce restarts the spring
  (campaign-2's applyNewDimensions → goBallistic observation). Anchor sampling
  at so' = 0 yields the first child at frac 0 — well-defined.
- **Bottom boundary + shrink-below-window.** When content shrinks under a
  resting-at-max viewport, the RANGE clamp wins: pixels move, the anchor is
  re-sampled from the displayed frame — pixels-primary at the boundary, never
  the anchor dragging the clamp (the O6 recalibration in step 10 already
  encodes this: removals at `pixels == max` legitimately change the top item).
  `align-start-fallback`'s `:attach-last` (window past the content end during
  a bounce) survives verbatim. Content shorter than the viewport: anchor =
  index 0 frac 0, lead exact 0, no corrections — degenerate and stable.
- **Empty collection.** Nothing attached ⇒ `viewAnchor` keeps its last value
  (supports remove-all + undo), geometry zero, and the LATCH is invalidated
  (nothing measurable). On repopulation the resolve clamps `idx` to `n−1`; if
  the key is gone and pixels sit at 0, reassignment lands index 0 — restoring
  to top is the honest answer when identity is lost.
- **Infinite collections (`item-count` nil).** The anchor model is MORE natural
  here than offset-primary (Compose's own case): lead works from checkpoints /
  `:approx-offset` (gallop inversion already handles unbounded), trailing extent
  is infinity (`flow-total-extent` unchanged), no bottom boundary exists, and
  exact landing at index 0 still works. No special-casing needed — the design
  must simply never require `item-count` outside trailing extrapolation (it
  doesn't today).
- **Indicator honesty (§2.4 lineage).** `scrollExtent = lead + window +
  trailing-estimate`. The rules carry over verbatim: while lead is inexact the
  reported extent is raised toward the honest `:approx-offset` estimate, never
  below laid content (O5's horizon guards it); estimate refinements reach the
  indicator through Δ-paired corrections, so pixels and extent move together
  and the thumb fraction stays continuous; residual thumb drift during
  refinement is accepted — honest and monotone-ish beats frozen.
- **The `:settled` divergence (step 11's handoff) — dissolved, deliberately.**
  Step 11 preserved `:origin-refine?` OFF in `:settled`, so a jump-to-top taken
  mid-segment defers its exact landing until the host clock ends. Under
  anchor-primary there IS no separate origin capability — landing is the
  Δ-epilogue, and `:settled` already has `:emit-rebase?` true. This design owns
  that behavior change: exactness lands during the playing clock (a domain-
  exited segment runs the resting drivers anyway; no set-point competes — step
  11 verified that). `:overscan?`/`:prune-commit?` stay OFF in `:settled` —
  they are unrelated to the offset model and remain deferred conservatism.

## 4. Equivalence gate for step 13

- **Per stage**: default suite green (`-x "known-red || fuzz"`), O7 cold-vs-warm
  run explicitly, `bin/check` clean — unchanged from steps 9–11.
- **Byte-comparison protocol (step 11 demonstrated it works)**: before each
  stage, run the full campaign (220 episodes, four profiles, the fixed seed
  ranges) at the pre-stage commit and record the failing set as
  (seed, op-index, oracle, :why) tuples; after the stage, the set must be
  **byte-identical** — same seeds, same op indices, same oracles — EXCEPT for a
  stage that claims a residual family: there the diff must be *exactly the
  claimed seeds turning green* (strict subset; any new failure or any
  unexplained flip = stop). Record both sets in the stage's commit message.
- **New oracles anchor-primary needs** (added in stage 0, so the baseline is
  sharpest before the engine moves):
  - **O10 anchor-rest stability**: at rest, `viewAnchor`'s (key, frac) is
    bit-stable across idle pumps — the anchor never teleports at rest.
  - **O11 truth equation**: after every settled pass,
    `|anchor.off + frac·extent − so'| ≤ eps` (harness reads the engine map;
    today it holds by sampling construction, post-migration it is the model).
  - **O12 extent quiescence**: at rest with no mutations, `maxScrollExtent` is
    constant across pumps (the campaign-2 spring-restart observation as an
    oracle) — catches Δ storms and EMA-driven extent oscillation.
  - **O9 per-pass WORK probe** — the deferred design lands HERE (stage 0):
    bound `pass-materialized` per pass by `materialization-budget` (already the
    engine's own tripwire arithmetic, exposed through the engine map), and
    RELAX the mid-segment retention denominators (`cache-n`) that the landing-
    collapsed attached window miscalibrates. This closes the deliberately-open
    O9 calibration and should resolve the residual `:o9` families into either
    calibrated green or a sharp red before any engine change.

## 5. Step-13 stages

Each stage: default-green + O7 + byte-compared per §4; each deletes more state
than it adds where possible.

0. **Oracles first, no engine change**: O9 work-probe split, O10/O11/O12;
   re-baseline the campaign. Expected diff: only `:o9` reclassifications (the
   documented denominator family); the `:o6` and `:o5` residuals must survive
   as the sharp targets for stages 2/4. *(Optional, budget-permitting: one
   BouncingScrollPhysics batch — the harness is clamping-only today and §6's
   bounce rules are otherwise untested.)*
1. **Anchor-seeded pass, behavior-preserving**: the pass resolves the anchor
   first; `seed-cache!`'s branches re-expressed as resolve/reassign;
   `align-start!`/`walk!` seeded from the anchor's frame. At rest the anchor
   and the cache-head bsearch coincide, so the fuzz set must be byte-identical.
   Flips: none by design.
2. **The Δ-epilogue**: backfill producers accumulate lead deltas instead of
   shifting; one translate+emit site; delete `shift-attached-from!`-as-producer,
   the `wsx` loop, `reanchor-band`, the landing cache-drop/`:init` dance,
   `underflow-decision`, `pendingRebase` and `landing-reseed-decision`;
   `anchor-before`/`anchor-delta` die. **Claimed family: the `:o6` above-window
   count-change seeds (151, 111, 17, 205, 52, 55, 34)** — this stage's success
   criterion; if any survive, they are a different root and the stage gate is
   "no new failures + the family strictly shrinks", with the survivors becoming
   the next triage before stage 3 proceeds.
3. **`anchoredTo0` → `leadEstimate`, capability collapse, `:settled` exactness
   ON**: delete the flag's write sites and `landing-decision`;
   `:refine?`/`:origin-refine?` → `:lead-emit?`; new green test: a jump-to-top
   taken mid-segment lands index 0 exactly before the host clock ends (the
   step-11 divergence, now a pinned behavior). Fuzz: byte-identical except any
   seed whose vector hits the settled-landing path (enumerate by trace before
   flipping).
4. **The `:o5` remove-extent family (seeds 23, 322)**: triage under the derived
   extent (`lead + window + trailing`); the remove-segment extent producer is
   expected to reduce to the same derivation or a one-line fix in it.
5. **Cleanup handoff to step 14**: dead kernels out, docs/CollectionRectAnimator
   §9/§2.4 updated to the truth equation, memory updated; field-count audit
   (net ≤ −1).

## 6. Risks

- **Ballistic scroll during estimate drift.** Corrections during a ballistic
  restart simulations when they change dimensions (NEW-6's evidence:
  `applyNewDimensions → goBallistic` on every extent move). Mitigations are
  structural: the latching invariant (§2.1) means NO Δ exists without a
  measurement event; mid-band Δ defers under the existing velocity suppression;
  exactness Δ is one-shot by nature (lead can become exact once per tenure).
  Residual risk: a long backfill measuring one seam per frame during a slow
  backward drag emits a small Δ per frame — same cadence as today's telescoped
  seam corrections, so no regression, but O12 now watches it.
- **BouncingScrollPhysics.** iOS bounce + boundary Δ = spring restarts. The
  §3 rules (no estimate Δ while overscrolled; exactness Δ allowed) are designed
  for it, but the harness is clamping-only — the stage-0 bouncing batch is the
  cheap insurance; without it, this ships on reasoning alone.
- **maxScrollExtent oscillation.** Extent now has three terms; if EMA drift
  leaked into `lead` the extent would breathe at rest. The latching invariant
  forbids it; O12 enforces it. Trailing-estimate drift at rest is already
  impossible (EMA over a constant attached window converges).
- **Velocity-sensor interplay.** The epilogue MUST feed `refineLastOffset`'s
  correction compensation (as every emitter does today), or Δ reads as user
  motion, flips the suppression gate, and the gate then defers the next Δ — a
  latch/suppress oscillation. One line, but load-bearing; pin with a pure test.
- **`flying-flow!` interplay.** Flying bypasses the epilogue entirely (it
  writes geometry directly and never emits); on landing, reassignment picks the
  anchor from the attached window — the existing contract. Risk is only that
  stage 1 accidentally routes flying through the anchor resolve; keep its
  entry short-circuit explicit.
- **Translation cost.** The epilogue's translate is O(window + checkpoints) per
  emitting pass — the same order as today's `shift-attached-from!` +
  cache-rewrite paths, and emitting passes are rare (suppressed under motion).
  The materialization tripwire and the new O9 work probe bound regressions.

## 7. Open questions (with leanings)

1. **Keep a `leadExact` bool vs deriving exactness from the checkpoint store?**
   Leaning: keep the bool inside `leadEstimate` — one writer (epilogue +
   reassignment), debug-assertable, and the checkpoint store prunes too
   aggressively to be the source of record. The point is single ownership, not
   field elimination.
2. **`:settled` exactness ON (stage 3)** — is the mid-clock landing acceptable
   while dying cells play out? Leaning: yes — the domain-exited segment already
   runs resting drivers and the set-point is gone; deferring was conservatism,
   not contract. The new pinned test decides empirically; if a visible artifact
   shows up, fall back to "exactness Δ waits for `:resting`" (one capability
   row, not a design change).
3. **Δ hysteresis** — new threshold or today's `1e-3`? Leaning: `1e-3`; the
   latching invariant removes the storm source, so no new tunable. Revisit only
   if O12 flags residual breathing.
4. **Bouncing-physics batch in stage 0** — in or out? Leaning: in (one batch,
   ~10 s of campaign time, covers §3's two bounce rules); out only if the
   campaign budget objects.
5. **Does stage 2 fully close the o6 count-change family?** Leaning: yes — the
   drift mechanism (window moved by leading repairs while the anchor's offset
   is held) is structurally impossible once window layout is anchor-relative
   pre-translation. If seeds survive, the stage gate degrades gracefully
   (strict shrink + survivors become the next triage) rather than blocking.
6. **Should reassignment after a far jump preserve `:frac` when the same key
   lands in-window by luck?** Leaning: no — a far jump is user intent to leave
   the anchor; frac 0 at the inverted index is honest and matches today. Not
   worth the lookup.
