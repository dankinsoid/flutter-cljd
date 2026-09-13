# Shader demo research

Research summary for a demo proving Flutter can deliver high-end GPU animation.
Collected 2026-09-12/13 by five parallel agents; key claims verified by hand.

Interactive catalogue of all 41 effects, with live links and filters:
<https://claude.ai/code/artifact/3ade0560-03c9-4a9b-9466-a7c2c5e1617a>

---

## 1. The central finding

Most of what we picked is not a set of separate effects. It is **one shader** with
several entry points, all reducible to the same chain:

```
height field → normal → what you do with the normal
```

| Entry | Height comes from | Normal is used to | Gives |
|---|---|---|---|
| SDF + `smin` | squircle profile from the shape border | — | shape, merging, morphing |
| glass | same profile | refract the backdrop | liquid glass |
| water | decaying rings from touch points | refract the backdrop | ripple across glass |
| metal | SDF blob | reflect a procedural matcap | liquid chrome |
| rubber | falloff around the grab point | specular | edge pull |
| membrane | sum of eigenmodes | refract + specular | drum, trampoline |
| paper | Voronoi cell edges | flat-shade the facets | crumple |

**Architectural consequence:** keep the height field in *screen* coordinates and shared
across the scene, never local to a single panel. Then "ripples inside an element" and
"a ripple travelling across the background into the glass" are the same architecture with
different parameters, not a rewrite.

**Where the system ends.** Cloth and irreversible crumpling need per-frame state — that is
`drawVertices` plus CPU physics, a separate path. Particles are also outside it: they
describe a set of objects, not a surface. The bridge is the **flow field** — one field
either drives particles or displaces surface UVs.

### Shared building blocks

Three pieces get reused across many effects, so build them once:

- **Mip pyramid** — progressive blur, glass frost, bloom. Also makes expensive things cheap:
  blur at half resolution costs a quarter.
- **FBM noise** — aurora, caustics, smoke reveal, dissolve, plus the flow field.
- **Flow field** — particle motion and surface displacement, same data.

---

## 2. Verified facts

Agent reports contained claims worth checking. These were opened and read directly.

| Subject | Claim | Verdict |
|---|---|---|
| iPhone Duo | Announced 2026-09-09 | **Confirmed** — [Apple Newsroom](https://www.apple.com/newsroom/2026/09/apple-unveils-iphone-duo/). 7.6"/5.4", matching aspect ratio, $1999, Oct 23 |
| Flutter Liquid Glass | "Good support announced yesterday" | **Not confirmed.** [flutter.dev/blog](https://flutter.dev/blog) has one September post: Material/Cupertino decoupling (Sep 9). Liquid Glass is *promised*, no API, no date |
| `ImageFilter.shader` | Backdrop reaches the shader on-GPU | **Confirmed** — [API docs](https://api.flutter.dev/flutter/dart-ui/ImageFilter/ImageFilter.shader.html). First uniform `vec2` = input size, first `sampler2D` filled by the engine. Impeller-only; Y axis flipped on GLES |
| `drawVertices` on Impeller | Broken on iOS (#127486) | **Closed**, fix landed. Mesh path is open |
| `liquid_glass_renderer` | Suggested as the reference | Correct SDFs and a cached geometry pass; last commit 2026-04, `0.2.0-dev.4`, 16-shape cap. Blur is a separate clipped layer while the shader displaces UVs, so edge samples read outside the clip → **own implementation planned** |
| `duo_motion` | Suggested as a 1:1 basis | Exists (Apache 2.0, four `.frag` files, no external users yet). Useful for its shader maths; its constants are unverified against the original |
| `riveo_page_curl` | Page curl works as a fragment shader | **Confirmed** — [340 stars](https://github.com/Rahiche/riveo_page_curl), `page_curl.frag` over a snapshot |
| gl-transitions | Named "liquid" transitions | **Partly.** `crosswarp`, `WaterDrop`, `Dreamy`, `Fold`, `Rolls`, `BookFlip` exist. `ripple` and `undulatingBurnOut` **do not** |
| `flutter_tearable_cloth` | Cloth already solved | Physics yes (80×60, Verlet, tearing, MIT). Rendering draws lines and points via CustomPainter, so a widget is not textured onto the mesh |

---

## 3. Platform constraints

**Fragment shaders only.** Vertex shaders live in `flutter_gpu`, preview on master.
Mesh deformation therefore goes through `Canvas.drawVertices` with a CPU grid.

**No shader composition.** Each `FragmentProgram` is a separately compiled artifact;
the engine cannot merge two programs into one pass. Layers on the same resolution must be
merged into **one fat shader** with branches, not stacked. Only passes that genuinely change
resolution (the mip pyramid) stay separate. `BackdropGroup` merges backdrop *captures*, not
shader code.

**Uniform limits.** Arrays are indexed by constants only; Impeller's uniform buffer is small
— `liquid_glass_renderer` cut shapes from 64 to 16 for this reason. This repo's
`shaders/mesh_gradient.frag` already solves it via a texture path (`uUseTexture`, data in an
`Nx2` texture) — reuse that for SDF shapes.

**Cost is reads, not area.** GPUs parallelise per pixel, but every pixel's texture reads cost
bandwidth, and scattered reads (refraction) miss cache.

| Effect | Reads/pixel | Full-screen? |
|---|---|---|
| doppler, aberration, dither | 1–3 | yes |
| holographic, aurora | 0–1 (procedural) | yes |
| refraction without blur | 1–3 | yes |
| glass with frost blur | tens | local only |
| fluid | tens × ~30 passes | no |

This is why iOS Liquid Glass stutters on large components: area × reads × poor locality,
and `BackdropFilter` processes the whole screen regardless of panel size.

**Snapshots in transitions.** `AnimatedSampler` rebuilds the texture every frame, so a "live"
widget is not protection against a mid-animation jump — it *guarantees* one. Transitions
should **freeze** snapshots A and B. Apple's fold does exactly this: the live UI is never
rendered mid-transition (the two states have different layouts anyway).

---

## 4. The two iPhone Duo targets

### Fold

Content stays "locked in space" while the phone reorients around the software.

1. Two half-planes around the hinge, true pinhole perspective — not a 2.5D skew
2. UI projected in **screen** coordinates, not attached to the rotating geometry — this is the whole illusion
3. Progressive defocus: sharp at the hinge, blurred toward the free edge
4. Cross-fade between two snapshots
5. Inner-corner darkening plus a specular sheen sliding across the glass

**Progress is an input parameter, not a duration** — driven by hinge angle, settled with a
spring. In Flutter: `SpringSimulation`.

Deformation is analytic, so this is a fragment shader over a snapshot — no mesh. Note the
trap: `drawVertices` cannot be combined with a *custom* fragment shader, since interpolated
UVs never reach the runtime shader ([#173835](https://github.com/flutter/flutter/issues/173835)).

### Tab bar

The gradient-blur hunch was right, and Apple says it verbatim: the iOS 26 scroll edge effect
is "variable blur combined with a soft gradient".

1. Scroll edge effect: σ grows toward the screen edge over a ~80–120 lpx gradient mask
2. Liquid Glass pill on top: refraction via SDF capsule normals
3. Specular rim, tint, saturation boost
4. Collapse on scroll, tied to scroll offset rather than time

Two independent shader layers — but the blur's mip chain feeds the glass frost.

**Why naive blur fails:** `blur + alpha mask` is a cross-fade between sharp and blurred, so
high frequencies survive as ghosting. A genuine varying radius is required — mip pyramid with
lod-lerp, or a two-pass variable gaussian. The CSSWG
[explicitly rejects](https://github.com/w3c/csswg-drafts/issues/13285) the mask workaround as
equivalent.

**Apple's glass profile is a squircle**, not a sphere: `h(x) = (1-(1-x)^4)^(1/4)`. It keeps the
refraction gradient smooth when stretched into a rectangle. This repo's `smooth_corner.cljd`
(810 lines of squircle geometry) is the foundation.

---

## 5. Priorities

### High — the core

1. **Foundation** — three shader primitives: self-contained, over a child snapshot, over the
   backdrop (with a non-Impeller fallback). Today `shader.cljd` is 55 lines hard-coding two
   specific shaders.
2. **Holographic card** — accelerometer-driven. Zero infrastructure, best wow-to-effort ratio,
   and the only effect that *requires* a phone in hand — you cannot fake it with a screenshot.
   Good pipeline shakedown before the backdrop machinery.
3. **Tab bar** — progressive blur + glass pill. Also the skeleton of water UI.
4. **Water UI** — glass + water + metaballs in one height field.
5. **Fold** — the iPhone Duo target.

Sensor pipeline for the holographic card matters more than the shader: One Euro filter
(adaptive — smooths at rest, stays responsive in motion), gravity vector separate from linear
acceleration, and a spring on the input so the card trails the tilt with slight inertia. Glass
specular reuses the same light vector.

### Medium

- **Fluted glass** — first in this group: one height function on top of working glass gives a
  visually distinct second material, proving the architecture. Watch for moire — sample from a
  mip level based on rib width in screen pixels.
- **Liquid metal** — reflect a matcap instead of refracting. An evening on top of glass.
- **Rubber-band pull (2D)** — analytic; does not exist as a Flutter package.
- **Disintegration** — one shader, three modes: Thanos (hash-based shards), smoke reveal
  (FBM threshold), blow-away (flow source at the gesture point). All are
  `snapshot sample displaced along a field + mask`.
- **Genie** — analytic, like curl. The two-phase timing is what makes it read as suction:
  the leading edge leaves before the trailing one.
- **Displacement transitions** — frozen A/B snapshots.
- **God rays behind text** — text is an ideal occluder (many thin gaps), and with live input the
  rays animate themselves. **Not for glass**: rays mean *blocked* light, which contradicts
  transparency. Glass glows via caustics or edge bloom instead.
- **Caustics** — divergence of the existing normal field, no extra samples. Best use: a glass
  element's shadow that focuses light instead of being a blurred rectangle. More convincing
  from water than from static glass, because it moves.
- **Staggered spring scroll** — the only non-shader item, see §6.
- **Liquid wave transitions** — port `crosswarp` / `WaterDrop` from gl-transitions.

### Open / flagship

- **Cloth** — take the physics from `flutter_tearable_cloth`, write the rendering: `drawVertices`
  with snapshot UVs, plus light baked into vertex colors (`drawVertices` has no lighting, and
  without it drapery reads flat — folds show only as texture distortion, never as shading).
- **GPU fluid** — the only effect needing true ping-pong. Applying it to *widgets* appears to be
  unexplored. Try **advected noise** first (analytic velocity field, one pass); only go to real
  Navier–Stokes if that falls short. Bloom pairs well: it reveals flow structure, and reuses the
  mip pyramid.
- **Real GPU particles** — blocked on missing vertex shaders. A project, not an evening.

---

## 6. Repository split

**Separate repo:** shaders and materials. They are self-contained — a `.frag` plus a thin
binding — and this keeps them portable for a possible contribution to
[fluttershaders.com](https://fluttershaders.com/) (9 shaders there, 7 of them post-process
filters over an image; glass, materials and composed height fields are not covered). Keep all maths in the `.frag` and only
parameter passing in the binding, so a shader detaches as one file. That site's shaders are deliberately isolated,
one technique each — our strength is composition, which would need breaking apart before it
made a useful contribution there.

**This repo:** staggered spring scroll, and bottom-sheet stretch as its 1D case. It wedges
between the layout's target position and the painted position, so it belongs inside the
collection engine and cannot be extracted. Watch the existing rule — inner fixed duration
outranks outer — so per-index delay must not stretch each element's spring, or the last cards
move slowly instead of moving later.

Existing assets worth reusing: `shaders/mesh_gradient.frag` (Coons patch, Newton inversion,
Catmull-Rom, OKLab, dithering — and the texture path around the uniform limit) makes the best
possible background for glass, since refraction is only visible over detail;
`smooth_corner.cljd` for squircle SDF; `animations.cljd` and `curves.cljd` for springs.

---

## 7. Reference notes

**Shape morphing.** The user's own approach in VDAnimation — resample both contours to equal
vertex counts, match by polar angle around the centroid — is a named technique (angular
correspondence). It is correct exactly when both shapes are **star-shaped** about the centroid;
it breaks on crescents and horseshoes, where a ray crosses the contour three times, or where the
centroid falls outside the shape. Two cheap patches: pick the cyclic shift minimising total pair
distance (otherwise rotation always goes one way), and resample by arc length rather than angle.
Versus SDF interpolation (`mix(sdA, sdB, t)`), which needs no correspondence at all but may split
a shape in two: the contour method guarantees topology, which matters more for UI shapes.

**Noise families.** FBM (summed octaves) drives fire, clouds, aurora — the noise *is* the pixel.
Minecraft-style terrain uses Perlin to generate *geometry*, computed once and stored, not per
frame. Voronoi/Worley is different in kind — distance to nearest point, giving sharp cell edges;
that is what makes crumple read as paper. FBM instead of Voronoi makes it read as cloth: the
classic mistake.

**Chromatic aberration vs doppler.** Aberration dephases RGB channels, so it colours *edges* and
leaves flat fills grey. A doppler shift moves the whole spectrum, colouring everything. The
latter is the user's idea and is not in the catalogue anywhere — scroll up reddens, scroll down
blues, stronger at the edges than at centre. Cheap (no extra samples), though on photographic
content it may read as a broken screen.
