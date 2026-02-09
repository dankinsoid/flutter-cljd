#version 460 core

#include <flutter/runtime_effect.glsl>

precision highp float;

// --- Uniforms ---
// 0-1
uniform vec2 uSize;
// 2
uniform float uType;      // 0=linear, 1=radial, 2=sweep
// 3
uniform float uTileMode;  // 0=clamp, 1=repeat, 2=mirror, 3=decal
// 4
uniform float uColorCount;
// 5
uniform float uUseTexture; // 0=uniforms, 1=texture
// 6-9
uniform vec2 uParam0;     // linear: begin; radial/sweep: center
uniform vec2 uParam1;     // linear: end; radial: (radius, _); sweep: (startAngle, endAngle)
// 10-137
uniform vec4 uColors[32];
// 138-169
uniform float uStops[32];

uniform sampler2D uData;   // Nx2 data texture (row 0: RGB+A=255, row 1: stop+alpha+A=255)

out vec4 fragColor;

// --- sRGB <-> Linear ---

float srgbToLinear(float c) {
    return c > 0.04045
        ? pow((c + 0.055) / 1.055, 2.4)
        : c / 12.92;
}

float linearToSrgb(float c) {
    return c > 0.0031308
        ? 1.055 * pow(c, 1.0 / 2.4) - 0.055
        : 12.92 * c;
}

vec3 srgbToLinearV(vec3 c) {
    return vec3(srgbToLinear(c.r), srgbToLinear(c.g), srgbToLinear(c.b));
}

vec3 linearToSrgbV(vec3 c) {
    return vec3(linearToSrgb(c.r), linearToSrgb(c.g), linearToSrgb(c.b));
}

// --- Linear RGB <-> OKLab ---

vec3 linearToOklab(vec3 rgb) {
    float l = pow(max(0.4122214708 * rgb.r + 0.5363325363 * rgb.g + 0.0514459929 * rgb.b, 0.0), 1.0 / 3.0);
    float m = pow(max(0.2119034982 * rgb.r + 0.6806995451 * rgb.g + 0.1073969566 * rgb.b, 0.0), 1.0 / 3.0);
    float s = pow(max(0.0883024619 * rgb.r + 0.2817188376 * rgb.g + 0.6299787005 * rgb.b, 0.0), 1.0 / 3.0);
    return vec3(
        0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
        1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
        0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s
    );
}

vec3 oklabToLinear(vec3 lab) {
    float l_ = lab.x + 0.3963377774 * lab.y + 0.2158037573 * lab.z;
    float m_ = lab.x - 0.1055613458 * lab.y - 0.0638541728 * lab.z;
    float s_ = lab.x - 0.0894841775 * lab.y - 1.2914855480 * lab.z;
    float l = l_ * l_ * l_;
    float m = m_ * m_ * m_;
    float s = s_ * s_ * s_;
    return vec3(
         4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
        -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
        -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s
    );
}

// --- Data access (uniform or texture) ---

vec4 getColor(int i) {
    if (uUseTexture < 0.5) {
        return uColors[i];
    }
    float u = (float(i) + 0.5) / uColorCount;
    // Row 0: RGB color (A=255 in texture to avoid premul issues)
    vec3 rgb = texture(uData, vec2(u, 0.25)).rgb;
    // Row 1: R=stop, G=alpha
    float alpha = texture(uData, vec2(u, 0.75)).g;
    return vec4(rgb, alpha);
}

float getStop(int i) {
    if (uUseTexture < 0.5) {
        return uStops[i];
    }
    float u = (float(i) + 0.5) / uColorCount;
    return texture(uData, vec2(u, 0.75)).r;
}

// --- Gradient position ---

float computeT(vec2 pos) {
    if (uType < 0.5) {
        // Linear: project onto begin→end line
        vec2 delta = uParam1 - uParam0;
        float len2 = dot(delta, delta);
        if (len2 < 0.0001) return 0.0;
        return dot(pos - uParam0, delta) / len2;
    } else if (uType < 1.5) {
        // Radial: distance from center / radius
        float radius = uParam1.x;
        if (radius < 0.0001) return 0.0;
        return length(pos - uParam0) / radius;
    } else {
        // Sweep: angle from center, normalized to [0,1]
        vec2 d = pos - uParam0;
        float angle = atan(d.y, d.x);
        if (angle < 0.0) angle += 6.28318530718;
        float startAngle = uParam1.x;
        float endAngle = uParam1.y;
        float range = endAngle - startAngle;
        if (abs(range) < 0.0001) return 0.0;
        return (angle - startAngle) / range;
    }
}

// --- Tile mode ---

float applyTileMode(float t) {
    if (uTileMode < 0.5) {
        // Clamp
        return clamp(t, 0.0, 1.0);
    } else if (uTileMode < 1.5) {
        // Repeat
        return fract(t);
    } else if (uTileMode < 2.5) {
        // Mirror
        float mt = mod(t, 2.0);
        return mt > 1.0 ? 2.0 - mt : mt;
    } else {
        // Decal: -1 signals "outside"
        return (t < 0.0 || t > 1.0) ? -1.0 : t;
    }
}

// --- Main ---

void main() {
    vec2 pos = FlutterFragCoord().xy;
    float t = computeT(pos);
    t = applyTileMode(t);

    // Decal: transparent outside gradient
    if (t < 0.0) {
        fragColor = vec4(0.0);
        return;
    }

    int count = int(uColorCount);

    // Find the segment containing t
    int idx = 0;
    for (int i = 0; i < count - 1; i++) {
        if (t >= getStop(i)) idx = i;
    }

    // Local interpolation factor within the segment
    float sA = getStop(idx);
    float sB = getStop(idx + 1);
    float range = sB - sA;
    float localT = range > 0.0001 ? clamp((t - sA) / range, 0.0, 1.0) : 0.0;

    // Source colors (sRGB 0-1)
    vec4 cA = getColor(idx);
    vec4 cB = getColor(idx + 1);

    // sRGB → linear → OKLab
    vec3 labA = linearToOklab(srgbToLinearV(cA.rgb));
    vec3 labB = linearToOklab(srgbToLinearV(cB.rgb));

    // Interpolate in OKLab
    vec3 labMix = mix(labA, labB, localT);
    float alphaMix = mix(cA.a, cB.a, localT);

    // OKLab → linear → sRGB, clamped
    vec3 rgb = clamp(linearToSrgbV(oklabToLinear(labMix)), 0.0, 1.0);

    // Premultiplied alpha
    fragColor = vec4(rgb * alphaMix, alphaMix);
}
