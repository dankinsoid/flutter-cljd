#version 460 core

#include <flutter/runtime_effect.glsl>

precision highp float;

// --- Uniforms ---
// 0-1
uniform vec2 uSize;
// 2
uniform float uGridWidth;
// 3
uniform float uGridHeight;
// 4
uniform float uUseTexture; // 0=uniforms, 1=texture
// 5-132
uniform vec2 uPositions[64];
// 133-388
uniform vec4 uColors[64];

uniform sampler2D uData;   // Nx2 data texture (row 0: pos x,y + alpha, row 1: RGB)

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

vec2 getPosition(int i) {
    if (uUseTexture < 0.5) {
        return uPositions[i];
    }
    float n = uGridWidth * uGridHeight;
    float u = (float(i) + 0.5) / n;
    // Row 0: R=x, G=y, B=alpha (8-bit each)
    return texture(uData, vec2(u, 0.25)).rg;
}

vec4 getMeshColor(int i) {
    if (uUseTexture < 0.5) {
        return uColors[i];
    }
    float n = uGridWidth * uGridHeight;
    float u = (float(i) + 0.5) / n;
    // Row 0 B channel = alpha
    float alpha = texture(uData, vec2(u, 0.25)).b;
    // Row 1: RGB color
    vec3 rgb = texture(uData, vec2(u, 0.75)).rgb;
    return vec4(rgb, alpha);
}

// --- Inverse bilinear interpolation ---
// Based on Inigo Quilez's method.
// Given a point p and a quad (a, b, c, d) where:
//   a = top-left, b = top-right, c = bottom-right, d = bottom-left
// Returns (u, v) in [0,1] if p is inside the quad, else (-1, -1).

float cross2d(vec2 a, vec2 b) {
    return a.x * b.y - a.y * b.x;
}

vec2 invBilinear(vec2 p, vec2 a, vec2 b, vec2 c, vec2 d) {
    vec2 e = b - a;
    vec2 f = d - a;
    vec2 g = a - b + c - d;
    vec2 h = p - a;

    float k2 = cross2d(g, f);
    float k1 = cross2d(e, f) + cross2d(h, g);
    float k0 = cross2d(h, e);

    float u, v;

    if (abs(k2) < 0.0001) {
        // Near-linear case (parallel edges)
        if (abs(k1) < 0.0001) return vec2(-1.0);
        v = -k0 / k1;
        float dx = e.x + g.x * v;
        float dy = e.y + g.y * v;
        u = (abs(dx) > abs(dy))
            ? (h.x - f.x * v) / dx
            : (h.y - f.y * v) / dy;
    } else {
        float w = k1 * k1 - 4.0 * k0 * k2;
        if (w < 0.0) return vec2(-1.0);
        w = sqrt(w);

        v = (-k1 - w) / (2.0 * k2);
        float dx = e.x + g.x * v;
        float dy = e.y + g.y * v;
        u = (abs(dx) > abs(dy))
            ? (h.x - f.x * v) / dx
            : (h.y - f.y * v) / dy;

        if (u < -0.001 || u > 1.001 || v < -0.001 || v > 1.001) {
            v = (-k1 + w) / (2.0 * k2);
            dx = e.x + g.x * v;
            dy = e.y + g.y * v;
            u = (abs(dx) > abs(dy))
                ? (h.x - f.x * v) / dx
                : (h.y - f.y * v) / dy;
        }
    }

    if (u < -0.001 || u > 1.001 || v < -0.001 || v > 1.001)
        return vec2(-1.0);

    return clamp(vec2(u, v), 0.0, 1.0);
}

// --- Main ---

void main() {
    vec2 pos = FlutterFragCoord().xy / uSize;
    int gw = int(uGridWidth);
    int gh = int(uGridHeight);

    for (int row = 0; row < gh - 1; row++) {
        for (int col = 0; col < gw - 1; col++) {
            int i00 = row * gw + col;
            int i10 = i00 + 1;
            int i01 = i00 + gw;
            int i11 = i01 + 1;

            vec2 p00 = getPosition(i00);
            vec2 p10 = getPosition(i10);
            vec2 p01 = getPosition(i01);
            vec2 p11 = getPosition(i11);

            // AABB early reject
            vec2 lo = min(min(p00, p10), min(p01, p11));
            vec2 hi = max(max(p00, p10), max(p01, p11));
            if (pos.x < lo.x || pos.x > hi.x || pos.y < lo.y || pos.y > hi.y)
                continue;

            // a=top-left, b=top-right, c=bottom-right, d=bottom-left
            vec2 uv = invBilinear(pos, p00, p10, p11, p01);
            if (uv.x < 0.0) continue;

            float u = uv.x;
            float v = uv.y;

            // Fetch corner colors
            vec4 c00 = getMeshColor(i00);
            vec4 c10 = getMeshColor(i10);
            vec4 c01 = getMeshColor(i01);
            vec4 c11 = getMeshColor(i11);

            // sRGB -> linear -> OKLab
            vec3 lab00 = linearToOklab(srgbToLinearV(c00.rgb));
            vec3 lab10 = linearToOklab(srgbToLinearV(c10.rgb));
            vec3 lab01 = linearToOklab(srgbToLinearV(c01.rgb));
            vec3 lab11 = linearToOklab(srgbToLinearV(c11.rgb));

            // Bilinear interpolation in OKLab
            vec3 labTop = mix(lab00, lab10, u);
            vec3 labBot = mix(lab01, lab11, u);
            vec3 labFinal = mix(labTop, labBot, v);

            // Alpha interpolation
            float aTop = mix(c00.a, c10.a, u);
            float aBot = mix(c01.a, c11.a, u);
            float alpha = mix(aTop, aBot, v);

            // OKLab -> linear -> sRGB, clamped
            vec3 rgb = clamp(linearToSrgbV(oklabToLinear(labFinal)), 0.0, 1.0);

            // Premultiplied alpha
            fragColor = vec4(rgb * alpha, alpha);
            return;
        }
    }

    // Outside all cells
    fragColor = vec4(0.0);
}
