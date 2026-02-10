#version 460 core

#include <flutter/runtime_effect.glsl>

precision highp float;

// --- Uniforms ---
// 0-1
uniform vec2 uSize;
// 2-3
uniform vec2 uOffset;     // rect.topLeft — FragCoord is absolute canvas coords
// 4
uniform float uGridWidth;
// 5
uniform float uGridHeight;
// 6
uniform float uUseTexture; // 0=uniforms, 1=texture
// 7-134
uniform vec2 uPositions[64];
// 135-390
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

// --- Catmull-Rom bicubic interpolation ---
// Weights for t: samples at positions -1, 0, 1, 2

vec4 catmullRomWeights(float t) {
    float t2 = t * t;
    float t3 = t2 * t;
    return vec4(
        -0.5*t3 + t2 - 0.5*t,
         1.5*t3 - 2.5*t2 + 1.0,
        -1.5*t3 + 2.0*t2 + 0.5*t,
         0.5*t3 - 0.5*t2
    );
}

// Get OKLab + alpha at grid position, clamped to grid bounds
vec4 gridColorLab(int row, int col, int gw, int gh) {
    int r = clamp(row, 0, gh - 1);
    int c = clamp(col, 0, gw - 1);
    vec4 srgb = getMeshColor(r * gw + c);
    return vec4(linearToOklab(srgbToLinearV(srgb.rgb)), srgb.a);
}

// --- Coons patch (cubic Bezier boundaries) ---

vec2 gridPos(int row, int col, int gw, int gh) {
    return getPosition(clamp(row, 0, gh - 1) * gw + clamp(col, 0, gw - 1));
}

// Evaluate cubic Bezier B(t) and derivative B'(t)
vec2 cubicBezier(vec2 B0, vec2 B1, vec2 B2, vec2 B3, float t, out vec2 dBdt) {
    float s = 1.0 - t;
    float s2 = s * s;
    float t2 = t * t;
    dBdt = 3.0 * (s2*(B1-B0) + 2.0*s*t*(B2-B1) + t2*(B3-B2));
    return s2*s*B0 + 3.0*s2*t*B1 + 3.0*s*t2*B2 + t2*t*B3;
}

// Evaluate Coons patch S(u,v) and Jacobian at cell (row,col)
// Boundary curves are cubic Bezier with tangents from Catmull-Rom
vec2 evalCoons(float u, float v, int row, int col, int gw, int gh,
               out vec2 dSdu, out vec2 dSdv) {
    // Corners
    vec2 c00 = gridPos(row, col, gw, gh);
    vec2 c10 = gridPos(row, col+1, gw, gh);
    vec2 c01 = gridPos(row+1, col, gw, gh);
    vec2 c11 = gridPos(row+1, col+1, gw, gh);

    // Bezier control points from Catmull-Rom tangents: offset = (next - prev) / 6
    // Top edge (horizontal, row)
    vec2 tB1 = c00 + (c10 - gridPos(row, col-1, gw, gh)) / 6.0;
    vec2 tB2 = c10 - (gridPos(row, col+2, gw, gh) - c00) / 6.0;
    // Bottom edge (horizontal, row+1)
    vec2 bB1 = c01 + (c11 - gridPos(row+1, col-1, gw, gh)) / 6.0;
    vec2 bB2 = c11 - (gridPos(row+1, col+2, gw, gh) - c01) / 6.0;
    // Left edge (vertical, col)
    vec2 lB1 = c00 + (c01 - gridPos(row-1, col, gw, gh)) / 6.0;
    vec2 lB2 = c01 - (gridPos(row+2, col, gw, gh) - c00) / 6.0;
    // Right edge (vertical, col+1)
    vec2 rB1 = c10 + (c11 - gridPos(row-1, col+1, gw, gh)) / 6.0;
    vec2 rB2 = c11 - (gridPos(row+2, col+1, gw, gh) - c10) / 6.0;

    // Evaluate boundary curves
    vec2 topD, botD, leftD, rightD;
    vec2 top   = cubicBezier(c00, tB1, tB2, c10, u, topD);
    vec2 bot   = cubicBezier(c01, bB1, bB2, c11, u, botD);
    vec2 left  = cubicBezier(c00, lB1, lB2, c01, v, leftD);
    vec2 right = cubicBezier(c10, rB1, rB2, c11, v, rightD);

    // Coons patch = ruled columns + ruled rows - bilinear corners
    float su = 1.0 - u, sv = 1.0 - v;
    vec2 bilin = su*sv*c00 + u*sv*c10 + su*v*c01 + u*v*c11;
    vec2 S = sv*top + v*bot + su*left + u*right - bilin;

    // Jacobian
    dSdu = sv*topD + v*botD - left + right - (sv*(c10-c00) + v*(c11-c01));
    dSdv = -top + bot + su*leftD + u*rightD - (su*(c01-c00) + u*(c11-c10));

    return S;
}

// --- Main ---

void main() {
    vec2 pos = (FlutterFragCoord().xy - uOffset) / uSize;
    int gw = int(uGridWidth);
    int gh = int(uGridHeight);

    for (int row = 0; row < 15; row++) {
        if (row >= gh - 1) break;
        for (int col = 0; col < 15; col++) {
            if (col >= gw - 1) break;

            // Corner positions (reused for AABB and initial guess)
            vec2 c00 = getPosition(row * gw + col);
            vec2 c10 = getPosition(row * gw + col + 1);
            vec2 c01 = getPosition((row + 1) * gw + col);
            vec2 c11 = getPosition((row + 1) * gw + col + 1);

            // AABB with margin for Bezier overshoot
            vec2 lo = min(min(c00, c10), min(c01, c11));
            vec2 hi = max(max(c00, c10), max(c01, c11));
            vec2 margin = max(hi - lo, vec2(0.05));
            lo -= margin;
            hi += margin;
            if (pos.x < lo.x || pos.x > hi.x || pos.y < lo.y || pos.y > hi.y)
                continue;

            // Initial guess: linear projection onto cell parallelogram
            vec2 d = pos - c00;
            vec2 ex = c10 - c00;
            vec2 ey = c01 - c00;
            float det0 = ex.x * ey.y - ex.y * ey.x;
            float u = abs(det0) > 1e-6
                ? clamp((d.x * ey.y - d.y * ey.x) / det0, 0.0, 1.0) : 0.5;
            float v = abs(det0) > 1e-6
                ? clamp((ex.x * d.y - ex.y * d.x) / det0, 0.0, 1.0) : 0.5;

            // Newton-Raphson on Coons patch
            // Multiple starting points for robustness with non-monotonic grids
            bool accepted = false;
            float bestDist = 1e6;
            float bestU = 0.5, bestV = 0.5;

            for (int attempt = 0; attempt < 4; attempt++) {
                if (attempt == 1) { u = 0.5; v = 0.5; }
                else if (attempt == 2) { u = 0.25; v = 0.25; }
                else if (attempt == 3) { u = 0.75; v = 0.75; }

                for (int iter = 0; iter < 10; iter++) {
                    vec2 dSdu, dSdv;
                    vec2 S = evalCoons(u, v, row, col, gw, gh, dSdu, dSdv);
                    vec2 err = pos - S;

                    if (dot(err, err) < 1e-12) break;

                    float det = dSdu.x * dSdv.y - dSdu.y * dSdv.x;
                    if (abs(det) < 1e-12) break;

                    float du = ( err.x * dSdv.y - err.y * dSdv.x) / det;
                    float dv = (-err.x * dSdu.y + err.y * dSdu.x) / det;

                    // Damp large steps to prevent divergence near singular Jacobian
                    float stepLen = max(abs(du), abs(dv));
                    if (stepLen > 1.0) { du /= stepLen; dv /= stepLen; }

                    u += du;
                    v += dv;
                }

                float cu = clamp(u, 0.0, 1.0);
                float cv = clamp(v, 0.0, 1.0);
                vec2 dSdu_c, dSdv_c;
                vec2 Sc = evalCoons(cu, cv, row, col, gw, gh, dSdu_c, dSdv_c);
                float dist = length(pos - Sc);

                if (dist < bestDist) {
                    bestDist = dist;
                    bestU = cu;
                    bestV = cv;
                }

                if (dist < 0.002) {
                    u = cu; v = cv;
                    accepted = true;
                    break;
                }
            }

            // Fallback: use best solution if reasonably close
            if (!accepted && bestDist < 0.008) {
                u = bestU; v = bestV;
                accepted = true;
            }

            if (!accepted) continue;

            // Bicubic Catmull-Rom color interpolation in OKLab
            vec4 wu = catmullRomWeights(u);
            vec4 wv = catmullRomWeights(v);

            vec3 labFinal = vec3(0.0);
            float alpha = 0.0;

            for (int j = 0; j < 4; j++) {
                vec3 rowLab = vec3(0.0);
                float rowA = 0.0;
                for (int i = 0; i < 4; i++) {
                    vec4 labA = gridColorLab(row + j - 1, col + i - 1, gw, gh);
                    rowLab += labA.xyz * wu[i];
                    rowA += labA.w * wu[i];
                }
                labFinal += rowLab * wv[j];
                alpha += rowA * wv[j];
            }

            alpha = clamp(alpha, 0.0, 1.0);
            vec3 rgb = clamp(linearToSrgbV(oklabToLinear(labFinal)), 0.0, 1.0);

            fragColor = vec4(rgb * alpha, alpha);
            return;
        }
    }

    fragColor = vec4(0.0);
}
