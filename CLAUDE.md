# ClojureDart Coding Rules

## Type Annotations
Always annotate Dart interop types for efficient compiled code:
- Function return types: `(defn ^w/Widget my-fn ...)`
- deftype fields: `^w/Key? k`, `^:mutable ^int count`
- Local bindings with Dart casts: `(let [^Button widget (.-widget this)] ...)`
- Generic types: `^#/(w/State MyWidget)`, `^#/(a/Tween double)`
- Optional types use `?`: `^w/Widget?`, `^dc/Duration?`

## Build & Verify
`clojure -M:cljd compile` only translates ClojureDart → Dart; it does not check
that the emitted Dart type-checks. To also statically analyze the generated Dart,
run:

```sh
bin/check
```

It compiles, then runs `dart analyze lib/cljd-out`, failing only on
error-severity issues (warnings/infos are the generated-code baseline).

Scope: `dart analyze` catches malformed emitted Dart, broken imports, and type
mismatches on statically-resolved targets. It does NOT catch calls to unresolved
members — ClojureDart emits those as `(x as dynamic).m()`. Dynamic dispatch is
idiomatic here (~240 `DYNAMIC WARNING:` lines on a clean build), so those
warnings are baseline noise, not an error gate — don't try to read them.
