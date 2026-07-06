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

Caveat: ClojureDart emits unresolved member accesses as `(x as dynamic).m()`, so
those are NOT caught by `dart analyze` — the compiler reports them as
`DYNAMIC WARNING:` lines during the compile step instead. Watch both.
