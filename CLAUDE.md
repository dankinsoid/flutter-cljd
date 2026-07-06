# ClojureDart Coding Rules

## Type Annotations
Always annotate Dart interop types for efficient compiled code:
- Function return types: `(defn ^w/Widget my-fn ...)`
- deftype fields: `^w/Key? k`, `^:mutable ^int count`
- Local bindings with Dart casts: `(let [^Button widget (.-widget this)] ...)`
- Generic types: `^#/(w/State MyWidget)`, `^#/(a/Tween double)`
- Optional types use `?`: `^w/Widget?`, `^dc/Duration?`

## Build & Verify
`clojure -M:cljd compile` only translates ClojureDart → Dart; it never runs the
Dart front-end. An internally-inconsistent build still reports success and only
fails when `flutter run` compiles it at launch. To catch that ahead of time:

```sh
bin/check           # compile example + dart analyze (fast, matches dev loop)
bin/check --clean   # wipe the ClojureDart cache first
```

It compiles the **example app** (which pulls the library in as a dep, so library
errors surface too) and runs `dart analyze` on the generated Dart, failing only
on error-severity issues. Warnings/infos are the generated-code baseline.

Most common launch failure: `Type '...IFnMixin_...' not found` — a stale
incremental build where hashed mixin references got out of sync between files.
`dart analyze` reports it as `mixin_of_non_class`; the fix is `bin/check --clean`
(or `rm -rf example/.clojuredart example/lib/cljd-out` then recompile).

What `dart analyze` does NOT catch: calls to members ClojureDart couldn't resolve
— those are emitted as `(x as dynamic).m()` and slip past both the analyzer and
the real compiler, failing only at runtime. Dynamic dispatch is idiomatic here
(~240 `DYNAMIC WARNING:` lines on a clean build), so those compiler warnings are
baseline noise, not an error gate — don't try to read them.
