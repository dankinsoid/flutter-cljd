# ClojureDart Coding Rules

## Type Annotations
Always annotate Dart interop types for efficient compiled code:
- Function return types: `(defn ^w/Widget my-fn ...)`
- deftype fields: `^w/Key? k`, `^:mutable ^int count`
- Local bindings with Dart casts: `(let [^Button widget (.-widget this)] ...)`
- Generic types: `^#/(w/State MyWidget)`, `^#/(a/Tween double)`
- Optional types use `?`: `^w/Widget?`, `^dc/Duration?`
