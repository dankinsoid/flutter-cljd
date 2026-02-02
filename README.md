# flutter-cljd

![Alpha Status](https://img.shields.io/badge/status-alpha-red)

ClojureDart wrapper for Flutter Material widgets, designed to simplify and compact Flutter development in ClojureDart. It provides concise, Clojure-like syntax to work with Flutter’s Material components and types, making code more readable and expressive for Clojure developers building Flutter apps.

## Main Goals

### Provide a more concise and readable UI syntax
The library focuses on simplifying the syntax for building Flutter UIs, making it more compact, intuitive, and aligned with Clojure’s functional style.

### Use Clojure data structures for better consistency and flexibility
The API is designed around pure Clojure types instead of Dart’s, offering a more seamless and consistent experience for Clojure developers while increasing code flexibility.

### Streamline and enhance Dart APIs
The library simplifies certain Dart APIs, making them easier to use and more expressive, improving the overall developer experience.

## Extensions
Besides wrappers the library contains some extensions to Flutter
- [Animations](./docs/Animations.md)
- [Button](./docs/Button.md)

## Examples

```clojure
;; Basic button with styling
(->> (text "Click me!")
     (with-text-style :color :blue, :size 16)
     (padding :h 16 :v 8)
     (button #(println "Clicked!")))
```
```clojure
;; Card with multiple elements
(->> (row :spacing 10
       (text "Title" :size 20, :weight :bold)
       (text "Subtitle" :color :grey))
     (padding 16)
     (card :elevation 2 :radius 8))
```

### Styling and Layout
```clojure
;; Applying styles and layouts
(->> (text "Styled Text")
     (with-text-style :color :blue
                      :size 20
                      :weight :bold)
     (padding 16)
     (center))

;; Responsive layouts
(->> (column :spacing 8
       (text "Header")
       (expanded
         (list-view
           (for [i (range 10)]
             (text (str "Item " i)))))
       (text "Footer"))
     (container :color :white))
```

### Interactive Components
```clojure
;; Button with feedback
(->> (text "Click Me!")
     (button #(println "Clicked!")
             {:on-long-press #(println "Long pressed!")
              :on-hover #(println "Hover: " %)}))

;; Form elements
(let [controller (atom "")]
  (->> (text-field
         {:controller controller
          :on-changed #(reset! controller %)
          :decoration {:label "Enter text"}})
       (padding 16)))
```

### Complete Example
```clojure
(ns readme.example
  (:require [flutter-cljd.widgets :as ui]))

;; User profile card
(defn profile-card [& {:keys [name role avatar]}]
  (->> (ui/row
         ;; Avatar section
         (->> avatar
              (ui/circle :size 40)
              (ui/padding :right 12))
         
         ;; Text content
         (ui/column :spacing 10
           (ui/text name :size 18, :weight :bold)
           (ui/text role :size 14)))
       (ui/with-text-style :color :grey)
       (ui/padding 16)
       (ui/card :elevation 4 :radius 12)
       (ui/center)))

;; Usage example
(def user-profile
  (profile-card
     :name "John Doe"
     :role "Senior Developer"
     :avatar (ui/image "path/to/avatar.png")))
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

### Guidelines

- Follow the existing code style
- Add tests for new features
- Update documentation as needed
- Keep commits focused and atomic
- Write clear commit messages (I recommend Aider to generate commit messages)

## License

This project is licensed under the MIT License - see the LICENSE file for details.
