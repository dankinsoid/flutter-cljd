# flutter-cljd dev tasks

# Device to run on (override: make repl DEVICE=chrome)
DEVICE ?= macos

.PHONY: repl check check-clean compile test

# Run the example app as a hot-reload ClojureDart REPL (main = repl.cljd).
# Must run from example/ — running from the repo root packages the library and fails.
repl:
	cd example && clojure -M:cljd flutter -d $(DEVICE) --enable-impeller

# Compile the example + dart analyze the generated Dart (fast, matches dev loop).
check:
	bin/check

# Same as check but wipes the ClojureDart cache first (fixes stale builds).
check-clean:
	bin/check --clean

# Compile the library (ClojureDart -> Dart) only.
compile:
	clojure -M:cljd compile

# Run the test suite.
test:
	clojure -M:cljd test
