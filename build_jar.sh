#!/usr/bin/env bash
# Build a self-contained, runnable PepQueryMHC jar and name it after the
# version declared in Constants.java (e.g. PepQueryMHC.v1.0.7.jar).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CONSTANTS_FILE="src/main/java/progistar/scan/data/Constants.java"
VERSION=$(grep -m1 'VERSION' "$CONSTANTS_FILE" | sed -E 's/.*"(v[0-9]+\.[0-9]+\.[0-9]+)".*/\1/')

if [[ -z "$VERSION" ]]; then
  echo "Failed to read VERSION from $CONSTANTS_FILE" >&2
  exit 1
fi

echo "Building PepQueryMHC $VERSION ..."
mvn -q clean package -DskipTests

BUILT_JAR="target/PepQueryMHC.jar"
if [[ ! -f "$BUILT_JAR" ]]; then
  echo "Expected build output not found: $BUILT_JAR" >&2
  exit 1
fi

OUTPUT_JAR="target/PepQueryMHC.${VERSION}.jar"
mv "$BUILT_JAR" "$OUTPUT_JAR"

echo "Built: $OUTPUT_JAR"
