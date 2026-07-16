#!/usr/bin/env bash
set -eu

if [ -t 1 ] && [ "${NO_COLOR:-}" = "" ]; then
  COLOR_GREEN=$'\033[32m'
  COLOR_RED=$'\033[31m'
  COLOR_RESET=$'\033[0m'
else
  COLOR_GREEN=""
  COLOR_RED=""
  COLOR_RESET=""
fi

print_result_line() {
  local name="$1"
  local status="$2"
  local color="$3"
  printf "%-68s %s%s%s\n" "$name" "$color" "$status" "$COLOR_RESET"
}

run_test() {
  local name="$1"
  shift

  local log
  log="$(mktemp "${TMPDIR:-/tmp}/objcryst-test.XXXXXX.log")"

  if "$@" >"$log" 2>&1; then
    print_result_line "$name" "PASS" "$COLOR_GREEN"
    rm -f "$log"
  else
    print_result_line "$name" "FAIL" "$COLOR_RED"
    echo "---- Captured output for: $name ----"
    cat "$log"
    echo "---- End captured output ----"
    rm -f "$log"
    return 1
  fi
}
