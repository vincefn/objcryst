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
  local timeout_bin=""
  local test_timeout="${TEST_TIMEOUT:-}"
  local live_output="${TEST_LIVE_OUTPUT:-0}"
  local status=0
  log="$(mktemp "${TMPDIR:-/tmp}/objcryst-test.XXXXXX")"

  if command -v timeout >/dev/null 2>&1; then
    timeout_bin="timeout"
  elif command -v gtimeout >/dev/null 2>&1; then
    timeout_bin="gtimeout"
  fi

  if [ "$live_output" = "1" ]; then
    set +e
    if [ -n "$test_timeout" ] && [ -n "$timeout_bin" ]; then
      "$timeout_bin" --signal=TERM --kill-after=20s "$test_timeout" "$@" 2>&1 | tee "$log"
      status=${PIPESTATUS[0]}
    else
      "$@" 2>&1 | tee "$log"
      status=${PIPESTATUS[0]}
    fi
    set -e
  else
    set +e
    if [ -n "$test_timeout" ] && [ -n "$timeout_bin" ]; then
      "$timeout_bin" --signal=TERM --kill-after=20s "$test_timeout" "$@" >"$log" 2>&1
      status=$?
    else
      "$@" >"$log" 2>&1
      status=$?
    fi
    set -e
  fi

  if [ "$status" -eq 0 ]; then
    print_result_line "$name" "PASS" "$COLOR_GREEN"
    rm -f "$log"
    return 0
  fi

  print_result_line "$name" "FAIL" "$COLOR_RED"
  echo "---- Captured output for: $name ----"
  cat "$log"
  echo "---- End captured output ----"
  rm -f "$log"
  return 1
}
