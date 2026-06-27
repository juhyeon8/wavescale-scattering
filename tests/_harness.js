"use strict";
let pass = 0, fail = 0;
function test(name, fn) {
  try { fn(); pass++; console.log("  ok   -", name); }
  catch (e) { fail++; console.error("  FAIL -", name, "\n   ", e.message); }
}
function approx(actual, expected, tol, msg) {
  if (Math.abs(actual - expected) > tol)
    throw new Error((msg || "approx") + `: got ${actual}, expected ${expected} ±${tol}`);
}
function done() {
  console.log(`\n${pass} passed, ${fail} failed`);
  if (fail) process.exit(1);
}
module.exports = { test, approx, done };
