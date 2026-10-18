# numeric

Root finding, integration and differential equations for
[Tessel](https://tessel-lang.syncster.dev) programs: the methods
engineering calculations keep needing, each saying how its answer came out
(whether it converged, in how many steps, how accurate it is) rather than
quietly returning `NaN`.

```tessel
import numeric

fn main() {
    // The depth (m) at which a 2 m wide channel carries 3 m³/s
    // (Manning's n 0.013, slope 0.001).
    let r = numeric.brent(low: 0.01, high: 10.0) { y in
        let area = 2.0 * y
        let radius = area / (2.0 + 2.0 * y)
        area * pow(radius, 2.0 / 3.0) * sqrt(0.001) / 0.013 - 3.0
    }
    if r.converged {
        print("normal depth: {r.value} m, in {r.iterations} steps")
    }
}
```

```text
normal depth: 0.9841196687128626 m, in 8 steps
```

## Using it

Add it to a project (Tessel 0.8 or later):

```sh
tessel add numeric
```

then `import numeric` and call its functions as `numeric.brent(…)` (or,
with Tessel 0.9 or later, `import numeric as num` and `num.brent(…)`). The
function each method works on comes last, so it's written as a block after
the call, with its parameters before `in`.

## Roots: where f(x) = 0

Write the problem as f(x) = 0; for "f(x) equals a target", use
`f(x) - target`.

| Function | Needs | Use it when |
|---|---|---|
| `brent(low:high:) { x in … }` | f(low) and f(high) of opposite signs | Almost always: it can't miss a root in the range, and is fast. |
| `bisection(low:high:) { x in … }` | The same | You want the simplest method, and speed doesn't matter. |
| `newton(start:derivative:) { x in … }` | A first guess, and f's derivative | The derivative is known and the guess is good. |
| `secant(start:next:) { x in … }` | Two first guesses | No derivative, and the guesses are good. |
| `bracket(start:step:) { x in … }` | A starting point | Finding a range for `brent`, widening outwards. |
| `scan(low:high:steps:) { x in … }` | A range to look in | Finding a range for `brent` where x must stay within limits (a depth above zero), or the root nearest `low`. |

Each root finder takes `tolerance:` (in x, default `1e-12`) and
`maxIterations:`, and returns a `Root`:

| Field | |
|---|---|
| `value` | The root: `NaN` when there's no answer (no sign change). |
| `residual` | f(value), to check how close to zero it really is. |
| `iterations` | How many steps it took. |
| `converged` | Whether it met the tolerance. When `false`, `value` isn't an answer. |

`bracket` and `scan` return a `Bracket`: `low`, `high` and `found`.

```tessel
let range = numeric.scan(low: 0.0, high: 10.0, steps: 100) { x in sin(x) - 0.5 }
if range.found {
    let r = numeric.brent(low: range.low, high: range.high) { x in sin(x) - 0.5 }
    print(r.value)   // 0.5235987755982989 (pi / 6)
}
```

## Integration: the area under a curve

| Function | |
|---|---|
| `integrate(from:to:) { x in … }` | To a tolerance (`tolerance:`, default `1e-10`), taking small steps only where f bends sharply. Returns an `Integral`. |
| `simpson(from:to:steps:) { x in … }` | Simpson's rule with equal steps; exact for cubics. |
| `trapezoid(from:to:steps:) { x in … }` | The trapezoid rule with equal steps; exact for straight lines. |
| `trapezoidPoints(x:y:)` | Measured points rather than a function: the x's needn't be evenly spaced. |

An `Integral` has the `value`, an `error` estimate, the number of
`evaluations` of f, and whether it `converged` (false when f has a jump or
spike that `maxDepth:` halvings can't resolve: the value is then still the
best estimate).

```tessel
let r = numeric.integrate(from: 0.0, to: pi) { x in sin(x) }
print(r.value)   // 2.0, to about 1e-10 (r.error says how close)

// Flows (m³/s) read every hour: the volume that passed, in m³.
let volume = numeric.trapezoidPoints(x: [0.0, 3600.0, 7200.0, 10800.0], y: [10.0, 40.0, 25.0, 12.0])
print(volume)   // 273600.0
```

## Differential equations: dy/dx = f(x, y)

| Function | |
|---|---|
| `solve(x:y:to:) { x, y in … }` | To a tolerance (`tolerance:`, default `1e-8`), choosing its own steps (Dormand-Prince). Returns a `Solution`. |
| `rk4(x:y:to:steps:) { x, y in … }` | Fourth-order Runge-Kutta with equal steps: a `[Sample]` at regular x's. |
| `rk4System(x:y:to:steps:) { x, y in … }` | The same for several quantities that change together, `y` being a `[Float]`: a `[State]`. |

All of them go backwards when `to` is less than `x`, as a water surface
profile computed upstream from a control does. A `Sample` is an `x` and a
`y`; a `State` an `x` and a `[Float]` `y`. A `Solution` has the `samples`,
the number of `steps`, and whether it `converged` (reached `to`): when it
didn't, the samples end where it stopped, often where the solution blows
up or f gives `NaN`.

```tessel
// A reservoir's level h (m) through a day: 50 m³/s in, Q = 30 h^1.5 out
// over the spillway, 2 km² of surface.
let day = numeric.solve(x: 0.0, y: 0.5, to: 86400.0) { t, h in
    (50.0 - 30.0 * pow(max(h, 0.0), 1.5)) / 2000000.0
}
print(day.samples[day.samples.count - 1].y)   // about 1.30

// A second-order equation as a system: a pendulum (small angles),
// angle'' = -angle, with y[0] the angle and y[1] its rate.
let swing = numeric.rk4System(x: 0.0, y: [1.0, 0.0], to: pi, steps: 100) { t, s in [s[1], -s[0]] }
print(swing[100].y[0])   // about -1.0
```

## Tests

```sh
tessel run tests
```

checks every method against answers known exactly, and every example in
the documentation, and exits with 1 if any check fails.

## What's where

| File | |
|---|---|
| `roots.tsl` | `brent`, `bisection`, `newton`, `secant`, `bracket`, `scan`, `Root` and `Bracket`. |
| `integrate.tsl` | `integrate`, `simpson`, `trapezoid`, `trapezoidPoints` and `Integral`. |
| `ode.tsl` | `solve`, `rk4`, `rk4System`, `Sample`, `State` and `Solution`. |
| `tests/main.tsl` | The tests. |
| `tessel.toml` | The library's name, version, description and licence. |

Every function's comment in the code says what it does, its parameters,
what it returns, and gives an example.

## Licence

MIT: see [LICENSE](LICENSE).
