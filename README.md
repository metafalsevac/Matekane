# Matekane
Numerical integration of the differential state-space equations of a simulated aircraft taking off
from the Matekane airstrip, using 4th-order Runge-Kutta (plus optional adaptive and geometric
integrators). Two models are available:

* **2D** – the original vertical-plane point-mass model (flight-path angle formulation).
* **3D** – a full 6-degree-of-freedom rigid-body model with unit-quaternion attitude.

The scenario: the aircraft leaves the end of the strip at 7649.63 ft with V = 83 ft/s and
γ = 0, over a valley whose floor drops to ~6048 ft and rises again to level ground (7480.58 ft)
at x = 6233.58 ft. For each lift coefficient in a sweep the program reports **CRASH** (terrain
contact, with interpolated touchdown point) or **SURVIVE** (clears x = 10000 ft). Lift
coefficients above CLmax = 2.195 are skipped (stalled wing).

## Build and run

Requires only a C++17 compiler (`g++`/`clang++`) and `make`.

```sh
make                      # builds ./matekane
make test                 # builds and runs the test suite

./matekane                                    # 2D, RK4, prompts for CL sweep (original behaviour)
./matekane --cl-min 1.0 --cl-max 2.1 --cl-step 0.1
./matekane --mode 3d --cl 1.5                 # 6-DOF, quaternion attitude
./matekane --mode 3d --cl 1.8 --aileron 2 --integrator lie-rk4 --stats
./matekane --integrator rk45 --rtol 1e-10 --atol 1e-10 --cl 1.5
./matekane --help
```

| option | meaning |
|---|---|
| `--mode 2d\|3d` | model (default `2d`) |
| `--integrator rk4\|rk45\|lie-rk4\|lie-rk45` | time integrator (default `rk4`; `lie-*` require `--mode 3d`) |
| `--dt S` | RK4 step and output interval for all integrators (default 0.1 s) |
| `--rtol X`, `--atol X` | RK45 tolerances (default 1e-8) |
| `--cl X` / `--cl-min --cl-max --cl-step` | lift-coefficient value / sweep |
| `--tmax S` | simulated-time limit per run (default 600 s) |
| `--aileron --elevator --rudder DEG` | constant control deflections (3D) |
| `--no-renorm` | disable quaternion renormalisation after `rk4`/`rk45` steps (3D) |
| `--out FILE`, `--stats` | output file, print step / function-evaluation counts |

Output: `output.txt` (2D; columns `x h V ground gamma Veas CL`, same as the original program plus a
`#` header) or `output3d.csv` (3D; time, NED position, altitude, terrain, body velocity, quaternion,
body rates, Euler angles, α, β, V, CL, |q|−1). Optional plots (needs `numpy` + `matplotlib`):

```sh
python3 scripts/plot.py output.txt      # altitude vs distance per CL, with terrain
python3 scripts/plot.py output3d.csv    # 3D trajectory over terrain + attitude/alpha/beta/|q|-1 vs time
```

## Units, frames and conventions

US customary units throughout (ft, slug, lbf, s, °R); angles in radians internally (the CLI takes
control deflections in degrees, the CSV reports Euler angles/α/β in degrees).

* **Inertial frame: NED** (flat, non-rotating Earth): x North = down-valley distance, y East,
  z Down. Altitude h = −z. The terrain varies only with x.
* **Body frame: FRD** – x out the nose, y out the right wing, z down; origin at the CG.
* **Quaternion: scalar-first, Hamilton** `q = [w, x, y, z]`, unit norm, representing the body
  attitude relative to NED. It rotates body vectors into NED: `v_ned = q ⊗ [0, v_b] ⊗ q* = R(q) v_b`.
  Kinematics: `q̇ = ½ q ⊗ [0, ω_b]` with ω_b = (p, q, r) the body angular rates.
  Euler angles (3-2-1 yaw ψ, pitch θ, roll φ) are derived from q for output only.
* Wind axes: α = atan2(w, u), β = asin(v / V).
* Control signs: +aileron rolls right, +elevator pitches nose down, +rudder yaws nose left.

## 2D model (original)

State `y = [V, γ, h, x]`, CL constant (the swept parameter), thrust along the velocity vector:

```
V̇ = (T − D)/m − g sin γ          T = σ (T0SL − aSL V²),  σ = ρ/ρ0
γ̇ = L/(m V) − g cos γ / V        L = ½ ρ V² S CL,  D = ½ ρ V² S (0.025 + 0.056 CL²)
ḣ = V sin γ
ẋ = V cos γ
```

ISA troposphere: T = T0 + a h, p = p0 (1 + a h/T0)^5.2561, ρ = p/(R T).
W = 2712.5 lbf, m = W/g, S = 195.75 ft², T0SL = 714.9375 lbf, aSL = 0.0334 lbf·s²/ft².

## 3D model (6-DOF, quaternion attitude)

State (13): `[p_ned(3) | v_b = (u, v, w) | q = (w, x, y, z) | ω_b = (p, q, r)]`

```
ṗ_ned = R(q) v_b
v̇_b   = F_b / m − ω_b × v_b
q̇     = ½ q ⊗ [0, ω_b]                 (renormalised after each rk4/rk45 step)
ω̇_b   = I⁻¹ (M_b − ω_b × I ω_b)        I = [[Ixx,0,−Ixz],[0,Iyy,0],[−Ixz,0,Izz]]
F_b   = R(q)ᵀ [0, 0, m g] + [T, 0, 0] + C_bw [−D, Y, −L]
```

Thrust uses the same T(ρ, V) law as 2D and acts along body x. Aerodynamics are a linear
stability-derivative generalisation of the 2D model (q̄ = ½ρV², p̂ = pb/2V, q̂ = qc̄/2V, r̂ = rb/2V):

```
CL = CL_cmd + CLα α + CLq q̂ + CLδe δe        CD = 0.025 + 0.056 CL²   (2D polar)
CY = CYβ β + CYδr δr
Cl = Clβ β + Clp p̂ + Clr r̂ + Clδa δa + Clδr δr
Cm = Cm0 + Cmα α + Cmq q̂ + Cmδe δe
Cn = Cnβ β + Cnp p̂ + Cnr r̂ + Cnδa δa + Cnδr δr
```

`CL_cmd` is the swept CL, interpreted as the lift coefficient of the (high-lift) configuration at
zero body angle of attack. With Cm0 = 0 the aircraft is trimmed at α = 0, where CL = CL_cmd and
thrust is aligned with the velocity – exactly the 2D assumptions. Default coefficients
(`src/sim3d.hpp`, `Params3D`) are representative Cessna-172-class values; they are **not** data for
the original assignment aircraft, which only specified the 2D quantities:

| | | | |
|---|---|---|---|
| b = 36 ft | c̄ = S/b = 5.4375 ft | Ixx = 948, Iyy = 1346, Izz = 1967 slug·ft² | Ixz = 0 |
| CLα = 4.6 | CLq = 3.9 | CLδe = 0.43 | |
| CYβ = −0.31 | CYδr = 0.187 | | |
| Clβ = −0.089 | Clp = −0.47 | Clr = 0.096 | Clδa = 0.178, Clδr = 0.0147 |
| Cm0 = 0 | Cmα = −0.89 | Cmq = −12.4 | Cmδe = −1.28 |
| Cnβ = 0.065 | Cnp = −0.03 | Cnr = −0.099 | Cnδa = −0.053, Cnδr = −0.0657 |

Because the 3D aircraft must rotate to follow its flight path (short-period dynamics, pitch
damping), its results differ somewhat from 2D for the same CL (e.g. CL = 1.0 touches down at
x ≈ 6073 ft instead of 6109 ft). With a very stiff pitch axis (|Cmα| → ∞, CLq = 0) the 3D model
converges to the 2D one; this is checked by the tests. Not modelled: stall/post-stall (linear lift
curve; CL > CLmax only skipped at input), ground effect, ground roll, wind, propeller torque,
Earth rotation/curvature, terrain variation in the East direction.

## Integrators

| name | type | attitude handling | notes |
|---|---|---|---|
| `rk4` (default) | classic RK4, fixed step `--dt` | q̇ integrated as a 4-vector, then renormalised | the original scheme; simple, predictable cost |
| `rk45` | Dormand–Prince 5(4), adaptive, FSAL, 4th-order dense output | as rk4 | error-controlled (`--rtol/--atol`); output at the `--dt` grid via dense output |
| `lie-rk4` | Runge–Kutta–Munthe-Kaas with the RK4 tableau | `q ← q ⊗ exp(Θ)`, Θ̇ = dexp⁻¹(Θ) ω | 3D only; |q| = 1 to rounding with no renormalisation |
| `lie-rk45` | RKMK with the Dormand–Prince tableau, adaptive | as lie-rk4 | 3D only |

RKMK: in each step the attitude is written `q = q_n ⊗ exp(Θ)` with Θ ∈ ℝ³ (rotation vector),
Θ(t_n) = 0. The RK tableau is applied to `[p, v_b, Θ, ω_b]` with
`Θ̇ = Jr⁻¹(Θ) ω = ω + ½ Θ×ω + (1 − (|Θ|/2) cot(|Θ|/2))/|Θ|² · Θ×(Θ×ω)` (exact inverse right
Jacobian of SO(3), so the method keeps the full order of the tableau), and `q_{n+1} = q_n ⊗ exp(Θ_{n+1})`.

Trade-offs:

* **rk4** – cheapest per step (4 f-evals), fixed cost, no error estimate. Renormalisation is a
  projection: harmless for accuracy at these step sizes but not "structure preserving".
* **rk45** – 6 f-evals per accepted step (FSAL) but takes large steps where the solution is smooth:
  for the 2D model at the same accuracy (~2e-8) it needs 135 steps / 811 evaluations vs 600 / 2400
  for RK4. Step size adapts to transients (e.g. short-period oscillations in 3D). Event detection
  (terrain contact) is done on the `--dt` output grid using dense output, and touchdown is
  interpolated between grid points.
* **lie-rk4 / lie-rk45** – the quaternion never leaves the unit sphere (no projection, no drift
  even over 10⁵ steps of a tumbling body), at the cost of a few trig evaluations per stage.
  Accuracy order is the same as the corresponding classical method. Preferred for long
  simulations or aggressive rotation; for this takeoff problem all four agree to ~5 significant digits.

## Tests

`make test` runs `tests/test_main.cpp` (no external framework):

* RK4 observed order ≈ 4 (scalar ODE with exact solution, and the 2D aircraft model);
* the refactored 2D model reproduces the original program's output;
* Dormand–Prince tableau: 5th-order global convergence, dense output local error O(h⁵);
* RK45 meets its tolerance vs a fine-step reference and is cheaper than RK4 at equal accuracy;
* |q| stays 1 during rolling/turning flight for all integrators;
* pure rotation with constant ω gives `q0 ⊗ exp(ωT)` (and yaw-only rotation gives ψ = rT);
* geometric integrators keep |q| = 1 to rounding **without** renormalisation, are 4th order,
  and conserve angular momentum/energy of a torque-free tumbling body;
* free fall / ballistic flight (with and without body rotation) matches the analytic parabola;
* 3D constrained to the pitch plane reproduces 2D (lateral states stay exactly zero; differences
  shrink ∝ 1/|Cmα|).
