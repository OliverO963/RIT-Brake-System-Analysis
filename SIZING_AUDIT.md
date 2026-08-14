# Per-Component Sizing Audit — BrakeCoeffOptimizer.m & padSpecificPower.m

Audit of every physical "sizing" parameter in both scripts, checking that each is
**declared** and **used** on a strictly per-single-component basis:

- pad area = **one pad** (×2 → one corner, ×4 → one axle)
- rotor surface area = **one rotor**
- rotor mass = **one rotor**
- rotational inertia = **one corner**
- piston count = **one caliper**

---

## Summary

| Parameter | Declared basis | Used basis (before) | Verdict |
|---|---|---|---|
| `VehicleMass` | whole vehicle | whole vehicle | ✅ already correct |
| `RotorMass_front/rear` | one rotor | one rotor | ✅ already correct |
| `RotorArea_front/rear` | one rotor | one rotor | ✅ already correct |
| `I` (inertia) | one corner | one corner (×4 → vehicle) | ✅ correct, but aggregation fixed (D) |
| `A_pad_front/rear` | one pad | **one corner (2 pads)** | ❌ **fixed (A)** |
| `*_piston_count` | one caliper | **one caliper used as one side** | ❌ **fixed (B)** |
| rotor radius as lever arm | — | rotor **outer** radius | ❌ **fixed (C)** |

Two hard 2× errors were found, both in `padSpecificPower.m`. Four adjacent
correctness issues were also found and fixed.

---

## A. Pad energy divided by one pad's area, but computed for a whole corner — **2× error**

**Where:** `padSpecificPower.m`, `simulate_pad_power`

The energy chain is `whole vehicle → ×Tbias → one axle → ×corner_split → ONE CORNER`.
A corner has **one rotor** but **two pads**. The rotor branch divided corner energy
by one `RotorMass` (correct); the pad branch divided the same corner energy by
`A_pad_cm2`, which is **one pad's** area.

```matlab
% BEFORE - CorrectedEnergyPad is the energy into TWO pads
CorrectedEnergyPad = friction_energy(i) * corner_split * Tbias(i) * ... ;
pad_energy(i) = CorrectedEnergyPad;
q_inst(i)     = CorrectedEnergyPad / tbrake / A_pad_cm2;   % <- one pad's area

% AFTER
CorrectedEnergyPad = friction_energy(i) * PadFrac * BrakeFrac * CalibrationFactor;
pad_energy(i) = CorrectedEnergyPad / pads_per_corner;      % now ONE pad
q_inst(i)     = pad_energy(i) / tbrake / A_pad_cm2;
```

The code comments claimed the correct basis while the math did not — `% J, per
single pad, per step` on a variable holding two pads' worth of energy, and the file
header's `[applied per-pad]`.

**Impact:** all specific-power outputs (`q_inst`, `q_lag`, `avg_q`, fade-onset flux,
every Output 1–7 plot) were **2× overstated**.

### A2. Coupled change in the regen-failure margin (Output 8)

`P_front_peak` was **per corner** (`× TbiasF × 0.5`), which was *accidentally
self-consistent* with the corner-based `qOnset` — so the **margin factor was
correct** even though `q_inst` was 2× high (pad area cancels out of the ratio).

Halving `q_inst` alone would have silently halved the reported margin. Both were
moved to per-pad together:

```matlab
P_front_peak = P_total_peak * TbiasF * 0.5 / pads_per_corner;   % per SINGLE front pad
```

---

## B. Piston count is a per-caliper total but was used as single-side area — **2× error**

**Where:** `compute_derived_quantities` in **both** scripts

```matlab
% BEFORE - uses ALL 6 front / 4 rear pistons as the clamping area
front_piston_area = front_piston_count * pi * (front_piston_dia/2)^2;
fl_Tbrake = -2 * mu_front .* fl_clamp_force .* front_rotor_radius;
```

The `2` in the torque identity `T = 2·μ·F_clamp·r` is the **two friction faces**.
For that to be right, `F_clamp` must be the **single-side** force — a caliper is a
closed force loop, so the far side reacts the load rather than adding to it. With
opposed calipers (confirmed: 3+3 front, 2+2 rear) only half the pistons generate
clamp force.

Now handled once, in `loadVehicleParams.m`, driven by an explicit
`front/rear_caliper_opposed` flag so a floating caliper can be configured without
editing physics code:

```matlab
frontSides = 1 + (VP.front_caliper_opposed ~= 0);
VP.front_piston_area_per_side = (VP.front_piston_count / frontSides) * pi * (VP.front_piston_dia/2)^2;
```

**Impact:**
- `BrakeCoeffOptimizer.m`: **none numerically.** These feed only `Tbias_brake`, a
  ratio in which the factor cancels. *Verified:* mean `Tbias_brake` = **0.856 before
  and after**.
- `padSpecificPower.m`: halves `T_predicted`, so **`mu_ratio` roughly doubles.**
  This is why the bulk of `mu_ratio` previously sat well below 1.0. The
  `min_pressure_muratio_psi = 150` gate and the Tukey outlier fence were both
  calibrated against the old biased distribution and **should be re-checked.**

---

## C. Rotor outer radius used as the friction lever arm (~27% overstatement)

`front_rotor_radius = front_rotor_dia / 2` is the rotor **OD** radius; the correct
moment arm is the **mean pad radius**. Now a first-class documented spreadsheet
input (`front_pad_mean_radius` / `rear_pad_mean_radius`) rather than being derived
from rotor diameter. Rotor diameter remains in the sheet as reference geometry.

Same cancels-in-`Tbias` / doesn't-cancel-in-`mu_ratio` split as B.

> ⚠️ **Seeded at the current value** (rotor OD/2 = 0.09271 m front, 0.09398 m rear)
> so nothing changed silently. **Replace these with the documented mean pad radii** —
> this is the single highest-impact number still outstanding.

---

## D. Rotational KE aggregated per axle instead of per corner

```matlab
% BEFORE - mean of ONE axle's two omegas, squared, scaled to all four wheels
omegaP  = (omega_wheel_L(i-1) + omega_wheel_R(i-1)) / 2;
Energy2 = 4 * (0.5 * I * (omegaP^2 - omegaN^2));

% AFTER - each corner contributes using its OWN omega
Energy2 = sum(0.5 * I * (omega_all(i-1,:).^2 - omega_all(i,:).^2));
```

The old form made the front and rear calls compute *different* whole-vehicle
energies for the same timestep, and relied on `mean(ω)² == mean(ω²)`.

---

## E. Per-corner regen was pooled vehicle-wide before the corner split

Regen is measured at all four corners, but was summed and subtracted from vehicle
KE **before** the `×Tbias ×corner_split` step — so front regen reduced **rear** pad
heat in proportion to brake bias.

```matlab
% AFTER - split to this corner first, then subtract this corner's own regen
Energy_corner      = Energy * Tbias(i) * corner_split * (1 - AeroFrac);
regen_energy(i)    = abs(axle_regen_power(i)) * tbrake * corner_split;
friction_energy(i) = max(Energy_corner - regen_energy(i), 0);
```

`compute_derived_quantities` now exposes `front_axle_regen_power` /
`rear_axle_regen_power` instead of a single pooled `total_regen_power`.

---

## F. `padSpecificPower` applied the fitted coefficients outside their fit conditions

`BrakeCoeffOptimizer` derives a **per-dataset** ambient from each file's first rotor
readings and seeds the integration from the **measured** temperature.
`padSpecificPower` used a hardcoded 25 °C for both. Now mirrored, so the
coefficients are applied under the conditions they were fit in.

---

## G. Dead parameter removed

`WheelR` was threaded through ~9 call sites in both scripts and appeared in **zero**
expressions (flagged `%#ok<INUSD>`). Removed from all signatures; retained in the
spreadsheet as a real vehicle parameter.

---

## Centralized parameter spreadsheet

`VehicleParameters.xlsx` (sheets `Parameters` + `MuTable`) is now the single source
for both scripts, read via the new shared `loadVehicleParams.m`. Every row carries
an explicit **Basis** column — undocumented basis is what caused every bug above.

38 parameters centralized, including values that were previously duplicated
literals in both files (pressure/temperature ADC calibration, `BrakeFrac`,
`CalibrationFactor`, `min_pressure`, aero polynomial, mu table). `BrakeCoeffOptimizer`
had four of these as *anonymous literals* (`0.82`, `0.8`, `0.5`, `5`) — editing the
named copy in `padSpecificPower` used to silently desynchronize the two.

The loader validates that every expected parameter is present and numeric, and
errors clearly on a renamed/deleted row rather than propagating an undefined field
into the physics.

**Values were transcribed exactly as they existed — nothing was rescaled.**

---

## ⚠️ Consequence: the logged padFrac coefficients are now stale

Fixes **D, E, and F** change the energy/temperature model itself, so the six
`h_w`/PadFrac coefficient sets currently hardcoded in `padSpecificPower.m` no longer
correspond to the physics that produced them. `BrakeCoeffOptimizer` must be re-run
and the model list updated. (Fixes A/B/C do **not** affect the fit — they cancel out
of `Tbias`.)

---

## Values to verify against the spec sheet

| Value | Concern |
|---|---|
| `front/rear_pad_mean_radius` | **Currently seeded at rotor OD/2.** Replace with documented mean pad radius (~27% overstatement while unreplaced). |
| `RotorArea_rear` = 0.0226 m² | 1.68× smaller than front despite near-equal rotor diameters — front/rear appear to use different geometric conventions. |
| `rear_rotor_dia` = 0.18796 m | Rear rotor is currently **larger** than front (0.18542 m). Confirm not transposed. |
| `VehicleMass` = 259 kg | Confirm whether driver mass is included (used for both KE and the regen-failure case). |
| `front/rear_caliper_opposed` = 1 | Assumes opposed 3+3 / 2+2. Set to 0 if floating. |
