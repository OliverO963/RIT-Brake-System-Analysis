### 8/18/2026 (later)

Add regen fraction vs charge current limit



\- Add a figure sweeping the charge current limit (10-70 A) and plotting the

&#x20; modeled regen share of dissipated energy for each drive type, with the

&#x20; measured share at the recorded 35 A limit overlaid for comparison. The

&#x20; denominator is held fixed across the sweep: raising the limit moves energy

&#x20; from the friction brakes to the motors, it does not change how much energy

&#x20; the lap required the car to shed

\- Change regen energy vs deceleration rate to regen SHARE vs deceleration

&#x20; rate. A harder stop sheds more of everything, so absolute regen energy

&#x20; against decel rate mostly restated that the stop was big; the share shows

&#x20; how well regen keeps up as braking gets harder. Event duration keeps

&#x20; absolute energy, where the magnitude is the point

\- Factor the regen torque model into predict\_regen\_power so the swept curve

&#x20; and the single-limit console figure cannot drift apart (verified equal to

&#x20; 5 decimal places at the recorded limit)



### 8/18/2026

Add regenerative braking energy contribution analysis



\- Add regenAnalysis.m: measures how much of the energy dissipated in a braking

&#x20; event is recovered by regen, across as many datasets as the user selects

\- Combine BrakeDataAnalysis.m's regen physics with the multi-dataset front end

&#x20; of padSpecificPower.m and BrakeCoeffOptimizer.m (multi-select picker,

&#x20; 21/22/27-column auto-detection, drive-type detection, event segmentation)

\- Report regen contribution two ways, energy-weighted and mean-per-event,

&#x20; broken down by autocross / endurance / other-misc

\- Define a braking event as ANY deceleration event, so pure-regen stops are

&#x20; counted; excluding them would bias the regen share downward

\- Exclude datasets recorded without regen so they cannot dilute the averages.

&#x20; Both 4-17 endurance sessions recover ~0.2% of dissipated energy against

&#x20; 10-48% for every later session, and are dropped

\- Add the regen parameters (current limit, torque/power limits, max regen curve

&#x20; coefficients, efficiency, speed and energy thresholds) to

&#x20; VehicleParameters.xlsx, so they are shared rather than duplicated per script



Notes on data format:

\- Pack current exists only in the raw logger files. Preprocessed exports from

&#x20; F34BrakeDataPreprocessor.m drop it, and the column in that position holds

&#x20; Steering\_pct instead, so the channel is located by name and both the regen

&#x20; index and the no-regen dataset gate fall back to regen torque/energy

\- Mechanical braking energy is half what BrakeDataAnalysis.m reports for the

&#x20; same file. That script put all 6 front pistons into the clamp force; the

&#x20; caliper is opposed with 3 per side (confirmed), so clamp force uses one

&#x20; side's 3 and the factor of 2 in T = 2\*mu\*F\*r is the two friction faces.

&#x20; Regen and aero energies match that script exactly. An energy budget from

&#x20; vehicle dynamics corroborates: the 6-piston figure claims 126% of the

&#x20; energy available to dissipate on AutoX 5-27, the 3-per-side figure 84%



Next: validate the regen model coefficients against a wider set of sessions

### 7/30/2026

Refine data preprocessing and run coefficient optimizer on brake data



\- Process 12 recorded sessions, yielding 7 usable files after preprocessing

\- Improve session boundary detection: apply adaptive start/end offsets based on

&#x20; inter-session gaps to better capture cooling behavior and ambient conditions

\- Add/refine peak filtering to remove outliers in brake pressure data

\- Run coefficient optimizer on preprocessed data to explore padFrac relationships

\- Confirm linear model (temperature + pressure + interaction term) provides

&#x20; best fit for padFrac coefficients

\- Generate braking event specific power script (untested; pending validated fits)



Known issues to address:

\- Some non-session records still passing preprocessing filters

\- Coefficient optimizer shows inconsistent padFrac outputs (e.g., zeros above

&#x20; 400psi); likely resolves with curated input dataset

\- Peak filtering behavior needs careful validation



Next: curate input data, validate preprocessing filters, troubleshoot optimizer

