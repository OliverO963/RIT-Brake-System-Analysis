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

