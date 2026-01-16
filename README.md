# Anti-Aliasing Signal Recovery Toolkit

<p align="center">
  <img src="docs/figures/aliasing_demo.png" alt="Aliasing Demo" width="600">
</p>

Signal recovery methods under sub-Nyquist sampling conditions.

---

## 🎯 Overview

When the sampling rate $F_s$ is below the Nyquist rate ($F_s < 2 f_{max}$), high-frequency components **alias** to lower frequencies, making direct recovery impossible.

This toolkit demonstrates **6 different methods** to overcome aliasing:

| Method | Key Idea | When to Use |
|--------|----------|-------------|
| [Bandpass Sampling](methods/01_bandpass_sampling/) | Narrowband signals can be sampled below Nyquist | Known signal bandwidth |
| [CRT Multi-rate](methods/02_crt_multirate/) | Multiple sample rates + Chinese Remainder Theorem | Multiple ADCs available |
| [Staggered PRF](methods/03_staggered_prf/) | Alternating pulse intervals | Pulsed radar/sonar |
| [Non-uniform Sampling](methods/04_nonuniform_sampling/) | Irregular sampling breaks aliasing periodicity | Flexible sampling times |
| [Compressed Sensing](methods/05_compressed_sensing/) | Exploit signal sparsity | Sparse spectra |
| [Parametric Methods](methods/06_parametric_methods/) | MUSIC/ESPRIT assume signal model | Known number of sinusoids |

---

## 📂 Repository Structure

```
anti-aliasing-recovery/
├── README.md                          # This file
├── docs/
│   ├── theory_overview.md            # Mathematical foundations
│   └── figures/                      # Diagrams and results
│
├── methods/
│   ├── 01_bandpass_sampling/
│   │   ├── README.md
│   │   ├── demo_1d_rf_signal.m       # 1D: RF signal undersampling
│   │   ├── demo_2d_ultrasound.m      # 2D: Ultrasound imaging
│   │   └── app_sdr_receiver.m        # Application: SDR
│   │
│   ├── 02_crt_multirate/
│   │   ├── README.md
│   │   ├── demo_1d_frequency_est.m   # 1D: Dual-ADC frequency estimation
│   │   ├── demo_2d_phase_unwrap.m    # 2D: InSAR phase unwrapping
│   │   └── app_tof_ranging.m         # Application: ToF camera
│   │
│   ├── 03_staggered_prf/
│   │   ├── README.md
│   │   ├── demo_1d_doppler.m         # 1D: Doppler velocity
│   │   ├── demo_2d_weather_radar.m   # 2D: Weather radar
│   │   └── app_velocity_extension.m  # Application: Automotive radar
│   │
│   ├── 04_nonuniform_sampling/
│   │   ├── README.md
│   │   ├── demo_1d_irregular.m       # 1D: Irregular sampling
│   │   ├── demo_2d_mri_kspace.m      # 2D: MRI k-space
│   │   └── app_astronomical.m        # Application: Astronomy
│   │
│   ├── 05_compressed_sensing/
│   │   ├── README.md
│   │   ├── demo_1d_sparse_recovery.m # 1D: Sparse spectrum recovery
│   │   ├── demo_2d_image_cs.m        # 2D: Image CS
│   │   └── app_radar_imaging.m       # Application: SAR/ISAR
│   │
│   └── 06_parametric_methods/
│       ├── README.md
│       ├── demo_1d_music_esprit.m    # 1D: Super-resolution
│       ├── demo_2d_doa.m             # 2D: 2D DOA
│       └── app_array_processing.m    # Application: Array signal
│
└── utils/
    ├── alias_freq.m                  # Compute aliased frequency
    ├── crt_solver.m                  # CRT recovery algorithm
    └── plot_spectrum.m               # Spectrum visualization
```

---

## 🚀 Quick Start

```matlab
% Add paths
addpath(genpath('utils'));
addpath(genpath('methods'));

% Run CRT frequency estimation demo
cd methods/02_crt_multirate
demo_1d_frequency_est
```

---

## 📖 Why Does Aliasing Happen?

<p align="center">
  <img src="docs/figures/aliasing_wheel.gif" alt="Wagon Wheel Effect" width="400">
</p>

**The Wagon Wheel Effect**: A wheel spinning at 10 rev/s filmed at 6 fps appears to spin at 4 rev/s backward!

Mathematically:
$$f_{alias} = ((f_{true} + F_s/2) \mod F_s) - F_s/2$$

Multiple true frequencies map to the **same** aliased frequency:
```
f_true:  50, 150, 250, 350, ... Hz
Fs:      100 Hz
f_alias: 50,  50,  50,  50, ... Hz  (all look the same!)
```

---

## 🔬 Method Comparison

| Method | Complexity | Hardware Req. | Best For |
|--------|------------|---------------|----------|
| Bandpass | Low | Single ADC | Narrowband RF |
| CRT | Medium | 2+ ADCs | General |
| Staggered PRF | Medium | Variable timing | Pulsed systems |
| Non-uniform | High | Flexible timing | Astronomy |
| Compressed Sensing | High | Random sampling | Sparse signals |
| MUSIC/ESPRIT | Medium | Single ADC | Few sinusoids |

---

## 🤔 How to Choose?

```
                    ┌─ Is signal narrowband? ─── Yes ──► Bandpass Sampling
                    │
Start ──► What do   ├─ Multiple ADCs available? ─ Yes ──► CRT Multi-rate
          you have? │
                    ├─ Can vary sample timing? ─── Yes ──► Staggered PRF
                    │                                      or Non-uniform
                    │
                    ├─ Is spectrum sparse? ────── Yes ──► Compressed Sensing
                    │
                    └─ Know # of sinusoids? ───── Yes ──► MUSIC/ESPRIT
```

---

## 📚 References

1. **Bandpass Sampling**: Vaughan et al., "The Theory of Bandpass Sampling," IEEE TSP, 1991.
2. **CRT**: Xia & Wang, "Phase Unwrapping and a Robust Chinese Remainder Theorem," IEEE SPL, 2007.
3. **Compressed Sensing**: Candès et al., "Robust Uncertainty Principles," IEEE TIT, 2006.
4. **MUSIC**: Schmidt, "Multiple Emitter Location and Signal Parameter Estimation," IEEE TAP, 1986.

---

## 📝 License

MIT License - Feel free to use for education and research.


