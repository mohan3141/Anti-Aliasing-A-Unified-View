# Compressed Sensing for Anti-Aliasing Recovery

Compressed Sensing (CS) enables signal recovery from far fewer samples than Nyquist rate requires, by exploiting **sparsity** in some transform domain.

---

## 1. Theory Background

### 1.1 The Core Problem

Traditional sampling (Nyquist-Shannon):
$$F_s \geq 2 f_{max}$$

Compressed Sensing breakthrough: If signal is **K-sparse** in some basis, we only need:
$$M \geq C \cdot K \cdot \log(N/K)$$

measurements, where $M \ll N$.

### 1.2 Mathematical Framework

**Signal Model:**
$$\mathbf{x} = \mathbf{\Psi} \mathbf{s}$$

where:
- $\mathbf{x} \in \mathbb{R}^N$: Original signal
- $\mathbf{\Psi} \in \mathbb{R}^{N \times N}$: Sparsifying basis (e.g., Fourier, DCT, Wavelet)
- $\mathbf{s} \in \mathbb{R}^N$: Sparse coefficient vector ($K$ non-zeros)

**Measurement Process:**
$$\mathbf{y} = \mathbf{\Phi} \mathbf{x} = \mathbf{\Phi} \mathbf{\Psi} \mathbf{s} = \mathbf{A} \mathbf{s}$$

where:
- $\mathbf{y} \in \mathbb{R}^M$: Measurements ($M \ll N$)
- $\mathbf{\Phi} \in \mathbb{R}^{M \times N}$: Measurement matrix
- $\mathbf{A} = \mathbf{\Phi} \mathbf{\Psi}$: Sensing matrix

**Recovery Problem:**
$$\min_{\mathbf{s}} \|\mathbf{s}\|_1 \quad \text{s.t.} \quad \mathbf{A}\mathbf{s} = \mathbf{y}$$

or with noise:
$$\min_{\mathbf{s}} \|\mathbf{s}\|_1 \quad \text{s.t.} \quad \|\mathbf{A}\mathbf{s} - \mathbf{y}\|_2 \leq \epsilon$$

### 1.3 Key Conditions

#### Restricted Isometry Property (RIP)

Matrix $\mathbf{A}$ satisfies RIP of order $K$ with constant $\delta_K$ if:
$$(1 - \delta_K) \|\mathbf{s}\|_2^2 \leq \|\mathbf{A}\mathbf{s}\|_2^2 \leq (1 + \delta_K) \|\mathbf{s}\|_2^2$$

for all $K$-sparse vectors $\mathbf{s}$.

**Theorem:** If $\delta_{2K} < \sqrt{2} - 1 \approx 0.41$, then L1-minimization exactly recovers $\mathbf{s}$.

#### Incoherence

The measurement basis $\mathbf{\Phi}$ should be "incoherent" with sparsity basis $\mathbf{\Psi}$:
$$\mu(\mathbf{\Phi}, \mathbf{\Psi}) = \sqrt{N} \cdot \max_{i,j} |\langle \phi_i, \psi_j \rangle|$$

Lower $\mu$ → fewer measurements needed.

| Basis Pair | Coherence |
|------------|-----------|
| Time ↔ Fourier | $\mu = 1$ (maximally incoherent) |
| Random ↔ Any | $\mu \approx \sqrt{2 \log N}$ |
| DCT ↔ Wavelet | Low |

---

## 2. Recovery Algorithms

### 2.1 L1-Minimization (Basis Pursuit)

$$\min_{\mathbf{s}} \|\mathbf{s}\|_1 \quad \text{s.t.} \quad \mathbf{A}\mathbf{s} = \mathbf{y}$$

**Algorithms:**
- Interior Point Methods
- ADMM (Alternating Direction Method of Multipliers)
- Proximal Gradient Methods (ISTA, FISTA)

**Pros:** Provable recovery guarantees  
**Cons:** Computationally expensive for large problems

### 2.2 Greedy Algorithms

#### Orthogonal Matching Pursuit (OMP)

```
1. Initialize: residual r = y, support S = {}
2. For i = 1 to K:
   a. Find: j* = argmax_j |<r, a_j>|
   b. Update: S = S ∪ {j*}
   c. Project: s_S = A_S† y
   d. Update residual: r = y - A_S s_S
3. Return: s
```

**Pros:** Fast, simple implementation  
**Cons:** May fail for highly correlated dictionaries

#### Other Greedy Methods

| Algorithm | Complexity | Notes |
|-----------|------------|-------|
| OMP | O(KMN) | Most common |
| CoSaMP | O(MN log N) | Adds/removes elements |
| IHT | O(MN) | Iterative hard thresholding |
| SP | O(KMN) | Subspace pursuit |

### 2.3 LASSO (L1-Regularized Least Squares)

$$\min_{\mathbf{s}} \frac{1}{2}\|\mathbf{A}\mathbf{s} - \mathbf{y}\|_2^2 + \lambda \|\mathbf{s}\|_1$$

**Algorithms:**
- Coordinate Descent
- LARS (Least Angle Regression)
- ADMM

---

## 3. Measurement Matrix Design

### 3.1 Random Matrices

| Type | Construction | RIP? |
|------|--------------|------|
| Gaussian | $\Phi_{ij} \sim \mathcal{N}(0, 1/M)$ | Yes, w.h.p. |
| Bernoulli | $\Phi_{ij} \in \{-1/\sqrt{M}, +1/\sqrt{M}\}$ | Yes, w.h.p. |
| Subsampled Fourier | Random rows of DFT | Yes, w.h.p. |

### 3.2 Structured Matrices

For hardware implementation:
- Sparse binary matrices
- Toeplitz/Circulant matrices
- Hadamard matrices

### 3.3 Minimum Measurements

| Sparsity K | Signal length N | Min M (theory) | M (practical) |
|------------|-----------------|----------------|---------------|
| 5 | 256 | ~40 | ~50-60 |
| 10 | 256 | ~80 | ~100-120 |
| 20 | 1024 | ~180 | ~200-250 |

Rule of thumb: $M \approx 4K \log(N/K)$

---

## 4. Demo Files

| File | Dimension | Scenario | Key Concepts |
|------|-----------|----------|--------------|
| `demo_1d_sparse_recovery.m` | 1D | Sparse spectrum recovery | OMP, L1-minimization, RIP |
| `demo_2d_image_cs.m` | 2D | Image reconstruction | Random sampling, TV regularization |
| `app_mri_acceleration.m` | Application | MRI k-space undersampling | Variable density sampling |
| `app_radar_imaging.m` | Application | ISAR sparse aperture | Range-Doppler imaging |

---

## 5. Quick Start

```matlab
% Add paths
addpath(genpath('../../utils'));

% Run 1D demo: Recover sparse frequency spectrum
demo_1d_sparse_recovery

% Run 2D demo: Image CS reconstruction
demo_2d_image_cs

% Run MRI acceleration application
app_mri_acceleration

% Run radar imaging application
app_radar_imaging
```

---

## 6. When to Use Compressed Sensing?

### ✅ Good Scenarios

- Signal is sparse in known domain (frequency, wavelet, gradient)
- Physical sampling is expensive (MRI scan time, radar aperture)
- Power/bandwidth limited (IoT sensors)
- Random access to samples is possible

### ❌ Bad Scenarios

- Signal is not sparse
- Real-time processing required (recovery is slow)
- Measurement matrix must be deterministic with bad coherence
- Very high SNR requirements

---

## 7. Comparison with Other Methods

| Aspect | Compressed Sensing | CRT Multi-rate | MUSIC/ESPRIT |
|--------|-------------------|----------------|--------------|
| **Assumption** | Sparsity | Multiple rates | Finite sinusoids |
| **Hardware** | Random sampler | Multiple ADCs | Single ADC |
| **Computation** | High (optimization) | Low (CRT) | Medium (SVD) |
| **Robustness** | Moderate | Sensitive to phase | Very sensitive |
| **Max sparsity** | ~M/(4 log N) | N/A | ~N/3 |

---

## 8. References

1. Candès, E. J., Romberg, J., & Tao, T. (2006). "Robust uncertainty principles: exact signal reconstruction from highly incomplete frequency information." *IEEE TIT*.

2. Donoho, D. L. (2006). "Compressed sensing." *IEEE TIT*.

3. Baraniuk, R. G. (2007). "Compressive sensing." *IEEE Signal Processing Magazine*.

4. Tropp, J. A., & Gilbert, A. C. (2007). "Signal recovery from random measurements via orthogonal matching pursuit." *IEEE TIT*.

5. Lustig, M., Donoho, D., & Pauly, J. M. (2007). "Sparse MRI: The application of compressed sensing for rapid MR imaging." *Magnetic Resonance in Medicine*.

---

## 9. Further Reading

- **Total Variation (TV) regularization** for piecewise constant signals
- **Dictionary Learning** for adaptive sparsifying bases
- **Structured Sparsity** (group sparsity, block sparsity)
- **One-bit Compressed Sensing** for extreme quantization
- **Deep Unfolding** combining CS with neural networks
