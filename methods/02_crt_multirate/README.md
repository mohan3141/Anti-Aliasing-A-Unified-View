# Chinese Remainder Theorem for Anti-Aliasing

中国剩余定理（CRT）在欠采样信号恢复中的应用

## 1. 问题背景

当采样率 $F_s$ 不满足Nyquist条件时，真实频率 $f_{true}$ 会混叠到：

$$f_{alias} = ((f_{true} + F_s/2) \mod F_s) - F_s/2$$

**核心问题**：从 $f_{alias}$ 能否恢复 $f_{true}$？

单一采样率下：**不能**（无穷多个频率混叠到同一位置）

多采样率下：**可以**（利用CRT）

## 2. 数学原理

### 2.1 中国剩余定理

若 $m_1, m_2$ 互质，则对于任意整数 $a_1, a_2$，同余方程组：

$$x \equiv a_1 \pmod{m_1}$$
$$x \equiv a_2 \pmod{m_2}$$

在 $[0, m_1 \cdot m_2)$ 内有**唯一解**。

### 2.2 应用于频率估计

设两个采样率 $F_{s1}, F_{s2}$，观测到的混叠频率为 $f_1, f_2$：

$$f_{true} \equiv f_1 \pmod{F_{s1}}$$
$$f_{true} \equiv f_2 \pmod{F_{s2}}$$

**唯一恢复条件**：$f_{true} < \text{LCM}(F_{s1}, F_{s2})$

### 2.3 最优采样率选择

为最大化恢复范围，选择**互质**的采样率：

| $F_{s1}$ | $F_{s2}$ | LCM | 可恢复范围 |
|----------|----------|-----|-----------|
| 100 Hz   | 110 Hz   | 1100 Hz | 0-1100 Hz |
| 100 Hz   | 111 Hz   | 11100 Hz | 0-11100 Hz |
| 97 Hz    | 101 Hz   | 9797 Hz | 0-9797 Hz |

**技巧**：选择相邻质数作为采样率比例

## 3. Demo 列表

| 文件 | 维度 | 场景 |
|------|------|------|
| `demo_1d_frequency_est.m` | 1D | 双ADC频率估计 |
| `demo_2d_phase_unwrap.m` | 2D | InSAR相位解缠绕 |
| `app_tof_ranging.m` | 应用 | ToF相机测距解模糊 |

## 4. 运行示例

```matlab
% 1D频率估计
demo_1d_frequency_est

% 2D相位解缠绕
demo_2d_phase_unwrap

% ToF测距应用
app_tof_ranging
```

## 5. 参考文献

1. Xia, X. G., & Wang, G. (2007). Phase unwrapping and a robust Chinese remainder theorem. *IEEE Signal Processing Letters*.
2. Li, X., et al. (2019). Generalized CRT for frequency estimation. *IEEE Trans. Signal Processing*.
