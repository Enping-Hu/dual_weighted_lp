# Dual Weighted $L_p$ — Low-Light Image Enhancement

Official source code for the paper:

> **Low-Light Image Enhancement Using a Retinex-Based Variational Model With Weighted $L_p$ Norm Constraint**
> *IEEE Transactions on Multimedia*, vol. 28, pp. 152–166, 2026.
> DOI: [10.1109/TMM.2025.3618567](https://doi.org/10.1109/TMM.2025.3618567)

[English](#english) · [中文](#中文)

---

## English

### Overview

Images captured under low-light conditions usually suffer from limited visibility, low
contrast, and severe noise, which degrades downstream computer-vision tasks. Most
Retinex-based variational methods constrain the illumination and reflectance with
*integer-order* norms ($L_0$, $L_1$, or $L_2$), which are tied to specific statistical
distributions (e.g. $L_2$↔Gaussian, $L_1$↔Laplacian) and often fail to match the
complexity of real scenes, leading to imperfect Retinex decomposition.

This work proposes a **unified Retinex-based variational model with flexible weighted
$L_p$ norm constraints**. By allowing $p$ to vary, the model simultaneously enforces
piece-wise smoothness on the illumination and rich textures on the reflectance, while an
explicit noise term suppresses hidden noise without amplification. The result is a
**structure-aware illumination** and a **detail-revealed reflectance**.

### Key contributions

- **Flexible $L_p$ constraints on both layers.** A unified model applies $L_p$ norms with
  different $p$ ranges to illumination ($0 \le p_L \le 1$, for piece-wise smoothness) and
  reflectance ($1 \le p_R \le 2$, for abundant texture), going beyond fixed integer-order
  norms.
- **Two pixel-wise weight matrices.** A variance-based weight $W_L$ preserves the
  structural edges of the illumination, and a gradient-based weight $W_R$ retains fine
  details in the reflectance.
- **Explicit noise estimation.** An $L_2$ term models the overall noise level, so the
  method predicts a noise map while jointly estimating a smooth illumination and a
  textured reflectance — avoiding noise amplification.
- **Strong results.** Competitive qualitative and quantitative performance on real and
  synthetic datasets, with the model also extending to image dehazing and underwater
  enhancement.

### Method

Under Retinex theory, the captured image $S$ is decomposed into illumination $L$,
reflectance $R$, and a noise term $N$:

$$ S = L \circ R + N $$

where $\circ$ denotes element-wise multiplication. The enhancement is formulated as the
energy minimization (Eq. 3 in the paper):

$$
\arg\min_{L,R,N}\; \lVert L\circ R + N - S\rVert_2^2
\;+\; \lambda_L\,\lVert W_L\circ\nabla L\rVert_{p_L}
\;+\; \lambda_R\,\lVert W_R\circ\nabla R\rVert_{p_R}
\;+\; \alpha\,\lVert N\rVert_2^2
\;+\; \beta\,\lVert L - B\rVert_2^2
$$

- **Fidelity** $\lVert L\circ R + N - S\rVert_2^2$ — keeps the recomposed image close to the input.
- **Illumination prior** $\lambda_L\lVert W_L\circ\nabla L\rVert_{p_L}$ — weighted $L_{p_L}$ norm for a piece-wise smooth, edge-preserving illumination.
- **Reflectance prior** $\lambda_R\lVert W_R\circ\nabla R\rVert_{p_R}$ — weighted $L_{p_R}$ norm to reveal textures and details.
- **Noise term** $\alpha\lVert N\rVert_2^2$ — controls overall noise intensity.
- **Bright-channel term** $\beta\lVert L - B\rVert_2^2$ — drives $L$ toward an initial illumination $B$ (bright-channel prior, max over a $4\times4$ window).

The weight functions are:

$$ W_L(x) = e^{-\delta\,\lvert g_{r,\sigma}(S^2(x)) - (g_{r,\sigma}(S(x)))^2 \rvert} $$

$$ W_R(x) = \left(\max\!\left(\Big(\nabla S(x) - \tfrac{1}{\lvert\Omega\rvert}\textstyle\sum_\Omega \nabla S(x)\Big)^2 + \Big\lvert \tfrac{1}{\lvert\Omega\rvert}\textstyle\sum_\Omega \nabla S(x)\Big\rvert^{0.5},\, \varepsilon\right)\right)^{-1} $$

where $g_{r,\sigma}$ is a Gaussian filter and $\varepsilon$ is a small constant.

**Solver.** Because the $L_p$ regularizers make the objective non-convex, it is optimized
with **block coordinate descent** combined with an **iteratively reweighted least squares
(IRLS)** reformulation of the $L_p$ terms. The $L$, $R$, and $N$ sub-problems have
closed-form / fast solutions (the $L$ sub-problem uses a preconditioned conjugate
gradient solver). For the special case $p_L = 0$, the $L$ sub-problem is solved via
**ADMM with FFT** and an $L_0$ shrinkage step.

**Final enhancement.** The algorithm runs on the **V channel in HSV color space**, then
applies Gamma correction ($\gamma = 2.2$) to the estimated illumination and recombines:

$$ \hat{S} = \hat{L}^{1/\gamma} \circ \hat{R} $$

### Repository structure

```
dual_weighted_lp/
├── demo.m                 # Entry-point demo script
├── dual_weighted_lp.p     # Core algorithm (P-coded / protected MATLAB)
└── Third_codes/           # Third-party helper functions
    ├── boxfilter.m
    ├── guidedfilter.m     # Guided filter
    ├── convBox.m
    ├── convMax.m
    ├── convConst.cpp      # Fast convolution (Piotr's Toolbox)
    ├── convConst.mexa64   # Compiled MEX (Linux)
    ├── convConst.mexw64   # Compiled MEX (Windows)
    ├── sse.hpp
    └── wrappers.hpp
```

> **Note:** the core method `dual_weighted_lp.p` is a P-coded MATLAB file; it can be
> called but its source is not readable.

### Requirements

- MATLAB (developed with **R2021a**; Image Processing Toolbox required).
- The repository ships pre-compiled MEX binaries (`convConst.mexa64`, `convConst.mexw64`).
  If neither matches your platform, recompile from `convConst.cpp`:

  ```matlab
  mex Third_codes/convConst.cpp
  ```

### Usage

```matlab
% in demo.m
filename  = 'path/to/your/low_light_image.png';   % change to your input
outputDir = 'path/to/output';                      % change to your output folder
im = im2double(imread(filename));

addpath(genpath('Third_codes'));

% Retinex decomposition: illumination L, reflectance R, noise N
[L, R, N] = dual_weighted_lp(im);

% Enhancement on the V channel of HSV with gamma correction
gamma   = 2.2;
hsv     = rgb2hsv(im);
L_gamma = L.^(1/gamma);
S_gamma = R .* L_gamma;
hsv(:,:,3) = real(S_gamma);
enhance = hsv2rgb(hsv);

imwrite(enhance, fullfile(outputDir, 'E.png'));
```

> The original `demo.m` uses hard-coded Windows paths
> (e.g. `D:\our\dataset\LOL\1.png`). Edit `filename` and `outputDir` to your own paths
> before running.

### Default parameters

| Parameter | Symbol | Default | Role |
|-----------|--------|---------|------|
| Illumination weight | $\lambda_L$ | `0.1` ($10^{-1}$) | Smoothness strength of illumination |
| Reflectance weight | $\lambda_R$ | `1e-6` ($10^{-6}$) | Detail strength of reflectance |
| Noise weight | $\alpha$ | `0.1` ($10^{-1}$) | Overall noise level control |
| Bright-channel weight | $\beta$ | `2` | Fidelity to initial illumination |
| Illumination norm | $p_L$ | `0.5` | $L_{p_L}$ order for illumination |
| Reflectance norm | $p_R$ | `1.5` | $L_{p_R}$ order for reflectance |
| Weight sharpness | $\delta$ | `100` | Sharpness in $W_L$ |
| Max iterations | $K$ | `15` | Stopping criterion |

Smaller $\lambda_L$ → smoother illumination; smaller $\lambda_R$ → richer reflectance
detail; larger $\alpha$ → more noise retained in the noise map. The defaults work well
across datasets without per-image tuning.

### Datasets & evaluation

Evaluated on **DICM** (64), **LIME** (10), **MEF** (17), **NPE** (85), and **LOL** (500).
Quality is measured with three no-reference metrics — **NIQE**, **ARISMC**, and
**BRISQUE** (lower is better). The method ranks first on NIQE across three real-world
datasets and is consistently top-two on ARISMC and BRISQUE, compared against MSRCR,
NPEA, WVM, LIME, JIEP, PnPR, RMSP, FMSTD, Zero-DCE, RUAS, PIE, and Low-FaceNet.

### Citation

```bibtex
@article{hu2026lowlight,
  title   = {Low-Light Image Enhancement Using a Retinex-Based Variational Model With Weighted $L_p$ Norm Constraint},
  author  = {Hu, Enping and Liu, Yun and Wang, Anzhi and Shiri, Babak and Ren, Wenqi and Lin, Weisi},
  journal = {IEEE Transactions on Multimedia},
  volume  = {28},
  pages   = {152--166},
  year    = {2026},
  doi     = {10.1109/TMM.2025.3618567}
}
```

---

## 中文

### 简介

低照度条件下拍摄的图像往往可见度低、对比度差、噪声严重,会影响后续的计算机视觉任务。大多数
基于 Retinex 的变分方法使用**整数阶范数**($L_0$、$L_1$、$L_2$)来约束光照层和反射层,而
这些范数各自对应特定的统计分布(如 $L_2$↔高斯、$L_1$↔拉普拉斯),难以匹配真实场景的复杂
性,从而无法得到理想的 Retinex 分解。

本工作提出一个**统一的、带可变加权 $L_p$ 范数约束的 Retinex 变分模型**。通过让 $p$ 取不同
的值,模型可同时让光照层保持分片平滑、让反射层呈现丰富纹理;同时引入显式噪声项抑制隐藏噪声、
避免噪声放大,最终得到**结构感知的光照层**与**细节丰富的反射层**。

### 主要贡献

- **对两个分量都使用灵活的 $L_p$ 约束。** 统一模型对光照层($0 \le p_L \le 1$,保证分片平
  滑)和反射层($1 \le p_R \le 2$,保留丰富纹理)分别施加不同 $p$ 区间的 $L_p$ 范数,突破
  了固定整数阶范数的局限。
- **两个像素级权重矩阵。** 基于方差的权重 $W_L$ 用于保持光照层的结构边缘;基于梯度的权重
  $W_R$ 用于保留反射层的精细细节。
- **显式噪声估计。** 用 $L_2$ 项对整体噪声水平建模,在估计平滑光照与纹理反射的同时预测噪声
  图,从而避免噪声放大。
- **效果优越。** 在真实与合成数据集上均取得有竞争力的主客观表现,模型还可推广到图像去雾与水
  下图像增强。

### 方法

在 Retinex 理论下,采集图像 $S$ 被分解为光照层 $L$、反射层 $R$ 与噪声项 $N$:

$$ S = L \circ R + N $$

其中 $\circ$ 表示逐元素相乘。增强被建模为如下能量最小化问题(论文式 3):

$$
\arg\min_{L,R,N}\; \lVert L\circ R + N - S\rVert_2^2
\;+\; \lambda_L\,\lVert W_L\circ\nabla L\rVert_{p_L}
\;+\; \lambda_R\,\lVert W_R\circ\nabla R\rVert_{p_R}
\;+\; \alpha\,\lVert N\rVert_2^2
\;+\; \beta\,\lVert L - B\rVert_2^2
$$

- **保真项** $\lVert L\circ R + N - S\rVert_2^2$ —— 让重组图像贴近输入。
- **光照先验** $\lambda_L\lVert W_L\circ\nabla L\rVert_{p_L}$ —— 加权 $L_{p_L}$ 范数,得到分片平滑、保边的光照层。
- **反射先验** $\lambda_R\lVert W_R\circ\nabla R\rVert_{p_R}$ —— 加权 $L_{p_R}$ 范数,凸显纹理与细节。
- **噪声项** $\alpha\lVert N\rVert_2^2$ —— 控制整体噪声强度。
- **亮通道项** $\beta\lVert L - B\rVert_2^2$ —— 让 $L$ 逼近初始光照 $B$($B$ 由亮通道先验得到,在 $4\times4$ 窗口内取最大值)。

权重函数为:

$$ W_L(x) = e^{-\delta\,\lvert g_{r,\sigma}(S^2(x)) - (g_{r,\sigma}(S(x)))^2 \rvert} $$

$$ W_R(x) = \left(\max\!\left(\Big(\nabla S(x) - \tfrac{1}{\lvert\Omega\rvert}\textstyle\sum_\Omega \nabla S(x)\Big)^2 + \Big\lvert \tfrac{1}{\lvert\Omega\rvert}\textstyle\sum_\Omega \nabla S(x)\Big\rvert^{0.5},\, \varepsilon\right)\right)^{-1} $$

其中 $g_{r,\sigma}$ 为高斯滤波器,$\varepsilon$ 为防止除零的小常数。

**求解。** 由于 $L_p$ 正则项使目标函数非凸,采用**块坐标下降**结合对 $L_p$ 项的**迭代重加权
最小二乘(IRLS)**重构来优化。$L$、$R$、$N$ 三个子问题均有闭式解或快速解($L$ 子问题用预
条件共轭梯度求解)。当 $p_L = 0$ 这一特殊情形时,$L$ 子问题改用 **ADMM + FFT** 并配合
$L_0$ 收缩(shrinkage)步骤求解。

**最终增强。** 算法在 **HSV 色彩空间的 V 通道**上进行,对估计出的光照层做 $\gamma = 2.2$
的 Gamma 校正后重组:

$$ \hat{S} = \hat{L}^{1/\gamma} \circ \hat{R} $$

### 仓库结构

```
dual_weighted_lp/
├── demo.m                 # 演示入口脚本
├── dual_weighted_lp.p     # 核心算法(P-code 加密的 MATLAB 文件)
└── Third_codes/           # 第三方辅助函数
    ├── boxfilter.m
    ├── guidedfilter.m     # 导向滤波
    ├── convBox.m
    ├── convMax.m
    ├── convConst.cpp      # 快速卷积(Piotr's Toolbox)
    ├── convConst.mexa64   # 已编译 MEX(Linux)
    ├── convConst.mexw64   # 已编译 MEX(Windows)
    ├── sse.hpp
    └── wrappers.hpp
```

> **说明:** 核心方法 `dual_weighted_lp.p` 是 P-code 加密文件,可以调用但源码不可读。

### 运行环境

- MATLAB(开发环境为 **R2021a**;需要 Image Processing Toolbox)。
- 仓库已附带预编译的 MEX 二进制(`convConst.mexa64`、`convConst.mexw64`)。若都与你的平台不
  匹配,可从 `convConst.cpp` 重新编译:

  ```matlab
  mex Third_codes/convConst.cpp
  ```

### 用法

```matlab
% demo.m 中
filename  = 'path/to/your/low_light_image.png';   % 改成你的输入图像
outputDir = 'path/to/output';                      % 改成你的输出目录
im = im2double(imread(filename));

addpath(genpath('Third_codes'));

% Retinex 分解:光照 L、反射 R、噪声 N
[L, R, N] = dual_weighted_lp(im);

% 在 HSV 的 V 通道上做 gamma 校正增强
gamma   = 2.2;
hsv     = rgb2hsv(im);
L_gamma = L.^(1/gamma);
S_gamma = R .* L_gamma;
hsv(:,:,3) = real(S_gamma);
enhance = hsv2rgb(hsv);

imwrite(enhance, fullfile(outputDir, 'E.png'));
```

> 原始 `demo.m` 使用了硬编码的 Windows 路径(如 `D:\our\dataset\LOL\1.png`)。运行前请把
> `filename` 与 `outputDir` 改成你自己的路径。

### 默认参数

| 参数 | 符号 | 默认值 | 作用 |
|------|------|--------|------|
| 光照权重 | $\lambda_L$ | `0.1`($10^{-1}$) | 光照层平滑强度 |
| 反射权重 | $\lambda_R$ | `1e-6`($10^{-6}$) | 反射层细节强度 |
| 噪声权重 | $\alpha$ | `0.1`($10^{-1}$) | 整体噪声水平控制 |
| 亮通道权重 | $\beta$ | `2` | 对初始光照的保真度 |
| 光照范数 | $p_L$ | `0.5` | 光照层的 $L_{p_L}$ 阶数 |
| 反射范数 | $p_R$ | `1.5` | 反射层的 $L_{p_R}$ 阶数 |
| 权重锐度 | $\delta$ | `100` | $W_L$ 中的锐度系数 |
| 最大迭代 | $K$ | `15` | 停止条件 |

$\lambda_L$ 越小光照越平滑;$\lambda_R$ 越小反射细节越丰富;$\alpha$ 越大噪声图中保留的噪声
越多。默认参数无需逐图调节即可在各数据集上取得良好效果。

### 数据集与评测

在 **DICM**(64 张)、**LIME**(10 张)、**MEF**(17 张)、**NPE**(85 张)、**LOL**(500
张)上评测,采用三个无参考指标 **NIQE**、**ARISMC**、**BRISQUE**(数值越小越好)。在三个真
实数据集上 NIQE 排名第一,ARISMC 与 BRISQUE 多数情况下位列前二;对比方法包括 MSRCR、NPEA、
WVM、LIME、JIEP、PnPR、RMSP、FMSTD、Zero-DCE、RUAS、PIE、Low-FaceNet。

### 引用

```bibtex
@article{hu2026lowlight,
  title   = {Low-Light Image Enhancement Using a Retinex-Based Variational Model With Weighted $L_p$ Norm Constraint},
  author  = {Hu, Enping and Liu, Yun and Wang, Anzhi and Shiri, Babak and Ren, Wenqi and Lin, Weisi},
  journal = {IEEE Transactions on Multimedia},
  volume  = {28},
  pages   = {152--166},
  year    = {2026},
  doi     = {10.1109/TMM.2025.3618567}
}
```
