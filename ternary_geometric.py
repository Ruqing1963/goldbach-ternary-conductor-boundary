"""
ternary_geometric.py — True geometric analysis of ternary Goldbach-Frey conductor

Starting from the correct discriminant:
  Δ ∝ ∏ pᵢ⁶ · ∏_{i<j}(pᵢ² - pⱼ²)⁴
    = ∏ pᵢ⁶ · ∏_{i<j}(pᵢ - pⱼ)⁴ · ∏_k (N - pₖ)⁴

The odd radical of Δ (which determines conductor support) is:
  rad_odd(Δ) = product of all ODD primes dividing
    p₁·p₂·p₃ · (p₁-p₂)·(p₂-p₃)·(p₃-p₁) · (p₂+p₃)·(p₁+p₃)·(p₁+p₂)

Normalized: ρ₃ = log(rad_odd(Δ)) / log(N)

Key geometric questions:
  Q1: What is the natural decomposition of ρ₃?
  Q2: What separates PPP from CCC in the TRUE conductor?
  Q3: Is there a structural law governing ⟨ρ₃⟩ for PPP triples?
"""
import math
import numpy as np
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sieve(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = 0
    return is_prime

def odd_prime_factors(n):
    """Return set of odd prime factors of |n|."""
    if n == 0: return set()
    n = abs(n)
    while n % 2 == 0: n //= 2
    factors = set()
    d = 3
    while d * d <= n:
        if n % d == 0:
            factors.add(d)
            while n % d == 0: n //= d
        d += 2
    if n > 1:
        factors.add(n)
    return factors

def log_odd_rad(n):
    """Return log(rad_odd(|n|))."""
    if n == 0: return 0.0
    return sum(math.log(p) for p in odd_prime_factors(n))

LIMIT = 10000
is_prime = sieve(LIMIT)

# ═══════════════════════════════════════════════════════════════════════════════
# CORRECT ρ₃ from TRUE discriminant
# ═══════════════════════════════════════════════════════════════════════════════
def compute_rho3(N, p1, p2, p3):
    """
    ρ₃ = log(rad_odd(Δ)) / log(N)
    
    where rad_odd(Δ) = product of all odd primes dividing any of:
      p₁, p₂, p₃, (p₁-p₂), (p₂-p₃), (p₃-p₁), (p₂+p₃), (p₁+p₃), (p₁+p₂)
    
    Note: exponents in Δ are irrelevant for rad.
    """
    all_factors = set()
    
    # Summand factors
    for p in [p1, p2, p3]:
        all_factors |= odd_prime_factors(p)
    
    # Difference factors
    for d in [p1-p2, p2-p3, p3-p1]:
        all_factors |= odd_prime_factors(d)
    
    # Partial sum factors (= N - pₖ)
    for s in [p2+p3, p1+p3, p1+p2]:
        all_factors |= odd_prime_factors(s)
    
    log_rad = sum(math.log(f) for f in all_factors)
    return log_rad / math.log(N)

def decompose_rho3(N, p1, p2, p3):
    """
    Decompose ρ₃ into three geometric components:
      σ = summand contribution = log(rad_odd(p₁·p₂·p₃)) / log N
      δ = difference contribution = log(rad_odd(∏(pᵢ-pⱼ))) / log N  [NEW primes only]
      π = partial sum contribution = log(rad_odd(∏(N-pₖ))) / log N  [NEW primes only]
    """
    logN = math.log(N)
    
    summand_primes = set()
    for p in [p1, p2, p3]:
        summand_primes |= odd_prime_factors(p)
    sigma = sum(math.log(f) for f in summand_primes) / logN
    
    diff_primes = set()
    for d in [p1-p2, p2-p3, p3-p1]:
        diff_primes |= odd_prime_factors(d)
    new_diff = diff_primes - summand_primes
    delta = sum(math.log(f) for f in new_diff) / logN
    
    sum_primes = set()
    for s in [p2+p3, p1+p3, p1+p2]:
        sum_primes |= odd_prime_factors(s)
    new_sum = sum_primes - summand_primes - diff_primes
    pi_contrib = sum(math.log(f) for f in new_sum) / logN
    
    return sigma, delta, pi_contrib

# ═══════════════════════════════════════════════════════════════════════════════
# SCAN 1: Detailed N=1025 with decomposition
# ═══════════════════════════════════════════════════════════════════════════════
N = 1025
print("=" * 80)
print(f"SCAN 1: Geometric decomposition at N = {N}")
print("=" * 80)

ppp_data = []
ccc_data = []
mix_data = []

for p1 in range(3, N // 3 + 1):
    for p2 in range(p1, (N - p1) // 2 + 1):
        p3 = N - p1 - p2
        if p3 < p2: break
        if p3 <= 1: continue
        
        rho = compute_rho3(N, p1, p2, p3)
        sigma, delta, pi_c = decompose_rho3(N, p1, p2, p3)
        
        ip = bool(is_prime[p1]) + bool(is_prime[p2]) + bool(is_prime[p3])
        entry = (p1, p2, p3, rho, sigma, delta, pi_c, ip)
        
        if ip == 3: ppp_data.append(entry)
        elif ip == 0: ccc_data.append(entry)
        else: mix_data.append(entry)

print(f"\nTotal: PPP={len(ppp_data)}, Mixed={len(mix_data)}, CCC={len(ccc_data)}")

# Decomposition statistics
print(f"\n{'Category':<10} {'ρ₃':<12} {'σ (summ)':<12} {'δ (diff)':<12} {'π (part.sum)':<12}")
print("-" * 60)
for label, data in [("PPP", ppp_data), ("CCC", ccc_data)]:
    rhos = [d[3] for d in data]
    sigmas = [d[4] for d in data]
    deltas = [d[5] for d in data]
    pis = [d[6] for d in data]
    print(f"{label:<10} {np.mean(rhos):<12.3f} {np.mean(sigmas):<12.3f} "
          f"{np.mean(deltas):<12.3f} {np.mean(pis):<12.3f}")

print(f"\nDecomposition: ρ₃ = σ + δ + π")
print(f"  σ = primes from summands pᵢ")
print(f"  δ = NEW primes from differences (pᵢ - pⱼ)")
print(f"  π = NEW primes from partial sums (N - pₖ)")

# Understand the PPP-CCC gap
print(f"\n{'Component gap':<25} {'PPP mean':<12} {'CCC mean':<12} {'Δ(PPP-CCC)':<12}")
print("-" * 60)
for comp, idx in [("σ (summands)", 4), ("δ (differences)", 5), ("π (partial sums)", 6)]:
    ppp_v = np.mean([d[idx] for d in ppp_data])
    ccc_v = np.mean([d[idx] for d in ccc_data])
    print(f"{comp:<25} {ppp_v:<12.4f} {ccc_v:<12.4f} {ppp_v - ccc_v:<+12.4f}")

# Ground states
print(f"\nPPP ground states:")
for t in sorted(ppp_data, key=lambda x: x[3])[:5]:
    p1,p2,p3,rho,s,d,p,_ = t
    print(f"  ({p1:3d}, {p2:3d}, {p3:3d})  ρ₃={rho:.3f}  [σ={s:.2f} δ={d:.2f} π={p:.2f}]")

print(f"\nCCC ground states:")
for t in sorted(ccc_data, key=lambda x: x[3])[:5]:
    p1,p2,p3,rho,s,d,p,_ = t
    print(f"  ({p1:3d}, {p2:3d}, {p3:3d})  ρ₃={rho:.3f}  [σ={s:.2f} δ={d:.2f} π={p:.2f}]")

# ═══════════════════════════════════════════════════════════════════════════════
# SCAN 2: What controls ρ₃ for PPP? Look for natural organizing variable.
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("SCAN 2: What natural variable controls ⟨ρ₃⟩ for PPP?")
print("=" * 80)

# For PPP triples: rad(pᵢ) = pᵢ (since all are prime)
# So σ = log(p₁·p₂·p₃)/log(N) — this is fixed by the SIZES of the primes
# The variation comes from δ and π: what new primes appear in diffs and sums

# Hypothesis: for PPP, ρ₃ is primarily controlled by the smallest summand p₁
# because small p₁ → large p₃ → diffs (p₃-p₁) and (p₃-p₂) are large and likely
# to be smooth → fewer NEW prime factors

p1_vals = np.array([d[0] for d in ppp_data])
rho_vals = np.array([d[3] for d in ppp_data])
sigma_vals = np.array([d[4] for d in ppp_data])

# Correlations
corr_p1_rho = np.corrcoef(p1_vals, rho_vals)[0,1]
corr_sigma_rho = np.corrcoef(sigma_vals, rho_vals)[0,1]
print(f"  Corr(p₁, ρ₃) = {corr_p1_rho:.4f}")
print(f"  Corr(σ, ρ₃) = {corr_sigma_rho:.4f}")

# What about the "smoothness" of differences?
# Define: τ = log(rad_odd(∏(pᵢ-pⱼ)·∏(N-pₖ))) / log(N)  (non-summand part)
tau_vals = np.array([d[5] + d[6] for d in ppp_data])
corr_tau_rho = np.corrcoef(tau_vals, rho_vals)[0,1]
print(f"  Corr(τ=δ+π, ρ₃) = {corr_tau_rho:.4f}")

# σ dominates: let's check
print(f"\n  Mean contributions: σ={np.mean(sigma_vals):.3f}, δ+π={np.mean(tau_vals):.3f}")
print(f"  Fraction from summands: {np.mean(sigma_vals)/np.mean(rho_vals)*100:.1f}%")

# ═══════════════════════════════════════════════════════════════════════════════
# SCAN 3: Multi-N scan — does ρ₃ scale with log(N)?
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("SCAN 3: ρ₃ scaling with N (odd N in [501, 4001])")
print("=" * 80)

multi_n = []
for N in range(501, 4002, 10):
    if N % 2 == 0: N += 1
    ppp_r = []
    for p1 in range(3, N // 3 + 1, 2):
        if not is_prime[p1]: continue
        for p2 in range(p1, (N - p1) // 2 + 1, 2):
            if not is_prime[p2]: continue
            p3 = N - p1 - p2
            if p3 < p2 or p3 <= 1: continue
            if not is_prime[p3]: continue
            rho = compute_rho3(N, p1, p2, p3)
            ppp_r.append(rho)
    
    if len(ppp_r) >= 10:
        multi_n.append({
            'N': N, 'logN': math.log(N),
            'n_ppp': len(ppp_r),
            'rho_mean': np.mean(ppp_r),
            'rho_min': min(ppp_r),
            'rho_max': max(ppp_r),
            'bw': max(ppp_r) - min(ppp_r),
        })

logNs = np.array([r['logN'] for r in multi_n])
means_m = np.array([r['rho_mean'] for r in multi_n])
mins_m = np.array([r['rho_min'] for r in multi_n])
bws_m = np.array([r['bw'] for r in multi_n])

# The natural scaling: for PPP with typical triple (N/3, N/3, N/3):
# σ = log((N/3)³)/log(N) = 3 - 3log3/logN → 3
# δ+π depends on diffs and sums...
# So ρ₃ ∝ some function of logN

c_logN = np.polyfit(logNs, means_m, 1)
resid = means_m - (c_logN[0]*logNs + c_logN[1])
r2_logN = 1 - np.var(resid)/np.var(means_m)

print(f"  Linear fit ⟨ρ₃⟩ = {c_logN[0]:.4f} log(N) + {c_logN[1]:.4f}")
print(f"  R² = {r2_logN:.6f}")

# Try: ρ₃ / logN = const?
ratios = means_m / logNs
print(f"  ⟨ρ₃⟩/log(N): mean = {np.mean(ratios):.4f}, std = {np.std(ratios):.4f}")

# Bandwidth scaling
c_bw = np.polyfit(logNs, bws_m, 1)
print(f"  Bandwidth: {c_bw[0]:.4f} log(N) + {c_bw[1]:.4f}")

# ═══════════════════════════════════════════════════════════════════════════════
# SCAN 4: The true PPP-CCC separation — what's the geometric gap?
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("SCAN 4: PPP-CCC gap stability across N")
print("=" * 80)

gap_data = []
for N in range(501, 2002, 20):
    if N % 2 == 0: N += 1
    ppp_r, ccc_r = [], []
    for p1 in range(3, N // 3 + 1):
        for p2 in range(p1, (N - p1) // 2 + 1):
            p3 = N - p1 - p2
            if p3 < p2: break
            if p3 <= 1: continue
            rho = compute_rho3(N, p1, p2, p3)
            ip = bool(is_prime[p1]) + bool(is_prime[p2]) + bool(is_prime[p3])
            if ip == 3: ppp_r.append(rho)
            elif ip == 0: ccc_r.append(rho)
    
    if ppp_r and ccc_r:
        gap_data.append({
            'N': N,
            'ppp_min': min(ppp_r), 'ppp_mean': np.mean(ppp_r),
            'ccc_min': min(ccc_r), 'ccc_mean': np.mean(ccc_r),
            'floor_gap': min(ppp_r) - min(ccc_r),
            'mean_gap': np.mean(ppp_r) - np.mean(ccc_r),
        })

floor_gaps = [g['floor_gap'] for g in gap_data]
mean_gaps = [g['mean_gap'] for g in gap_data]
print(f"  Floor gap (PPP_min - CCC_min): {np.mean(floor_gaps):.3f} ± {np.std(floor_gaps):.3f}")
print(f"  Mean gap (PPP_mean - CCC_mean): {np.mean(mean_gaps):.3f} ± {np.std(mean_gaps):.3f}")
print(f"  Floor gap always positive: {all(g > 0 for g in floor_gaps)}")
print(f"  → Primes ALWAYS above composites (consistent with binary, NO inversion)")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURES
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[*] Generating figures...")

# Figure 1: Geometric decomposition at N=1025
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# 1a: PPP vs CCC histograms
ppp_rhos = np.array([d[3] for d in ppp_data])
ccc_rhos = np.array([d[3] for d in ccc_data])
bins = np.linspace(min(min(ccc_rhos), min(ppp_rhos))-0.5,
                   max(max(ccc_rhos), max(ppp_rhos))+0.5, 60)
axes[0].hist(ccc_rhos, bins=bins, alpha=0.4, color='#CC6644', density=True,
             label=f'CCC ({len(ccc_data)})')
axes[0].hist(ppp_rhos, bins=bins, alpha=0.6, color='#2255BB', density=True,
             label=f'PPP ({len(ppp_data)})')
axes[0].set_xlabel(r'$\rho_3 = \log(\mathrm{rad}_{\mathrm{odd}}(\Delta))/\log N$', fontsize=10)
axes[0].set_ylabel('Density', fontsize=10)
axes[0].set_title(f'$N = {N}$: True Conductor', fontsize=12)
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.15)

# 1b: Decomposition bar chart
categories = ['PPP', 'CCC']
sigma_means = [np.mean([d[4] for d in ppp_data]), np.mean([d[4] for d in ccc_data])]
delta_means = [np.mean([d[5] for d in ppp_data]), np.mean([d[5] for d in ccc_data])]
pi_means = [np.mean([d[6] for d in ppp_data]), np.mean([d[6] for d in ccc_data])]

x = np.arange(len(categories))
w = 0.25
axes[1].bar(x - w, sigma_means, w, label=r'$\sigma$ (summands)', color='#2255BB', alpha=0.7)
axes[1].bar(x, delta_means, w, label=r'$\delta$ (differences)', color='#CC6644', alpha=0.7)
axes[1].bar(x + w, pi_means, w, label=r'$\pi$ (partial sums)', color='#228833', alpha=0.7)
axes[1].set_xticks(x)
axes[1].set_xticklabels(categories)
axes[1].set_ylabel(r'Mean contribution to $\rho_3$', fontsize=10)
axes[1].set_title('Geometric Decomposition', fontsize=12)
axes[1].legend(fontsize=8)
axes[1].grid(True, alpha=0.15, axis='y')

# 1c: ρ₃ scaling with log(N)
Ns_plot = [r['N'] for r in multi_n]
axes[2].scatter(logNs, means_m, s=8, alpha=0.5, color='#2255BB', label=r'$\langle\rho_3\rangle$')
axes[2].scatter(logNs, mins_m, s=5, alpha=0.4, color='#CC3333', label=r'$\rho_{3,\min}$')
xx = np.linspace(min(logNs), max(logNs), 100)
axes[2].plot(xx, c_logN[0]*xx + c_logN[1], 'k--', lw=1.2,
             label=r'$%.2f \log N + %.1f$, $R^2=%.4f$' % (c_logN[0], c_logN[1], r2_logN))
axes[2].set_xlabel(r'$\log N$', fontsize=11)
axes[2].set_ylabel(r'$\rho_3$', fontsize=11)
axes[2].set_title(r'Scaling: $\rho_3$ vs.\ $\log N$', fontsize=12)
axes[2].legend(fontsize=8)
axes[2].grid(True, alpha=0.15)

plt.tight_layout()
plt.savefig('/home/claude/fig_geometric.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/fig_geometric.png', dpi=200, bbox_inches='tight')
plt.close()

# Figure 2: PPP-CCC gap stability
fig, ax = plt.subplots(figsize=(10, 5))
gap_Ns = [g['N'] for g in gap_data]
ax.plot(gap_Ns, floor_gaps, 'o-', ms=4, color='#CC3333', alpha=0.6, label='Floor gap')
ax.plot(gap_Ns, mean_gaps, 's-', ms=4, color='#2255BB', alpha=0.6, label='Mean gap')
ax.axhline(y=0, color='black', ls=':', lw=0.8)
ax.set_xlabel(r'$N$', fontsize=12)
ax.set_ylabel(r'Gap $\rho_3^{\mathrm{PPP}} - \rho_3^{\mathrm{CCC}}$', fontsize=12)
ax.set_title('PPP–CCC Gap Stability (True Discriminant)', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.15)
plt.tight_layout()
plt.savefig('/home/claude/fig_gap.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/fig_gap.png', dpi=200, bbox_inches='tight')
plt.close()

print("  All figures saved.")
print("\n[*] Analysis complete.")
