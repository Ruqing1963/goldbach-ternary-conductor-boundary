"""
ternary_corrected.py — GSp(6) conductor based on the TRUE discriminant

The genus-3 hyperelliptic curve C: y² = x(x²-p₁²)(x²-p₂²)(x²-p₃²)
has discriminant:
  Δ ∝ ∏ pᵢ⁶ · ∏_{i<j} (pᵢ-pⱼ)⁴ · ∏_k (N-pₖ)⁴

where (N - pₖ) = pᵢ + pⱼ (the complementary partial sum).

The conductor proxy is therefore:
  N_proxy = rad_odd(∏ pᵢ⁶ · ∏(pᵢ-pⱼ)⁴ · ∏(N-pₖ)⁴)

Since rad strips exponents:
  rad_odd(Δ) = rad_odd(p₁·p₂·p₃ · (p₁-p₂)(p₂-p₃)(p₃-p₁) · (N-p₁)(N-p₂)(N-p₃))

Key difference from Paper #11: NO independent factor of N.
N only enters via (N-pₖ) = partial sums.
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sieve_and_radicals(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    rads = np.ones(limit + 1, dtype=np.int64)
    for i in range(2, limit + 1):
        if is_prime[i]:
            for j in range(i, limit + 1, i):
                if j != i: is_prime[j] = 0
                rads[j] *= i
    return is_prime, rads

def odd_rad(x, rads):
    if x == 0: return 0
    x = abs(x)
    r = int(rads[x])
    while r % 2 == 0: r //= 2
    return r if r > 0 else 1

LIMIT = 20000
is_prime, rads = sieve_and_radicals(LIMIT)

# ═══════════════════════════════════════════════════════════════════════════════
# CORRECTED CONDUCTOR PROXY
# ═══════════════════════════════════════════════════════════════════════════════
def true_rho3(N, p1, p2, p3, rads):
    """
    Compute ρ₃ from the TRUE discriminant of the genus-3 Frey curve.
    
    Δ ∝ (p1·p2·p3)⁶ · ((p1-p2)(p2-p3)(p3-p1))⁴ · ((N-p1)(N-p2)(N-p3))⁴
    
    rad_odd(Δ) = rad_odd(p1·p2·p3 · (p1-p2)·(p2-p3)·(p3-p1) · (N-p1)·(N-p2)·(N-p3))
    
    Note: (N-pₖ) = sum of other two summands.
    """
    logN = math.log(N)
    
    # Summand radicals
    r1 = odd_rad(p1, rads)
    r2 = odd_rad(p2, rads)
    r3 = odd_rad(p3, rads)
    
    # Difference radicals
    d12 = odd_rad(abs(p1 - p2), rads)
    d23 = odd_rad(abs(p2 - p3), rads)
    d31 = odd_rad(abs(p3 - p1), rads)
    
    # Partial sum radicals (N - pₖ)
    s1 = odd_rad(N - p1, rads)  # = p2 + p3
    s2 = odd_rad(N - p2, rads)  # = p1 + p3
    s3 = odd_rad(N - p3, rads)  # = p1 + p2
    
    # Full radical of discriminant (exponents stripped by rad)
    # Exponents from Δ: pᵢ⁶, diffs⁴, sums⁴
    # rad ignores exponents, so:
    # log(rad_odd(Δ)) = log(rad_odd(p1·p2·p3·d12·d23·d31·s1·s2·s3))
    # But rad of a product ≤ product of rads (with shared factors counted once)
    # For the conductor, we want the actual radical of the full product
    
    # Compute rad_odd of the combined product via prime factorization
    # Efficient: collect all prime factors
    factors = set()
    for val in [r1, r2, r3, d12, d23, d31, s1, s2, s3]:
        if val == 0: continue
        v = val
        for pp in range(3, int(math.isqrt(v)) + 2, 2):
            if v % pp == 0:
                factors.add(pp)
                while v % pp == 0: v //= pp
            if v <= 1: break
            if pp * pp > v:
                if v > 1: factors.add(v)
                break
        if v > 1: factors.add(v)
    
    log_rad = sum(math.log(f) for f in factors) if factors else 0.0
    
    # ρ₃ based on discriminant exponents:
    # log(conductor) ~ 6·log(rad(p's)) + 4·log(rad(diffs)) + 4·log(rad(sums))
    # But since rad strips exponents, we use the actual exponents from Δ
    
    # More precisely: conductor exponent at prime ℓ is determined by ℓ-adic valuation
    # For our proxy, use: log N_proxy = 6·log r_p + 4·log r_d + 4·log r_s
    # where r_p = rad_odd(p1·p2·p3), r_d = rad_odd(diffs), r_s = rad_odd(sums)
    
    log_rp = math.log(r1) + math.log(r2) + math.log(r3) if (r1*r2*r3 > 0) else 0
    log_rd = 0
    for d in [d12, d23, d31]:
        if d > 0: log_rd += math.log(d)
    log_rs = 0
    for s in [s1, s2, s3]:
        if s > 0: log_rs += math.log(s)
    
    log_cond = 6 * log_rp + 4 * log_rd + 4 * log_rs
    rho = log_cond / logN if logN > 0 else 0
    
    return rho

# ═══════════════════════════════════════════════════════════════════════════════
# SCAN: N = 1025 with CORRECTED proxy
# ═══════════════════════════════════════════════════════════════════════════════
N = 1025
print(f"N = {N}")
print(f"Using TRUE discriminant: Δ ∝ ∏pᵢ⁶ · ∏(pᵢ-pⱼ)⁴ · ∏(N-pₖ)⁴")
print(f"NO artificial rad_odd(N) factor.\n")

ppp_rhos, ccc_rhos, mix_rhos = [], [], []
ppp_data, ccc_data = [], []

for p1 in range(3, N // 3 + 1):
    for p2 in range(p1, (N - p1) // 2 + 1):
        p3 = N - p1 - p2
        if p3 < p2: break
        if p3 <= 1: continue
        
        rho = true_rho3(N, p1, p2, p3, rads)
        
        ip1 = bool(is_prime[p1])
        ip2 = bool(is_prime[p2])
        ip3 = bool(is_prime[p3])
        np_ = ip1 + ip2 + ip3
        
        if np_ == 3:
            ppp_rhos.append(rho)
            ppp_data.append((p1, p2, p3, rho))
        elif np_ == 0:
            ccc_rhos.append(rho)
            ccc_data.append((p1, p2, p3, rho))
        else:
            mix_rhos.append(rho)

ppp_rhos = np.array(ppp_rhos)
ccc_rhos = np.array(ccc_rhos)
mix_rhos = np.array(mix_rhos)

print(f"{'Category':<16} {'Count':<8} {'ρ_min':<10} {'⟨ρ⟩':<10} {'ρ_max':<10} {'BW':<10}")
print("-" * 64)
for label, arr in [("PPP", ppp_rhos), ("Mixed", mix_rhos), ("CCC", ccc_rhos)]:
    if len(arr) > 0:
        print(f"{label:<16} {len(arr):<8} {min(arr):<10.4f} {np.mean(arr):<10.4f} {max(arr):<10.4f} {max(arr)-min(arr):<10.4f}")

if len(ppp_rhos) > 0 and len(ccc_rhos) > 0:
    print(f"\nPPP floor: {min(ppp_rhos):.4f}")
    print(f"CCC floor: {min(ccc_rhos):.4f}")
    print(f"Gap (PPP_min - CCC_min): {min(ppp_rhos) - min(ccc_rhos):.4f}")

print("\nPPP ground states (lowest ρ₃):")
sorted_ppp = sorted(ppp_data, key=lambda t: t[3])
for t in sorted_ppp[:5]:
    p1,p2,p3,rho = t
    print(f"  ({p1}, {p2}, {p3})  ρ₃ = {rho:.4f}")
    # Show components
    d12 = abs(p1-p2); d23 = abs(p2-p3); d31 = abs(p3-p1)
    s1 = N-p1; s2 = N-p2; s3 = N-p3
    print(f"    diffs: {d12},{d23},{d31}  sums: {s1},{s2},{s3}")

print("\nCCC ground states (lowest ρ₃):")
sorted_ccc = sorted(ccc_data, key=lambda t: t[3])
for t in sorted_ccc[:5]:
    p1,p2,p3,rho = t
    print(f"  ({p1}, {p2}, {p3})  ρ₃ = {rho:.4f}")
    d12 = abs(p1-p2); d23 = abs(p2-p3); d31 = abs(p3-p1)
    s1 = N-p1; s2 = N-p2; s3 = N-p3
    print(f"    diffs: {d12},{d23},{d31}  sums: {s1},{s2},{s3}")

# ═══════════════════════════════════════════════════════════════════════════════
# BSL TEST with CORRECTED proxy over odd N in [1001, 2001]
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("BSL TEST: Corrected proxy, odd N in [1001, 2001]")
print("=" * 80)

# What is the "natural" static variable now?
# In the TRUE discriminant, N only enters via (N - pₖ).
# There is NO independent N-factor to define a universal ξ.
# Instead, let's test: does ⟨ρ₃⟩ correlate with rad_odd(N)?

bsl = []
for N_scan in range(1001, 2002, 2):
    rN = odd_rad(N_scan, rads)
    logN = math.log(N_scan)
    xi_old = 2 * math.log(rN) / logN if rN > 1 else 0  # the OLD (wrong) variable
    
    ppp_r = []
    for p1 in range(3, N_scan // 3 + 1, 2):
        if not is_prime[p1]: continue
        for p2 in range(p1, (N_scan - p1) // 2 + 1, 2):
            if not is_prime[p2]: continue
            p3 = N_scan - p1 - p2
            if p3 < p2 or p3 <= 1: continue
            if not is_prime[p3]: continue
            rho = true_rho3(N_scan, p1, p2, p3, rads)
            ppp_r.append(rho)
    
    if ppp_r:
        bsl.append({
            'N': N_scan, 'rad_N': rN, 'xi_old': xi_old,
            'n_ppp': len(ppp_r),
            'rho_mean': np.mean(ppp_r),
            'rho_min': min(ppp_r),
            'rho_max': max(ppp_r),
            'bw': max(ppp_r) - min(ppp_r),
        })

xis_old = np.array([r['xi_old'] for r in bsl])
means = np.array([r['rho_mean'] for r in bsl])
mins = np.array([r['rho_min'] for r in bsl])
bws = np.array([r['bw'] for r in bsl])

# Regression against the OLD ξ (which is now just a correlate, not a cause)
c = np.polyfit(xis_old, means, 1)
resid = means - (c[0]*xis_old + c[1])
r2 = 1 - np.var(resid)/np.var(means)

print(f"  Regression ⟨ρ₃⟩ vs ξ_old (= 2log(rad_odd(N))/logN):")
print(f"    slope = {c[0]:.4f}, intercept = {c[1]:.4f}, R² = {r2:.6f}")
print(f"  Bandwidth: {np.mean(bws):.4f} ± {np.std(bws):.4f}")
print(f"  PPP count range: [{min([r['n_ppp'] for r in bsl])}, {max([r['n_ppp'] for r in bsl])}]")

# Also try: regression against log(N) itself
logNs = np.array([math.log(r['N']) for r in bsl])
c2 = np.polyfit(logNs, means, 1)
resid2 = means - (c2[0]*logNs + c2[1])
r2_logN = 1 - np.var(resid2)/np.var(means)
print(f"\n  Regression ⟨ρ₃⟩ vs log(N):")
print(f"    slope = {c2[0]:.4f}, intercept = {c2[1]:.4f}, R² = {r2_logN:.6f}")

# Try: regression against nothing (is ρ₃ just roughly constant for PPP?)
print(f"\n  ⟨ρ₃⟩ statistics: mean = {np.mean(means):.4f}, std = {np.std(means):.4f}")
print(f"  ρ₃_min statistics: mean = {np.mean(mins):.4f}, std = {np.std(mins):.4f}")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE: Corrected distributions
# ═══════════════════════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

bins = np.linspace(min(min(ccc_rhos), min(ppp_rhos)) - 0.5,
                   max(max(ccc_rhos), max(ppp_rhos)) + 0.5, 80)
ax1.hist(ccc_rhos, bins=bins, alpha=0.4, color='#CC6644', density=True, label=f'CCC (n={len(ccc_rhos)})')
ax1.hist(mix_rhos, bins=bins, alpha=0.2, color='#888888', density=True, label=f'Mixed (n={len(mix_rhos)})')
ax1.hist(ppp_rhos, bins=bins, alpha=0.6, color='#2255BB', density=True, label=f'PPP (n={len(ppp_rhos)})')
ax1.axvline(x=min(ppp_rhos), color='#2255BB', ls='--', lw=1.5)
ax1.axvline(x=min(ccc_rhos), color='#CC6644', ls='--', lw=1.5)
ax1.set_xlabel(r"Corrected $\rho_3$ (true discriminant)", fontsize=11)
ax1.set_ylabel('Density', fontsize=11)
ax1.set_title(f'$N = 1025$: Corrected Conductor Distribution', fontsize=13)
ax1.legend(fontsize=8, loc='upper left')
ax1.grid(True, alpha=0.15)

ax2.scatter(xis_old, means, s=6, alpha=0.4, color='#2255BB', label=r'$\langle\rho_3\rangle$')
ax2.scatter(xis_old, mins, s=4, alpha=0.3, color='#CC3333', label=r'$\rho_{3,\min}$')
xx = np.linspace(0, max(xis_old)*1.05, 100)
ax2.plot(xx, c[0]*xx + c[1], 'k--', lw=1.2,
         label=r'$R^2 = %.4f$' % r2)
ax2.set_xlabel(r'$\xi = 2\log(\mathrm{rad}_{\mathrm{odd}}(N))/\log N$ (NOT in discriminant)', fontsize=9)
ax2.set_ylabel(r'Corrected $\rho_3$', fontsize=11)
ax2.set_title('BSL Test: Corrected Proxy', fontsize=13)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.15)

plt.tight_layout()
plt.savefig('/home/claude/ternary_corrected.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/ternary_corrected.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved ternary_corrected.pdf/.png")
