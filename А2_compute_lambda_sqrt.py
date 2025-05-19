# A2 – Визначення λ для заданого ω

# Мета:
# Знайти λ ∈ 𝔽_q таке, що:     φ(P) = (ωx, y) = [λ]P
# Це λ буде власним значенням ендоморфізму φ, яке узгоджено з ω, знайденим у A1.

# Метод
# - Оскільки φ³ = id, то λ повинно задовольняти:
#       λ² + λ + 1 ≡ 0 mod q
# - Розв'язки цього рівняння:
#       λ = (–1 ± √–3)/2 mod q
# - Після перевірки обидвох λ, залишили те, для якого:
#       φ(G) = [λ]G для G = (1, 2)


from sympy.ntheory import sqrt_mod

p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
q = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
omega = 21888242871839275220042445260109153167277707414472061641714758635765020556616
G = (1, 2)
phiG = (omega, 2)

# Обчислення коренів з –3 mod q
neg3 = (-3) % q
roots = sqrt_mod(neg3, q, all_roots=True)

# Кандидати на λ = (–1 ± √–3)/2 mod q
lambdas = [(( -1 + r ) * pow(2, -1, q)) % q for r in roots]

# Еліптична арифметика
def point_add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and y1 != y2: return None
    if P == Q:
        m = (3 * x1 * x1) * pow(2 * y1, -1, p) % p
    else:
        m = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    return (x3, y3)

def scalar_mult(k, P):
    R = None
    while k:
        if k & 1:
            R = point_add(R, P)
        P = point_add(P, P)
        k >>= 1
    return R

# Перевірка: φ(G) == [λ]G ?
for idx in range(len(lambdas)):
    lam = lambdas[idx]
    result = scalar_mult(lam, G)
    print("λ =", lam)
    print("φ(G) =", phiG)
    print("[λ]G =", result)
    print("φ(G) == [λ]G:", result == phiG)
    print("---")
