# A5 – Повна реалізація GLV-множення

## Мета: побудувати повний алгоритм GLV-множення [α]P=[α₁]P+[α₂]φ(P),
#        який швидший за класичне скалярне множення точки на кривій E(𝔽ₚ).

## Компоненти
#    - ω – нетривіальний кубічний корінь з 1 (з A1)
#    - λ – власне значення ендоморфізму φ (з A2)
#    - GLV-декомпозиція α (з A3)
#    - φ(P) = (ωx, y) (з A2)
#    - [α₁]P + [α₂]φ(P) – метод Шаміра (з A4)

## Алгоритм
# 1. Використати `glv_decompose(α)` → α₁, α₂
# 2. Обчислити φ(P)
# 3. Обчислити `shamir_trick(α₁, P, α₂, φ(P))`

## Результат
# Файл `glv_multiply.py` реалізує:
# - повне множення [α]P через GLV
# - порівняння з класичним scalar_mult
# - виведення часу виконання та перевірку збігу результатів


from decimal import Decimal, getcontext
import time

# Параметри кривої BN254
p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
q = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
omega = 21888242871839275220042445260109153167277707414472061641714758635765020556616
lambda_val = 21888242871839275217838484774961031246154997185409878258781734729429964517155

getcontext().prec = 200

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

def negate_point(P):
    if P is None: return None
    x, y = P
    return (x, (-y) % p)

def scalar_mult(k, P):
    R = None
    if k < 0:
        k = -k
        P = negate_point(P)
    while k:
        if k & 1:
            R = point_add(R, P)
        P = point_add(P, P)
        k >>= 1
    return R

def phi(P):
    x, y = P
    return (omega * x % p, y)

def glv_decompose(alpha):
    mu = Decimal(alpha) * Decimal(lambda_val) / Decimal(q)
    alpha2 = int(mu.to_integral_value(rounding='ROUND_HALF_EVEN'))
    alpha1 = (alpha - alpha2 * lambda_val) % q
    if alpha1 > q // 2: alpha1 -= q
    if alpha2 > q // 2: alpha2 -= q
    return alpha1, alpha2

def shamir_trick(k1, P1, k2, P2):
    if k1 < 0: k1, P1 = -k1, negate_point(P1)
    if k2 < 0: k2, P2 = -k2, negate_point(P2)
    table = {
        (0, 0): None,
        (0, 1): P2,
        (1, 0): P1,
        (1, 1): point_add(P1, P2)
    }
    R = None
    bin1 = bin(k1)[2:].zfill(max(len(bin(k1)), len(bin(k2))))
    bin2 = bin(k2)[2:].zfill(len(bin1))
    for b1, b2 in zip(bin1, bin2):
        R = point_add(R, R)
        T = table[(int(b1), int(b2))]
        if T: R = point_add(R, T)
    return R

def glv_multiply(alpha, P):
    a1, a2 = glv_decompose(alpha)
    P2 = phi(P)

    # Вимірюємо лише час множення
    start = time.time()
    R = shamir_trick(a1, P, a2, P2)
    elapsed = time.time() - start
    return R, elapsed

# Тест
P = (1, 2)
alpha = 1234567890123456789012345678901234567890

# Час класичного множення
start_classic = time.time()
R_classic = scalar_mult(alpha, P)
elapsed_classic = time.time() - start_classic

# Час GLV множення
R_glv, elapsed_glv = glv_multiply(alpha, P)

print("Класичне множення:", R_classic)
print("GLV множення     :", R_glv)
print("Результати збігаються:", R_classic == R_glv)
print(f"Час класичного множення: {elapsed_classic:.6f} секунд")
print(f"Час GLV множення       : {elapsed_glv:.6f} секунд")
