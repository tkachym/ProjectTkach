# A3 – GLV-декомпозиція скаляра α

# Цей модуль реалізує Алгоритм 3.74 з книги [1]:
# 'Guide to Elliptic Curve Cryptography' (Hankerson, Menezes, Vanstone)

# Мета:
# Знайти представлення α = α₁ + α₂λ (mod q), де:
# - α₁, α₂ мають приблизно вдвічі менший бітовий розмір, ніж α
# - λ — власне значення ендоморфізму φ, знайдене у A2
# - q — порядок групи точок на кривій E(𝔽ₚ)

# Формула:
#     μ = α ⋅ λ / q
#     α₂ = round(μ)
#     α₁ = α − α₂ ⋅ λ (mod q)

# Використовується точна арифметика Decimal для контролю точності.


from decimal import Decimal, getcontext

# Висока точність для роботи з дробами
getcontext().prec = 200

# GLV-декомпозиція: α = α₁ + α₂λ (mod q), де α, q, λ ∈ ℤ
def glv_decompose(alpha, q, lam):
    mu = Decimal(alpha) * Decimal(lam) / Decimal(q)
    alpha2 = int(mu.to_integral_value(rounding='ROUND_HALF_EVEN'))
    alpha1 = (alpha - alpha2 * lam) % q

    if alpha1 > q // 2:
        alpha1 -= q
    if alpha2 > q // 2:
        alpha2 -= q

    return alpha1, alpha2

# Приклад:
q = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
lambda_val = 21888242871839275217838484774961031246154997185409878258781734729429964517155
alpha = 1234567890123456789012345678901234567890

a1, a2 = glv_decompose(alpha, q, lambda_val)
print("α₁ =", a1)
print("α₂ =", a2)
