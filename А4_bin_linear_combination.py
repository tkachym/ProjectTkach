# A4: Ефективне обчислення [α₁]P₁+[α₂]P₂ (реалізація Алгоритму 3.48 - Shamir’s Trick)

## Мета: реалізувати обчислення [α₁]P₁+[α₂]P₂
# у контексті GLV-множення, де:
#    - P₁ = P
#    - P₂ = φ(P)
#     - α₁, α₂ — отримані з GLV-декомпозиції (A3)

## Метод
# Застосовується Алгоритм 3.48 з книги [1] (Hankerson, Menezes, Vanstone), також відомий як метод Шаміра (Shamir’s Trick).

### Основна ідея:
# - Побудувати таблицю додавання для 4 можливих комбінацій бітів:
#     (0,0) → 𝒪  
#     (0,1) → P₂  
#     (1,0) → P₁  
#     (1,1) → P₁ + P₂
# - Далі, побітно обробляються двійкові подання α₁ та α₂.

## Результат
# Модуль `shamir_trick.py` обчислює [α₁]P+[α₂]φ(P) ефективніше, ніж окреме множення з наступним додаванням.


p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
b = 3

def point_add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and y1 != y2:
        return None
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

def shamir_trick(k1, P1, k2, P2):
    if k1 < 0:
        k1 = -k1
        P1 = negate_point(P1)
    if k2 < 0:
        k2 = -k2
        P2 = negate_point(P2)

    # Побудова таблиці
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
        if T:
            R = point_add(R, T)
    return R

# Приклад
P = (1, 2)
omega = 21888242871839275220042445260109153167277707414472061641714758635765020556616
phiP = (omega * P[0] % p, P[1])

alpha1 = -8373803423573080002309799091858454828798290858057793217301015255930102072719
alpha2 = 1234567890123456788763724639560935449249

R = shamir_trick(alpha1, P, alpha2, phiP)
print("Результат [α₁]P + [α₂]φ(P):", R)
