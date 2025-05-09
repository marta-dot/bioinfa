import random
import sys


def generate_dna_sequence(n):
    """Generuje losową sekwencję DNA o długości n."""
    return ''.join(random.choice('ACGT') for _ in range(n))


def introduce_errors(original_lmers, num_negative, num_positive, l):
    """Tworzy zbiór S z błędami negatywnymi i pozytywnymi."""
    # Błędy negatywne: usuwamy część oryginalnych l-merów
    S = list(original_lmers)
    random.shuffle(S)
    S = set(S[:len(S) - num_negative])  # Usuwamy `num_negative` l-merów

    # Błędy pozytywne: dodajemy losowe l-mery nieobecne w oryginale
    nucleotides = ['A', 'C', 'G', 'T']
    while len(S) < (len(original_lmers) - num_negative + num_positive):
        new_mer = ''.join(random.choice(nucleotides) for _ in range(l))
        if new_mer not in original_lmers:
            S.add(new_mer)

    return list(S)


# Parametry
n = 200  # Długość oryginalnej sekwencji
l = 10  # Długość oligonukleotydów
num_negative = 40  # Liczba brakujących l-merów (błędy negatywne)
num_positive = 60  # Liczba fałszywych l-merów (błędy pozytywne)

# Generowanie danych
original_seq = generate_dna_sequence(n)
original_lmers = [original_seq[i:i + l] for i in range(n - l + 1)]
S = introduce_errors(original_lmers, num_negative, num_positive, l)

# Zapisz do pliku
filename = f"{n}+{num_positive}+{num_negative}.txt"
with open(filename, 'w') as f:
    for mer in S:
        f.write(mer + '\n')

# Statystyki
print(f"Oryginalna sekwencja: {original_seq[:50]}... (długość {n})")
print(f"Liczba l-merów w idealnym spektrum: {len(original_lmers)}")
print(f"Liczba l-merów w S: {len(S)}")
print(f"  - Brakujące l-mery (błędy negatywne): {num_negative}")
print(f"  - Fałszywe l-mery (błędy pozytywne): {num_positive}")
print(f"Plik wyjściowy: {filename}")