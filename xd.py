import re
from collections import defaultdict, deque
import heapq
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


def solve_dna_sequencing(S, l, n):
    if not S or l <= 0 or n < l:
        return ""

    # Remove duplicate oligos
    S = list(set(mer for mer in S if len(mer) == l))

    # Budowanie grafu De Bruijna i jednoczesne obliczanie stopni wyjściowych dla sufiksów
    de_bruijn = defaultdict(int)
    suffix_out_degree = defaultdict(int)

    for mer in S:
        prefix = mer[:l - 1]
        suffix = mer[1:]
        de_bruijn[prefix] += 1
        suffix_out_degree[suffix] += 1

    # Preprocess prefix maps for overlaps
    prefix_maps = {}
    for o in range(1, l):
        prefix_map = defaultdict(list)
        for mer in S:
            prefix = mer[:o]
            prefix_map[prefix].append(mer)
        prefix_maps[o] = prefix_map

    max_count = 0
    best_sequence = ""

    # Try each oligo as starting point
    for start_mer in S:
        used = {start_mer}
        current_seq = start_mer
        current_len = l
        current_count = 1

        while True:
            found_next = False
            max_o = min(l - 1, len(current_seq))
            # Check overlaps from largest to smallest
            for o in range(max_o, 0, -1):
                suffix = current_seq[-o:]
                candidates = prefix_maps[o].get(suffix, [])
                # Sort candidates by suffix out-degree (descending) to prioritize better extensions
                candidates.sort(key=lambda x: (-suffix_out_degree.get(x[1:], 0), x))
                for cand in candidates:
                    if cand in used:
                        continue
                    new_len = current_len + (l - o)
                    if new_len > n:
                        continue
                    # Use this candidate
                    used.add(cand)
                    current_seq += cand[o:]
                    current_len = new_len
                    current_count += 1
                    found_next = True
                    break  # Move to next step
                if found_next:
                    break
            if not found_next:
                break  # No further extensions

        # Update best sequence
        if (current_count > max_count or
                (current_count == max_count and len(current_seq) < len(best_sequence))):
            max_count = current_count
            best_sequence = current_seq

    return best_sequence[:n]


def verify_solution(result, S, l, n):
    """Weryfikuje poprawność zrekonstruowanej sekwencji DNA."""
    print(f"Długość sekwencji: {len(result)}")

    # Wyodrębnij l-mery z sekwencji wynikowej
    extracted_lmers = [result[i:i + l] for i in range(len(result) - l + 1)]

    # Policz, ile l-merów z S znajduje się w sekwencji
    S_set = set(S)
    used_count = sum(1 for lm in extracted_lmers if lm in S_set)
    total_S = len(S)

    print(f"Liczba l-merów z S w sekwencji: {used_count}/{total_S}")
    print(f"Procent wykorzystania S: {100 * used_count / total_S:.2f}%")

    # W funkcji verify_solution
    for i in range(len(result) - l):
        current = result[i:i + l]
        next_lmer = result[i + 1:i + l + 1]
        if current[1:] != next_lmer[:-1]:
            print(f"Błąd ciągłości: {current} → {next_lmer}")
            return False

    return used_count


# -----------------------------------------------------------------------------------
# Example usage
# Parametry
n = 209  # Długość oryginalnej sekwencji
l = 10  # Długość oligonukleotydów
num_negative = 40  # Liczba brakujących l-merów (błędy negatywne)
num_positive = 30  # Liczba fałszywych l-merów (błędy pozytywne)

# Generowanie danych
original_seq = generate_dna_sequence(n)
original_lmers = [original_seq[i:i + l] for i in range(n - l + 1)]
S = introduce_errors(original_lmers, num_negative, num_positive, l)

# Zapisz do pliku
filename = f"{n}+{num_positive}.txt"
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

# -----------------------------------------------------------------------------------
file_name = "200+60+40.txt"

with open(file_name, 'r') as file:
    oligos = [line.strip() for line in file if line.strip()]

S = oligos

# Parse n and l from filename (e.g., "200-40.txt" -> n=200, l=10)
# match = re.match(r"(\d+)-\d+\.txt", file_name)
# match = re.match(r"(\d+)\+\d+\.txt", file_name)
# match = re.match(r"(\d+)\+\d+\.txt", file_name)
# match = re.match(r"(\d+)[+-]\d+\.txt", file_name)
match = re.match(r"(\d+)([+-]\d+)+\.txt", file_name)
n = int(match.group(1))
l = 10

result = solve_dna_sequencing(S, l, n)
print(f"Reconstructed sequence: {result}")
print(f"Sequence length: {len(result)}")

verify_solution(result, S, l, n)
