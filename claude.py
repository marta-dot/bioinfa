import csv
import os
import re
from collections import defaultdict, deque
import heapq
import random
import sys
from time import time


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


def build_de_bruijn_graph(S, l):
    """Buduje graf de Bruijna z oligo w S."""
    # Budujemy graf jako słownik: prefix -> lista (suffix, cały oligo) par
    graph = defaultdict(list)
    # Potrzebujemy również zbioru wszystkich wierzchołków dla znalezienia potencjalnych startów
    all_nodes = set()

    for oligo in S:
        prefix = oligo[:l - 1]
        suffix = oligo[1:]
        graph[prefix].append((suffix, oligo))
        all_nodes.add(prefix)
        all_nodes.add(suffix)

    return graph, all_nodes


def calculate_node_metrics(graph, all_nodes):
    """Oblicza metryki dla każdego węzła (stopnie wejściowe i wyjściowe)."""
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)

    for prefix, edges in graph.items():
        out_degree[prefix] = len(edges)
        for suffix, _ in edges:
            in_degree[suffix] += 1

    # Znajdź potencjalne wierzchołki startowe i końcowe
    potential_starts = []
    potential_ends = []

    for node in all_nodes:
        # Potencjalne starty: węzły z przewagą wyjść nad wejściami
        if out_degree[node] > in_degree[node]:
            potential_starts.append((node, out_degree[node] - in_degree[node]))

        # Potencjalne końce: węzły z przewagą wejść nad wyjściami
        if in_degree[node] > out_degree[node]:
            potential_ends.append((node, in_degree[node] - out_degree[node]))

    # Sortuj wg różnicy stopni (malejąco)
    potential_starts.sort(key=lambda x: x[1], reverse=True)
    potential_ends.sort(key=lambda x: x[1], reverse=True)

    return in_degree, out_degree, potential_starts, potential_ends


def find_probable_path(graph, start_node, used_oligos, n, l):
    """Znajduje prawdopodobną ścieżkę w grafie używając heurystyki pokrycia."""
    path = []
    current = start_node
    current_seq = start_node

    while len(current_seq) < n and current in graph:
        # Sortuj kandydatów wg heurystyki pokrycia
        candidates = []
        for next_node, oligo in graph[current]:
            if oligo in used_oligos:
                continue

            # Heurystyka - liczba nowych oligo, które możemy potencjalnie pokryć
            # plus stopień wyjściowy następnego węzła
            coverage_score = sum(1 for o in graph.get(next_node, [])
                                 if o[1] not in used_oligos)

            # Preferujemy oligo, które prowadzą do większej liczby niewykorzystanych oligo
            candidates.append((coverage_score, next_node, oligo))

        if not candidates:
            break

        # Wybierz najlepszego kandydata wg heurystyki
        candidates.sort(reverse=True)  # Sortuj malejąco wg coverage_score
        _, next_node, oligo = candidates[0]

        path.append(oligo)
        used_oligos.add(oligo)
        current = next_node
        current_seq += next_node[-1]  # Dodaj ostatni nukleotyd

    return path, current_seq


def solve_dna_sequencing_improved(S, l, n):
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
    print(f"Długość sekwencji: {len(result)}")

    extracted_lmers = {result[i:i + l] for i in range(len(result) - l + 1)}
    S_set = set(S)

    used = extracted_lmers & S_set  # Przecięcie zbiorów
    used_count = len(used)
    total_S = len(S_set)

    print(f"Liczba unikalnych l-merów z S w sekwencji: {used_count}/{total_S}")
    print(f"Procent wykorzystania S: {100 * used_count / total_S:.2f}%")

    # Sprawdź ciągłość sekwencji
    for i in range(len(result) - l):
        current = result[i:i + l]
        next_lmer = result[i + 1:i + l + 1]
        if current[1:] != next_lmer[:-1]:
            print(f"Błąd ciągłości: {current} → {next_lmer}")
            return False

    return used_count


# -----------------------------------------------------------------------------------
# Przykład użycia
if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        file_name = sys.argv[1]
    else:
        # Domyślny plik
        file_name = "9.200-80.txt"

    with open(file_name, 'r') as file:
        oligos = [line.strip() for line in file if line.strip()]

    S = oligos

    # Parse n and l from filename (e.g., "200+40.txt" -> n=200, l=10)
    match = re.match(r".*?(\d+)\.(\d+)([+-])(\d+)", file_name)
    if match:
        n = int(match.group(2))
        l = 10  # Zakładamy l=10 jak w oryginalnym kodzie
        print("n: ",n)
    else:
        print("Nie można odczytać parametrów n i l z nazwy pliku. Użycie domyślnych wartości.")
        n = 209
        l = 10

    # Testowanie oryginalnego algorytmu
    # print("\n=== Oryginalny algorytm ===")

    # Testowanie ulepszonego algorytmu
    print("\n=== Ulepszony algorytm ===")
    start_time = time()
    improved_result = solve_dna_sequencing_improved(S, l, n)
    improved_time = time() - start_time
    print(f"Czas wykonania: {improved_time:.4f} s")
    improved_coverage = verify_solution(improved_result, S, l, n)

    print(f"Pokrycie - ulepszony: {improved_coverage}/{len(S)} ({100 * improved_coverage / len(S):.2f}%)")


    def process_all_files_in_directory(directory="."):
        results = []
        for file_name in os.listdir(directory):
            if not file_name.endswith(".txt"):
                continue

            file_path = os.path.join(directory, file_name)

            with open(file_path, 'r') as file:
                oligos = [line.strip() for line in file if line.strip()]
            S = oligos

            # Domyślne wartości
            n, l = None, 10

            # Spróbuj wyciągnąć n z nazwy pliku
            match = re.match(r".*?(\d+)\.(\d+)([+-])(\d+)", file_name)
            if match:
                n = int(match.group(2)) + 9
            else:
                print(f"Ostrzeżenie: Nie można odczytać n z {file_name}, używam n=200")
                n = 200

            print(f"Przetwarzanie pliku: {file_name} (n={n}, l={l})")

            # Uruchomienie ulepszonego algorytmu
            start_time = time()
            result = solve_dna_sequencing_improved(S, l, n)
            elapsed_time = time() - start_time
            coverage = verify_solution(result, S, l, n)

            results.append({
                "file_name": file_name,
                "n": n,
                "l": l,
                "time": round(elapsed_time, 4),
                "coverage_abs": coverage,
                "coverage_percent": round(100 * coverage / len(S), 2)
            })

        # Zapisz wyniki do CSV
        output_csv = os.path.join(directory, "results.csv")
        with open(output_csv, "w", newline="") as csvfile:
            fieldnames = ["file_name", "n", "l", "time", "coverage_abs", "coverage_percent"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        print(f"\nZapisano wyniki do {output_csv}")


    def process_all_files(input_dir="test_zbiory", output_csv="wyniki.csv", l=10):
        rows = []

        for file_name in sorted(os.listdir(input_dir)):
            if not file_name.endswith(".txt"):
                continue

            match = re.match(r".*?(\d+)\.(\d+)(?:\+(\d+))?(?:-(\d+))?", file_name)
            if not match:
                print(f"❌ Pomijam niepasujący plik: {file_name}")
                continue

            base = int(match.group(2))
            num_pos = int(match.group(3)) if match.group(3) else 0
            num_neg = int(match.group(4)) if match.group(4) else 0
            n = base + l - 1

            with open(os.path.join(input_dir, file_name)) as f:
                S = [line.strip() for line in f if line.strip()]

            print(f"▶ Przetwarzam {file_name} (n={n}, l={l}, +{num_pos}, -{num_neg})")

            start = time()
            result = solve_dna_sequencing_improved(S, l, n)
            elapsed = time() - start

            coverage = verify_solution(result, S, l, n)

            rows.append({
                "plik": file_name,
                "n": n,
                "l": l,
                "błędy_pozytywne": num_pos,
                "błędy_negatywne": num_neg,
                "czas_s": round(elapsed, 4),
                "pokrycie": coverage,
                "rozmiar_S": len(S),
                "pokrycie_%": round(100 * coverage / len(S), 2)
            })

        # Zapis do CSV
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

        print(f"\n📄 Wyniki zapisane do: {output_csv}")


    if __name__ == "__main__":
        process_all_files_in_directory(".")
        process_all_files()