import random
import os

def generate_lmers_with_errors(n, l, num_negative, num_positive, seed=None):
    """Tworzy listę l-merów z błędami pozytywnymi i negatywnymi."""
    if seed is not None:
        random.seed(seed)

    dna = ''.join(random.choice('ACGT') for _ in range(n))
    original_lmers = {dna[i:i + l] for i in range(n - l + 1)}

    # Usuń część poprawnych l-merów
    lmers = list(original_lmers)
    random.shuffle(lmers)
    lmers = lmers[:-num_negative] if num_negative < len(lmers) else []

    # Dodaj losowe błędne l-mery
    while len(lmers) < len(original_lmers) - num_negative + num_positive:
        candidate = ''.join(random.choice('ACGT') for _ in range(l))
        if candidate not in original_lmers and candidate not in lmers:
            lmers.append(candidate)

    random.shuffle(lmers)
    return lmers

def generate_test_instances(output_dir="test_instances", instances=5, n=200, l=10, seed_base=42):
    """Tworzy wiele instancji testowych do katalogu `output_dir`."""
    os.makedirs(output_dir, exist_ok=True)

    for i in range(instances):
        n = (i+2) * 100 + 9
        num_negative = random.randint(1, 10) * 10
        num_positive = random.randint(1, 10) * 10
        seed = seed_base + i

        lmers = generate_lmers_with_errors(n, l, num_negative, num_positive, seed=seed)

        # liczba oryginalnych l-merów: n - l + 1
        base = n - l + 1
        file_name = f"9.{base}+{num_positive}-{num_negative}.txt"
        path = os.path.join(output_dir, file_name)

        with open(path, "w") as f:
            for oligo in lmers:
                f.write(oligo + "\n")

        print(f"✔ Zapisano: {path}")


if __name__ == "__main__":
    generate_test_instances(
        output_dir="test_zbiory",
        instances=5,
        n=200,
        l=10,
        seed_base=100
    )
