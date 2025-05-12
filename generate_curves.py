import json
import time
from math import ceil, log
from pathlib import Path
from sage.all import *

def get_mode():
    print("Select generation mode:")
    print("1 - Shamir curves")
    print("2 - ECDLP curves")
    while True:
        choice = input("Enter choice (1 or 2): ").strip()
        if choice in ["1", "2"]:
            return choice
        print("Invalid input. Please enter 1 or 2.")

def get_ecdlp_solver_type():
    print("Select what ECDLP analysis should generated curves be for:")
    print("1 - Naive search")
    print("2 - Pollard's Rho or BSGS")
    while True:
        choice = input("Enter choice (1 or 2): ").strip()
        if choice in ["1", "2"]:
            return choice
        print("Invalid input. Please enter 1 or 2.")

def get_bit_size_range():
    while True:
        try:
            user_input = input("Enter bit size range with optional step (e.g., 220-250,5): ").strip()
            if ',' in user_input:
                range_part, step_part = user_input.split(',')
                step = int(step_part)
            else:
                range_part = user_input
                step = 1
            start, end = map(int, range_part.split('-'))
            if start > end or start < 1 or step < 1:
                raise ValueError
            return list(range(start, end + 1, step))
        except Exception:
            print("Invalid input. Use format like 220-250,5 or 220-250.")

def get_iterations():
    while True:
        try:
            iterations = int(input("Enter number of instances per bit size: "))
            if iterations <= 0:
                raise ValueError
            return iterations
        except Exception:
            print("Invalid input. Enter a positive integer.")

def generate_shamir_curves():
    bit_sizes = get_bit_size_range()
    iterations_per_size = get_iterations()
    output_path = Path("shamircurves_extra.json")
    results = []

    for p_size in bit_sizes:
        print(f"\n[---] Generating data for {p_size} bits ---")
        for i in range(iterations_per_size):
            print(f"[*] Iteration {i+1}/{iterations_per_size} for {p_size} bits")
            try:
                start_time = time.time()
                p = random_prime(2**p_size - 1, False, 2**(p_size - 1))
                a = randint(1, p - 1)
                b = randint(1, p - 1)
                E = EllipticCurve(GF(p), [a, b])
                G = E.random_element()
                G2 = E.random_element()
                results.append({
                    "bit_size": str(p_size),
                    "prime": str(p),
                    "curve": {
                        "a": str(a),
                        "b": str(b)
                    },
                    "curve_order": str(E.order()),
                    "P": {
                        "x": str(G.xy()[0]),
                        "y": str(G.xy()[1])
                    },
                    "Q": {
                        "x": str(G2.xy()[0]),
                        "y": str(G2.xy()[1])
                    },
                })
                print(f"[OK] Done in {round(time.time() - start_time, 2)}s")
            except Exception as e:
                print(f"[WRONG] Error: {e}")
                continue

    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)

    print(f"\n[OK] Script completed. Total curves generated: {len(results)}")

def generate_ecdlp_curves():
    type = get_ecdlp_solver_type()
    bit_sizes = get_bit_size_range()
    iterations_per_size = get_iterations()
    
    if type == "1":
        output_path = Path("ecdlpcurvesnaive_extra.json")
    else:
        output_path = Path("ecdlpcurves_extra.json")
        
    results = []

    for p_size in bit_sizes:
        print(f"[...] Trying for {p_size} bits...")
        for _ in range(iterations_per_size):
            while True:
                try:
                    p = random_prime(2**p_size - 1, False, 2**(p_size - 1))
                    a = randint(1, p - 1)
                    b = randint(1, p - 1)
                    E = EllipticCurve(GF(p), [a, b])
                    E_order = E.order()
                    for _ in range(100):
                        G = E.random_element()
                        if G.order() == E_order:
                            break
                    else:
                        continue
                    G_order = G.order()
                    if is_prime(G_order):
                        results.append({
                            "bit_size": p_size,
                            "prime": str(p),
                            "curve": {
                                "a": str(a),
                                "b": str(b)
                            },
                            "curve_order": str(E_order),
                            "generator": {
                                "x": str(G.xy()[0]),
                                "y": str(G.xy()[1])
                            },
                            "generator_order": str(G_order),
                            "log2_order": str(ceil(log(G_order) / log(2)))
                        })
                        print(f"[OK] {p_size} bits.")
                        break
                except Exception as e:
                    print(f"[RETRY] {p_size} bits: {e}")
                    continue

    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)

    print(f"\n[OK] Script completed. Total curves generated: {len(results)}")

def main():
    mode = get_mode()
    if mode == "1":
        generate_shamir_curves()
    else:
        generate_ecdlp_curves()

main()
