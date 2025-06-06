# Bachelor's Thesis – Elliptic Curves Focused Console Application

Welcome! 👋
This repository contains the codebase for my Bachelor's thesis project — a C++ console application designed to solve mathematical problems on elliptic curves. It particularly focuses on the **Elliptic Curve Discrete Logarithm Problem (ECDLP)** and **Shamir's trick** for calculating linear combinations of two points.

## 🔍 What It Does

The application is built to:
1.  Test and demonstrate solving ECDLP instances using various algorithms.
2.  Implement and analyze the performance of Shamir's Trick for multiple instances with varying bit sizes.

When using the app, you can choose which algorithm to apply for solving the ECDLP:

-   **Baby-Step Giant-Step (BSGS)**
-   **Pollard’s Rho Algorithm**
-   **Naive Search**

Alternatively, you can choose to start analyzing the execution of **Shamir's trick** or it's alternative - **Naive scalar point multiplication**.

## ⚙️ Technologies Used

-   **C++**
-   **NTL (Number Theory Library)** for number theory and elliptic curve operations
    -   ➤ [Download NTL](https://libntl.org/download.html)
    -   ➤ [Unix installation guide](https://libntl.org/doc/tour-unix.html)
-   **Python 3** for generating new input data

## 📂 Input Data

The application uses structured `.json` files for its input. These files define the elliptic curve parameters, base points, and specific instances of the problems to be solved.

-   **Standard Input:** Files like `ecdlpcurves.json` and `shamircurves.json` contain the default problem sets. The correct `.json` file is automatically selected based on the functionality triggered in the app.
-   **Custom Input:** If you generate new test cases using the provided Python script (see below), the newly generated files will have an `_extra` suffix (e.g., `ecdlpcurves_extra.json`). The application will automatically prioritize using these `_extra` files if they exist.

## 📊 Output Data

The application generates `.csv` files to store the results and performance metrics of the executed algorithms.

-   **Generated Results:** When you run the application, it outputs results into `.csv` files (e.g., `ECDLP_BSGS.csv`, `shamir_results.csv`).
-   **Thesis Results:** Files included in the repository that end with `-test.csv` (e.g., `ECDLP_BSGS-test.csv`, `ECDLP_PollardRho-test.csv`) contain the specific results that were tested and are presented as part of the bachelor's thesis documentation. The application will generate corresponding files *without* the `-test` suffix during its execution.

## 🚀 How to Compile and Run

1.  **Prerequisites:** Ensure you have a C++ compiler (like `g++`), the NTL library, and its dependency GMP installed on your system.
2.  **Compile:** Navigate to the repository's root directory in your terminal and run the following command:
    ```bash
    g++ thesis.cpp -o thesis -lntl -lgmp -g -march=native -O2
    ```
3.  **Run:** Execute the compiled application:
    ```bash
    ./thesis
    ```
    The application will then present a menu for you to choose the desired operation.
    

## 🧪 Generating New Test Cases

Included in the repository is a Python script: `generate_curves.py`.

-   You can run this script using the command `python3 generate_curves.py`.
-   The script will guide you through the process of creating new datasets (elliptic curves and problem instances).
-   These new datasets are saved as `.json` files with the `_extra` suffix (e.g., `ecdlpcurves_extra.json`, `shamircurves_extra.json`).
-   As mentioned above, the C++ application will automatically detect and use these `_extra.json` files if they are present, overriding the default input files.

---
