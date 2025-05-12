# Bachelor's Thesis – Elliptic Curves focused console application

Welcome! 👋  
This repository contains the codebase for my Bachelor's thesis project — a C++ console application designed to solve mathematical problems on elliptic curves, particularly focused on the **Elliptic Curve Discrete Logarithm Problem (ECDLP)** as well as **Shamir's trick for calculating linear combination of two points**..

## 🔍 What It Does

The application is built to test and demonstrate solving of ECDLP instances, as well as implement Shamir's Trick for multiple instances with ranging bit sizes.

When using the app, the user can choose which algorithm to apply for solving ECDLP:

- **Baby-Step Giant-Step (BSGS)**
- **Pollard’s Rho Algorithm**
- **Naive search**

Or user can chose to start analyzing **Shamir's tricks** execution

This provides flexibility for testing and comparing the performance and behavior of these algorithms.

## ⚙️ Technologies Used

- **C++**
- **NTL (Number Theory Library)** for number theory and elliptic curve operations  
  ➤ [Download NTL](https://libntl.org/download.html)  
  ➤ [Unix installation guide](https://libntl.org/doc/tour-unix.html)
- **Python** for generating new input data

## 📂 Input Data

The app uses structured `.json` files for its input. These files define the elliptic curve parameters, base points, and instances of the problems being solved.

The correct `.json` file is automatically selected based on the functionality the user triggers in the app.

## 🧪 Generating New Test Cases

Included in the repo is a Python script: `generate_curves.py`.

- Users can run the script using `python3 generate_curves.py` command and it will guide them through the process of creating new datasets in form of `.json` files for the application to test. 
- If generated, the application will automatically prioritize these new files over the original ones.

---
