# Gosper's Summation Algorithm
An implementation of Gosper's algorithm for computing indefinite sums of hypergeometric functions. 

## Overview
This project contains a Maple implementation and a Python implementation of Gosper's algorithm.

## Installation

### Maple

**Note:**
The Maple implementation is currently broken.

The Maple implementation requries Maple 2022. It might work on other Maple versions but this is untested.

### Python
The Python implementation requires Python 3.9.6, Sympy 1.11.1 and any compatible version of pytest. It might work on other versions but this is untested. To run the test suite, execute
```
pytest .
```
in the root directory of the project.

## Execution

### Maple

**Note:**
The Maple implementation is currently broken.

### Python

The `sympy` folder contains a sample `main.py` file for playing with Gosper's algorithm. 

## References

Gosper, R. William Jr. "Decision Procedure for Indefinite Hypergeometric Summation." (1977).

Graham, R. L., Knuth, D. E., Patashnik, O. "Concrete Mathematics: A Foundation for Computer Science." (1989).

Hayden, Michael B. and Lamanga, Edmund A. "Gosper's Summation Algorithm and Some Applications." (1988).

Koepf, Wolfram and Roberto Pirastu. “Summation in Maple.” (1996).
