# Truss

Input truss is:

![Alt text](truss.jpeg?raw=true "Truss")

# Features

Geometric nonlinearity.

1. Deformations are added to the node positions.
   1. `pos_i+1 = pos_i + dis_i`
1. Convergence criterion.
   1. A threshold for horizontal translation of the top node.
   1. `( dis_i+1 - dis_i ) / dis_i * 100 < 0.01%`

# Installation

Octave can be installed on Linux by:

https://wiki.octave.org/Octave_for_openSUSE

# Run

Program can be run by:

```bash
octave truss.m
```
