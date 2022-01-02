# Truss

Input truss is:

![Alt text](truss.jpeg?raw=true "Truss")

# TODO

Geometric nonlinearity.

1. Add deformations to the node positions.
1. Convergence criterion.
   1. If deformation is less than a threshold.
   1. `pos_i+1` = `pos_i` + `dis_i`
   1. Horizontal translation of top node.
   1. `( dis_i+1 - dis_i ) / dis_i * 100 < 0.01%`
