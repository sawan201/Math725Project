# csf_trees.sage

SageMath code accompanying the Math 725 independent project report,
*On Distinguishing Trees by Their Chromatic Symmetric Functions: A
Report on Martin–Morin–Wagner*.

## Requirements

SageMath 9.0 or later. No dependencies beyond the standard SageMath
distribution.

## Running

From the project directory:

    sage csf_trees.sage

This runs three demos in sequence:

1. The full power-sum expansion of the chromatic symmetric function
   of the spider T_{(2,1,1)}, reproducing Example 6.1 of the report
   and verifying each part of Proposition 4.1.

2. Reconstruction of the caterpillar with leaf-number sequence
   (1, 2, 3) from its CSF coefficient data, reproducing Example 6.2.

3. A brute-force CSF-collision search over all positive caterpillars
   on up to 11 vertices, providing the empirical evidence for
   Conjecture 7.1.

The collision search up to n = 11 takes roughly a minute on a laptop.
Pushing it to n = 12 is feasible but slower; beyond that the
2^|E| edge-subset enumeration starts to bite.

## Main functions

`chromatic_symmetric_function(G)` returns X_G as an element of the
power-sum basis of `Sym`, computed by enumerating all edge subsets of
G and applying Stanley's expansion.

`csf_coefficient_dict(G)` returns the same data as a
`{Partition: int}` dictionary with zero coefficients dropped. Handier
when you want to look up specific c_λ.

`caterpillar(leaf_numbers)` builds the caterpillar with the given
leaf-number sequence as a Sage `Graph`.

`reconstruct_caterpillar(coeffs)` takes the coefficient dictionary of
a caterpillar with strictly positive, pairwise distinct leaf numbers
and returns the leaf-number sequence (canonicalised so the smaller
endpoint comes first). This is Theorem 5.1 of the report implemented
step-for-step. Raises `ValueError` when the input is inconsistent
with those hypotheses.

`test_caterpillar_conjecture(n_max)` runs the empirical
CSF-collision search over all positive caterpillars on at most
`n_max` vertices.

## Example session

    sage: load("csf_trees.sage")
    sage: T = caterpillar((1, 2, 3))
    sage: coeffs = csf_coefficient_dict(T)
    sage: reconstruct_caterpillar(coeffs)
    (1, 2, 3)

## Notes

The CSF is stored in the power-sum basis throughout because Stanley's
expansion gives the p-coefficients directly. Conversion to another
basis is one call away if needed:

    XT = chromatic_symmetric_function(T)
    XT_e = SymmetricFunctions(QQ).elementary()(XT)

Edge-subset enumeration is 2^|E|, so this code is meant for small
examples and for verifying the claims in the report, not for trees on
more than around 25 edges. The brute-force conjecture test is the
main bottleneck; everything else runs comfortably for the sizes
considered in the report.
