# generated_modes
generated_modes is a repository in which the twining modes in different crystal structures can be predicted.

# Code Description

## Functions

### equivalent.m

The function will generate `3x3` integer matrices with determinant $\pm 1$, where the largest absolute value is fixed. 

 + Input parameters: **n** -  The maximum absolute value for the integer in the matrices generated.
 + Output -
   * **mo**: A cell with all the GL(3,Z) matrices.
   * **count**: Number of matrices.

## Source Codes

### generated_modes11072017.m

+ Creates all the twinning modes for GL(3,Z) matrices with maximum integer 4. Variables `k_1` and `eta_1` describe the vectors needed to define the twinning mode.
+ Calculates the shear associated with each twinning mode. This is in the variable `shear`.

