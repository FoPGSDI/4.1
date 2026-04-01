# V. METHODS FOR CALCULATING THE COMPLEX OUTGOING NORMAL MODES

In §IV we have outlined the key role played in real, physical pulsations by the outgoing normal modes. We next turn to the question of, given a general-relativistic stellar model, how might one go about calculating numerically the eigenfrequencies and functions of the outgoing modes.

In our present state of ignorance about analytic properties of the normal modes, the only computational method possible is iterative trial and error. One specifies once and for all the spherical-harmonic index, $l$, and the pulsation amplitude, $B$ (eq. [15a]), at the center of the star. One then specifies successive trial values for the complex frequency, $\omega$; and for each trial frequency one integrates the eigenequations (14), subject to the boundary conditions (15), to obtain that unique eigenfunction which is well-behaved at both the center and the surface of the star. The resultant eigenfunction for each value of $\omega$ has a particular amplitude, $C_\omega^{(1)}$, for incoming gravitational waves. The purely outgoing normal modes for which one is searching are those for which the ratio

$$\frac{C_\omega^{(1)}}{B} = \left(\frac{\text{incoming wave amplitude}}{\text{central pulsation amplitude}}\right) \tag{35}$$

vanishes. Hence, one chooses one trial value of $\omega$ after another in search of the zeros of $C_\omega^{(1)}/B$.

In practice, this technique should not require very many iterations. One can use, as initial trial frequency, values obtained from the linearized theory of gravitation (see Wheeler 1962; Zee and Wheeler 1967; Chau 1967) and from the theory of radial pulsations (RSSD, chaps. iv and vii); and one can devise efficient predictor-corrector techniques for choosing successive trial frequencies in the search for zeros of $C_\omega^{(1)}/B$.
