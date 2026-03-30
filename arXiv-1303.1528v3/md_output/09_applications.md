# Applications

The I-Love-Q relations have 3 immediate applications to observational astrophysics, GWs and fundamental physics. Let us look at each application in turn.

## Observational Astrophysics

On an astrophysical front, a measurement of any single member of the I-Love-Q trio automatically provides information about the other two, even when measuring the other two directly might not be possible with current observations. For example, one might be able to measure $\bar{I}$ within 10% accuracy by measuring the orbits of binary pulsars sufficiently accurately, so as to extract the spin-orbit coupling effect in the advance rate of the periastron of the double binary pulsar J0737-3039 [lattimer-schutz,kramer-wex]. If such measurement is accomplished, one can then automatically obtain the quadrupole moment and tidal Love number of the primary pulsar by using the I-Love-Q relations. Similarly, if an equal-mass NS binary within 300Mpc is about to coalesce, one might be able to determine the $\bar{\lambda}^\mathrm{(tid)}$ of the constituents with second-generation, ground-based GW interferometers [flanagan-hinderer-love,read-love,hinderer-lackey-lang-read,lackey,damour-nagar-villain]. Then, from the I-Love-Q relations, one may obtain the moment of inertia and quadrupole moment of the binary constituents, which again, would be very difficult to measure with GWs.

Let us stress that the I-Love-Q relations cannot be used to measure the equation of state, but rather to infer two members in the I-Love-Q trio when the third is measured. The main result in this paper is, in fact, that the I-Love-Q relations seem to be rather insensitive to the EOS. Inferring the quadrupole moment and the Love number would provide important information about the properties of NSs. The quadrupole moment would tell us how much a NS can be quadrupolarly deformed (squeezed at the poles), while the Love number would tell us how much it can be deformed, for example, in the presence of a companion.

A small caveat should be presented here. The I-Love-Q relations hold for the dimensionless (barred) moment of inertia, quadrupole moment and Love number, which are normalized by the NS mass and spin. In particular, the observed NS mass differs from the mass used to normalize the I-Love-Q relations by factors of ${\cal{O}}(\chi^{2})$, ie. $M_{\rm obs} =  M_{*} \left(1 + \chi^{2} \delta M\right)$, where $\delta M = {\cal{O}}(0.3)$ [Berti:2004ny,berti-iyer-will]. For stars spinning with $\chi \lesssim 0.1$, this induces differences between $M_{\rm obs}$ and $M_{*}$ of ${\cal{O}}(10^{-3})$, which would spoil the I-Love-Q universality. However, this non-universality is much smaller than the accuracy to which $M_{*}$ can be observationally determined, and thus, it does not spoil the use of the I-Love-Q relations in observational astrophysics.

Of course, these application assumes that the universality of the I-Love-Q relations holds, which rests on the assumptions of uniform and slow-rotation, small tidal perturbations and that GR is the correct theory. Therefore, this technique cannot be applied to, for example, newly-born NSs that are differentially and rapidly rotating. NSs that source GWs in the sensitivity band of ground based detectors, however, are expected to be old, and thus uniformly rotating with large spin periods (they should have spun down by the time they are visible by GW detectors [bildsten-cutler]), so that the slow-rotation approximation is well-justified. The primary NS in the double binary pulsar has a period of $22 \; {\rm{ms}}$, which implies a $\chi \sim 0.018$, small enough that the slow-rotation approximation is again well-justified.

Millisecond binary pulsars with short-periods, ie. periods below $1 \; {\rm{ms}}$, would be spinning too fast for the above relations to be directly applicable. However, we expect I-Love-Q type universality with respect to the EoS to still hold in this case, except that now the coefficients in Table [ref:table:coeff] will also depend somewhat on the spin angular frequency (or the spin period). One can correct the universal I-Love-Q relations for non-negligible spins by considering rapidly rotating NSs [berti-stergioulas,berti-white,benhar,pappas-apostolatos], but we leave this to future work.

## Gravitational Wave Astrophysics

Another application of the I-Love-Q relations is to GW astrophysics, as a means to break the degeneracy between individual spins and the quadrupole moments of NSs in the GWs emitted during binary NS inspirals. Let us first discuss gravitational waveforms of spinning, tidally-deformed NS binaries, and then, carry out a back-of-the-envelope parameter estimation study using Fisher theory. The latter will allow us to determine the degree to which degeneracies are broken through the I-Love-Q relations and the projected accuracy to which individual NS spins could be measured given a GW detection.

### Waveforms

The sky-averaged gravitational waveform (in the Fourier domain) generated by a compact NS binary in a quasi-circular orbit with masses $m_1$ and $m_2$ and at distance $D_L$ is given by [cutlerflanagan] $\tilde{h}(f) = A(f) \exp [i \Psi(f)]$, with[^1]

$$\begin{aligned}
A(f) &=& \frac{1}{\sqrt{30} \pi^{2/3}} \frac{\mathcal{M}^{5/6}}{D_L} f^{-7/6}\,, \\
\Psi(f) &=& \Psi_\mathrm{tp}(f) + \Psi_{\bar{Q}} (f) +\Psi_{\bar{\lambda}} (f)\,.
\end{aligned}$$

Here, $f$ is the GW frequency, $\mathcal{M} =  m \eta^{3/5}$ is the chirp mass, $\eta = {m_1 m_2}/{m^2}$ is the symmetric mass ratio and $m =  m_1 + m_2$ is the total mass. The quantity $\Psi_\mathrm{tp}(f)$ is the gravitational waveform phase in the test-particle limit while $\Psi_{\bar{Q}}$ and $\Psi_{\bar{\lambda}}$ represent terms that deviate from this limit, where the former corresponds to a quadrupole moment deformation, while the latter depends on the tidal Love number.

The test-particle term, to 3.5 PN order, is given by [arun35PN,arunbuonanno,blanchet3PN]

\allowdisplaybreaks

$$\begin{aligned}
\Psi_\mathrm{tp} (f) &=& 2\pi f t_c - \phi_c - \frac{\pi}{4} +\frac{3}{128} (\pi \mathcal{M} f)^{-5/3} \Bigg \{ 1 + \left( \frac{3715}{756} + \frac{55}{9} \eta \right) x - (16\pi -4 \beta ) x^{3/2}  \\
& & + \left( \frac{15293365}{508032} + \frac{27145}{504} \eta + \frac{3085}{72} \eta^2 - 10\sigma \right) x^2 + \left( \frac{38645}{756} \pi - \frac{65}{9}\pi \eta - \gamma \right) (1+3\log{v}) x^{5/2}  \\
& & +\Bigg[ \frac{11583231236531}{4694215680} - \frac{640 \pi^2}{3} - \frac{6848}{21} \gamma_{\text{E}}  - \left(\frac{15737765635}{3048192} - \frac{2255}{12} \pi^2 \right) \eta + \frac{76055}{1728} \eta^2 - \frac{127825}{1296} \eta^3   \\
& & - \frac{6848}{21} \log (4v) + \alpha \Bigg] x^3 + \left( \frac{77096675}{254016} + \frac{1014115}{3024} \eta - \frac{36865}{378} \eta^2 \right) \pi x^{7/2} \Bigg \}\,,

\end{aligned}$$

where $x\equiv v^2 = (\pi m f)^{2/3}$ and $(t_c,\phi_c)$ correspond to the time and phase at coalescence, respectively, with $\gamma_{\text{E}}$ the Euler constant. The spin parameters $\beta$ and $\sigma$ [kidder-spin,vasuth-spinspin,arunbuonanno], $\gamma$ [arunbuonanno] and $\alpha$ [blanchet3PN] are given by[^2]

$$\begin{aligned}

\beta &=&  \left( \frac{113}{12} - \frac{19}{3}\eta \right)   \left(\hat{\boldsymbol{L}} \cdot \boldsymbol{\chi}_s \right) + \frac{113}{12} \delta_m (\boldsymbol{\chi}_s \cdot \boldsymbol{\chi}_a) \,, \\
\sigma &=&\frac{719}{48} \delta_m   \left(\hat{\boldsymbol{L}} \cdot \boldsymbol{\chi}_s \right)  \left(\hat{\boldsymbol{L}} \cdot \boldsymbol{\chi}_a \right) -\frac{233}{48} \delta_m (\boldsymbol{\chi}_s \cdot \boldsymbol{\chi}_a)  \\
& & + \left( \frac{719}{96}+\frac{1}{24} \eta \right) \left(\hat{\boldsymbol{L}} \cdot \boldsymbol{\chi}_s \right) ^2 + \left(\frac{719}{96} -30 \eta \right)  \left(\hat{\boldsymbol{L}} \cdot \boldsymbol{\chi}_a \right) ^2  \\
& & - \left( \frac{233}{96}+\frac{7}{24} \eta \right) \chi_s^2 - \left( \frac{233}{96}-10 \eta  \right) \chi_a^2\,,  \\
\gamma &=& \left(\frac{732985}{2268}-\frac{24260}{81} \eta-\frac{340}{9} \eta^2 \right) \left( \hat{\boldsymbol{L}} \cdot \boldsymbol{\chi}_s \right)   \\
& & + \left( \frac{732985}{2268}+\frac{140}{9} \eta \right) \delta_m \left( \hat{\boldsymbol{L}} \cdot \boldsymbol{\chi}_a \right)\,, \\
\alpha &=& \frac{2270 \pi}{3} \left[ \left(1-\frac{227}{156} \eta \right) \left( \hat{\boldsymbol{L}} \cdot \boldsymbol{\chi}_s \right)+ \delta_m \left(  \hat{\boldsymbol{L}} \cdot \boldsymbol{\chi}_a \right) \right]\,,  \\
\end{aligned}$$

where $\hat{\boldsymbol{L}}$ is the unit orbital angular momentum, $\delta_{m} \equiv (m_{1} - m_{2})/m$ is the dimensionless mass difference, $\boldsymbol{\chi}_s \equiv (\boldsymbol{\chi}_1 + \boldsymbol{\chi}_2)/2$ and $\boldsymbol{\chi}_a \equiv (\boldsymbol{\chi}_1 - \boldsymbol{\chi}_2)/2$ with $\boldsymbol{\chi}_i \equiv \boldsymbol{S}_i/m_i^2$ denoting the dimensionless spin vector of the $i$-th body. Notice that we are here referring to the individual NS spin vectors by $\boldsymbol{S}_{i}$

The quadrupole moment contribution correction to the test-particle limit in the GW phase enters at 2PN order and it is given by [poisson-quadrupole,vasuth-spinspin]

$$\begin{aligned}
\Psi_{\bar{Q}} (f) &=& \frac{3}{128} \frac{x^{-5/2}}{\eta} \left\{ -50  \left[  \left( \frac{m_1^2}{m^2} \chi_1^2 + \frac{m_2^2}{m^2} \chi_2^2 \right) (\bar{Q}_s -1) \right. \right.  \\
&& \left. \left. +\left( \frac{m_1^2}{m^2} \chi_1^2 - \frac{m_2^2}{m^2} \chi_2^2 \right) \bar{Q}_a \right] x^2 \right\}\,,
\end{aligned}$$

where

$$
\bar{Q}_s \equiv \frac{\bar{Q}_1 + \bar{Q}_2}{2}\,, \quad \bar{Q}_a \equiv \frac{\bar{Q}_1 - \bar{Q}_2}{2}\,.
$$

$\bar{Q}_s$ is strongly correlated with $\sigma$, which enters at the same PN order as $\bar{Q}_s$.

The leading-order contribution of $\Psi_{\bar{\lambda}} (f)$ to the GW phase enters at 5PN order through [flanagan-hinderer-love]

$$\begin{aligned}
\Psi_{\bar{\lambda}}^\mathrm{5PN} (f) &=& -\frac{3}{128} \frac{x^{-5/2}}{\eta}24 \left[ (1+7\eta -31 \eta^2 ) \bar{\lambda}_s \right.   \\
& & \left. + (1+9 \eta -11 \eta^2 ) \bar{\lambda}_a \delta_m \right] x^{5}\,,
\end{aligned}$$

where

$$
\bar{\lambda}_s \equiv \frac{\bar{\lambda}^\mathrm{(tid)}_1+\bar{\lambda}^\mathrm{(tid)}_2}{2}\,, \quad \bar{\lambda}_a \equiv \frac{\bar{\lambda}^\mathrm{(tid)}_1-\bar{\lambda}^\mathrm{(tid)}_2}{2}\,.
$$

Higher PN contributions to $\Psi_{\bar{\lambda}} (f)$ can be found in [vines1,vines2,damour-nagar-villain]. When carrying out parameter estimation studies, as explained below, we will use $\Psi_{\bar{\lambda}} (f)$ as given in [damour-nagar-villain], which includes up to 2.5PN order corrections relative to $\Psi_{\bar{\lambda}}^\mathrm{5PN}$[^3].

[^1]: $A(f)$ needs to be multiplied by $\sqrt{3}/2$ when calculating the Fisher matrix for LISA and DECIGO/BBO.
[^2]: $\sigma$ includes the quadrupole-monopole interaction in the test-particle limit.
[^3]: Tidal effects on the gravitational waveform phase have been calculated to 1.5PN order relative to the leading 5PN contribution, in addition to tail effects at 2.5PN order. Ref. [damour-nagar-villain] estimated that currently unknown terms should be subdominant, at least for an equal-mass binary.

### Parameter Estimation

For stationary and Gaussian detector noise, the measurement accuracy of parameters $\theta^a$ can be estimated as

$$
\Delta \theta^a = \sqrt{\frac{(\Gamma^{-1}){}^{aa}}{N}}\,,
$$

where $N$ is the number of effective interferometers and

$$
\Gamma_{ab} \equiv 4 \; \mathrm{Re} \int^{f_\mathrm{max}}_{f_\mathrm{min}} \frac{\partial_a \tilde{h}(f) \partial_b \tilde{h}(f)}{S_n(f)} df
$$

is the Fisher matrix, where the partial derivatives are with respect to the parameters $\theta^a$. The noise spectral density $S_n(f)$ is given in Refs. [cornish-PPE,mishra,yagi:brane] for Adv. LIGO, ET and DECIGO/BBO, respectively. We take the lower cutoff frequencies $f_\mathrm{min}$ to be 10Hz for Adv. LIGO, 1Hz for ET and the frequency 1yr before coalescence for DECIGO/BBO. For Adv. LIGO and ET, we take the higher cutoff frequency to be that of the innermost stable circular orbit (ISCO), $f_\mathrm{max}=f_\mathrm{ISCO}=1/(6^{3/2} \pi m)$, while we set $f_\mathrm{max}=100$Hz for DECIGO/BBO. $N$ is the number of effective interferometers, which we take to be 5 for 2nd-generation ground-based detectors (corresponding to 2 Adv. LIGO, Adv. VIRGO, KAGRA and INLIGO), 2 for ET (like LISA [cutler1998]) and 8 for DECIGO/BBO [yagi:brane].

We focus here on GWs emitted during the quasi-circular inspiral of NSs with aligned spins, since this is a realistic astrophysical scenario [kesden-berti]. Given that the NS masses are expected to be approximately the same, we will not include $\bar{Q}_a$ and $\bar{\lambda}_a$ in the parameter vector, as this must be close to zero. We will consider the case of slightly unequal NS masses. We choose two parameterizations of the waveform. Parameterization $A$ uses the parameter set [damour-nagar-villain]

$$
\{ \theta_A^i \} =  (\ln \mathcal{M}, \ln\eta, \beta, D_L, t_c, \phi_c, \bar{\lambda}_s)

$$

with the priors $|\eta| < 0.25$ and $|\beta| < 0.8$. We do not include $\sigma$ in this set because the NS spins at the time of coalescence are expected to be small [bildsten-cutler][^1].
Parametrization $B$ uses the parameter set

$$
\{ \theta_B^i \} =  (\ln \mathcal{M}, \delta_m, \chi_s, \chi_a, D_L, t_c, \phi_c, \bar{Q}_s(\bar{\lambda}_s), \bar{\lambda}_s)\,,

$$

with the priors $|\delta_m| < 1/3$, $|\chi_s| < 0.1$ and $|\chi_a| < 0.1$. Moreover, we use the Q-Love relation to express $\bar{Q}_s$ in terms of $\bar{\lambda}_s$, and thus, partially break the degeneracy between $\bar{Q}_s$ and $\chi_s$.

Figure [ref:fig:spin] shows the measurement accuracies of spin parameters using second-generation, ground-based detectors. We assume that the detected GW was emitted by a source at $D_L=100$Mpc with ${\rm{SNR}} \sim 30$. We consider 3 different systems:  (i) $(m_1,m_2)=(1.45,1.35)M_\odot$, $\chi_1=\chi_2$, (ii) $(m_1,m_2)=(1.45,1.35)M_\odot$, $\chi_1=2 \chi_2$ and (iii) $(m_1,m_2)=(1.4,1.35)M_\odot$, $\chi_1=\chi_2$.

Observe, that the averaged spin $\chi_s$ can be measured to $\mathcal{O}(10)%$. Such an accuracy on $\chi_s$ is inaccessible without the Q-Love relation.

We can understand this enhanced accuracy in the extraction of $\chi_{s}$ as follows. First, notice that, for an equal-mass and spin-aligned binary, $\beta \sim \mathcal{O}(10) \chi_s$ [see e.g. Eq. (Eq. beta)]. Given that the measurement accuracy of $\beta$ is $\Delta \beta = \mathcal{O}(0.1)$, this implies a measurement accuracy of $\chi_{s}$ of $\Delta \chi_s \approx 0.01$, which corresponds to $\Delta \ln \chi_s \approx 0.1$ for $\chi_s \approx 0.1$.  The measurement accuracy of $\chi_s$ for system (ii) increases as $\chi_1 \to 0.1$, as shown in Fig. [ref:fig:spin]. This is because $\partial \tilde{h}/\partial \delta_m \approx 0 \approx \partial \tilde{h}/\partial \chi_a$ as $\chi_1 \approx 0.1$ and the priors lead to $\delta_m$, $\chi_a$ and other parameters being effectively uncorrelated. In such a case, however, the assumptions that underlie the Fisher approximation may be violated [vallisneri-fisher], and hence, one requires a Bayesian analysis [cornish-PPE] to confirm these results. Finally, we have checked that the Q-Love relation does not improve the measurement accuracy of $\bar{\lambda}_s$.

![](Deltam-lambda-PRD.eps)

Up until now, we have considered equal-mass NS binaries, but realistic systems may not have identical masses. If the NS masses are not equal, one must then take into account the parameters $\bar{Q}_a$ and $\bar{\lambda}_a$, whose inclusion could in principle degrade the accuracy to which other parameters are extracted. Let us then study the range of masses for which neglecting $\bar{Q}_a$ and $\bar{\lambda}_a$ is a good approximation. A rough estimate of this range can be obtained by investigating the systems for which the accumulated GW phase induced by terms proportional to $\bar{Q}_a$ and $\bar{\lambda}_a$ is less than one radian. Let us then define the latter by $\Psi_{\bar{Q}_a}$ and $\Psi_{\bar{\lambda}_a}$ respectively, where

$$\begin{aligned}
\Psi_{\bar{Q}_a} (f) &=& - \frac{75}{64} \frac{1}{\eta} \left[ \left( \frac{m_1^2}{m^2} \chi_1^2 - \frac{m_2^2}{m^2} \chi_2^2 \right) \bar{Q}_a \right] x^{-1/2}\,,
 \\
\Psi_{\bar{\lambda}_a}^\mathrm{5PN} (f) &=&  -\frac{9}{16} \frac{1}{\eta} (1+9 \eta -11 \eta^2 ) \bar{\lambda}_a \delta_m  x^{5/2}\,,
\end{aligned}$$

to leading PN order.

Figure [ref:fig:deltam] shows the range of masses for which $\Psi_{\bar{\lambda}_a} = 1$ for various realistic EoSs. Systems above these lines would lead to  $\Psi_{\bar{\lambda}_a} > 1$, while those below this line lead to  $\Psi_{\bar{\lambda}_a} < 1$. In particular, for systems that satisfy the latter inequality, we can in principle neglect $\bar{\lambda}_a$ if the SNR is $\mathcal{O}(10)$. This figure implies that for second-generation, ground-based detectors, the parameter estimation study presented above is probably valid, even for unequal-mass systems provided, for example, that $\bar{m}=1.4 M_\odot$ and $\Delta m \lesssim \mathcal{O}(0.1)M_\odot$, where $\bar{m} \equiv (m_1 + m_2)/2$ is the averaged mass of the binary and $\Delta m \equiv |m_1 - m_2|$ is the mass difference. This same conclusion also applies to neglecting $\bar{Q}_{a}$, provided $|\chi| < 0.1$.

[^1]: Damour *et al.* [damour-nagar-villain] estimated that at the time of coalescence, $|\beta| < 0.2$ and $|\sigma | < 10^{-4}$. We use a conservative prior $|\beta| < 0.8$ which corresponds to $|\chi|<0.1$.

## Fundamental Physics

The independent measurement of any 2 members of the I-Love-Q trio would allow us to perform model-independent and EoS-independent tests of GR.  For example, if one can measure $\bar{I}$ and $\bar{\lambda}^\mathrm{(tid)}$ independently, one can plot a point in the I-Love plane with an error box. If the I-Love relation in GR crosses the error box, then GR is consistent with the observations. Otherwise, one would have found model-independent evidence for some type of departure from GR. Moreover, one can constrain non-GR theories by requiring that the I-Love relation in that theory crosses the error box.

The accuracy of such a test depends, of course, in how accurately two elements in the I-Love-Q set can be measured. One way to measure $\bar{I}$ would be to look for a spin-orbit correction to the rate of advance of the periastron of a binary system. Future double binary pulsar observations may measure $\bar{I}$ with an accuracy of roughly 10% [lattimer-schutz,kramer-wex]. Probably, the best way to measure $\bar{\lambda}^\mathrm{(tid)}$ and $\bar{Q}$ would be to use GW observations.

In what follows, we first discuss the possibility of measuring $\bar{Q}$ and $\bar{\lambda}^\mathrm{(tid)}$ simultaneously, given GW observations. Then, we discuss how well GR tests can be carried out by combining GW observations with binary pulsar observations. For concreteness, we apply all of this to a specific modified gravity theory (dynamical CS gravity [CSreview]).

### Redundancy Tests with GW Observations Only

{\renewcommand{1.2}{1.2}
| \noalign{\smallskip} |  | \multicolumn{4}{c}{$\bar{\lambda}^\mathrm{(tid)}=400.0$} |  | \multicolumn{3}{c}{$M_*=1.3382M_\odot$} |
|---|---|---|---|---|
| EoS |  | \multicolumn{1}{c}{$M_*$} | \multicolumn{1}{c}{$\mathcal{R}_*$} | \multicolumn{1}{c}{$\bar{Q}$} | \multicolumn{1}{c}{$f_\mathrm{spin}$} |  | \multicolumn{1}{c}{$\mathcal{R}_*$} | \multicolumn{1}{c}{$\bar{I}$} | \multicolumn{1}{c}{$\bar{\lambda}^\mathrm{(tid)}$} |
|  |  | ($M_\odot$) | (km) |  | (Hz) |  | (km) |  |  |
| APR |  | 1.40 | 12.2 | 5.52 | 194 |  | 12.2 | 13.3 | 520 |
| SLy |  | 1.32 | 11.6 | 5.54 | 206 |  | 11.6 | 12.2 | 375 |
| LS220 |  | 1.38 | 13.5 | 5.56 | 198 |  | 13.6 | 13.1 | 506 |
| Shen |  | 1.55 | 14.6 | 5.54 | 176 |  | 15.0 | 15.9 | 1012 |
| \noalign{\smallskip} |
}
Let us estimate how accurately future ground-based detectors may determine $\bar{\lambda}^\mathrm{(tid)}$ and $\bar{Q}$ *simultaneously*.  The measurability of $\bar{\lambda}^\mathrm{(tid)}$ has been discussed extensively in [flanagan-hinderer-love,read-love,hinderer-lackey-lang-read,kiuchi,kyutoku,hotokezaka,vines1,vines2,lackey,damour-nagar-villain], while the effect of the quadrupole moment on compact binary GWs has been estimated in [bildsten-cutler,lai-rasio-shapiro,poisson-quadrupole], but so far no study has been performed to study their simultaneous extraction. Let us consider an equal-mass NS binary with $\bar{\lambda}^\mathrm{(tid)} = 400$, which corresponds to NSs with $M_*=1.4M_\odot$ with the APR EoS (other parameters are shown in Table [ref:table:parameters]).

![](error-dup.eps)

Figure [ref:fig:error] shows the measurement accuracy of $\bar{Q}_s$ and $\bar{\lambda}_s$ with second-generation, ground-based detectors, ET and DECIGO/BBO. We assume an equal-mass, spin-aligned NS/NS binary with $(1.4, 1.4)M_\odot$, $\chi_1 < 0.1$ and $\chi_2 =0$ at $D_L=100$Mpc, and we also assume that the APR EoS is the correct one. The measurement accuracy of $\Delta \bar{Q}_s$ with second-generation, ground-based detectors is $\Delta \ln \bar{Q}_s \approx 30$, which would increase to  $\Delta \ln \bar{Q}_s \approx 2 $ using future detectors, such as ET or DECIGO/BBO. These results imply that it may be difficult to measure $\Delta \bar{Q}_s$ due to its strong degeneracies with spin parameters. Notice, however, that even though one may not be able to detect $\bar{Q}_s$, one can still place an upper and lower bound on $\bar{Q}_s$. Such a bound would be sufficient to perform model-independent GR tests. The measurement accuracy of $\bar{\lambda}_s$ with second-generation, ground-based detectors is $\Delta\bar{\lambda}_s \approx 0.8$, which would increase by roughly an order of magnitude using ET. Although the error bars are large with current detectors, it may be possible to measure $\bar{\lambda}_s$ with future GW observations.

Given the above measurement errors, we can now simulate a GR test. Figure [ref:fig:QLoveerror] presents the Q-Love relation for realistic EoS, together with a fiducial GW measurement of the pair $(\bar{Q},\bar{\lambda}^\mathrm{(tid)})$ and its estimated errors. Notice that the error in $\bar{Q}$ is larger than the value about which the error is centered. This implies that a GW measurement would not be able to measure $\bar{Q}$, but it would be able to say the region of allowed $\bar{Q}$ that is consistent with the GW detection. Therefore, such a GW detection would automatically constitute a model-independent test of GR; for GR to be consistent with these measurements, the GR Q-Love curve must cross the GW error box in $(\bar{Q},\bar{\lambda}^\mathrm{(tid)})$.

![](Q-Love-error-lambda400-PRD.eps)

The above test is quite robust. First, although we employ a uniform and slow-rotation approximation, the NS/NS binaries that ground-based detectors will observe will have spun down by the time they enter the detector’s sensitivity band, and thus, the slow-rotation approximation should be excellent. Second, the error box of Fig. [ref:fig:QLoveerror] depends on how accurately $(\bar{\lambda}_s, \bar{Q}_s)$ can be measured, which in turn depends on whether $(\bar{\lambda}_{a}, \bar{Q}_{a})$ need to be included in the parameter set. This would be the case if the binary system detected were not an equal-mass one. As shown in Fig. [ref:fig:deltam], however, there is a wide range of mass ratios for which these parameters can be neglected, even outside of the equal-mass point; thus, the discussion presented above should be robust.

### Joint Tests with GW and Electromagnetic Observations

Since $\bar{Q}$ is a quantity that is difficult to measure with GW observations, let us consider model-independent and EoS-independent tests of GR with the I-Love relation that uses a combination of GW and double binary pulsar observations. Let us then assume that $\bar{I}$ has been measured to $10%$ by future double binary pulsar (J0737-3039) observations [lattimer-schutz,kramer-wex], and that $\bar{\lambda}^\mathrm{(tid)}$ has been measured to $40%$ with future GW observations. The latter assumes an ET detection of an equal-mass, non-spinning NS/NS binary at 3Gpc, with the individual NS masses exactly equal to that of the primary pulsar in J0737-3039, $M_*=1.3382M_\odot$, assuming the Shen EoS[^1]. All of this is shown in Fig. [ref:fig:CS], together with the fiducial measurement of $(\bar{I},\bar{\lambda}^\mathrm{(tid)})$ as a big black cross. As shown in that figure, one can constrain modified theories of gravity, such as dynamical CS gravity, by requiring that the I-Love curve crosses the error box.

The test described above has one major problem: the NS mass $m_\mathrm{pulsar}$ of the primary pulsar in J0737-3039 and the individual NS masses in the binary system that generated the detected GW will all in principle be different from each other. As explained in Sec. [ref:sec:parameter-estimation], the accuracy to which $\bar{\lambda}_{s}$ can be measured assumed that $\bar{Q}_{a}$ and $\bar{\lambda}_{a}$ could be neglected, which holds for a certain range of mass differences $\Delta m$, shown in Fig. [ref:fig:deltam]. In general, the typical maximum mass difference would need to be $\Delta m = \mathcal{O}(0.1)M_\odot$, assuming observations with $\mathrm{SNR} \approx 10$. Given current event rate estimates, one expects to detect NS/NS binaries with such similar masses, and thus, this is not in principle a problem. The test above, however, also requires that $m_\mathrm{pulsar} \approx \bar{m}$, since after all the I-Love-Q relations assume one is investigating NSs with the same mass.

Let us then estimate how much the I-Love-Q relations would change if $m_\mathrm{pulsar} \neq \bar{m}_\text{GW}$. The top panel of Fig. [ref:fig:I-Love-massratio] shows the I-Love relation with $m_\mathrm{pulsar}/\bar{m}_\text{GW}=0.9$, $1.0$ and $1.1$ for realistic EoSs, while the bottom panel shows the relative fractional difference between the I-Love curves for different EoSs and the APR EoS as a reference. The relative fractional difference when $m_\mathrm{pulsar}/\bar{m}_\text{GW}=0.9$ (not shown in this figure) is similar to that of $m_\mathrm{pulsar}/\bar{m}_\text{GW}=1.1$. One sees that the dependence on the EoS becomes stronger as the mass difference $m_\mathrm{pulsar}$ and $\bar{m}_\text{GW}$ increases. However, this dependence is still weak if the mass difference is sufficiently small (or order $0.1 M_{\odot}$), and most importantly, the loss of universality (the difference between curves with different EoS) is much, much smaller than the observational error in measuring either the moment of inertia or the Love number. Therefore, one can perform the GR test described above, even when $m_\mathrm{pulsar} \neq \bar{m}_\text{GW}$.

![](I-Love-mass-ratio-Dup-PRD.eps)

Of course, the test described here assumes that the uniform and slow-rotation approximation used to derive the I-Love-Q relation holds for binary pulsars. This is indeed the case, provided the period is sufficiently long, such that each binary component is spinning slowly. However, the approximation might break down for (currently unobserved) sub-millisecond pulsars, ie. those with periods shorter than $1 \; {\rm{ms}}$. For such systems, the I-Love-Q relations will also now depend on the spin frequency. A cursory analysis, however, suggests that the spin-frequency effect breaks universality at the $10%$ level [berti-stergioulas,berti-white,benhar,pappas-apostolatos]. Therefore, the difference in I-Love-Q relations for different EoSs will be rather small, and in particular smaller than the errors in the first binary pulsar measurements of the moment of inertia.

Binary pulsar observations may measure $\bar{I}$ within 10% accuracy, but this amazing measurement will be difficult to accomplish in the *near* future [kramer-wex]. This is because the effect of the moment of inertia (or equivalently, the spin-orbit coupling) in the motion of the binary is of $\mathcal{O}(10^{-5})$ relative to the leading-order contribution. This means that one needs to measure at least 3 post-Keplerian parameters to this accuracy, in order to determine the two masses and the moment of inertia. Of course, such a measurement is a big challenge, although not out of the question, as binary pulsar observations improve within the next ten years.

![](Love-C-error-2.eps)

Alternatively, one could use the NS compactness $C$ instead of $\bar{I}$ to perform model-independent GR tests. Currently, $C$ has been measured to $\mathcal{O}(10)%$ accuracy with low-mass X-ray binary observations [steiner-lattimer-brown,ozel-baym-guver,ozel-review,guver]. Figure [ref:fig:Love-C] shows the $\bar{\lambda}^\mathrm{(tid)}$--$C$ relation for realistic EoSs and for the $n=1$ polytrope. We also show in this figure a fiducial measurement of $(C,\bar{\lambda}^\mathrm{(tid)})$ (big black cross), as well as projected measurement accuracies (dashed lines). For the latter, we assume $\Delta C =0.05$ from electromagnetic observations of a NS with $M_*=1.4 M_\odot$ and a roughly $70%$ measurement of $\bar{\lambda}^\mathrm{(tid)}$ from GW detection (equal-mass, non-spinning NS/NS binary with the individual NS mass of $M_* = 1.4M_\odot$ at $D_L=100$Mpc) with second-generation, ground-based detectors[^2]. Although the dependence on the EoS is relatively large compared to the universality of the I-Love-Q relations, the measurement errors are {*larger*} than the uncertainties due to the EoS. This shows that one might be able to use the Love-C relation to perform model-independent tests of gravity.

[^1]: Adv. LIGO is expected to detect NS/NS binaries out to $D_L \approx 300$Mpc with the detection rate of $\mathcal{O}(10)$/yr [abadie]. Therefore, if we consider ET detecting GW signals from a NS/NS binary at 3Gpc, the expected detection rate would be $\mathcal{O}(10)/\mathrm{yr}\times 10^3 \sim \mathcal{O}(10^4)$/yr. With this detection rate, we may detect an equal-mass NS/NS binary, with the individual masses very close to that of the primary pulsar of J0737-3039, $M_*=1.3382M_\odot$
[^2]: The measurement error of $\Delta \bar{\lambda}^\mathrm{(tid)} \approx 0.7$ is slightly better than that shown in Fig. [ref:fig:QLoveerror]. This is because we assumed parameter set A [Eq. (Eq. parameter-A)] instead of B [Eq (Eq. parameter-B)].

### Example: Dynamical CS Gravity

We now apply the results obtained above to see how testing a specific theory of gravity would go about. As an example, we choose dynamical CS gravity [jackiw,CSreview], which is well-motivated from the Standard Model, superstring theory [polchinski1,polchinski2], loop quantum gravity [alexandergates,taveras,calcagni] and inflation [weinberg-CS]. Dynamical CS gravity is a parity-violating, quadratic-curvature theory, where the Einstein-Hilbert action is modified through the Pontryagin density (the contraction of the Riemann tensor and its dual), coupled to a dynamical scalar field. This theory has a characteristic length $\xi^{1/4}$, which has been constrained by Solar System experiments, using Gravity Probe B [GPB] and LAGEOS [LAGEOS], to $\xi^{1/4} < \mathcal{O}(10^8)$km [alihaimoud-chen]. Dynamical CS gravity should be treated as an {*effective*} theory, and thus, one should work to leading-order in a small coupling expansion, ie. to leading order in the dimensionless coupling constant $\zeta \equiv \xi M_*^2/\mathcal{R}_*^6$ [kent-CSNS].

NSs in dynamical CS gravity have been studied before. Reference [kent-CSNS] found that it would be difficult to meaningfully constrain this theory with binary pulsar observations in the standard fashion. This is because the largest CS correction appears in the rate of change of periastron advance at 1PN order, and thus it is suppressed by the ratio of the binary’s mass to its separation (for J0737-3039, this is of ${\cal{O}}(10^{-6})$). The CS correction to the NS moment of inertia
was calculated in [yunespsaltis,alihaimoud-chen], while the CS correction to the NS quadrupole moment was obtained in [kent-CSNS]. In the small coupling approximation, the CS corrections to $\bar{I}$ and $\bar{Q}$ scale linearly with $\zeta$. The leading-order CS correction to tidal effects enters through the gravitomagnetic tidal tensor (because of the parity) [quadratic], and hence $\bar{\lambda}^\mathrm{(tid)}$ is not affected at leading-order.

Figure [ref:fig:CS] shows the I-Love relation in dynamical CS gravity with a fixed value of the coupling constant $\xi = 1.85 \times 10^4 M_*^4$. Observe that the dependence on the EoS is stronger than that of the GR I-Love relation. We believe that this is because a compact object in dynamical CS gravity depends on the scalar-dipole charge, which encodes information on the internal structure of the body [kent-CSNS]. Given this reasoning, we expect that the I-Love-Q relations should be more sensitive to the NSs’ internal structure in dynamical CS gravity than in GR. The bottom panel of Fig. [ref:fig:CS] shows that there is indeed a loss of universality, but the latter is still preserved to a few % level.

With this in hand, let us estimate the projected bound that one could place on dynamical CS gravity using the I-Love relation. With $\xi=1.85 \times 10^4 M_*^4$, the I-Love curve in dynamical CS gravity barely crosses the error box. Since the larger $\xi$, the higher the CS curves, such an I-Love observation would automatically constrain $\xi < 1.85 \times 10^{4} M_{*}^{4}$, which corresponds to $\zeta = 0.0977$. Converting back to dimensional quantities, such a test would impose the constraint

$$
\xi^{1/4} < \mathcal{O}(50) \mathrm{km}\,.

$$

NS observations would then allow us to probe the theory within NS length scales, like the NS radius. Notice that the above bound is stronger than Solar System [alihaimoud-chen] and table-top [kent-CSBH] ones by more than six orders of magnitude. Notice, however, that this bound is slightly weaker than the proposed projected bound with GW observations of BH/BH binaries [kent-CSGW]. This is because the “radius” of a BH is smaller than that of a typical NS, and thus, with the former, we can probe shorter length scales. Notice also that the bound given above is dominated by the measurement error on the NS moment of inertia. This is reasonable because the tidal Love number is unaffected in dynamical CS gravity to the order of approximation considered here.

The measurement accuracy shown in Fig. [ref:fig:CS] is obtained by assuming that the Shen EoS is the correct one, which gives the weakest bound on the theory among the realistic EoSs considered in this paper. This is because, with the NS masses fixed to $M_*=1.3382M_\odot$, the Shen EoS gives the largest $\bar{\lambda}^\mathrm{(tid)}$ (see Table [ref:table:parameters]). Figure [ref:fig:CS] shows that the deviation away from GR becomes larger for smaller $\bar{\lambda}^\mathrm{(tid)}$. This is because the compactness becomes larger for smaller $\bar{\lambda}^\mathrm{(tid)}$ (or smaller $\bar{I}$ and $\bar{Q}$) which allows us to probe stronger gravity. These studies suggest that the I-Love-Q relations can be very powerful in testing GR in the strong-field regime.