   



**An Introduction to CLEED**

... and (almost) everything else  
you need to know about carrying out  
structural analysis by LEED  

**Georg Held**

**2025-10-02**



## Preface

This manual is still very much work in progress. It is a compilation of notes, lectures, and instructions that I have written over the years and I have been adding and correcting whenever I find a little time. Some of the non-essential introductory chapters are still incomplete or missing. This is information that can be found in most textbooks that cover LEED. The description of the input files for the CLEED package, however, as well as the instructions for installation, should be complete and accurate. If you spot any errors or omissions, please let me know.

   
Georg Held (georg.held.@.diamond.ac.uk)

# Introduction: History and Principles

## History

### General

When Davisson and Germer reported in 1927 that the elastic scatting of low-energy electrons from well ordered surfaces leads to diffraction spots similar to those observed in X-ray diffraction , this was the first experimental proof of the wave nature of electrons. A few years before, in 1923, De Broglie had postulated that electrons have a wave length in Å of

``` math
\begin{equation}
\lambda_e =  \frac{h}{m_e v} = \sqrt{\frac{150 \mbox{eV}}{E_{kin}}}
\end{equation}
```

and a corresponding wave vector of the length

``` math
\begin{equation}
k = \frac{2 \pi}{\lambda_e}
\end{equation}
```

where $`h`$ is Planck’s constant, $`m_e`$ the electron mass, $`v`$ the velocity, and $`E_{kin}`$ the kinetic energy of the electron. Already Davisson and Germer realised that the diffraction of low-energy electrons (LEED) can be used to determine the structure of single crystal surfaces in analogy to X-ray diffraction if their kinetic energy is between 40 and 500 eV, i.e. their wave length ranges between 0.5 and 2 Å. Due to their small inelastic mean free path (IMFP) of typically less than 10 Å electrons in this energy range sample only the top-most atomic layers of a surface and are, therefore, better suited for the analysis of surface geometries than X-ray photons which have a much larger mean free path (typically a few $`\mu`$m). However, unlike for photon diffraction, multiple scattering plays an important role in the diffraction process of electrons at solid surfaces. Therefore, the analysis of LEED data with respect to the exact positions of atoms at the surface is somewhat more complicated and requires fully dynamical quantum mechanical scattering calculations.

The use of LEED for surface analysis became important when large enough single crystals became available for surface studies. At first the technique was only used for qualitative characterisation of surface ordering and quantitative determination of the two-dimensional surface lattice parameters (e.g. superstructures, see below). The information about the positions of the atoms in the surface is hidden in the energy-dependence of the diffraction spot intensities, the so-called LEED I-V, or $`I(E)`$, curves. In the late 1960’s computer programs became available which could perform fully dynamical scattering calculations for simple surface geometries. The comparison of such theoretical I-V curves for a set of model geometries with experimental data allows the determination of the atomic positions within the surface by trial and error. With the immense growth of available computer power and speed since then, LEED I-V structure determination could be applied to a large number of more and more complex surface geometries. This has made LEED a standard technique of modern surface crystallography.

To date, the structure analysis by LEED I-V is the most accurate and reliable way of determining the atomic positions at surfaces. Especially the fact that LEED is sensitive to the atomic positions not only in the top-most layer but down to several layers below the surface makes it an ideal tool for studying spontaneous or adsorbate-induced surface reconstructions. General restrictions are that LEED usually requires ordered and conducting surfaces, otherwise charging would distort the LEED pattern. The time needed for calculating the I-V curves and the number of trial geometries are factors limiting the complexity of accessible surface structures on the computational side; the density of LEED spots on the fluorescent screen is an experimental factor limiting the size of unit cells that can be studied.

Since LEED theory was initially developed for close packed clean metal surfaces, these are the most reliably determined surface structures with theory-experiment agreement comparable to that between two experimental sets of I-V curvesand error bars for the atomic coordinates as small as 0.01 Å. A good overview over state-of-the-art LEED structure determinations of clean metal surfaces and further references are found in two recent articles by Heinz et al. . For more open adsorbate covered and/or reconstructed surfaces certain approximations used in the standard programs are less accurate which leads to higher RP factor values. For simple superstructures of mono-atomic adsorbates or small molecules on metal the agreement is usually not as good and worse for large molecules, which often adsorb in complex superstructures and semiconductor surfaces. Since the error margins scale with the agreement between experimental data and theoretical calculations, the accuracy of atomic coordinates is these latter cases is smaller with error bars up to 0.1 Å. A review of structure determinations of molecular adsorbates was published by Over ; the structure determination of semiconductor surfaces was reviewed by Kahn . In general the accuracy is higher for coordinates perpendicular to the surface than for lateral coordinates. For further details about the history, experimental setup, and theoretical approaches of LEED refer to the books by Pendry, , Van Hove and Tong , Van Hove, Weinberg and Chan , and Clarke . This chapter relies extensively on these books.

### CLEED package

The CLEED project was started by G. Held in 1994 when working in D. A. King’s group in Cambridge, UK, out of a certain degree of frustration. Then all available LEED codes were FORTRAN based, which did at the time not allow dynamical memory allocation, and, hence, the codes had to be edited and re-compiled for every new surface structure. The aim was to provide a more user-friendly program package that can handle most single crystal surface structures without the necessity of re-compiling and would determine most non-geometrical parameters, such as the beams to be used in the calculation, matrix size, distribution of atoms into layers, etc., automatically, thus reducing the input to - more or less - only the surface geometry.

It was decided to write the code from scratch, based on the well-established layer doubling and combined space methods by Pendry and Van Hove & Tong, respectively . The C programming language was used for several reasons:

- C allows dynamical memory allocation and, thus, diffraction matrices can be adapted easily to the surface geometry.

- C allows ”system calls” of other programs from within a master program and, thus, a modular and very flexible program structure.

- C compilers are an intrinsic part of any Linux distribution and, thus, freely available on every PC workstation.

Initially, the name ”CLEED” was just a working title originating from the programming language used. Later it was decided to that ”C” should stand for Cambridge, however major contributions to the package were added since by W. Braun in Erlangen, Germany, (1998 - 2001), P. de Andres and M. Blanco-Rey in Madrid, Spain, (2000 - 2005) and Z. V. Zheleva in Reading, UK, (2007 - 2011).

## LEED pattern

### Surface Lattice

Low-energy electrons do not penetrate into the crystal bulk far enough to experience its three-dimensional periodicity, therefore the diffraction pattern is determined by the two-dimensional surface periodicity described by the lattice vectors $`\vec{a}_1`$ and $`\vec{a}_2`$, which are parallel to the surface plane. A general lattice point within the surface plane is an integer multiple of these lattice vectors:

``` math
\begin{equation}
\vec{R} = n_1 \vec{a}_1 + n_2 \vec{a}_2
\end{equation}
```

Note, for every periodic lattice there is an infinite number of pairs of lattice vectors that fulfil the above condition and for many lattices the choice is not obvious. Usually the following conventions are used:

- The angle $`\gamma`$ between $`\vec{a}_1`$ and $`\vec{a}_2`$ is between $`90^\circ`$ and $`180^\circ`$ and as close to $`90^\circ`$ as possible.

- $`|\vec{a}_1| \geq |\vec{a}_2 |`$

- $`\vec{a}_1`$ and $`\vec{a}_2`$ form a right-handed system, i.e. going anticlockwise from $`\vec{a}_1`$ to $`\vec{a}_2`$ is a smaller angle than going clockwise (like the angle between thumb and index finger of the right hand when facing the palm)

The left-hand column of Figure [fig_int_LEED_patt] shows the conventional lattice vectors for the most common square, rectangular, and hexagonal lattices. The vectors $`\vec{a}_1`$ and $`\vec{a}_2`$ can be combined to form the lattice matrix:
``` math
\begin{equation}
\cal A
=  \left(
\begin{array}{c}
\vec{a}_1 \\ \vec{a}_2
\end{array} \right)
= \left(
\begin{array}{cc}
a_{1x} & a_{1y} \\
a_{2x} & a_{2y}
\end{array} \right)
\end{equation}
```
The determinante $`\det {\cal A}  = \vec{a}_1 \times \vec{a}_2  =  a_{1x}a_{2y} - a_{1y}a_{2x}`$ is equal to the area $`A_A`$ of the unit cell described by $`\vec{a}_1`$ and $`\vec{a}_2`$. $`\det \cal A`$ is positive if $`(\vec{a}_1, \vec{a}_2)`$ form a right-handed system.

<figure id="fig_int_LEED_patt" data-latex-placement="hb">
<img src="docs/Fig_LEED_Patt_1.jpg" style="height:30.0%" />
<img src="docs/Fig_LEED_Patt_2.jpg" style="height:30.0%" />
<figcaption> Examples of LEED patterns for common surfaces: fcc{100}, fcc{110}, fcc{111}, fcc{100}-<span class="math inline"><em>p</em>(2 × 1)</span>. </figcaption>
</figure>

### Reciprocal Lattice

The two-dimensional Bragg condition leads to the definition of reciprocal lattice vectors $`\vec{a}_1^\ast`$ and $`\vec{a}_2^\ast`$ which fulfil the following set of equations:
``` math
\begin{equation}
\begin{array}{lclcl}
\vec{a}_1 \cdot \vec{a}_1^\ast $\; =\; $ \vec{a}_2 \cdot \vec{a}_2^\ast $\; =\; $ 2\pi \\
\vec{a}_1 \cdot \vec{a}_2^\ast $\; =\; $ \vec{a}_2 \cdot \vec{a}_1^\ast $\; =\; $ 0 \\
\end{array}
\end{equation}
```
The Bragg condition can be written as a matrix equation in a more compact form
``` math
\begin{equation}
{\cal A}^\dagger \cdot {\cal A}^\ast =
2\pi \left( \begin{array}{cc} 1&0\\0&1 \end{array} \right)
\label{eq_int_Bragg}
\end{equation}
```
by using the real and reciprocal lattice matrices $`\cal A`$ (see above) and $`\cal A^\ast`$.
``` math
\begin{equation}
\cal A^\ast
=  \left(
\begin{array}{c}
\vec{a}_1^\ast \\ \vec{a}_2^\ast
\end{array} \right)
= \left(
\begin{array}{cc}
a_{1x}^\ast & a_{1y}^\ast \\
a_{2x}^\ast & a_{2y}^\ast
\end{array} \right)
\end{equation}
```

From Equation ([eq_int_Bragg]) we get the following expressions for $`\cal{A}^\ast`$ and, hence, for the reciprocal lattice vectors $`\vec{a}_1^\ast`$ and $`\vec{a}_2^\ast`$:
``` math
\begin{eqnarray}
\cal{A}^\ast & = &
2\pi \cdot {\cal A}^{\dagger -1} =
\frac{ 2\pi }{\det (\cal{A})}
\left( \begin{array}{rr} a_{2y} & -a_{2x} \\ -a_{1y} & a_{1x}
       \end{array} \right) \\
       \vec{a_1^\ast} &=& \frac{2\pi }{a_{1x}a_{2y} - a_{1y}a_{2x}}
       \left( \begin{array}{r}
a_{2y} \\ -a_{2x}
\end{array} \right)
\\
\vec{a_2^\ast} &=& \frac{2\pi }{a_{1x}a_{2y} - a_{1y}a_{2x}}
\left( \begin{array}{r}
-a_{1y} \\ a_{1x}
\end{array} \right)
\end{eqnarray}
```

These reciprocal lattice vectors have units of Å$`^{-1}`$ and are also parallel to the surface. In analogy to the real lattice vectors, they define the LEED pattern in $`k`$-space. Each diffraction spot corresponds to the sum of integer multiples of $`\vec{a}_1^\ast`$ and $`\vec{a}_2^\ast`$.

``` math
\begin{equation}
\vec{g}(n_1, n_2) = n_1 \vec{a}_1^\ast + n_2 \vec{a}_2^\ast
\end{equation}
```

The pair of integer numbers $`(n_1,n_2)`$ are used as indices to label the spot. The parallel component of the corresponding wave vector is:

``` math
\begin{equation}
\vec{k}_\parallel(n_1, n_2) = \vec{k}_{\parallel,0} + \vec{g}(n_1, n_2)
\end{equation}
```

where $`\vec{k}_{\parallel,0}`$ is the parallel component of the wave vector, $`\vec{k}_0`$, of the incoming electron beam. The length of $`\vec k`$ is $`|\vec k| = \sqrt{2 m_e E_{kin} / \hbar^2 }`$ for both the incoming and elastically scattered electrons. The vertical component, $`k_z`$, of the electrons back-diffracted into spot $`(n_1,n_2)`$ is defined by energy conservation:

``` math
\begin{equation}
k_z(n_1, n_2) = \sqrt{\frac{2 m_e E_{kin}}{ \hbar^2 } - |\vec{k}_\parallel(n_1, n_2)|^2 }
\end{equation}
```

This equation also limits the number of observable LEED spots since the argument of the square root must be greater than zero. With increasing electron energy the number of LEED spots increases while the polar emission angle with respect to the specular spot (0,0) decreases for each spot. The specular spot does not change its position if the angle of incidence is kept constant. Figure [fig_int_LEED_patt] shows examples of common surface unit cells and the corresponding LEED patterns.

### Superstructures

In many cases the periodicity of the surface lattice is larger than expected for the bulk-truncated surface of a given crystal due to adsorption or reconstruction. This leads to additional (superstructure) spots in the LEED pattern, for which fractional indices are used. There are two nomenclatures commonly used to describe such superstructures, the matrix notation and the Wood notation , which will be described in the following. As far as the input for CLEED is concerned, only the matrix notation is used as this is more general.

#### Matrix Notation

The lattice vectors $`\vec{b}_1`$ and $`\vec{b}_2`$ of such superstructures can be expressed as multiples of the $`(1 \times 1)`$ lattice vectors $`\vec{a}_1`$ and $`\vec{a}_2`$:

``` math
\begin{equation}
\begin{array}{lcl}
\vec{b}_1 $\; =\; $ m_{11} \vec{a}_1  +  m_{12} \vec{a}_2  \\
\vec{b}_2 $\; =\; $ m_{21} \vec{a}_1  +  m_{22} \vec{a}_2  \\
\end{array}
\end{equation}
```

where the numbers $`m_{ij}`$ are the coefficients of the superstructure matrix , which offers a straightforward way of characterising the superstructure.
``` math
\begin{eqnarray}
\cal M &=&
\left( \begin{array}{cc} m_{11} & m_{12} \\ m_{21} & m_{22} \end{array} \right)
\nonumber \\
\left( \begin{array}{c} \vec{b}_1\\ \vec{b}_2 \end{array} \right)
&=& \cal B = \cal M \cdot \cal A
\label{eq_int_M}
\end{eqnarray}
```

The area of the superstructure unit cell, $`A_B`$, is given by
``` math
\begin{equation}
A_B = \left | \vec{b}_1 \times \vec{b}_2 \right | = \det {\cal B} = \det {\cal M} \cdot \det {\cal A}
    = \det {\cal M} \cdot A_A =
\end{equation}
```
where $`A_A = \det \cal{A}`$ is the area of the $`(1 \times 1)`$ unit cell (see above) and the determinante of the superstructure matrix is
``` math
\begin{equation}
\det {\cal M} = (m_{11}m_{22} - m_{12}m_{21})
\end{equation}
```

The positions and indices of the additional LEED spots can be calculated directly from the Bragg condition analogue to Equation ([eq_int_Bragg]) when we use a matrix $`\cal M^\ast`$ which correlates the reciprocal vectors of the superstructure to those of the substrate:
``` math
\begin{equation}
\cal B^\ast = \cal M^\ast \cdot \cal A^\ast
\label{eq_int_Mrec}
\end{equation}
```
Then
``` math
\begin{equation}
{\cal B}^\dagger \cdot {\cal B}^\ast
= ({\cal M A})^\dagger \cdot ({\cal M^\ast A^\ast})
= {\cal A}^\dagger \cdot ({\cal M}^\dagger {\cal M}^\ast ) \cdot {\cal A}^\ast
= 2\pi \left( \begin{array}{cc} 1&0\\0&1 \end{array} \right)
\label{eq_int_BraggB}
\end{equation}
```
Since $`{\cal A}`$ and $`{\cal A}^\ast`$ fulfil the Bragg condition already, we find that the matrix product $`{\cal M}^\dagger {\cal M}^\ast`$ must be unity, and hence
``` math
\begin{equation}
{\cal M}^\ast =
{\cal M} ^{{\dagger} -1} =
\frac{1}{m_{11}m_{22} - m_{12}m_{21}}
\left( \begin{array}{rr}
m_{22} & -m_{21} \\
-m_{12} & m_{11}
\end{array} \right)
\end{equation}
```
Re-substituting this result into Equation ([eq_int_Mrec]) we arrive at the following set of formulae for $`\vec{b}_1^\ast`$ and $`\vec{b}_2^\ast`$:

``` math
\begin{equation}
\begin{array}{lcl}
\vec{b}_1^\ast $\; =\; $ (m_{11}m_{22} - m_{12}m_{21})^{-1} \cdot ( +m_{22} \vec{a}_1^\ast  -  m_{21} \vec{a}_2^\ast )  \\
\vec{b}_1^\ast $\; =\; $ (m_{11}m_{22} - m_{12}m_{21})^{-1} \cdot ( -m_{12} \vec{a}_1^\ast  -  m_{11} \vec{a}_2^\ast )  \\
\end{array}
\end{equation}
```

The fractional indices of the superstructure spots are multiples of the pre-factors of $`\vec{a}_1^\ast`$ and $`\vec{a}_2^\ast`$ in the above equations. The area of the superstructure unit cell, A, in units of the $`(1 \times 1)`$ unit cell area can also easily be calculated from the determinante of the superstructure matrix:

``` math
\begin{equation}
A = (m_{11}m_{22} - m_{12}m_{21})
\end{equation}
```

#### Wood Notation

An alternative, less general notation according to Wood specifies the lengths of the vectors $`\vec{b}_1`$ and $`\vec{b}_2`$ in units of $`\vec{a}_1`$ and $`\vec{a}_2`$ together with the rotation angle $`\alpha`$ between $`\vec{b}_1`$ and $`\vec{a}_1`$ (only specified if non-zero):

``` math
\begin{equation}
p/c \left (  \frac{|b_1|}{|a_1|} \times \frac{|b_2|}{|a_2|}  \right ) R \alpha
\end{equation}
```

”p” indicates a ”primitive” and ”c” a ”centred” surface unit cell. Examples are ”$`p(2 \times 2)`$”, ”$`p(\sqrt{3} \times \sqrt{3}) R 30^\circ`$”, and ”$`p(2 \times 2)`$”. This notation is not applicable in all cases but it is more frequently used than the matrix notation because it is shorter. Figure [fig_int_LEED_patt] shows an example of a $`p(2 \times 1)`$ superstructure with the corresponding matrix and Wood notations. Other examples are shown in Figure [fig_int_LEED_domains]

### Domains

<figure id="fig_int_LEED_domains" data-latex-placement="h">
<img src="docs/Fig_int_Domains_r3xr3.jpg" style="height:30.0%" />
<img src="docs/Fig_int_Domains_c2x4.jpg" style="height:30.0%" />
<img src="docs/Fig_int_Domains_r7xr7.jpg" style="height:30.0%" />
<figcaption> Experimental LEED patterns formed by CO adsorbed on Ni{111} (left) and corresponding real-space unit cells (right): <span class="math inline">$p(\sqrt{3} \times \sqrt{3}) R30^\circ$</span> (top, only one domain) <span class="math inline"><em>c</em>(2 × 4)</span> (middle, 3 rotational domains) and <span class="math inline">$p(\sqrt{7} \times \sqrt{7}) R19^\circ$</span> (bottom, 2 mirror domains) <span class="citation" data-cites="Held98a"></span>. Note that real space diagrams are rotated by about <span class="math inline">30<sup>∘</sup></span> with respect to the crystal orientation of the experiment. </figcaption>
</figure>

When the symmetry of the superstructure lattice is lower than the symmetry of the substrate and the adsorption sites one expects domains of different superstructure orientations to be present at the surface in equal amounts. These orientations are found by applying the missing symmetry operations to the superstructure lattice. As the electron beam of a conventional LEED system is orders of magnitude larger than the typical domain size (100 - 1000 Å), the observed LEED pattern is normally a superpositions of all rotational and/or mirror domains. Hence the LEED pattern usually shows the same symmetry as the substrate, even if the superlattice symmetry is lower. For examples see Figure [fig_int_LEED_domains]. Important exceptions are the superstructures of enantiopure chiral molecules, which themselves do not have any mirror symmetry. For these structures mirror domains are not observed since applying a mirror operation would turn the molecule into the opposite enantiomer. For the same enantiomer applying a mirror operation to the lattice would alter the adsorption site geometry of the molecule.

Each domain is characterised by its own superstructure matrix. They are related by the transformation matrices of the corresponding symmetry operations $`\cal T`$.
``` math
\begin{equation}
{\cal M}_{T} =
{\cal M} \cdot ({\cal A} \cdot {\cal T}^{\dagger} \cdot {\cal A}^{-1} )
\label{eq_int_dom}
\end{equation}
```
The matrices for common symmetry operations, rotation by an angle $`\varphi`$, $`{\cal R}_{\varphi}`$, and reflection with respect to the $`x`$ or $`y`$ axis, $`{\cal S}_{x,y}`$ are given below
``` math
\begin{equation}
{\cal R}_{\varphi} =
\left ( \begin{array}{rr}
\cos \varphi  & - \sin \varphi \\  \sin \varphi & \cos \varphi \end{array}
\right) ;\;
{\cal S}_x =
\left ( \begin{array}{rr} 1 & 0 \\ 0 & -1 \end{array} \right ) ;\;
{\cal S}_y =
\left ( \begin{array}{rr} -1 & 0 \\ 0 & 1 \end{array} \right ) .
\label{eq_int_sym}
\end{equation}
```

## Principles of LEED-IV analysis

Because of the strong contribution of multiple scattering in low energy electron diffraction, a direct conversion of LEED-IV data into structural parameters by Fourier Transformation, similar to the Patterson function used in X-ray diffraction, is not reliable. Several attempts have been made of devising direct methods for the conversion of LEED data into atom coordinates but so far they do not deliver the same precision and as conventional ”trial and error” analysis.

The latter is based on full dynamical (multiple) scattering calculations of the LEED intensities over a range of energies for trial surface geometries. These $`I(E)`$ or IV curves are compared with the experimental data, whereby the agreement is quantified through a reliability or $`R`$ factor. The trial geometries are optimised until a minimum of the $`R`$ factor is found.

The following chapters will introduce the basic theories behind the intensity calculations, the definition of $`R`$ factors and the principles of some common search algorithms used in the CLEED package.

# Data Acquisition

## Instrumentation

A sufficiently large experimental data set of good quality is as essential for a successful structure analysis, as are the careful execution of IV calculations and the structural search. We therefore discuss briefly the components of the LEED experiment and the key steps in the data treatment. Especially when Pendry’s R factor is used in the structure optimisation it is important to reduce experimental noise as much as possible.

### LEED Optics

<figure id="fig_LEED_opt" data-latex-placement="h">
<img src="docs/Fig_LEED_Schem.jpg" style="width:80.0%" />
<figcaption> Schematic diagram of a conventional LEED system. </figcaption>
</figure>

The standard modern LEED optics is of the "rear view" type which is schematically shown in Figure [fig_LEED_opt]. The electron beam is emitted from a hot filament, accelerated by the potential $`V_0`$ and collimated by an array of electrostatic lenses (not shown in Figure [fig_LEED_opt]). The electron gun is behind the hemispherical fluorescent screen, which is made of glass and has a hole in the middle for the electrons to pass through. Typically, the electron beam has a current of around 1 $`\mu`$A and illuminates an area of 1 mm$`^2`$. The surface is in the centre of the hemisphere so that all back-diffracted electrons travel towards the LEED screen on radial trajectories. Before the electrons hit the screen they have to pass a retarding field energy analyser. It typically consists of four (or three) hemispherical grids concentric with the screen, each containing a central hole, through which the electron gun is inserted. The first grid (nearest to the sample) is connected to earth ground, as is the sample, to provide an essentially field-free region between the sample and the first grid. This minimises an undesirable electrostatic deflection of diffracted electrons. A suitable negative potential $`-(V_0 - \Delta V)`$ is applied to the second and third (only second) grid, the so-called suppressor grids, to allow only electrons to be transmitted to the fluorescent screen whose energy is within e?V of the elastically scattered electrons. The fourth (third) grid is usually grounded in order to reduce field penetration of the screen voltage to the suppressor grids. The screen voltage is of the order of 5-7 kV; it provides the electrons with enough energy to make the diffraction pattern visible on the fluorescent screen. Older designs have opaque screens, which only allow the observation from behind the sample, where the sample holder obstructs the view. All modern instruments have transparent fluorescent screens allowing the pattern to be observed through a view-port from behind the screen. With this design only the electron gun assembly (diameter $`< 15`$ mm) limits the view slightly. So-called MCP-LEED systems with position sensitive single or double-stack ”multi channel plate” (MCP) electron multipliers between the grids and the fluorescent screen have become commercially available in recent years for applications that require low incident electron currents, either to avoid beam damage (see below) or charging of insulating samples. These systems can be operated with electron currents as low as 1 nA. Several groups have also developed position sensitive electron counting devices, which use primary currents in the pA range . Such extremely low currents are, however, more imposed by the maximum count rates of the detectors, than buy the need of keeping the electron-induced damage within reasonable limits. Typical LEED systems have diameters of around 140 mm so that they fit into a 150 mm (ID) flange, which leads to a screen diameter of around 100 mm, but also smaller versions for 100 mm flanges are available.

The standard way of recording LEED patterns is by using a video camera, either at standard video rate or with variable slow scan CCD systems. Still images are recorded for different kinetic energies, typically between 50 and $`> 300`$ eV in steps of 1 eV. The spot intensities are then extracted from these images using suitable image processing software (see Chapter [chap_exp_Proc]).

The $`k`$-space resolution of LEED systems is characterised by the so-called transfer width. Features in the spot profiles of the LEED pattern (spot splitting, spot broadening) correspond to typical distances on the surface (step separation, island size). The largest resolvable length is the transfer width . Typical values are around 150 Å for conventional rear view LEED systems. Specialised spot profile analysis (SPA) LEED systems can achieve transfer widths up to 1000 Å .

### Quantum efficiency of LEED detectors

It is worthwhile to estimate the limits of video LEED systems in general as opposed to position-sensitive counting detectors. Let us assume the same retarding field optics and therefore the same transmission in front of both detector types. About every second electron arriving at the channel plates of a counting device will be detected . In a video LEED system the electron will hit the fluorescent screen and excite, on average, a number of photons $`N_p`$ which is proportional to a small power ($`1.5 \leq n \leq 2`$) of the applied acceleration voltage $`V_{acc}`$ :
``` math
\begin{equation}
N_p = \frac{e_{phos}}{V_0^n} \cdot V_{acc}^n
\end{equation}
```
$`N_p`$ is about 50 for a commerical LEED system with P-22G fluorescent coating and an acceleration voltage of 7 kV . For a MCP LEED system $`N_p`$ has to be multiplied with the MCP gain factor, typically $`g_{MCP} = 10^4`$ for a single and $`10^7`$ for a double MCP . Since the photons are emitted isotropically the fraction arriving at the camera lens and eventually recorded by the CCD array is determined by the acceptance angle $`\Delta \Omega`$ of the camera
``` math
\begin{equation}
\Delta \Omega = \frac{\pi r_l^2}{d_l^2}
\end{equation}
```
($`r_l`$ is the effective radius of the lens aperture, $`d_l`$ the distance between screen and lens). For a typical system, as described in , with $`r_l`$ = f/2.4 = 10 mm and $`d_l`$ = 260 mm we calculate $`\Delta \Omega = 4.6\times 10^{-3}`$ sr). The overall optical transmission factor $`t_{opt}`$ of screen, view ports, filters and camera lens can be estimated to about 50%. A photon eventually arriving at the CCD will excite an electron with a probability $`e_{ccd}`$ of about 30% for the green light emitted from phosphors of standard LEED systems. Hence, the probability for an electron arriving at the LEED screen to cause an electron–hole pair in the CCD to be excited is:
``` math
\begin{equation}
e_{tot} = e_{ccd} \cdot t_{opt} \cdot \frac{\Delta \Omega}{4\pi} \cdot N_p \cdot g_{MCP}
\end{equation}
```
or
``` math
\begin{equation}
e_{tot} = \frac{e_{ccd} \cdot t_{opt} \cdot e_{phos}}{4} \cdot
          \frac{r_l^2}{d_l^2} \cdot \frac{V_{acc}^n}{V_0^n} \cdot g_{MCP}
\label{eq_tot2}
\end{equation}
```
We then get a total efficiency of $`e_{tot} = g_{MCP} \cdot 3 \times 10^{-3}`$ for the system described in this paper. This marks the limit achievable using commercial LEED optics and camera components.

Equation ([eq_tot2]) shows that a single MCP system already leads to a total efficiency greater than one, i.e. every electron hitting the screen will be detected. But also with conventional LEED systems there is still great potential for increasing $`e_{tot}`$ by improving the optical system. Especially $`\Delta \Omega`$ can be dramatically increased by using a lens with larger aperture, which requires a larger focal length for practical reasons and therefore a larger CCD chip. A commercially available photographic lens with 55 mm focal length and a relative aperture of 1.2 could already increase $`e_{tot}`$ by a factor of about 20, at the expense of an about three times larger CCD. An even more dramatic improvement can be achieved by coupling the LEED screen and the CCD directly through fiber-optical system which is, however, not feasible with a conventional LEED system . A higher acceleration voltage $`V_{acc}`$ will also increase $`e_{tot}`$ about quadratically. 20 kV instead of the 7 kV used conventional LEED systems would therefore lead to a total efficiency about eight times higher. Together with a camera lens, as described before, a conversion factor of $`e_{tot} = 20 \cdot 8 \cdot 3 \times 10^{-3} = 0.5`$ could be achieved.

The above considerations show that commercial video LEED detectors, if operated under the right conditions, are effectively counting devices. This enables us to estimate the ideal operating conditions in terms of primary current, $`I_p`$. From LEED IV calculations we know that relative backscattered intensities for energies above 50 eV are between $`I_{rel} = 10^{-6}`$ and $`10^{-3}`$. For narrow diffraction spots this intensity is distributed over the area of the beam diameter, typically 1 mm$`^2`$. When a LEED screen with diameter 100 mm is imaged onto a 1 Mpixel CCD this corresponds to around $`10 \times 10 = 100`$ pixels per mm$`^2`$. hence, the current per pixel is

``` math
\begin{equation}
I_{pix} = I_{p} \cdot I_{rel} \cdot e_{tot} \cdot 10^{-2}
\label{eq_Ipix}
\end{equation}
```

For several reasons (surface contamination, settling times of LEED electronics, etc) the optimum data acquisition time is around 1-2 s per energy point. In order to reduce noise it is best to accumulate charge on a slow-scan CCD and read out only once per energy as opposed to reading out at video frequency and averaging several images in the data acquisition computer. The former procedure is limited by the fact that the maximum charge per pixel is limited to around $`10^5`$ e. Hence, the maximum primary current is limited by the fact that this charge should not be exceeded when irradiated for 1 s by photons from an intense spot, i.e. $`I_{pix, max} = 10^5`$ e s$`^{-1}`$.
``` math
\begin{equation}
I_{p,max} = \frac{I_{pix, max} }{I_{rel,max} \cdot e_{tot}} \cdot 10^{2} = \frac{10^{10}}{e_{tot}} \mbox{ e s}^{-1}
\label{eq_Ipmax}
\end{equation}
```
This primary current would lead to a signal of $`10^2`$ e s$`^{-1}`$ per pixel or $`10^4`$ e s$`^{-1}`$ per spot for weak intensities, i.e. enough for an acceptable signal-to-noise ratio. Thus, depending on $`e_{tot}`$ the optimum primary current is between around 1 and 500 nA.

If charging of the sample is a problem lower currents and correspondingly longer acquisition times may be required.

### Beam damage

Although LEED is considered a non-destructive method, electron-induced damage to the sample can be a problem, in particular when organic molecules are studied. As a rule of thumb, a conventional LEED system with a beam current of 1 $`\mu`$A and a spot cross section of 1 mm$`^2`$ leads to an exposure of the order of 1 electron per molecule and second. If I-V curves are recorded over an energy range of 250 eV in steps of 1 eV with a data collection time of 1 s per step this leads to a total exposure of 250 electrons per molecule. There are cases in the literature, however, where an exposure of the order of one electron per molecule can cause significant damage to a molecular layer . Beam damage can be reduced by using lower currents in the 1-100 nA range in connection with fast video-LEED systems and by scanning the sample during data acquisition. Short data acquisition times also minimise contamination effects from the residual gas. Such experimental set-ups enable detailed LEED studies of systems which are highly sensitive to beam damage and/or to contamination. Examples are ordered structures of water , benzene or rare gases on metal surfaces .

The minimum beam currents required for video-LEED systems with MCP-intensified LEED screens (about 1 nA for substrate spots only, about 10-100 nA for weak superstructures) are significantly higher than beam currents used in counting devices (pA range). Even considering the much longer data acquisition times necessary for these systems the total electron doses are still smaller. However, such extreme conditions are imposed by the maximum count rates of their detectors, and are not needed in order to keep electron–induced damage within reasonable limits. One can test the stability of a particular adsorbate system against electron bombardment by simply recording the intensity of a superstructure beam versus electron dose (i.e. time of irradiation). As a typical example, the intensity of the ($`\frac{1}{7},\frac{2}{7}`$) beam of the $`(\sqrt{7} \times \sqrt{7})`$ R19$`^\circ`$ superstructure of benzene on Ru{0001} at an energy of 45 eV and a beam current of 130 nA stayed roughly constant with irradiation time for 10 minutes and decayed slowly afterwards; after 50 minutes of irradiation it had decreased by about 30%.

The most beam-sensitive overlayer we have analysed up to now is $`(\sqrt{3} \times \sqrt{3})`$ R30$`^\circ`$ H$`_2`$O/Ru(001) , which is already degraded significantly by about 60 s of irradiation at 200 nA (intensity change of about 10%). In this case, data acquisition had to be split up into segments of 50 eV (steps of 1 eV, 1 s per energy point) in order to achieve reliable results. However, the IV curves of the majority of surface structures investigated in our laboratory so far could be collected in one run without compromising the quality of the data, since a full energy scan (30 eV - 400 eV) can be performed within a total exposure time which is well below the onset of significant sample degradation (e.g. 370 s for benzene on Ru(001) with 1 s irradiation per energy point and steps of 1 eV).

## Data Preparation

The raw data set is a series of still images recorded over the desired kinetic energy range, typically in steps of 1 eV. The spot intensities are then extracted from these images and plotted as function of the kinetic energy for each diffraction spot. The final step before the actual LEED-IV analysis consists of smoothing the raw IV curves and selecting the energy range to be entered into the structure determination.

### Extraction of LEED-IV curves

We have developed a computer program which allows to extract the IV curves of all beams simultaneously from the series of LEED images taken during the energy scan. The only input required are the superstructure matrix and the image positions and indices of at least three spots at an arbitrary energy. The fact that all data processing is done off-line on a workstation allows detailed analysis of the LEED pattern. Dense LEED patterns or large differences in the intensities of neighboring spots can be handled reliably within the data analysis program. Having the entire image information available at a relatively high spatial resolution (about 0.2$`^\circ`$/ per pixel for a $`512 \times 512`$ pixel image) also allows spot profile analysis to some extent , which is mainly limited by the Moiré pattern of the grids.

#### Determine Spot Position and Intensity

Within the program the spot positions at any given energy are then calculated from the reciprocal lattice vectors and the position of the (0,0) beam determined from the spot positions at the preceding energy and scaled according to the ratio between last and current electron momentum, i.e. the square root of the energy ratio.

In order to eliminate minor distortions introduced by the electron and light optics the final spot positions are determined by searches first for the maximum intensity in the vicinity of the calculated position and then for the local first moment of the intensity in the vicinity of the maximum intensity, defined as:
``` math
\begin{equation}
\vec{x}_{i,\vec{g}} =
\sum _{|\vec{x}_i - \vec{x}_{i,max}| < r}
       I(\vec{x}_i)\cdot \vec{x}_i
\left /
\sum _{|\vec{x}_i - \vec{x}_{i,max}| < r} I(\vec{x}_i)
\right .
\end{equation}
```
($`\vec{x}_{i}`$ are image coordinates, $`r`$ is the radius for the local search). This latter search is particularly important if spot profiles are split or otherwise distorted by long range order effects.

The intensity is then integrated over an elliptical or circular area - depending on the spot shape - centered at this final spot position $`\vec{x}_{\vec{g},i}`$. The half axes (or radii) of the integration areas, as well as the radii for the maximum searches, are kept constant with respect to the reciprocal lattice vectors (typically 5% to 10% of their lengths) and therefore scaled with 1/$`\sqrt{E}`$. This way the integration and search areas always cover the same section of $`k_\parallel`$ when the energy is varied. A lower limit for the radius can be defined to account for the finite spot size of the incoming electron beam. Because the background of stray light and diffusely scattered electrons cannot be expected to be uniform, the local background is determined separately for each spot from a ring of about equal area surrounding the integration area and subtracted from the intensity inside. Before background subtraction the signal–to–background ratio is evaluated; if it is below a certain value (typically 10%), the maximum searches for this particular spot are repeated within a smaller range to avoid being misled by tails of closeby intense spots. The calculated spot position is used if an intensity maximum cannot be found within this range. These provisions have proved to be absolutely necessary for the extraction of IV curves from very dense LEED patterns such as $`p(\sqrt{7} \times \sqrt{7})`$ or $`p(2\sqrt{3} \times 2\sqrt{3})`$ superstructures at higher energies (up to 400 eV). Also, this ensures that beams which are temporarily invisible due to total extinction can be followed very easily since other spots are used to determine the reciprocal lattice. Equally, the collection of diffuse IV curves at any $`k_\parallel`$–position (defined by an appropriate superstructure matrix) is straightforward since the reciprocal lattice is defined through the integral order beams (in this case, of course, the background subtraction must be omitted).

Finally the background–corrected integrated intensity is normalized with respect to the beam current at the relevant energy and the exposure time. The positions of the most intense spots are used to recalibrate the reciprocal lattice vectors and the position of the (0,0) beam (which is back–reflected into the electron gun for data collection at normal incidence) in order to calculate the spot positions for the next energy point. While the IV curves are extracted, the current image with the two integration areas for each spot (signal and background) is displayed on the terminal screen of the workstation in order to enable the user to monitor the data processing. Depending on the number of spots and the number of energy points (i.e. images), extracting a set of IV curves from the stored images takes between 15 min. and several hours on a DEC 5000 workstation.

#### Off-normal Incidence

See MKIV manual

#### Image Distortion

##### Mapping on a Flat Screen (MCP-LEED)

See MKIV Manual.

##### Finite Distance of Camera (Curved Screen)

In addition, a geometrical correction has to be included which accounts for the distortion introduced by the finite distance between LEED screen and camera lens:
``` math
\begin{equation}
\vec{x}_{i,real} =
\frac{M \cdot \vec{k_\parallel}}{1 + \frac{r_s}{d_l}(1-|k_\parallel|/|k|)}
\label{eq_xreal}
\end{equation}
```
($`\vec{x}_{i,real}`$ are the recorded image coordinates, $`M`$ is the magnification factor of the optical system, $`r_s`$ = 66 mm is the screen radius, $`d_l`$ = 260 mm the distance between screen and camera lens; see also Figure [fig_LEED_screen]).

<figure id="fig_LEED_screen">
<embed src="docs/Fig_LEED_Screen.eps" style="width:40.0%" />
<figcaption>Distortion of the LEED pattern due to the finite distance <span class="math inline"><em>d</em><sub><em>l</em></sub></span> between camera lens and screen: For a very large distance <span class="math inline"><em>d</em><sub><em>l</em></sub></span> between screen and camera lens the curvature of the screen <span class="math inline"><em>r</em><sub><em>s</em></sub></span> can be ignored and the spot coordinates in the image plane <span class="math inline"><em>x⃗</em><sub><em>i</em>, <em>i</em><em>d</em><em>e</em><em>a</em><em>l</em></sub></span> are an undistorted map of the reciprocal lattice (i.e. uniformly proportional to the parallel components of the outgoing wave vectors). For finite distances the screen curvature leads to image coordinates <span class="math inline"><em>x⃗</em><sub><em>i</em>, <em>r</em><em>e</em><em>a</em><em>l</em></sub></span> that are shifted towards the (0,0) beam with respect to <span class="math inline"><em>x⃗</em><sub><em>i</em>, <em>i</em><em>d</em><em>e</em><em>a</em><em>l</em></sub></span>. The shift is bigger with increasing <span class="math inline"><em>k</em><sub>∥</sub>/<em>k</em></span> (see equation ([eq_xreal])). </figcaption>
</figure>

### Smoothing of IV curves, Selection of energy range

# Principles of Data Analysis: Theory and Algorithms

## LEED-IV Calculation - Program structure and Algorithms

### Overview over Program structure

In the **CLEED** code fully dynamical scattering theory has been implemented along the lines of the algorithms described by Pendry (layer doubling for multiple scattering between successive layers of atoms ) and Van Hove/Tong (combined space method for layers with more than one atom per unit cell ). The extensive use of dynamical memory allocation - an intrinsic feature of the C programming language - allows the memory requirements, even for large surface unit cells, to be kept small. It also allows a very flexible input format which does not impose any restrictions on the number of bulk layers and overlayers nor on the number of atoms therein (other than the physically available memory). The input from the optimization program to the LEED program is simply a set of atomic coordinates from which the program creates its own set of Bravais and/or composite layers. For these layers scattering matrices are calculated and used to evaluate the amplitudes for multiple scattering between the layers.

The main emphasis in the first development stage of the code had been put on reducing the time needed for the giant matrix inversion which is part of the combined space method. In the second stage the symmetries among beams were exploited in order to gain additional reductions in computing times. It turns out, however that the savings in computer time are not big enough to justify the additional restrictions that have to be imposed on the search and the possibility of non-normal incidence angles. Therefore, this path was abandoned and to date only one the non–symmetrised version `cleed_nsym` is supported. It can treat most cases of surface geometries, provided, there is at least one bulk inter–layer distance greater than `MIN_DIST` = 1.0 Å.

`cleed_nsym` does not include any linear approximations such as ”Tensor LEED” , which are in general incompatible with global search strategies such as simulated annealing or the genetic algorithm, and can lead to errors even in downhill–oriented optimizations when the search path enters regimes of the parameter space where the approximation is inaccurate. It is, however, planned to implement this as an option for locally confined searches in the future version.

The input information required from the user is kept as little as possible; the program creates most of its parameters from the geometry input provided in the input files. The distribution of atoms over different layers is performed by the program based on the atom coordinates and the minimum possible vertical distance (`MIN_DIST`) between two layers. The list of beams (basis set for plane wave representation) is created on the basis of the final energy and the minimum distance between the layers.

### Stacking of layers

### Multiple scattering in Bravais layers

### Multiple scattering in Composite layers

### Description of scatterers

#### Muffin-tin potential, phase shifts

#### Isotropic thermal vibrations

#### Anisotropic thermal vibrations - non-spherical scattering potentials

#### Molecular T-matrix

### Domain Mixing

### General Considerations

## R-factors

### General Considerations

LEED Structure determination is based on comparing experimental IV curves with those calculated for a trial surface geometry. Ultimately, the structure determination is the process of searching for the model geometry with the best agreement. In order to be able to automatise the search procedure, one needs to quantify the level of agreement between two curves, which is done by ”reliability” or $`R`$ factors Due to several approximations made in the model calculations, the agreement will never be perfect. In particular, the absolute intensity of the back-scattered electron beams is not well reproduced by the LEED calculations. Different $`R`$ factors emphasize certain properties of the IV curves more than others. In order to make the right choice for a particular set of data, it is important to be aware of these criteria, which will be discussed in the following sections.

#### "Distance between two mathematical functions"

The experimental and all calculated IV curves of a spot $`\vec{g}`$ (or simply $`g`$) can be described in mathematical terms as a set of positive limited functions of $`E`$ defined in the overlap interval $`[E_{i},E_{f}]`$.
``` math
\begin{equation}
{\cal I}  = \{I_{g,x}(E): E \in [E_{i},E_{f}]\}\;,
\end{equation}
```
the index $`x`$ is either $`e`$ (experimental) or $`t`$ (theoretical). Mathematical topology enables calculating ”distances” between such functions by defining a metric, which essentially serve the purpose of an $`R`$ factor. The general properties of a metric are as follows: :
``` math
\begin{equation}
\begin{array}{clr}
(0)   &  R(I_{g,e},I_{g,t}) \geq 0 & \\
(i)   &  R(I_{g,e},I_{g,t}) = 0  \Leftrightarrow I_{g,e}=I_{g,t} & \\
(ii)  &  R(I_{g,e},I_{g,t}) = R(I_{g,t},I_{g,e})  & \\
(iii) &  R(I_{g,e},I_{g,t}) \leq R(I_{g,e},I_{g,t'}) + R(I_{g,t'},I_{g,t})
      &  \mbox{(triangle inequality)} \\
\end{array}
\label{eq_r_def_metric}
\end{equation}
```
Note that only condition (iii), the ”triangle inequality”, enables us to talk about large and small distances between functions. In general, the scaling of experimental and theoretical IV curves will be very different. Therefore, only relative intensities, i.e. the curve shape, should be compared . Hence conditions (i) and (ii) of ([eq_r_def_metric]) must be modified:
``` math
\begin{equation}
\begin{array}{cl}
(ia)  & R(I_{g,e},I_{g,t}) = 0 \Leftrightarrow I_{g,e}  = a \cdot I_{g,t}  \\
(iia)  & R(I_{g,e}, I_{g,t}) = R(a\cdot I_{g,t},b\cdot I_{g,e})\\
      & (a, b  \in \cal{R}) \\
\end{array}
\label{eq_r_def_rfa}
\end{equation}
```
((iia) is a consequence of (ia) and (iii)). Consequently, before calculating the actual distance the IV curves need to go through a suitable transformation in order to eliminate differences due to differences in scaling. In mathematical terms, this is a projection onto a sub-space of reduced dimensionality.

#### Normalisation of R factors

The final modification of the general definition of a metric to generate a suitable $`R`$ factor is introducing normalisation. In order to put meaning to the absolute values of $`R`$ factors, we add condition (ib) as follows:
``` math
\begin{equation}
\begin{array}{cl}
(ib) & R(I_{g,e},I_{g,t}) = 1 \Leftrightarrow I_{g,e} \mbox{ and } I_{g,t}
       \mbox{ are completely un-correlated}\\
\end{array}
%\label{eq_r_def_rfb}
\end{equation}
```

(ib) can be achieved dividing an un-normalised metric by the value found for two un-correlated functions Through this final condition $`R`$ factors will mostly assume values between 0 and 1. $`R`$ factors greater than 1 correspond to ”anti-correlated” curves where maxima of one curve coincide mostly with minima of the other. The above conditions still allow a multitude of different R factor definitions but if they fulfil condition (ib) they should lead to similar values.

The $`R`$ factor for a certain model geometry is the weighted average of all available pairs of IV curves:
``` math
\begin{equation}
R_{avg} =  \frac{ \sum_{g}(E_{f}-E_{i})_{g} \cdot R[I_{g,e}(E),I_{g,t}(E)]}
                { \sum_{g}(E_{f}-E_{i})_{g} }
\label{eq_r_avg_rf}
\end{equation}
```

There is a number of mathematical transformations that can be used to construct an $`R`$ factor metric. Some of them are used as standard $`R`$ factors for LEED, which will be discussed in the following sections. We also introduce a few new $`R`$ factors which could be used to emphasize additional curve properties. Not all of the standard $`R`$ factors comply with the triangle inequality (iii) and are therefore actually concave functions of a true metric. As a consequence, they lead to sharper contrast, which can speed up the search for minima.

### $`R_1`$ and $`R_2`$

### Pendry’s R factor $`R_P`$

Logarithmic derivatives have been suggested by Pendry as another way of mapping IV curves onto a subspace which is independent of scaling :

``` math
\begin{equation}
I(E) \rightarrow  L[I(E)] = \frac{\frac{\partial I_{x}(E)}{\partial E}}
                   {I(E)}
\label{eq_r_rp_proj}
\end{equation}
```

Since differentiation is a linear operation, the logarithmic derivative of a function $`I(E)`$ is the same as that of $`a \cdot I(E)`$. But also features of different intensities within the same curve have the same weight. For instance, a model IV curve which is composed of Lorentzian functions of equal width $`2 V_{0i}`$ but different intensities $`a_{n}`$
``` math
\begin{equation}
I(E) = \sum_{n} \frac {a_{n}}{(E-E_{n})^{2} + V_{0i}^2},
\label{eq_r_rp_modell}
\end{equation}
```
has a logarithmic derivative which is independent of the individual intensities $`a_{n}`$, provided the Lorentzians are well separated. The logarithmic derivative is:
``` math
\begin{equation}
L[I(E)] = \sum_{n} \frac {-2(E-E_{n})}{(E-E_{n})^{2} + V_{0i}^2}
\end{equation}
```

One drawback of the logarithmic derivative ([eq_r_rp_proj]) is that it has singularities at the points where I(E) passes through zero. Such singularities cannot be excluded completely, therefore Pendry suggested the so-called $`Y`$-function whose denominator $`(1 + L^2 \cdot V_{0i}^2)`$ avoids such singularities.
``` math
\begin{equation}
Y(E) = \frac {L} {1 + L^2 \cdot V_{0i}^2}
\label{eq_r_rp_Y}
\end{equation}
```
$`V_{0i}`$ should be equal to the imaginary part of the optical potential used in the IV calculations, typically around 4 eV. For Lorentzians of the form ([eq_r_rp_modell]) Y(E) oscillates between $`\pm \frac{1}{2|V_{0i}|}`$. These extrema occur at the positions $`E_{n} \pm |V_{0i}|`$ and are independent of the intensity of the IV curve. More generally, the maxima of Y(E) mark the rising slopes and the minima the falling slopes of the corresponding IV curve. $`Y(E)`$ passes through zero at the points where the IV curve has a maximum or minimum.

<figure id="fig_RFAC_IY" data-latex-placement="h">
<embed src="docs/FIG_RU001_IV.eps" style="width:60.0%" />
<embed src="docs/FIG_RU001_YR.eps" style="width:66.0%" />
<figcaption> Top: Experimental (solid) and calculated (dotted) IV curve of the (1, 0) spot of the clean Ru{0001} surface. Bottom: corresponding <span class="math inline"><em>Y</em></span> functions. </figcaption>
</figure>

Figure [fig_RFAC_IY] shows experimental and calculated IV curves of the (1, 0) spot of the clean Ru{0001} surface (top) together with the corresponding $`Y`$ functions (bottom: $`V_{0i} = -6\mbox{ eV} \Rightarrow \max |Y(E)| = 0.83\mbox{ eV}^{-1}`$).

Pendry’s $`R`$ factor, $`R_{P}`$, is defined as the normalised mean square deviation of the two functions $`Y_{g,t}(E)`$:
``` math
\begin{equation}
R_{P} = \frac{ \int _{E_i}^{E_f} (Y_{g,t} - Y_{g,e} ) ^2 }
             { \int _{E_i}^{E_f} (Y_{g,t}^2 + Y_{g,e}^2) }
\label{eq_r_rp_def}
\end{equation}
```

For any pair of uncorrelated functions $`Y_{g,t}`$ and $`Y_{g,e}`$ the mean of their product is equal to the product of the means:
``` math
\begin{equation}
\langle Y_{g,t} \cdot Y_{g,e} \rangle =
\langle Y_{g,t} \rangle \cdot \langle Y_{g,e} \rangle \; ,
\label{eq_r_rp_unkorr}
\end{equation}
```
Since the $`Y`$ function oscillates around zero, their mean will be zero
``` math
\begin{equation}
\langle Y_{g,t} \rangle = \langle Y_{g,e} \rangle = 0 \; ,
\end{equation}
```
and the denominator in Equation ([eq_r_rp_def]) is equal to the mean square deviation of two un-correlated $`Y`$ functions:
``` math
\begin{equation}
\langle (Y_{g,t} - Y_{g,e} ) ^2 \rangle _{unkorr} =
\langle Y_{g,t}^2 \rangle + \langle Y_{g,e}^2 \rangle.
\end{equation}
```
Hence $`R_P = 1`$ for uncorrelated IV curves, which fulfils the normalisation condition (ib).

Due to the properties of the $`Y`$ function $`R_P`$ is particularly sensitive to differences in the peak positions between the two curves. Also, small features, e.g. shoulders and small peaks, have a much greater effect on $`R_P`$ than on $`R_1`$ or $`R_2`$. On the other hand, $`R_P`$ is less sensitive to difference in the absolute peak heights. Therefore non-linearities in the experiment and slow variations in the primary current have little effect. In the example of Figure [fig_RFAC_IY] there is a minimum and peak in the calculated curve between 60 eV and 80 eV, which makes a large contribution to the value of $`R_P`$, as can be seen from the large difference in the $`Y`$ function in this energy range. The relatively large differences in the peak heights at 90 eV and around 245 eV, however, have very little effect. As the $`Y`$ function is largely insensitive to the height of a peak it is essential to filter out small artificial structures due to experimental noise as much as possible. Otherwise they would be treated with similar weight as true peaks. The strong oscillations of the experimental $`Y`$ function in the region between 108 eV and 125 eV, where the intensity is essentially zero, demonstrates the effect of experimental noise very well. In order to avoid such artifacts, both the experimental and theoretical data are convoluted with a Lorentzian function of width $`V_{0i}`$ . Due to the attenuation of elastically scattered electrons features in the IV curves are expected to be at least $`2 V_{0i}`$ wide. Therefore, this procedure does not lead to a significant loss of information. It is applied in `crfac` also when calculating $`R_1`$ and $`R_2`$. Despite such smoothing, it is still not recommended to use $`R_P`$ for the comparison of IV curves where the intensity is close to zero over a large energy range. In such cases $`R_1`$ or $`R_2`$ usually perform better .

## Structure Optimisation

### Overview over Program structure

LEED IV structural analysis is based on the ”trial and error” approach. IV curves and the R factor are calculated for a start geometry. This geometry is varied until a minimum in the R factor is found. A number of established optimisation algorithms have been proposed for LEED, including the Simplex method , Powell’s algorithm, the Marquard’s algorithm, Simmulated annealing and genetic algorithms . Several of these are implemented in the CLEED package. `csearch` is the master program performing the structural optimisation. It iteratively calls sub-programs for the calculation of IV curves and the R factor. Following one of several optimisation algorithms chosen by the user, it creates new trial geometries based on the history of the search and and terminates the optimisation when a R factor minimum is found.

#### Data structure

The search parameters are not necessarily identical with the Cartesian coordinates of the atoms. Iternally the program uses a matrix $`{\cal P}`$ which relates the set of search parameters $`p_1`$, $`p_2`$, $`p_3`$, …, $`p_{Np}`$ to the Cartesian coordinates $`x_1`$, $`y_1`$, $`z_1`$, …, $`x_{Na}`$, $`y_{Na}`$, $`z_{Na}`$ of the $`N_a`$ atoms included in the search and the angles of incidence, $`\theta, \varphi`$ (if the angles are optimised):

``` math
\begin{equation}
\left( \begin{array}{c}
x_1 \\ y_1 \\ z_1 \\ \vdots  \\ z_{Na} \\ x_{ang} \\ y_{ang}
\end{array} \right)
= 
\left( \begin{array}{c}
x_{1,0} \\ y_{1,0} \\ z_{1,0} \\ \vdots  \\ z_{Na,0} \\ x_{ang,0} \\ y_{ang,0}
\end{array} \right)
+ {\cal P} \cdot
\left( \begin{array}{c}
p_1 \\ p_2 \\ p_3 \\ \vdots \\ p_{Np} \\
\end{array} \right)
\end{equation}
```
$`N_a`$ and $`N_p`$, respectively, are the numbers of atoms and search parameters. $`x_{a,0}`$, $`y_{a,0}`$, $`z_{a,0}`$, $`\cdots`$ are the coordinates of the start geometry. The $`3N_a`$ coordinate displacements of the atoms from their start values are defined by the $`j`$th row of $`{\cal P}`$:
``` math
\begin{equation}
\Delta x_i = \sum_{j=1}^{N_P} {\cal P}_{i,j} \cdot p_j
\end{equation}
```
the matrix can either be supplied directly via the input file (`spp: - x/y/z` parameter specifier, see Section [sec_CSEARCH_in]) or be generated by the program using the symmetry of the surface structure.

If atoms are related through a rotation axis or a mirror plane perpendicular to the surface, their $`z`$ coordinates have to be identical. In addition, their lateral coordinates are linked such that there are only three search parameters linked to the entire group of symmetry-related atoms. Through the direct input of $`{\cal P}`$ more complicated symmetry operations can be realised, such as glide lines. Certain atoms can be optimised in $`x,y,z`$, while others are only optimised in $`z`$ or not at all; or groups of atoms can be kept as a rigid unit only optimising their centre of mass position.

Angles of incidence $`\theta, \varphi`$ are not optimised directly but through
``` math
\begin{equation}
\begin{array}{rcl}
x_{ang,0} &=& \cos \varphi \sin \theta \\
y_{ang,0} &=& \sin \varphi \sin \theta \\
\end{array}
\end{equation}
```
and
``` math
\begin{equation}
\begin{array}{rcl}
\varphi &=& \arctan \left ( \frac{x_{ang,0} + \Delta x_{ang}}{y_{ang,0} + \Delta y_{ang}} \right ) \\
\theta  &=& \arcsin \left ( \sqrt{(x_{ang,0} + \Delta x_{ang})^2 + (y_{ang,0} + \Delta y_{ang})^2} \right ) \\
\end{array}
\end{equation}
```
This way $`\theta = 0^\circ`$ can be used as start value, which would otherwise cause a degeneracy in $`\varphi`$. The surface geometry can be optimised using data sets recorded at different angles of incidence. Normally, the angles will have to be optimised together with the structural parameters. Each pair of incidence angles increases the number of search parameters by two.

### Search Algorithms

Several optimisation algorithms have been proposed for LEED, including the Simplex method , Powell’s algorithm, the Marquard’s algorithm, Simmulated annealing and genetic algorithms . The latter two are based on random generation of trial geometries and should, in principle, find the global R factor minimum if the parameters of the search are chosen appropriately. The other algorithms are guided by the shape of the R factor hypersurface and follow the direction of steepest decent, which means that they can be trapped in local R factor minima. The choice of the right search algorithm is, in general, a trade-off between large numbers of calculations necessary for random searches vs the risk of being trapped in a local minimum.

#### Simplex Method

#### Simulated Annealing

# Manual for the CLEED Package

## CLEED

### General Description

#### Syntax

The general calling syntax of the LEED program is:

    cleed_nsym -i <parameter file> -b <bulk parameter file> -o <results file>
               -r <storage file>   -w <storage file>

The first argument, `-i <parameter file>`, specifying the parameter input file is the only mandatory. The file contains all the geometric and non–geometric parameters needed for the LEED calculations. A sample file is shown in Table [tab_prg_inp_par]. Alternatively the input can be split into two files, the parameter file and a bulk parameter file. The latter file (`-b <bulk parameter file>`) contains all the parameters which are not varied during the optimisation. Consequently, the search program has to produce only the parameter file containing the optimised atom positions of the overlayer in each iteration step of an automated search.

`-o <results file>`  
specifies the file name where the beam intensities will be written to. If this argument is not used, the results will be written to the file `leed.res`. Any existing files with the same name as the results file will be overwritten.

`-r <storage file>` and `-w <storage file>`  
specify files where the bulk reflexion matrices can be read from (`-r`) or written to (`-w`). This feature can be used in order to save computer time. It is, however, only recomended for very time consuming bulk matrix calculations, i.e. when the bulk unit cells contain more than one atom.

#### UNIX–Environment

The program has been developed in a UNIX environment and uses a few UNIX specific features. The most important is the use of environment variables:

<div class="description">

`CLEED_PHASE`: directory where the files are stored which contain the energy dependent atomic phase shifts.



These variables have to be set using the export or setenv UNIX commands, respectively before the program is called for the first time.

### Input files

#### Input of non-geometrical Parameters

Table [tab_leed_inp_par] shows an example for a single parameter input file for a $`p(\sqrt 3 \times \sqrt 3)`$ H$`_2`$O / Ru(0001) overlayer structure.



    # Sample input file for
    c: Ru(0001) + p(r3xr3)-2O
    # bulk unit cell parameters
    a1:       1.3525       -2.3426    0.0000
    a2:       1.3525        2.3426    0.0000
    a3:       0.0000        0.0000   -4.2800
    # superstructure matrix
    m1:  2. 1.
    m2: -1. 1.
    #
    # Overlayer Atoms:
    # O atoms on hcp site
    po: O_H2O                     2.7050  0.0000  4.3100  dr3 0.05 0.05 0.05
    po: O_H2O                     1.3525  2.3426  4.1900  dr3 0.05 0.05 0.05
    # Ru atoms of the first layer
    po: /home/CLEED/PHASE/Ru.phs  0.0000  0.0000  2.0400  dr1 0.0707
    po: /home/CLEED/PHASE/Ru.phs  2.7050  0.0000  2.0500  dr1 0.0707
    po: /home/CLEED/PHASE/Ru.phs  1.3525  2.3426  2.1100  dr1 0.0707
    #
    # Bulk layer Atoms:
    pb: /home/CLEED/PHASE/Ru.phs  0.0000 -1.5617  0.0000  dmt 400. 101. 200.
    pb: /home/CLEED/PHASE/Ru.phs  0.0000  0.0000 -2.1400  dmt 400. 101. 200.
    # NON-GEOMETRIC PARAMETERS:
    # optical potentials
    vr:   -13.00
    vi:     4.50
    # energies
    ei: 32.
    ef: 260.1
    es: 4.
    # angles of incidence
    it: 0.
    ip: 0.
    # epsilon, lmax
    ep: 1.e-4
    lm: 8

Each parameter or set of parameters, respectively, is specified by two letters and a colon at the beginning of a line. There can only be one parameter specifier per line. Comments can either be indicated by a leading ”`#`” or by ”`c:`”. In the first case the comment is ignored by the LEED program, in the latter case it is stored by the program and written to the output file.

The meaning and the syntax of the parameter specifiers is (`n`: integer number, `f`: floating point number `c`: character):

<div class="description">

`a1:`, `a2:`, `a3:` `f f f`  
Lattice vectors of the three–dimensional bulk unit cell (x, y, z in Å). a1 and a2 are parallel to the surface plane, i.e. they define the two–dimensional $`(1\times 1)`$ unit cell. a3 must contain a component perpendicular to the surface. For the symmetrised program version a3 must not have any parallel components. In some cases (e.g. fcc lattice) a non–primitive bulk unit cell must be used in order to fulfill this condition.

`m1:`, `m2:` `f f`  
Super structure matrix defining the relationship between the superstructure lattice vectors b1 and b2 and the $`(1\times 1)`$ lattice vectors a1 and a2:
``` math
\begin{array}{lcl}
\vec b_1 &=& m1(1) \cdot \vec a_1 + m1(2) \cdot \vec a_2\\
\vec b_2 &=& m2(1) \cdot \vec a_1 + m2(2) \cdot \vec a_2\\
\end{array}
```

`vr:` `f`  
Real part of optical potential (in eV).

`vi:` `f`  
Imaginary part of optical potential (in eV).

`po:` `<phasestring> f f f ccc f {f f} `  
Atom parameters in overlayer (super structure unit cell):  
`<phasestring>`: Name of phase shift file. It can either be specified by the full path (starting with a leading slash ’/’, as in `/home/CLEED/PHASE/Ru.phs`) or by the file name body (no leading ’/’, e.g. `O_H2O`) which will then be expanded into `CLEED_PHASE/<phasestring>.phs` by using the environment variable `CLEED_PHASE`. If `CLEED_PHASE = /home/CLEED/PHASE` in the current example, the input of `O_H2O` is equivalent to `/home/CLEED/PHASE/O_H2O.phs`.  
`f f f`: atom position x, y, z (in Å);  
`ccc f {f f}`: three character specifier for input of vibrational displacements:

<div class="description">

`dmt f f f`: input of Debye temperature ($`\Theta_D`$ in K), temperature ($`T`$ in K), mass (in amu). From these values the isotropic radial root mean square displacement
``` math
\sqrt{<\Delta r^2>} = \frac{9}{2  m \; k_B \Theta_D } \cdot
                        \sqrt{\frac{1}{16} + \frac{T^2}{\Theta_D^2}}
```
will be calculated and used inside the program.

`dr1 f`: input of radial root mean square displacement in Å.
``` math
\sqrt{<\Delta r^2>} = \sqrt{<\Delta x^2> + <\Delta y^2> + <\Delta z^2> }
```

`dr3 f f f`: separate input of root mean square displacements along the coordinates $`\sqrt{<\Delta x^2}>`$, $`<\sqrt{\Delta y^2}>`$, and $`<\sqrt{\Delta z^2}>`$ in Å. From these values the average radial root mean square displacement $`\sqrt{<\Delta r^2>}`$ (see above) will be calculated and used as isotropic displacement inside the program, thus the inputs ”`dr1 0.0866`” and ”`dr3 0.05 0.05 0.05`” are equivalent.

`nd3 f f f`: separate input of root mean square displacements along the coordinates $`\sqrt{<\Delta x^2}>`$, $`<\sqrt{\Delta y^2}>`$, and $`<\sqrt{\Delta z^2}>`$ in Å (Only non–symmetric version). The displacements are used in order to set up a non–diagonal atomic $`t`$–matrix which is used in the program.



`pb:` `<phasestring> f f f ccc f {f f} `  
Atom parameters in bulk layers ($`(1\times 1)`$ unit cell); for syntax see `po:`.

`ei:` `f`  
First energy in energy loop for which intensities are calculated (in eV).

`ef:` `f`  
Last energy in loop (in eV).

`es:` `f`  
Energy step in loop (in eV).

`it:` `f`  
Polar angle of incidence with respect to surface normal ($`0^\circ \leq \vartheta < 90^\circ`$ ).

`ip:` `f`  
Azimuthal angle of incidence with respect to x–axis ($`0^\circ \leq \varphi < 360^\circ`$).

`ep:` `f`  
Epsilon: criterion for convergence of lattice sums and layer doubling.

`lm:` `n`  
Maximum angular momentum quantum number ($`l_{max}`$).



The number of overlayer atoms specified by `po:` must be exactly the same as the number of atoms within one two–dimensional overlayer unit cell given by the overlayer matrix. However, they can lie in different unit cells. The bulk atoms specified by `pb:` must be exactly those within the topmost three–dimensional bulk unit cell specified by `a1`, `a2`, and `a3`. All overlayer atoms must have larger $`z`$ coordinates than the top–most bulk atom. The program will produce unreliable results if the vertical distance between the top–most bulk atom and the bottom–most overlayer atom is shorter than `MIN_DIST` = 1.0 Å (Note, the value of `MIN_DIST` can be changed by editing the `leed_def.h` header file and re-compiling the program).

#### Separate bulk and overlayer parameter input

If the parameter input is split into two files, all parameters except for the overlayer atom parameters are read from the bulk parameter input file specified by the `-b` option. The overlayer atom parameters (`po:`) are read from the parameter input file (option `-i`). The two files corresponding to the above example of $`p(\sqrt 3 \times \sqrt 3)`$ (H$`_2`$)O / Ru(0001) are shown in Tables [tab_leed_inp_bul] and [tab_leed_inp_ovl].



    # sample bulk geometry input file
    c: Ru(0001) + p(r3xr3)
    # bulk unit cell parameters
    a1:       1.3525       -2.3426    0.0000
    a2:       1.3525        2.3426    0.0000
    a3:       0.0000        0.0000   -4.2800
    # superstructure matrix
    m1:  2. 1.
    m2: -1. 1.
    #
    # bulk layers:
    pb: /home/CLEED/PHASE/Ru.phs  0.0000 -1.5617  0.0000  dmt 400. 101. 200.
    pb: /home/CLEED/PHASE/Ru.phs  0.0000  0.0000 -2.1400  dmt 400. 101. 200.
    #
    # NON-GEOMETRIC PARAMETERS:
    # optical potentials
    vr:   -13.00
    vi:     4.50
    # energies
    ei: 32.
    ef: 260.1
    es: 4.
    # angles of incidence
    it: 0.
    ip: 0.
    # epsilon, lmax
    ep: 1.e-4
    lm: 8



    # sample overlayer geometry input file
    c: Ru(0001) + p(r3xr3)-2O
    # Overlayers:
    # O atoms on hcp site
    po: O_H2O                     2.7050  0.0000  4.3100  dr3 0.05 0.05 0.05
    po: O_H2O                     1.3525  2.3426  4.1900  dr3 0.05 0.05 0.05
    # Ru atoms of the first layer
    po: /home/CLEED/PHASE/Ru.phs  0.0000  0.0000  2.0400  dr1 0.0707
    po: /home/CLEED/PHASE/Ru.phs  2.7050  0.0000  2.0500  dr1 0.0707
    po: /home/CLEED/PHASE/Ru.phs  1.3525  2.3426  2.1100  dr1 0.0707

This way of input is used within the automated search where only the optimised `po:` parameters are supplied by the search program whereas the unchanged parameters are read from a separate bulk file provided by the user.

#### Surfaces with small inter-layer distances ($`< 1`$Å)

#### Phase shifts input

The energy dependent atomic phase shifts $`\delta_l(E)`$ are read from the files specified in the atom parameter sets (`po:` and `pb:`). These files have the format as shown in Table [tab_leed_inp_phase]



    40 9 neng, lmax (Ru, m= 101)
     0.5000
     2.2748-0.3054-0.5527 0.0180 0.0007 0.0000 0.0000 0.0000 0.0000 0.0000
     1.0000
     1.8705-0.5465-0.6358 0.1689 0.0115 0.0009 0.0001 0.0000 0.0000 0.0000
     1.5000
     1.5527-0.7538-0.7190 0.5972 0.0497 0.0056 0.0005 0.0000 0.0000 0.0000
     2.0000
     1.2909-0.9420-0.7769 1.1747 0.1270 0.0190 0.0024 0.0003 0.0000 0.0000

                        .....................................

    19.0000
    -1.2839-3.0896-1.8887 2.6390 1.6041 1.0384 0.6948 0.4697 0.3216 0.2165
    19.5000
    -1.3173-3.1180-1.9057 2.6399 1.6129 1.0490 0.7059 0.4795 0.3300 0.2243
    20.0000
    -1.3499-3.1458-1.9223 2.6407 1.6212 1.0593 0.7166 0.4891 0.3382 0.2319

The first line specifies the number of energy values and maximum angular momentum quantum number $`l_{max}`$ for which phase shifts are stored in the file. The rest of the first line is treated as comment and ignored by the program. The rest of the file consists of pairs of lines containing the energy value (in Hartree) in the first and the phase shifts $`\delta_l`$ (in the order $`l`$ = 0, 1, 2, ..., $`l_{max}`$) for this energy in the second line. The phase shifts are either separated by blanks or by ’–’ which is interpreted as sign of the following number. Inside the program the phase shifts are interpolated for each energy of the I–V curve or extrapolated, respectively, when the I–V curve surpasses the energy range of the input file. It is, however not extrapolated to energies below the lowest energy value in the input list.

### Output file

#### Results - IV curves

Table [tab_leed_out_res] shows a sample output file for a $`p(\sqrt 3 \times \sqrt 3)`$ H$`_2`$O / Ru{0001} overlayer structure.

    # ####################################### #
    #            output from CLEED            #
    # ####################################### #
    #vn 0.20 (version GH/30.08.97)
    #ts Thu Jul  2 07:18:21 1998
    #
    #en 58 32.000000 260.100000 4.000000
    #bn 109
    #bi 0 0.000000 0.000000 0
    #bi 1 -1.000000 0.000000 0
    #bi 2 -1.000000 1.000000 0
    #bi 3 0.000000 -1.000000 0
    #bi 4 0.000000 1.000000 0
    #bi 5 1.000000 -1.000000 0
    #bi 6 1.000000 0.000000 0
    #bi 7 -2.000000 1.000000 0
    #bi 8 -1.000000 -1.000000 0
    #bi 9 -1.000000 2.000000 0
    #bi 10 1.000000 -2.000000 0

    .......................

    #bi 104 -3.333333 2.666667 2
    #bi 105 0.666667 -3.333333 2
    #bi 106 0.666667 2.666667 2
    #bi 107 2.666667 -3.333333 2
    #bi 108 2.666667 0.666667 2
    32.00 1.758462e-02 1.557583e-03 6.515479e-04 6.515850e-04 1.557595e-03
    1.557651e-03 6.515733e-04 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
    ............ 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00

                           ..............

    260.00 1.303590e-02 7.845590e-05 5.554994e-03 5.555435e-03 7.845676e-05
    7.851732e-05 5.555437e-03 8.115542e-04 7.950326e-04 7.949291e-04 8.117350e-04
    ............ 1.905631e-05 1.905916e-05 1.286778e-05 1.286476e-05 1.906064e-05

The header of the file consits of lines beginning with ’`#`’. They contain information needed e.g. by the R factor program:

<div class="description">

`vn` program version.

`ts` creation date and time.

`en` number of energy points, first energy, last energg, energy step.

`bn` total number of beams in output files.

`bi` beam counter, 1st beam index, 2nd beam index, beam set (for each beam).



The header is followed by lines containing the energy (in eV) and the relative beam intensities in the same order as the beam list separated by blanks.

### Special Functions and Libraries

#### Matrix Library

##### Matrix type

In contrast to the conventional definition of matrices in C as arrays of pointers to arrays of numbers (integer, float, etc.), in our programs matrices are structures `struct mat_str` (defined in `"mat.h"`, see table [tab_mat_str]). This has the advantage that all relevant information such as matrix dimensions (`rows,cols`) and the type of matrix elements (`num_type`) can be exchanged between functions together with the matrix elements themselves through a single variable of type `mat`, which is a pointer to this structure (`typedef struct mat_str * mat` in `"mat.h"`).

    struct mat_str     /* real or complex matrix */
    {
     int mag_no;       /* magic number */
     int mat_type;     /* type of matrix (square, etc.) */
     int num_type;     /* type of matrix elements (read, complex) */
     int rows;         /* 1st dimension of matrix (number of rows) */
     int cols;         /* 2nd dimension of matrix (number of columns) */
     real *rel;        /* pointer to real parts of matrix elements */
     real *iel;        /* pointer to imag. parts of matrix elements */
    };

In order to make the matrix type independent of the type of the matrix elements, complex matrices are effectively stored as two matrices, one containing the real part of the elements (`*rel`) and one containing the imaginary part (`*iel`). For real matrices `iel` is a null pointer, in addition, the structure element `num_type` is an integer number representing the type of the matrix elements.

The structure element `mat_type` is another integer number representing the matrix type (rectangular, square, diagonal) which decides upon the storage scheme of the matrix elements and various other things. For example, the matrix inversion routine checks first of all, if the matrix to be inverted is square at all, and sends an error message if not.

The storage scheme for rectangular and square matrices is:
``` math
\begin{equation}
\Re[M(m,n)] \rightarrow \mathtt{*(rel + (m-1)*rows + n)}
\end{equation}
```
``` math
\begin{equation}
\Im[M(m,n)] \rightarrow \mathtt{*(iel + (m-1)*rows + n)}
\end{equation}
```
diagonal matrices are stored as:
``` math
\begin{equation}
\Re[M(m,m)] \rightarrow \mathtt{*(rel + m)}
\end{equation}
```
``` math
\begin{equation}
\Im[M(m,m)] \rightarrow \mathtt{*(iel + m)}
\end{equation}
```
Note that the first element in the arrays `rel` and `iel` is never used, therefore these arrays have the length $`[(n\cdot m + 1) \cdot \mathtt{sizeof(real)}]`$ for a rectangular or square $`n\times m`$ matrix and $`[(n + 1) \cdot \mathtt{sizeof(real)}]`$ for a diagonal matrix.

##### Matrix operations

The general format for matrix functions which have a matrix as results is:



` (mat) result = matfunction( (mat) storage, operand_1, operand_2, ...) `



The first parameter of the function call (`storage`) is usually a matrix where the result will be written to (type `mat`, i.e. a pointer to a matrix structure). If the dimensions or matrix type of this matrix do not match with those of the result, the necessary adjustments in memory will be done automatically. In particular, if `storage` is the null–pointer (`NULL`), a new matrix will be created. In any case, the return value `result` is a pointer to this updated matrix which is not necessarily identical with `storage`. The types of the other parameters (operands) depend on the purpose of the particular matrix function. They can be either matrices or single numbers.

##### Basic Matrix Operations

The LEED-IV calculations spend most of the time ($`> 90`$%) multiplying or inverting large matrices with dimensions often greater than 1000.

**Matrix multiplication** is performed through the function



` mat matmul( mat Mr, mat M1, mat M2 )`



where `M1` and `M2` and the operands and `Mr` is the result. The function uses a ”textbook style” algorithm:
``` math
\begin{equation}
\mathtt{Mr}(m,n) = \sum _i \mathtt{M1}(m,i) \cdot \mathtt{M2}(i,n)
\end{equation}
```
with no optimisation other than discriminating between multiplication of real/complex and diagonal/non-diagonal.

 

**Matrix Inversion** is performed by the function



` mat matinv( mat A_1, mat A)`



where `A` is the input matrix to be inverted and `A_1` the inverse (result). The function discriminates between real and complex and then uses ”LU decomposition” to perform the inversion (functions `c_ludcmp` and `c_luinv` in `matclu.c`, which are based on the algorithms described in ”Numerical Receipes” ). Matrix inversion is a key operation in, both layer doubling and the calculation of scattering matrices for composite layers (see Chapters [chap_the LEED] and [sec_CLEED_ms]). Several optimisations are used but they all ultimately call `matinv`.

 

In a recent addition by Michael Fink (Univ. of Innsbruck) the core functionalities of `matmul` and `matinv` have been replaced by the corresponding **BLAS** (`*gemm`) and **LAPACK** (`*getrf and *getri`) functions. This speeds up the calculation very significantly, by a factor of more than 2 (depending on the type of structure) when a single processor is used and up to an order of magnitude if the local BLAS/LAPACK library can access multiple processors (e.g. cygwin).

#####  Matrix allocation functions

The matrix functions described in this section are mainly called by other ”high level” matrix functions in order to check, allocate, free, etc. memory space used to store a matrix. They are listed in table [tab_mat_bas].

<div id="tab_mat_bas">

|                                    |                            |             |
|:-----------------------------------|:---------------------------|:------------|
| function name                      | description                | source file |
| `mat matalloc(mat, int, int, int)` | allocate memory for matrix | matalloc.c  |
| `int matcheck(mat)`                | check validity of pointer  | matcheck.c  |
| `mat matcop(mat, mat)`             | copy matrices              | matcop.c    |
| `int matfree(mat)`                 | free memory                | matfree.c   |

Basic matrix functions



##### Display and Control Functions for Matrices

The matrix functions described in this section are not thought to be used in the ”production state” of a LEED program since they produce very large outputs. During the development of LEED (or other electron scattering) programs it is however useful, to display the contents of matrices in the right format. Two functions are available for this purpose. `int matshow(mat)` displays the real and imaginary part of complex matrix elements while `int matshowabs(mat)` diplays only the modulus of complex matrix elements and therefore reduces the amount of output. For example the matrix:
``` math
\left (
\begin{array}{cc} 1 & 1+i \\ 1-i & i \\ \end{array}
\right )
```
would be displyed by `matshow` as:

    (2 rows) x (2 columns):
    (1.00, 0.00)   (1.00, 1.00)
    (1.00,-1.00)   (0.00, 1.00)

while `matshowabs` would produce the following output:

    *** abs ***  (2 rows) x (2 columns):
    1.00   1.41
    1.41   1.00

Note that both functions also display the matrix dimensions which is particularly useful, when the rows are too long to be diplayed in one line of the screen.

#### Clebsh–Gordan Coefficients

|                                                                |
|:---------------------------------------------------------------|
| `(qmcgc.c:)`                                                   |
| `int mk_cg_coef(int l_max)`                                    |
| `double cg(int l1, int m1, int l2, int m2, int l3, int m3)`    |
| `double blm(int l1, int m1, int l2, int m2, int l3, int m3)`   |
| `double gaunt(int l1, int m1, int l2, int m2, int l3, int m3)` |

All necessary Clebsh–Gordan coefficients (or Gaunt’s coefficients) up to $`l_{max}`$ are calculated once by the function `mk_cg_coef`. The coefficients are stored in a list and can be extracted from there by using the functions `cg`, `blm`, or `gaunt`. These functions read from the same list but have return different values:

<div class="description">

`double cg(int l1, int m1, int l2, int m2, int l3, int m3)`  
returns the value of the integral described in Slater’s book:
``` math
(-1)^{m_2}
\int Y_{l_1 m_1}(\Omega)  Y_{l_2 m_2}^\ast(\Omega)  Y_{l_3 m_3}(\Omega) d\Omega
```
Allowed (non–zero) values:  
$`0 \leq l_1     \leq 2 \cdot l_{max}, \;\; 0 \leq l_{2,3} \leq l_{max}`$;  
$`m_1 =  m_2 + m_3 \;\; \Leftrightarrow \;\; m_1 -  m_2 - m3 = 0`$,  
$`|l_2 - l_3| \leq l_1 \leq l_2 + l_3`$;  
Transformations:  
$`cg(l_1, m_1, l_2, m_2, l_3, m_3) =
 cg(l_1, m_1, l_3, m_3, l_2, m_2) =
 cg(l_2, m_2, l_1, m_1, l_3,-m_3)`$

`double blm(int l1, int m1, int l2, int m2, int l3, int m3)`  
returns the value of the integral used in Pendry’s book:
``` math
\int Y_{l_1 m_1}(\Omega)  Y_{l_2 m_2}(\Omega)  Y_{l_3 m_3}(\Omega) d\Omega
```
Allowed (non–zero) values:  
$`0 \leq l_2     \leq 2 \cdot l_{max}, \;\; 0 \leq l_{1,3} \leq l_{max}`$;  
$`m_1 +  m_2 + m_3 = 0`$,  
$`|l_1 - l_3| \leq l_2 \leq l_1 + l_3`$;  
Transformations:  
$`blm(l_1, m_1, l_2, m_2, l_3, m_3) =
 cg(l_2, -m_2, l_1, m_1, l_3, m_3)`$

`double gaunt(int l1, int m1, int l2, int m2, int l3, int m3)`  
returns Gaunt’s integral:
``` math
(-1)^{m_3}
\int Y_{l_1 m_1}(\Omega)  Y_{l_2 m_2}(\Omega)  Y_{l_3 m_3}^\ast(\Omega) d\Omega
```
Allowed (non–zero) values:  
$`0 \leq l_2     \leq 2 \cdot l_{max}, \;\; 0 \leq l_{1,3} \leq l_{max}`$;  
$`m_3 = m_1 + m_2 \;\; \Leftrightarrow \;\; m_1 + m_2 - m_3 = 0`$,  
$`|l_2 - l_3| \leq l_1 \leq l_2 + l_3`$;  
Transformations:  
$`gaunt(l_1, m_1, l_2, m_2, l_3, m_3) =
 cg(l_2, m_2, l_1, -m_1, l_3,m_3)`$



The memory requirements of the list created by `mk_cg_coef` for a given value of $`l_{max}`$ are:
``` math
(2 \cdot l_{max} + 1) (2 \cdot l_{max} + 2)/2 \cdot
(l_{max} + 1)^2 (l_{max}/2 + 1) \cdot
\mbox{\texttt{ sizeof(double)}}
```



| $`l_{max}`$ | memory (bytes) |
|------------:|---------------:|
|           6 |        142,688 |
|           8 |        495,720 |
|          10 |      1,341,648 |
|          12 |      3,075,800 |
|          14 |      6,264,000 |
|          16 |     11,673,288 |



#### Spherical Harmonics and Hankle Functions

##### Spherical Harmonics

|                                                                 |
|:----------------------------------------------------------------|
| `(qmylm.c:)`                                                    |
| `mat r_ylm( mat Ylm, real x, real phi, int l_max )`             |
| `mat c_ylm( mat Ylm, real z_r, real z_i, real phi, int l_max )` |
| `(lmsymat.c:)`                                                  |
| `(lmsypy.c:)`                                                   |

##### Spherical Hankle Functions

|                                                         |
|:--------------------------------------------------------|
| `(qmhank.c:)`                                           |
| `mat r_hank1 ( mat Hl, real x, int l_max )`             |
| `mat c_hank1 ( mat Hl, real z_r, real z_i, int l_max )` |

### Multiple Scattering Functions

#### Multiple scattering Matrix: Single Bravais Layer

|  |  |
|:---|:---|
| `(lmsbravl.c:)` |  |
| `int ms_bravl` | `( mat `$`\ast`$`p_Tpp, mat `$`\ast`$`p_Rpm,` |
|  | ` struct var_str `$`\ast`$` v_par, struct layer_str `$`\ast`$` layer,` |
|  | ` struct beam_str `$`\ast`$` beams) ` |
| `(lmsbravlnd.c:)` |  |
| `int ms_bravl_nd` | `( mat `$`\ast`$`p_Tpp, mat `$`\ast`$`p_Tmm mat `$`\ast`$`p_Rpm, mat `$`\ast`$`p_Rmp,` |
|  | ` struct var_str `$`\ast`$` v_par, struct layer_str `$`\ast`$` layer,` |
|  | ` struct beam_str `$`\ast`$` beams) ` |
| `(lmslsumii.c:)` |  |
| `mat ms_lsum_ii` | `( mat Llm, real k_r, real k_i, real `$`\ast`$`k_in, real `$`\ast`$`a,` |
|  | ` int l_max, real epsilon ) ` |
| `(lmstmatii.c:)` |  |
| `mat ms_tmat_ii` | `( mat Tii, mat Llm, mat Tl, int l_max) ` |
| `(lmstmatndii.c:)` |  |
| `mat ms_tmat_nd_ii` | `( mat Tii, mat Llm, mat Tlm_in, int l_max) ` |
| `(lmsymat.c:)` |  |
| `mat ms_ymat` | `( mat Ymat, int l_max,` |
|  | ` struct beam_str `$`\ast`$`beams, int n_beams) ` |
| `(lmsypy.c:)` |  |
| `mat ms_yp_yxp` | `( mat Yxmat, mat Ymat) ` |
| `mat ms_yp_yxm` | `( mat Yxmat, mat Ymat) ` |
| `mat ms_yp_ym` | `( mat Ymmat, mat Ypmat) ` |

The function `ms_bravl` calculates the multiple scattering matrix for a single Bravais layer of scatterers. According to equation (??) in Ref. . This equation can be written in Matrix notation as:

``` math
\begin{equation}
 M = \frac{1}{\kappa \cdot A} \;\;
     Y^{out\;\pm} \left \{ (1 - \tau \cdot G)^{-1} \tau \right \} Y^{in\;\pm}
\end{equation}
```

The matrices $`Y^{in}`$ and $`Y^{out}`$ perform the transformation from the plane wave representation (used between the layers) into the spherical wave representation (used inside the layers). They are defined as:
``` math
\begin{equation}
 Y^{in\;\pm}_{lm,\vec k} = Y^\ast_{lm}(\vec k^{\pm})
\end{equation}
```
computed by the functions `ms_ymat` and `ms_yp_ym` and
``` math
\begin{equation}
 Y^{out\;\pm}_{\vec k,lm} = \frac{8i}{k_z} \; Y_{lm}(\vec k^{\pm})
\end{equation}
```
computed by the functions `ms_yp_yxp` and `ms_yp_yxm`.

The expression embraced by the $`Y`$ matrices describes the multiple scattering process inside the layer. $`G`$ is the propagator matrix
``` math
\begin{equation}
G_{l_1 m_1, l_2 m_2} = \sum _{l_3 m_3}
  i^{l_1 - l_2} (-1)^{(l_2 - l_1 - l_3) / 2 - m_1)} \cdot
  C_{l_1 m_1, l_3 m_3, l_2 m_2}       \cdot
  L_{l_3 m_3}
\end{equation}
```
with $`C_{l_1 m_1, l_3 m_3, l_2 m_2}`$ being a Clebsh–Gordan–coefficient and $`L_{l_3 m_3}`$ the lattice sum:
``` math
\begin{equation}
L_{l_3 m_3} = (-1)^{m_3} 4 \pi \cdot Y_{l_3 m_3}(\cos \frac{\pi}{2}, 0) \cdot
  \sum_{\vec R_j}
       \exp [i \vec k_{in} \vec R_j - i \varphi(\vec R_j) ] \cdot
       h^{(1)}_{l_3}(\kappa \cdot |\vec R_j|)
\end{equation}
```
computed in the function `ms_lsum_ii`. $`\vec R_j`$ is the position of the $`j`$th atom inside the layer; $`h^{(1)}_{l_3}`$ the Hankle function of the first kind, and $`\kappa= \sqrt{2 E}`$ the length of the wave vector.

$`\tau`$ is the atomic scattering matrix $`t`$ times $`-\kappa`$. The functions `ms_bravl` and `ms_bravl_sym` assume isotropic scatterers with a diagonal scattering matrix which only depends on $`l`$ and not on the $`m`$ quantum numbers.
``` math
\begin{equation}
 \tau_{l_1 m_1, l_2 m_2} =
 -\kappa \cdot t_{l_1 m_1, l_2 m_2} =
 e^{i \delta_{l_1}} \sin (\delta _{l_1}) \;\;
 \delta _{l_1 l_2} \delta _{m_1 m_2}
\end{equation}
```
The diagonality of $`\tau`$ is explicitly used in the function `ms_tmat_ii` called by `ms_bravl` and `ms_bravl_sym` which computes the matrix expression
``` math
\begin{equation}
 (1 - \tau \cdot G)^{-1} \tau
\end{equation}
```

The function `ms_tmat_nd_ii` and the calling function `ms_bravl_nd` treat the more general case of anisotropic scatterers (e.g. anisotropic thermal vibrations) with non–diagonal scattering matrices which can also depend on $`m`$ (see below).

Before performing any calculations, the function checks whether some quantities needed in the calculation can be reused from the last call of `ms_bravl`. In the nex steps the function calculate the lattice sum (if needed) by calling `ls_lsum_ii`, $`(1 - t \cdot G)^{-1} t`$ by calling `ls_tmat_ii`, and the transformation matrices $`Y^\pm`$ by calling `ms_ymat` first and then `ms_yp_yxp` and `ms_yp_yxm`, respectively. After multiplying the outgoing wavefunctions ($`Y^+`$) with $`i 8 \pi / (|k| A k ^\prime _\perp)`$, the reflection and transmission are put together.

#### Multiple scattering Matrix: Composite Layer

#### Multiple scattering Matrix: Layer Doubling

### CPARA: Parallel Calculation of IV curves on Multi-processor Machines

LEED calculations can be easily distributed over multiple processors by assigning each processor a different energy point. These calculations are completely independent of each other and can therefore be performed by clones of the same program call.

The wrapper program `c_para` reads the input `*.bul` file for `cleed` and generates `n_para` copies,`*.bul.p1`, `*.bul.p2`, etc., each with a different sequence of energies, as illustrated in Figure [fig_cleed_para_eng]

<figure id="fig_cleed_para_eng" data-latex-placement="h">
<img src="docs/Fig_CPARA_Energies.jpg" style="width:80.0%" />
<figcaption> Diagram showing the sequence of energies calculated if <code>c_para</code> is used. <span id="fig_cleed_para_eng" data-label="fig_cleed_para_eng"></span> </figcaption>
</figure>

<figure id="fig_cleed_para_flow" data-latex-placement="h">
<img src="docs/Fig_CPARA_FlowChart.jpg" style="width:50.0%" />
<figcaption> Diagram illustrating the call sequence of <code>cleed</code> from <code>c_para</code>. <span id="fig_cleed_para_flow" data-label="fig_cleed_para_flow"></span> </figcaption>
</figure>

The program then launches `n_para` system calls executing `cleed` simultaneously, one call for each copy of `*.bul.p*` with identical `*.par` input files. After all `cleed` calls have finished, `c_para` reads all `*.res.p*` output files and combines the IV data into one result file, `*.res`, which contains the intensities for all energies in ascending order. This is illustrated in the flow chart of Figure  [fig_cleed_para_flow].

 

The general calling syntax of `c_para` is compatible to `cleed`:

    c_para -i <parameter file> -b <bulk parameter file> -o <results file> -n <n_para>

with the extra `-n` option to provide `n_para`, the number of parallel processes.

Alternatively, `n_para` can also be specified through an additional line:

    np: <n_para>

in the `*.bul` input file. This line will be ignored by `cleed` and produced a warning message. If a value is provided through the `-n` option in the command-line this will overwrite any input from the `*.bul` input file. The default value of `n_para` is 1. This is used if no input is provided or if the input is less than 1. Also, the maximum value of `n_para` is half the total number of energy steps. It will be automatically reduced if the input value is higher.

The program to be executed for the actual calculation of IV curves is specified by the environment variable `CPARA_LEED`.

Tests show that significant savings in overall times can be achieved, albeit not the same factor as `n_para`. Values of 4 or 8 for `n_para` seem to be the optimum for a quad-processor machine running cygwin, however if the cygwin LAPACK library is used for matrix operations the savings are relatively small, as this makes already use of multiple processors.

## CRFAC

### General Description

#### Program description

The R factor program performs the comparison of calculated and experimental IV curves. The agreement is quantified in terms of the R(eliability) factor. It offers the choice between four different R factors that can be used for the search: $`R_1`$, $`R_2`$ , $`R_P`$ , and $`R_B`$ . The output R factor value is the optimum achieved by shifting the energy axes of experimental and theoretical I–V curves with respect to each other. This shift acts as a correction for any non–optimum value of the optical potential, $`V_{0r}`$, in the LEED calculations, which need therefore not be optimized by the search program. In this way, one dimension is eliminated from the search parameter space on which the search program operates, hence reducing the number of LEED calculations to be performed. The assignment of experimental and theoretical I–V curves (or average of curves) for comparison, and of the relative weight that a particular I–V curve has in the overall R factor, is performed by the user prior to the search.

#### Syntax

The general calling syntax of the R factor program is:

    crfac -a <ID flag>  -c <control file>   -h    -o <output file>
          -r <R factor> -s <shift1,shift2,shift3> -t <theoretical file>
          -v <optical potential>  -w <IV output file>

The options `-c <control file>` and `-t <theoretical file>` are mandatory for normal use of the program.

<div class="description">

`-a <ID flag>`  
defines whether only the average R factor is calculated (argument `average`, default) or partial R factors for each subset of IV curves sharing a common ID number (argument `all`). Only the first two characters of the argument are significant.

`-c <control file>`   
specifies the control file which defines the correlation between experimental and theoretical IV curves. A sample control file is shown in Table [tab_rfac_ctr].

`-h`  
causes the program to show a short list of arguments.

`-o <output file>`   
specifies the output file where the R factor values are written to (default: standard output).

`-r <R factor>`   
specifies the R factor to be calculated. Valid arguments are: `r1` ($`R_1`$), `r2` ($`R_2`$), `rb` ($`R_{B1}`$ and $`R_{B2}`$), and `rp` ($`R_P`$, default)

`-s <shift1,shift2,shift3>`   
defines the range (arguments `shift1` and `shift2`) and step width (`shift3`) of energy shifts between experimental and theoretical IV curves.

`-t <theoretical file>`   
specifies the file containing the theoretical IV curves, i.e. the results file from the LEED program.

`-v <optical potential>`   
specifies the value of the optical potential $`V_{0i}`$ used in the evaluation of Pendry’s R–factor (default: 4 eV).

`-w <IV output file>`   
causes the program to write all normalised IV curves as energy/intensity pairs to separate files so that they can be plotted. `<IV output file>` specifies the body of the file names to which the letters `e` (experimental) or `t` (theoretical) and the number of the pair of curves is added.



#### UNIX–Environment

The program has been developed in a UNIX environment and uses a few UNIX specific features. The most important is the use of environment variables:

<div class="description">

`RF_HELP_FILE`: file to be shown, when the `-h` option is chosen.



This variable has to be set using the `export` or `setenv` UNIX commands, respectively before the program is called for the first time.

### CRFAC Input / Output files

#### Control file: assignment of experimental data

The control file is used for the assignment of experimental data files to beams in the output from the IV calculation. Table [tab_rfac_ctr] shows an example R factor control file for a $`p(1 \times 1)`$ Er / Si(111) structure.

    #  Si(111)-Er (1x1)
    #  (control file for R factor program CRFAC)
    #
    #  Theor. indices are listed for 3-fold rot. symmetry output
    #  and a calculation up to 250 eV. Delete if no experimental
    #  data are available for certain spots.
    #
    #
    #  ef=<experimental input file>
    #  ti=<corresponding indices in theoretical input file>
    #  id=<group ID> (does not really matter)
    #  wt=<weight>
    #  : = separator
    #
    #
    ef=Si111_Er_expt_01:ti=(0.00,1.00):id=01:wt=1.
    ef=Si111_Er_expt_02:ti=(1.00,0.00):id=02:wt=1.
    ef=Si111_Er_expt_03:ti=(1.00,1.00)+(2.00,-1.00):id=03:wt=1.
    ef=Si111_Er_expt_04:ti=(0.00,2.00):id=04:wt=1.
    ef=Si111_Er_expt_05:ti=(2.00,0.00):id=05:wt=1.
    ef=Si111_Er_expt_06:ti=(-1.00,3.00)+(1.00,2.00):id=06:wt=1.
    ef=Si111_Er_expt_07:ti=(3.00,-1.00)+(2.00,1.00):id=07:wt=1.
    ef=Si111_Er_expt_08:ti=(0.00,3.00):id=08:wt=1.
    ef=Si111_Er_expt_09:ti=(3.00,0.00):id=09:wt=1.
    ef=Si111_Er_expt_10:ti=(2.00,2.00)+(4.00,-2.00):id=10:wt=1.
    ef=Si111_Er_expt_11:ti=(-1.00,4.00)+(1.00,3.00):id=11:wt=1.
    ef=Si111_Er_expt_12:ti=(4.00,-1.00)+(3.00,1.00):id=12:wt=1.
    ef=Si111_Er_expt_13:ti=(0.00,4.00):id=13:wt=1.
    ef=Si111_Er_expt_14:ti=(4.00,0.00):id=14:wt=1.

Each line defines the correlation of one pair of theoretical and experimental IV curves using the following syntax:

    ef=<expt. file>:ti=<index list>:id=<ID number>:wt=<weight>

<div class="description">

`<expt. file>` contains one experimental IV curve in a two-column format (see below for details)

`<index list>` is a list of indices of beams to be averaged and compared with the experimental IV curve. The two indices of each beam are in round brackets and can be preceeded by a weighting factor followed by $`\ast`$. Different beams are connected through `+`. The general form is:

    {<weight>}*(<index_1>,<index_2>){+{<weight>}*(<ind_1>,<ind_2>)}

`<ID number>` is a integer number attached to this pair of IV curves. When the option `-a all` is specified, the partial R factors will be calculated for all groups of IV curves sharing the same ID numbers.

`<weight>` is a positive real number defining the weight of this pair of IV curves when calculating the average R factor. This weight (default: 1) will be multiplied with the fraction of the total energy range covered by this pair of IV curves.



#### Input of experimental and theoretical IV curves

The input file containing the calculated IV curves (specified by the `-t` option in the command line) is expected to have the same format at the output files from `cleed`.

 

The input files containing the experimental IV curves are specified in the control file (see above and Table [tab_rfac_ctr]). They are ASCII files with two columns, for energies (in eV) and intensities (arb. units), respectively. For an example, see Table [tab_rfac_exp]. Lines starting with `#` at any point within the file are ignored and can be used for comments. Both experimental and theoretical IV curves must have equal intervals between the energy points, however these intervals do not have to be the same for calculated and experimental curves.

    # Sin Fourier Smooth: version 1.1
    # cutoff: 0.500, tailoff: 10.000 => k range: 0.974
    # average over beams:
    #   (-0.50,1.33)
    #   (0.50,1.33)
    #   (-0.50,-1.33)
    #   (0.50,-1.33)
    8.400000e+01 4.440016e+00
    8.600000e+01 3.741519e+00
    8.800000e+01 2.960896e+00
    ......
    2.440000e+02 2.295693e-01
    2.460000e+02 3.104824e-01
    2.480000e+02 3.515357e-01
    # End of Data

#### R factor Output

A typical output line is shown below:

    0.318127 0.185645 3.00 928.50           #  Rp  RR  shift  range

The first value is the actual R factor, followed by Pendry’s $`RR`$ factor defining the statistical error. The third number is the energy shift between the experimental and theoretical IV curves, the last value is the total overlaqp (in eV) between experimental and theoretical IV curves. These values are followed by a comment specifying them.

## CSEARCH

### CSEARCH Input files

#### Parameter Input

##### Automatic assignment of search parameters to coordinates

The name of the parameter input file should have the format `<project>.inp`. Through this the name of the project is defined and `csearch` will also need to be supplied with input files called `<project>.bul` (additional input for `cleed`) and `<project>.ctr` (additional input for `crfac`). Table [tab_search_inp1] shows an example input file for a $`p(1 \times 1)`$ Er / Si(111) overlayer structure.

    # input file for SEARCH
    #  Si(111)-Er (1x1)
    #  atomic positions according to Spence et al. Phys. Rev. B 61 (2000) 5707.
    #
    #
    a1:      1.9162240  -3.3189974   0.0
    a2:      1.9162240   3.3189974   0.0
    #
    m1:  1. 0.
    m2:  0. 1.
    #
    # bulk: (see *.bul)
    # pb: Si   0.0         2.2126649  -0.7822951  dr3 0.05  0.05  0.05
    # pb: Si   0.0         0.0         0.0        dr3 0.05  0.05  0.05
    #
    # overlayer:
    #
    po: Si   0.0         0.0        +2.353      dr3 0.05  0.05  0.05
    po: Si   0.0        -2.2126649  +3.153      dr3 0.05  0.05  0.05
    po: Ru   0.0         0.0        +5.453      dr3 0.05  0.05  0.05
    po: Si   0.0        +2.2126649  +7.233      dr3 0.05  0.05  0.05
    po: Si   0.0        -2.2126649  +8.053      dr3 0.05  0.05  0.05
    #
    # minimum radii:
    #
    rm: Si 1.00
    rm: Ru 1.00
    # z range
    zr: 1.00  9.00
    # sz: 0 (xyz search), 1 (z only)
    sz: 1
    # sr: rotational axis
    sr: 3 0.0 0.0

Each parameter or set of parameters, respectively, is specified by two letters and a colon at the beginning of a line. There can only be one parameter specifier per line. Comments can either be indicated by a leading ”`#`” or by ”`c:`”. In the first case the comment is ignored by the LEED program, in the latter case it is stored by the program and written to the output file. (The meaning and the syntax of the parameter specifiers is `n`: integer number, `f`: floating point number `c`: character):

<div class="description">

`a1:`, `a2:` `f f f`  
Lattice vectors of the two–dimensional $`(1\times 1)`$ unit cell (x, y, z in Å). a1 and a2 have to be parallel to the surface plane. Note `a1:` and `a2:` can also be specified in the `<project>.bul` input file. If they are found in both, the input from `<project>.inp` overwrites the input from `<project>.bul`.

`m1:`, `m2:` `f f`  
Super structure matrix defining the relationship between the superstructure lattice vectors b1 and b2 and the $`(1\times 1)`$ lattice vectors a1 and a2:
``` math
\begin{array}{lcl}
\vec b_1 &=& m1(1) \cdot \vec a_1 + m1(2) \cdot \vec a_2\\
\vec b_2 &=& m2(1) \cdot \vec a_1 + m2(2) \cdot \vec a_2\\
\end{array}
```
Note `m1:` and `m2:` can also be specified in the `<project>.bul` input file. If they are found in both, the input from `<project>.inp` overwrites the input from `<project>.bul`.

`po:` `<phasestring> f f f ccc f {f f} `  
Atom parameters in overlayer (super structure unit cell): same syntax as for the LEED program, see chapter [sec_leed_inp].

`rm:` `<phasestring> f`  
Minimum radius of atoms specified by `<phasestring>`. If the distance between two atoms is smaller than the sum of their miniumum radii, an additional number is added to the actual R factor in order to repell the search from this un–physical geometry.

`zr:` `f f`  
$`z`$–range. The atoms are forced to stay within this range of $`z`$ coordinates. If an atom exceeds this range, an additional number is added to the actual R factor in order to repell the search from this un–physical geometry.

`sz:` `n`  
Variation of vertical parameters (1) or vertical and lateral parameters (0) in agreement with the specified symmetry.

`sr:` `n f f`  
Rotational symmetry: degree (n–fold axis), position of axis (x, y in Å).

`sa:` `n`  
Switch angle search on (1)/off (0).

`it:` `f`  
Polar angle of incidence (start value). Only used if `sa: 1`;will overwrite `it:` parameter in `*.bul` file. Note `it:` can also be specified in the `<project>.bul` input file. If it is found in both, the input from `<project>.inp` overwrites the input from `<project>.bul`.

`ip:` `f`  
Azimuthal angle of incidence (start value). will overwrite `ip:` parameter in `*.bul` file. Note `it:` can also be specified in the `<project>.bul` input file. If it is found in both, the input from `<project>.inp` overwrites the input from `<project>.bul`.



##### Manual assignment of search parameters to coordinates

The above notation only allows imposing relatively simple symmetry constraints and these have to be imposed on all optimised atoms in the same way, for instance `sz:` does not allow to specify a subgroup of atoms to be optimised in $`x,y,z`$ only and another subgroup in $`z`$ only.

A much more flexible way of introducing constraints to the search is shown in Table [tab_search_inp2]. Here two additional types of parameter specifiers are used:

<div class="description">

`spn:` `n`  
defines the number of parameters to be optimised in the search. The relationship between the parameters and the coordinates is defined through the parameters following `spp: - <x/y/z>`. `spn:` must be specified before the list of atoms.

`spp: - <x/y/z>` `f f f .... ` (number of parameters defined by `spn:`  
defines the relationship between the search parameters being optimised and the coordinates of the last atom defined in the input file.



Note, if `spn:` and `spp:` are used, `sz:` `sr` will be ignored.

Within `csearch` the relationship between search parameters and coordinates is defined by a matrix $`{\cal P}`$:
``` math
\begin{equation}
\left( \begin{array}{c}
x_1 \\ y_1 \\ z_1 \\ \vdots \\ x_{Na} \\ y_{Na} \\ z_{Na} \\
\end{array} \right)
= {\cal P} \cdot
\left( \begin{array}{c}
p_1 \\ p_2 \\ p_3 \\ \vdots \\ p_{Np} \\
\end{array} \right)
\end{equation}
```

$`N_a`$ and $`N_P`$ are the numbers of atoms and search parameters, $`p_i`$, as defined by `spn` and the list of atoms in the input file, respectively. The positions of the $`N_a`$ atoms are defined by $`3N_a`$ coordinates; hence the $`[3(i-1)+1]`$th row of $`{\cal P}`$ defines the dependence on the search parameters of the $`x`$ coordinate of the $`i`$th atom through:
``` math
\begin{equation}
x_i = \sum_{j=1}^{N_P} {\cal P}_{i,j} \cdot p_j
\end{equation}
```
this row is defined through the line of numbers following the `spp: - x` parameter specifier after the first atom. Using this approach one can specify a subgroup of atoms to be optimised in $`x,y,z`$ (e.g. the top 4 layers of Cu atoms in Table  [tab_search_inp2]) and another subgroup in $`z`$ only (e.g. the 5th and 6th layer of Cu atoms in Table  [tab_search_inp2]). Note that if `spp: - <x/y/z>` is not defined for an atom, this coordinate will not be optimised, i.e this is equivalent to a line of only zeros after `spp: - <x/y/z>`. The atoms in layer 7 to 16 are not optimised. In addition, more complicated symmetry operations, such as glide lines, can be realised, or groups of atoms can be kept as a rigid unit only optimising their centre of mass position.

    # Cu{531} - geometry input file for the LEED program.
    a1:       5.7116618       0.000        0.0000
    a2:      -2.8557407       3.3789890    0.0000
    #
    m1:  1. 0.
    m2:  0. 1.
    # matrix format
    spn: 14
    # 1st layer - x, y, z optimisation
    po: Cu_met -0.012000 0.028268 -0.142097 dr1 0.110000
    spp: - z       1.0 0.0 0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0
    spp: - y       0.0 0.0 0.0 0.0 0.0 0.0   1.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0
    spp: - x       0.0 0.0 0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0   1.0 0.0 0.0 0.0
    # 2nd layer - x, y, z optimisation
    po: Cu_met 2.262675 1.011614 -0.621688  dr1 0.100000
    spp: - z       0.0 1.0 0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0
    spp: - y       0.0 0.0 0.0 0.0 0.0 0.0   0.0 1.0 0.0 0.0   0.0 0.0 0.0 0.0
    spp: - x       0.0 0.0 0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0   0.0 1.0 0.0 0.0
    # 3rd layer - x, y, z optimisation
    po: Cu_met -1.191446 2.000753 -1.237364 dr1 0.090000
    spp: - z       0.0 0.0 1.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0
    spp: - y       0.0 0.0 0.0 0.0 0.0 0.0   0.0 0.0 1.0 0.0   0.0 0.0 0.0 0.0
    spp: - x       0.0 0.0 0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 1.0 0.0
    # 4th layer - x, y, z optimisation
    po: Cu_met 1.137943 2.921576 -1.756463  dr1 0.080000
    spp: - z       0.0 0.0 0.0 1.0 0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0
    spp: - y       0.0 0.0 0.0 0.0 0.0 0.0   0.0 0.0 0.0 1.0   0.0 0.0 0.0 0.0
    spp: - x       0.0 0.0 0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 1.0
    # 5th layer - only z optimisation
    po: Cu_met 0.535428 0.452284 -2.454718  dr1 0.06500
    spp: - z       0.0 0.0 0.0 0.0 1.0 0.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0
    # 6th layer - only z optimisation
    po: Cu_met 2.911708 1.544562 -3.101747  dr1 0.06500
    spp: - z       0.0 0.0 0.0 0.0 0.0 1.0   0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0
    # 7th layer - no optimisation
    po: Cu_met  -0.5715   2.4133  -3.6640  dr1 0.065

      ......................

    # minimum radii:
    rm: Cu_met 1.00
    # z-range
    zr: -18.00 4.00
    # sa: 1 (angle search on), 0 (angle search off)
    sa: 1
    it: 4.320
    ip: 210.23

#### Bulk parameters

`<project>.bul`. See chapter [sec_CLEED_inp]

#### R factor control file

`<project>.ctr`. See chapter [sec_CRFAC_io]

### CSEARCH Output files

#### Log file

A protocol of the search is stored in the file `<project>.log`, where `<project>` is the name of the parameter input file without extension. The log file starts with a list of atom positions in the start geometry followed by a listing of the parameters and the corresponding R factor values of each search step.

#### LEED files

The input parameter file for the LEED program is called `<project>.par`. In order to call the LEED program successfully, there must also exist a file called `<project>.bul` containing the bulk parameters, i.e. the geometrical and non–geometrical parameters which are not varied within the structure search. This will be copied into `<project>.bsr`, which is the bulk parameter file used during the search. If the angles of incidence are optimised as part of the search, these will differ between `<project>.bul` and `<project>.bsr`; all other parameters remain unchanged. The calculated IV curves will be stored in `<project>.res`. `<project>.out` contains the control output from the LEED program. The input parameter file and the calculated IV curves of the current best fit geometry are copied to `<project>.pmin`, `<project>.bmin` and `<project>.rmin`, respectively.

#### R factor files

The output from the R factor program is written to the file `<project>.dum`. Only the last line will be used, which contains the current R factor value.

#### Simplex Vertex file

If the simplex method is used, the vertices of the simplex are stored in the file `<project>.ver` after each change of the simplex. The file `<project>.vbk` contains a backup of the iteration step before. If the search had to be stopped for any reason (e.g. a computer shutdown) before the R factor minimum was found, it can be re-started at the last iteration of the search the command-line option `-v <project>.ver`.

### Modifications to input and output files for CDOM and CAOI

As explained in section [chap_CSEARCH_GEN], the input files for `csearch` do not change when multiple domains and/or multiple angles of incidence are considered. All additional parameters are supplied either through environment variables (`CDOM_LEED`, `CAOI_LEED`, `CAOI_LEED`) or through the `<project>.bul` input file.

#### CDOM

The environment variable `CDOM_LEED` specifies the program which is used to calculate the IV curves for a single domain, i.e. `cleed_nsym` (See also Table [tab_csearch_env]). The only extra input required for `c_dom`, in addition to the input for `cleed_nsym`, is the symmetry of the substrate. More precisely, we need to specify the rotation and mirror symmetry operations that generate new domains on the surface. These are specified through `dr:` and `dm:` in `<project>.bul`.

<div class="description">

`dr:` `n`  
Degree of rotational symmetry. `n` can have the values 1, 2, 3, 4, or 6. Any other value will generate an error message. The default value is 1 (if no rotation symmetry is specified).

`dm:` `c/cc`  
Mirror plane. The character string following `dm:` can either be `n` (no mirror plane), `x` (mirror plane along the $`x`$-axis), `y` (mirror plane along the $`y`$-axis), or `xy` (mirror plane $`45^\circ`$ from the the $`x`$-axis). The default value is `n` (no mirror plane)



`dr:` and `dm:` do not have to be specified if there is no rotational and/or mirror symmetry. In this case `c_dom` calls `cleed_nsym` only once and generates the same output as a direct call to `cleed_nsym` would generate. Note, the symmetry operations specified here do not reflect the symmetry of the trial overlayer geometry (as specified through `sr:` and `sm:` in the parameter input file for `csearch`). In fact, they have to be different. If the overlayer geometry has the same symmetry as the substrate, there cannot be rotation or mirror domains.

`c_dom` generates a number of intermediate output files:

<div class="description">

`<project>.bul.d1`, `<project>.bul.d2`, ... (one for each domain),

`<project>.res.d1`, `<project>.res.d2`, ... (one for each domain),

`<project>.res.tmp`, and

`<project>.res.out` (containing the control output of the current LEED calculation).



The only output file used by `csearch`/`crfac` is `<project>.res`, which contains correctly averaged IV curves of the spots of all domains generated by the above symmetry operations.

Note, the `<project>.ctr` input file may have to be modified compared to a single domain as the number of spots and their averaging may be different.

#### CAOI

The environment variable `CAOI_LEED` specifies the program which is used to calculate the IV curves for a single angle of incidence single. In case of a single domain this is `cleed_nsym`; for multiple domains this is `c_dom`. `CAOI_RFAC` specifies the R factor program, i.e. `crfac` (See also Table [tab_csearch_env]).

Additional information that needs to be specified in the `<project>.bul` and/or `<project>.inp` files are the number of incidence angles and first guesses. This is done through `sa:` and `itp:`

<div class="description">

`sa:` `n`  
The number of angles of incidence for which data sets are used in the analysis.

`itp:` ` n f f`  
An integer number associated with each data set (1, 2, up to the number specified through `sa:`) taken at certain angle of incidence; the polar angle $`\theta`$ (with respect to the surface normal) and the azimuthal angle $`\varphi`$ (with respect to the $`x`$-axis) of this set.



## Setup and Run the CLEED Package

### Compile Programs under Linux

The standard distribution (`CLEED_DIS_YYMM`) contains is running versions (compiled under LINUX) of all programs (`cleed_nsym, crfac, csearch`) in the directory `.../CLEED_DIS/BIN`. It is, however, recommended to recompile the programs on your own computer. The source codes and Makefiles are in the directories `LEED_NSYM`, `RFAC`, and `SEARCH`. If the directory tree is kept as it is set up in the `*.tar` file usually only two variables need to be updated in the Makefile of each directory: `MYPRG`, the absolute path to the `CLEED_DIS` directory, and `CCOMP`, the name of the C compiler installed on the system.

    #============================================================================
    #  GH/21.06.2006
    #  makefile for the nonsymmetry LEED program
    #============================================================================
    .SUFFIXES: .o .c .h
    #=============================================================================
    MYPRG = /home/leeduser/CLEED_DIS
    INCLUDEDIR = $(MYPRG)/INCLUDE
    LIBDIR = .

    #
    # OBJECTDIR  and .c.o:  ... are defined in Makeconf
    #==flags=====================================================================
    CCOMP = cc
    OPT = O2
    CFLAGSSUB = -c -$(OPT) -I$(INCLUDEDIR) -L$(LIBDIR)
    FFLAGSSUB = -c
    CFLAGS = -$(OPT) -I$(INCLUDEDIR) -L$(LIBDIR)
    LDFLAGS = -lm
    #============================================================================

In order to compile the programs use the following steps in each of the three directories:

- Edit to Makefile to update `MYPRG` and `CCOMP`.

- Delete all object ($`\ast`$.o) files.

- Run `make`; this should produce executables called `testrfac`, `test_search`, and `test_nsym`, respectively without any error messages. Warning messages may appear depending on the compiler.

- Copy the executables into `.../CLEED_DIS/BIN` and rename them to `crfac`, `csearch`, `cleed_nsym`.

On other operating systems check the compatibility of the trigonometric functions and data type used for real in `real.h` (in directory `INCLUDE`)

### Set up a Search

A set of sample input files for Er/Si{111} (no expt. IV curves), Benzene on Ru{0001} and O and Cu/Ni{111} can be found in `EXAMPLES/SIER`, `RUBENZ`, `NIO` and `NICU`. There are three files for each search (in addition to phase shifts and experimental data) which have to be created by the user:

- `*.bul`: bulk geometry and non-geometric parameters

- `*.inp`: start geometry for the search and search parameters (symmetry etc.)

- `*.ctr`: correspondence between theoretical and experimental IV curves.

Detailed descriptions of the file formats can be found in this manual (Sections [sec_leed_inp], [sec_rfac_io], [sec_search_io], but modifying the sample input files are probably the best way to start a new search. All other files are created by the programs during the search. The search can be monitored through the `*.log` file; at each point of the search the current best fit geometries and IV curves can be found in `*.rmin`, `*.pmin`, and `*.bmin` (for a description see the manual).

The phase shift input is described in Section [sec_leed_inp]. It must have the same format as for the VanHove/Tong programs with a first line added that contains the number of energies and $`l_{max}`$. For each type of atoms there must be a separate phase shift file in the directory `.../CLEED_DIS/PHASE`, called `<atom>.phs`. The directory is passed to the LEED program through the environment variable `CLEED_PHASE`. Before starting the search at least 3 environment variables have to be set:

- `CLEED_PHASE`: name of the phase shift directory (used by `cleed_nsym`).

- `CSEARCH_LEED`: name of the executable LEED-IV program, normally `cleed_nsym` (used by `csearch`).

- `CSEARCH_RFAC`: name of the executable R factor program, normally `crfac` (used by `csearch`).

a shell script performing these assignments can be found in `.../CLEED_DIS/BIN/set_env` and is shown in Table [tab_setup_set_env]. Note, the directory path of the parent directory (`CLEED_HOME`) must be changed accordingly.

            export CLEED_HOME=/home/leeduser/CLEED_DIS
            export CLEED_PHASE=$CLEED_HOME/PHASE
            export CSEARCH_LEED=$CLEED_HOME/BIN/cleed_nsym
            export CSEARCH_RFAC=$CLEED_HOME/BIN/crfac

### Start a Search

We use the example given in `CLEED_DIS/EXAMPLES/NIO`. The directory contains all input files, `Ni111_2x2O.inp`, `Ni111_2x2O.bul`, `Ni111_2x2O.ctr`, as well as the data files `nio_*.fsm`. The search is started from the directory that contains the input files by

`.../CLEED_DIS/BIN/csearch -i Ni111_2x2O.inp`

`csearch` is the master program that calls the LEED-IV program (`cleed_nsym`) and the R factor program (`crfac`) to calculate IV curves and R factors for a given trial structure and optimises the geometry parameters of the surface geometry in order to minimise the R factor.

The LEED-IV program can be called from outside the search by

`...CLEED_DIS/BIN/cleed_nsym -i Ni111_2x2O.par -b Ni111_2x2O.bul`  
`-o Ni111_2x2O.res 1> Ni111_2x2O.out`

`Ni111_2x2O.par` is usually created by the search program. It contains only the positions of the overlayer atoms that are optimised during the search. The IV curves are written to `Ni111_2x2O.res` (indicated by the ”-o” option). In addition, the LEED-IV program produces a large amounts of control output, which is written to ”stdout”. It is best to re-directed ”stdout” to a file using the ”`1>`” command (in bsh).

After significant alterations to the input file the LEED-IV program should always be checked separately first to ensure it produces sensible output. The R factor program is called by

`...CLEED_DIS/BIN/crfac -c Ni111_2x2O.ctr -t Ni111_2x2O.res`

Once the user is satisfied that all parts of the IV calculation run smoothly, the search can be started. It is recommended to use the ’`nohup`’ command together with ’`&`’. This runs the search in the background:

`nohup .../CLEED_DIS/BIN/csearch -i Ni111_2x2O.inp 1> out 2> err & `

In the above example each iteration takes about 90-100s (2.4 GHz Linux PC), the final R factor after convergence should be around 0.1322 and it takes around 183 iteration to get there from the start geometry specified in `Ni111_2x2O.inp` (see $`\ast`$.log file, which is also included). The example in NICU leads to an R factor of 0.0633 after 63 iterations (around 15 s cpu time per iteration).

All information necessary to restart the search at the current position is stored in $`\ast`$`.ver`.If the search stops for any reason, e.g. because it has reached the limit of iterations (currently 2000), it can be restarted by:

`nohup  /CLEED_DIS/BIN/csearch -i Ni111_2x2O.inp -v Ni111_2x2O.ver 1> out 2> err &`

# Appendix

## Other Useful Programs

In the following we describe a few little programmes that are quite useful in the preparation of a LEED-IV analysis. They can be found in the directory "MISC", which is included in the CLEED package.

### Smoothing of Data

#### Fourier Smoothing

The program `ftsmooth` performs a Fourier transformation of the original data. The high-frequency part of the Fourier spectrum (noise) is cut off and back-transformed. The energy values of the output spectrum are the same as in the input spectrum. When the program is compiled the mathematical functions library has to be linked (`cc ftsmooth.c -lm`).

##### Command line

    ftsmooth -i <input> -o <output>
    -m <s(n): (no) smooth> -c <cutoff(0.5)> -t <tailoff(10)>\n");

##### Parameters

<div class="description">

`-m s` smooth; `-m n`: no smooth.

`-c` cutoff in Fourier transform; default: 0.5 (corresponds to 4 eV)

`-t` tail-off in Fourier transform; default: 10



#### Three-Point Smoothing

The program `smooth` performs a three-point smooth of the original data: $`I_{sm}(i) = \frac{1}{4} [I(i-1) + 2 \cdot I(i) + I(i+1]`$. The energy values of the output spectrum are the same as in the input spectrum.

##### Command line

    smooth <input> <output>

### LEED Pattern (patt)

The program `pattern` prints the LEED pattern for a given superstructure specified in the input file. The output file is in EPS format.

##### Command line

    pattern -i <input> -o <output> -ni -rs <spot size> -rg <spot size>

##### Parameters

<div class="description">

`-i <input>`: specify input file name.

`-o <output>`: specify output file name.

`-ni` (no argument): do not print spot indices.

`-ns` (no argument): do not use different symbols for superstructure spots.

`-rs <spot size>`: spot size of SS spots in points (PS units).

`-rg <spot size>`: spot size of GS spots in points (PS units).



##### Format of input file

An Example of an input file or two Domains of a $`p(\sqrt{7} \times \sqrt{7}) R19^\circ`$ structure on a hexagonal surface is given in Table [tab_OTHER_patt]:

    c 2 Domains of (3 2 / -2 1) = (r7 x r7)
    1.0  1.732    (a1  <lattice vector>)
    1.0 -1.732    (a2  <lattice vector>)
    1.5           (radius <Radius = longest rec. lattice vector * radius>)
    2             (number of domains)
    #M1
      3   2  M1   (<first domain>)
     -2   1  M1
    #M2
    Sy            (<second domain: mirror image of first domain w/r y-axis>)

”`#`” specifies a comment which does not appear in the plot. ”`c`” specifies a comment which appears as title of plot. There two different ways of specifying the superstructure matrix. First by specifying the actual matrix:

    m11 m12  (actual matrix)
    m21 m22

The matrix of the first domain has to be specified like this. The following domains can be specified by symmetry operations:

<div class="description">

`R<phi>`: rotation of previous matrix by `<phi>` degrees

`Sx`: mirror previous matrix w/r x-axis

`Sy`: mirror previous matrix w/r y-axis



### Surface Coordinates (latt)

The program `latt` produces an output file in XYZ format which contains the atom coordinates for a given surface.

##### Command line

    latt -h <h> -k <k> -l <l> -a <latt const. a> -c <latt const. c>
         -n <name of atom> -t <lattice type>
         -i <input of basis vectors> -o <output> (default is stdout)

##### Parameters

<div class="description">

`-i <input>`: input of file containing basis vectors for a general lattice. If an input file is specified, the options `-a, -c, -t, -n` are ignored.

`-o <output>`: specify output file name (default is ”stdout”). Two output files are created, the file with the name specified by -o contains the XYZ data, another file with the extra extension ”`.inf`” contains additional information.

`-h <h> -k <k> -l <l>`: Miller indices.

`-a <latt const. a> `: lattice constant for cubic lattices.

`-c <latt const. c>`: lattice constant c for hcp lattice.

`-n <name of atom>`: name of atoms.

`-t <lattice type>`: lattice type; implemented types are: `fcc, bcc, hcp, dia`.



##### Format of input file

Table [tab_OTHER_latt] shows the input file for hexagonal ice.

    # Hexagonal ice
    a1:     2.3469     -4.0650     0.0000
    a2:     2.3469      4.0650     0.0000
    a3:     0.0000      0.0000     7.3600
    #
    pb:  O  0.0000      0.0000     0.0000
    pb:  H  0.8660     -0.5000    -0.3580
    pb:  H  0.0000      0.0000     1.0000
    #
    pb:  O  0.0000     +2.7100    -0.9706
    pb:  H  0.8660     +3.2100    -0.6107
    pb:  H  0.0000     +1.7100    -0.6107
    #
    pb:  O  0.0000     +2.7100    -3.6806
    pb:  H  0.8660     +3.2000    -4.0386
    pb:  H  0.0000     +2.7100    -2.6806
    #
    pb:  O  0.0000      0.0000    -4.6500
    pb:  H  0.8660     -0.5000    -4.2920
    pb:  H  0.0000     +1.0000    -4.2920
    #
    il: 20.
    nl: 2

`il` and `nl` specify the number of unit cells in lateral and vertical direction, respectively.

### Calculation of Phase Shifts

For the calculation of phase shifts we recommend the Barbieri/Van Hove phase shift package, which can be downloaded from:

    http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/leedpack.html
