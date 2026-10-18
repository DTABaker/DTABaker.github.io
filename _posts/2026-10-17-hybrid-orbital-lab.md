---
layout: post
title: "Orbital Observatory: How Orbitals Mix"
short_title: "How orbitals mix"
date: 2026-10-17 12:00:00 -0500
permalink: /hybrid-orbitals/
topic: "Chemistry & teaching"
tags: ["Teaching & outreach"]
description: "Continuously vary s/p/d orbital character, compare wavefunctions and densities, and explore complete sp, sp², sp³, dsp², sp³d, and sp³d² constructions."
---

What changes when we mix an s orbital with a p orbital? The answer depends on what we add. Adding **wavefunctions** allows constructive and destructive interference. Adding their **probability densities** leaves out that interference and gives a different distribution.

The **Hybrid Orbital Lab** makes that distinction interactive. Move the s-character slider, reverse the relative sign, and compare three linked views of the same mathematical ingredients. Then explore how one hybrid fits into a complete sp, sp², or sp³ set. A further section brings d orbitals into the basis through square-planar dsp², trigonal-bipyramidal sp³d, and octahedral sp³d² constructions, with linked sliders to vary one orbital’s s, p, and d character.

This continues the [Orbital Observatory]({{ '/2026/09/07/orbital-observatory.html' | relative_url }}) series, alongside [From Wavefunction to Orbital]({{ '/wavefunction-to-orbital/' | relative_url }}) and [A Sense of Scale]({{ '/orbital-scale/' | relative_url }}).

<p><a href="{{ '/assets/hybrid-orbital-lab.html' | relative_url }}" target="_blank" rel="noopener">Open the Hybrid Orbital Lab in a full window ↗</a></p>

<iframe
  id="hybrid-orbital-embed"
  src="{{ '/assets/hybrid-orbital-lab.html' | relative_url }}"
  title="Hybrid Orbital Lab: s/p character, interference, and complete s/p/d hybrid constructions"
  width="100%" height="1600"
  style="display:block; width:100%; border:1px solid #d8dfdf; border-radius:10px;"
  loading="lazy">
</iframe>

Scroll within the model to reach the complete s/p sets, the section including d orbitals, and the model notes. The full-window link gives the comparisons more room, particularly on a phone.

## Character is a squared coefficient

For the normalized, orthogonal 2s and 2p basis functions in the slider section, one hybrid is

<p><strong>h = √f φ<sub>s</sub> + σ√(1 − f) φ<sub>p</sub>.</strong></p>

Here f is the s-character fraction, σ is +1 or −1, and φ<sub>p</sub> is the selected 2p function. The coefficient of the s function is √f, so **25% s character means an s coefficient of 0.5**, not 0.25. These percentages describe components in the chosen orbital basis, not separate pieces of an electron occupying different lobes.

Squaring the combined wavefunction gives the hybrid probability density. For these real functions,

<p><strong>ρ = f&#124;φ<sub>s</sub>&#124;² + (1 − f)&#124;φ<sub>p</sub>&#124;² + 2σ√[f(1 − f)]φ<sub>s</sub>φ<sub>p</sub>.</strong></p>

The last term is the interference term. It can increase or decrease the density at a point; the total density remains nonnegative. The density-only mixture omits this cross term. It represents an incoherent mixture of the two basis states, rather than the same pure hybrid state drawn another way.

## Try these comparisons

1. **Reverse the relative sign.** Start at 50% s character. Changing the p contribution from + to − reflects the hybrid density along the selected p axis. The density-only mixture stays unchanged. Reversing the overall sign of the entire wavefunction, by contrast, would leave its density unchanged.
2. **Inspect a point.** Move the shared probe across a slice and compare hybrid density, mixture density, and their difference. At a node of either component, the cross term vanishes.
3. **Visit both endpoints.** At pure s or pure p, one coefficient is zero. There is no interference between s and p, so the two density views agree everywhere.
4. **Build the whole set.** The sp construction contains two hybrids and two unused p orbitals; sp² has three hybrids and one unused p orbital; sp³ has four hybrids. Select individual orbitals when their surfaces overlap.

The equivalent sp, sp², and sp³ hybrids have 50%, 33⅓%, and 25% s character, respectively. Their direction vectors form the familiar 180°, 120°, and approximately 109.47° angles. These are ideal orbital constructions, not predictions of every molecular bond angle. A single orbital with 33⅓% s character does not by itself specify a trigonal arrangement.

The complete four-orbital bases remain orthonormal. Summing all four orbital densities, including unused p functions, gives the same unweighted sum as the original 2s/2p basis. That identity does not specify a molecule's occupied orbitals or electron density.

## Bringing d orbitals into the basis

The d-containing section uses nine normalized hydrogenic functions: **3s, three 3p functions, and five real 3d functions**, all at Z = 1. Choose a construction, inspect one orbital at a time, or overlay the set. Three side-by-side panels show its **wavefunction, hybrid density, and density-only mixture**, linked to the character sliders. They share a camera, spatial scale, and density cutoff. The wavefunction and hybrid-density panels share exactly the same boundary geometry. Adding the unused orbitals to an overlay shows separate functions, not their summed wavefunction, summed density, or a molecular electron density.

| Construction | Character of each hybrid: s / p / d | Unused functions in the nine-function basis |
| --- | --- | --- |
| dsp²: four square-planar hybrids | 25% / 50% / 25% | Five: p<sub>z</sub>, d<sub>xy</sub>, d<sub>xz</sub>, d<sub>yz</sub>, d<sub>z²</sub> |
| sp³d: three equatorial and two axial hybrids | Equatorial: 33⅓% / 66⅔% / 0%<br>Axial: 0% / 50% / 50% | Four: d<sub>xy</sub>, d<sub>xz</sub>, d<sub>yz</sub>, d<sub>x²−y²</sub> |
| sp³d²: six octahedral hybrids | 16⅔% / 50% / 33⅓% | Three: d<sub>xy</sub>, d<sub>xz</sub>, d<sub>yz</sub> |

The sp³d set is particularly instructive: **its five orbitals are not equivalent**. In the orthonormal convention used here, the three equatorial functions are sp² combinations, while the axial pair mixes p<sub>z</sub> with d<sub>z²</sub>. The set name counts the input functions; it does not mean every member contains every orbital type. The 20% s / 60% p / 20% d composition is the average across all five, not the composition of each one. The coefficient convention follows [Wannier90's explicit hybrid table](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/projections/#hybrids); the display uses hydrogenic n = 3 radial functions rather than calculating Wannier functions.

The three linked character sliders let you change a selected orbital continuously while keeping the total at 100%. Moving one proportionally rebalances the other two; if both are zero, they share the remainder equally. **Restore selected orbital** returns to the exact preset. The sliders edit one orbital, so an overlay must first be changed to an individual selection.

During a scan, the selected reference supplies fixed normalized s, p, and d components, including their internal signs. Their weights are the square roots of the character fractions. If a component is absent from the reference, the editor explicitly supplies positive 3s, 3p<sub>z</sub>, or 3d<sub>z²</sub>, as appropriate. A custom orbital remains normalized, but substituting it into the reference set would generally break orthogonality. The original complete sets remain available unchanged.

Calling those normalized components S, P, and D, the displayed hybrid is h = √f<sub>s</sub>S + √f<sub>p</sub>P + √f<sub>d</sub>D, with f<sub>s</sub> + f<sub>p</sub> + f<sub>d</sub> = 1. The density-only mixture is

<p><strong>ρ<sub>mix</sub> = f<sub>s</sub>&#124;S&#124;² + f<sub>p</sub>&#124;P&#124;² + f<sub>d</sub>&#124;D&#124;².</strong></p>

It removes interference **between the s, p, and d components**, while keeping the chosen shape within each group. In particular, P or D may itself combine several p or d basis functions whose internal interference remains. The difference between the hybrid density and this mixture consists of the three cross terms: 2√(f<sub>s</sub>f<sub>p</sub>)SP + 2√(f<sub>s</sub>f<sub>d</sub>)SD + 2√(f<sub>p</sub>f<sub>d</sub>)PD. Both densities are nonnegative and integrate to one for a normalized selected orbital, but they describe different states. Try moving a character slider to 100%: the two density panels agree everywhere, including for an oriented p or d combination.

The **linked cross-sections and shared probe** let students check this distinction at the same spatial point. Choose the xy, xz, or yz plane through the nucleus, then point to a slice, tap it, use its arrow keys, or enter coordinates. Each slice spans −24 to +24a₀ along both displayed axes. Changing the plane resets the probe to the nucleus; changing character keeps its coordinates. The readout shows signed amplitude, both densities, their difference, and the individual s–p, s–d, and p–d interference terms. Select one orbital to use the slices and probe; an overlay has no single combined wavefunction to sample.

These are values of the analytic functions on a plane, **not projections or integrals through the orbital**. The 3D cutoff, cutaway, and camera do not change them. The two density slices share a fixed square-root color scale that saturates at 0.0004a₀⁻³; the signed-amplitude scale saturates at ±0.02a₀<sup>−3/2</sup>. Saturation affects colors only: the probe reports the underlying values without clipping. A slice lying in an angular nodal plane can be blank even though the orbital is nonzero elsewhere—try another plane. Density and interference readouts have units a₀⁻³; signed amplitude has units a₀<sup>−3/2</sup>.

Character remains a sum of squared coefficients. For example, an octahedral hybrid can contain both d<sub>z²</sub> and d<sub>x²−y²</sub>; their squared coefficients add to its total d character. For each reference preset, including all unused functions completes an orthonormal nine-function basis, and its unweighted sum of orbital densities equals that of the original nine functions. This does not assign nine electrons to the basis.

These constructions illustrate symmetry and interference. **A molecular shape does not establish substantial d hybridization.** The historical sp³d and sp³d² descriptions of main-group molecules such as PF₅ and SF₆ should not be taken as their bonding mechanism; ab initio studies support much smaller d-orbital contributions and descriptions involving polar or multicenter bonding. See [Noury, Silvi, and Gillespie's study of hypervalent bonding](https://doi.org/10.1021/ic011003v). In transition-metal compounds, the metal–ligand interactions, orbital energies, and occupations must be considered together, for example in [molecular-orbital and ligand-field treatments](https://doi.org/10.1103/PhysRevB.85.165113).

The n = 3 basis also retains radial structure: 3s has two radial nodes, 3p has one at r = 6a₀, and 3d has none. Direction labels describe the symmetry of the coefficient construction; radial phase changes mean they need not identify the largest lobe at every radius. This section has its own scale bar and lower density cutoffs. It does not model a transition metal's 3d/4s/4p shell ordering, ligand orbitals, or spin state.

## What this model represents

The s/p slider and first complete-set section use normalized hydrogenic **2s, 2pₓ, 2pᵧ, and 2p<sub>z</sub> at Z = 1**; the d-containing section uses the separate n = 3 basis described above. Both use a fixed-nucleus, nonrelativistic one-electron model. Functions within a given n shell are degenerate in this model, which does not establish equal orbital energies in a many-electron atom or molecule. The slider scans coefficients; it is not a movie of an atom forming a bond, an excitation, or a chemical reaction. The display illustrates the mathematics of hybridization and does not calculate carbon bonding orbitals, molecular energies, or bond strengths.

The 2s basis function has a radial node at 2a₀, where a₀ is the Bohr radius. A hybrid's nodes change through interference and need not retain that spherical node. This radial structure is one reason these surfaces can look different from simplified textbook hybrid lobes. A cutaway and the analytic slices in the s/p comparison help reveal the interior structure.

Colors in the wavefunction view indicate positive and negative amplitude, not electrical charge. The wavefunction and hybrid-density panels use **the same boundary geometry**, colored differently. Each three-panel comparison shares one physical scale internally. The n = 3 comparison has its own 5a₀ scale bars and lower cutoffs, distinct from the n = 2 views.

Surfaces are drawn at a selected fixed probability-density cutoff in a₀⁻³. They do **not** enclose a specified percentage of probability, and they are not hard orbital edges. Small features near a cutoff may be unresolved by the finite mesh. In From Wavefunction to Orbital and A Sense of Scale, each orbital uses its own cutoff chosen for approximately 90% probability. Here, the comparison uses one selected density cutoff without fixing the enclosed probability.

For the underlying mathematics, see [Oregon State University's hydrogen radial functions and spherical harmonics](https://books.physics.oregonstate.edu/GMM/hydroform.html) and [TU Braunschweig's normalized hybrid-orbital constructions](https://www.pci.tu-bs.de/aggericke/PC4e/Kap_II/Hybride/Konstruktion.htm). The bundled Three.js renderer retains its [MIT license]({{ '/assets/THREE-LICENSE.txt' | relative_url }}).
