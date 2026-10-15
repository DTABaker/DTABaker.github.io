---
layout: post
title: "Orbital Observatory: A Sense of Scale"
short_title: "Orbitals on the same ruler"
date: 2026-10-14 00:00:00 -0500
permalink: /orbital-scale/
topic: "Chemistry & teaching"
tags: ["Teaching & outreach"]
description: "Compare sixteen hydrogenic orbitals on one physical scale, explore nuclear-charge contraction, and connect three-dimensional probability surfaces to radial probability distributions."
---

Orbital drawings often make every shape about the same size. That helps us compare their angular patterns, but hides their relative spatial extent. **A Sense of Scale** puts the sixteen Orbital Observatory states on one shared ruler, with synchronized rotation and zoom.

This is a companion to [the original Orbital Observatory]({{ '/2026/09/07/orbital-observatory.html' | relative_url }}) and [From Wavefunction to Orbital]({{ '/wavefunction-to-orbital/' | relative_url }}). The original viewer emphasizes shape; the wavefunction model connects mathematical factors to a spatial distribution. This page makes physical size the comparison.

<p><a href="{{ '/assets/orbitals-to-scale.html' | relative_url }}" target="_blank" rel="noopener">Open A Sense of Scale in a full window ↗</a></p>

<iframe
  id="orbital-scale-embed"
  src="{{ '/assets/orbitals-to-scale.html' | relative_url }}"
  title="Orbital Observatory: A Sense of Scale, sixteen orbitals with one shared physical ruler"
  width="100%" height="1600"
  style="display:block; width:100%; border:1px solid #d8dfdf; border-radius:10px;"
  loading="lazy">
</iframe>

<noscript><p>The interactive model needs JavaScript and WebGL. The explanation and reference values below remain available without JavaScript.</p></noscript>

Scroll within the model to reach the radial graphs and assumptions. The full-window link gives the comparisons more room, particularly on a phone.

## Does increasing nuclear charge really make an orbital smaller?

**Yes, for a fixed orbital state in this one-electron model.** Nuclear charge Z is the number of protons. Increasing it strengthens the attraction between the nucleus and the same single electron. In the nonrelativistic, fixed-nucleus Coulomb model, every characteristic length contracts in proportion to **1/Z**.

For the same state, He⁺ (Z = 2) therefore has half the characteristic size of H (Z = 1); Li²⁺ (Z = 3) has one-third, and Be³⁺ (Z = 4) has one-quarter. The camera scale stays fixed when you change Z, so the contraction is visible rather than canceled by automatically fitting the orbital to its panel.

The mean electron–nucleus distance is

<p class="orbital-equation"><strong>⟨r⟩ = [3n² − ℓ(ℓ + 1)]a₀ / (2Z).</strong></p>

Here n is the principal quantum number, ℓ is the orbital angular-momentum quantum number, and a₀ is the Bohr radius, approximately 0.529177 Å. For 1s, ⟨r⟩ is 1.5a₀ in H and 0.75a₀ in He⁺. The normalized three-dimensional probability density increases as the distribution contracts: ρ<sub>Z</sub>(**r**) = Z³ρ<sub>1</sub>(Z**r**), keeping its integral equal to one.

**He⁺ is not neutral helium.** Neutral many-electron atoms also involve electron repulsion and shielding. Adding a proton and an electron changes more than Z in a one-electron equation. This viewer should not be read as a quantitative sequence of neutral-atom sizes or as an exact rule for all periodic trends.

## Try these comparisons

1. **Keep the ruler fixed.** Start at Z = 1 and compare 1s, 2p, 3d, and 4f. Their mean radii are 1.5, 5, 10.5, and 18a₀. These are expectation values, not the outer edges of the displayed surfaces.
2. **Change only Z.** Switch from H to He⁺. The surfaces and mean radii contract by one-half. Return to H and notice that the common viewing window is unchanged.
3. **Change units.** Switch between Bohr radii and ångströms. The labels and ruler convert; the physical orbital sizes do not change.
4. **Compare angular forms.** Use the orbital selectors, or show all sixteen orbitals. States in the same subshell share the same radial distribution and mean radius, although their angular shapes and farthest surface points can differ. Not every real d or f angular form is simply a rotation of another.
5. **Read the radial graph.** Its curves show P(r) = r²&#124;R(r)&#124;², the probability per unit radial distance. For these n = ℓ + 1 states, the radial-probability peak is at n²a₀/Z. A graph peak, mean radius, and isodensity boundary answer different questions.

Drag any panel to rotate all of them, use the front/top/3D buttons, or use arrow keys on a focused view. Zoom is shared across panels. The radial chart's horizontal distance axis stays fixed as Z changes, while the labeled probability-density axis rescales.

## What the surfaces mean

The collection contains **1s, three 2p, five 3d, and seven 4f real orbitals**. These states have no radial nodes. They are alternative states of a one-electron system, not sixteen electrons occupying one atom, and they do not include every orbital through n = 4.

Each displayed surface approximates the boundary of an orbital-specific isodensity region containing 90% of that state's probability, summed over all lobes and both signs. Probability extends outside the region, and there is no hard atomic edge. The chosen density threshold is different for different orbitals. Blue and coral represent wavefunction sign, not opposite charge.

Every panel uses the same orthographic viewing width and physical coordinate system. The small nucleus marker and axes are annotations, not physical-size models of a nucleus. The meshes are finite-resolution approximations; the underlying analytic regions and coordinate conversions have independent numerical checks. These hydrogenic shapes do not provide quantitative sizes for orbitals in neutral many-electron atoms.

For the mathematical background, see [Oregon State University's hydrogen radial functions and spherical harmonics](https://books.physics.oregonstate.edu/GMM/hydroform.html) and Richard Fitzpatrick's [Hydrogen Atom](https://farside.ph.utexas.edu/teaching/qmech/lectures/node82.html). The Bohr-radius conversion uses [NIST's physical constants](https://physics.nist.gov/cuu/Constants/Table/allascii.txt). The bundled Three.js renderer retains its [MIT license]({{ '/assets/THREE-LICENSE.txt' | relative_url }}).
