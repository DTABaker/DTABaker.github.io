---
layout: post
title: "Orbital Observatory: How Orbitals Mix"
short_title: "How orbitals mix"
date: 2026-10-17 12:00:00 -0500
permalink: /hybrid-orbitals/
topic: "Chemistry & teaching"
tags: ["Teaching & outreach"]
description: "Change s and p character, compare wavefunctions with probability densities, and explore complete sp, sp², and sp³ orbital sets."
---

What changes when we mix an s orbital with a p orbital? The answer depends on what we add. Adding **wavefunctions** allows constructive and destructive interference. Adding their **probability densities** leaves out that interference and gives a different distribution.

The **Hybrid Orbital Lab** makes that distinction interactive. Move the s-character slider, reverse the relative sign, and compare three linked views of the same mathematical ingredients. Then explore how one hybrid fits into a complete sp, sp², or sp³ set.

This continues the [Orbital Observatory]({{ '/2026/09/07/orbital-observatory.html' | relative_url }}) series, alongside [From Wavefunction to Orbital]({{ '/wavefunction-to-orbital/' | relative_url }}) and [A Sense of Scale]({{ '/orbital-scale/' | relative_url }}).

<p><a href="{{ '/assets/hybrid-orbital-lab.html' | relative_url }}" target="_blank" rel="noopener">Open the Hybrid Orbital Lab in a full window ↗</a></p>

<iframe
  id="hybrid-orbital-embed"
  src="{{ '/assets/hybrid-orbital-lab.html' | relative_url }}"
  title="Hybrid Orbital Lab: s and p character, interference, and complete hybrid sets"
  width="100%" height="1600"
  style="display:block; width:100%; border:1px solid #d8dfdf; border-radius:10px;"
  loading="lazy">
</iframe>

Scroll within the model to reach the complete sets and model notes. The full-window link gives the comparisons more room, particularly on a phone.

## Character is a squared coefficient

For the normalized, orthogonal s and p basis functions used here, one hybrid is

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

The equivalent hybrids have 50%, 33⅓%, and 25% s character, respectively. Their direction vectors form the familiar 180°, 120°, and approximately 109.47° angles. These are ideal orbital constructions, not predictions of every molecular bond angle. A single orbital with 33⅓% s character does not by itself specify a trigonal arrangement.

The complete four-orbital bases remain orthonormal. Summing all four orbital densities, including unused p functions, gives the same unweighted sum as the original 2s/2p basis. That identity does not specify a molecule's occupied orbitals or electron density.

## What this model represents

The functions are normalized hydrogenic **2s, 2pₓ, 2pᵧ, and 2p<sub>z</sub> at Z = 1**, in a fixed-nucleus, nonrelativistic one-electron model. The slider scans coefficients; it is not a movie of an atom forming a bond, an excitation, or a chemical reaction. The display illustrates the mathematics of hybridization and does not calculate carbon bonding orbitals, molecular energies, or bond strengths.

The 2s basis function has a radial node at 2a₀, where a₀ is the Bohr radius. A hybrid's nodes change through interference and need not retain that spherical node. This radial structure is one reason these surfaces can look different from simplified textbook hybrid lobes. A cutaway and the analytic slices help reveal the interior structure.

Colors in the wavefunction view indicate positive and negative amplitude, not electrical charge. The wavefunction and hybrid-density panels use **the same boundary geometry**, colored differently. All three comparison panels share one physical scale. The complete-set panel has its own scale bar.

Surfaces are drawn at a selected fixed probability-density cutoff in a₀⁻³. They do **not** enclose a specified percentage of probability, and they are not hard orbital edges. Small features near a cutoff may be unresolved by the finite mesh. In From Wavefunction to Orbital and A Sense of Scale, each orbital uses its own cutoff chosen for approximately 90% probability. Here, the comparison uses one selected density cutoff without fixing the enclosed probability.

For the underlying mathematics, see [Oregon State University's hydrogen radial functions and spherical harmonics](https://books.physics.oregonstate.edu/GMM/hydroform.html) and [TU Braunschweig's normalized hybrid-orbital constructions](https://www.pci.tu-bs.de/aggericke/PC4e/Kap_II/Hybride/Konstruktion.htm). The bundled Three.js renderer retains its [MIT license]({{ '/assets/THREE-LICENSE.txt' | relative_url }}).
