> **Research input — Chapter 1 content research.** Raw output of an Opus research agent (2026-09-10), kept for deep dives. It has been verified by Fable agents and condensed into `notes/chapters/ch01-what-to-build.md`, which is the authoritative, corrected version. Verifier verdict: consistent (minor nits only — e.g. a few internal cross-references and the 1.37/1.38 "projectile" wording, both corrected in the final report). No text from Taylor's book is reproduced; problems are paraphrased.

---

# Taylor, *Classical Mechanics* (2005) — Chapter 1 Research Report
### For a student building a C++ force-based particle simulator alongside the course

**Tagging convention used here (strict):**
- **[WEB]** = I actually fetched the source this session (URL given inline). This includes the 51 Chapter-1 solution PDFs from stemjock.com, which I downloaded and text-extracted locally, so their **problem statements** are directly verified.
- **[MEM]** = my own knowledge of the book, *not* independently confirmed by a page I fetched. Some [MEM] items were corroborated only by a search-engine summary (never a fetched page) — those are called out explicitly and repeated in "Could not confirm."

No text, problem statement, or figure from Taylor is reproduced here. Every problem is paraphrased in my own words. I did not fetch, and deliberately avoided, the several pirated full-text copies that appeared in search results (archive.org full-text dump, pdfcoffee, silo.tips, slideshare, studylib, scribd, zlib.pub).

I read `/Users/jackery/Github/newtonian-library/docs/research/taylor-ch1-7-map.md` (lines 38–122) as a starting hypothesis. **It is mostly right but contains at least two factual errors about Chapter 1** — see §3.6 and §6 below.

---

## 1. Section-by-section map of Chapter 1

Chapter structure (seven numbered sections, then a summary page, then the problems):

| § | Title | Start page |
|---|---|---|
| 1.1 | Classical Mechanics | 3 |
| 1.2 | Space and Time | 4 |
| 1.3 | Mass and Force | 9 |
| 1.4 | Newton's First and Second Laws; Inertial Frames | 13 |
| 1.5 | The Third Law and Conservation of Momentum | 17 |
| 1.6 | Newton's Second Law in Cartesian Coordinates | 23 |
| 1.7 | Two-Dimensional Polar Coordinates | 26 |
| — | Principal Definitions and Equations of Chapter 1 | 33 |
| — | Problems for Chapter 1 | 34 |

**[MEM]** for the section titles and page numbers as a table. Search summaries agreed on the titles and on 1.7 → p. 26 and problems → p. 34, but the page I would have fetched to confirm (`henry.pha.jhu.edu/chapters.html`) returned HTTP 403, and the other TOC hits were pirated copies I did not open. **[WEB]** for one page number: the stemjock solution to problem 1.50 states that equation (1.51) is on page 31, which is consistent with §1.7 running from 26 into the low 30s — https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p50.pdf

The end-of-chapter problems run **1.1 through 1.51 — exactly 51 problems** **[WEB]**, confirmed independently by two fetched sources: the stemjock index (https://stemjock.com/taylorcm.htm) and StudySoup's Chapter 1 listing (https://studysoup.com/tsg/physics/90/classical-mechanics/chapter/1765/1). Taylor groups them by section, so the problem numbers themselves tell you which section a topic belongs to (see §3).

---

### 1.1 Classical Mechanics — the domain of validity

Key idea: Newtonian mechanics is the theory of macroscopic bodies moving slowly compared with light. It is superseded by relativity at high speed and by quantum mechanics at small scale, but remains exact enough for essentially all engineering and most astronomy. **[MEM]** — corroborated in spirit by an independent fetched university lecture-notes PDF (Brod, U. Cincinnati) which makes the same framing point about relativity and quantum mechanics refining Newton's laws: https://homepages.uc.edu/~brodjm/teaching/lecture_notes/mechanics.pdf

No equations. Nothing to implement.

---

### 1.2 Space and Time — vectors and vector calculus

This is where Taylor sets up the entire vector toolkit the rest of the book uses. This section is the direct source of your `Vec2`/`Vec3` requirements.

Key equations (plain text):

```
r         = (x, y, z)                        position vector, Cartesian components
r         = x x_hat + y y_hat + z z_hat      unit-vector form
r         = sqrt(x^2 + y^2 + z^2)            magnitude (Pythagoras applied twice)

r + s     = (r_x + s_x, r_y + s_y, r_z + s_z)
c r       = (c r_x, c r_y, c r_z)

r . s     = r s cos(theta)                   (1.6)  geometric definition
r . s     = sum_i r_i s_i                    (1.7)  component definition
                                             1.6 and 1.7 are equal - proved in problem 1.7

r x s     = (r_y s_z - r_z s_y,
             r_z s_x - r_x s_z,
             r_x s_y - r_y s_x)              (1.9)  component definition
|r x s|   = r s sin(theta),  direction by right-hand rule

d/dt (r . s) = r . (ds/dt) + (dr/dt) . s     product rule, dot
d/dt (r x s) = r x (ds/dt) + (dr/dt) x s     product rule, cross - ORDER MATTERS

v = dr/dt = r_dot = (x_dot, y_dot, z_dot)
a = dv/dt = r_ddot
```

**[WEB]** for the equation-number bindings (1.6) = the `rs cos θ` form, (1.7) = the component sum, (1.9) = the component cross product: these numbers are cited verbatim inside the fetched statements of problems 1.7, 1.8, 1.15, 1.16, and 1.17 — e.g. https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p7.pdf and .../TaylorCMCh1p15.pdf

**[MEM]** that Taylor also introduces the "BAC−CAB" triple-product rule in or near this section — though the fetched statement of problem 1.33 explicitly instructs the student to use it, which confirms the book expects the reader to have it: **[WEB]** https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p33.pdf

---

### 1.3 Mass and Force

Key ideas: mass as the measure of inertia (operationally defined by comparing accelerations under the same force); force as the vector quantity that appears on the right of Newton II; forces add as vectors, so the net force is the vector sum. Weight near the Earth's surface is `w = m g`. **[MEM]** — I found no fetched page confirming the specifics of this section's treatment.

```
F_net = sum over all forces acting                (superposition of forces)
w     = m g                                       (weight)
```

Library implication: forces must be **summed**, never replaced. That is the single most important API consequence of §1.3.

---

### 1.4 Newton's First and Second Laws; Inertial Frames

Key ideas:
- **First law**: a body subject to zero net force moves with constant velocity (possibly zero). The first law is not a special case of the second — it *defines* the class of frames in which the second law holds. **[WEB]** — the fetched statements of problems 1.26 and 1.27 both open by stating exactly this as the defining property of an inertial frame: https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p26.pdf
- **Second law**: `F = dp/dt`, reducing to `F = m a` when mass is constant.
- **Inertial frames**: a frame that is unaccelerated and non-rotating relative to a known inertial frame is itself inertial. A frame moving at constant velocity relative to an inertial frame is inertial (Galilean boost); a uniformly *accelerating* frame is not; a *rotating* frame is not. **[WEB]** — this is precisely what problems 1.26 (boost vs. accelerating frame) and 1.27/1.46 (rotating turntable) are constructed to demonstrate; the fetched stemjock solutions state that S′ moving at constant velocity is inertial "because it neither accelerates nor rotates," while S″ (constant acceleration) is not.

```
p         = m v                                  momentum
F         = dp/dt                                Newton II, general form
F         = m a                                  Newton II, constant mass
F = 0  =>  v = const                             Newton I

Galilean transformation between inertial frames S and S':
   r'  = r - V t          (S' moves at constant V relative to S)
   v'  = v - V
   a'  = a                 -> acceleration, hence F = m a, is Galilean invariant

Uniformly accelerating frame S'' (NOT inertial):
   r'' = r - (1/2) A t^2
   a'' = a - A             -> a fictitious force -mA appears
```

**[WEB]** for the boost/accelerating-frame results as they appear in the fetched solution to problem 1.26 (there written out in components: `x' = -Vt`, `x'' = -½At²` for a puck kicked along `y`).

---

### 1.5 The Third Law and Conservation of Momentum

Key ideas:
- **Third law (weak form)**: `F_12 = -F_21`. Equal and opposite.
- **Third law (strong form)**: additionally, the interaction force lies along the line joining the two particles (a *central* force). The strong form is what you need for angular-momentum conservation. **[WEB]** — the fetched statement of problem 1.34 explicitly says you need, in addition to the third law, the assumption that interparticle forces are central, and it forward-references Chapter 3 for the full angular-momentum discussion: https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p34.pdf
- **Momentum conservation** for an N-particle system is derived from the third law: internal forces cancel pairwise. Equations (1.25) through (1.29) are the steps of that derivation. **[WEB]** — the fetched statement of problem 1.28 names exactly that equation range as "the proof of conservation of momentum": https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p28.pdf
- **The third law is not universal.** Taylor shows (with a figure, Figure 1.8) that two charges moving with velocities perpendicular to each other exert magnetic forces that violate `F_12 = -F_21`. Momentum conservation is rescued only by giving the electromagnetic field its own momentum. The violation is of relative order `v1 v2 / c^2`, hence negligible in mechanics. **[WEB]** — problem 1.32's fetched statement references Figure 1.8, gives the E and B fields of a slowly moving charge, and asks the student to show `F_mag_12 <= (v1 v2 / c^2) F_el_12`: https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p32.pdf

```
F_12      = -F_21                                Newton III, weak form
F_12 parallel to (r_1 - r_2)                     Newton III, strong form (central)

P         = sum_alpha p_alpha                    total momentum of a system
dP/dt     = F_ext                                only EXTERNAL forces change P
F_ext = 0 =>  P = const                          conservation of momentum

Breakdown scale for Newton III (magnetic forces between moving charges):
   |F_mag_12| <= (v_1 v_2 / c^2) |F_el_12|
```

---

### 1.6 Newton's Second Law in Cartesian Coordinates

Key idea: in Cartesian coordinates and *only* in Cartesian coordinates, the vector equation `F = m r̈` resolves into three **independent** scalar second-order ODEs, because `x̂`, `ŷ`, `ẑ` are constant in time. This is the whole reason Cartesian coordinates are the natural internal representation for a simulator. The equations are only *uncoupled* if the force components happen not to mix the coordinates; they are always *separable into three scalar statements*.

```
m x_ddot = F_x(r, v, t)
m y_ddot = F_y(r, v, t)                          equation (1.35) in Taylor
m z_ddot = F_z(r, v, t)
```

**[WEB]** for the binding of equation number (1.35) to Newton II in Cartesian components used for projectile motion — the fetched statement of problem 1.35 (golf ball) instructs the student to use Newton's second law "(1.35)": https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p35.pdf

A second, subtler idea in this section: **you may tilt your axes**. For a body on an incline it is far easier to put one axis along the slope and one normal to it than to use horizontal/vertical. Gravity then has components `(-mg sin φ, -mg cos φ, 0)`. **[WEB]** — this is exactly the axis choice dictated by the fetched statements of problems 1.37, 1.38, and 1.39, and the fetched stemjock solution to 1.39 works out `W_x = -mg sin φ`, `W_y = -mg cos φ`: https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p39.pdf

---

### 1.7 Two-Dimensional Polar Coordinates

Key idea: in polar coordinates the basis vectors rotate with the particle, so `d(r̂)/dt ≠ 0`. Differentiating the position vector therefore produces "extra" terms that look like forces but are pure kinematics.

```
Transformations (Taylor's equation set (1.37)):
   x   = r cos(phi)
   y   = r sin(phi)
   r   = sqrt(x^2 + y^2)
   phi = arctan(y / x)      <-- INCOMPLETE; needs quadrant resolution (atan2)

Unit vectors (equation (1.59), stated inside problem 1.43):
   r_hat   =  x_hat cos(phi) + y_hat sin(phi)
   phi_hat = -x_hat sin(phi) + y_hat cos(phi)

Time derivatives of the rotating basis:
   d(r_hat)/dt   =  phi_dot * phi_hat        equation (1.42)
   d(phi_hat)/dt = -phi_dot * r_hat          equation (1.46)

Position, velocity, acceleration:
   r_vec = r * r_hat
   v     = r_dot * r_hat + r * phi_dot * phi_hat
   a     = (r_ddot - r phi_dot^2) r_hat + (r phi_ddot + 2 r_dot phi_dot) phi_hat

Newton II in polar coordinates - equation (1.48):
   F_r   = m (r_ddot - r phi_dot^2)          -r*phi_dot^2 is the centripetal term
   F_phi = m (r phi_ddot + 2 r_dot phi_dot)  2*r_dot*phi_dot is the Coriolis-like term
```

**[WEB]** for the polar velocity and acceleration formulas and for `d(r̂)/dt = φ̇ φ̂` — independently confirmed on a fetched third-party notes page: https://dandavison.github.io/taylor-classical-mechanics.html

**[WEB]** for the equation-number bindings: (1.37) = the four rectangular↔polar transformation equations, and the note that the `φ` equation is not quite complete (problem 1.42's fetched statement asks precisely this — it is the `atan2` quadrant problem); (1.42) and (1.46) = the two unit-vector time derivatives, and (1.59) = the `r̂ = x̂ cos φ + ŷ sin φ` expression (problem 1.43's fetched statement); (1.48) = Newton II in polar coordinates (problem 1.41's fetched statement).

The section ends by applying this machinery to a curved-track problem that yields a nonlinear pendulum equation, and then linearising it:

```
Equation (1.51), on page 31:
   phi_ddot = -(g/R) sin(phi)                 exact, nonlinear, NOT solvable in
                                              elementary functions

Small-angle limit, equation (1.55):
   phi_ddot = -omega^2 phi,   with omega = sqrt(g/R)

General solution, equation (1.56):
   phi(t) = A sin(omega t) + B cos(omega t)

Equation (1.57) is that solution specialised to release from rest at phi_0:
   phi(t) = phi_0 cos( sqrt(g/R) t )
   period T = 2 pi sqrt(R/g)                  independent of amplitude (small-angle only)
```

**[WEB]** for all of these: equation (1.51) with its page number and its exact form is quoted as an equation reference in the fetched stemjock solutions to problems 1.50 and 1.51; (1.55) and (1.56) are named in the fetched statement of problem 1.44; (1.57) is named in the fetched statement of problem 1.50. https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p44.pdf and .../TaylorCMCh1p50.pdf

---

## 2. Worked examples in Chapter 1

Taylor's Chapter 1 is unusually light on numbered worked examples — it is a setup chapter. I can confirm one by number and one by strong inference.

### Example 1.2 — the oscillating skateboard in a curved half-pipe (§1.7) — **[WEB]**

A skateboard (treated as a point mass) is released from rest on a frictionless semicircular trough of radius `R` and slides back and forth. Working in polar coordinates centred on the trough's centre of curvature: `r = R` is fixed, so `ṙ = r̈ = 0`, and the polar decomposition collapses to

```
radial:      -F_r_net = m ( -R phi_dot^2 )   ->  gives the normal force
tangential:  -m g sin(phi) = m R phi_ddot
             =>  phi_ddot = -(g/R) sin(phi)          [Taylor eq. (1.51)]
```

Physics illustrated: (a) how the rotating polar basis turns a constrained 2-D problem into a *single* scalar ODE in `φ`; (b) that the exact equation is nonlinear and has no elementary closed-form solution — Taylor says so outright, which is the book's first explicit argument for numerical methods; (c) small-angle linearisation to simple harmonic motion with period `2π sqrt(R/g)`; (d) that the "pendulum" here is a constraint problem, not a string-tension problem, but has identical mathematics.

**[WEB]** that this is *Example 1.2* and that it is the source of equation (1.51): the fetched statement of problem 1.50 names "the skateboard of Example 1.2" and its equation (1.51). **[WEB]** for the numbers `R = 5 m`, `g = 9.8 m/s²`, from the fetched statements of 1.50/1.51 and from the fetched dandavison notes page.

### Example 1.1 — a block sliding down an incline (§1.6) — **[MEM]**

A block slides on an inclined plane; the point is the choice of tilted axes (one along the slope, one normal to it), the resolution of gravity into `mg sin θ` down-slope and `mg cos θ` into the surface, the normal force from `a_normal = 0`, and (with friction) `f = μ N`. Physics illustrated: free-body diagrams, non-obvious axis choice, and constraint forces determined by a kinematic condition rather than a force law.

**Honesty note:** a search-engine summary asserted that Example 1.1 is titled with a block sliding down an incline and sits in §1.6, but **I never fetched a page that says so**, so I am tagging this [MEM]. What *is* [WEB] is that the incline-with-tilted-axes technique is the backbone of problems 1.37, 1.38 and 1.39, which is only sensible if the text worked such an example first.

### Other set-piece derivations in the running text (not numbered examples)

- The N-particle momentum-conservation derivation, equations (1.25)–(1.29). **[WEB]** (named in problem 1.28's fetched statement.)
- The two-moving-charges counterexample to the third law, Figure 1.8. **[WEB]** (named in problem 1.32's fetched statement.)

**Not present in Chapter 1, contrary to the local map:** I found no Atwood's machine and no conical pendulum anywhere in Chapter 1's problems — see §3.6.

---

## 3. Problem archetypes, 1.1–1.51

I downloaded and text-extracted all 51 Chapter-1 solution PDFs from stemjock.com, so **every problem statement below is [WEB]-verified** (base URL: `https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p<N>.pdf`). All descriptions are my own paraphrases.

The numbering tracks the sections: **1.1–1.25 ≈ §1.2 (vectors/calculus)**, **1.26–1.34 ≈ §1.4–1.5 (frames, third law, momentum)**, **1.35–1.40 ≈ §1.6 (Cartesian Newton II)**, **1.41–1.51 ≈ §1.7 (polar/cylindrical)**.

### 3.1 Vector algebra and vector products (1.1–1.23) — 23 problems, ~45% of the chapter

- **1.1, 1.2** — Given two concrete vectors, compute sum, a linear combination, dot product, cross product. Pure `Vec3` arithmetic drill.
- **1.3** — Derive `r² = x² + y² + z²` by applying Pythagoras twice.
- **1.4** — Extract the angle between two given vectors from their dot product.
- **1.5** — Angle between a cube's body diagonal and a face diagonal (the classic `arccos(sqrt(2/3)) ≈ 35.26°` construction).
- **1.6** — Find the scalar parameter making two parameterised vectors orthogonal.
- **1.7** — Prove the geometric and component definitions of the dot product agree.
- **1.8** — (a) dot product is distributive; (b) the product rule holds for `d/dt (r·s)`.
- **1.9** — Law of cosines falls straight out of expanding `|a+b|²`.
- **1.13** — Prove `b² = (û·b)² + (û×b)²` (parallel/perpendicular decomposition).
- **1.14** — Triangle inequality `|a+b| ≤ |a| + |b|`.
- **1.15** — Prove the component cross product equals the geometric (perpendicular, `rs sin θ`, right-hand-rule) definition.
- **1.16** — (a) `|r| = sqrt(r·r)`; (b) the dot product is invariant under change of orthogonal axes.
- **1.17** — (a) cross product is distributive; (b) product rule for `d/dt (r×s)`, minding operand order.
- **1.18** — Triangle area as `½|a×b|` in three equivalent ways; deduce the law of sines.
- **1.19** — Show `d/dt [a·(v×r)] = ȧ·(v×r)` (the other two terms vanish because of repeated vectors).
- **1.20** — Area of a triangle from three position vectors: `½|B×C + C×A + A×B|`.
- **1.21** — Parallelepiped volume as the scalar triple product `|a·(b×c)|`.
- **1.22** — Derive the cosine and sine addition formulas by evaluating dot and cross products two ways.
- **1.23** — Invert the system `b·v = λ`, `b×v = c` for the unknown vector `v`.

**Simulator relevance:** 1.1–1.9 and 1.13–1.23 are a *ready-made unit-test suite for your `Vec2`/`Vec3` class*. Every one is a closed-form identity you can assert to machine epsilon. 1.17(b) in particular catches the classic bug of writing `d/dt(r×s)` with the operands swapped.

### 3.2 Kinematics from a given trajectory (1.10–1.12) — no dynamics needed

- **1.10** — Show that `r(t) = x̂ R cos(ωt) + ŷ R sin(ωt)` describes uniform circular motion; differentiate to get `v` and `a`; identify `|v| = Rω`, `|a| = Rω²` directed toward the centre.
- **1.11** — Identify the orbit of `r(t) = x̂ b cos(ωt) + ŷ c sin(ωt)` as an ellipse `x²/b² + y²/c² = 1`.
- **1.12** — Same plus a `ẑ v₀ t` term: an elliptical helix.

**Simulator relevance:** these are exact analytic trajectories with no force model at all. They are the cheapest possible tests of your differentiation/interpolation and of any dense-output machinery.

### 3.3 Differential-equation warm-ups (1.24, 1.25, 1.44)

- **1.24** — Solve `df/dt = f` by separation; observe that an nth-order ODE's general solution carries n arbitrary constants.
- **1.25** — Same for `df/dt = -3f`.
- **1.44** — Verify by direct substitution that `A sin(ωt) + B cos(ωt)` solves `φ̈ = -ω²φ`.

Taylor is explicitly teaching ODEs from scratch here — he says he will introduce the necessary ideas as they are needed. This matters for you: **Chapter 1 assumes zero prior ODE knowledge**, so anything your library does beyond "integrate `ẏ = f(y,t)`" is ahead of the course.

### 3.4 Inertial frames, Galilean transformation, third law, momentum (1.26–1.34)

- **1.26** — A frictionless puck is kicked across a floor. Describe its coordinates in the original inertial frame, in a frame boosted at constant velocity, and in a frame under constant acceleration. Which frames are inertial? *(Straight Galilean-transformation drill; the accelerating frame gives a parabola.)*
- **1.27** — The same puck slid straight across a **rotating turntable**, passing through the centre. Describe the path a turntable-riding observer sees. Qualitative; the stemjock solution notes the puck's apparent speed is `sqrt(v₀² + R²ω²)` at the rim and `v₀` at the centre.
- **1.28, 1.29** — Write out the momentum-conservation derivation explicitly for N = 3 and N = 4 with all sums expanded.
- **1.30** — Perfectly inelastic collision: mass `m₁` at velocity `v` strikes stationary `m₂` and they stick. `v' = m₁ v / (m₁ + m₂)`.
- **1.31** — The converse of the text's result: if momentum is conserved for every subgroup of particles, the interparticle forces must obey the third law.
- **1.32** — Electromagnetism: show the magnetic force between two slowly moving charges is smaller than the electric force by `v₁v₂/c²`, so the third law is *approximately* valid at low speed.
- **1.33** — Show the magnetic forces between two *steady current loops* do obey the third law exactly (Biot–Savart plus the BAC−CAB rule).
- **1.34** — Prove total angular momentum `L = Σ r_α × p_α` is conserved for an isolated system, *given* central internal forces. Taylor forward-references Chapter 3.

### 3.5 Newton II in Cartesian components — projectiles and inclines (1.35–1.40)

All six are vacuum (drag-free) projectile problems that differ in the choice of axes.

- **1.35** — Golf ball hit at speed `v₀`, angle `θ`, in 3-D axes (east/north/up). Find `r(t)`, time of flight, range. Results: `T = 2v₀ sin θ / g`, range `= v₀² sin(2θ)/g`.
- **1.36** — A plane flying horizontally at height `h`, speed `v₀`, drops a bundle to a raft. Find `r(t)`, the horizontal lead distance `v₀ sqrt(2h/g)`, and the timing tolerance `Δt = Δx / v₀`. With `v₀ = 50 m/s`, `h = 100 m`, `g ≈ 10 m/s²`: lead ≈ 224 m, tolerance ≈ ±0.2 s for ±10 m.
- **1.37** — Frictionless puck kicked straight up an incline of angle `θ`. Find `x(t)` and the return time (`2v₀ / (g sin θ)`).
- **1.38** — Puck kicked across a **tilted board** with a general 2-component initial velocity in the board's plane. Axes: along the bottom edge, up the slope, normal to the board. Find the time to reach floor level and the distance from the origin.
- **1.39** — Ball thrown up an inclined plane at angle `θ` above the *slope*, the slope at `φ` above horizontal. Show `R = 2v₀² sin θ cos(θ+φ) / (g cos²φ)` and `R_max = v₀² / [g(1 + sin φ)]` (attained at `2θ + φ = π/2`).
- **1.40** — Cannonball at angle `θ`. Find the largest `θ` for which the *distance from the launch point* `r(t)` increases monotonically throughout the flight. Answer: `θ_max = arcsin(sqrt(8/9)) ≈ 1.23 rad ≈ 70.5°`.

**All the numeric answers above are [WEB]** — read out of the fetched stemjock solution PDFs for 1.35, 1.36, 1.37, 1.39, 1.40.

A recurring student trap, worth encoding as a comment in your examples: in 1.39 there is **no normal force** during flight, because the ball is in the air, not on the plane — the tilted axes are a bookkeeping convenience only. **[WEB]** — this is the entire subject of a fetched Physics Forums thread on the problem: https://www.physicsforums.com/threads/classical-mechanics-taylor-1-39.432209/

### 3.6 Polar and cylindrical coordinates (1.41–1.51)

- **1.41** — An astronaut in gravity-free space twirls a mass on a string of length `R` at constant `ω`. Write Newton II in polar form (eq. 1.48) and get the tension: `T = m R ω²`. The `φ` equation is identically satisfied.
- **1.42** — Prove the four rectangular↔polar transformations (1.37) and fix the incomplete `arctan` formula for `φ`. *This is the `atan2` problem.*
- **1.43** — Derive `r̂ = x̂ cos φ + ŷ sin φ` and the matching `φ̂`; differentiate to re-derive equations (1.42) and (1.46).
- **1.44** — Verify the SHM solution satisfies `φ̈ = -ω²φ`.
- **1.45** — Prove: a time-dependent vector has constant magnitude **iff** its derivative is orthogonal to it. Taylor points out two payoffs: it explains why `d(r̂)/dt` must point along `φ̂`, and it shows a charged particle's speed is constant in a magnetic field.
- **1.46** — The quantitative version of 1.27: the turntable puck's `(r, φ)` in the ground frame and `(r', φ')` in the rotating frame. The rotating-frame path is a spiral in and back out. Is the rotating frame inertial? (No.)
- **1.47** — Define **cylindrical polar coordinates** `(ρ, φ, z)`; express them in Cartesian terms; note the notational hazard of writing `r` for `ρ`.
- **1.48** — Express `ρ̂`, `φ̂`, `ẑ` in Cartesian terms and differentiate them with respect to time.
- **1.49** — A frictionless puck confined between two nearly-coincident concentric vertical cylinders of radius `R`, under gravity. Solve in cylindrical polars. Result: `ρ = R`, `φ(t) = ωt + φ₀`, `z(t) = -½gt² + v₀z t + z₀` — an **elongating helix** — with normal force `N = m R ω²`.
- **1.50 — [Computer]** — Numerically solve the skateboard equation (1.51) with `φ₀ = 20°`, `R = 5 m`, `g = 9.8 m/s²`; plot `φ(t)` over two or three periods; overlay the small-angle solution (1.57) at the same `φ₀`; comment. The ODE reduces to `φ̈ = -1.96 sin φ`, `φ(0) = π/9`, `φ̇(0) = 0`. Taylor's stated framing is that the ODE has no elementary closed-form solution but is easy to solve numerically; the problem names Mathematica/Maple/Matlab as example tools.
- **1.51 — [Computer]** — Repeat 1.50 with `φ₀ = π/2`. The published solution's verdict: the two curves do not overlap, so the small-angle approximation fails badly at 90°.

### The [Computer] problems in Chapter 1 — **exactly two: 1.50 and 1.51** — **[WEB]**

I verified this exhaustively, not by sampling: I text-extracted all 51 statements and searched for Taylor's `[Computer]` tag. It appears in **1.50 and 1.51 only**. (A naive full-document grep also hits every file, but only inside the embedded "Computer Modern" font-licence strings — a false positive I checked and excluded.)

This is the single most important fact in this report for your library. **Chapter 1 asks the student to do exactly one numerical thing:** integrate a second-order scalar ODE `φ̈ = -(g/R) sin φ` from rest at two different amplitudes, and plot it against the linearised solution. If your library can do that cleanly, it has fully discharged Chapter 1.

Small-angle reference values for those two problems: `ω = sqrt(9.8/5) = 1.4 rad/s` exactly, so the linearised period is `2π/1.4 ≈ 4.488 s`, and 1.50's plot window of 0–14 s covers just over three periods (1.51's is 0–21 s). **[WEB]** from the fetched 1.50 and 1.51 solutions.

### 3.7 Corrections to `docs/research/taylor-ch1-7-map.md`

Two errors and several omissions in the existing map's Chapter 1 section:

1. **"Atwood's machine" is listed as a Chapter 1 canonical example. It is not in Chapter 1's problems.** I searched all 51 extracted statements and solutions for "Atwood" and "pulley" — **zero hits** **[WEB]**. My recollection is that Taylor's Atwood machine is an *Example in Chapter 7* (the Lagrangian showcase), which the map itself says in its own Chapter 7 section. Whether an Atwood machine appears in Chapter 1's *running text* I could not confirm either way; the evidence says it does not.
2. **"Conical pendulum / uniform circular motion" is listed as a Chapter 1 example.** Uniform circular motion is certainly there (problems 1.10, 1.41). A *conical* pendulum is not — problem 1.41's rotating mass is explicitly in gravity-free space, so there is no cone. **[WEB]**
3. **Omission: cylindrical polar coordinates.** Problems 1.47–1.49 introduce `(ρ, φ, z)` and use them for real dynamics. The map defers cylindrical coordinates to Chapter 4. Chapter 1 already needs them. **[WEB]**
4. **Omission: the Galilean transformation and non-inertial frames.** Problems 1.26, 1.27, 1.46 are a coherent three-problem thread that the map's Chapter 1 section doesn't mention, but which maps directly to a `Frame`/`observe_from` feature.
5. **Omission: the vector-identity problems as a test suite.** 1.13–1.23 are nineteen free unit tests for your vector type.
6. The map says "this is the chapter where [all the integrators] should appear." That is a design choice, not something Chapter 1 demands. Chapter 1 justifies exactly one integrator.

---

## 4. Simulation-relevant ideas from Chapter 1 **only**

What follows is what a force-based particle simulator can legitimately do using nothing beyond Chapter 1 physics. Anything requiring a later chapter is flagged in red-letter terms so you don't build it by accident and then find you can't validate it.

### 4.1 In scope with Chapter 1 alone

**(a) The core state and the core loop.**
```
state y = (r, v)  in R^6
dy/dt = f(y, t) = ( v , F(r, v, t) / m )
```
Chapter 1 gives you `F = m a` and nothing else; that is enough for the whole reduction. §1.6's point — that Cartesian components separate cleanly — is the justification for storing state in Cartesian coordinates internally, always, and converting to polar only for reporting.

**(b) Constant forces.** Uniform gravity `F = (0, 0, -mg)`. Every one of problems 1.35–1.40 is this force with different axes and initial conditions. Fully validated by closed form (§5).

**(c) Tilted / rotated coordinate frames.** Problems 1.37, 1.38, 1.39 all reduce to "rotate the gravity vector, then run the constant-force integrator." A `Frame` type with a rotation and an origin, plus `to_world()` / `from_world()`, is directly motivated by Chapter 1 and lets you reproduce 1.37–1.39 without new physics.

**(d) Two-body Newton III and momentum conservation.** A pairwise force interface `Vec3 f(const Particle& a, const Particle& b)` that is *by construction* applied as `+f` to `a` and `-f` to `b` makes the weak third law structurally true rather than a runtime hope. Add a `is_central()` trait (force parallel to `r_a - r_b`) because problem 1.34 shows the strong form is what angular momentum needs. A `total_momentum()` diagnostic that must be constant to roundoff under internal-only forces is your best invariant test in the whole chapter.

Caution worth encoding: Chapter 1 also shows the third law *fails* for magnetic interactions between moving charges (problem 1.32). If you ever add a velocity-dependent pairwise force, your "momentum is conserved" assertion is no longer guaranteed. Document that boundary now.

**(e) Inertial vs. non-inertial observation — the Galilean transformation.** Problems 1.26, 1.27, 1.46 are *pure post-processing*. You simulate once in the inertial frame and then re-express the trajectory:
```
boost:      r'(t)  = r(t) - V t
accelerate: r''(t) = r(t) - (1/2) A t^2
rotate:     r'(t)  = R(-omega t) r(t)      i.e. phi' = phi - omega t,  r' = r
```
This is the *right* Chapter-1 design: an `observe_from(Frame)` view over a recorded trajectory, **not** fictitious forces added to the equation of motion. Chapter 1 never adds a centrifugal or Coriolis *force*; it only asks what the moving observer sees.

**(f) Polar decomposition of velocity and acceleration.** Implement
```
polar_kinematics(r, r_dot, r_ddot, phi, phi_dot, phi_ddot) -> (a_r, a_phi)
```
and its Cartesian inverse. Then write a test that runs a Cartesian simulation and verifies numerically that `a_r == r̈ - r φ̇²` and `a_φ == r φ̈ + 2 ṙ φ̇` along the trajectory. This is both a correctness test and the single best teaching artifact you can build from §1.7. Include a `to_polar` that uses `atan2`, and cite problem 1.42 in the comment — Taylor's own problem is about exactly that bug.

**(g) Uniform circular motion.** Problems 1.10, 1.41, 1.49 give you exact circular and helical solutions with `|v| = Rω`, `|a| = Rω²`, `T = mRω²`, `N = mRω²`. Perfect closed-form targets.

**(h) The half-pipe / pendulum problem in polar form — the chapter's flagship numeric.** Support integrating a *constrained scalar* second-order ODE directly:
```
phi_ddot = -(g/R) sin(phi),  phi(0) = phi_0,  phi_dot(0) = 0
```
Two ways to build it, both defensible with Chapter 1 physics: (i) a dedicated 1-DOF `AngularState {phi, phi_dot}` integrator, which is what Taylor literally asks for; or (ii) a 2-D Cartesian simulation with a rigid-radius constraint whose normal force you compute from the radial equation `N = m(g cos φ + R φ̇²)`. Option (i) is honest to Chapter 1; option (ii) is a stronger test because it exercises your vector path, and Chapter 1 gives you the radial equation you need to make it work. Doing both and comparing them is a genuinely good exercise.

**(i) Constraint problems of the "fixed radius" type.** Problem 1.49 (puck between coaxial cylinders) is the cleanest one: `ρ` pinned, `φ` and `z` free, gravity along `-z`, giving a helix that stretches. It needs cylindrical coordinates and one algebraic normal force. It is a good second constraint test after the half-pipe because the answer is closed-form.

**(j) Atwood-type constraint problems — buildable, but not from a Chapter 1 problem.** The physics an Atwood machine needs is entirely Chapter 1 (Newton II on each mass, third law for the string tension, and the constraint `a₁ = -a₂`), and the result `a = (m₁ - m₂)g/(m₁ + m₂)` is a good closed-form validation target. But **Taylor does not pose an Atwood machine in Chapter 1** — see §3.6. So implement it if you want the test; just don't tell yourself you're following the book here.

**(k) A `Vec2`/`Vec3` type with a real test suite.** Problems 1.1–1.9 and 1.13–1.23 hand you roughly twenty exact identities. Use them.

**(l) Trajectory recording and plotting.** Both [Computer] problems ask for a *plot* of `φ(t)` over two to three periods with a second analytic curve overlaid. A `Recorder` that samples at fixed output intervals (decoupled from the integrator's internal step) plus CSV output is the minimum deliverable that lets you reproduce 1.50 and 1.51 exactly.

### 4.2 Natural extensions that REQUIRE later chapters — do not build these yet

| Tempting feature | Chapter it actually needs | Why you can't validate it from Ch. 1 |
|---|---|---|
| **Air resistance / drag**, linear `-bv` or quadratic `-cv²` | **Chapter 2** ("Projectiles and Charged Particles") **[MEM]**, though a fetched UCSD homework set confirms Taylor 2.14 is a velocity-dependent-force ODE problem **[WEB]** (https://courses.physics.ucsd.edu/2010/Fall/physics110a/HOMEWORK/hwsolutions1_F10.pdf) | Ch. 1 has no velocity-dependent force at all. Every Ch. 1 projectile is explicitly stated to neglect air resistance. You'd have no analytic reference and no terminal-velocity concept. |
| **Lorentz force / charged particle in a B field** | **Chapter 2** **[MEM]** | Ch. 1 touches `E` and `B` only in problems 1.32/1.33, and only to discuss the third law — not to integrate motion. (Problem 1.45's remark that a magnetic field leaves speed constant is a nice hint, but that's all it is.) |
| **Rocket / variable-mass propulsion** | **Chapter 3** **[WEB]** — the fetched UCSD set works Taylor 3.11 and 3.13, both rocket problems | `F = dp/dt` with `ṁ ≠ 0` is beyond §1.4's constant-mass reduction. |
| **Angular momentum as a first-class diagnostic; centre of mass** | **Chapter 3** **[WEB]** — problem 1.34's own statement forward-references Ch. 3 | Ch. 1 gives you `L` conservation as a one-off exercise under the strong third law, not the machinery. |
| **Work, kinetic/potential energy, energy conservation as an integrator quality metric** | **Chapter 4** **[MEM]**, corroborated by the fetched UCSD set working Taylor 4.2 / 4.23 / 4.36 **[WEB]** | This is the big temptation. "Symplectic integrators conserve energy" is the standard motivation for Verlet — **but Chapter 1 has no concept of energy**, so you cannot state that motivation yet. Use *momentum* conservation and *exactness on constant force* as your Ch. 1 integrator metrics instead. |
| **Damped and driven oscillators, resonance, Q factor** | **Chapter 5** **[MEM]** | Ch. 1 gives you undamped SHM only as a small-angle *approximation* to one specific problem. Don't generalise it into an `Oscillator` module yet. |
| **Centrifugal and Coriolis forces as forces in a rotating frame** | **Chapter 9** ("Mechanics in Noninertial Frames") **[MEM]**, corroborated by an independent fetched lecture-notes PDF deriving `m r̈ = F - m ω̇×r - m ω×(ω×r) - 2m ω×ṙ` for rotating frames **[WEB]** | Ch. 1's `2 ṙ φ̇` term is *kinematic*, not a force. Problems 1.27 and 1.46 deliberately ask only what the rotating observer *sees*. Adding fictitious forces to your integrator now would be jumping seven chapters. |
| **Lagrangian / generalised coordinates / constraint solvers** | **Chapters 6–7** **[MEM]** | The half-pipe is a constraint problem solved *by hand* in Ch. 1. A general constraint engine is Ch. 7 territory. |
| **Central-force orbits, Kepler** | **Chapter 8** **[MEM]** | Ch. 1's ellipse (problem 1.11) is a *prescribed kinematic* ellipse, not a dynamical orbit. |

**One design decision Chapter 1 forces on you early.** Velocity Verlet and the symplectic family assume `F = F(r)` only. Every Chapter 1 force satisfies that. Chapter 2's drag and Chapter 5's damping do not. So add the force-trait flag `depends_on_velocity()` *now*, while everything returns `false`, and make the integrator selection check it — rather than discovering the problem when you hit Chapter 2. The existing map already makes this point and it is correct.

---

## 5. Analytic closed-form results from Chapter 1 alone, for validating numerics

Every entry is derivable from Chapter 1 physics. Confirmed values are marked [WEB] where I read them out of a fetched solution PDF.

### 5.1 Constant force — the sharpest first test

```
r(t) = r0 + v0 t + (1/2) (F/m) t^2
v(t) = v0 + (F/m) t
```
A degree-2 polynomial trajectory. **RK4, RK2, velocity Verlet and leapfrog must all reproduce this to roundoff for any step size**, while explicit Euler must not. That one test discriminates a correct integrator implementation from a broken one before you write anything else. **[MEM]** (elementary; no source needed).

### 5.2 Vacuum projectile, flat ground, launch from origin at speed `v0`, angle `θ`

```
x(t)   = v0 cos(theta) t
y(t)   = v0 sin(theta) t - (1/2) g t^2
T      = 2 v0 sin(theta) / g            time of flight        [WEB, problem 1.35]
R      = v0^2 sin(2 theta) / g          range                 [WEB, problem 1.35]
h_apex = v0^2 sin^2(theta) / (2 g)      apex height           [MEM]
t_apex = v0 sin(theta) / g                                    [MEM]
R is maximised at theta = 45 degrees                          [MEM]
```

### 5.3 Vacuum projectile, horizontal launch from height `h` (problem 1.36) — **[WEB]**

```
t_fall = sqrt(2 h / g)
x_lead = v0 sqrt(2 h / g)
Delta_t tolerance for a horizontal miss of Delta_x:  Delta_t = Delta_x / v0

Numeric check: v0 = 50 m/s, h = 100 m, g = 10 m/s^2
   x_lead  = 100 sqrt(5) ~ 224 m
   +/-10 m  <->  +/-0.2 s
```

### 5.4 Projectile up an inclined plane (problem 1.39) — **[WEB]**

Slope at `φ` above horizontal, launch at `θ` above the slope:
```
R     = 2 v0^2 sin(theta) cos(theta + phi) / (g cos^2(phi))
R_max = v0^2 / [ g (1 + sin(phi)) ]        attained at 2 theta + phi = pi/2
```
Excellent test of a rotated-frame implementation, because setting `φ = 0` must collapse it to `v0² sin(2θ)/g`.

### 5.5 Puck up a frictionless incline (problem 1.37) — **[WEB]**

```
x(t)      = v0 t - (1/2) g sin(theta) t^2
t_return  = 2 v0 / (g sin(theta))
```

### 5.6 Monotonic-range condition (problem 1.40) — **[WEB]**

```
theta_max = arcsin( sqrt(8/9) ) ~ 1.23 rad ~ 70.5 degrees
```
The distance `|r(t)|` from the launch point increases throughout the flight iff `θ ≤ θ_max`. A lovely non-obvious test: your simulator should find `d|r|/dt` dips negative just above 70.5° and never below it.

### 5.7 Uniform circular motion (problems 1.10, 1.41) — **[WEB]**

```
r(t)   = ( R cos(omega t), R sin(omega t) )
|v|    = R omega                          constant speed
|a|    = R omega^2 = |v|^2 / R            directed toward the centre
T_period = 2 pi / omega = 2 pi R / |v|
Tension (or any central force holding the circle):  F = m R omega^2 = m v^2 / R
```
Also gives you two *invariant* assertions to run every step: `|r| = R` and `a · v = 0` (the latter is exactly problem 1.45's theorem).

### 5.8 Constrained helix — puck between coaxial cylinders (problem 1.49) — **[WEB]**

```
rho(t) = R                      (constant)
phi(t) = omega t + phi_0
z(t)   = -(1/2) g t^2 + v0z t + z0
N      = m R omega^2
```
This one is valuable because it is a *three-dimensional constrained* motion with a fully closed-form answer — rare, and ideal for testing a cylindrical-coordinate path.

### 5.9 Two-body momentum conservation

```
Internal forces only:  P(t) = P(0) exactly, for all t
Perfectly inelastic collision (problem 1.30):  v' = m1 v / (m1 + m2)     [WEB]
```
For a pairwise-symmetric force implementation, `P` should be conserved to *roundoff*, not to truncation error — because the cancellation is structural. If your momentum drift grows with step size, your force application is asymmetric. That is a much sharper test than energy drift, and it is the one Chapter 1 licenses you to use.

### 5.10 Atwood's machine

```
a = (m1 - m2) g / (m1 + m2)
T = 2 m1 m2 g / (m1 + m2)
```
Derivable from Chapter 1 (Newton II per mass + third law for the tension + the constraint `a₁ = -a₂`), and a good closed-form target — but **not** a Chapter 1 problem in Taylor. **[MEM]** for the formulas (elementary); **[WEB]** for the absence from Chapter 1.

### 5.11 Small-angle half-pipe / pendulum period (Example 1.2, problems 1.50–1.51) — **[WEB]**

```
Exact:        phi_ddot = -(g/R) sin(phi)          no elementary closed form
Small angle:  phi(t) = phi_0 cos( sqrt(g/R) t )
              T = 2 pi sqrt(R / g)                amplitude-independent

For Taylor's numbers R = 5 m, g = 9.8 m/s^2:
   omega = sqrt(9.8/5) = 1.4 rad/s exactly
   T     = 2 pi / 1.4 ~ 4.488 s
   phi_ddot = -1.96 sin(phi)
```
**How to use this as a graded numerical test**, mirroring what 1.50 and 1.51 actually ask:
- At `φ₀ = 20°`, the numerical solution and `φ₀ cos(sqrt(g/R) t)` should track closely over 2–3 periods; the true period is larger than `2π sqrt(R/g)` by roughly `1 + φ₀²/16 ≈ 0.76%`, i.e. about 4.522 s versus 4.488 s. (The `φ₀²/16` correction is **[MEM]**, standard pendulum theory, *not* in Taylor Chapter 1 — treat it as a tolerance guide, not as something the book gives you.)
- At `φ₀ = π/2`, the two curves visibly diverge — the fetched 1.51 solution's own conclusion is that the small-angle approximation is not good there. That divergence is the *expected* result, so your test should assert it, not treat it as failure.
- Bonus invariant Chapter 1 does **not** give you: the pendulum's energy. Resist using it. Use instead the fact that the motion is exactly periodic and symmetric — `φ(t)` must return to `+φ₀` with `φ̇ = 0`, and the turning points must not drift. Turning-point drift over many periods is a legitimate, energy-free integrator-quality metric available from Chapter 1 alone.

---

## Sources

Pages and files I actually fetched this session:

- Solutions index, Taylor *Classical Mechanics* — https://stemjock.com/taylorcm.htm (confirmed problems 1.1–1.51 exist)
- All 51 Chapter 1 solution PDFs, `https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p1.pdf` … `TaylorCMCh1p51.pdf` — downloaded and text-extracted locally; the primary evidence for every problem statement, every equation-number binding, and every closed-form answer in §3 and §5. Local copies: `/private/tmp/claude-501/-Users-jackery-Github-newtonian-library/e1090fdc-4da5-488d-aae7-cac525791aa8/scratchpad/ch1/`, with extracted statements in `.../scratchpad/ch1_statements.txt`
- Dan Davison, study notes on Taylor's *Classical Mechanics* — https://dandavison.github.io/taylor-classical-mechanics.html (polar velocity/acceleration, `d(r̂)/dt = φ̇ φ̂`, half-pipe with `R = 5 m` and `mRφ̈ = -mg sin φ`)
- StudySoup, Classical Mechanics Chapter 1 listing — https://studysoup.com/tsg/physics/90/classical-mechanics/chapter/1765/1 (independent confirmation of 51 problems and of the 1.1–1.23 topics)
- U.S. Naval Academy SP333 course page — https://www.usna.edu/Users/physics/cmorgan/SP333.php (Chapter 1 covered as sections 1–7 in one week; assigned problems 1.6, 1.12, 1.18, 1.25, 1.27, 1.32, 1.41, 1.44, 1.47)
- UCSD Physics 110A, HW 1 solutions — https://courses.physics.ucsd.edu/2010/Fall/physics110a/HOMEWORK/hwsolutions1_F10.pdf (Taylor 1.46 assigned; also Taylor 2.14, 3.11, 3.13, 4.2, 4.23, 4.36, confirming which topics live in which later chapters)
- Rhett Allain, Classical Mechanics course materials — https://rhettallain.com/classical-mechanics-course-material-and-videos/ (Chapter 1 taught as Newton II + polar coordinates, with Python/Jupyter numerical solutions)
- Physics Forums thread on Taylor 1.39 — https://www.physicsforums.com/threads/classical-mechanics-taylor-1-39.432209/ (the no-normal-force-in-flight confusion)
- Joachim Brod, *Classical Mechanics* lecture notes, University of Cincinnati — https://homepages.uc.edu/~brodjm/teaching/lecture_notes/mechanics.pdf (independent corroboration of inertial frames / Galilean invariance, and of the rotating-frame fictitious-force equation that belongs to Taylor Ch. 9, not Ch. 1)
- Local starting hypothesis, read but independently checked — `/Users/jackery/Github/newtonian-library/docs/research/taylor-ch1-7-map.md`

Sources that appeared in search results and that I deliberately did **not** open, because they are unauthorised full-text copies of a copyrighted book: the archive.org full-text dump, pdfcoffee, silo.tips, slideshare, studylib, scribd, zlib.pub, sweetstudy, and the Internet Archive copy of the instructor's solution manual.

Attempted but unavailable: `henry.pha.jhu.edu/chapters.html` (HTTP 403), the MIT Press book page (HTTP 403), and the University Science Books pages including the preface PDF (302-redirected to a publisher landing page after the July 2025 transfer to MIT Press).

---

## Could not confirm — check these against your physical copy

1. **The exact titles of sections 1.1–1.7 and their page numbers.** Search summaries agreed with the list in the task prompt (and gave 1.2 → p. 4, 1.3 → p. 9, 1.5 → p. 17, 1.7 → p. 26, summary page → p. 33, problems → p. 34), but **no page I fetched displayed them**. Only "equation (1.51) is on page 31" is [WEB]. Verify the TOC.
2. **That Example 1.1 is the block on an incline, and that it sits in §1.6.** [MEM]; a search summary agreed, but I never fetched a page confirming it. Also verify the block's details — whether Taylor's version includes friction, and whether he solves for the normal force.
3. **How many numbered examples Chapter 1 contains.** I confirmed Example 1.2 exists [WEB]. I do not know whether the chapter stops at 1.2 or continues to 1.3+. Check.
4. **Whether an Atwood's machine appears anywhere in Chapter 1's running text.** It is definitively absent from the 51 problems [WEB]. If it is in the text, tell me and the library's Atwood test gets a legitimate Chapter 1 citation.
5. **Whether the chapter contains a conical pendulum.** Absent from the problems [WEB]; unknown for the text.
6. **Taylor's difficulty-star convention** (single/double/triple star markings on problems). I know Taylor uses one [MEM] but I confirmed nothing about which Chapter 1 problems carry which rating — the stemjock reproductions do not preserve the stars. If you want to grade your examples by difficulty, read the stars off the book.
7. **Whether Taylor's `[Computer]` convention in Chapter 1 extends to any problem besides 1.50 and 1.51.** My evidence is strong (all 51 statements scanned, [WEB]) but it depends on stemjock having transcribed the tag faithfully everywhere, not just twice. A 30-second flip through pages 34–41 will settle it.
8. **The precise wording of equation (1.35).** I confirmed the *number* is Newton II in Cartesian form used for projectiles [WEB], but not whether Taylor writes it as three scalar equations, as a single vector equation, or both.
9. **Whether §1.2 explicitly presents the BAC−CAB triple-product identity, or leaves it to the problems.** Problem 1.33 assumes the reader has it [WEB]; where it is introduced, I don't know.
10. **The `φ₀²/16` amplitude correction to the pendulum period** used as a tolerance guide in §5.11 is standard physics [MEM] but is **not** Chapter 1 material. Don't cite it as coming from Taylor Ch. 1.