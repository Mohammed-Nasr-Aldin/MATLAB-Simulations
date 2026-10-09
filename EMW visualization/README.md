<h1 align="center">Electromagnetic Waves (EMW) · MATLAB Visualization</h1>

<p align="center">
  <b>Visualization and Analysis of Vector Fields and Plane Electromagnetic Waves</b><br>
  Divergence, Curl and Gradient of the Electric and Magnetic Fields: A MATLAB Study
</p>

<p align="center">
  <img alt="Tool" src="https://img.shields.io/badge/Tool-MATLAB%20R2024a-2A55E5?style=flat-square">
  <img alt="Course" src="https://img.shields.io/badge/Course-Electromagnetic%20Waves-2A55E5?style=flat-square">
  <img alt="Type" src="https://img.shields.io/badge/Type-Report-C8102E?style=flat-square">
</p>

<p align="center">
  <a href="report/EMW-Matlab_visualization.pdf">
    <img alt="Open the report" src="https://img.shields.io/badge/REPORT-Open%20the%20PDF%20%E2%86%92-2A55E5?style=for-the-badge&labelColor=0B1F4D">
  </a>
</p>

---

## 📖 Overview

This repository contains my report for the **Electromagnetic Waves (EMW)** course. The report uses MATLAB to visualize and analyze the three differential operators that Maxwell's equations are built on: **divergence**, **curl** and **gradient**.

Each field is first drawn in 3D, then on the three coordinate planes (XY, XZ and YZ), so you can see what each operator does and what a 2D slice hides. The results are checked analytically and numerically, including Stokes' theorem and the divergence theorem.

---

## 🔍 What is covered

<table align="center">
  <thead>
    <tr>
      <th align="center">Part</th>
      <th align="center">Topic</th>
      <th align="center">Main idea</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center">1</td>
      <td align="center">🌀 <b>Rotational field</b></td>
      <td align="center">No sources or sinks, but constant rotation. Circulation on the unit circle is checked with Stokes' theorem</td>
    </tr>
    <tr>
      <td align="center">2</td>
      <td align="center">🪑 <b>Saddle field</b></td>
      <td align="center">Zero divergence and zero curl, same magnitude as the rotational field but a completely different behaviour. It derives from a scalar potential</td>
    </tr>
    <tr>
      <td align="center">3</td>
      <td align="center">🌊 <b>Uniform plane wave</b></td>
      <td align="center">A wave in free space travelling in an arbitrary direction. E and H are examined one by one with their curl, divergence and gradient of magnitude</td>
    </tr>
  </tbody>
</table>

For the plane wave, ten figures show E and H, their curls, their divergences and the gradients of their magnitudes, each in 3D and on the three coordinate planes. Two more figures show the wave itself in 3D and the plane transverse to the direction of propagation.

---

## ✅ Key results

<table align="center">
  <thead>
    <tr>
      <th align="center">Field</th>
      <th align="center">Divergence</th>
      <th align="center">Curl</th>
      <th align="center">Conclusion</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center">Rotational field</td>
      <td align="center">Zero</td>
      <td align="center">Constant, along +z</td>
      <td align="center">Solenoidal but not irrotational. Circulation equals the flux of the curl</td>
    </tr>
    <tr>
      <td align="center">Saddle field</td>
      <td align="center">Zero</td>
      <td align="center">Zero</td>
      <td align="center">Both solenoidal and irrotational. Circulation and net flux vanish</td>
    </tr>
    <tr>
      <td align="center">Plane wave (E and H)</td>
      <td align="center">Zero</td>
      <td align="center">In quadrature with the fields</td>
      <td align="center">Agrees with Faraday's and Ampère's laws. E, H and the direction of propagation form a right-handed orthogonal triad</td>
    </tr>
  </tbody>
</table>

A few more takeaways:

- 📏 The magnitude of a field alone says nothing about its divergence or curl.
- 🧭 The gradient of the magnitude of E and H points along the direction of propagation.
- ⚡ The Poynting vector always points along the direction of propagation.
- 🔎 A 2D plot hides the component normal to its plane, so it has to be read with care.

---

## 🧪 Wave parameters used

<table align="center">
  <thead>
    <tr>
      <th align="center">Parameter</th>
      <th align="center">Value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center">Polar angle of the direction of propagation</td>
      <td align="center">50°</td>
    </tr>
    <tr>
      <td align="center">Azimuth of the direction of propagation</td>
      <td align="center">35°</td>
    </tr>
    <tr>
      <td align="center">Polarisation angle</td>
      <td align="center">30°</td>
    </tr>
    <tr>
      <td align="center">Electric field amplitude</td>
      <td align="center">1 V/m</td>
    </tr>
    <tr>
      <td align="center">Plotted region</td>
      <td align="center">±0.75 wavelength on each axis</td>
    </tr>
  </tbody>
</table>

These are set at the top of each script, so the same code can be reused for other directions, polarisations or phases.

---

## 🧰 Tools

<table align="center">
  <thead>
    <tr>
      <th align="center">Purpose</th>
      <th align="center">Tool</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center">🖥️ Computation and plotting</td>
      <td align="center"><b>MATLAB R2024a</b></td>
    </tr>
    <tr>
      <td align="center">🧮 Operators</td>
      <td align="center"><code>divergence</code>, <code>curl</code></td>
    </tr>
    <tr>
      <td align="center">🏹 Vector plots</td>
      <td align="center"><code>quiver</code>, <code>quiver3</code></td>
    </tr>
    <tr>
      <td align="center">🖼️ Figure export</td>
      <td align="center"><code>exportgraphics</code> (400 to 600 dpi)</td>
    </tr>
  </tbody>
</table>

---

## 📚 References

- W. H. Hayt Jr. and J. A. Buck, *Engineering Electromagnetics*, McGraw-Hill.
- D. K. Cheng, *Field and Wave Electromagnetics*, Pearson.
- C. A. Balanis, *Advanced Engineering Electromagnetics*, Wiley.
- D. J. Griffiths, *Introduction to Electrodynamics*, Cambridge University Press.
- The MathWorks, Inc., MATLAB documentation.

---

## 📬 Contact

Reach out for collaborations, discussions, or opportunities. 🤝

<table align="center">
  <tr>
    <td align="center" width="120">
      <a href="https://github.com/Mohammed-Nasr-Aldin">
        <img src="https://cdn.simpleicons.org/github/FFFFFF#gh-dark-mode-only" alt="GitHub" width="44" height="44">
        <br><sub><b>GitHub</b></sub>
      </a>
    </td>
    <td align="center" width="120">
      <a href="https://www.linkedin.com/in/mohammed-nasreldin">
        <img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/linkedin/linkedin-original.svg" alt="LinkedIn" width="44" height="44">
        <br><sub><b>LinkedIn</b></sub>
      </a>
    </td>
    <td align="center" width="120">
      <a href="https://wa.me/201156108363">
        <img src="https://cdn.simpleicons.org/whatsapp/25D366" alt="WhatsApp" width="44" height="44">
        <br><sub><b>WhatsApp</b></sub>
      </a>
    </td>
    <td align="center" width="120">
      <a href="mailto:mohammednasrsmail@gmail.com">
        <img src="https://cdn.simpleicons.org/gmail/EA4335" alt="Email" width="44" height="44">
        <br><sub><b>Email</b></sub>
      </a>
    </td>
  </tr>
</table>

<p align="center"><sub>© Mohammed Nasr Eldin</sub></p>
