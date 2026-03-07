# Hydrogen Quantum Orbital Visualizer (Web Edition)

## Overview

This project is a **web-based refactoring and reimplementation** of the original **Hydrogen Quantum Orbital Visualizer** by **Kavan Anderson**.

The original project is written in **C++ / OpenGL** and provides multiple renderers (raytracer, realtime, 2D) for visualizing hydrogen atom quantum orbitals derived from the Schrödinger equation.

**Original project (author & reference)**  
- Repository: https://github.com/kavan010/Atoms  
- Original web demo: https://www.kavang.com/atom  

This repository **does not replace** the original work.  
It is a **conceptual and technical adaptation for the web platform**, rewritten from scratch in **HTML, CSS and JavaScript**, targeting modern browsers.

---

## What this version does

- Reimplements **hydrogen orbital probability visualization** for the web
- Uses the same **physical and mathematical foundations**:
  - Quantum numbers *(n, l, m)*
  - Schrödinger equation–based orbital shapes
  - Probability density–driven rendering
- Runs **entirely in the browser**
- Requires **no compilation or native dependencies**
- Focuses on:
  - Real-time interaction
  - Educational clarity
  - Portability and accessibility

---

## Key differences from the original project

| Aspect | Original (C++) | This project (Web) |
|------|----------------|--------------------|
| Language | C++17 | JavaScript (ES6+) |
| Rendering | OpenGL (native) | WebGL / Canvas |
| Build system | CMake + vcpkg | None |
| Platform | Desktop | Browser |
| Codebase | Original implementation | Full rewrite (no shared source code) |

⚠️ This is **not a direct port** and **does not reuse the original C++ code**.  
It is a **ground-up reimplementation** inspired by the same scientific model.

---

## Technologies used

- HTML5  
- CSS3  
- JavaScript (ES6+)  
- WebGL / Canvas  

---

## How to run

Open `index.html` in a modern browser.  
For best WebGL compatibility, serve it via a local HTTP server:

```bash
npx serve .