<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" class="logo" width="120"/>

# Dynamic LEC with Voronoi Diagrams

## Overview

This project implements an interactive visualization and computation of the **Largest Empty Circle (LEC) with dynamic obstacles** using Voronoi diagrams. The main focus of this research project is the development and demonstration of efficient algorithms for computing the LEC in environments where obstacles can move along user-defined paths. The system leverages computational geometry techniques and real-time animation to provide an intuitive exploration of LEC behavior in dynamic settings.

Other sub-projects in this repository provide foundational algorithms and tools (such as static LEC computation, Voronoi diagram generation, and convex hull algorithms) that underpin and support the advanced features of this dynamic LEC project.

## Features

- **Interactive GUI:**
Built with Tkinter and Matplotlib, allowing users to:
    - Add Voronoi sites by clicking on the canvas
    - Create and animate obstacles with custom paths
    - Toggle visualization overlays (Voronoi diagram, convex hull, LEC)
    - Control animation speed and reset simulation time
- **Dynamic Obstacles:**
Obstacles can follow user-defined paths, and their positions are updated in real-time. The LEC is recalculated dynamically as obstacles move.
- **Efficient LEC Computation:**
Utilizes optimized geometric algorithms (including incremental and sampled Delaunay triangulation, candidate sampling, and pruning) to efficiently compute the largest empty circle that avoids all sites and moving obstacles.
- **Visualization:**
    - Voronoi diagram edges (blue)
    - Convex hull of sites (green)
    - Largest Empty Circle (orange)
    - Obstacles and their paths
- **Statistics \& Controls:**
Real-time display of the number of sites, obstacles, and simulation time, with controls for clearing, pausing, and resuming the simulation.


## Main Focus

> **This project—Dynamic LEC with Voronoi Diagrams—is the central focus of the research on computing the Largest Empty Circle in the presence of dynamic obstacles.**
> All other sub-projects in this repository serve as foundational work, providing essential algorithms and components that enable the advanced capabilities demonstrated here.

## Getting Started

### Requirements

- Python 3.x
- Tkinter (usually included with Python)
- matplotlib
- numpy
- scipy

Install required packages (if not already installed):

```bash
pip install matplotlib numpy scipy
```


### Running the Application

```bash
python dynamic_lec_voronoi.py
```


## Usage

- **Add Sites:**
Select "Add Sites" mode and click on the canvas to place Voronoi sites (red dots).
- **Add Obstacles:**
Switch to "Add Obstacles" mode, click to define the obstacle's path. Each obstacle can have a custom radius and speed.
- **Animate:**
Use the Play/Pause button to start or stop obstacle movement. Adjust speed as needed.
- **Visualization Toggles:**
Enable or disable overlays for Voronoi diagram, convex hull, or LEC.
- **Clear All:**
Remove all sites and obstacles to start a new simulation.


## Repository Structure

- **dynamic_lec_voronoi.py**
Main application for dynamic LEC computation and visualization (this project).
- **[Other sub-projects]**
    - Static LEC computation
    - Static Voronoi diagram generation
    - Convex hull algorithms
    - Utility modules for geometric primitives

These sub-projects are designed to provide the underlying logic and algorithms that make the dynamic LEC project possible.

## Research Context

This project explores the problem of finding the **Largest Empty Circle** in environments with both static and dynamic obstacles, an important challenge in computational geometry with applications in robotics, path planning, and spatial analysis. The dynamic aspect—handling moving obstacles in real time—distinguishes this work and is the primary contribution of this research.

## Acknowledgments

This project builds upon classic computational geometry algorithms and leverages open-source libraries such as scipy, numpy, and matplotlib for efficient computation and visualization.

**For more details on the algorithms and optimizations used, please refer to the in-code documentation and comments.**

