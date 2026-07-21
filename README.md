# persprot: Perspective Rotation Kinematics

**`persprot`** is a lightweight Python package designed to calculate and apply perspective rotation corrections for Milky Way satellite galaxies and globular clusters. 

Perspective rotation is an apparent radial velocity gradient across a stellar system that arises purely from the object's systemic proper motion and finite angular extent on the sky. This package provides tools to assess the magnitude of this effect and correct measured line-of-sight velocities ($v_{los}$) for individual stars, directly implementing the framework presented in **Hicks & Geha (2026)**.

The package seamlessly integrates with the [Local Volume Database (LVDB)](https://github.com/apace7/local_volume_database) to automatically fetch systemic proper motions, distances, and centers. All calculations utilize `astropy` native unit equivalencies to ensure strict numerical precision.

---

## Installation

You can install `persprot` directly from GitHub using `pip`:

```bash
pip install git+[https://github.com/YourUsername/persprot.git](https://github.com/YourUsername/persprot.git)