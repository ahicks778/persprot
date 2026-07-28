# persprot: Perspective Rotation Kinematics

**`persprot`** is a lightweight Python package designed to calculate and apply perspective rotation corrections for Milky Way satellite galaxies and globular clusters. 

Perspective rotation is an apparent radial velocity gradient across a stellar system that arises purely from the object's systemic proper motion and finite angular extent on the sky. This package provides tools to assess the magnitude of this effect and correct measured line-of-sight velocities ($v_{\mathrm{los}}$) for individual stars, directly implementing the framework presented in **Hicks & Geha (2026)**.

The package seamlessly integrates with the [Local Volume Database (LVDB)](https://github.com/apace7/local_volume_database) to automatically fetch systemic proper motions, distances, and centers. All calculations utilize `astropy` native unit equivalencies to ensure strict numerical precision.

---

## Installation

You can install `persprot` directly via PyPI using `pip`:

```bash
pip install persprot
```

To install the most up-to-date development version directly from GitHub:

```bash
pip install "git+[https://github.com/ahicks778/persprot.git](https://github.com/ahicks778/persprot.git)"
```

> **Note for macOS / zsh users:** Always wrap the git URL in double quotes as shown above to avoid shell syntax errors (`zsh: unknown file attribute`).

---

## API Reference & Core Functions

### `pcorr`

Applies perspective rotation corrections directly to measured line-of-sight velocities ($v_{\mathrm{los}}$). This is the primary function for reducing kinematic datasets.

> **Important:** You must specify **either** a `system_name` (to automatically fetch parameters from the LVDB) **or** provide the manual systemic parameters (`RA0`, `DEC0`, `pmra`, `pmdec`, `D`).

```python
pcorr(RA, DEC, v_los, e_v_los=None, RA0=None, DEC0=None, system_name=None,
      pmra=None, pmdec=None, D=None, e_pmra=None, e_pmdec=None, e_D=None, verbose=True)
```

#### Parameters
* **`RA`**, **`DEC`** (*astropy.units.Quantity*): Celestial coordinates of target stars.
* **`v_los`** (*astropy.units.Quantity*): Measured line-of-sight velocities (e.g., `u.km/u.s`).
* **`e_v_los`** (*astropy.units.Quantity, optional*): Velocity measurement uncertainties (e.g., `u.km/u.s`).
* **`RA0`**, **`DEC0`** (*astropy.units.Quantity, optional*): System center coordinates. Automatically fetched from LVDB if omitted.
* **`system_name`** (*str, optional*): Target system name for automatic LVDB query.
* **`pmra`**, **`pmdec`**, **`D`** (*astropy.units.Quantity, optional*): Systemic proper motions and distance.
* **`e_pmra`**, **`e_pmdec`**, **`e_D`** (*Quantity, float, or list, optional*):  Uncertainties associated with proper motions and distance (supports symmetric or asymmetric `[lower, upper]` inputs).
* **`verbose`** (*bool, optional*): If `True`, prints confirmation output. Defaults to `True`.

#### Returns
* **`v_pc`** (*astropy.units.Quantity*): Corrected line-of-sight velocities ($v_{\mathrm{los}} - \Delta p$) in `km/s`.
* **`e_v_pc`** (*astropy.units.Quantity*): Total velocity uncertainties added in quadrature ($\sqrt{e_{v_{\mathrm{los}}}^2 + \sigma_{\Delta p}^2}$) in `km/s`.

---

### `pcorr_values`

Computes the perspective rotation velocity bias ($\Delta p$) and associated uncertainty ($\sigma_{\Delta p}$) for target stars without modifying velocity arrays.

> **Important:** You must specify **either** a `system_name` (to automatically fetch parameters from the LVDB) **or** provide the manual systemic parameters (`RA0`, `DEC0`, `pmra`, `pmdec`, `D`).

```python
pcorr_values(RA, DEC, RA0=None, DEC0=None, system_name=None, pmra=None, pmdec=None, D=None, 
             e_pmra=None, e_pmdec=None, e_D=None, verbose=True)
```

#### Parameters
* **`RA`**, **`DEC`** (*astropy.units.Quantity*): Right Ascension and Declination of individual stars (e.g., `u.deg`).
* **`RA0`**, **`DEC0`** (*astropy.units.Quantity, optional*): System center coordinates (e.g., `u.deg`). Automatically fetched from LVDB if omitted.
* **`system_name`** (*str, optional*): Target name to query system center, proper motion, and distance from LVDB.
* **`pmra`**, **`pmdec`** (*astropy.units.Quantity, optional*): Proper motion vector components (e.g., `u.mas/u.yr`).
* **`D`** (*astropy.units.Quantity, optional*): Distance (e.g., `u.kpc`).
* **`e_pmra`**, **`e_pmdec`**, **`e_D`** (*Quantity, float, or list, optional*):  Uncertainties associated with proper motions and distance (supports symmetric or asymmetric `[lower, upper]` inputs).
* **`verbose`** (*bool, optional*): If `True`, prints calculation details to the console. Defaults to `True`.

#### Returns
* **`delta_p`** (*astropy.units.Quantity*): Perspective velocity correction array for each star in `km/s`.
* **`sigma_dp`** (*astropy.units.Quantity*): Propagated uncertainty array for the velocity corrections in `km/s`.

---

### `pgrad`

Calculates the physical perspective gradient magnitude and Position Angle (PA) across a stellar system.

> **Important:** You must specify **either** a `system_name` (to automatically fetch parameters from the LVDB) **or** provide the manual systemic parameters (`pmra`, `pmdec`, `rh`, `D`).

```python
pgrad(system_name=None, pmra=None, pmdec=None, rh=None, D=None, 
      e_pmra=None, e_pmdec=None, e_rh=None, e_D=None, scale_unit='rh', verbose=True)
```

#### Parameters
* **`system_name`** (*str, optional*): Name of the system in the LVDB. If provided, missing systemic parameters will be fetched automatically.
* **`pmra`**, **`pmdec`** (*astropy.units.Quantity, optional*): Systemic proper motions (e.g., `u.mas/u.yr`).
* **`rh`** (*astropy.units.Quantity, optional*): Half-light radius as an angle (e.g., `u.arcmin`) or physical length (e.g., `u.pc`).
* **`D`** (*astropy.units.Quantity, optional*): Distance to the system (e.g., `u.kpc`).
* **`e_pmra`**, **`e_pmdec`**, **`e_rh`**, **`e_D`** (*Quantity, float, or list, optional*): Uncertainties associated with proper motions, half-light radius, and distance (supports symmetric or asymmetric `[lower, upper]` inputs).
* **`scale_unit`** (*str or astropy.units.Quantity, optional*): Spatial scale over which to calculate the gradient. Defaults to `'rh'`.
* **`verbose`** (*bool, optional*): If `True`, prints calculated values to the console. Defaults to `True`.

#### Returns
* **`gradient`** (*tuple*): `(val, e_val)` where `val` is the gradient magnitude in `km/s` and `e_val` is its uncertainty (or asymmetric uncertainties `(el_val, eu_val)`).
* **`pa`** (*tuple*): `(pa_val, e_pa)` where `pa_val` is the position angle in degrees and `e_pa` is its uncertainty.