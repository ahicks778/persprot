import numpy as np
import astropy.units as u
import astropy.table as table
from astropy.utils.data import download_file 

_LVDB_CACHE = None

def _get_lvdb():
    """Downloads and caches the unified Local Volume Database table dynamically."""
    global _LVDB_CACHE
    if _LVDB_CACHE is None:
        print("Connecting to Local Volume Database (LVDB) public repository...")
        try:
            version_url = 'https://raw.githubusercontent.com/apace7/local_volume_database/main/code/release_version.txt'
            version_num = table.Table.read(version_url, format='ascii.fast_no_header')['col1'][0]
            csv_url = f'https://github.com/apace7/local_volume_database/releases/download/{version_num}/comb_all.csv'
            local_path = download_file(csv_url, cache=True)
            
            _LVDB_CACHE = table.Table.read(local_path, format='csv')
            _LVDB_CACHE['key_lower'] = [str(x).strip().lower() for x in _LVDB_CACHE['key']]
        except Exception as e:
            raise RuntimeError(f"Could not load the Local Volume Database: {e}")
    return _LVDB_CACHE

def _parse_err(err, default_unit):
    if err is None:
        return (0.0 * default_unit, 0.0 * default_unit, False)
    if isinstance(err, list):
        err = u.Quantity(err)
    try:
        if len(err) == 2:
            return (np.abs(err[0]).to(default_unit), np.abs(err[1]).to(default_unit), True)
    except TypeError:
        pass 
    return (np.abs(err).to(default_unit), np.abs(err).to(default_unit), False)

def _safe_db_err(em, ep, unit):
    if np.isnan(em) or np.isnan(ep):
        return None
    return [em * unit, ep * unit]

def _ensure_quantity(val, expected_physical_types, name):
    if not isinstance(val, u.Quantity):
        raise TypeError(f"'{name}' must be an astropy Quantity (e.g., {val} * u.mas/u.yr).")
    if val.unit.physical_type not in expected_physical_types:
        raise u.UnitsError(f"'{name}' has physical type '{val.unit.physical_type}', expected {expected_physical_types}.")
    return val


def pgrad(system_name=None, pmra=None, pmdec=None, rh=None, D=None, 
          e_pmra=None, e_pmdec=None, e_rh=None, e_D=None, scale_unit='rh', verbose=True):
    """Calculates physical perspective gradient magnitude and Position Angle."""
    
    if system_name is None:
        if any(param is None for param in (pmra, pmdec, rh, D)):
            raise ValueError(
                "You must provide a 'system_name' to fetch missing parameters from the database, "
                "or provide all parameters (pmra, pmdec, rh, D) manually."
            )
        system_name = "Custom System"

    sources = {"pmra": "User Input", "pmdec": "User Input", "rh": "User Input", "D": "User Input"}
    
    ## Database Lookup
    if any(param is None for param in (pmra, pmdec, rh, D)):
        lvdb = _get_lvdb()
        match_idx = np.where(lvdb['key_lower'] == system_name.strip().lower())[0]
        if len(match_idx) == 0:
            raise ValueError(f"System '{system_name}' not found in the Local Volume Database.")
        row = lvdb[match_idx[0]]
        
        if pmra is None:
            pmra = row['pmra'] * u.mas / u.yr
            e_pmra = _safe_db_err(row['pmra_em'], row['pmra_ep'], u.mas/u.yr)
            sources["pmra"] = "LVDB"
            
        if pmdec is None:
            pmdec = row['pmdec'] * u.mas / u.yr
            e_pmdec = _safe_db_err(row['pmdec_em'], row['pmdec_ep'], u.mas/u.yr)
            sources["pmdec"] = "LVDB"
            
        if rh is None:
            if not np.isnan(row['rhalf']):
                rh = row['rhalf'] * u.arcmin
                e_rh = _safe_db_err(row['rhalf_em'], row['rhalf_ep'], u.arcmin)
                sources["rh"] = "LVDB (Angular)"
            elif 'rhalf_physical' in row.colnames and not np.isnan(row['rhalf_physical']):
                rh = row['rhalf_physical'] * u.pc
                e_rh = _safe_db_err(row['rhalf_physical_em'], row['rhalf_physical_ep'], u.pc)
                sources["rh"] = "LVDB (Physical)"
                
        if D is None:
            D = row['distance'] * u.kpc
            e_D = _safe_db_err(row['distance_em'], row['distance_ep'], u.kpc)
            sources["D"] = "LVDB"

    pmra = _ensure_quantity(pmra, ['angular speed'], 'pmra').to(u.mas/u.yr)
    pmdec = _ensure_quantity(pmdec, ['angular speed'], 'pmdec').to(u.mas/u.yr)
    rh = _ensure_quantity(rh, ['length', 'angle'], 'rh')
    D = _ensure_quantity(D, ['length'], 'D').to(u.kpc)
    b_line = rh if scale_unit == 'rh' else _ensure_quantity(scale_unit, ['length', 'angle'], 'scale_unit')

    ## Perspective Gradient Calculation (Equation 11)
    b_is_length = (b_line.unit.physical_type == 'length')
    pmt = np.hypot(pmra, pmdec)
    raw_result = (b_line * pmt) if b_is_length else (b_line * D * pmt)
    val = raw_result.to(u.km / u.s, equivalencies=u.dimensionless_angles())
    pa = (0.0 * u.deg) if pmt.value == 0 else (np.arctan2(pmra, pmdec).to(u.deg)) % (360.0 * u.deg)
    
    el_pmra, eu_pmra, asym_pmra = _parse_err(e_pmra, pmra.unit)
    el_pmdec, eu_pmdec, asym_pmdec = _parse_err(e_pmdec, pmdec.unit)
    el_rh, eu_rh, asym_rh = _parse_err(e_rh, rh.unit)
    el_D, eu_D, asym_D = _parse_err(e_D, D.unit)
    
    el_b = el_rh if scale_unit == 'rh' else 0.0 * b_line.unit
    eu_b = eu_rh if scale_unit == 'rh' else 0.0 * b_line.unit
    
    if pmt.value != 0:
        el_pmt = np.sqrt((pmra * el_pmra)**2 + (pmdec * el_pmdec)**2) / pmt
        eu_pmt = np.sqrt((pmra * eu_pmra)**2 + (pmdec * eu_pmdec)**2) / pmt
        el_pa = (np.sqrt((pmdec * el_pmra)**2 + (pmra * el_pmdec)**2) / (pmt**2) * u.rad).to(u.deg)
        eu_pa = (np.sqrt((pmdec * eu_pmra)**2 + (pmra * eu_pmdec)**2) / (pmt**2) * u.rad).to(u.deg)
    else:
        el_pmt, eu_pmt, el_pa, eu_pa = 0*pmra.unit, 0*pmra.unit, 0*u.deg, 0*u.deg

    if b_is_length:
        el_val = val * np.sqrt((el_b/b_line)**2 + (el_pmt/pmt)**2)
        eu_val = val * np.sqrt((eu_b/b_line)**2 + (eu_pmt/pmt)**2)
    else:
        el_val = val * np.sqrt((el_b/b_line)**2 + (el_D/D)**2 + (el_pmt/pmt)**2)
        eu_val = val * np.sqrt((eu_b/b_line)**2 + (eu_D/D)**2 + (eu_pmt/pmt)**2)

    has_asym = any([asym_pmra, asym_pmdec, asym_rh if scale_unit=='rh' else False, asym_D if not b_is_length else False])

    if verbose:
        print(f"\n=== Gradient Results for {system_name} ===")
        print(f"Gradient: {val.value:.3f} km/s. Position Angle: {pa.value:.1f} deg")
        
    return ((val, (el_val, eu_val)), (pa, (el_pa, eu_pa))) if has_asym else ((val, eu_val), (pa, eu_pa))


def pcorr_values(RA, DEC, RA0=None, DEC0=None, system_name=None, pmra=None, pmdec=None, D=None, 
                 e_pmra=None, e_pmdec=None, e_D=None, verbose=True):
    """Calculates perspective rotation bias (delta_p) and error for specific stars."""
    
    ##Database Lookup 
    if system_name is not None and any(p is None for p in (RA0, DEC0, pmra, pmdec, D)):
        lvdb = _get_lvdb()
        idx = np.where(lvdb['key_lower'] == system_name.strip().lower())[0]
        if len(idx) == 0: raise ValueError(f"System '{system_name}' not found.")
        row = lvdb[idx[0]]
        
        RA0 = (row['ra'] * u.deg) if RA0 is None else RA0
        DEC0 = (row['dec'] * u.deg) if DEC0 is None else DEC0
        pmra = (row['pmra'] * u.mas/u.yr) if pmra is None else pmra
        pmdec = (row['pmdec'] * u.mas/u.yr) if pmdec is None else pmdec
        D = (row['distance'] * u.kpc) if D is None else D
        e_pmra = _safe_db_err(row['pmra_em'], row['pmra_ep'], u.mas/u.yr) if e_pmra is None else e_pmra
        e_pmdec = _safe_db_err(row['pmdec_em'], row['pmdec_ep'], u.mas/u.yr) if e_pmdec is None else e_pmdec
        e_D = _safe_db_err(row['distance_em'], row['distance_ep'], u.kpc) if e_D is None else e_D

    RA = _ensure_quantity(RA, ['angle'], 'RA').to(u.deg)
    DEC = _ensure_quantity(DEC, ['angle'], 'DEC').to(u.deg)
    RA0 = _ensure_quantity(RA0, ['angle'], 'RA0').to(u.deg)
    DEC0 = _ensure_quantity(DEC0, ['angle'], 'DEC0').to(u.deg)
    pmra = _ensure_quantity(pmra, ['angular speed'], 'pmra').to(u.mas/u.yr)
    pmdec = _ensure_quantity(pmdec, ['angular speed'], 'pmdec').to(u.mas/u.yr)
    D = _ensure_quantity(D, ['length'], 'D').to(u.kpc)

    ## BIAS CALCULATION (Equations 8 & 9) 
    delta_alpha = RA - RA0
    delta_delta = DEC - DEC0
    delta_alpha_star = delta_alpha * np.cos(DEC0.to(u.rad)) 
    
    # Kinematic term has units of [angular velocity * angle]
    kinematic_term = (pmra * delta_alpha_star) + (pmdec * delta_delta)
    raw_delta_p = D * kinematic_term
    
    # Convert exactly using astropy's dimensionless angles equivalency
    delta_p = raw_delta_p.to(u.km / u.s, equivalencies=u.dimensionless_angles())

    # Error Calculation (Equation 9)
    if all(e is not None for e in (e_pmra, e_pmdec, e_D)):
        el_p, eu_p, _ = _parse_err(e_pmra, pmra.unit)
        el_d, eu_d, _ = _parse_err(e_pmdec, pmdec.unit)
        el_D, eu_D, _ = _parse_err(e_D, D.unit)
        
        sigma_pmra = 0.5 * (el_p + eu_p)
        sigma_pmdec = 0.5 * (el_d + eu_d)
        sigma_D = 0.5 * (el_D + eu_D)

        term1 = (sigma_D**2) * (kinematic_term**2)
        term2 = (D**2) * (delta_alpha_star**2) * (sigma_pmra**2)
        term3 = (D**2) * (delta_delta**2) * (sigma_pmdec**2)

        var_dp = term1 + term2 + term3
        sigma_dp = np.sqrt(var_dp).to(u.km / u.s, equivalencies=u.dimensionless_angles())
    else:
        sigma_dp = np.zeros_like(delta_p.value) * (u.km / u.s)

    if verbose:
        print(f"Computed {len(np.atleast_1d(RA))} targets. Max correction: {np.max(np.abs(delta_p)):.3f}")

    return delta_p, sigma_dp


def pcorr(RA, DEC, v_los, e_v_los=None, RA0=None, DEC0=None, system_name=None,
          pmra=None, pmdec=None, D=None, e_pmra=None, e_pmdec=None, e_D=None, verbose=True):
    """Corrects radial velocities for perspective rotation."""
    
    v_los = _ensure_quantity(v_los, ['speed'], 'v_los').to(u.km / u.s)
    
    # Calculate the raw bias array using the dedicated function
    delta_p, sigma_dp = pcorr_values(RA, DEC, RA0, DEC0, system_name, pmra, pmdec, 
                                     D, e_pmra, e_pmdec, e_D, verbose=False)
    
    # Equation 7: Correcting the LOS velocity
    v_pc = v_los - delta_p

    if e_v_los is not None:
        e_v_los = _ensure_quantity(e_v_los, ['speed'], 'e_v_los').to(u.km / u.s)
        # Quadrature addition of uncertainties
        e_v_pc = np.sqrt(e_v_los**2 + sigma_dp**2)
    else:
        e_v_pc = sigma_dp

    if verbose:
        print(f"\n=== Perspective Correction Applied for {system_name or 'Custom'} ===")
        print(f"Corrected {len(np.atleast_1d(v_pc))} velocities.")
        
    return v_pc, e_v_pc