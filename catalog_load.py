import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.cosmology import Planck18 as cosmo

path_halo = '/sps/lsst/groups/clusters/cluster_comparison_project/after_matching/cosmoDC2_pywazp.DC2.tpz.T500k.pointEstimateMags/zband/member_matching/fshare_0.1_pref_more_massive/cat1.fits'

path_amas = '/sps/lsst/groups/clusters/cluster_comparison_project/after_matching/cosmoDC2_pywazp.DC2.tpz.T500k.pointEstimateMags/zband/member_matching/fshare_0.1_pref_more_massive/cat2.fits'

def table(path):
    hdulist = fits.open(path)
    data = hdulist[1].data
    hdulist.close()
    return Table(data)


#halo
# on a : id,         ra,       dec,             z,          M
#     'halo_id', 'ra_true', 'dec_true', 'redshift_true', 'm200c'

#cluster
# on a : id,         ra,       dec,             z,          lmbda
#       'id',       'ra',     'dec',           'zp',        'n200'

t_halo = table(path_halo)
t_amas = table(path_amas)
#print(t_galaxies.colnames)
#print(t_amas[:3])

#print(t_halo[0]['ra_true'])
#print(t_amas[0])

#transforme les coordonnées sphériques en carthésiennes
def spherical_to_cartesian(d, ra_deg, dec_deg):
    ra_rad = np.deg2rad(ra_deg)
    dec_rad = np.deg2rad(dec_deg)
    x = d * np.cos(dec_rad) * np.cos(ra_rad)
    y = d * np.cos(dec_rad) * np.sin(ra_rad)
    z = d * np.sin(dec_rad)
    return np.array([x, y, z])

#donne la distance entre deux points dans le ciel repérés par ra, dec et un z
def distance(ra1, dec1, z1, ra2, dec2, z2):
    # Conversion RA/DEC en coordonnées célestes
    coord1 = SkyCoord(ra1 * u.deg, dec1 * u.deg)
    coord2 = SkyCoord(ra2 * u.deg, dec2 * u.deg)
    # Distances comobiles (en Mpc)
    d1 = cosmo.comoving_distance(z1).value
    d2 = cosmo.comoving_distance(z2).value
    #vecteurs position
    vec1 = spherical_to_cartesian(d1, ra1, dec1)
    vec2 = spherical_to_cartesian(d2, ra2, dec2)
    #la norme de vec1-vec2 est la distance estimée entre les deux points
    return np.linalg.norm(vec1 - vec2)

# donne la distance entre un halo et un amas
def L(halo, amas):
    ra1, dec1, z1 = halo["ra_true"], halo["dec_true"], halo["redshift_true"]
    ra2, dec2, z2 = amas["ra"], amas["dec"], amas["zp"]
    return distance(ra1, dec1, z1, ra2, dec2, z2)

#listes contenant au rang i les propriétés du i-ème amas
pop, M, z = [], [], []
#nopmbre d'amas et de halos
n, m = len(t_amas["n200"]), len(t_halo["m200c"])
#ittération sur les amas et sur les halos
for i in range(n):
    mini = (0, L(halo = t_halo[0], amas = t_amas[0]))
    for j in range(1, m):# recherche du halo le plus proche de l'amas i
        D = L(halo = t_halo[j], amas = t_amas[i])
        if mini[1] > D:
            mini = (j, D)
    # récupération des populations, masses et z de l'amas i
    pop.append(t_amas[i]["n200"])
    M.append(t_halo[mini[0]]["m200c"])
    z.append((t_halo[mini[0]]["redshift_true"]+t_amas[i]["zp"])/2)
# création d'une table contenant à la fois la population, la masse et z
t_amas_final = Table()
t_amas_final["l"] = pop 
t_amas_final["M"] = M
t_amas_final["z"] = z

print(t_amas_final)
















