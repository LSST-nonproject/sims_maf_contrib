import numpy as np

__all__ = ['findTelescopes']


def findTelescopes(minSize=3.):
    """
    Return an array of telescopes that are above minSize limit
    """
    ## Aperture  Name Location http://astro.nineplanets.org/bigeyes.html
    telescopes =[
        [10.4, 'Gran Canarias','La Palma'],
        [10.0, 'Keck','Mauna Kea'],
        [10.0, 'Keck II', 'Mauna Kea'],
        [9.2, 'SALT', 'South African Astronomical Observatory'],
        [9.2,'Hobby-Eberly', 'Mt. Fowlkes'],
        [8.4,'Large Binocular Telescope','Mt. Graham'],
        [8.3,'Subaru','Mauna Kea'],
        [8.2,'Antu','Cerro Paranal'],
        [8.2,'Kueyen','Cerro Paranal'],
        [8.2,'Melipal','Cerro Paranal'],
        [8.2,'Yepun','Cerro Paranal'],
        [8.1,'Gemini North','Mauna Kea'],
        [8.1,'Gemini South','Cerro Pachon'],
        [6.5,'MMT','Mt. Hopkins'],
        [6.5,'Walter Baade','La Serena'],
        [6.5,'Landon Clay','La Serena'],
        [6.0,'Bolshoi Teleskop Azimutalnyi','Nizhny Arkhyz'],
        [6.0,'LZT','British Columbia'],
        [5.0,'Hale','Palomar Mountain'],
        [4.3,'Dicovery Channel','Lowell Observatory'],
        [4.2,'William Herschel','La Palma'],
        [4.2,'SOAR','Cerro Pachon'],
        [4.2,'LAMOST','Xinglong Station'],
        [4.0,'Victor Blanco','Cerro Tololo'],
        [4.0,'Vista','Cerro Paranal'],
        [3.9,'Anglo-Australian','Coonabarabran'],
        [3.8,'Mayall','Kitt Peak'],
        [3.8,'UKIRT','Mauna Kea'],
        [3.6,'360','Cerro La Silla'],
        [3.6,'Canada-France-Hawaii','Mauna Kea'],
        [3.6,'Telescopio Nazionale Galileo','La Palma'],
        [3.5,'MPI-CAHA','Calar Alto'],
        [3.5,'New Technology','Cerro La Silla'],
        [3.5,'ARC','Apache Point'],
        [3.5,'WIYN','Kitt Peak'],
        [3.0,'Shane','Mount Hamilton'],
        [3.0,'NASA IRTF','Mauna Kea'],
    ]

    scopes = np.zeros(len(telescopes), dtype = zip(
        ['apperture','name','lat','lon'],[float, '|S38', float, float]))


    # name, lat (S negative), lon (W negative)
    observatories = [
        ['Cerro Paranal', -24, 38, -70, 24],
        ['Nizhny Arkhyz', 43, 39, 41, 26],
        ['Cerro La Silla', -29, 15, -70, 44],
        ['Lowell Observatory', 35, 12, -111, 40],
        ['Apache Point', 32, 47, -105, 49],
        ['Mount Hamilton', 37, 21, -121, 38],
        ['South African Astronomical Observatory', -32, 23, 20 ,49],
        ['Cerro Pachon',  -30, 20, -70, 59],
        ['Coonabarabran',-31, 17, 149, 04],
        ['Mt. Fowlkes', 30,40,-104,1],
        ['La Palma', 28, 46, -17, 53],
        ['Mt. Graham',32, 42, -109, 53],
        ['Calar Alto', 37,13,-2,33],
        ['British Columbia',49,17, -122,34],
        ['Kitt Peak',31,57,-111,37],
        ['La Serena',-30,10,-70,48],
        ['Palomar Mountain',33,21,-116,52],
        ['Xinglong Station',40,23,105,50],
        ['Mt. Hopkins', 31,41, -110,53],
        ['Cerro Tololo',-30,10,-70,49],
        ['Mauna Kea', 19,50,-155,28]
    ]

    # Make a nice little dict to look up the observatory positions
    obs = {}
    for i,ob in enumerate(observatories):
        obs[ob[0]] = [(np.abs(ob[1])+ob[2]/60.)*(ob[1]/np.abs(ob[1])),
                      (np.abs(ob[3])+ob[4]/60.)*(ob[3]/np.abs(ob[3]))]

    for i,telescope in enumerate(telescopes):
        scopes['apperture'][i] = telescope[0]
        scopes['name'][i] = telescope[1]
        scopes['lat'][i],scopes['lon'][i] = obs[telescope[2]]

    scopes = scopes[np.where(scopes['apperture'] >= minSize)]

    return scopes
