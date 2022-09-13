"""
Paramters used in script_parallel_acquisition.py

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""

# Correspondences between files from the 2nd RSA (with old acquisition script)
# and files from the 1st RSA (using ViBR)
correspondences = {'NW_843_fil_5_spec_0001':'Fil 5_2',
                   'NW_843_fil_5_spec_0002':'Fil 5_3',
                   'NW_843_fil_5_spec_0003':'Fil 5_4',
                   'NW_843_fil_5_spec_0004':'Fil 5_5',
                   'NW_843_fil_5_spec_0005':'Fil 5_6',
                   'NW_843_fil_5_spec_0006':'Fil 5_7',
                   'NW_843_fil_5_spec_0007':'Fil 5_8',
                   'NW_843_fil_5_spec_0008':'Fil 5_10',
                   'NW_843_fil_5_spec_0009':'Fil 5_11',

                   'NW_843_fil_6_spec_0011':'Fil 6_1',
                   'NW_843_fil_6_spec_0012':'Fil 6_2',
                   'NW_843_fil_6_spec_0013':'Fil 6_3',
                   'NW_843_fil_6_spec_0014':'Fil 6_4',

                   'NW_843_fil_6_spec_0016':'Fil 6_5',
                   'NW_843_fil_6_spec_0017':'Fil 6_6',
                   'NW_843_fil_6_spec_0018':'Fil 6_7',
                   'NW_843_fil_6_spec_0019':'Fil 6_8',
                   }

# Direction of measurements
dd_sens = {'NW_843_fil_5_spec_0001':'10',
           'NW_843_fil_5_spec_0002':'01',
           'NW_843_fil_5_spec_0003':'10',
           'NW_843_fil_5_spec_0004':'10',
           'NW_843_fil_5_spec_0005':'01',
           'NW_843_fil_5_spec_0006':'10',
           'NW_843_fil_5_spec_0007':'01',
           'NW_843_fil_5_spec_0008':'10',
           'NW_843_fil_5_spec_0009':'01',

           'NW_843_fil_6_spec_0011':'01',
           'NW_843_fil_6_spec_0012':'10',
           'NW_843_fil_6_spec_0013':'01',
           'NW_843_fil_6_spec_0014':'10',

           'NW_843_fil_6_spec_0016':'10',
           'NW_843_fil_6_spec_0017':'01',
           'NW_843_fil_6_spec_0018':'01',
           'NW_843_fil_6_spec_0019':'10',
           }

# Threshold for NW edges detection
# Default value if not in this dictionary : -60
dd_threshold_NW_edges = {'NW_843_fil_5_spec_0001':-62,
                         'NW_843_fil_5_spec_0002':-62,
                         'NW_843_fil_5_spec_0003':-62,
                         'NW_843_fil_5_spec_0004':-62.5,
                         'NW_843_fil_5_spec_0005':-62,
                         'NW_843_fil_5_spec_0006':-62.5,
                         'NW_843_fil_5_spec_0007':-62,
                         'NW_843_fil_5_spec_0009':-65,

                         'NW_843_fil_6_spec_0013':-60,

                         'NW_843_fil_6_spec_0016':-77,
                         'NW_843_fil_6_spec_0017':-77,
                         'NW_843_fil_6_spec_0018':-75.5,
                         'NW_843_fil_6_spec_0019':-80,
                         }

# Nb of dB above noise level for peak detection
# Default value if not in this dictionary : 12.5
dd_threshold_peak = {'NW_843_fil_5_spec_0001':12.5,
                     'NW_843_fil_5_spec_0009':12.,

                     'NW_843_fil_6_spec_0011':12.6,
                     'NW_843_fil_6_spec_0013':12.5,

                     'NW_843_fil_6_spec_0016':12.5,
                     'NW_843_fil_6_spec_0017':12.4,
                     'NW_843_fil_6_spec_0018':13.5,
                     'NW_843_fil_6_spec_0019':12.9,
                     }

# Minimum frequency for mode 2.1
dd_freq_min1 = {'NW_843_fil_5_spec_0001':None,
                'NW_843_fil_5_spec_0008':5.06e6,
                'NW_843_fil_5_spec_0009':5.09e6,

                'NW_843_fil_6_spec_0011':3.06e6,
                'NW_843_fil_6_spec_0012':3.1e6,
                'NW_843_fil_6_spec_0013':3.15e6,
                'NW_843_fil_6_spec_0014':3.15e6,

                'NW_843_fil_6_spec_0016':560e3,
                'NW_843_fil_6_spec_0017':585e3,
                'NW_843_fil_6_spec_0018':590e3,
                'NW_843_fil_6_spec_0019':605e3,
                }

# Frequency between modes 2.1 and 2.2
dd_freq_middle = {'NW_843_fil_5_spec_0001':None,
                  'NW_843_fil_5_spec_0008':5.08e6,
                  'NW_843_fil_5_spec_0009':5.13e6,

                  'NW_843_fil_6_spec_0011':lambda t: 3.11e6 if t<1500 else 3.15e6,
                  'NW_843_fil_6_spec_0012':3.18e6,
                  'NW_843_fil_6_spec_0013':3.22e6,
                  'NW_843_fil_6_spec_0014':3.24e6,

                  'NW_843_fil_6_spec_0016':570e3,
                  'NW_843_fil_6_spec_0017':595e3,
                  'NW_843_fil_6_spec_0018':610e3,
                  'NW_843_fil_6_spec_0019':615e3,
                  }

# Maximum frequency for mode 2.2
dd_freq_max2 = {'NW_843_fil_5_spec_0001':None,
                'NW_843_fil_5_spec_0008':None,
                'NW_843_fil_5_spec_0009':None,

                'NW_843_fil_6_spec_0011':3.25e6,
                'NW_843_fil_6_spec_0012':3.25e6,
                'NW_843_fil_6_spec_0013':3.30e6,
                'NW_843_fil_6_spec_0014':3.33e6,

                'NW_843_fil_6_spec_0016':585e3,
                'NW_843_fil_6_spec_0017':605e3,
                'NW_843_fil_6_spec_0018':620e3,
                'NW_843_fil_6_spec_0019':630e3,
                }
