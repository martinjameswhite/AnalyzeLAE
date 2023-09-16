# Contains "fiducial" distances and interloper fractions for
# the various sources.  This way all scripts can import this
# information without the need to synchronize.
#
# The distances are from make_table while the interloper
# fractions are from VI, assuming redshift failures have
# the same z-distribution as non-failures.

# Comoving distance to field center.  Used to convert angles
# to distances.
chi_dict = {}
chi_dict['N419'] = 3941.0
chi_dict['N501'] = 4448.0
chi_dict['N673'] = 5160.0

# Interloper fractions.
fint_dict = {}
fint_dict['N501'] = {'s0':0.1066,'s1':0.0307,'s2':0.0280,'s3':0.0360}
