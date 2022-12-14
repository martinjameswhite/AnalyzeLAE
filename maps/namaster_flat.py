#  - Lx and Ly: the size of the patch in the x and y dimensions (in radians)
Lx,Ly = 20. * np.pi/180 ,  40. * np.pi/180
#  - Nx and Ny: the number of pixels in the x and y dimensions
Nx,Ny = 600,400

# Define a mask.
# You can also apodize it in the same way you do for full-sky masks:
mask = nmt.mask_apodization_flat(mask, Lx, Ly, aposize=2., apotype="C1")
# Make the "map", or in NaMaster language the "field".
f0 = nmt.NmtFieldFlat(Lx, Ly, mask, [scalarmap])
# Set up your favorite binning scheme:
l0_bins = np.arange(Nx/8)      * 8 * np.pi/Lx
lf_bins = (np.arange(Nx/8)+1) * 8 * np.pi/Lx
b = nmt.NmtBinFlat(l0_bins, lf_bins)
# The effective ells for these bandpowers can be obtained calling:
ells_uncoupled = b.get_effective_ells()

#then you simply call compute_full_master_flat:
#
#https://namaster.readthedocs.io/en/latest/pymaster.html?highlight=compute_full_master#pymaster.workspaces.compute_full_master_flat

#If you want the mode coupling matrix (normally I do) then it's

w00 = nmt.NmtWorkspaceFlat()
w00.compute_coupling_matrix(f0, f0, b)
