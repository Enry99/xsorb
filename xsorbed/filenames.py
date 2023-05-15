VERSION                  = '1.0'  #program version

rydbergtoev = 13.605703976

hybrid_screening_thresholds = [5e-3, 5e-2]
N_relax_default             = 5

#just a collection of the preset (non-user-defined) filenames
pw_files_prefix             = ''
labels_filename             = 'site_labels.csv'
screening_energies_filename = 'screening_energies.csv'
relax_energies_filename     = 'relax_energies.csv'
images_dirname              = 'images'   
screening_outdir            = 'screening_outdirs'
relax_outdir                = 'relax_outdirs'
jobscript_filename          = 'jobscript' #standard name used in the copied version inside the outdirs

#user-defined atomic colors, with format [atomic number, [r,g,b]]:
#USER_COLORS_DEFS = [
#    [1,  [1.00, 1.00, 1.00]],
#    [6,  [0.35, 0.35, 0.40]],
#    [7,  [0.50, 0.75, 1.00]],
#    [8,  [0.70, 0.00, 0.00]],
#    [14, [1.00, 1.00, 0.22]],
#    [16, [0.90, 0.95, 0.30]],
#    [26, [0.40, 0.40, 1.00]],
#    [42, [0.70, 0.70, 0.70]],
#]

RADIUS_DEFAULT = 0.8  #radius for bondpairs in povray

