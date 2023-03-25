from ase.visualize import view
from ase.io.pov import get_bondpairs
from ase.io import read, write
import glob, sys, os
from slab import Slab
from settings import Settings
from filenames import *

def plot_adsorption_sites():

    settings = Settings()
    slab = Slab(settings.slab_filename)
  
    slab.find_adsorption_sites(* {
        "distance":0, 
        'symm_reduce':settings.symm_reduce, 
        'near_reduce':settings.near_reduce, 
        'no_obtuse_hollow':True}.values(),
         save_image=True)


def config_images(which : str, povray = False):

    if which == 'scf':
        prefix = pwi_prefix
        pw = 'pwi'
    elif which == 'finalrelax':
        prefix = pwo_prefix
        pw = 'pwo'
    print('Reading files...')
    pw_list=glob.glob(prefix+'_'+which+'_*.'+pw)
    configs = [(read(file) if pw == 'pwi' else read(file, results_required=False)) for file in pw_list]
    print('All files read.')

    print('Saving images...')
    if(not os.path.exists(which+'_'+images_dirname)):
        os.mkdir(which+'_'+images_dirname)
    for i, config in enumerate(configs):
        label = pw_list[i].split('.'+pw)[0].split('_')[-1]

        if(povray):
            os.chdir(which+'_'+images_dirname)    

            write(prefix+'_'+which+'_{0}_pov.pov'.format(label), 
                config, 
                format='pov',
                radii = 0.65, 
                rotation='-10z,-80x', 
                povray_settings=dict(canvas_width=4000, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(config, radius=1))
                #camera_type='perspective'
            ).render()
            os.remove(prefix+'_'+which+'_{0}_pov.pov'.format(label))
            os.remove(prefix+'_'+which+'_{0}_pov.ini'.format(label))

            os.chdir('..')
        else:
            write(which+'_'+images_dirname+'/'+prefix+'_'+which+'_{0}.png'.format(label), config, rotation='-10z,-80x', scale = 100)
            write(which+'_'+images_dirname+'/'+prefix+'_'+which+'_{0}_top.png'.format(label), config, scale = 100)
    print('All images saved to {0}.'.format(which+'_'+images_dirname))


def view_config(which : str, index : int):
    if which == 'scf':
        prefix = pwi_prefix
        pw = 'pwi'
    elif which == 'finalrelax':
        prefix = pwo_prefix
        pw = 'pwo'
    file = prefix+'_'+which+'_{0}.'.format(index)+pw
    try:
        config = read(file, index=':')
        view(config)
    except:
        config =  read(file) if pw == 'pwi' else read(file, results_required=False)
        print('Not possible to load the full relaxation history because of the ase bug for Espresso > 6.7. Showing only the last configuration.')
        view(config)        


def relax_animations(index = -1):

    print('Reading files...')
    pwo_list=glob.glob(pwo_prefix+'_finalrelax_*.pwo')
    configs = [read(file, index=':') for file in pwo_list]
    labels = [pwo.split('.pwo')[0].split('_')[-1] for pwo in pwo_list]
    print('All files read.')
    
    print('Generating animation(s)...')
    if index == -1:
        for i, config in enumerate(configs):
            write('finalrelax_'+images_dirname+'/'+pwo_prefix+'_finalrelax_{0}.gif'.format(labels[i]), config, rotation='-10z,-80x', interval=200, scale = 100, save_count=None)
    else:
        if index not in labels:
            print('Selected configuration not found. Quitting.')
            sys.exit(1)
        else:
            config = configs[labels.index(index)]
            write('finalrelax_'+images_dirname+'/'+pwo_prefix+'_finalrelax_{0}.gif'.format(index), config, rotation='-10z,-80x', interval=200, scale = 100, save_count=None)
    print('All animations saved to {0}.'.format('finalrelax_'+images_dirname))