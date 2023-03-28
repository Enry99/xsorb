from ase.visualize import view
from ase.io.pov import get_bondpairs
from ase.io import read, write
import glob, sys, os, shutil
from slab import Slab
from settings import Settings
from filenames import *

def plot_adsorption_sites():

    settings = Settings()
    slab = Slab(settings.slab_filename, surface_sites_height=settings.surface_height)
  
    slab.find_adsorption_sites(* {
        "distance":0, 
        'symm_reduce':settings.symm_reduce, 
        'near_reduce':settings.near_reduce, 
        'no_obtuse_hollow':True}.values(),
         save_image=True)


def config_images(which : str, povray = False, witdth_res=3000):

    if witdth_res is None and povray: witdth_res = 3000 
    if which == 'scf':
        prefix = pwi_prefix
        pw = 'pwi'
    elif which == 'relax':
        prefix = pwo_prefix
        pw = 'pwo'
    print('Reading files...')
    pw_list=glob.glob(prefix+which+'_*.'+pw)
    configs = [(read(file) if pw == 'pwi' else read(file, results_required=False)) for file in pw_list]
    print('All files read.')

    print('Saving images...')
    if(not os.path.exists(which+'_'+images_dirname)):
        os.mkdir(which+'_'+images_dirname)
    os.chdir(which+'_'+images_dirname) 

    for i, config in enumerate(configs):
        label = pw_list[i].split('.'+pw)[0].split('_')[-1]

        if(povray):
            config_copy = config.copy()
            config_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas
            write(prefix+which+'_{0}_pov.pov'.format(label), 
                config, 
                format='pov',
                radii = 0.65, 
                rotation='-10z,-80x', 
                povray_settings=dict(canvas_width=witdth_res, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(config_copy, radius=0.95))
                #camera_type='perspective'
            ).render()
            os.remove(prefix+which+'_{0}_pov.pov'.format(label))
            os.remove(prefix+which+'_{0}_pov.ini'.format(label))

        else:
            write(prefix+which+'_{0}.png'.format(label), config, rotation='-10z,-80x', scale = 100)
            write(prefix+which+'_{0}_top.png'.format(label), config, scale = 100)
    print('All images saved to {0}.'.format(which+'_'+images_dirname))


def view_config(which : str, index : int):
    if which == 's':
        which = 'scf'
        prefix = pwi_prefix
        pw = 'pwi'
    elif which == 'r':
        which = 'relax'
        prefix = pwo_prefix
        pw = 'pwo'
    file = prefix+which+'_{0}.'.format(index)+pw
    try:
        config = read(file, index=':')
        view(config)
    except:
        config =  read(file) if pw == 'pwi' else read(file, results_required=False)
        print('Not possible to load the full relaxation history because of the ase bug for Espresso > 6.7. Showing only the last configuration.')
        view(config)        


def relax_animations(povray = False, witdth_res=3000):


    print('Reading files...')
    pwo_list=glob.glob(pwo_prefix+'relax_*.pwo')
    configs = [read(file, index=':') for file in pwo_list]
    labels = [pwo.split('.pwo')[0].split('_')[-1] for pwo in pwo_list]
    print('All files read.')
    
    print('Generating animation(s)...')
    if(not os.path.exists('relax_'+images_dirname)):
        os.mkdir('relax_'+images_dirname)
    os.chdir('relax_'+images_dirname)

    if(povray):
        if witdth_res is None: witdth_res = 3000 
        for i, config in enumerate(configs):
            if os.path.exists('temp'):
                shutil.rmtree('temp')
            os.mkdir('temp')  

            os.chdir('temp')
            for j, step in enumerate(config):
                
                step_copy = step.copy()
                step_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas
                write(pwo_prefix+'step_{:04d}.pov'.format(j), 
                    step, 
                    format='pov',
                    radii = 0.65, 
                    rotation='-10z,-80x', 
                    povray_settings=dict(canvas_width=witdth_res, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(step_copy, radius=0.95))
                    #camera_type='perspective'
                ).render()

            os.chdir('..')
            #os.system('convert -delay 20 temp/step_*.png '+pwo_prefix+'relax_{0}_pov.gif'.format(labels[i]))
            os.system('ffmpeg -framerate 10 -i temp/step_%04d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p '+pwo_prefix+'relax_{0}_pov.mp4'.format(labels[i]))
            shutil.rmtree('temp')

    else:
        for i, config in enumerate(configs):
            write(pwo_prefix+'relax_{0}.gif'.format(labels[i]), config, rotation='-10z,-80x', interval=150, scale = 100, save_count=None)

    print('All animations saved to {0}.'.format('relax_'+images_dirname))