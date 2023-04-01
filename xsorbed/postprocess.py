from ase.visualize import view
from ase.io.pov import get_bondpairs
from ase.io import read, write
import glob, sys, os, shutil
from natsort import natsorted
from slab import Slab
from settings import Settings
from filenames import *

def read_energy(filename: str, Eslab = 0, Emol = 0):
    with open(filename, 'r') as f:
        pwo = f.readlines()

    toten = 0
    scf_terminated = False
    for line in pwo: #make sure to get the last one (useful in relaxations)
        if '!' in line: 
            toten = line.split()[4]
        if 'End of self-consistent calculation' in line:
            scf_terminated = True
    if(scf_terminated and toten != 0): 
        return (float(toten) - Eslab - Emol)*rydbergtoev
    else: return None


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
    from natsort import natsorted
    pw_list=natsorted(glob.glob(prefix+which+'_*.'+pw))
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
                povray_settings=dict(canvas_width=witdth_res, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(config_copy, radius=0.7))
                #camera_type='perspective'
            ).render()
            os.remove(prefix+which+'_{0}_pov.pov'.format(label))
            os.remove(prefix+which+'_{0}_pov.ini'.format(label))

        else:
            write(prefix+which+'_{0}.png'.format(label), config, rotation='-10z,-80x', scale = 100)
            write(prefix+which+'_{0}_top.png'.format(label), config, scale = 100)


    if(which=='relax'):
        os.chdir('..')
        import numpy as np
        from matplotlib import pyplot as plt
        import matplotlib.image as mpimg

        rows_fig = max(int(np.ceil(len(configs)/5)), 1)
        cols_fig = max(int(np.ceil(len(configs)/rows_fig)), 1)
        fig = plt.figure(figsize=(1.5 * cols_fig, 5 * rows_fig / 3))
        fig.subplots_adjust(wspace=0.001)
        axes = [fig.add_subplot(rows_fig,cols_fig,i) for i in range(1,len(configs) + 1)]

        print(rows_fig, cols_fig)
        try:
            E_slab_mol = Settings().E_slab_mol
        except: #option in case settings.in or some other input file is no more present.
            E_slab_mol = []
    
        energies = [read_energy(file, *E_slab_mol) for file in pw_list]

        for i, conf in enumerate(configs):
            label = int(pw_list[i].split('.'+pw)[0].split('_')[-1])
            img = mpimg.imread('relax_'+images_dirname+'/'+prefix+which+'_{0}{1}.png'.format(label, '_pov' if povray else ''))
            axes[i].imshow(img)
            axes[i].set_title('{0}: {1:.3f} eV'.format(label, energies[i]), fontsize = 7, pad=1)
            axes[i].set_xticks([])
            axes[i].set_yticks([])
        fig.savefig('relax_'+images_dirname+'/'+'relax_overview.png', dpi=1500, bbox_inches='tight')

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
                    povray_settings=dict(canvas_width=witdth_res, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(step_copy, radius=0.7))
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


def plot_energy_evolution():

    settings = Settings()
    Eslab, Emol = (settings.E_slab_mol[0], settings.E_slab_mol[1])

    print('Reading files...')
    pwo_list=natsorted(glob.glob(pwo_prefix+'relax_*.pwo'))
    labels = [int(pwo.split('.pwo')[0].split('_')[-1]) for pwo in pwo_list]

    totens = []
    relax_terminated = []
    for file in pwo_list:

        end = False
        totens.append([])

        with open(file, 'r') as f:
            pwo = f.readlines()

        for line in pwo: #make sure to get the last one (useful in relaxations)
            if '!' in line: 
                totens[-1].append( (float(line.split()[4]) - (Eslab+Emol)) * rydbergtoev )
            if 'Final energy' in line:
                end = True
        relax_terminated.append(end)
    print('All files read.')


    from matplotlib import pyplot as plt
    plt.axhline(y=0, linestyle='--', color='black', linewidth=1)
    for i, config_e in enumerate(totens):
        if(False):#len(config_e) > 10): #skip first 10 steps
            plt.plot([*range(10, len(config_e))], config_e[10:], '-', label=labels[i])
            plt.xlim(xmin=10)
        else:
            plt.plot(config_e, '-', label=labels[i])
            if (not relax_terminated[i]): plt.plot(len(config_e)-1, config_e[-1], 'x', color='black')
            #plt.xlim(xmin=0)
            
    
    plt.title('Relax energies')
    plt.xlabel('step')
    plt.ylabel('energy (eV)')
    plt.grid(linestyle='dotted')
    plt.legend(title="Config")
    plt.savefig(energy_plot_filename, dpi=300, bbox_inches='tight')

    print('plot saved in {0}'.format(energy_plot_filename))