from ase.visualize import view
from ase.io.pov import get_bondpairs
from ase.io import read, write
import glob, sys, os, shutil
from natsort import natsorted
from xsorbed.slab import Slab
from xsorbed.settings import Settings
from xsorbed.filenames import *


RADIUS = 0.8  #radius for bondpairs in povray

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


def plot_adsorption_sites(ALL = False):

    settings = Settings(read_energies=False)
    slab = Slab(settings.slab_filename, surface_sites_height=settings.surface_height)
  
    slab.find_adsorption_sites(* {
        "distance":0, 
        'symm_reduce': 0 if ALL else settings.symm_reduce, 
        'near_reduce':settings.near_reduce, 
        'no_obtuse_hollow':True}.values(),
        selected_sites= [] if ALL else settings.selected_sites,
        save_image=True)


def config_images(which : str, povray = False, witdth_res=3000, index : str = None, rotations : str = None):

    if witdth_res is None and povray: witdth_res = 3000 
    if which == 's':
        which = 'scf'
        prefix = pwi_prefix
        pw = 'pwi'
    elif which == 'r':
        which = 'relax'
        prefix = pwo_prefix
        pw = 'pwo'
    print('Reading files...')
    pw_list=natsorted(glob.glob(prefix+which+'_*.'+pw))
    if(not pw_list):
        print("Files not found. Quitting.")
        sys.exit(1)

    configs = []
    labels = []
    uncompleted = []
    for file in pw_list:
        with open(file, 'r') as f:
            lines = f.readlines()
        STOP = False
        for line in lines:
            if 'convergence NOT achieved' in line:
                print("Found 'convergence NOT achieved' in {0}, so relaxation was not completed. It will be skipped.".format(file))
                STOP = True
                uncompleted.append(file)
                break
        if(STOP): continue
        configs.append(read(file) if pw == 'pwi' else read(file, results_required=False))
        labels.append(file.split('.'+pw)[0].split('_')[-1])

    print('All files read.')

    if index is not None:
        configs = [configs[labels.index(index)]] #just one
        labels = [labels[labels.index(index)]]

    print('Saving images...')
    if(not os.path.exists(which+'_'+images_dirname)):
        os.mkdir(which+'_'+images_dirname)
    os.chdir(which+'_'+images_dirname) 


    from ase.data.colors import jmol_colors
    ATOM_COLORS = jmol_colors.copy()
    for color in USER_COLORS_DEFS:
        ATOM_COLORS[color[0]] = color[1]
    colors = [ATOM_COLORS[atom.number] for atom in configs[0]]

    try:
        settings = Settings()
        E_slab_mol = settings.E_slab_mol
        slab_filename = settings.slab_filename
        Nbulk = len(read(filename=slab_filename, results_required=False) if slab_filename.split('.')[-1]=='pwo' else read(filename=slab_filename))
    except Exception as e:
        print("Error while reading settings.in: ", e, "Energies will be given as total energies, and top view will not have faded surface atoms.")
        Nbulk = 0
        E_slab_mol = [0, 0]


    for i, config in enumerate(configs):
        label = labels[i]

        if(povray):
            config_copy = config.copy()
            config_copy.set_pbc([0,0,0]) #to avoid drawing bonds with invisible replicas

            if rotations is not None: #use specified rotation
                write(prefix+which+'_{0}_{1}_pov.pov'.format(label, rotations.replace(',','_')), 
                    config, 
                    format='pov',
                    radii = 0.65, 
                    rotation=rotations,
                    colors=colors,
                    povray_settings=dict(canvas_width=witdth_res, celllinewidth=0, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(config_copy, radius=RADIUS))
                    #camera_type='perspective'
                ).render()
                os.remove(prefix+which+'_{0}_{1}_pov.pov'.format(label, rotations.replace(',','_')))
                os.remove(prefix+which+'_{0}_{1}_pov.ini'.format(label, rotations.replace(',','_')))
            else: #front and top view
                #front view
                write(prefix+which+'_{0}_pov.pov'.format(label), 
                    config, 
                    format='pov',
                    radii = 0.65, 
                    rotation='-10z,-80x', 
                    colors=colors,
                    povray_settings=dict(canvas_width=witdth_res, celllinewidth=0, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(config_copy, radius=RADIUS))
                    #camera_type='perspective'
                ).render()
                os.remove(prefix+which+'_{0}_pov.pov'.format(label))
                os.remove(prefix+which+'_{0}_pov.ini'.format(label))

                #top view
                zmin = min([atom.z for atom in config])
                zmax = max([atom.z for atom in config])
                delta_z = zmax - zmin
                transmittances = [0]*len(config) #[1-(atom.z - zmin)/delta_z for atom in config] #linearization

                textures = ['pale'] * Nbulk + ['ase3'] * (len(config)-Nbulk)
                write(prefix+which+'_{0}_top_pov.pov'.format(label), 
                    config, 
                    format='pov',
                    radii = 0.65, 
                    colors=colors,
                    povray_settings=dict(canvas_width=witdth_res, celllinewidth=0, transparent=False, transmittances=transmittances, textures = textures,
                        camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(config_copy, radius=RADIUS))
                    #camera_type='perspective'
                ).render()
                os.remove(prefix+which+'_{0}_top_pov.pov'.format(label))
                os.remove(prefix+which+'_{0}_top_pov.ini'.format(label))

        else:
            if rotations is not None: #use specified rotation
                write(prefix+which+'_{0}_{1}.png'.format(label, rotations), config, rotation=rotations, scale = 100, colors=colors)
            else: #front and top view           
                write(prefix+which+'_{0}.png'.format(label), config, rotation='-10z,-80x', scale = 100, colors=colors) #front view
                write(prefix+which+'_{0}_top.png'.format(label), config, scale = 100, colors=colors) #top view


    if(which=='relax' and index is None):
        os.chdir('..')
        import numpy as np
        from matplotlib import pyplot as plt
        import matplotlib.image as mpimg

        rows_fig = max(int(np.ceil(len(configs)/5)), 1)
        cols_fig = max(int(np.ceil(len(configs)/rows_fig)), 1)
        fig = plt.figure(figsize=(1.5 * cols_fig, 5 * rows_fig / 3))
        fig.subplots_adjust(wspace=0.001)
        axes = [fig.add_subplot(rows_fig,cols_fig,i) for i in range(1,len(configs) + 1)]
    
        energies = [read_energy(file, *E_slab_mol) for file in pw_list if file not in uncompleted]

        for i, conf in enumerate(configs):
            label = labels[i]
            img = mpimg.imread('relax_'+images_dirname+'/'+prefix+which+'_{0}{1}.png'.format(label, '_pov' if povray else ''))
            axes[i].imshow(img)
            axes[i].set_title('{0:.3f} eV'.format(energies[i]), fontsize = 7, pad=1)
            axes[i].text(0.018, 0.983, label, bbox=dict(boxstyle='square', linewidth=0.5, fc="w", ec="k"),transform=axes[i].transAxes, 
                fontsize = 4, color="black",horizontalalignment="left", verticalalignment="top")
            axes[i].set_xticks([])
            axes[i].set_yticks([])
        fig.savefig('relax_'+images_dirname+'/'+'relax_overview{0}.png'.format('_pov' if povray else ''), dpi=1500, bbox_inches='tight')

    print('All images saved in {0}.'.format(which+'_'+images_dirname))


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


    from ase.data.colors import jmol_colors
    ATOM_COLORS = jmol_colors.copy()
    for color in USER_COLORS_DEFS:
        ATOM_COLORS[color[0]] = color[1]
    colors = [ATOM_COLORS[atom.number] for atom in configs[0]]


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
                    povray_settings=dict(canvas_width=witdth_res, transparent=False, camera_type='orthographic', camera_dist=50., bondatoms=get_bondpairs(step_copy, radius=RADIUS))
                    #camera_type='perspective'
                ).render()

            os.chdir('..')
            #os.system('convert -delay 20 temp/step_*.png '+pwo_prefix+'relax_{0}_pov.gif'.format(labels[i]))
            os.system('ffmpeg -framerate 10 -i temp/step_%04d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p '+pwo_prefix+'relax_{0}_pov.mp4'.format(labels[i]))
            shutil.rmtree('temp')

    else:
        for i, config in enumerate(configs):
            write(pwo_prefix+'relax_{0}.gif'.format(labels[i]), config, rotation='-10z,-80x', interval=150, scale = 100, colors=colors, save_count=None)

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
            if config_e:
                plt.plot(config_e, '-', label='{0}: {1:.2f}{2} eV'.format(labels[i], config_e[-1], '' if relax_terminated[i] else '*'))
                if (not relax_terminated[i]): plt.plot(len(config_e)-1, config_e[-1], 'x', color='black')
            else:
                print('Config. {0} job has not reached scf convergence. It will be skipped.'.format(i))
            #plt.xlim(xmin=0)
            
    
    plt.title('Relax energies')
    plt.xlabel('step')
    plt.ylabel('energy (eV)')
    plt.grid(linestyle='dotted')
    plt.legend(title="Config, energy")
    plt.savefig(energy_plot_filename, dpi=300, bbox_inches='tight')

    print('plot saved in {0}'.format(energy_plot_filename))