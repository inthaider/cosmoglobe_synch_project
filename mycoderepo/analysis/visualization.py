

from datetime import datetime
from pathlib import Path

import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
# mpl.rcParams['figure.titlesize'] = 12
# mpl.rcParams['axes.titlesize'] = 12
# mpl.rcParams['font.size'] = 8
# mpl.rcParams['figure.dpi'] = 300
# plt.rcParams['savefig.facecolor']='white'

now = datetime.now()



# data_path = Path("../data/")
# LWA1_source_path = Path(data_path/"LWA1_source/")
# commander_in_path = Path(data_path/"commander_in/")

# figs_path.mkdir(parents=True, exist_ok=True)
# LWA1_map_images_path.mkdir(parents=True, exist_ok=True)
# LWA1_source_path.mkdir(parents=True, exist_ok=True)
# commander_in_path.mkdir(parents=True, exist_ok=True)


visualization_path = Path("../visualization/")
figs_path = Path(visualization_path/"figs/")
LWA1_map_images_path = Path(visualization_path/"map-images/LWA1_source/")
figs_path.mkdir(parents=True, exist_ok=True)
LWA1_map_images_path.mkdir(parents=True, exist_ok=True)


#
# Display the maps to check if they're okay
#
def plot_maps(maps, cbrange_K)=None:
    """Plot the maps in a 3x3 grid.

    Args:
        maps (dict): Dictionary of maps to plot.
        cbrange_K (dict): Dictionary of colorbar ranges for each map. Defaults to None.

    Returns:
        None: Displays the maps in a 3x3 grid.

    Notes:
        The maps are plotted in a 3x3 grid, with the colorbar range for each map
        set to the values in the cbrange_K dictionary.

    """
    fig, axs = plt.subplots(3, 3, figsize=(7.3, 4.8), dpi=300)
    mpl.rcParams['figure.titlesize'] = 10
    mpl.rcParams['axes.titlesize'] = 6
    mpl.rcParams['font.size'] = 4
    # mpl.rcParams['legend.title_fontsize'] = 8
    # mpl.rcParams['legend.fontsize'] = 8
    # mpl.rcParams['xtick.labelsize'] = 8
    # mpl.rcParams['ytick.labelsize'] = 8
    # mpl.rcParams['axes.labelsize'] = 8

    # plt.tight_layout()
    count = 0
    for ax,freq,map in zip(axs.flat, maps.keys(), maps.values()):
        # print(ax)
        plt.axes(ax)
        mrgn = 1.0
        margins = (mrgn,mrgn,mrgn,mrgn)
        
        if cbrange_K is None:
            hp.mollview(map, title=f"{freq} MHz", xsize=2000, cmap='jet', unit='K', hold=True)
        else:
            hp.mollview(map, title=f"{freq} MHz", min=cbrange_K[freq][0], max=cbrange_K[freq][1], xsize=2000, cmap='jet', unit='K', hold=True)
        # hp.mollview(map, title=f"{freq} MHz", min=0.0, max=8000.0, xsize=2000, cmap='jet', unit='K', margins=margins, hold=True)
        hp.graticule(lw=0.5)

        # hp.mollview(celest_lwa180/1000, min=0.9, max=4.5, coord=['C', 'G'], cmap='gist_heat', unit='kK')
        # hp.mollview(map_lwa180/1000, min=0.9, max=4.5, coord=['C', 'G'], cmap='gist_heat')
        # hp.mollview(map_lwa180, coord=['C', 'G'])
        # hp.mollview(map_lwa180,  coord=['C', 'G'], norm="hist")

        ### Change font size of colorbar units text in hp.mollview ###
        f = plt.gcf().get_children() # accessing the current figure...
        # print(f)
        # print(len(f))
        # print(f[-1]) # ... then the colorbar's elements
        CbAx = f[-1].get_children() # ... then the colorbar's elements
        # print(CbAx)
        # print(len(CbAx))
        coord_text_obj = CbAx[2] # [1] corresponds to the particular label of the
                                            # colorbar, i.e. "Field value" in this case
        coord_text_obj.set_fontsize(6)
        coord_text_obj.set_y(-4.0)
        # plt.show()
        
        #
        # Save just the single frequency map (subplot)
        #
        # Pad the saved area by 10% in the x-direction and 10% in the y-direction
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

        LWA1_map_image_dir = Path(LWA1_map_images_path/f'LWA1-00{freq}/')
        LWA1_map_image_dir.mkdir(parents=True, exist_ok=True)
        LWA1_map_image_file = f"LWA1-00{freq}_map.png"
        fig.savefig(LWA1_map_image_dir/LWA1_map_image_file, bbox_inches=extent.expanded(1.01, 1.20), dpi=300)
        # fig.savefig(LWA1_map_images_path/f'LWA1-00{freq}/LWA1-00{freq}_map.png', bbox_inches=extent)

        count = count + 1

    fig.suptitle("LWA1 Survey Maps", y=0.95)
    plt.tight_layout()
    # fig.subplots_adjust(top=0.95)
    # fig.tight_layout()

    # 
    # Save the whole figure as a png image
    # 
    LWA1_allmaps_image_file = f"LWA1_allmaps.png"
    plt.savefig(LWA1_map_images_path/LWA1_allmaps_image_file, dpi=300, bbox_inches='tight')


# Create a function to plot the percentage error/noise maps in a 3x3 grid the same way as above
def plot_noise_maps(pct_errs):
    """Plot the noise maps in a 3x3 grid.
    
    Args:
        pct_errs (dict): Dictionary of percent noise maps to plot.
        
    Returns:
        None: Displays the maps in a 3x3 grid.

    Notes:
        The maps are plotted in a 3x3 grid, with the colorbar range for each map
        set to the values in the cbrange_K dictionary.

    """
    fig, axs = plt.subplots(3, 3, figsize=(7.3, 4.8), dpi=300)
    mpl.rcParams['figure.titlesize'] = 10
    mpl.rcParams['axes.titlesize'] = 6
    mpl.rcParams['font.size'] = 4

    for ax,freq,pct_err in zip(axs.flat, pct_errs.keys(), pct_errs.values()):
        plt.axes(ax)
        hp.mollview(pct_err, title=f"{freq} MHz", min=0.0, max=25.0, xsize=2000, cmap='gist_heat', unit='% Error', hold=True)
        # Draw graticule lines using hp.graticule() with line-width of 0.5:
        hp.graticule(lw=0.5)

        

        ### Change font size of colorbar units text in hp.mollview ###
        f = plt.gcf().get_children() # accessing the current figure...
        CbAx = f[-1].get_children() # ... then the colorbar's elements
        coord_text_obj = CbAx[2] # [1] corresponds to the particular label of the
                                            # colorbar, i.e. "Field value" in this case
        coord_text_obj.set_fontsize(6)
        coord_text_obj.set_y(-4.0)

        # Save just the single frequency map (subplot)
        # Pad the saved area by 10% in the x-direction and 10% in the y-direction
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        LWA1_pcterr_image_dir = Path(LWA1_map_images_path/f'LWA1-00{freq}/')
        LWA1_pcterr_image_dir.mkdir(parents=True, exist_ok=True)
        LWA1_pcterr_image_file = f"LWA1-00{freq}_pcterr.png"
        fig.savefig(LWA1_pcterr_image_dir/LWA1_pcterr_image_file, bbox_inches=extent.expanded(1.01, 1.20), dpi=300)
    
    fig.suptitle("LWA1 Percent Noise Maps", y=0.95)
    plt.tight_layout()

    # 
    # Save the whole figure as a png image
    # 
    LWA1_allpcterrs_image_file = f"LWA1_allpcterrs.png"
    plt.savefig(LWA1_map_images_path/LWA1_allpcterrs_image_file, dpi=300, bbox_inches='tight')

# plot_noise_maps(pct_errs_lwa1)
    