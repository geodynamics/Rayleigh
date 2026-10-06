#!/home/pgrad1/2831664s/anaconda3/bin/python

from plot_meridional import animate, VISCOUS_FORCE_ROW, CURL_VISCOUS_FORCE_ROW

if __name__ == '__main__':

    animate(
        'Meridional_Viscous_Force.mp4',
        row_funcs=(VISCOUS_FORCE_ROW, CURL_VISCOUS_FORCE_ROW),
        row_labels=(
            r'$F_V$',
            r'$\nabla\times F_V$'
        ),
        frame_stride=4  # longer wavelength / longer run here: 2000 slices -> 500 frames
    )
