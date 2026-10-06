#!/home/pgrad1/2831664s/anaconda3/bin/python

from plot_meridional import animate, force_residual, curl_force_residual

if __name__ == '__main__':

    animate(
        'Meridional_Residuals.mp4',
        row_funcs=(force_residual, curl_force_residual),
        row_labels=(
            'Force residual',
            'Curl-of-force residual'
        ),
        frame_stride=4  # longer wavelength / longer run here: 2000 slices -> 500 frames
    )
