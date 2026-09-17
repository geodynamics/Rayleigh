#!/home/pgrad1/2831664s/anaconda3/bin/python

from plot_meridional import animate, J_CROSS_B_ROW, CURL_J_CROSS_B_ROW

if __name__ == '__main__':

    animate(
        'Meridional_J_Cross_B.mp4',
        row_funcs=(J_CROSS_B_ROW, CURL_J_CROSS_B_ROW),
        row_labels=(
            r'$\mathbf{J}\times\mathbf{B}$',
            r'$\nabla\times(\mathbf{J}\times\mathbf{B})$'
        )
    )
