Kinetic Energy Equation
====================================================================

====================================================================================================================================================== ====== ======================== 
 :math:`-c_3\mathrm{f}_1\boldsymbol{v}\cdot\boldsymbol{\nabla}\left(\frac{P}{\mathrm{f}_1}\right)`                                                      1901    press\_work    
 :math:`-c_3\mathrm{f}_1\boldsymbol{v'}\cdot\boldsymbol{\nabla}\left(\frac{P'}{\mathrm{f}_1}\right)`                                                    1902    press\_work\_pp 
 :math:`-c_3\mathrm{f}_1\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\left(\frac{\overline{P}}{\mathrm{f}_1}\right)`                                1903    press\_work\_mm 
 :math:`c_2v_r\mathrm{f}_2\Theta`                                                                                                                       1904    buoy\_work    
 :math:`c_2v'_r\mathrm{f}_2\Theta'`                                                                                                                     1905    buoy\_work\_pp 
 :math:`c_2\overline{v_r}\mathrm{f}_2\overline{\Theta}`                                                                                                 1906    buoy\_work\_mm 
 :math:`c_5\boldsymbol{v}\cdot\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\right]`                                                            1907    visc\_work    
 :math:`c_5\boldsymbol{v'}\cdot\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D'}}\right]`                                                          1908    visc\_work\_pp 
 :math:`c_5\overline{\boldsymbol{v}}\cdot\left[\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal{D}}}\right]`                                      1909    visc\_work\_mm 
 :math:`\mathrm{f}_1\boldsymbol{v}\cdot\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]`                                               1910    advec\_work     
 :math:`\mathrm{f}_1\boldsymbol{v'}\cdot\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]`                                            1911    advec\_work\_ppp 
 :math:`\mathrm{f}_1\overline{\boldsymbol{v}}\cdot\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]`                                  1912    advec\_work\_mpp 
 :math:`\mathrm{f}_1\boldsymbol{v'}\cdot\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]`                                  1913    advec\_work\_pmp 
 :math:`\mathrm{f}_1\boldsymbol{v'}\cdot\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]`                                  1914    advec\_work\_ppm 
 :math:`\mathrm{f}_1\overline{\boldsymbol{v}}\cdot\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]`              1915    advec\_work\_mmm 
 :math:`c_4\boldsymbol{v}\cdot\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B}\right]`                                    1916    mag\_work     
 :math:`c_4\boldsymbol{v'}\cdot\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]`                                 1918    mag\_work\_ppp 
 :math:`c_4\overline{\boldsymbol{v}}\cdot\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]`                       1919    mag\_work\_mpp 
 :math:`c_4\boldsymbol{v'}\cdot\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\boldsymbol{B'}\right]`                       1920    mag\_work\_pmp 
 :math:`c_4\boldsymbol{v'}\cdot\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\overline{\boldsymbol{B}}\right]`                       1921    mag\_work\_ppm 
 :math:`c_4\overline{\boldsymbol{v}}\cdot\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\overline{\boldsymbol{B}}\right]`   1922    mag\_work\_mmm 
 :math:`\frac{1}{2}\mathrm{f}_1v_r\,v^2`                                                                                                                1923    ke\_flux\_radial 
 :math:`\frac{1}{2}\mathrm{f}_1v_\theta\,v^2`                                                                                                           1924    ke\_flux\_theta  
 :math:`\frac{1}{2}\mathrm{f}_1v_\phi\,v^2`                                                                                                             1925    ke\_flux\_phi    
 :math:`\frac{1}{2}\mathrm{f}_1\,\overline{v_r}\,\overline{v}^2`                                                                                        1926    mke\_mflux\_radial 
 :math:`\frac{1}{2}\mathrm{f}_1\,\overline{v_\theta}\,\overline{v}^2`                                                                                   1927    mke\_mflux\_theta  
 :math:`\frac{1}{2}\mathrm{f}_1\,\overline{v_\phi}\,\overline{v}^2`                                                                                     1928    mke\_mflux\_phi    
 :math:`\frac{1}{2}\mathrm{f}_1\,\overline{v_r}\,{v'}^2`                                                                                                1929    pke\_mflux\_radial  
 :math:`\frac{1}{2}\mathrm{f}_1\,\overline{v_\theta}\,{v'}^2`                                                                                           1930    pke\_mflux\_theta   
 :math:`\frac{1}{2}\mathrm{f}_1\,\overline{v_\phi}\,{v'}^2`                                                                                             1931    pke\_mflux\_phi     
 :math:`\frac{1}{2}\mathrm{f}_1\,v'_r\,{v'}^2`                                                                                                          1932    pke\_pflux\_radial  
 :math:`\frac{1}{2}\mathrm{f}_1\,v'_\theta\,{v'}^2`                                                                                                     1933    pke\_pflux\_theta   
 :math:`\frac{1}{2}\mathrm{f}_1\,v'_\phi\,{v'}^2`                                                                                                       1934    pke\_pflux\_phi     
 :math:`c_5\left[\boldsymbol{v}\cdot\boldsymbol{\mathcal{D}} \right]_r`                                                                                 1935    visc\_flux\_r     
 :math:`c_5\left[\boldsymbol{v}\cdot\boldsymbol{\mathcal{D}} \right]_\theta`                                                                            1936    visc\_flux\_theta 
 :math:`c_5\left[\boldsymbol{v}\cdot\boldsymbol{\mathcal{D}} \right]_\phi`                                                                              1937    visc\_flux\_phi   
 :math:`c_5\left[\boldsymbol{v'}\cdot\boldsymbol{\mathcal{D'}} \right]_r`                                                                               1938    visc\_fluxpp\_r     
 :math:`c_5\left[\boldsymbol{v'}\cdot\boldsymbol{\mathcal{D'}} \right]_\theta`                                                                          1939    visc\_fluxpp\_theta 
 :math:`c_5\left[\boldsymbol{v'}\cdot\boldsymbol{\mathcal{D'}} \right]_\phi`                                                                            1940    visc\_fluxpp\_phi   
 :math:`c_5\left[\boldsymbol{\overline{v}}\cdot\boldsymbol{\overline{\mathcal{D}}} \right]_r`                                                           1941    visc\_fluxmm\_r     
 :math:`c_5\left[\boldsymbol{\overline{v}}\cdot\boldsymbol{\overline{\mathcal{D}}} \right]_\theta`                                                      1942    visc\_fluxmm\_theta 
 :math:`c_5\left[\boldsymbol{\overline{v}}\cdot\boldsymbol{\overline{\mathcal{D}}} \right]_\phi`                                                        1943    visc\_fluxmm\_phi   
 :math:`-c_3v_r P`                                                                                                                                      1944    press\_flux\_r     
 :math:`-c_3v_\theta P`                                                                                                                                 1945    press\_flux\_theta 
 :math:`-c_3v_\phi P`                                                                                                                                   1946    press\_flux\_phi   
 :math:`-c_3v'_r P'`                                                                                                                                    1947    press\_fluxpp\_r     
 :math:`-c_3v'_\theta P'`                                                                                                                               1948    press\_fluxpp\_theta 
 :math:`-c_3v'_\phi P'`                                                                                                                                 1949    press\_fluxpp\_phi   
 :math:`-c_3\overline{v_r}\, \overline{P}`                                                                                                              1950    press\_fluxmm\_r     
 :math:`-c_3\overline{v_\theta}\, \overline{P}`                                                                                                         1951    press\_fluxmm\_theta 
 :math:`-c_3\overline{v_\phi}\, \overline{P}`                                                                                                           1952    press\_fluxmm\_phi   
 :math:`--`                                                                                                                                             1953    production\_shear\_ke  
 :math:`--`                                                                                                                                             1954    production\_shear\_pke 
 :math:`--`                                                                                                                                             1955    production\_shear\_mke 
====================================================================================================================================================== ====== ======================== 
