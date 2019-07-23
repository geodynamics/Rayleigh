Momentum Equation
====================================================================

=============================================================================================================================== ====== ========================== 
 :math:`\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_r`                                         1201    v\_grad\_v\_r       
 :math:`\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\theta`                                    1202    v\_grad\_v\_theta   
 :math:`\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\phi`                                      1203    v\_grad\_v\_phi     
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_r`                             1204    vp\_grad\_vm\_r     
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\theta`                        1205    vp\_grad\_vm\_theta 
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi`                          1206    vp\_grad\_vm\_phi   
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_r`                             1207    vm\_grad\_vp\_r     
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\theta`                        1208    vm\_grad\_vp\_theta 
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi`                          1209    vm\_grad\_vp\_phi   
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_r`                                       1210    vp\_grad\_vp\_r     
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\theta`                                  1211    vp\_grad\_vp\_theta 
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi`                                    1212    vp\_grad\_vp\_phi   
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_r`                   1213    vm\_grad\_vm\_r     
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\theta`              1214    vm\_grad\_vm\_theta 
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi`                1215    vm\_grad\_vm\_phi   
 :math:`c_2\mathrm{f}_2\Theta`                                                                                                   1216    buoyancy\_force  
 :math:`c_2\mathrm{f}_2\Theta'`                                                                                                  1217    buoyancy\_pforce 
 :math:`c_2\mathrm{f}_2\overline{\Theta}`                                                                                        1218    buoyancy\_mforce 
 :math:`-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v}\right]_r`                                                 1219    Coriolis\_Force\_r      
 :math:`-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v}\right]_\theta`                                            1220    Coriolis\_Force\_theta  
 :math:`-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v}\right]_\phi`                                              1221    Coriolis\_Force\_phi    
 :math:`-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v'}\right]_r`                                                1222    Coriolis\_pForce\_r     
 :math:`-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v'}\right]_\theta`                                           1223    Coriolis\_pForce\_theta 
 :math:`-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v'}\right]_\phi`                                             1224    Coriolis\_pForce\_phi   
 :math:`-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\overline{\boldsymbol{v}}\right]_r`                                      1225    Coriolis\_mForce\_r     
 :math:`-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\overline{\boldsymbol{v}}\right]_\theta`                                 1226    Coriolis\_mForce\_theta 
 :math:`-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\overline{\boldsymbol{v}}\right]_\phi`                                   1227    Coriolis\_mForce\_phi   
 :math:`c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\right]_r`                                                      1228    viscous\_Force\_r       
 :math:`c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\right]_\theta`                                                 1229    viscous\_Force\_theta   
 :math:`c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\right]_\phi`                                                   1230    viscous\_Force\_phi     
 :math:`c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D'}}\right]_r`                                                     1231    viscous\_pForce\_r      
 :math:`c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D'}}\right]_\theta`                                                1232    viscous\_pForce\_theta  
 :math:`c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D'}}\right]_\phi`                                                  1233    viscous\_pForce\_phi    
 :math:`c_5\left[\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal{D}}}\right]_r`                                           1234    viscous\_mForce\_r      
 :math:`c_5\left[\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal{D}}}\right]_\theta`                                      1235    viscous\_mForce\_theta  
 :math:`c_5\left[\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal{D}}}\right]_\phi`                                        1236    viscous\_mForce\_phi    
 :math:`-c_3\mathrm{f}_1\frac{\partial}{\partial r}\left(\frac{P}{\mathrm{f}_1} \right)`                                         1237    pressure\_Force\_r       
 :math:`-c_3\frac{1}{r}\frac{\partial P}{\partial \theta}`                                                                       1238    pressure\_Force\_theta   
 :math:`-c_3\frac{1}{r\mathrm{sin}\theta}\frac{\partial P}{\partial \phi}`                                                       1239    pressure\_Force\_phi     
 :math:`-c_3\mathrm{f}_1\frac{\partial}{\partial r}\left(\frac{P'}{\mathrm{f}_1} \right)`                                        1240    pressure\_pForce\_r      
 :math:`-c_3\frac{1}{r}\frac{\partial P'}{\partial \theta}`                                                                      1241    pressure\_pForce\_theta  
 :math:`-c_3\frac{1}{r\mathrm{sin}\theta}\frac{\partial P'}{\partial \phi}`                                                      1242    pressure\_pForce\_phi    
 :math:`-c_3\mathrm{f}_1\frac{\partial}{\partial r}\left(\frac{\overline{P}}{\mathrm{f}_1} \right)`                              1243    pressure\_mForce\_r      
 :math:`-c_3\frac{1}{r}\frac{\partial \overline{P}}{\partial \theta}`                                                            1244    pressure\_mForce\_theta  
 :math:`-c_3\frac{1}{r\mathrm{sin}\theta}\frac{\partial \overline{P}}{\partial \phi}`                                            1245    pressure\_mForce\_phi    
 :math:`c_2\mathrm{f}_2\Theta_{00}`                                                                                              1246    buoyancy\_force\_ell0 
 :math:`-c_3\mathrm{f}_1\frac{\partial}{\partial r}\left(\frac{P_{00}}{\mathrm{f}_1} \right)`                                    1247    pressure\_force\_ell0\_r 
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B}\right]_r`                              1248    j\_cross\_b\_r       
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B}\right]_\theta`                         1249    j\_cross\_b\_theta   
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B}\right]_\phi`                           1250    j\_cross\_b\_phi     
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\overline{\boldsymbol{B}}\right]_r`                  1251    jp\_cross\_bm\_r     
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\overline{\boldsymbol{B}}\right]_\theta`             1252    jp\_cross\_bm\_theta 
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\overline{\boldsymbol{B}}\right]_\phi`               1253    jp\_cross\_bm\_phi   
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\boldsymbol{B'}\right]_r`                  1254    jm\_cross\_bp\_r     
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\boldsymbol{B'}\right]_\theta`             1255    jm\_cross\_bp\_theta 
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\boldsymbol{B'}\right]_\phi`               1256    jm\_cross\_bp\_phi   
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\overline{\boldsymbol{B}}\right]_r`        1257    jm\_cross\_bm\_r     
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\overline{\boldsymbol{B}}\right]_\theta`   1258    jm\_cross\_bm\_theta 
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\overline{\boldsymbol{B}}\right]_\phi`     1259    jm\_cross\_bm\_phi   
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_r`                            1260    jp\_cross\_bp\_r     
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_\theta`                       1261    jp\_cross\_bp\_theta 
 :math:`c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_\phi`                         1262    jp\_cross\_bp\_phi   
=============================================================================================================================== ====== ========================== 
