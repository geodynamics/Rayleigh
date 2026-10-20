Curl Momentum Equation
====================================================================

===================================================================================================================================================================== ====== =============================== 
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right]_r`                                                      1301    curl\_v\_grad\_v\_r       
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right]_\theta`                                                 1302    curl\_v\_grad\_v\_theta   
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right]_\phi`                                                   1303    curl\_v\_grad\_v\_phi     
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_r`                                          1304    curl\_vp\_grad\_vm\_r     
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\theta`                                     1305    curl\_vp\_grad\_vm\_theta 
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\phi`                                       1306    curl\_vp\_grad\_vm\_phi   
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_r`                                          1307    curl\_vm\_grad\_vp\_r     
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\theta`                                     1308    curl\_vm\_grad\_vp\_theta 
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\phi`                                       1309    curl\_vm\_grad\_vp\_phi   
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_r`                                                    1310    curl\_vp\_grad\_vp\_r     
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\theta`                                               1311    curl\_vp\_grad\_vp\_theta 
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\phi`                                                 1312    curl\_vp\_grad\_vp\_phi   
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_r`                                1313    curl\_vm\_grad\_vm\_r     
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\theta`                           1314    curl\_vm\_grad\_vm\_theta 
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\phi`                             1315    curl\_vm\_grad\_vm\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\theta`                                                                               1316    curl\_buoyancy\_force\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\phi`                                                                                 1317    curl\_buoyancy\_force\_phi  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta'\right)\right]_\theta`                                                                              1318    curl\_buoyancy\_pforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta'\right)\right]_\phi`                                                                                1319    curl\_buoyancy\_pforce\_phi  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\overline{\Theta}\right)\right]_\phi`                                                                      1321    curl\_buoyancy\_mforce\_phi  
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_r`                                               1322    curl\_coriolis\_force\_r      
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\theta`                                          1323    curl\_coriolis\_force\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\phi`                                            1324    curl\_coriolis\_force\_phi    
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_r`                                              1325    curl\_coriolis\_pforce\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_\theta`                                         1326    curl\_coriolis\_pforce\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_\phi`                                           1327    curl\_coriolis\_pforce\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_r`                              1328    curl\_coriolis\_mforce\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_\theta`                         1329    curl\_coriolis\_mforce\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_\phi`                           1330    curl\_coriolis\_mforce\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_r`                                                   1331    curl\_viscous\_force\_r       
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\theta`                                              1332    curl\_viscous\_force\_theta   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\phi`                                                1333    curl\_viscous\_force\_phi     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_r`                                                  1334    curl\_viscous\_pforce\_r      
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_\theta`                                             1335    curl\_viscous\_pforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_\phi`                                               1336    curl\_viscous\_pforce\_phi    
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_r`                                        1337    curl\_viscous\_mforce\_r      
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_\theta`                                   1338    curl\_viscous\_mforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_\phi`                                     1339    curl\_viscous\_mforce\_phi    
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_3\mathrm{f}_1\boldsymbol{\nabla}\left(\frac{P}{\mathrm{f}_1}\right)\right)\right]_\theta`                              1341    curl\_pressure\_force\_theta   
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_3\mathrm{f}_1\boldsymbol{\nabla}\left(\frac{P}{\mathrm{f}_1}\right)\right)\right]_\phi`                                1342    curl\_pressure\_force\_phi     
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_3\mathrm{f}_1\boldsymbol{\nabla}\left(\frac{P'}{\mathrm{f}_1}\right)\right)\right]_\theta`                             1344    curl\_pressure\_pforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_3\mathrm{f}_1\boldsymbol{\nabla}\left(\frac{P'}{\mathrm{f}_1}\right)\right)\right]_\phi`                               1345    curl\_pressure\_pforce\_phi    
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_3\mathrm{f}_1\boldsymbol{\nabla}\left(\frac{\overline{P}}{\mathrm{f}_1}\right)\right)\right]_\phi`                     1348    curl\_pressure\_mforce\_phi    
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_r`                              1349    curl\_j\_cross\_b\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\theta`                         1350    curl\_j\_cross\_b\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\phi`                           1351    curl\_j\_cross\_b\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_r`                  1352    curl\_jp\_cross\_bm\_r       
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_\theta`             1353    curl\_jp\_cross\_bm\_theta   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_\phi`               1354    curl\_jp\_cross\_bm\_phi     
 :math:`\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_r`                              1355    curl\_jm\_cross\_bp\_r     
 :math:`\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_\theta`                         1356    curl\_jm\_cross\_bp\_theta 
 :math:`\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_\phi`                           1357    curl\_jm\_cross\_bp\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_r`        1358    curl\_jm\_cross\_bm\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_\theta`   1359    curl\_jm\_cross\_bm\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_\phi`     1360    curl\_jm\_cross\_bm\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_r`                            1361    curl\_jp\_cross\_bp\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_\theta`                       1362    curl\_jp\_cross\_bp\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_\phi`                         1363    curl\_jp\_cross\_bp\_phi   
===================================================================================================================================================================== ====== =============================== 
