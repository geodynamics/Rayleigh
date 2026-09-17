Curl Momentum Equation
====================================================================

===================================================================================================================================================================== ====== ======================================= 
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right]_r`                                                      1301    curl\_v\_grad\_v\_r       
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right]_\theta`                                                 1302    curl\_v\_grad\_v\_theta   
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right]_\phi`                                                   1303    curl\_v\_grad\_v\_phi     
 :math:`\left(\left[\nabla\times\left(\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right)\right]_r\right)^2`                          1304    curl\_v\_grad\_v\_r\_squared       
 :math:`\left(\left[\nabla\times\left(\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right)\right]_\theta\right)^2`                     1305    curl\_v\_grad\_v\_theta\_squared   
 :math:`\left(\left[\nabla\times\left(\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right)\right]_\phi\right)^2`                       1306    curl\_v\_grad\_v\_phi\_squared     
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_r`                                          1307    curl\_vp\_grad\_vm\_r     
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\theta`                                     1308    curl\_vp\_grad\_vm\_theta 
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\phi`                                       1309    curl\_vp\_grad\_vm\_phi   
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_r`                                          1310    curl\_vm\_grad\_vp\_r     
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\theta`                                     1311    curl\_vm\_grad\_vp\_theta 
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\phi`                                       1312    curl\_vm\_grad\_vp\_phi   
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_r`                                                    1313    curl\_vp\_grad\_vp\_r     
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\theta`                                               1314    curl\_vp\_grad\_vp\_theta 
 :math:`\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\phi`                                                 1315    curl\_vp\_grad\_vp\_phi   
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_r`                                1316    curl\_vm\_grad\_vm\_r     
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\theta`                           1317    curl\_vm\_grad\_vm\_theta 
 :math:`\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\phi`                             1318    curl\_vm\_grad\_vm\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\theta`                                                                               1319    curl\_buoyancy\_force\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\phi`                                                                                 1320    curl\_buoyancy\_force\_phi  
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\theta\right)^2`                                                                1321    curl\_buoyancy\_force\_theta\_squared  
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\phi\right)^2`                                                                  1322    curl\_buoyancy\_force\_phi\_squared  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta'\right)\right]_\theta`                                                                              1323    curl\_buoyancy\_pforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta'\right)\right]_\phi`                                                                                1324    curl\_buoyancy\_pforce\_phi  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\overline{\Theta}\right)\right]_\phi`                                                                      1326    curl\_buoyancy\_mforce\_phi  
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_r`                                               1327    curl\_coriolis\_force\_r      
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\theta`                                          1328    curl\_coriolis\_force\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\phi`                                            1329    curl\_coriolis\_force\_phi    
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_r\right)^2`                                1330    curl\_coriolis\_force\_r\_squared      
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\theta\right)^2`                           1331    curl\_coriolis\_force\_theta\_squared  
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\phi\right)^2`                             1332    curl\_coriolis\_force\_phi\_squared    
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_r`                                              1333    curl\_coriolis\_pforce\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_\theta`                                         1334    curl\_coriolis\_pforce\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_\phi`                                           1335    curl\_coriolis\_pforce\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_r`                              1336    curl\_coriolis\_mforce\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_\theta`                         1337    curl\_coriolis\_mforce\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_\phi`                           1338    curl\_coriolis\_mforce\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_r`                                                   1339    curl\_viscous\_force\_r       
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\theta`                                              1340    curl\_viscous\_force\_theta   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\phi`                                                1341    curl\_viscous\_force\_phi     
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_r\right)^2`                                    1342    curl\_viscous\_force\_r\_squared       
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\theta\right)^2`                               1343    curl\_viscous\_force\_theta\_squared   
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\phi\right)^2`                                 1344    curl\_viscous\_force\_phi\_squared     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_r`                                                  1345    curl\_viscous\_pforce\_r      
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_\theta`                                             1346    curl\_viscous\_pforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_\phi`                                               1347    curl\_viscous\_pforce\_phi    
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_r`                                        1351    curl\_viscous\_mforce\_r      
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_\theta`                                   1352    curl\_viscous\_mforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_\phi`                                     1353    curl\_viscous\_mforce\_phi    
 :math:`\left[\nabla\times\nabla P\right]_\theta`                                                                                                                      1358    curl\_pressure\_force\_theta   
 :math:`\left[\nabla\times\nabla P\right]_\phi`                                                                                                                        1359    curl\_pressure\_force\_phi     
 :math:`\left(\left[\nabla\times\nabla P\right]_\theta\right)^2`                                                                                                       1361    curl\_pressure\_force\_theta\_squared   
 :math:`\left(\left[\nabla\times\nabla P\right]_\phi\right)^2`                                                                                                         1362    curl\_pressure\_force\_phi\_squared     
 :math:`\left[\nabla\times\nabla P'\right]_\theta`                                                                                                                     1364    curl\_pressure\_pforce\_theta  
 :math:`\left[\nabla\times\nabla P'\right]_\phi`                                                                                                                       1365    curl\_pressure\_pforce\_phi    
 :math:`\left[\nabla\times\nabla \overline{P}\right]_\phi`                                                                                                             1368    curl\_pressure\_mforce\_phi    
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_r`                              1369    curl\_j\_cross\_b\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\theta`                         1370    curl\_j\_cross\_b\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\phi`                           1371    curl\_j\_cross\_b\_phi   
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_r\right)^2`               1372    curl\_j\_cross\_b\_r\_squared     
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\theta\right)^2`          1373    curl\_j\_cross\_b\_theta\_squared 
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\phi\right)^2`            1374    curl\_j\_cross\_b\_phi\_squared   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_r`                  1375    curl\_jp\_cross\_bm\_r       
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_\theta`             1376    curl\_jp\_cross\_bm\_theta   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_\phi`               1377    curl\_jp\_cross\_bm\_phi     
 :math:`\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_r`                              1378    curl\_jm\_cross\_bp\_r     
 :math:`\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_\theta`                         1379    curl\_jm\_cross\_bp\_theta 
 :math:`\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_\phi`                           1380    curl\_jm\_cross\_bp\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_r`        1381    curl\_jm\_cross\_bm\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_\theta`   1382    curl\_jm\_cross\_bm\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_\phi`     1383    curl\_jm\_cross\_bm\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_r`                            1384    curl\_jp\_cross\_bp\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_\theta`                       1385    curl\_jp\_cross\_bp\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_\phi`                         1386    curl\_jp\_cross\_bp\_phi   
 :math:`\left|\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right|`                                                                                      1387    curl\_buoyancy\_force\_abs 
 :math:`\left|\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right|`                                                 1388    curl\_coriolis\_force\_abs 
 :math:`\left|\nabla\times\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right|`                                                        1389    curl\_v\_grad\_v\_abs       
 :math:`\left|\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right|`                                1390    curl\_j\_cross\_b\_abs      
 :math:`\left|\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right|`                                                     1391    curl\_viscous\_force\_abs  
===================================================================================================================================================================== ====== ======================================= 
