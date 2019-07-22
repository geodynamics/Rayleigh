Induction Equation
====================================================================

======================================================================================================================================= ====== =========================== 
 :math:`\left[\boldsymbol{B}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_r`                                                             1601    induct\_shear\_r          
 :math:`-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v} \right)B_r`                                                                        1602    induct\_comp\_r           
 :math:`-\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{B}\right]_r`                                                            1603    induct\_advec\_r          
 :math:`\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v}\times\boldsymbol{B}\right)\right]_r`                                         1604    induct\_r                
 :math:`-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B}\right)\right]_r`                   1605    induct\_diff\_r           
 :math:`\left[\boldsymbol{B}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\theta`                                                        1606    induct\_shear\_theta      
 :math:`-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v} \right)B_\theta`                                                                   1607    induct\_comp\_theta       
 :math:`-\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{B}\right]_\theta`                                                       1608    induct\_advec\_theta      
 :math:`\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v}\times\boldsymbol{B}\right)\right]_\theta`                                    1609    induct\_theta            
 :math:`-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B}\right)\right]_\theta`              1610    induct\_diff\_theta       
 :math:`\left[\boldsymbol{B}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\phi`                                                          1611    induct\_shear\_phi        
 :math:`-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v} \right)B_\phi`                                                                     1612    induct\_comp\_phi         
 :math:`-\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{B}\right]_\phi`                                                         1613    induct\_advec\_phi        
 :math:`\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v}\times\boldsymbol{B}\right)\right]_\phi`                                      1614    induct\_phi              
 :math:`-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B}\right)\right]_\phi`                1615    induct\_diff\_phi         
 :math:`\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_r`                                       1616    induct\_shear\_vmbm\_r     
 :math:`-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)\overline{B_r}`                                       1617    induct\_comp\_vmbm\_r      
 :math:`-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_r`                                      1618    induct\_advec\_vmbm\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\overline{\boldsymbol{B}}\right)\right]_r`                   1619    induct\_vmbm\_r           
 :math:`-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\right]_r`        1620    induct\_diff\_bm\_r        
 :math:`\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\theta`                                  1621    induct\_shear\_vmbm\_theta 
 :math:`-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)\overline{B_\theta}`                                  1622    induct\_comp\_vmbm\_theta  
 :math:`-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_\theta`                                 1623    induct\_advec\_vmbm\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\overline{\boldsymbol{B}}\right)\right]_\theta`              1624    induct\_vmbm\_theta       
 :math:`-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\right]_\theta`   1625    induct\_diff\_bm\_theta    
 :math:`\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi`                                    1626    induct\_shear\_vmbm\_phi   
 :math:`-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)\overline{B_\phi}`                                    1627    induct\_comp\_vmbm\_phi    
 :math:`-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_\phi`                                   1628    induct\_advec\_vmbm\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\overline{\boldsymbol{B}}\right)\right]_\phi`                1629    induct\_vmbm\_phi         
 :math:`-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\right]_\phi`     1630    induct\_diff\_bm\_phi      
 :math:`\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_r`                                                 1631    induct\_shear\_vmbp\_r     
 :math:`-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)B'_r`                                                 1632    induct\_comp\_vmbp\_r      
 :math:`-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_r`                                                1633    induct\_advec\_vmbp\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\boldsymbol{B'}\right)\right]_r`                             1634    induct\_vmbp\_r           
 :math:`-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\right]_r`                  1635    induct\_diff\_bp\_r        
 :math:`\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\theta`                                            1636    induct\_shear\_vmbp\_theta 
 :math:`-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)B'_\theta`                                            1637    induct\_comp\_vmbp\_theta  
 :math:`-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_\theta`                                           1638    induct\_advec\_vmbp\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\boldsymbol{B'}\right)\right]_\theta`                        1639    induct\_vmbp\_theta       
 :math:`-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\right]_\theta`             1640    induct\_diff\_bp\_theta    
 :math:`\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi`                                              1641    induct\_shear\_vmbp\_phi   
 :math:`-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)B'_\phi`                                              1642    induct\_comp\_vmbp\_phi    
 :math:`-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_\phi`                                             1643    induct\_advec\_vmbp\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\boldsymbol{B'}\right)\right]_\phi`                          1644    induct\_vmbp\_phi         
 :math:`-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\right]_\phi`               1645    induct\_diff\_bp\_phi      
 :math:`\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_r`                                                 1646    induct\_shear\_vpbm\_r     
 :math:`-\left(\overline{\boldsymbol{\nabla}}\cdot\boldsymbol{v'} \right)\overline{B_r}`                                                 1647    induct\_comp\_vpbm\_r      
 :math:`-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_r`                                                1648    induct\_advec\_vpbm\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\overline{\boldsymbol{B}}\right)\right]_r`                             1649    induct\_vpbm\_r           
 :math:`\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\theta`                                            1650    induct\_shear\_vpbm\_theta 
 :math:`-\left(\overline{\boldsymbol{\nabla}}\cdot\boldsymbol{v'} \right)\overline{B_\theta}`                                            1651    induct\_comp\_vpbm\_theta  
 :math:`-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_\theta`                                           1652    induct\_advec\_vpbm\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\overline{\boldsymbol{B}}\right)\right]_\theta`                        1653    induct\_vpbm\_theta       
 :math:`\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi`                                              1654    induct\_shear\_vpbm\_phi   
 :math:`-\left(\overline{\boldsymbol{\nabla}}\cdot\boldsymbol{v'} \right)\overline{B_\phi}`                                              1655    induct\_comp\_vpbm\_phi    
 :math:`-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_\phi`                                             1656    induct\_advec\_vpbm\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\overline{\boldsymbol{B}}\right)\right]_\phi`                          1657    induct\_vpbm\_phi         
 :math:`\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_r`                                                           1658    induct\_shear\_vpbp\_r     
 :math:`-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v'} \right)B'_r`                                                                      1659    induct\_comp\_vpbp\_r      
 :math:`-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_r`                                                          1660    induct\_advec\_vpbp\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\right]_r`                                       1661    induct\_vpbp\_r           
 :math:`\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\theta`                                                      1662    induct\_shear\_vpbp\_theta 
 :math:`-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v'} \right)B'_\theta`                                                                 1663    induct\_comp\_vpbp\_theta  
 :math:`-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_\theta`                                                     1664    induct\_advec\_vpbp\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\right]_\theta`                                  1665    induct\_vpbp\_theta       
 :math:`\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi`                                                        1666    induct\_shear\_vpbp\_phi   
 :math:`-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v'} \right)B'_\phi`                                                                   1667    induct\_comp\_vpbp\_phi    
 :math:`-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_\phi`                                                       1668    induct\_advec\_vpbp\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\right]_\phi`                                    1669    induct\_vpbp\_phi         
======================================================================================================================================= ====== =========================== 
