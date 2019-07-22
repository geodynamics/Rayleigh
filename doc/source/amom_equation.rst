Angular Momentum Equation
====================================================================

=================================================================================================================================================== ====== ====================== 
 :math:`r\,\mathrm{sin}\theta\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi`                                   1801    samom\_advec\_pp 
 :math:`r\,\mathrm{sin}\theta\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi`               1802    samom\_advec\_mm 
 :math:`-c_1\mathrm{f}_1r\mathrm{sin}\theta \left(\mathrm{cos}\theta\, v_\theta + \mathrm{sin}\theta\, v_r\right)`                                   1803    samom\_coriolis 
 :math:`r\,\mathrm{sin}\theta\left[\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal{D}}}\right]_\phi`                                          1804    samom\_diffusion 
 :math:`r\,\mathrm{sin}\theta c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\overline{\boldsymbol{B}}\right]_\phi`   1805    samom\_lorentz\_mm   
 :math:`r\,\mathrm{sin}\theta c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_\phi`                       1806    samom\_lorentz\_pp   
 :math:`\mathrm{f}_1r\mathrm{sin}\theta v'_r v'_\phi`                                                                                                1807    famom\_fluct\_r     
 :math:`\mathrm{f}_1r\mathrm{sin}\theta v'_\theta v'_\phi`                                                                                           1808    famom\_fluct\_theta 
 :math:`\mathrm{f}_1r\mathrm{sin}\theta\, \overline{v_r}\,\, \overline{v_\phi}`                                                                      1809    famom\_dr\_r        
 :math:`\mathrm{f}_1r\mathrm{sin}\theta\, \overline{v_\theta}\,\, \overline{v_\phi}`                                                                 1810    famom\_dr\_theta    
 :math:`\frac{c_1}{2}\mathrm{f}_1r^2\mathrm{sin}^2\theta\, \overline{v_r}`                                                                           1811    famom\_mean\_r      
 :math:`\frac{c_1}{2}\mathrm{f}_1r^2\mathrm{sin}^2\theta\, \overline{v_\theta}`                                                                      1812    famom\_mean\_theta  
 :math:`\mathrm{f}_1\nu\mathrm{sin}\theta\left(v_\phi-r\frac{\partial\,\overline{v_\phi}}{\partial r}   \right)`                                     1813    famom\_diff\_r       
 :math:`\mathrm{f}_1\nu\left(\mathrm{cos}\theta\, \overline{v_\phi}-\mathrm{sin}\theta\frac{\partial\, \overline{v_\phi}}{\partial \theta}\right)`   1814    famom\_diff\_theta   
 :math:`-r\mathrm{sin}\theta c_4\,B'_r\,B'_\phi`                                                                                                     1815    famom\_maxstr\_r     
 :math:`-r\mathrm{sin}\theta c_4\,B'_\theta\,B'_\phi`                                                                                                1816    famom\_maxstr\_theta 
 :math:`-r\mathrm{sin}\theta c_4\,\overline{B_r}\,\overline{B_\phi}`                                                                                 1817    famom\_magtor\_r     
 :math:`-r\mathrm{sin}\theta c_4\,\overline{B_\theta}\,\overline{B_\phi}`                                                                            1818    famom\_magtor\_theta 
 :math:`\mathrm{f}_1r\mathrm{sin}\theta v_\phi`                                                                                                      1819    amom\_z 
 :math:`\mathrm{f}_1r(-\mathrm{sin}\theta v_\phi - \mathrm{cos}\phi v_\theta)`                                                                       1820    amom\_x 
 :math:`\mathrm{f}_1r(-\mathrm{cos}\theta v_\phi + \mathrm{cos}\phi v_\theta)`                                                                       1821    amom\_y 
 :math:`\mathrm{f}_1r\mathrm{sin}\theta v_\phi'`                                                                                                     1822    amomp\_z 
 :math:`\mathrm{f}_1r(-\mathrm{sin}\theta v_\phi' - \mathrm{cos}\phi v_\theta')`                                                                     1823    amomp\_x 
 :math:`\mathrm{f}_1r(-\mathrm{cos}\theta v_\phi' + \mathrm{cos}\phi v_\theta')`                                                                     1824    amomp\_y 
 :math:`\mathrm{f}_1r\mathrm{sin}\theta \overline{v_\phi}`                                                                                           1825    amomm\_z 
 :math:`\mathrm{f}_1r(-\mathrm{sin}\theta \overline{v_\phi} - \mathrm{cos}\phi \overline{v_\theta'}`                                                 1826    amomm\_x 
 :math:`\mathrm{f}_1r(-\mathrm{cos}\theta\overline{v_\phi} + \mathrm{cos}\phi \overline{v_\theta})`                                                  1827    amomm\_y 
=================================================================================================================================================== ====== ====================== 
