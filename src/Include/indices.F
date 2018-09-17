#define RHSP wsp%p3b
#define DO_IDX Do t = my_theta%min, my_theta%max;    Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define END_DO enddo; enddo; enddo
#define IDX k,r,t
#define FIELDSP wsp%p3a
#define IDXM r,t
#define DO_IDX2 Do mp = my_mp%min, my_mp%max; m = m_values(mp); Do imi = 1, 2; Do r = my_r%min, my_r%max
#define IDX2 m:l_max,r,imi
#define SBUFFA wsp%s2a(mp)%data
#define SBUFFB wsp%s2b(mp)%data
#define ASBUFFA cobuffer%s2a(mp)%data
#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define PSI k,r,t
#define PSI2 r,t
#define DDBUFF d2buffer%p3a

