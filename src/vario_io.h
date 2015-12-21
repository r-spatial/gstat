double sem_cov_ab(VARIOGRAM *v, DPOINT *a, DPOINT *b, int sem);
/* covariance: */
#define COVARIANCE(v,a,b) ((IS_POINT(a) && IS_POINT(b) && !gl_longlat) ? \
	(get_covariance(v,a->x - b->x,a->y - b->y, a->z - b->z)) : \
	sem_cov_ab(v,a,b,0))
/* generalized covariance: */
#define GCV(v,a,b) ((IS_POINT(a) && IS_POINT(b) && !gl_longlat) ? \
	(v->max_val - get_semivariance(v,a->x - b->x,a->y - b->y, a->z - b->z)) : \
	(v->max_val - sem_cov_ab(v,a,b,1)))
/* 
 * CME is the measurement error-adjustment to GCV or COVARIANCE:
 * see Cressie, Statistics for Spatial Data, revised ed. 1993,
 * eq. 3.2.25-3.2.27, and page 379 
 */
#define CME(v,a,b,dist) ((IS_POINT(a) && IS_POINT(b) && \
	(a == b || dist(a, b) == 0.0)) ? v->measurement_error : 0.0)
#define GCV0(v,a,b,dist)        (GCV(v,a,b)        - CME(v,a,b,dist))
#define COVARIANCE0(v,a,b,dist) (COVARIANCE(v,a,b) - CME(v,a,b,dist))
