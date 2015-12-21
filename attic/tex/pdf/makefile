#
# TransFig makefile
#

all = badlos.eps \
	boxtest.eps \
	closed.eps \
	difflos.eps \
	dist.eps \
	los.eps \
	open.eps \
	spec_los.eps \
	decis1.eps \
	decis2.eps \
	msteps.eps \
	bes.eps \
	cir.eps \
	cov.eps \
	ell.eps \
	example.eps \
	exp.eps \
	gau.eps \
	lin.eps \
	log.eps \
	map.eps \
	nbh.eps \
	nug.eps \
	pen.eps \
	per.eps \
	powa.eps \
	powb.eps \
	scat.eps \
	sph.eps \
	vgm.eps

default:	$(all)

zn_map.eps:
	sh make.jgr

# translation into ps
%eps:  %fig
	fig2dev -L ps $< > $@

%eps:  %jgr
	jgraph $< > $@

clean:
	rm $(all)
