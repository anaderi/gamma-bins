TEV_URL=https://gamma-cat.readthedocs.io/output/gammacat.fits.gz
GEV_URL=https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4yr_catalog/gll_psc_v16.fit
all:	create_dir	add_data
create_dir	: 
	if	not	exist	"data"	mkdir	data
	if	not	exist	"markeddata"	mkdir	markeddata
	if	not	exist	"spectrums"	mkdir	spectrums
add_data	:
	if	not	exist	"data/gammacat.fits.gz"	curl	$(TEV_URL)	--output	data/gammacat.fits.gz
	if	not	exist	"data/gll_psc_v16.fit"	curl	$(GEV_URL)	--output	data/gll_psc_v16.fit