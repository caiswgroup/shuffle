all:

	##################################################
	###               kissat_mab                   ###
	##################################################
	cd kissat_mab && ./configure
	+ $(MAKE) -C kissat_mab

	##################################################
	###                 PaInleSS                   ###
	##################################################
	+ $(MAKE) -C painless-src
	mv painless-src/painless painless-mcomsps

clean:
	##################################################
	###                 PaInleSS                   ###
	##################################################
	+ $(MAKE) clean -C painless-src
	rm -f painless-mcomsps
