ANT = ant
DOWNLOAD_FROM_WEB = ./download_svn.py

CODEDIR = deploy

all: ${CODEDIR}
	${ANT}

$(CODEDIR):
ifdef DOWNLOAD_CODES
	$(RM) -Rf deploy
	$(DOWNLOAD_FROM_WEB)
else
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. Ibis will not be downloaded and build."
	@echo "If you do want Ibis, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif

clean:
	$(ANT) clean
	
distclean:
	$(ANT) distclean

