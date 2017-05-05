MODULE_NAME=test_ii_shadow_build
TEST_ANALYSIS=Test_ii_shadow_build


test:
	(\
		mkdir -pv tmp-work; cd tmp-work;  \
	 	source /home/isdc/savchenk/osa11_deployment/deployment/init.sh && \
	 	PYTHONPATH=..:$(PYTHONPATH) rundda.py $(TEST_ANALYSIS) -m ddosa -m git://ddosa11 -m $(MODULE_NAME) -a 'ddosa.ScWData(input_scwid="066500220010.001")'\
	)
