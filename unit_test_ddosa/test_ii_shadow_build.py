import ddosa
from astropy.io import fits

class BinEventsImage(ddosa.BinEventsImage):
    input_osatools=ddosa.get_OSA_tools("ii_shadow_build")

class Test_sane_shadows_ii_shadow_build(ddosa.DataAnalysis):
    input_ii_shadow_build=BinEventsImage

    def main(self):
        f_e=fits.open(self.input_ii_shadow_build.shadow_efficiency.get_path())
        f_d=fits.open(self.input_ii_shadow_build.shadow_detector.get_path())
        print "opened",f_e,f_d

class Test_noise_ii_shadow_build(ddosa.DataAnalysis):
    input_ii_shadow_build=BinEventsImage

    def main(self):
        pass
        #raise Exception("fail test!")

class Test_ii_shadow_build(ddosa.DataAnalysis):
    input_t1=Test_sane_shadows_ii_shadow_build
    input_t2=Test_noise_ii_shadow_build

