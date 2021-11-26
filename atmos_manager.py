import os


class AtmosManager:

    def __init__(self, teff, logg, feh, vmicro=None, cfe=0., alphafe=0.):
        # Need to decide if we are going to have c and alpha input as flags
        # Or if we want to just input values and just force add set C and alpha
        self.teff = teff
        self.logg = logg
        self.feh = feh
        if vmicro is None:
            self.vmicro = self.vmicro_calc()
        else:
            self.vmicro = vmicro

        self.cfe = cfe
        self.alphafe = alphafe

    def vmicro_calc(self):
        if self.logg < 3.5:
            if self.teff < 5250.:
                vmicro = 1.15 + 2e-4 * (self.teff - 5500.) + 3.95e-7 * (
                    self.teff - 5500.)**2. - 0.13 * (
                    self.logg - 4.0) + 0.13 * (self.logg - 4.0) ** 2.
            else:
                vmicro = 1.15 + 2e-4 * (5250. - 5500.) + 3.95e-7 * (
                    5250. - 5500.)**2. - 0.13 * (
                    self.logg - 4.0) + 0.13 * (self.logg - 4.0) ** 2.
        else:
            vmicro = 0.94 + 2.2E-5 * (self.teff - 5500) - 0.5E-7 * (
                self.teff - 5500)**2 - 0.1 * (
                self.logg - 4.0) + 0.04 * (
                self.logg - 4.0)**2 - 0.37 * self.feh - 0.07 * self.feh**2

        return vmicro

    def interp_atmos(self, star=None, file_format='ts'):
        if file_format.lower() in ['moog', 'm']:
            format_code = 'MOOG'
        elif file_format.lower() in ['turbospectrum',
                                     'turbo', 'turbospec', 'ts']:
            format_code = 'TurboSpectrum'

        if star is None:
            self.filename = (f'{self.teff:.0f}g{self.logg:.2f}m1.0z'
                             f'{self.feh:.2f}.int')
        else:
            self.filename = (f'{self.teff:.0f}g{self.logg:.2f}m1.0z'
                             f'{self.feh:.2f}_{star}.int')

        atmos_param_arg = (f'{self.teff:.0f} {self.logg:.2f} {self.feh:.2f} '
                           f'{self.vmicro:.2f} {self.cfe:.2f} '
                           f'{self.alphafe:.2f} {format_code} {self.filename}')

        os.system(f'./load_atmos_param.com {atmos_param_arg}')

    def moog_abund(self, abund_dict):
        # To be implemented
        pass

    def move(self, path):
        # To be implemented
        pass
